#!/usr/bin/env python3
"""Export per-cell SpoQC data without editing the input SpatialData."""

import argparse
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import zarr

BINS = 1001
BATCH_PIXELS = 1_048_576
BLOCK_PIXELS = 16_777_216
MAP_PREFIX = "mask_smoothed_raw_output_"


def parquet_parts(path):
    if path.is_file():
        return [path]
    parts = []
    for part in path.glob("*.parquet"):
        match = re.fullmatch(r"part[.-](\d+)\.parquet", part.name)
        if not match:
            raise ValueError(f"Unknown Parquet part naming: {part}")
        parts.append((int(match[1]), part))
    if not parts or len({i for i, _ in parts}) != len(parts):
        raise ValueError(f"No unique numbered Parquet parts in {path}")
    return [part for _, part in sorted(parts)]


def probability_batches(path, column, expected_pixels):
    parts = parquet_parts(path)
    count = 0
    for part in parts:
        with pq.ParquetFile(part) as reader:
            count += reader.metadata.num_rows
    if count != expected_pixels:
        raise ValueError(f"{path} contains {count} pixels, expected {expected_pixels}")
    for part in parts:
        with pq.ParquetFile(part) as reader:
            for batch in reader.iter_batches(batch_size=BATCH_PIXELS, columns=[column]):
                yield batch.column(0).to_numpy(zero_copy_only=False).astype(np.float32, copy=False)


def map_blocks(batches, height, width):
    remaining = np.empty(0, dtype=np.float32)
    rows_per_block = max(1, BLOCK_PIXELS // width)
    for y in range(0, height, rows_per_block):
        end = min(y + rows_per_block, height)
        pixels = (end - y) * width
        values = np.empty(pixels, dtype=np.float32)
        start = 0
        while start < pixels:
            if not len(remaining):
                remaining = next(batches)
            take = min(pixels - start, len(remaining))
            values[start:start + take] = remaining[:take]
            remaining = remaining[take:]
            start += take
        yield y, end, values


def summarize_map(labels, ids, path, column):
    height, width = labels.shape
    order = ids.argsort()
    sorted_ids = ids[order]
    count = np.zeros(len(ids), dtype=np.int64)
    sums = np.zeros(len(ids), dtype=np.float64)
    hist = np.zeros((len(ids), BINS), dtype=np.uint32)
    batches = probability_batches(path, column, height * width)

    # Each map is stored bottom-to-top. Read the corresponding label rows in that order.
    for start, end, values in map_blocks(batches, height, width):
        block = np.asarray(labels[height - end:height - start, :])[::-1].ravel()
        positive = block > 0
        candidates = block[positive]
        positions = np.searchsorted(sorted_ids, candidates)
        in_range = positions < len(sorted_ids)
        matched = in_range.copy()
        matched[in_range] = sorted_ids[positions[in_range]] == candidates[in_range]
        rows = order[positions[matched]]
        probabilities = np.clip(values[positive][matched], 0, 1)
        if not np.isfinite(probabilities).all():
            raise ValueError(f"{path} has non-finite probabilities inside cell labels")
        count += np.bincount(rows, minlength=len(ids))
        sums += np.bincount(rows, weights=probabilities, minlength=len(ids))
        scores = np.rint(probabilities * (BINS - 1)).astype(np.intp)
        np.add.at(hist, (rows, scores), 1)

    mean = np.divide(sums, count, out=np.full(len(ids), np.nan), where=count > 0)
    median = np.full(len(ids), np.nan)
    for start in range(0, len(ids), 4096):
        end = min(start + 4096, len(ids))
        cumulative = hist[start:end].cumsum(axis=1, dtype=np.uint32)
        target = (count[start:end] + 1) // 2
        values = (cumulative >= target[:, None]).argmax(axis=1) / (BINS - 1)
        median[start:end] = np.where(target > 0, values, np.nan)
    return mean.astype(np.float32), median.astype(np.float32), count


def collect_outputs(directories):
    metrics = {}
    maps = {}
    for directory in directories:
        for tmp in (directory, directory / "spoQC_tmp"):
            for path in sorted(tmp.glob(f"{MAP_PREFIX}*")):
                name = path.name.removeprefix(MAP_PREFIX).removesuffix(".parquet")
                if re.fullmatch(r"hqpr_\d+|hqtr|hqcr", name):
                    previous = maps.setdefault(name, path)
                    if previous.resolve() != path.resolve():
                        raise ValueError(f"Multiple smoothed maps for {name}: {previous} and {path}")
            for path in sorted((tmp / "metrices/hqcr").glob("*_output_hqcr.parquet")):
                metrics.setdefault(path.name, path)
            traffic = tmp / "traffic_light_output_hqcr.parquet"
            if traffic.is_file():
                metrics.setdefault(traffic.name, traffic)
    return maps, metrics


def attach_metrics(obs, metrics):
    for filename, path in sorted(metrics.items()):
        df = pd.read_parquet(path)
        df.index = df.index.astype(str)
        if not df.index.is_unique or not df.index.isin(obs.index).all():
            raise ValueError(f"{path}: metric rows do not match unique input table observation IDs")
        df = df.reindex(obs.index)
        source = filename.removesuffix("_output_hqcr.parquet").removesuffix(".parquet")
        for name, column in df.items():
            if column.dtype == object and column.map(lambda value: isinstance(value, (np.ndarray, list))).any():
                column = column.map(lambda value: json.dumps(np.asarray(value).tolist())
                                    if isinstance(value, (np.ndarray, list)) else value)
            output_name = f"spoqc_{source}_{name}" if name in obs.columns else name
            if output_name in obs.columns:
                raise ValueError(f"Metric column {output_name} already exists")
            obs[output_name] = column


def export(args):
    source = zarr.open_group(args.spatialdata, mode="r")
    table = source["tables/table"]
    labels = source["labels/cell_labels/0"]
    if labels.ndim != 2:
        raise ValueError(f"Expected a 2D cell-label image, found {labels.shape}")
    obs = ad.io.read_elem(table["obs"])
    if not obs.index.is_unique:
        raise ValueError("Input table observation IDs must be unique")
    id_column = "cell_labels" if "cell_labels" in obs.columns else "cell_id"
    numeric_ids = pd.to_numeric(obs[id_column], errors="raise").to_numpy()
    ids = numeric_ids.astype(np.int64)
    if not np.array_equal(numeric_ids, ids) or np.any(ids <= 0) or len(np.unique(ids)) != len(ids):
        raise ValueError(f"obs[{id_column!r}] must contain unique positive integer label IDs")

    maps, metrics = collect_outputs(args.spoqc_output)
    expected = {"hqtr", "hqcr", *(f"hqpr_{staining}" for staining in args.staining)}
    missing = expected - maps.keys()
    if missing:
        raise ValueError("Missing smoothed SpoQC maps: " + ", ".join(sorted(missing)))
    attach_metrics(obs, metrics)
    for name, path in sorted(maps.items()):
        mean, median, count = summarize_map(labels, ids, path, f"{name}_beliefs_smoothed")
        obs[f"spoqc_{name}_mean"] = mean
        obs[f"spoqc_{name}_median"] = median
        obs["spoqc_n_pixels"] = count

    # Keep the spatial coordinates and gene metadata. The expression matrix and layers
    # can be large, and are not needed to display per-cell QC scores in the browser.
    var = ad.io.read_elem(table["var"])
    obsm = ad.io.read_elem(table["obsm"]) if "obsm" in table else {}
    result = ad.AnnData(obs=obs, var=var, obsm=obsm)
    with ad.settings.override(zarr_write_format=2):
        result.write_zarr(args.output)
    if args.h5ad:
        result.write_h5ad(args.h5ad)
    print(f"Exported {result.n_obs} cells and {result.obs.shape[1]} observation columns to {args.output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spatialdata", type=Path, required=True)
    parser.add_argument("--spoqc-output", type=Path, nargs="+", required=True)
    parser.add_argument("--staining", action="append", default=[])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--h5ad", type=Path)
    export(parser.parse_args())


if __name__ == "__main__":
    main()
