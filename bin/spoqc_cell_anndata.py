#!/usr/bin/env python3
"""Add SpoQC rasters and observation columns to a SpatialData store."""

import argparse
import json
import re
import shutil
import tempfile
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


def map_blocks(batches, height, width, block_pixels=BLOCK_PIXELS):
    remaining = np.empty(0, dtype=np.float32)
    rows_per_block = max(1, block_pixels // width)
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


def map_grid(map_path, height, width):
    pixels = 0
    for part in parquet_parts(map_path):
        with pq.ParquetFile(part) as reader:
            pixels += reader.metadata.num_rows
    for level in range(5):
        scale = 2 ** level
        shape = (height // scale, width // scale)
        if shape[0] * shape[1] == pixels:
            return shape, scale
    raise ValueError(f"{map_path} pixel count does not match a morphology grid")


def summarize_map(labels, ids, path, column):
    height, width = labels.shape
    (rows, cols), scale = map_grid(path, height, width)
    order = ids.argsort()
    sorted_ids = ids[order]
    count = np.zeros(len(ids), dtype=np.int64)
    observed = np.zeros(len(ids), dtype=np.int64)
    sums = np.zeros(len(ids), dtype=np.float64)
    hist = np.zeros((len(ids), BINS), dtype=np.uint32)
    batches = probability_batches(path, column, rows * cols)

    def match_labels(block):
        positive = block > 0
        candidates = block[positive]
        positions = np.searchsorted(sorted_ids, candidates)
        in_range = positions < len(sorted_ids)
        matched = in_range.copy()
        matched[in_range] = sorted_ids[positions[in_range]] == candidates[in_range]
        return positive, matched, order[positions[matched]]

    # Each map is stored bottom-to-top. Read the corresponding label rows in that order.
    for start, end, values in map_blocks(batches, rows, cols, BLOCK_PIXELS // scale**2):
        block = np.asarray(labels[(rows - end) * scale:(rows - start) * scale, :])[::-1].ravel()
        positive, matched, indices = match_labels(block)
        count += np.bincount(indices, minlength=len(ids))
        if scale != 1:
            expanded = values.reshape(end - start, cols).repeat(scale, axis=0).repeat(scale, axis=1)
            values = np.pad(expanded, ((0, 0), (0, width - cols * scale)), constant_values=np.nan).ravel()
        probabilities = values[positive][matched]
        valid = np.isfinite(probabilities)
        indices, probabilities = indices[valid], np.clip(probabilities[valid], 0, 1)
        observed += np.bincount(indices, minlength=len(ids))
        sums += np.bincount(indices, weights=probabilities, minlength=len(ids))
        scores = np.rint(probabilities * (BINS - 1)).astype(np.intp)
        np.add.at(hist, (indices, scores), 1)

    if rows * scale < height:
        _, _, indices = match_labels(np.asarray(labels[rows * scale:, :]).ravel())
        count += np.bincount(indices, minlength=len(ids))
    mean = np.divide(sums, observed, out=np.full(len(ids), np.nan), where=observed > 0)
    median = np.full(len(ids), np.nan)
    for start in range(0, len(ids), 4096):
        end = min(start + 4096, len(ids))
        cumulative = hist[start:end].cumsum(axis=1, dtype=np.uint32)
        target = (observed[start:end] + 1) // 2
        values = (cumulative >= target[:, None]).argmax(axis=1) / (BINS - 1)
        median[start:end] = np.where(target > 0, values, np.nan)
    return mean.astype(np.float32), median.astype(np.float32), count


def write_qc_spatialdata(path, maps, obs, input_obs, obs_backup):
    import dask
    import dask.array as da
    import spatialdata as sd
    from multiscale_spatial_image.to_multiscale import Methods
    from spatialdata.models import Image2DModel, Labels2DModel
    from spatialdata.transformations import Identity, Scale

    source = sd.read_zarr(str(path), selection=("images", "labels"))
    morphology = source.images["morphology_focus"]
    full = morphology["scale0"].image if hasattr(morphology, "keys") else morphology
    height, width = full.sizes["y"], full.sizes["x"]
    stainings = [str(name) for name in full.coords["c"].values]
    ordered = sorted(maps, key=lambda name: (0, int(name[5:])) if name.startswith("hqpr_")
                     else (1 if name == "hqtr" else 2, 0))
    names = []
    for key in ordered:
        if key.startswith("hqpr_"):
            index = int(key[5:])
            if index >= len(stainings):
                raise ValueError(f"Staining {index} is missing from morphology_focus")
            label = re.sub(r"[^a-z0-9]+", "_", stainings[index].lower()).strip("_")
            if not label:
                raise ValueError(f"Staining {index} has no usable name")
            names.append(f"hqpr_{label}")
        else:
            names.append(key)
    if len(names) != len(set(names)):
        raise ValueError("Staining names must be unique after normalizing")

    additions = {"spoqc_maps", "spoqc_mask", *(f"spoqc_mask_{name}" for name in names)}
    existing = {element_path.split("/")[-1] for element_path in source.elements_paths_on_disk()}
    conflicts = additions & existing
    if conflicts:
        raise ValueError(f"{path} already contains QC elements: {', '.join(sorted(conflicts))}")

    grids = {key: map_grid(maps[key], height, width) for key in ordered}
    if len(set(grids.values())) != 1:
        raise ValueError("SpoQC channels occupy different grids and cannot be stacked without resampling")
    (rows, cols), scale = grids[ordered[0]]
    transform = Identity() if scale == 1 else Scale([scale, scale], axes=("y", "x"))
    with tempfile.TemporaryDirectory(prefix="spoqc-rasters-", dir=".") as scratch:
        store = zarr.open_group(scratch, mode="w")
        channels = []
        labels = {}
        mask_arrays = []
        for key, name in zip(ordered, names):
            stripe = max(1, BLOCK_PIXELS // cols)
            probability = store.create_array(f"p_{name}", shape=(rows, cols), chunks=(stripe, 4096), dtype="u1")
            mask = store.create_array(f"m_{name}", shape=(rows, cols), chunks=(stripe, 4096), dtype="u1")
            missing_probability = missing_mask = 0
            for suffix, raster in (("beliefs", probability), ("mask", mask)):
                batches = probability_batches(maps[key], f"{key}_{suffix}_smoothed", rows * cols)
                for start, end, flat in map_blocks(batches, rows, cols):
                    valid = np.isfinite(flat)
                    if suffix == "beliefs":
                        missing_probability += int(np.count_nonzero(~valid))
                        encoded = np.full(flat.shape, 255, dtype=np.uint8)
                        encoded[valid] = np.rint(np.clip(flat[valid], 0, 1) * 254).astype(np.uint8)
                    else:
                        missing_mask += int(np.count_nonzero(~valid))
                        if np.any(valid & (flat != 0) & (flat != 1)):
                            raise ValueError(f"{key} contains a non-binary mask value")
                        if not valid.all():
                            continue
                        encoded = flat.astype(np.uint8)
                    raster[rows - end:rows - start, :] = encoded.reshape(end - start, cols)[::-1]
            print(f"{name}: {missing_probability} missing beliefs, {missing_mask} missing mask pixels "
                  f"out of {rows * cols}", flush=True)
            if missing_mask:
                raise ValueError(f"{key} has missing mask values that cannot be represented as binary uint8")
            channels.append(da.from_zarr(probability))
            mask_arrays.append(da.from_zarr(mask))
            labels[f"spoqc_mask_{name}"] = Labels2DModel.parse(
                mask_arrays[-1], dims=("y", "x"), transformations={"global": transform},
                scale_factors=[2, 2, 2, 2], chunks=(4096, 4096), method=Methods.DASK_IMAGE_NEAREST,
            )

        if len(mask_arrays) > 255:
            raise ValueError("Combined mask has too many channels for uint8")
        labels["spoqc_mask"] = Labels2DModel.parse(
            da.stack(mask_arrays).sum(axis=0, dtype=np.uint8), dims=("y", "x"),
            transformations={"global": transform}, scale_factors=[2, 2, 2, 2],
            chunks=(4096, 4096), method=Methods.DASK_IMAGE_NEAREST,
        )
        pixel_data = da.stack(channels)
        image = Image2DModel.parse(
            pixel_data, dims=("c", "y", "x"), c_coords=names,
            transformations={"global": transform}, scale_factors=[2, 2, 2, 2],
            chunks=(1, 4096, 4096), method=Methods.XARRAY_COARSEN,
        )
        # Keep the image pyramid metadata, but select actual encoded pixels.
        # Averaging the encoded bytes would mix missing value 255 into valid pixels.
        for index, level in enumerate(image.children.values()):
            factor = 2 ** index
            plane = level["image"]
            y, x = plane.sizes["y"], plane.sizes["x"]
            offset = factor // 2
            data = pixel_data[:, offset:y * factor:factor, offset:x * factor:factor]
            data = data.rechunk((1, 4096, 4096))
            level.dataset = level.to_dataset().assign(image=plane.copy(data=data))
        if obs_backup is not None:
            parent = obs_backup.parent
            parent.mkdir(parents=True, exist_ok=True)
            if obs_backup.exists():
                raise FileExistsError(f"Observation backup already exists: {obs_backup}")
            # Complete the backup before the first write to the existing store.
            with tempfile.TemporaryDirectory(prefix="obs-backup-", dir=parent) as staging:
                copied = Path(staging) / "obs"
                shutil.copytree(path / "tables/table/obs", copied)
                copied.rename(obs_backup)
            print(f"Backed up observations to {obs_backup}", flush=True)

        source.images["spoqc_maps"] = image
        for name, label in labels.items():
            source.labels[name] = label
        with dask.config.set(scheduler="threads", num_workers=2):
            source.write_element(["spoqc_maps", *labels], overwrite=False)

    group = zarr.open_group(path, mode="r+", use_consolidated=False)
    table = group["tables/table"]
    original = ad.io.read_elem(table["obs"])
    pd.testing.assert_frame_equal(original, input_obs)
    ad.io.write_elem(table, "obs", obs)
    if source.has_consolidated_metadata():
        source.write_consolidated_metadata()
    saved = ad.io.read_elem(zarr.open_group(path, mode="r")["tables/table/obs"])
    pd.testing.assert_frame_equal(saved, obs)
    print(f"Updated {path}: {len(obs)} cells, {len(obs.columns)} observation columns, "
          f"combined-mask maximum {len(mask_arrays)}", flush=True)


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
    input_obs = obs.copy(deep=True)
    attach_metrics(obs, metrics)
    for name, path in sorted(maps.items()):
        mean, median, count = summarize_map(labels, ids, path, f"{name}_beliefs_smoothed")
        obs[f"spoqc_{name}_mean"] = mean
        obs[f"spoqc_{name}_median"] = median
        obs["spoqc_n_pixels"] = count

    destination = args.spatialdata.resolve() if args.in_place else args.output.resolve()
    if not args.in_place:
        # Copy the complete bundle, including the expression matrix and any source
        # elements this exporter does not need to read.
        shutil.copytree(args.spatialdata.resolve(), destination)
    write_qc_spatialdata(destination, maps, obs, input_obs, args.obs_backup)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spatialdata", type=Path, required=True)
    parser.add_argument("--spoqc-output", type=Path, nargs="+", required=True)
    parser.add_argument("--staining", action="append", default=[])
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--output", type=Path, help="Write a full annotated copy")
    target.add_argument("--in-place", action="store_true", help="Update the input store")
    parser.add_argument("--obs-backup", type=Path, help="Backup directory for the original observations")
    args = parser.parse_args()
    if args.in_place != (args.obs_backup is not None):
        parser.error("--obs-backup is required for --in-place and only applies to that mode")
    if args.in_place and args.obs_backup.resolve().is_relative_to(args.spatialdata.resolve()):
        parser.error("--obs-backup must be outside the input SpatialData store")
    export(args)



if __name__ == "__main__":
    main()
