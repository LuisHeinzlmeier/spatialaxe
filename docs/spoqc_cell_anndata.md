# SpoQC SpatialData output

When SpoQC is enabled, the pipeline writes `results/<sample>/qc/<sample>_qc.zarr` for each sample. The export works in masks and report modes. It also writes the standalone AnnData Zarr v2 store `<sample>_qc.anndata.zarr`. To also write `<sample>_qc.h5ad`, set `spoqc_cell_h5ad: true` in the parameters file. To disable the export, set `spoqc_cell_anndata: false`.

The AnnData keeps the input cell IDs, observation order, original `obs`, gene metadata (`var`) and coordinates (`obsm`). It leaves out the expression matrix and layers to keep the browser-facing file small. The input SpatialData is read only.

The export joins SpoQC's per-cell metric tables by observation ID. This includes GeneralQC, CellQC, DoubletQC, VoidQC and the HQCR traffic light when available. When a metric name already exists in the input, the new column gets a `spoqc_<step>_` prefix. List-valued metric columns are encoded as JSON strings. AmbientQC's gene-level table and the transcript-level doublet table remain in the regular SpoQC outputs. SpoQC intermediates held only in memory, such as the original HQCR per-cell prior and report-only cell-cycle scores, are unavailable to this export. The HQCR mean and median below summarize the published refined map.

- `spoqc_hqpr_<staining>_mean` and `spoqc_hqpr_<staining>_median` for every staining.
- `spoqc_hqtr_mean`, `spoqc_hqtr_median`, `spoqc_hqcr_mean` and `spoqc_hqcr_median`.
- `spoqc_n_pixels` counts the full-resolution pixels of each cell's label. A cell with zero pixels has missing means and medians. Medians use probabilities rounded to 0.001, with the lower middle value for an even pixel count.

The map rows are flipped vertically to line up with the labels. The exporter reads the maps and label image in blocks and keeps histogram storage proportional to the number of cells.

## SpatialData layout

- `images/spoqc_maps` is a five-level uint8 pyramid with channels ordered by HQPR staining index, then HQTR and HQCR. HQPR channels use the staining names from `morphology_focus`, lowercased with punctuation replaced by underscores. Pixel values are the smoothed beliefs clipped to [0, 1], multiplied by 254 and rounded. Value 255 marks a missing belief.
- `labels/spoqc_mask_hqpr_<staining>`, `labels/spoqc_mask_hqtr` and `labels/spoqc_mask_hqcr` preserve SpoQC's smoothed binary labels. A value of 1 means high quality and 0 means low quality.
- `labels/spoqc_mask` sums the binary masks. Its range is 0 to the number of HQPR stainings plus two.
- `shapes/cell_circles` is copied from the input. `tables/table` is the per-cell table linked to it by `cell_id`. The table records the combined-mask maximum and missing pixel counts in `uns`.

The rasters are unflipped into the input image orientation and registered to `global`. Raster pyramid levels use nearest-neighbor sampling so mask values and the missing-belief marker remain intact. A missing binary mask pixel is reported and stops the export because it cannot be represented by either binary value.

A common coarser grid uses a scale transform. The export rejects samples whose masks have different grids because their combined pixel sum would require a resampling rule.
