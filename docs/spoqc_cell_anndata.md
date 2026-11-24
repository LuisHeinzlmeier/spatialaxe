# Per-cell SpoQC output

When SpoQC is enabled, the pipeline writes `results/<sample>/qc/<sample>_qc.zarr` for each sample. The export works in masks and report modes and uses AnnData's Zarr v2 format. To also write `<sample>_qc.h5ad`, set `spoqc_cell_h5ad: true` in the parameters file. To disable the export, set `spoqc_cell_anndata: false`.

The AnnData keeps the input cell IDs, observation order, original `obs`, gene metadata (`var`) and coordinates (`obsm`). It leaves out the expression matrix and layers to keep the browser-facing file small. The input SpatialData is read only.

The export joins SpoQC's per-cell metric tables by observation ID. This includes GeneralQC, CellQC, DoubletQC, VoidQC and the HQCR traffic light when available. When a metric name already exists in the input, the new column gets a `spoqc_<step>_` prefix. List-valued metric columns are encoded as JSON strings. AmbientQC's gene-level table and the transcript-level doublet table remain in the regular SpoQC outputs. SpoQC intermediates held only in memory, such as the original HQCR per-cell prior and report-only cell-cycle scores, are unavailable to this export. The HQCR mean and median below summarize the published refined map.

- `spoqc_hqpr_<staining>_mean` and `spoqc_hqpr_<staining>_median` for every staining.
- `spoqc_hqtr_mean`, `spoqc_hqtr_median`, `spoqc_hqcr_mean` and `spoqc_hqcr_median`.
- `spoqc_n_pixels` counts the full-resolution pixels of each cell's label. A cell with zero pixels has missing means and medians. Medians use probabilities rounded to 0.001, with the lower middle value for an even pixel count.

The map rows are flipped vertically to line up with the labels. The exporter reads the maps and label image in blocks and keeps histogram storage proportional to the number of cells.
