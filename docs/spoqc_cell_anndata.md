# SpoQC output in SpatialData

When SpoQC is enabled, spatialaxe publishes the complete annotated SpatialData at `results/<sample>/qc/<sample>_qc.zarr`. This works in masks and report modes. Set `spoqc_cell_anndata: false` in the parameters file to disable this export.

The pipeline copies the input before adding QC. This lets Nextflow retry or resume without changing an upstream task's inputs or a reused SpatialData store. The copy retains the original images, labels, points, shapes and sole `table`, including its expression matrix, layers, gene metadata, spatial coordinates and table annotation. The output requires space for the original bundle plus the new rasters. Publishing with copy mode also retains a complete copy in the Nextflow work directory.

For an existing store that is safe to update, the same exporter accepts `--in-place --obs-backup <path-outside-store>`. It copies the original `tables/table/obs` group before making changes, adds each raster with `write_element` without overwriting an existing element, and replaces only `table.obs`. Run it when the store is free of competing writers and pipeline readers.

## QC elements

- `images/spoqc_maps` is a five-level uint8 pyramid with channels ordered by HQPR staining index, then HQTR and HQCR. HQPR channel names come from `morphology_focus`, in lowercase with punctuation replaced by underscores. Smoothed beliefs are clipped to [0, 1], multiplied by 254 and rounded. Value 255 marks a missing belief.
- `labels/spoqc_mask_hqpr_<staining>`, `labels/spoqc_mask_hqtr` and `labels/spoqc_mask_hqcr` preserve the smoothed binary labels. A value of 1 means high quality and 0 means low quality.
- `labels/spoqc_mask` sums the binary masks. Its maximum is the number of HQPR stainings plus two.
- `tables/table/obs` keeps the original IDs, order and columns. It adds `spoqc_hqpr_<staining>_mean` and `_median` for each staining, the matching `spoqc_hqtr_*` and `spoqc_hqcr_*` summaries, and `spoqc_n_pixels`. Cells with zero pixels have missing summaries. Medians use the lower middle of the probabilities rounded to 0.001.

The rasters are unflipped and registered to `global`. The export streams the source pixels in blocks and writes spatially chunked pyramids. Counts of missing beliefs and mask pixels are logged for each channel. Missing binary labels stop the export. Channels on a common coarser grid get a scale transform. Mixed native grids require a separate resampling rule.

The exporter joins saved GeneralQC, CellQC, DoubletQC, and VoidQC tables and the HQCR traffic light by observation ID. A metric that shares a name with an original column gets a `spoqc_<step>_` prefix. List-valued metrics are encoded as JSON strings. AmbientQC's saved table is per gene. HQCR's original per-cell prior and report-only intermediates are held in memory by SpoQC. The HQCR mean and median summarize the published refined map.
