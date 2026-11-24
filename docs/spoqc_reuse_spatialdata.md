# Reuse SpatialData for SpoQC

For a SpoQC-only rerun, pass `--reuse_spatialdata bundles.csv` alongside the normal samplesheet and SpoQC settings. This skips SpatialData conversion and passes the listed bundles to SpoQC. Use one row per sample, matching the samplesheet ID, and paths visible to the executor.

```csv
sample,spatialdata
example,/shared/example/raw_bundle
```
