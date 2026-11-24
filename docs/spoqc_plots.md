# SpoQC plots and HTML reports

Spatialaxe disables SpoQC plots, figures, and HTML reports by default in this fork. SpoQC still writes data outputs, including annotation tables, metrics, and masks. Use `--spoqc_plots true` to generate figures and reports again. With `--spoqc_until masks`, the pipeline still stops before the downstream analysis and report steps.

The default requires a SpoQC version that accepts `--no-plots` in the container or on `--spoqc_pythonpath`. Use SpoQC revision `e811fe94ab992265c47be5b0ce47e22da2c27d1c` or a compatible newer revision. With an older SpoQC build, enable plots to omit the new CLI flag.
