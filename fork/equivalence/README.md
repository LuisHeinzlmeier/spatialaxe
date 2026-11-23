# Check a rebuild

Run `fork/equivalence/check_tree.sh run/2026-11-19-spatialaxe integration/spoqc-fast` to check the pipeline files. If it reports differences, inspect them with `git diff` and record why they are expected.

Run the reference and candidate on the same inputs with the same parameters, image and SpoQC revision. Give each run its own work and output directories. For each sample, compare the two published `spoQC_tmp` folders with `python fork/equivalence/compare_outputs.py REF CANDIDATE --exact`. The script checks Parquet, CSV and TSV outputs. Row order must match.
