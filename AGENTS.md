# Working in the spatialaxe fork

Use nf-core/spatialaxe as `upstream`. Keep `fork/base` on upstream `dev` with only this file and the tools under `fork/`.

- Make one branch per change that could be reviewed together. Start from `fork/base` or the change it depends on. Record dependencies, fork-only changes and check results in `fork/branches.txt`.
- Rebuild `integration/spoqc-fast` with `fork/rebuild_integration.sh`. Never commit directly to integration. Tag code used for runs and compare rebuilt pipeline files with the run tag.
- To prepare an upstream branch, run `fork/make_pr_branch.sh <change>`. If a dependency has been merged, update the change against upstream dev first. Ask before opening any PR.
- Push from the laptop. Keep GitHub credentials off the cluster. Use a separate checkout and run directory for benchmarks. Keep EMBL routing in the run's EMBL config.
- For speed changes, record the reference and candidate commits, inputs, queue wait, wall time, allocated and used CPU-hours, and output comparison in the commit message. Compare with the same SpoQC revision on both sides.
