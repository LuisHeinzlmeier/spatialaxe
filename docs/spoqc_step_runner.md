# SpoQC step runner

Run with `--spoqc --spoqc_container IMAGE` and choose a container that provides the SpoQC 0.1 step interface. The runner supports container execution with Docker, Singularity or Podman. Independent steps can run at the same time when the executor has resources available.

Use `--spoqc_until masks` to stop after the HQCR and HQTR probability maps and the HQPR maps for any supplied staining channels. The default, `report`, also runs the downstream analysis and report steps.
