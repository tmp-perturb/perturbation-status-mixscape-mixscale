# perturbation-status-mixscape-mixscale

Omnibenchmark wrappers around the frozen, config-driven perturbation-status
execution system. `shared.py` creates the reusable CalcPerturbSig signature and
target shards. `method.py` runs the organically linked Mixscape and Mixscale
branches. Shared technical dependencies do not imply method ranking or workflow
preference.

The reviewed implementation scripts and dataset-specific configuration are
versioned under `implementation/`; no runtime access to the original benchmark
development directory is required. Reproducible Apptainer definitions for the
two method environments are under `environment/` and pin upstream Git commits.
