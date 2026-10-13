# perturbation-status-mixscape-mixscale

Omnibenchmark wrappers around the frozen, config-driven perturbation-status
execution system. `shared.py` creates the reusable CalcPerturbSig signature and
target shards. `method.py` runs the organically linked Mixscape and Mixscale
branches. Shared technical dependencies do not imply method ranking or workflow
preference.

The reviewed implementation scripts and dataset-specific configuration are
versioned under `implementation/`; no runtime access to the original benchmark
development directory is required. Public dependencies are supplied by the
plan's Conda YAML files. Exact Seurat, presto, and Mixscale revisions are
archived under `vendor/` and installed into a task-local R library at runtime,
without a personal R library or network checkout.
