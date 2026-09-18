# perturbation-status-methods

Omnibenchmark wrappers around the frozen, config-driven perturbation-status
execution system. `shared.py` creates the reusable CalcPerturbSig signature and
target shards. `method.py` runs Mixscape, Mixscale, PS, or the target-expression
baseline as equal benchmark branches. Shared technical dependencies do not imply
method ranking or workflow preference.
