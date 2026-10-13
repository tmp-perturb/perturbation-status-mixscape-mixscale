# Bundled R package sources

These source archives are installed offline into a task-local R library. They
replace image-build-time `install_github()` calls and never write to a user R
library.

| Archive | Upstream | Commit |
|---|---|---|
| `seurat-084d9e40.tar.gz` | `https://github.com/satijalab/seurat` | `084d9e40f0e669fe09a37cd0d20efb16c87d8e0c` |
| `presto-3c97180d.tar.gz` | `https://github.com/immunogenomics/presto` | `3c97180d90d330524b5b8353a342ab2f4b990d21` |
| `mixscale-3d3a2cd1.tar.gz` | `https://github.com/satijalab/Mixscale` | `3d3a2cd12e62382df01c056536d967534b2de3a9` |

The archives were produced with `git archive` from the exact commits above;
their upstream license files are included inside each archive.

SHA-256:

```text
3f783be5a76072ea7b6a7a4ab26436120d51cd88fd1ad4ac310d6ae42973f4a5  seurat-084d9e40.tar.gz
0293826e79eb905c16bfeeb5142b90f17b60b8c85f1441fa8be66700e6c75193  presto-3c97180d.tar.gz
3914d393b0f623be359b1b23a10bb36f61a979f4f60c83cd68dedcd9a6a15cd9  mixscale-3d3a2cd1.tar.gz
```
