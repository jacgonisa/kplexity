# Bioconda recipe

`recipe/` is the bioconda recipe for the `kplexity` package: the `kplex` Python CLI plus the `Kplex` binary, built
from `Kplex.c` of the [FASTK fork](https://github.com/jacgonisa/FASTK) pinned to one commit. Only Kplex is
compiled; `FastK`, `Histex` and `Fastrm`, which Kplex calls for every k, come from bioconda's `fastk` package (the fork
is FastK 1.2 plus `Kplex.c`, so the two match). Kplex and kplex find them on `PATH`, so a conda install needs no
`KPLEX_BIN`.

## Submitting a release to bioconda

1. Tag and release this repo: `git tag -a v0.4.0 -m "kplexity 0.4.0" && git push origin v0.4.0`, then
   `gh release create v0.4.0 --title "kplexity 0.4.0" --notes "..."`.
2. Put the release tarball's checksum into `recipe/meta.yaml`:
   `curl -sL https://github.com/jacgonisa/kplexity/archive/refs/tags/v0.4.0.tar.gz | sha256sum`.
   If the FASTK fork has moved on and Kplex changed, update `fastk_commit` and its `sha256` the same way
   (`https://github.com/jacgonisa/FASTK/archive/<commit>.tar.gz`).
3. Fork [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes), copy `recipe/` to
   `recipes/kplexity/` on a new branch, and open a pull request. The bioconda CI builds and tests it on linux-64,
   linux-aarch64 and macOS; once it passes, comment `@BiocondaBot please add label` to request review.
4. Later releases: bioconda's autobump bot opens the version-update pull request from new GitHub tags.

## Testing the recipe locally

Build against local tarballs (no release needed): replace the two `url`/`sha256` pairs with `file://` paths to
`git archive`-style tarballs of this repo and of the fork, then

```bash
conda create -p ./cb -c conda-forge conda-build
./cb/bin/conda-build recipe -c conda-forge -c bioconda --override-channels --python 3.12
```

The recipe tests run `kplex --version`, `kplex --help`, `Kplex`, `kplex schedule` and `tests/test_run.py`, an end-to-end
`kplex run` on a toy genome.
