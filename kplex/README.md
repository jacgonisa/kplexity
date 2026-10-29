# kplex — k-plexity CLI

Fraction of unique k-mers vs *k*: a fast genome repeat-landscape QC metric.
Three commands: **`run`** (FASTA → curve), **`fit`** (curve → model parameters), **`plot`** (curve → figure).

## Install

```bash
# from the repo root (contains pyproject.toml)
pip install -e .            # or: pipx install .
```

Pure-Python deps (numpy/pandas/scipy/matplotlib) install automatically. The k-mer
counter is an external binary — point `kplex run` at it once via `--kplex` or an env var:

```bash
export KPLEX_BIN=/path/to/FASTK/Kplex     # FASTK/Kplex (preferred), OR
conda install -c bioconda kmc             # KMC fallback (kplex run --tool kmc)
```

## Usage

```bash
# 1) compute the curve (auto-detects FASTK/Kplex or KMC)
kplex run genome.fa -o mygenome --chromosomes-only -T 20
#   -> mygenome.kplex.csv , mygenome.theoretical.csv
#   add --fit --plot to also fit + render in one go

# 2) fit the curve (the kplexity model; --model gamma adds one parameter)
kplex fit mygenome.theoretical.csv -o mygenome
#   -> mygenome.fit.json  (A, identical = 1-A, S, lam, divergence, k0, b, R², RMSE, AICc)

# 3) plot (optionally overlay the fit: chance collisions x A, near-copy gap, plateau)
kplex plot mygenome.theoretical.csv --fit mygenome.fit.json -o mygenome.png
```

## Curve schema
`run` writes `k, unique_kmers, total_kmers, fraction_unique` (and a `.theoretical.csv`
adding `fraction_unique_theoretical`, the length-based denominator, plus a truncation metric).
`fit`/`plot` prefer `fraction_unique_theoretical` when present.

## Parameters
`y(k) = B(k)·[A − S·e^(−λ(k−k0))]`, `B = (1 − e^(−L))/L`, `L = 2e^(−b(k−k0))`
- **A** — resolvable fraction (plateau); **identical** = 1 − A, copies no k separates
- **S** — near-copy share at k0 (k-mers in diverged repeat copies)
- **lam**, **divergence** = 1 − e^(−lam) — how fast and how diverged those copies are
- **k0** ≈ log4(genome size), **b** — position and steepness of the first rise (random sequence: b = ln 4)
- `--model gamma` adds **a** (spread of divergence across copies)
