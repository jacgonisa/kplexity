# kplex — k-plexity CLI

Fraction of unique k-mers vs *k*: a fast genome repeat-landscape QC metric.
Three commands: **`run`** (FASTA → curve), **`fit`** (curve → double-sigmoid params), **`plot`** (curve → figure).

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

# 2) fit the double sigmoid
kplex fit mygenome.theoretical.csv -o mygenome
#   -> mygenome.fit.json  (R², asymptote=L1+L2, L2_fraction=L2/(L1+L2), k0_1, k0_2)

# 3) plot (optionally overlay the fit + component decomposition)
kplex plot mygenome.theoretical.csv --fit mygenome.fit.json -o mygenome.png
```

## Curve schema
`run` writes `k, unique_kmers, total_kmers, fraction_unique` (and a `.theoretical.csv`
adding `fraction_unique_theoretical`, the length-based denominator, plus a truncation metric).
`fit`/`plot` prefer `fraction_unique_theoretical` when present.

## Parameters (shape of the curve)
- **asymptote** = L1 + L2 — the high-k plateau (overall uniqueness)
- **L2_fraction** = L2/(L1+L2) — relative height of the second component
- **k0_1, k0_2** — the two midpoints; **s1, s2** — the two steepnesses
