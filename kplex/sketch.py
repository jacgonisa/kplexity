"""Sketched k schedules: count a few k values instead of the full 5..151 ladder.

The size-anchored schedule spreads n - 6 k values evenly over [kmin, kmax] and puts 6 more across the first rise,
centred at k = 1.213 log4(G) - 1.95 (+/- 5), where G is the number of ACGT bases. Most of the curve's shape is in
that rise, and its position is set by genome size, so the schedule is known before anything is counted.
Fits to sketched curves are weighted by k-spacing (kplex.model.spacing_weights), which `kplex fit` does
automatically for unevenly spaced k.

Benchmark (326 annotated DToL genomes, kplexity model): with 15 k values the median combined drift of the
fitted parameters from the full-ladder fit is 1.8% (12.9% for 15 evenly spaced k); T2T-CHM13, 0.9%.
Below about 12 k values the fit becomes unreliable.
"""
from __future__ import annotations
import math

KMIN, KMAX = 5, 151
DEFAULT_N = 15
MIN_N = 10
N_BAND, HALF_WIDTH = 6, 5.0


def rise_centre(genome_bp: float) -> float:
    """Expected position of the first rise from genome size (fitted across 326 eukaryotes and 100 prokaryotes)."""
    return 1.213 * math.log(genome_bp) / math.log(4) - 1.95


def _fill(points, n, kmin, kmax):
    """Round, dedupe and clip; if rounding merged points, add k at the middle of the widest gaps until there are n."""
    ks = sorted({int(round(x)) for x in points if kmin <= round(x) <= kmax})
    while len(ks) < n:
        gaps = [(b - a, a, b) for a, b in zip(ks, ks[1:]) if b - a > 1]
        if not gaps:
            break
        _, a, b = max(gaps)
        ks = sorted(ks + [(a + b) // 2])
    return ks[:n]


def _linspace(a, b, n):
    return [a] if n == 1 else [a + (b - a) * i / (n - 1) for i in range(n)]


def size_anchored(genome_bp: float, n: int = DEFAULT_N, kmin: int = KMIN, kmax: int = KMAX) -> list[int]:
    """n distinct integer k values in [kmin, kmax], including both ends."""
    if n < MIN_N:
        raise ValueError(f"a sketch needs at least {MIN_N} k values (got {n})")
    if n >= kmax - kmin + 1:
        return list(range(kmin, kmax + 1))
    c = rise_centre(genome_bp)
    pts = _linspace(kmin, kmax, n - N_BAND) + _linspace(c - HALF_WIDTH, c + HALF_WIDTH, N_BAND)
    return _fill(pts, n, kmin, kmax)
