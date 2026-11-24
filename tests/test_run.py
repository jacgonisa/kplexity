"""End-to-end `kplex run` on a toy genome; skipped when Kplex/FastK are not installed."""
import random

import pandas as pd
import pytest

from kplex.cli import main
from kplex.kmers import find_counter, missing_fastk_tools


def _kplex_or_skip():
    try:
        tool, path = find_counter("fastk", None, None)
    except SystemExit:
        pytest.skip("Kplex not found (conda install -c bioconda kplexity, or set KPLEX_BIN)")
    if missing_fastk_tools(path):
        pytest.skip(f"FastK programs missing: {missing_fastk_tools(path)}")
    return path


def toy_fasta(path, seed=1):
    """Two chromosomes of random sequence; chr1 carries 20 identical copies of a 2-kb repeat and an N gap."""
    rng = random.Random(seed)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    rep = rand(2000)
    chr1 = rand(60_000) + rep * 20 + "N" * 500 + rand(60_000)
    chr2 = rand(100_000)
    with open(path, "w") as f:
        for name, s in (("chr1", chr1), ("chr2", chr2)):
            f.write(f">{name}\n")
            f.writelines(s[i:i + 80] + "\n" for i in range(0, len(s), 80))
    return [60_000 + 40_000, 60_000, 100_000]       # ACGT run lengths


def test_run_toy_genome(tmp_path):
    _kplex_or_skip()
    runs = toy_fasta(tmp_path / "toy.fa")
    prefix = tmp_path / "toy"
    assert main(["run", str(tmp_path / "toy.fa"), "-o", str(prefix), "-k", "7,13,17,21,31,51", "-T", "2"]) == 0
    c = pd.read_csv(f"{prefix}.kplex.csv")
    assert list(c.k) == [7, 13, 17, 21, 31, 51]
    # Kplex's total_kmers is the exact, N-free count of k-mer positions
    for k, t in zip(c.k, c.total_kmers):
        assert t == sum(L - k + 1 for L in runs if L >= k)
    f = c.set_index("k").fraction_unique
    assert f[7] < 0.2 and f[13] > 0.7                           # chance collisions resolve near log4(G) ~ 9
    # past the collisions, the 20 identical repeat copies keep all but one copy's k-mers non-distinct
    t51 = c.set_index("k").total_kmers[51]
    assert f[51] == pytest.approx(1 - (40_000 - 51 + 1 - 2000) / t51, abs=0.005)
