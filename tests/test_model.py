"""Tests for the curve models: `python -m pytest tests`."""
import json
from pathlib import Path

import numpy as np
import pytest

from kplex import model as M
from kplex.cli import main
from kplex import sketch as SK
from kplex.kmers import parse_k_range

DATA = Path(__file__).parent / "data"
K = np.arange(5, 152, dtype=float)


@pytest.mark.parametrize("truth", [
    dict(b=1.386, k0=15.0, A=0.95, S=0.10, lam=0.03),
    dict(b=1.2, k0=12.5, A=0.85, S=0.35, lam=0.06),
    dict(b=1.386, k0=16.5, A=0.99, S=0.02, lam=0.02),
])
def test_recovers_parameters(truth):
    y = M.kplexity_curve(K, **truth)
    f = M.fit_kplexity(K, y, genome_bp=4 ** truth["k0"])
    assert f["rmse"] < 1e-5
    for p in ("A", "S", "k0"):
        assert f[p] == pytest.approx(truth[p], abs=2e-3)
    assert f["divergence"] == pytest.approx(1 - np.exp(-truth["lam"]), rel=0.02)
    assert f["plateau_reached"]


def test_random_sequence_has_no_repeats():
    G = 1e8
    y = M.random_curve(K, G)
    f = M.fit_kplexity(K, y, genome_bp=G)
    assert f["A"] == pytest.approx(1.0, abs=1e-3)
    assert f["S"] < 2e-3
    assert f["k0"] == pytest.approx(np.log(G) / np.log(4), abs=0.05)
    assert f["b"] == pytest.approx(np.log(4), rel=0.02)


@pytest.mark.parametrize("name, A, S, k0", [
    ("chm13", 0.9262, 0.2301, None),
    ("tair12", 0.8651, 0.1095, None),
])
def test_reference_curves(name, A, S, k0):
    csv = DATA / f"{name}.kplex.csv"
    k, y, _ = M.load_curve(str(csv))
    G = M.genome_bp_from_csv(str(csv))
    f = M.fit_kplexity(k, y, G)
    assert f["A"] == pytest.approx(A, abs=1e-3)
    assert f["S"] == pytest.approx(S, abs=1e-3)
    assert f["k0"] == pytest.approx(np.log(G) / np.log(4), abs=0.5)
    assert f["r_squared"] > 0.999


def test_sketched_curve_close_to_full():
    csv = DATA / "chm13.kplex.csv"
    k, y, _ = M.load_curve(str(csv))
    G = M.genome_bp_from_csv(str(csv))
    full = M.fit_kplexity(k, y, G)
    sel = np.isin(k, [5, 10, 13, 15, 16, 17, 18, 20, 23, 27, 33, 41, 55, 80, 115, 151])
    sk = M.fit_kplexity(k[sel], y[sel], G)
    assert sk["weighted"]
    assert sk["A"] == pytest.approx(full["A"], abs=0.01)
    assert sk["S"] == pytest.approx(full["S"], abs=0.02)


def test_size_anchored_schedule():
    ks = SK.size_anchored(3.06e9, 15)                     # T2T-CHM13
    assert ks == [5, 12, 14, 16, 18, 20, 22, 23, 42, 60, 78, 96, 114, 133, 151]
    for G in (1e6, 1e8, 1e10):
        for n in (10, 15, 25, 40):
            ks = SK.size_anchored(G, n)
            assert len(ks) == n and ks[0] == 5 and ks[-1] == 151 and ks == sorted(set(ks))
    with pytest.raises(ValueError):
        SK.size_anchored(1e8, 5)


def test_size_anchored_sketch_recovers_full_fit():
    csv = DATA / "chm13.kplex.csv"
    k, y, _ = M.load_curve(str(csv))
    G = M.genome_bp_from_csv(str(csv))
    full = M.fit_kplexity(k, y, G)
    sel = np.isin(k, SK.size_anchored(G, 15))
    sk = M.fit_kplexity(k[sel], y[sel], G)
    for p in ("A", "S", "k0", "b"):
        assert sk[p] == pytest.approx(full[p], rel=0.03)


def test_parse_k_list():
    assert parse_k_range("5:11:3") == [5, 8, 11]
    assert parse_k_range("5,12,14,151") == [5, 12, 14, 151]


def test_cli_schedule(capsys):
    assert main(["schedule", "--genome-bp", "3.06e9"]) == 0
    assert capsys.readouterr().out.strip() == "5,12,14,16,18,20,22,23,42,60,78,96,114,133,151"


@pytest.mark.parametrize("model, key", [("kplexity", "S"), ("gamma", "a")])
def test_cli_fit_writes_json(tmp_path, model, key):
    out = tmp_path / "t"
    assert main(["fit", str(DATA / "tair12.kplex.csv"), "--model", model, "-o", str(out), "--plot"]) == 0
    fit = json.loads(Path(f"{out}.fit.json").read_text())
    assert key in fit and "error" not in fit
    assert Path(f"{out}.png").exists()
