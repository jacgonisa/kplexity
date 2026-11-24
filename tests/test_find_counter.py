"""Kplex discovery: the conda location wins, and a script named Kplex on PATH (what a case-insensitive file system
returns for the kplex entry point) is never taken for the binary."""
import kplex.kmers as km


def _fake(path, head):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(head + b"\n")
    path.chmod(0o755)
    return path


def test_script_on_path_is_skipped(tmp_path, monkeypatch):
    _fake(tmp_path / "bin" / "Kplex", b"#!/usr/bin/env python")
    monkeypatch.setenv("PATH", str(tmp_path / "bin"))
    monkeypatch.delenv("KPLEX_BIN", raising=False)
    monkeypatch.setattr(km, "CONDA_KPLEX", tmp_path / "libexec" / "kplexity" / "Kplex")
    try:
        km.find_counter("fastk", None, None)
        raise AssertionError("a script was accepted as Kplex")
    except SystemExit:
        pass


def test_conda_location_is_found(tmp_path, monkeypatch):
    kp = _fake(tmp_path / "libexec" / "kplexity" / "Kplex", b"\x7fELF")
    _fake(tmp_path / "bin" / "Kplex", b"#!/usr/bin/env python")
    monkeypatch.setenv("PATH", str(tmp_path / "bin"))
    monkeypatch.delenv("KPLEX_BIN", raising=False)
    monkeypatch.setattr(km, "CONDA_KPLEX", kp)
    assert km.find_counter("fastk", None, None) == ("fastk", str(kp.resolve()))
