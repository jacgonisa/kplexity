"""k-mer counting engine for `kplex run`.

Counts distinct vs total k-mers across a range of k using FASTK/Kplex (preferred)
or KMC, then computes the k-plexity curve (fraction_unique = distinct / total) and
the theoretical fraction (denominator = number of k-mer positions in the sequence).

Curve CSV schema (matches the original pipeline):
    k, unique_kmers, total_kmers, fraction_unique
Theoretical CSV adds:
    total_kmers_theoretical, fraction_unique_theoretical, truncation_fraction, ...
"""
from __future__ import annotations
import gzip
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

DEFAULT_CHROM_REGEX = r"(?i)chromosome|\bchr\b|\bLG\d|linkage group"


# ----------------------------------------------------------------------------- helpers
def parse_k_range(k_range: str) -> list[int]:
    a, b, c = k_range.split(":")
    return list(range(int(a), int(b) + 1, int(c)))


def find_counter(tool: str, kplex_bin: str | None, kmc_bin: str | None) -> tuple[str, str]:
    """Resolve which counter to use. Returns (tool, binary_path)."""
    def _resolve_kplex():
        for cand in (kplex_bin, os.environ.get("KPLEX_BIN"), shutil.which("Kplex")):
            if cand and Path(cand).exists():
                return str(Path(cand).resolve())
        return None

    def _resolve_kmc():
        for cand in (kmc_bin, os.environ.get("KMC_BIN"), shutil.which("kmc")):
            if cand and (Path(cand).exists() or shutil.which(str(cand))):
                return str(cand)
        return None

    if tool in ("auto", "fastk"):
        kp = _resolve_kplex()
        if kp:
            return "fastk", kp
        if tool == "fastk":
            raise SystemExit("Kplex not found. Pass --kplex /path/to/FASTK/Kplex or set KPLEX_BIN.")
    if tool in ("auto", "kmc"):
        km = _resolve_kmc()
        if km:
            return "kmc", km
        if tool == "kmc":
            raise SystemExit("kmc not found. Pass --kmc /path/to/kmc or `conda install -c bioconda kmc`.")
    raise SystemExit(
        "No k-mer counter found. Install FASTK (and pass --kplex .../Kplex or set KPLEX_BIN), "
        "or install KMC (`conda install -c bioconda kmc`)."
    )


def prepare_fasta(fasta: Path, workdir: Path) -> Path:
    """Decompress .gz and give FASTK a friendly extension (.fna is sometimes rejected)."""
    if fasta.suffix == ".gz":
        out = workdir / (fasta.stem if Path(fasta.stem).suffix else fasta.stem + ".fasta")
        with gzip.open(fasta, "rt") as fin, open(out, "wt") as fout:
            shutil.copyfileobj(fin, fout)
        fasta = out
    if fasta.suffix in {".fa", ".fasta"}:
        return fasta
    link = workdir / "genome.fasta"
    if link.exists():
        link.unlink()
    try:
        os.symlink(fasta.resolve(), link)
    except OSError:
        shutil.copyfile(fasta, link)
    return link


def filter_chromosomes(fasta: Path, out_fasta: Path, regex: str) -> dict:
    rx = re.compile(regex)
    stats = dict(contigs_all=0, contigs_chr=0, bases_all=0, bases_chr=0)
    with open(fasta, "rt", errors="ignore") as fin, open(out_fasta, "wt") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                stats["contigs_all"] += 1
                keep = bool(rx.search(line))
                if keep:
                    stats["contigs_chr"] += 1
                    fout.write(line)
                continue
            seq = line.strip()
            stats["bases_all"] += len(seq)
            if keep:
                stats["bases_chr"] += len(seq)
                fout.write(line)
    if stats["contigs_chr"] == 0:
        raise SystemExit(f"--chromosomes-only: no sequences matched /{regex}/")
    return stats


def fasta_lengths(fasta: Path) -> list[int]:
    lengths, cur = [], 0
    with open(fasta, "rt", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur:
        lengths.append(cur)
    return lengths


def theoretical_totals(lengths: list[int], k_values: list[int]) -> dict[int, int]:
    return {k: sum((L - k + 1) for L in lengths if L >= k) for k in k_values}


def add_theoretical(curve_csv: Path, out_csv: Path, t_theory: dict[int, int]) -> None:
    df = pd.read_csv(curve_csv)
    df["k"] = pd.to_numeric(df["k"], errors="coerce").astype("Int64")
    df["total_kmers_theoretical"] = df["k"].map(t_theory)
    df["total_kmers_observed"] = pd.to_numeric(df["total_kmers"], errors="coerce")
    df["fraction_unique_theoretical"] = df["unique_kmers"] / df["total_kmers_theoretical"]
    df["truncation_fraction"] = 1.0 - (df["total_kmers_observed"] / df["total_kmers_theoretical"])
    df["truncation_delta"] = df["total_kmers_theoretical"] - df["total_kmers_observed"]
    df.to_csv(out_csv, index=False)


# ----------------------------------------------------------------------------- counters
def _run_fastk(fasta: Path, kplex_bin: str, k_range: str, h_range: str,
               threads: int, out_csv: Path, workdir: Path, verbose: bool) -> None:
    env = os.environ.copy()
    env["PATH"] = str(Path(kplex_bin).resolve().parent) + os.pathsep + env.get("PATH", "")
    cmd = [kplex_bin, "-i", str(fasta), "-k", k_range, "-h", h_range,
           "-T", str(threads), "-o", str(workdir / "kplex_tmp"), "-c", str(out_csv)]
    if verbose:
        print("  [fastk]", " ".join(cmd))
    subprocess.check_call(cmd, env=env)
    if not out_csv.exists() or out_csv.stat().st_size == 0:
        raise RuntimeError("Kplex produced no output")


def _run_kmc(fasta: Path, kmc_bin: str, k_range: str, threads: int,
             out_csv: Path, workdir: Path, verbose: bool) -> None:
    kmc_tools = str(Path(kmc_bin).parent / "kmc_tools") if Path(kmc_bin).parent.name else "kmc_tools"
    rows = []
    for k in parse_k_range(k_range):
        db = workdir / f"kmc_k{k}"
        tmp = workdir / f"kmc_tmp_k{k}"
        tmp.mkdir(parents=True, exist_ok=True)
        subprocess.run([kmc_bin, f"-k{k}", f"-t{threads}", "-ci1", "-cs1000000000",
                        "-fm", str(fasta), str(db), str(tmp)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        info = subprocess.run([kmc_tools, "info", str(db)], capture_output=True, text=True, check=True)
        total = unique = 0
        for line in info.stdout.splitlines():
            if "Total no. of k-mers" in line:
                total = int(line.split(":")[1].strip())
            elif "Total no. of unique k-mers" in line:
                unique = int(line.split(":")[1].strip())
        rows.append(dict(k=k, unique_kmers=unique, total_kmers=total,
                         fraction_unique=(unique / total if total else 0.0)))
        for f in db.parent.glob(f"{db.name}*"):
            f.unlink()
        shutil.rmtree(tmp, ignore_errors=True)
        if verbose:
            print(f"  [kmc k={k}] unique={unique} total={total}")
    pd.DataFrame(rows).to_csv(out_csv, index=False)


# ----------------------------------------------------------------------------- public API
def compute_curve(fasta: str, out_prefix: str, k_range: str = "5:151:1",
                  tool: str = "auto", kplex_bin: str | None = None, kmc_bin: str | None = None,
                  threads: int = 4, h_range: str = "1:10000000",
                  chromosomes_only: bool = False, chrom_regex: str = DEFAULT_CHROM_REGEX,
                  theoretical: bool = True, verbose: bool = False) -> dict:
    """Run k-plexity on one FASTA. Returns paths of the outputs produced."""
    fasta = Path(fasta)
    if not fasta.exists():
        raise SystemExit(f"FASTA not found: {fasta}")
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tool, binpath = find_counter(tool, kplex_bin, kmc_bin)

    with tempfile.TemporaryDirectory(prefix="kplex_") as td:
        workdir = Path(td)
        seq = prepare_fasta(fasta, workdir)
        if chromosomes_only:
            filt = workdir / "chromosomes.fasta"
            stats = filter_chromosomes(seq, filt, chrom_regex)
            seq = filt
            if verbose:
                print(f"  [chr] kept {stats['contigs_chr']}/{stats['contigs_all']} sequences, "
                      f"{stats['bases_chr']/1e6:.1f} Mb")
        curve_csv = Path(f"{out_prefix}.kplex.csv")
        if tool == "fastk":
            _run_fastk(seq, binpath, k_range, h_range, threads, curve_csv, workdir, verbose)
        else:
            _run_kmc(seq, binpath, k_range, threads, curve_csv, workdir, verbose)

        out = {"tool": tool, "curve": str(curve_csv)}
        if theoretical:
            k_values = parse_k_range(k_range)
            t_theory = theoretical_totals(fasta_lengths(seq), k_values)
            theo_csv = Path(f"{out_prefix}.theoretical.csv")
            add_theoretical(curve_csv, theo_csv, t_theory)
            out["theoretical"] = str(theo_csv)
    return out
