"""kplex command-line interface: `kplex run|fit|plot`."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

from . import __version__

MODELS = ["kplexity", "gamma", "double", "single"]


def _default_prefix(fasta: str) -> str:
    stem = Path(fasta).name
    for ext in (".gz", ".fasta", ".fa", ".fna"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    return stem


def cmd_run(a) -> int:
    from .kmers import compute_curve, DEFAULT_CHROM_REGEX
    prefix = a.out or _default_prefix(a.fasta)
    out = compute_curve(
        fasta=a.fasta, out_prefix=prefix, k_range=a.k, tool=a.tool,
        kplex_bin=a.kplex, kmc_bin=a.kmc, threads=a.threads, h_range=a.h,
        chromosomes_only=a.chromosomes_only, chrom_regex=a.chrom_regex or DEFAULT_CHROM_REGEX,
        theoretical=not a.no_theoretical, verbose=a.verbose, sketch=a.sketch)
    print(f"[kplex run] tool={out['tool']}  curve -> {out['curve']}")
    if a.sketch:
        print(f"[kplex run] sketch, k = {out['k']}")
    if "theoretical" in out:
        print(f"[kplex run] theoretical -> {out['theoretical']}")
    curve_for_next = out.get("theoretical", out["curve"])
    if a.fit or a.plot:
        from .model import load_curve, fit_curve, genome_bp_from_csv
        x, y, _ = load_curve(curve_for_next)
        fit = fit_curve(x, y, a.model, genome_bp_from_csv(curve_for_next))
        fit_json = f"{prefix}.fit.json"
        Path(fit_json).write_text(json.dumps(fit, indent=2))
        print(f"[kplex run] fit -> {fit_json}")
        _summary("[kplex run]", fit)
        if a.plot:
            from .plotting import plot
            png = plot(curve_for_next, f"{prefix}.png", fit=fit)
            print(f"[kplex run] plot -> {png}")
    return 0


def _summary(tag: str, fit: dict) -> None:
    if "error" in fit:
        print(f"{tag} fit FAILED: {fit['error']}", file=sys.stderr)
        return
    print(f"  R²={fit['r_squared']:.4f}  RMSE={fit['rmse']:.4f}  AICc={fit.get('aicc', float('nan')):.1f}")
    if fit["model"].startswith("kplexity"):
        print(f"  A={fit['A']:.3f} (identical 1-A={fit['identical']:.3f})  S={fit['S']:.3f}  "
              f"divergence={100 * fit['divergence']:.1f}%  k0={fit['k0']:.1f}  b={fit['b']:.2f}"
              + (f"  a={fit['a']:.2f}" if "a" in fit else ""))
        if not fit["plateau_reached"]:
            print(f"  warning: the near-copy gap is still {fit['gap_at_kmax']:.3f} at the largest k, so A (and 1-A) "
                  "is an extrapolation; extend k or read y at the largest k instead", file=sys.stderr)
    elif fit["model"] == "double_sigmoid":
        print(f"  asymptote(L1+L2)={fit['asymptote']:.3f}  L2/(L1+L2)={fit['L2_fraction']:.3f}  "
              f"k0_1={fit['k0_1']:.1f}  k0_2={fit['k0_2']:.1f}")
    else:
        print(f"  L={fit['L']:.3f}  s={fit['s']:.3f}  k0={fit['k0']:.1f}")


def cmd_fit(a) -> int:
    from .model import load_curve, fit_curve, genome_bp_from_csv
    x, y, ycol = load_curve(a.csv)
    G = a.genome_bp if a.genome_bp else genome_bp_from_csv(a.csv)
    fit = fit_curve(x, y, a.model, G)
    prefix = a.out or _default_prefix(a.csv)
    Path(f"{prefix}.fit.json").write_text(json.dumps(fit, indent=2))
    if "error" in fit:
        print(f"[kplex fit] FAILED: {fit['error']}", file=sys.stderr)
        return 1
    print(f"[kplex fit] {fit['model']}  (y={ycol}, n={fit['n_points']})")
    _summary("[kplex fit]", fit)
    print(f"[kplex fit] -> {prefix}.fit.json")
    if a.plot:
        from .plotting import plot
        print(f"[kplex fit] plot -> {plot(a.csv, f'{prefix}.png', fit=fit)}")
    return 0


def cmd_schedule(a) -> int:
    from .sketch import size_anchored
    G = a.genome_bp
    if G is None:
        from .kmers import fasta_lengths, prepare_fasta
        import tempfile
        with tempfile.TemporaryDirectory(prefix="kplex_") as td:
            G = sum(fasta_lengths(prepare_fasta(Path(a.fasta), Path(td))))
    print(",".join(map(str, size_anchored(G, a.n))))
    return 0


def cmd_plot(a) -> int:
    from .plotting import plot
    fit = None
    if a.fit:
        fit = json.loads(Path(a.fit).read_text())
    out = a.out or (_default_prefix(a.csv) + ".png")
    print(f"[kplex plot] -> {plot(a.csv, out, fit=fit, title=a.title)}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="kplex", description="k-plexity: unique-k-mer genome QC")
    p.add_argument("--version", action="version", version=f"kplex {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="compute the k-plexity curve from a FASTA")
    r.add_argument("fasta", help="input genome FASTA (.fa/.fasta/.fna, optionally .gz)")
    r.add_argument("-o", "--out", help="output prefix (default: FASTA name)")
    r.add_argument("-k", "--k", default="5:151:1", help="k values: start:end:step (default 5:151:1) or a comma list")
    r.add_argument("--tool", choices=["auto", "fastk", "kmc"], default="auto")
    r.add_argument("--kplex", help="path to FASTK/Kplex (else $KPLEX_BIN or PATH)")
    r.add_argument("--kmc", help="path to kmc (else $KMC_BIN or PATH)")
    r.add_argument("-T", "--threads", type=int, default=4)
    r.add_argument("--h", default="1:10000000", help="FASTK histogram range (default 1:10000000)")
    r.add_argument("--chromosomes-only", action="store_true", help="keep only chromosome-scale sequences")
    r.add_argument("--chrom-regex", default="", help="regex for chromosome headers")
    r.add_argument("--no-theoretical", action="store_true", help="skip theoretical-fraction CSV")
    r.add_argument("--sketch", type=int, nargs="?", const=15, default=None, metavar="N",
                   help="count only N k values (default 15) on a size-anchored schedule instead of every k in -k; "
                        "about ten times faster, median parameter drift about 2%% at N = 15")
    r.add_argument("--fit", action="store_true", help="also fit the curve (see --model)")
    r.add_argument("--model", choices=MODELS, default="kplexity",
                   help="kplexity (default, 5 parameters), gamma (6), double (double sigmoid) or single")
    r.add_argument("--plot", action="store_true", help="also render a PNG (implies --fit)")
    r.add_argument("-v", "--verbose", action="store_true")
    r.set_defaults(func=cmd_run)

    f = sub.add_parser("fit", help="fit a model to a curve CSV (default: the kplexity model)")
    f.add_argument("csv", help="curve CSV from `kplex run` (.kplex.csv or .theoretical.csv)")
    f.add_argument("-o", "--out", help="output prefix (default: CSV name)")
    f.add_argument("--model", choices=MODELS, default="kplexity",
                   help="kplexity (default, 5 parameters), gamma (6), double (double sigmoid) or single")
    f.add_argument("--genome-bp", type=float, default=None,
                   help="k-mer positions (about the assembly size) to seed k0 at log4(G); "
                        "default: read from the curve's total_kmers column")
    f.add_argument("--plot", action="store_true", help="also render a PNG with the fit overlaid")
    f.set_defaults(func=cmd_fit)

    sc = sub.add_parser("schedule", help="print the size-anchored sketch k list (for Kplex -k or kplex run -k)")
    g = sc.add_mutually_exclusive_group(required=True)
    g.add_argument("--genome-bp", type=float, help="genome size in bp")
    g.add_argument("--fasta", help="FASTA to measure (ACGT bases)")
    sc.add_argument("-n", type=int, default=15, help="number of k values (default 15, minimum 10)")
    sc.set_defaults(func=cmd_schedule)

    pl = sub.add_parser("plot", help="plot a curve, optionally with a fit overlay")
    pl.add_argument("csv", help="curve CSV from `kplex run`")
    pl.add_argument("-o", "--out", help="output PNG (default: CSV name + .png)")
    pl.add_argument("--fit", help="fit JSON from `kplex fit` to overlay")
    pl.add_argument("--title", help="plot title")
    pl.set_defaults(func=cmd_plot)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
