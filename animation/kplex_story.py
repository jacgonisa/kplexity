"""Manim animation explaining the k-plexity CLI tool.

Render:
    mamba run -n manim manim -qh animation/kplex_story.py KplexStory
(uses real fitted double-sigmoid params from an Abax parallelepipedus genome)
"""
import csv
import os
import numpy as np
from manim import *

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def load_xy(name):
    """Read (k, fraction_unique) from a real .kplex.csv, sorted by k."""
    pts = []
    with open(os.path.join(DATA, f"{name}.kplex.csv")) as f:
        for row in csv.DictReader(f):
            try:
                pts.append((float(row["k"]), float(row["fraction_unique"])))
            except (ValueError, KeyError):
                continue
    pts.sort()
    return [p[0] for p in pts], [p[1] for p in pts]

# real double-sigmoid fit (kplex fit on Abax_parallelepipedus.chr.theoretical.csv)
L1, S1, K01 = 0.5821, 0.6973, 16.09
L2, S2, K02 = 0.3232, 0.0749, 42.24
ASYMPTOTE = L1 + L2

BG = "#0d1117"
CYAN = "#4dd0e1"
AMBER = "#ffca28"
CORAL = "#ef5350"
GREEN = "#66bb6a"
GREY = "#8899a6"


def sig(x, L, s, k0):
    return L / (1.0 + np.exp(-s * (x - k0)))


def double_sig(x):
    return sig(x, L1, S1, K01) + sig(x, L2, S2, K02)


class KplexStory(Scene):
    def construct(self):
        self.camera.background_color = BG
        self.title_card()
        self.concept()
        self.definition()
        self.curve_scene()
        self.decomposition()
        self.cli_scene()
        self.outro()

    # ------------------------------------------------------------------ 1. title
    def title_card(self):
        title = Text("k-plexity", font="sans-serif", weight=BOLD, color=CYAN).scale(1.9)
        sub = Text("the fraction of unique k-mers, across k",
                   font="sans-serif", color=GREY).scale(0.55)
        sub.next_to(title, DOWN, buff=0.35)
        self.play(Write(title), run_time=1.2)
        self.play(FadeIn(sub, shift=UP * 0.2))
        self.wait(0.8)
        self.play(VGroup(title, sub).animate.scale(0.42).to_edge(UP), run_time=0.9)
        self.wait(0.2)
        self.header = VGroup(title, sub)

    # ------------------------------------------------------------------ 2. concept
    def concept(self):
        q = Text("How unique are a genome's k-mers?",
                 font="sans-serif", color=WHITE).scale(0.6).next_to(self.header, DOWN, buff=0.5)
        self.play(FadeIn(q, shift=UP * 0.2))

        seq = "GATTACA GATTACA CCG".replace(" ", "")
        letters = VGroup(*[
            Text(c, font="monospace", weight=BOLD,
                 color={"G": GREEN, "A": AMBER, "T": CORAL, "C": CYAN}[c]).scale(0.75)
            for c in seq
        ]).arrange(RIGHT, buff=0.14).move_to(ORIGIN + DOWN * 0.4)
        self.play(LaggedStart(*[FadeIn(l) for l in letters], lag_ratio=0.05, run_time=1.2))

        k = 3
        win = SurroundingRectangle(letters[0:k], color=WHITE, buff=0.06, corner_radius=0.05)
        cap = Text("slide a k-mer window (k=3)", font="sans-serif", color=GREY).scale(0.45)
        cap.next_to(letters, DOWN, buff=0.5)
        self.play(Create(win), FadeIn(cap))

        # slide to the repeated GAT window and flag the repeat
        for start in range(1, 8):
            self.play(win.animate.move_to(letters[start:start + k].get_center()), run_time=0.22)
        rep1 = SurroundingRectangle(letters[0:k], color=CORAL, buff=0.06, corner_radius=0.05)
        rep2 = SurroundingRectangle(letters[7:7 + k], color=CORAL, buff=0.06, corner_radius=0.05)
        flag = Text("same k-mer seen twice → not unique",
                    font="sans-serif", color=CORAL).scale(0.5).next_to(letters, DOWN, buff=0.5)
        self.play(FadeOut(win), Create(rep1), Create(rep2), FadeTransform(cap, flag))
        self.wait(0.6)

        takeaway = VGroup(
            Text("small k  →  many shared k-mers  →  low uniqueness",
                 font="sans-serif", color=CORAL).scale(0.5),
            Text("large k  →  k-mers become unique  →  high uniqueness",
                 font="sans-serif", color=GREEN).scale(0.5),
        ).arrange(DOWN, buff=0.25).next_to(letters, DOWN, buff=0.5)
        self.play(FadeOut(flag), FadeIn(takeaway[0]))
        self.play(FadeIn(takeaway[1]))
        self.wait(1.0)
        self.play(FadeOut(letters, rep1, rep2, takeaway, q))

    # ------------------------------------------------------- 2b. how it's calculated
    def definition(self):
        # the formula, right under the header
        formula = Text("fraction unique  =  U / T", font="sans-serif",
                       weight=BOLD, color=WHITE).scale(0.8).next_to(self.header, DOWN, buff=0.45)
        self.play(Write(formula), run_time=1.0)
        defs = VGroup(
            Text("U  =  number of distinct k-mers", font="sans-serif", color=GREEN).scale(0.5),
            Text("T  =  total k-mers  (one per position, ≈ genome length)",
                 font="sans-serif", color=CYAN).scale(0.5),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2).next_to(formula, DOWN, buff=0.35)
        self.play(FadeIn(defs[0], shift=UP * 0.15))
        self.play(FadeIn(defs[1], shift=UP * 0.15))
        self.wait(0.7)
        self.play(FadeOut(defs))

        # tiny worked example at k = 3
        seq = "GATTACAGAT"
        cmap = {"G": GREEN, "A": AMBER, "T": CORAL, "C": CYAN}
        letters = VGroup(*[
            Text(c, font="monospace", weight=BOLD, color=cmap[c]).scale(0.62) for c in seq
        ]).arrange(RIGHT, buff=0.12).next_to(formula, DOWN, buff=0.6)
        exlab = Text("list every k-mer (k = 3):", font="sans-serif", color=GREY).scale(0.44)
        exlab.next_to(letters, DOWN, buff=0.3)
        self.play(LaggedStart(*[FadeIn(l) for l in letters], lag_ratio=0.05, run_time=0.9),
                  FadeIn(exlab))

        kmers = [seq[i:i + 3] for i in range(len(seq) - 2)]           # 8 windows
        tri = VGroup(*[Text(m, font="monospace", color=WHITE).scale(0.5) for m in kmers])
        tri.arrange(RIGHT, buff=0.32).next_to(exlab, DOWN, buff=0.35)
        self.play(LaggedStart(*[FadeIn(t) for t in tri], lag_ratio=0.1, run_time=1.3))

        # flag the one duplicate (GAT at both ends)
        dup = VGroup(
            SurroundingRectangle(tri[0], color=CORAL, buff=0.05, corner_radius=0.05),
            SurroundingRectangle(tri[7], color=CORAL, buff=0.05, corner_radius=0.05),
        )
        duplab = Text("GAT occurs twice → counts once toward U",
                      font="sans-serif", color=CORAL).scale(0.44).next_to(tri, DOWN, buff=0.35)
        self.play(Create(dup), FadeIn(duplab))
        self.wait(0.6)

        calc = VGroup(
            Text("T = 8", font="sans-serif", weight=BOLD, color=CYAN).scale(0.55),
            Text("U = 7", font="sans-serif", weight=BOLD, color=GREEN).scale(0.55),
            Text("→   U / T = 0.88", font="sans-serif", weight=BOLD, color=AMBER).scale(0.55),
        ).arrange(RIGHT, buff=0.7).next_to(duplab, DOWN, buff=0.4)
        self.play(FadeIn(calc, shift=UP * 0.15))
        self.wait(1.3)
        self.play(FadeOut(formula, letters, exlab, tri, dup, duplab, calc))

    # ------------------------------------------------------------------ 3. curve
    def curve_scene(self):
        ax = Axes(
            x_range=[0, 150, 30], y_range=[0, 1, 0.25],
            x_length=9.2, y_length=4.6,
            axis_config={"color": GREY, "stroke_width": 2},
            tips=False,
        ).to_edge(DOWN, buff=0.7)
        # manual Pango tick labels (avoids LaTeX)
        ticks = VGroup()
        for xv in (5, 50, 100, 150):
            t = Text(str(xv), font="sans-serif", color=GREY).scale(0.35)
            t.next_to(ax.c2p(xv, 0), DOWN, buff=0.12)
            ticks.add(t)
        for yv in (0.0, 0.5, 1.0):
            t = Text(f"{yv:.1f}", font="sans-serif", color=GREY).scale(0.35)
            t.next_to(ax.c2p(0, yv), LEFT, buff=0.12)
            ticks.add(t)
        xlab = Text("k-mer length  k", font="sans-serif", color=GREY).scale(0.5)
        xlab.next_to(ax.x_axis, DOWN, buff=0.5)
        ylab = Text("fraction unique  (U / T)", font="sans-serif", color=GREY).scale(0.5)
        ylab.rotate(PI / 2).next_to(ax.y_axis, LEFT, buff=0.5)
        self.play(Create(ax), FadeIn(ticks), FadeIn(xlab), FadeIn(ylab))
        self.ticks = ticks

        curve = ax.plot(double_sig, x_range=[5, 150, 0.5], color=CYAN, stroke_width=5)
        tracker = ValueTracker(5)
        dot = always_redraw(lambda: Dot(color=WHITE, radius=0.07).move_to(
            ax.c2p(tracker.get_value(), double_sig(tracker.get_value()))))
        self.add(dot)
        self.play(Create(curve), tracker.animate.set_value(150), run_time=2.6)
        self.play(FadeOut(dot))

        note = Text("one curve  —  one genome",
                    font="sans-serif", color=WHITE).scale(0.5)
        note.next_to(ax, UP, buff=0.2)
        self.play(FadeIn(note, shift=UP * 0.15))
        self.wait(0.8)
        self.play(FadeOut(note))
        self.ax, self.curve, self.axis_labels = ax, curve, VGroup(xlab, ylab)

    # ------------------------------------------------------------------ 4. fit / decomposition
    def decomposition(self):
        ax = self.ax
        title = Text("kplex fit  —  a double sigmoid", font="sans-serif",
                     color=AMBER, weight=BOLD).scale(0.55).next_to(ax, UP, buff=0.2)
        self.play(FadeIn(title))

        comp1 = ax.plot(lambda x: sig(x, L1, S1, K01), x_range=[5, 150, 0.5],
                        color=GREEN, stroke_width=3)
        comp2 = ax.plot(lambda x: sig(x, L2, S2, K02), x_range=[5, 150, 0.5],
                        color=CORAL, stroke_width=3)
        c1lab = Text("component 1  (L1)", font="sans-serif", color=GREEN).scale(0.38)
        c1lab.next_to(ax.c2p(150, sig(150, L1, S1, K01)), RIGHT, buff=0.1)
        c2lab = Text("component 2  (L2)", font="sans-serif", color=CORAL).scale(0.38)
        c2lab.next_to(ax.c2p(150, sig(150, L2, S2, K02)), RIGHT, buff=0.1)
        self.play(Create(comp1), FadeIn(c1lab))
        self.play(Create(comp2), FadeIn(c2lab))

        asym = DashedLine(ax.c2p(0, ASYMPTOTE), ax.c2p(150, ASYMPTOTE), color=GREY)
        asym_lab = Text(f"asymptote = {ASYMPTOTE:.2f}", font="sans-serif", color=GREY).scale(0.42)
        asym_lab.next_to(ax.c2p(0, ASYMPTOTE), UR, buff=0.08)
        self.play(Create(asym), FadeIn(asym_lab))
        self.wait(0.4)

        params = VGroup(
            Text(f"asymptote  L1 + L2  =  {ASYMPTOTE:.2f}",
                 font="sans-serif", color=WHITE).scale(0.42),
            Text(f"L1 = {L1:.2f}      L2 = {L2:.2f}",
                 font="sans-serif", color=WHITE).scale(0.42),
            Text(f"k0,1 = {K01:.0f}       k0,2 = {K02:.0f}",
                 font="sans-serif", color=WHITE).scale(0.42),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.18)
        bg = BackgroundRectangle(params, color=BG, fill_opacity=0.85, buff=0.28)
        box = SurroundingRectangle(params, color=AMBER, buff=0.22, corner_radius=0.1)
        panel = VGroup(bg, params, box).move_to(ax.c2p(80, 0.36))
        self.play(FadeOut(title), FadeIn(panel))
        self.wait(1.6)
        self.play(FadeOut(comp1, comp2, c1lab, c2lab, asym, asym_lab, panel,
                          self.ax, self.curve, self.axis_labels, self.ticks))

    # ------------------------------------------------------------------ 5. CLI
    def cli_scene(self):
        head = Text("three commands", font="sans-serif", color=WHITE).scale(0.6)
        head.next_to(self.header, DOWN, buff=0.6)
        self.play(FadeIn(head, shift=UP * 0.2))

        def cmd(cmd_txt, desc, col):
            c = Text(cmd_txt, font="monospace", weight=BOLD, color=col).scale(0.62)
            d = Text(desc, font="sans-serif", color=GREY).scale(0.44)
            d.next_to(c, RIGHT, buff=0.5)
            return VGroup(c, d)

        rows = VGroup(
            cmd("kplex run", "FASTA  →  k-plexity curve  (FASTK / KMC)", CYAN),
            cmd("kplex fit", "curve  →  double-sigmoid parameters", AMBER),
            cmd("kplex plot", "curve + fit  →  publication figure", GREEN),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.55).move_to(ORIGIN)
        # left-align the command column
        maxx = max(r[0].get_left()[0] for r in rows)
        for r in rows:
            r.shift(RIGHT * (rows[0][0].get_left()[0] - r[0].get_left()[0]))

        for r in rows:
            self.play(FadeIn(r[0], shift=RIGHT * 0.3), run_time=0.5)
            self.play(FadeIn(r[1]), run_time=0.4)
        self.wait(1.2)
        self.play(FadeOut(rows, head))

    # ------------------------------------------------------------------ 6. outro
    def outro(self):
        install = Text("pip install -e .", font="monospace", color=GREEN).scale(0.7)
        gh = Text("github.com/jacgonisa/kplexity", font="sans-serif", color=GREY).scale(0.5)
        grp = VGroup(install, gh).arrange(DOWN, buff=0.4).move_to(ORIGIN)
        self.play(FadeIn(install, scale=1.1))
        self.play(FadeIn(gh, shift=UP * 0.2))
        self.wait(1.5)
        self.play(FadeOut(grp, self.header))
        self.wait(0.3)


def _curve_end(curve, ax):
    """Return (x,y) of the current end of a partially-drawn curve, in data coords."""
    pt = curve.get_end()
    x, y = ax.p2c(pt)
    return x, y
