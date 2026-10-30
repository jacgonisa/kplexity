"""Manim animation explaining the k-plexity CLI tool.

Render:
    mamba run -n manim manim -qh animation/kplex_story.py KplexStory
(uses a real curve and its kplexity-model fit: Abax parallelepipedus, chromosomes only)
"""
import numpy as np
from manim import *

# real curve, odd k from 5 to 151 (Abax_parallelepipedus.chr.theoretical.csv, exact N-free T)
DATA_K = np.arange(5, 152, 2)
DATA_Y = np.array([
    0.0000, 0.0000, 0.0002, 0.0036, 0.0529, 0.2493, 0.4378, 0.5352, 0.5811, 0.6095, 0.6324, 0.6529,
    0.6715, 0.6887, 0.7044, 0.7190, 0.7324, 0.7449, 0.7564, 0.7670, 0.7769, 0.7861, 0.7945, 0.8024,
    0.8097, 0.8165, 0.8228, 0.8287, 0.8342, 0.8393, 0.8440, 0.8485, 0.8527, 0.8566, 0.8603, 0.8637,
    0.8670, 0.8701, 0.8730, 0.8757, 0.8783, 0.8808, 0.8832, 0.8854, 0.8876, 0.8896, 0.8915, 0.8934,
    0.8952, 0.8969, 0.8986, 0.9001, 0.9017, 0.9031, 0.9045, 0.9059, 0.9072, 0.9085, 0.9097, 0.9109,
    0.9120, 0.9131, 0.9142, 0.9153, 0.9163, 0.9173, 0.9182, 0.9192, 0.9201, 0.9209, 0.9218, 0.9226,
    0.9234, 0.9242])

# kplexity-model fit of that curve (kplex fit): y = B(k) * [A - S exp(-lam (k - k0))]
B_STEEP, K0, A, S, LAM = 0.7155, 14.88, 0.9219, 0.4025, 0.0331
LOG4G = 14.53                      # log4 of the 579.6 Mb of k-mer positions
DIVERGENCE = 1 - np.exp(-LAM)

BG = "#0d1117"
CYAN = "#4dd0e1"
AMBER = "#ffca28"
CORAL = "#ef5350"
GREEN = "#66bb6a"
GREY = "#8899a6"


def first_rise(x):
    """Chance-collision term B(k): (1 - e^-L) / L with L = 2 exp(-b (k - k0))."""
    L = 2.0 * np.exp(-B_STEEP * (x - K0))
    return -np.expm1(-L) / L


def collisions_plateau(x):
    return first_rise(x) * A


def kplexity(x):
    return first_rise(x) * (A - S * np.exp(-LAM * (x - K0)))


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
        ).to_edge(DOWN, buff=0.95)
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
        xlab.next_to(ax.x_axis, DOWN, buff=0.38)
        ylab = Text("fraction unique  (U / T)", font="sans-serif", color=GREY).scale(0.5)
        ylab.rotate(PI / 2).next_to(ax.y_axis, LEFT, buff=0.5)
        self.play(Create(ax), FadeIn(ticks), FadeIn(xlab), FadeIn(ylab))
        self.ticks = ticks

        dots = VGroup(*[Dot(ax.c2p(k, y), radius=0.035, color=WHITE)
                        for k, y in zip(DATA_K, DATA_Y)])
        self.play(LaggedStart(*[FadeIn(d, scale=0.5) for d in dots], lag_ratio=0.03, run_time=2.2))
        curve = ax.plot(kplexity, x_range=[5, 151, 0.25], color=CYAN, stroke_width=4)
        self.dots = dots

        note = Text("one genome  —  Abax parallelepipedus (ground beetle)",
                    font="sans-serif", color=WHITE).scale(0.5)
        note.next_to(ax, UP, buff=0.2)
        self.play(FadeIn(note, shift=UP * 0.15))
        self.wait(0.8)
        self.play(FadeOut(note))
        self.ax, self.curve, self.axis_labels = ax, curve, VGroup(xlab, ylab)

    # ------------------------------------------------------------------ 4. fit / the kplexity model
    def decomposition(self):
        ax = self.ax
        title = Text("kplex fit  —  the kplexity model", font="sans-serif",
                     color=AMBER, weight=BOLD).scale(0.5)
        formula = Text("y(k)  =  B(k) · [ A  −  S · e^(−λ (k − k₀)) ]",
                       font="sans-serif", color=WHITE).scale(0.5)
        top = VGroup(title, formula).arrange(DOWN, buff=0.15).next_to(ax, UP, buff=0.15)
        self.play(FadeIn(title), Create(self.curve), run_time=1.4)
        self.play(Write(formula), run_time=1.0)

        def caption(lines, col):
            g = VGroup(*[Text(t, font="sans-serif", color=col).scale(0.4) for t in lines])
            g.arrange(DOWN, aligned_edge=LEFT, buff=0.12)
            return g.move_to(ax.c2p(102, 0.36))

        # 1. chance collisions: the rise every genome shares, at k0 ~ log4 G
        rise = ax.plot(collisions_plateau, x_range=[5, 151, 0.25], color=GREEN, stroke_width=3)
        k0line = DashedLine(ax.c2p(K0, 0), ax.c2p(K0, A), color=GREEN, stroke_width=2)
        k0lab = Text(f"k₀ = {K0:.1f}", font="sans-serif", color=GREEN).scale(0.36)
        k0lab.next_to(ax.c2p(K0, 0.08), RIGHT, buff=0.12)
        cap1 = caption(["B(k): chance collisions, as in random sequence",
                        f"rise at k₀ ≈ log₄(genome size) = {LOG4G:.1f}"], GREEN)
        self.play(Create(rise), Create(k0line), FadeIn(k0lab), FadeIn(cap1), run_time=1.4)
        self.wait(1.2)

        # 2. near copies: the gap that closes as k grows
        gap = ax.get_area(rise, x_range=[K0, 151], bounded_graph=self.curve,
                          color=CORAL, opacity=0.45)
        cap2 = caption([f"S = {S:.2f}: near copies (repeats) at k₀",
                        "they become distinct as k grows,",
                        f"at rate λ  (divergence {100 * DIVERGENCE:.1f}%)"], CORAL)
        self.play(FadeIn(gap), FadeTransform(cap1, cap2), run_time=1.0)
        self.wait(1.6)

        # 3. identical copies: never resolved at any k
        aline = DashedLine(ax.c2p(5, A), ax.c2p(151, A), color=AMBER, stroke_width=2)
        band = Rectangle(width=ax.c2p(151, 0)[0] - ax.c2p(5, 0)[0],
                         height=ax.c2p(0, 1)[1] - ax.c2p(0, A)[1],
                         stroke_width=0, fill_color=AMBER, fill_opacity=0.3)
        band.move_to(ax.c2p(78, (1 + A) / 2))
        cap3 = caption([f"A = {A:.2f}: plateau, the share that can be resolved",
                        f"1 − A = {1 - A:.2f}: identical copies",
                        "(no k tells them apart)"], AMBER)
        self.play(Create(aline), FadeIn(band), FadeTransform(cap2, cap3), run_time=1.0)
        self.wait(1.6)

        params = VGroup(
            Text(f"A  =  {A:.2f}          1 − A  =  {1 - A:.2f}", font="sans-serif", color=AMBER).scale(0.4),
            Text(f"S  =  {S:.2f}          divergence  =  {100 * DIVERGENCE:.1f}%",
                 font="sans-serif", color=CORAL).scale(0.4),
            Text(f"k₀  =  {K0:.1f}          b  =  {B_STEEP:.2f}", font="sans-serif", color=GREEN).scale(0.4),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.18)
        bg = BackgroundRectangle(params, color=BG, fill_opacity=0.85, buff=0.25)
        box = SurroundingRectangle(params, color=WHITE, buff=0.2, corner_radius=0.1)
        panel = VGroup(bg, params, box).move_to(ax.c2p(102, 0.36))
        self.play(FadeOut(cap3), FadeIn(panel))
        self.wait(2.0)
        self.play(FadeOut(title, formula, rise, k0line, k0lab, gap, aline, band, panel, self.dots,
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
            cmd("kplex run", "FASTA  →  k-plexity curve  (FastK / KMC; --sketch: 15 k)", CYAN),
            cmd("kplex fit", "curve  →  A, S, divergence, k₀  (kplexity model)", AMBER),
            cmd("kplex plot", "curve + fit  →  figure", GREEN),
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

