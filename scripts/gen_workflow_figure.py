"""GVClass workflow "subway map" figure for the docs and README.

Standalone: no repo imports, only matplotlib. Regenerate with

    uv run --with matplotlib python scripts/gen_workflow_figure.py

which writes docs/assets/gvclass_workflow.png (300 dpi). The layout follows
the NeLLi workflow-figure style (phase bands, colored rails, station nodes,
legend). Every station states behavior verified in
docs/explanation/how-it-works.md, docs/reference/markers.md,
docs/reference/cli.md, and docs/reference/output.md.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle

OUTPUT = Path(__file__).resolve().parents[1] / "docs" / "assets" / "gvclass_workflow.png"

# Rail identity palette (Okabe-Ito based, CVD-safe).
TRACKS = {
    "genome": "#293845",
    "marker": "#0072B2",
    "placement": "#009E73",
    "quality": "#D55E00",
    "runner": "#6D747C",
}
TEXT = "#26323A"


def add_phase_band(ax, x0, x1, label, color):
    ax.add_patch(
        Rectangle(
            (x0, 0.06),
            x1 - x0,
            6.35,
            facecolor=color,
            edgecolor="white",
            linewidth=1.4,
            alpha=0.20,
            zorder=0,
        )
    )
    ax.text(
        (x0 + x1) / 2,
        6.55,
        label,
        ha="center",
        va="center",
        fontsize=10.5,
        fontweight="bold",
        color="#27323A",
    )


def add_track(ax, points, color, *, linewidth=7, dashed=False):
    xs, ys = zip(*points)
    ax.plot(
        xs,
        ys,
        color=color,
        linewidth=linewidth,
        solid_capstyle="round",
        linestyle=(0, (5, 4)) if dashed else "solid",
        alpha=0.92,
        zorder=2,
    )


def add_station(ax, x, y, label, color, *, dy=0.34, ha="center", transfer=False, fontsize=9.6):
    radius = 0.155 if not transfer else 0.2
    ax.add_patch(
        Circle((x, y), radius, facecolor="white", edgecolor=color, linewidth=2.2, zorder=5)
    )
    if transfer:
        ax.add_patch(
            Circle((x, y), radius + 0.08, facecolor="none", edgecolor=color, linewidth=1.7, zorder=4)
        )
    if label:
        va = "bottom" if dy >= 0 else "top"
        ax.text(
            x,
            y + dy,
            label,
            ha=ha,
            va=va,
            fontsize=fontsize,
            color=TEXT,
            linespacing=0.97,
            zorder=6,
            path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
        )


def add_note(ax, x, y, text, *, ha="center", fontsize=8.8):
    ax.text(
        x,
        y,
        text,
        ha=ha,
        va="top",
        fontsize=fontsize,
        style="italic",
        color="#3A4750",
        zorder=6,
        path_effects=[pe.withStroke(linewidth=2.2, foreground="white")],
    )


def add_transfer(ax, x, y0, y1, color="#A8B0B8"):
    ax.plot([x, x], [y0, y1], color=color, linewidth=1.7, linestyle=(0, (2, 3)), zorder=1)


def add_resource_box(ax):
    box = FancyBboxPatch(
        (4.20, 5.95),
        4.35,
        0.36,
        boxstyle="round,pad=0.05,rounding_size=0.1",
        facecolor="#F2F6FA",
        edgecolor="#8AA2B6",
        linewidth=1.1,
        zorder=3,
    )
    ax.add_patch(box)
    ax.text(
        6.375,
        6.13,
        "GVClass database v2.0.0 · HMMs, references, taxonomy",
        ha="center",
        va="center",
        fontsize=9.0,
        color="#28343C",
        zorder=4,
    )
    # Dashed drops: database feeds the HMM scan and the reference pull.
    ax.plot([4.35, 4.35], [5.95, 4.53], color="#8AA2B6", linewidth=1.1, linestyle=(0, (2, 3)))
    ax.plot([7.20, 7.20], [5.95, 3.26], color="#8AA2B6", linewidth=1.1, linestyle=(0, (2, 3)))


def build():
    mpl.rcParams["font.family"] = "DejaVu Sans"

    fig, ax = plt.subplots(figsize=(12.8, 6.4), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_xlim(0, 14.6)
    ax.set_ylim(-0.85, 6.85)
    ax.axis("off")
    ax.set_position([0.01, 0.02, 0.98, 0.96])

    phases = [
        (1.60, 3.80, "gene calling", "#6F92A8"),
        (3.80, 6.75, "marker detection", "#0072B2"),
        (6.75, 10.90, "per-marker placement", "#009E73"),
        (10.90, 13.95, "classify + quality", "#8E5EA2"),
    ]
    for x0, x1, label, color in phases:
        add_phase_band(ax, x0, x1, label, color)

    add_resource_box(ax)

    # Rails.
    add_track(ax, [(1.0, 5.35), (13.85, 5.35)], TRACKS["genome"], linewidth=7.5)
    add_track(ax, [(3.7, 4.35), (12.45, 4.35)], TRACKS["marker"], linewidth=7)
    add_track(ax, [(5.0, 3.88), (6.3, 3.88)], TRACKS["marker"], linewidth=4.2, dashed=True)
    add_track(ax, [(6.55, 3.08), (12.45, 3.08)], TRACKS["placement"], linewidth=7)
    add_track(ax, [(6.05, 1.80), (12.45, 1.80)], TRACKS["quality"], linewidth=7)
    add_track(ax, [(1.0, 0.55), (13.85, 0.55)], TRACKS["runner"], linewidth=4.2, dashed=True)

    # Transfers (information handoffs between rails). Colored where a rail
    # begins; neutral grey where evidence feeds across.
    add_transfer(ax, 3.7, 4.35, 5.35, color=TRACKS["marker"])
    add_transfer(ax, 5.0, 3.88, 4.35, color=TRACKS["marker"])
    add_transfer(ax, 6.3, 3.88, 4.35, color=TRACKS["marker"])
    add_transfer(ax, 6.55, 3.08, 4.35, color=TRACKS["placement"])
    add_transfer(ax, 6.05, 1.80, 4.35, color=TRACKS["quality"])
    add_transfer(ax, 9.5, 1.80, 4.35)
    add_transfer(ax, 12.45, 1.80, 5.35)

    stations = [
        # Genome / protein rail.
        (1.0, 5.35, "Query genomes\n.fna / .faa", TRACKS["genome"], 0.36, False, 9.6, "center"),
        (2.55, 5.35, "Gene calling\npyrodigal", TRACKS["genome"], 0.36, False, 9.6, "center"),
        (3.7, 5.35, "Query\nproteins", TRACKS["genome"], 0.36, True, 9.6, "center"),
        (12.45, 5.35, "Majority-vote taxonomy\n+ confidence", TRACKS["genome"], 0.40, True, 9.6, "center"),
        (13.72, 5.35, "gvclass_summary.tsv\n44 columns", TRACKS["genome"], -0.44, False, 9.6, "center"),
        # Marker evidence rail.
        (4.35, 4.35, "HMM marker scan\npyhmmer", TRACKS["marker"], 0.30, False, 9.2, "center"),
        (6.55, 4.35, "Marker\nproteins", TRACKS["marker"], 0.30, True, 9.2, "center"),
        (9.5, 4.35, "Copy counts\nper panel", TRACKS["marker"], 0.30, False, 9.2, "center"),
        (5.55, 3.88, "Optional: Extended panel\n576 order-level markers", TRACKS["marker"], -0.27, False, 8.8, "center"),
        # Placement rail.
        (7.2, 3.08, "Reference pull\npyswrd · top 100", TRACKS["placement"], 0.30, False, 9.2, "center"),
        (8.75, 3.08, "Align + trim\npyfamsa · pytrimal", TRACKS["placement"], -0.30, False, 9.2, "center"),
        (10.05, 3.08, "Gene tree\nper marker", TRACKS["placement"], 0.30, False, 9.2, "center"),
        (11.5, 3.08, "Nearest neighbor\n+ distance", TRACKS["placement"], 0.30, False, 9.2, "center"),
        # Quality rail.
        (6.4, 1.80, "Cellular flags\nBUSCO · UNI56", TRACKS["quality"], -0.34, False, 9.2, "center"),
        (10.5, 1.80, "Completeness\nnovelty-aware model", TRACKS["quality"], -0.34, False, 9.2, "center"),
        (12.0, 1.80, "Contamination\ntrained model", TRACKS["quality"], -0.34, False, 9.2, "center"),
        # Runner rail.
        (1.0, 0.55, "pixi run gvclass ·\ngvclass-a (Apptainer)", TRACKS["runner"], -0.22, False, 8.8, "center"),
        (5.5, 0.55, "Parallel per-query workers\n-j workers · -t threads", TRACKS["runner"], -0.22, False, 8.8, "center"),
        (9.3, 0.55, "--resume\nrun_status.json", TRACKS["runner"], -0.22, False, 8.8, "center"),
        (13.72, 0.55, "run.log ·\nper-query archives", TRACKS["runner"], -0.22, False, 8.8, "center"),
    ]
    for x, y, label, color, dy, transfer, fontsize, ha in stations:
        add_station(ax, x, y, label, color, dy=dy, transfer=transfer, fontsize=fontsize, ha=ha)

    # Junctions converging on the classification transfer.
    for y, color in ((4.35, TRACKS["marker"]), (3.08, TRACKS["placement"]), (1.80, TRACKS["quality"])):
        add_station(ax, 12.45, y, "", color, dy=0.0)

    # Italic notes: facts that belong to a station but not in its name.
    add_note(ax, 1.0, 5.02, "(.faa skips gene calling)")
    add_note(ax, 2.55, 4.88, "9 genetic codes · best by\nmarkers, score, or density")
    add_note(ax, 3.6, 4.12, "Nucleocytoviricota · Mirusviricota\nPreplasmiviricota · cellular")
    add_note(ax, 11.3, 2.60, "Optional: IQ-TREE (--tree-method) ·\nsupermatrix & species tree (--species-tree)")

    legend_items = [
        ("Genome / protein flow", TRACKS["genome"], "-"),
        ("Marker evidence (HMM)", TRACKS["marker"], "-"),
        ("Phylogenetic placement", TRACKS["placement"], "-"),
        ("Quality models", TRACKS["quality"], "-"),
        ("Execution / resume", TRACKS["runner"], (0, (5, 4))),
    ]
    handles = [
        Line2D([0], [0], color=color, lw=4, linestyle=ls, solid_capstyle="round", label=label)
        for label, color, ls in legend_items
    ]
    ax.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.01),
        ncol=3,
        frameon=False,
        fontsize=9.6,
        handlelength=2.2,
        columnspacing=1.7,
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=300, facecolor="white")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    mpl.use("Agg")
    build()
