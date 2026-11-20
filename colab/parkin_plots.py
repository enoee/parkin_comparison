"""Colours and matplotlib drawing for the Parkin MAVE comparison.

draw_domains / draw_ss / draw_zinc / draw_heatmap / draw_median_track each draw one track
on one Axes. heatmap_figure assembles them into an A4-width figure split into panels, from
the per-residual data returned by parkin_data.residual_data; update_heatmap_figure redraws
it in place for another residual.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import TwoSlopeNorm, to_rgb
from matplotlib.patches import Rectangle, FancyArrow
from matplotlib.lines import Line2D

import parkin_data as pdata

A4_WIDTH_IN = 7.0   # full text width of an A4 page (~178 mm)

# ------------------------------------------------------------------------- colours
HEATMAP_CMAP = "RdBu_r"
MISSING_COLOR = "0.35"
HOVER_COLOR = "black"      # explorer: variant under the mouse
SELECT_COLOR = "#00a651"   # explorer: selected residue

DOMAIN_COLORS = {  # as in the Parkin domain-architecture figure; ACT is not shown there and
    # is dark grey so it stands out from near-zero (white) residuals on the structure
    "UBL": "#ff0000", "ACT": "#4d4d4d", "RING0": "#7dc21e", "RING1": "#9ed9e3",
    "IBR": "#e01e8c", "REP": "#ffff00", "RING2": "#ff8000",
}


def text_color_on(background):
    """Black or white text, whichever reads better on the given background colour."""
    r, g, b = to_rgb(background)
    return "white" if 0.299 * r + 0.587 * g + 0.114 * b < 0.5 else "black"


def heatmap_norm(vlim=0.5):
    return TwoSlopeNorm(vcenter=0.0, vmin=-vlim, vmax=vlim)


def heatmap_cmap():
    cm = plt.get_cmap(HEATMAP_CMAP).copy()
    cm.set_bad(MISSING_COLOR)
    return cm


# ------------------------------------------------------------------ track drawing
# Tracks share the residue axis (x); the y-axis of annotation tracks runs over [0, 1].
def _setup_track(ax, start, end):
    ax.set_xlim(start - 0.5, end + 0.5)
    ax.set_ylim(0, 1)
    ax.set_axis_off()


def _in_range(s, e, start, end):
    return e >= start and s <= end


def draw_domains(ax, segments, start, end, colors=DOMAIN_COLORS, height=0.8, fontsize=6,
                 label_min_len=12):
    """Domain boxes with labels along the sequence; a thin line marks the backbone."""
    _setup_track(ax, start, end)
    ax.plot([start - 0.5, end + 0.5], [0.5, 0.5], color="0.5", lw=0.6, zorder=0)
    lo = 0.5 - height / 2
    for s, e, name in segments:
        if not _in_range(s, e, start, end):
            continue
        s0, e0 = max(s, start) - 0.5, min(e, end) + 0.5
        fc = colors.get(name, "0.85")
        ax.add_patch(Rectangle((s0, lo), e0 - s0, height, fc=fc, ec="0.2", lw=0.5))
        if e0 - s0 >= label_min_len:
            ax.text((s0 + e0) / 2, 0.5, name, ha="center", va="center", fontsize=fontsize,
                    color=text_color_on(fc))


def draw_ss(ax, segments, start, end, helix_color="#c0392b", strand_color="#2c6fbb",
            coil_color="0.4", helix_period=3.6, helix_amp=0.32, helix_lw=1.0,
            strand_width=0.45, head_width=0.8, head_len=1.5, coil_lw=0.6):
    """Secondary structure: coil as a line, helices as a coil spring, strands as arrows.

    helix_period: residues per turn of the drawn coil; helix_amp: amplitude (track units)
    strand_width / head_width: arrow body / head width (track units); head_len in residues
    """
    _setup_track(ax, start, end)
    ax.plot([start - 0.5, end + 0.5], [0.5, 0.5], color=coil_color, lw=coil_lw, zorder=0)
    for s, e, kind in segments:
        if not _in_range(s, e, start, end):
            continue
        s0, e0 = max(s, start) - 0.5, min(e, end) + 0.5
        if kind == "helix":
            p = np.linspace(s0, e0, max(int((e0 - s0) * 12), 20))
            off = 0.5 + helix_amp * np.sin(2 * np.pi * (p - s0) / helix_period)
            # white underlay hides the backbone line behind the spring
            ax.plot([s0, e0], [0.5, 0.5], color="white", lw=coil_lw * 3, zorder=1,
                    solid_capstyle="butt")
            ax.plot(p, off, color=helix_color, lw=helix_lw, zorder=2, solid_capstyle="butt")
        else:
            length = e0 - s0
            ax.add_patch(FancyArrow(s0, 0.5, length, 0.0, width=strand_width,
                                    head_width=head_width,
                                    head_length=min(head_len, length * 0.6),
                                    length_includes_head=True, fc=strand_color, ec="none",
                                    zorder=2))


def _nested_levels(intervals):
    """Depth level per interval: shortest intervals first, each on the lowest level free of
    overlap. Nested sites thus sit closer to the SS track than the site enclosing them, so
    the outer site's drop lines do not cross the inner site's line."""
    placed = []  # (level, start, end)
    levels = [0] * len(intervals)
    for i in sorted(range(len(intervals)), key=lambda i: intervals[i][1] - intervals[i][0]):
        s, e = intervals[i]
        lvl = 0
        while any(l == lvl and s <= pe + 1 and ps <= e + 1 for l, ps, pe in placed):
            lvl += 1
        levels[i] = lvl
        placed.append((lvl, s, e))
    return levels


def draw_zinc(ax, sites, start, end, color="0.25", lw=0.6, d_base=0.3, d_step=0.25):
    """Zinc coordination drawn as trees hanging from the SS track.

    Each site is one line spanning its first to last ligand residue at depth
    `d_base` + level * `d_step` below the top of the track, with a drop line from every
    ligand residue to it. Overlapping sites get different levels.
    """
    _setup_track(ax, start, end)
    levels = _nested_levels([(min(s), max(s)) for s in sites])
    for site, lvl in zip(sites, levels):
        if not _in_range(min(site), max(site), start, end):
            continue
        y = 1 - (d_base + lvl * d_step)
        kw = dict(color=color, lw=lw, solid_capstyle="projecting")
        for r in site:
            ax.plot([r, r], [1, y], **kw)
        ax.plot([min(site), max(site)], [y, y], **kw)


def _heatmap_values(mat, wt, start, end):
    """Panel slice as (n_aa, n_pos) masked array; wild-type cells set to 0 (white)."""
    sub = mat[start - 1:end].copy()
    sub[wt[start - 1:end] & np.isnan(sub)] = 0.0
    return np.ma.masked_invalid(sub.T)


def draw_heatmap(ax, mat, wt, start, end, norm=None, aa_order=pdata.AA_ORDER, fontsize=4.5,
                 wt_marker_size=1.2):
    """Variant heatmap drawn as vector cells (crisp at any zoom / in PDF).

    mat, wt: (n_res, n_aa) arrays from variant_matrix (row 0 = residue 1), columns matching
    `aa_order`. Wild-type cells are white with a black dot; missing cells are grey.
    Amino-acid labels are shown on both ends of the map.
    Returns the QuadMesh so its values can be updated with set_array.
    """
    norm = norm or heatmap_norm()
    x_edges = np.arange(start - 0.5, end + 1)
    y_edges = np.arange(-0.5, len(aa_order))
    mesh = ax.pcolormesh(x_edges, y_edges, _heatmap_values(mat, wt, start, end),
                         cmap=heatmap_cmap(), norm=norm, edgecolors="face", linewidth=0.05)
    rr, cc = np.nonzero(wt[start - 1:end])
    ax.scatter(rr + start, cc, s=wt_marker_size, c="k", lw=0)
    ax.set_ylim(len(aa_order) - 0.5, -0.5)
    ax.set_yticks(range(len(aa_order)), aa_order, fontsize=fontsize)
    ax.set_xticks([])
    ax.tick_params(length=1.5, width=0.4, pad=1, left=True, right=True, labelleft=True,
                   labelright=True)
    for s in ax.spines.values():
        s.set_linewidth(0.4)
    return mesh


def draw_median_track(ax, med, start, end, ylim=None, color="0.3", fontsize=5,
                      tick_every=10):
    """Per-residue median bars (NaN = no bar).

    This track also carries the residue-number ticks (every `tick_every` residues).
    """
    pos = np.arange(start, end + 1)
    m = med[start - 1:end]
    ok = np.isfinite(m)
    ax.bar(pos[ok], m[ok], width=0.9, lw=0, color=color)
    ylim = ylim or np.nanmax(np.abs(med)) * 1.05
    ax.set_ylim(-ylim, ylim)
    ticks = [t for t in range(0, end + 1, tick_every) if t >= start]
    ax.set_xticks(ticks, [str(t) for t in ticks], fontsize=fontsize)
    ax.set_yticks([-0.2, 0, 0.2], ["-0.2", "0", "0.2"])
    ax.tick_params(labelsize=fontsize, length=1.5, width=0.4, pad=1)
    ax.axhline(0, color="0.4", lw=0.4)
    for s in ax.spines.values():
        s.set_linewidth(0.4)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# ---------------------------------------------------------------- full figure
def heatmap_figure(rd, annot, zinc_sites, n_panels=2, vlim=0.5, ss_kwargs=None,
                   width=A4_WIDTH_IN, dpi=100, pyplot=True):
    """Residual heatmap split into `n_panels` stacked panels with domain, SS, zinc and
    per-residue median tracks.

    rd: parkin_data.residual_data(...) for the residual to show.
    ss_kwargs: ss_segments options (helix_codes, strand_codes, min_helix, min_strand) and
    draw_ss styling (helix_period, helix_amp, ...).
    pyplot=False builds a bare Figure (not shown by the notebook), e.g. for saving.
    Returns (fig, layout); layout holds the panel axes and artists needed to update the
    figure in place (see update_heatmap_figure).
    """
    ss_kwargs = ss_kwargs or {}
    seg_keys = ("helix_codes", "strand_codes", "min_helix", "min_strand")
    ss = pdata.ss_segments(annot["SS"],
                           **{k: v for k, v in ss_kwargs.items() if k in seg_keys})
    ss_draw = {k: v for k, v in ss_kwargs.items() if k not in seg_keys}
    doms = pdata.domain_segments(annot["domain"])
    n_res = len(annot["resi"])
    chunk = int(np.ceil(n_res / n_panels))

    # row heights in inches: domain, SS, zinc, heatmap, gap, median (+ gap between panels)
    sizes = [0.13, 0.13, 0.22, len(pdata.AA_ORDER) * 0.06, 0.03, 0.3]
    gap = 0.35
    ratios = []
    for p in range(n_panels):
        ratios += sizes + ([gap] if p < n_panels - 1 else [])
    height = sum(ratios) + 0.9
    fig = (plt.figure if pyplot else Figure)(figsize=(width, height), dpi=dpi)
    gs = fig.add_gridspec(len(ratios), 1, height_ratios=ratios, hspace=0, left=0.055,
                          right=0.975, top=1 - 0.15 / height, bottom=0.75 / height)
    stride = len(sizes) + 1
    norm = heatmap_norm(vlim)
    panels = []
    for p in range(n_panels):
        start, end = p * chunk + 1, min((p + 1) * chunk, n_res)
        full_end = start + chunk - 1  # identical cell width in every panel
        ax_dom, ax_ss, ax_zn, ax_hm = (fig.add_subplot(gs[p * stride + k, 0]) for k in range(4))
        ax_med = fig.add_subplot(gs[p * stride + 5, 0])
        draw_domains(ax_dom, doms, start, full_end)
        draw_ss(ax_ss, ss, start, full_end, **ss_draw)
        draw_zinc(ax_zn, zinc_sites, start, full_end)
        mesh = draw_heatmap(ax_hm, rd["mat"], rd["wt"], start, end, norm)
        ax_hm.set_xlim(start - 0.5, full_end + 0.5)
        panels.append(dict(start=start, end=end, full_end=full_end, ax_hm=ax_hm,
                           ax_med=ax_med, mesh=mesh))

    cax = fig.add_axes([0.3, 0.3 / height, 0.4, 0.07 / height])
    cb = fig.colorbar(panels[0]["mesh"], cax=cax, orientation="horizontal", extend="both")
    cb.ax.tick_params(labelsize=5, length=1.5, width=0.4, pad=1)
    cb.outline.set_linewidth(0.4)
    handles = [Rectangle((0, 0), 1, 1, fc=MISSING_COLOR, ec="none"),
               Line2D([], [], marker="o", ls="", color="k", ms=2.5)]
    fig.legend(handles, ["Missing", "Wild type"], fontsize=5, frameon=False, ncol=2,
               loc="lower left", bbox_to_anchor=(0.72, cax.get_position().y0 - 0.006),
               handlelength=1, handletextpad=0.4, columnspacing=1)
    layout = dict(panels=panels, colorbar=cb, norm=norm)
    update_heatmap_figure(fig, layout, rd)
    return fig, layout


def update_heatmap_figure(fig, layout, rd):
    """Redraw the residual-dependent parts (cells, median track, colour-bar label) for
    rd = parkin_data.residual_data(...)."""
    med = rd["med"]
    med_lim = np.nanmax(np.abs(med)) * 1.05
    for i, p in enumerate(layout["panels"]):
        p["mesh"].set_array(_heatmap_values(rd["mat"], rd["wt"], p["start"], p["end"]))
        ax = p["ax_med"]
        ax.clear()
        draw_median_track(ax, med, p["start"], p["end"], ylim=med_lim)
        ax.set_xlim(p["start"] - 0.5, p["full_end"] + 0.5)
        if i == 0:
            ax.set_ylabel("Median", fontsize=5, labelpad=1)
    layout["colorbar"].set_label(f"{rd['residual']} residual", fontsize=6, labelpad=1)
    fig.canvas.draw_idle()
