"""Linked scatter + heatmap + structure explorer (needs ipympl: `%matplotlib widget`).

Layout: controls, a scatter of the two MAVE scores and (optionally) the 3D structure on top;
the residual heatmap below. Points, cells and residues are linked in two ways:
* hover (black): moving over a scatter point outlines its heatmap cell and vice versa;
  the variant and its scores are shown under the controls.
* selected residue (green): set by clicking a hovered point or cell (clicking it again
  clears it) or by typing a residue number. It is boxed in the heatmap, its variants are
  ringed and labelled in the scatter and it is shown as sticks on the structure.
The structure is coloured by domain or by median residual (toggle): it starts with
domains, switches to residuals when a residue is selected, and "Reset markings" clears all
marks and returns to domains.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgba

import parkin_data as pdata
import parkin_plots as pp
import parkin_structure as ps

NO_RESIDUAL_COLOR = "0.6"   # scatter points without a residual


def _style_canvas(fig):
    fig.canvas.header_visible = False
    fig.canvas.footer_visible = False
    fig.canvas.toolbar_visible = "fade-in-fade-out"
    fig.canvas.resizable = False
    _keep_hidpi(fig)


def _keep_hidpi(fig):
    """Work around an ipympl race: on a fresh 'Run all' a figure can end up with the
    screen's device-pixel ratio set but its dpi not scaled, so it shows at half size and
    resolution. After each screen draw, restore dpi = original dpi * pixel ratio
    (uses matplotlib's private _original_dpi / _set_dpi, as matplotlib itself does)."""
    canvas = fig.canvas

    def check(event):
        if event.canvas is not canvas or getattr(canvas, "_is_saving", False):
            return
        expected = fig._original_dpi * canvas.device_pixel_ratio
        if fig.dpi != expected:
            fig._set_dpi(expected, forward=False)
            canvas.manager.resize(*fig.bbox.size, forward=False)  # tell the front-end
            canvas.draw_idle()

    canvas.mpl_connect("draw_event", check)


class LinkedExplorer:
    """Interactive view; display `LinkedExplorer(...).widget` in a notebook cell."""

    def __init__(self, data, annot, zinc_sites, residual=pdata.DEFAULT_RESIDUAL,
                 ss_kwargs=None, pdb_path=None, struct_vlim=0.3, dpi=140, pick_radius_px=6,
                 figure_dir="figures"):
        """data / annot / zinc_sites: from parkin_data (load_data, residue_table,
        zinc_sites_from_column). ss_kwargs: secondary-structure options for the heatmap.
        pdb_path: structure to show (None = no structure); struct_vlim: colour-scale limit
        (+-) for median residuals on the structure. figure_dir: where "Save" writes PDFs."""
        self.data, self.annot, self.zinc_sites = data, annot, zinc_sites
        self.ss_kwargs, self.dpi, self.pick_radius = ss_kwargs, dpi, pick_radius_px
        self.figure_dir, self.struct_vlim = figure_dir, struct_vlim
        self.pdb_str = ps.load_pdb(pdb_path) if pdb_path else None
        self.n_res = len(annot["resi"])
        self._index_variants()
        self._hover = None             # variant index under the mouse
        self._suspend_struct = False   # batch several changes into one structure rebuild

        self._build_controls(residual)
        self._build_structure_widgets()
        self.rd = pdata.residual_data(data, residual)
        with plt.ioff():  # figures are shown inside the widget layout, not on creation
            self.fig_hm, self.layout = pp.heatmap_figure(self.rd, annot, zinc_sites,
                                                         ss_kwargs=ss_kwargs, dpi=dpi)
            self._build_scatter()
        _style_canvas(self.fig_hm)
        self._build_heatmap_overlays()
        self._connect_events()
        self._update_residual()

        panels = [self.controls, self.fig_sc.canvas]
        if self.pdb_str is not None:
            panels.append(widgets.VBox([self.w_color_mode, self.w_struct,
                                        self.w_struct_caption]))
        top = widgets.HBox(panels, layout=widgets.Layout(align_items="flex-start",
                                                         flex_flow="row wrap"))
        self.widget = widgets.VBox([top, self.fig_hm.canvas])

    # ------------------------------------------------------------------ data
    def _index_variants(self):
        """Per-variant arrays plus a (residue, aa) -> variant index lookup."""
        d = self.data
        aa_col = {a: j for j, a in enumerate(pdata.AA_ORDER)}
        self.v_resi, self.v_wt, self.v_mut = d["resi"], d["wt"], d["mut_aa"]
        self.v_col = np.array([aa_col.get(m, -1) for m in self.v_mut])
        self.cell_index = np.full((self.n_res, len(pdata.AA_ORDER)), -1)
        for i, (r, j) in enumerate(zip(self.v_resi, self.v_col)):
            if j >= 0:
                self.cell_index[r - 1, j] = i

    def _variant_name(self, i):
        return f"{self.v_wt[i]}{self.v_resi[i]}{self.v_mut[i]}"

    def _residue_name(self, r):
        return f"{self.annot['wt'][r - 1]}{r}"

    # ----------------------------------------------------------------- build
    def _build_controls(self, residual):
        style = dict(description_width="70px")
        self.w_residual = widgets.Dropdown(options=list(pdata.RESIDUALS), value=residual,
                                           description="Residual", style=style,
                                           layout=widgets.Layout(width="290px"))
        self.w_residue = widgets.BoundedIntText(value=0, min=0, max=self.n_res,
                                                description="Residue", style=style,
                                                layout=widgets.Layout(width="160px"))
        self.w_reset = widgets.Button(description="Reset markings", icon="eraser",
                                      tooltip="Clear the selected residue and hover marks",
                                      layout=widgets.Layout(width="200px"))
        self.w_save = widgets.Button(description="Save heatmap PDF", icon="download",
                                     layout=widgets.Layout(width="200px"))
        self.w_info = widgets.HTML(layout=widgets.Layout(width="290px", min_height="90px"))
        self.w_status = widgets.HTML(layout=widgets.Layout(width="290px"))
        help_txt = widgets.HTML(
            "<small>Residue 0 = no selection. "
            f"<span style='color:{pp.SELECT_COLOR}'>Green</span> = selected residue, "
            "black = variant under the mouse. Click a point or cell to select its residue; "
            "click the same residue again to clear. Hover the structure to see residue "
            "labels.</small>", layout=widgets.Layout(width="290px"))

        self.w_residual.observe(lambda c: self._update_residual(), names="value")
        self.w_residue.observe(lambda c: self._update_selection(), names="value")
        self.w_reset.on_click(self._reset_markings)
        self.w_save.on_click(self._save_pdf)
        self.controls = widgets.VBox([self.w_residual, self.w_residue, help_txt, self.w_info,
                                      self.w_reset, self.w_save, self.w_status])

    def _build_structure_widgets(self):
        self.w_struct = widgets.HTML()
        self.w_struct_caption = widgets.HTML(layout=widgets.Layout(width="420px"))
        self.w_color_mode = widgets.ToggleButtons(
            options=["Domains", "Residuals"], value="Domains", description="Colour",
            style=dict(description_width="45px", button_width="90px"))
        self.w_color_mode.observe(lambda c: self._update_structure(), names="value")

    def _build_scatter(self):
        self.fig_sc, ax = plt.subplots(figsize=(3.6, 3.3), dpi=self.dpi)
        self.ax_sc = ax
        _style_canvas(self.fig_sc)
        # thin grey edge keeps near-zero (white) residuals visible on the white background
        self.sc = ax.scatter([], [], s=6, edgecolors="0.6", linewidths=0.15)
        self.sel_sc = ax.scatter([], [], s=28, facecolors="none", edgecolors=pp.SELECT_COLOR,
                                 lw=1.0, zorder=3)
        self.sel_labels = []
        (self.hover_pt,) = ax.plot([], [], "o", ms=7, mfc="none", mec=pp.HOVER_COLOR,
                                   mew=1.3, zorder=4)
        ax.tick_params(labelsize=6, length=2, width=0.5)
        for s in ax.spines.values():
            s.set_linewidth(0.5)
        self.fig_sc.subplots_adjust(left=0.16, right=0.97, top=0.95, bottom=0.14)

    def _build_heatmap_overlays(self):
        """Per panel: a hover outline (one cell) and a selection box (one column)."""
        self.hover_rects, self.sel_rects = [], []
        for p in self.layout["panels"]:
            hr = Rectangle((0, 0), 1, 1, fill=False, ec=pp.HOVER_COLOR, lw=1.2,
                           visible=False, zorder=5)
            sr = Rectangle((0, -0.5), 1, len(pdata.AA_ORDER), fill=False,
                           ec=pp.SELECT_COLOR, lw=1.2, visible=False, zorder=4)
            p["ax_hm"].add_patch(hr)
            p["ax_hm"].add_patch(sr)
            self.hover_rects.append(hr)
            self.sel_rects.append(sr)

    def _connect_events(self):
        self.fig_hm.canvas.mpl_connect("motion_notify_event", self._on_heatmap_motion)
        self.fig_sc.canvas.mpl_connect("motion_notify_event", self._on_scatter_motion)
        for fig in (self.fig_hm, self.fig_sc):
            fig.canvas.mpl_connect("figure_leave_event", lambda e: self._set_hover(None))
            fig.canvas.mpl_connect("axes_leave_event", lambda e: self._set_hover(None))
            fig.canvas.mpl_connect("button_press_event", self._on_click)

    def _redraw(self):
        self.fig_sc.canvas.draw_idle()
        self.fig_hm.canvas.draw_idle()

    # --------------------------------------------------------------- updates
    def _update_residual(self):
        """New residual: recompute its data once, then update scatter, heatmap,
        selection, hover info and structure."""
        rd = self.rd = pdata.residual_data(self.data, self.w_residual.value)
        ok = np.isfinite(rd["x"]) & np.isfinite(rd["y"])
        order = np.nonzero(ok)[0]   # strongest residuals drawn last (on top)
        order = order[np.argsort(np.nan_to_num(np.abs(rd["res"][order]), nan=-1))]
        self.sc_index = order       # scatter point k -> variant index
        self.sc_xy = np.column_stack([rd["x"][order], rd["y"][order]])

        r = rd["res"][order]
        colors = pp.heatmap_cmap()(self.layout["norm"](np.nan_to_num(r)))
        colors[np.isnan(r)] = to_rgba(NO_RESIDUAL_COLOR)
        self.sc.set_offsets(self.sc_xy)
        self.sc.set_facecolors(colors)
        lo, hi = self.sc_xy.min(0), self.sc_xy.max(0)
        pad = 0.04 * (hi - lo)
        self.ax_sc.set_xlim(lo[0] - pad[0], hi[0] + pad[0])
        self.ax_sc.set_ylim(lo[1] - pad[1], hi[1] + pad[1])
        label = pdata.KIND_LABEL[rd["kind"]]
        self.ax_sc.set_xlabel(f"MAVE 1 {label} ({pdata.column_label(rd['xcol'])})", fontsize=7)
        self.ax_sc.set_ylabel(f"MAVE 2 {label} ({pdata.column_label(rd['ycol'])})", fontsize=7)

        pp.update_heatmap_figure(self.fig_hm, self.layout, rd)
        self._update_selection(redraw=False)
        self._set_hover(self._hover, force=True)

    def _update_selection(self, redraw=True):
        """Mark the selected residue (0 = none) in the scatter, heatmap and structure."""
        r = self.w_residue.value
        for t in self.sel_labels:
            t.remove()
        self.sel_labels = []
        for p, rect in zip(self.layout["panels"], self.sel_rects):
            rect.set_visible(bool(r) and p["start"] <= r <= p["end"])
            rect.set_x(r - 0.5)
        if r:
            k = np.nonzero(self.v_resi[self.sc_index] == r)[0]
            self.sel_sc.set_offsets(self.sc_xy[k] if len(k) else np.empty((0, 2)))
            for (px, py), i in zip(self.sc_xy[k], self.sc_index[k]):
                self.sel_labels.append(self.ax_sc.annotate(
                    self.v_mut[i], (px, py), xytext=(2, 2), textcoords="offset points",
                    fontsize=5, color=pp.SELECT_COLOR, zorder=5))
            med = self.rd["med"][r - 1]
            title = f"{self._residue_name(r)}: {len(k)} variants"
            if np.isfinite(med):
                title += f", median residual {med:.2f}"
            self.ax_sc.set_title(title, fontsize=6, color=pp.SELECT_COLOR)
            if self.w_color_mode.value != "Residuals":
                self.w_color_mode.value = "Residuals"   # its observer rebuilds the structure
            else:
                self._update_structure()
        else:
            self.sel_sc.set_offsets(np.empty((0, 2)))
            self.ax_sc.set_title("")
            self._update_structure()
        if redraw:
            self._redraw()

    def _reset_markings(self, _button=None):
        """Clear the selection and hover marks; structure back to domain colours
        (rebuilt once)."""
        self._suspend_struct = True
        try:
            self.w_residue.value = 0   # triggers _update_selection
            self.w_color_mode.value = "Domains"
        finally:
            self._suspend_struct = False
        self._update_structure()
        self._set_hover(None, force=True)

    # ------------------------------------------------------------- structure
    def _update_structure(self):
        """Rebuild the 3D view and its caption for the current colour mode and selection."""
        if self.pdb_str is None or self._suspend_struct:
            return
        selected = ("Selected residue: side-chain sticks in its backbone colour, marked by "
                    f"a <span style='color:{pp.SELECT_COLOR}'>green</span> sphere and label.")
        if self.w_color_mode.value == "Domains":
            colors = ps.domain_colors(self.annot["domain"])
            key = " ".join(self._domain_swatch(d)
                           for d in dict.fromkeys(d for d in self.annot["domain"] if d))
            self.w_struct_caption.value = f"<small>Domains: {key}<br>{selected}</small>"
        else:
            colors = ps.median_colors(self.rd["med"], self.struct_vlim)
            self.w_struct_caption.value = (
                f"<small>{self.rd['residual']} median residual (colour scale "
                f"±{self.struct_vlim}); grey: fewer than {self.rd['min_n']} variants."
                f"<br>{selected}</small>")
        r = self.w_residue.value or None
        self.w_struct.value = ps.structure_html(
            self.pdb_str, colors, selected=r,
            selected_label=self._residue_name(r) if r else None,
            hover_text=self._structure_hover_text())

    @staticmethod
    def _domain_swatch(d):
        c = pp.DOMAIN_COLORS.get(d, ps.LINKER_COLOR)
        return (f"<span style='background:{c};color:{pp.text_color_on(c)};"
                f"padding:0 4px;border:1px solid #888'>{d}</span>")

    def _structure_hover_text(self):
        """Per-residue hover label, e.g. 'W403 REP: -0.01' (median residual, or NA)."""
        text = {}
        for r in range(1, self.n_res + 1):
            m = self.rd["med"][r - 1]
            dom = self.annot["domain"][r - 1] or "linker"
            text[r] = f"{self._residue_name(r)} {dom}: {f'{m:+.2f}' if np.isfinite(m) else 'NA'}"
        return text

    # ----------------------------------------------------------------- mouse
    def _on_click(self, event):
        """Left click on a hovered point/cell selects its residue (again = clear).
        Ignored while the toolbar's pan/zoom mode is active."""
        toolbar = getattr(event.canvas, "toolbar", None)
        if event.button != 1 or self._hover is None or (toolbar is not None and toolbar.mode):
            return
        r = int(self.v_resi[self._hover])
        self.w_residue.value = 0 if self.w_residue.value == r else r

    def _on_scatter_motion(self, event):
        if event.inaxes is not self.ax_sc or not len(self.sc_xy):
            return self._set_hover(None)
        pix = self.ax_sc.transData.transform(self.sc_xy)
        d2 = ((pix - (event.x, event.y)) ** 2).sum(1)
        k = int(np.argmin(d2))
        dpr = getattr(self.fig_sc.canvas, "device_pixel_ratio", 1)
        hit = d2[k] <= (self.pick_radius * dpr) ** 2
        self._set_hover(int(self.sc_index[k]) if hit else None)

    def _on_heatmap_motion(self, event):
        for p in self.layout["panels"]:
            if event.inaxes is p["ax_hm"] and event.xdata is not None:
                r, j = int(round(event.xdata)), int(round(event.ydata))
                if p["start"] <= r <= p["end"] and 0 <= j < len(pdata.AA_ORDER):
                    i = self.cell_index[r - 1, j]
                    return self._set_hover(int(i) if i >= 0 else None)
        self._set_hover(None)

    def _set_hover(self, i, force=False):
        """Show variant i (or nothing) as hover marks in both figures and in the info box."""
        if i == self._hover and not force:
            return
        self._hover = i
        for rect in self.hover_rects:
            rect.set_visible(False)
        self.hover_pt.set_data([], [])
        self.w_info.value = ""
        if i is not None:
            rd, r, j = self.rd, self.v_resi[i], self.v_col[i]
            for p, rect in zip(self.layout["panels"], self.hover_rects):
                if p["start"] <= r <= p["end"] and j >= 0:
                    rect.set_xy((r - 0.5, j - 0.5))
                    rect.set_visible(True)
            if np.isfinite(rd["x"][i]) and np.isfinite(rd["y"][i]):
                self.hover_pt.set_data([rd["x"][i]], [rd["y"][i]])
            fmt = lambda v: "NA" if not np.isfinite(v) else f"{v:.3f}"
            self.w_info.value = (
                f"<b>{self._variant_name(i)}</b> ({self.annot['domain'][r - 1] or 'linker'})<br>"
                f"{pdata.column_label(rd['xcol'])}: {fmt(rd['x'][i])}<br>"
                f"{pdata.column_label(rd['ycol'])}: {fmt(rd['y'][i])}<br>"
                f"residual: {fmt(rd['res'][i])}")
        self._redraw()

    # ------------------------------------------------------------------ save
    def _save_pdf(self, _button):
        fig, _ = pp.heatmap_figure(self.rd, self.annot, self.zinc_sites,
                                   ss_kwargs=self.ss_kwargs, pyplot=False)
        os.makedirs(self.figure_dir, exist_ok=True)
        kind, method = pdata.RESIDUALS[self.rd["residual"]]
        path = os.path.join(self.figure_dir, f"heatmap_{kind}_{method}.pdf")
        fig.savefig(path)
        self.w_status.value = f"<small>Saved {path}</small>"
