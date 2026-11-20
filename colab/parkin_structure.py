"""3D structure view (py3Dmol), coloured per residue (by domain or median residual).

The view is returned as an <iframe srcdoc> HTML string for an ipywidgets.HTML widget
(same approach as the HPRT notebook): a nested document always runs its own scripts.
- The 3Dmol.js library is embedded from js/3Dmol-min.js, so no internet connection is
  needed (falls back to the CDN if that file is missing).
- The iframe is rebuilt on every update, so the camera orientation is stored in the
  browser's localStorage and restored: rotating and then selecting a residue keeps the view.
- Hovering a residue shows a label (3Dmol hover callbacks, inside the viewer only).
"""

import html
import json
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import py3Dmol
from matplotlib.colors import to_hex

import parkin_plots as pp

NO_DATA_COLOR = "#bdbdbd"      # median NaN (too few variants)
LINKER_COLOR = "#f2f2f2"       # residues outside annotated domains

# Local copy of the 3Dmol.js library (same version py3Dmol requests from the CDN)
LOCAL_3DMOL = Path(__file__).resolve().parent / "js" / "3Dmol-min.js"
_CDN_LOADER = re.compile(r"loadScriptAsync\('[^']*3Dmol-min\.js'\)")


@lru_cache(maxsize=1)
def _local_3dmol_js():
    return LOCAL_3DMOL.read_text() if LOCAL_3DMOL.exists() else None


def _use_local_3dmol(page):
    """Swap py3Dmol's CDN loader for the embedded local library (if present)."""
    lib = _local_3dmol_js()
    if lib is None or not _CDN_LOADER.search(page):
        return page  # fall back to loading from the CDN
    return f"<script>{lib}</script>\n" + _CDN_LOADER.sub("Promise.resolve()", page)


def load_pdb(path):
    with open(path) as fh:
        return fh.read()


# ---------------------------------------------------------------- residue colours
def domain_colors(domain, colors=pp.DOMAIN_COLORS, linker_color=LINKER_COLOR):
    """Per-residue hex colours {resi: '#rrggbb'} from the per-residue domain labels."""
    return {i + 1: to_hex(colors.get(d, linker_color)) if d else linker_color
            for i, d in enumerate(domain)}


def median_colors(med, vlim=0.3):
    """Per-residue hex colours {resi: '#rrggbb'} on the heatmap colour map (+-vlim);
    NaN medians are grey."""
    cmap, norm = pp.heatmap_cmap(), pp.heatmap_norm(vlim)
    return {i + 1: to_hex(cmap(norm(m))) if np.isfinite(m) else NO_DATA_COLOR
            for i, m in enumerate(med)}


# ------------------------------------------------------------------- javascript
# restores the stored camera and keeps storing it whenever the view changes
_VIEW_JS = """
<script>
(function() {
  var key = %(key)s, name = %(viewer)s, tries = 0;
  function hook() {
    var v = window[name];
    if (!v) { if (tries++ < 200) setTimeout(hook, 50); return; }
    try {
      var saved = window.localStorage.getItem(key);
      if (saved) { v.setView(JSON.parse(saved)); v.render(); }
      v.setViewChangeCallback(function(view) {
        try { window.localStorage.setItem(key, JSON.stringify(view)); } catch (e) {}
      });
    } catch (e) {}
  }
  hook();
})();
</script>
"""


def _hover_js(hover_text):
    """3Dmol hover/unhover callbacks showing one label for the residue under the mouse."""
    lookup = json.dumps({str(k): v for k, v in hover_text.items()})
    hover = (
        "function(atom, viewer) {"
        f"  var txt = {lookup}[String(atom.resi)] || (atom.resn + ' ' + atom.resi);"
        "   if (viewer._hoverLabel) viewer.removeLabel(viewer._hoverLabel);"
        "   viewer._hoverLabel = viewer.addLabel(txt, {"
        "       position: {x: atom.x, y: atom.y, z: atom.z},"
        "       backgroundColor: '#222', fontColor: 'white', fontSize: 12,"
        "       backgroundOpacity: 0.85, inFront: true, alignment: 'bottomCenter'});"
        "   viewer.render();"
        "}")
    unhover = (
        "function(atom, viewer) {"
        "   if (viewer._hoverLabel) { viewer.removeLabel(viewer._hoverLabel);"
        "       viewer._hoverLabel = null; viewer.render(); }"
        "}")
    return hover, unhover


# ------------------------------------------------------------------------- view
def structure_html(pdb_str, colors, selected=None, selected_label=None, hover_text=None,
                   width=420, height=460, select_color=pp.SELECT_COLOR,
                   view_key="parkin_structure_view"):
    """py3Dmol cartoon coloured per residue, returned as an <iframe srcdoc=...> string.

    colors: {resi: hex} (domain_colors or median_colors).
    selected: residue number or None; drawn as side-chain sticks in its backbone colour
    plus a CA sphere and a label (`selected_label`) in `select_color`.
    hover_text: {resi: label} shown when hovering a residue (default: name and number).
    """
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_str, "pdb")
    view.setStyle({}, {"cartoon": {"color": NO_DATA_COLOR}})
    by_color = {}
    for r, c in colors.items():
        by_color.setdefault(c, []).append(r)
    for c, res in by_color.items():
        view.setStyle({"resi": res}, {"cartoon": {"color": c}})
    if selected:
        view.addStyle({"resi": [selected], "not": {"atom": ["N", "C", "O"]}},
                      {"stick": {"color": colors.get(selected, NO_DATA_COLOR), "radius": 0.3}})
        # CA sphere so the position is visible also for glycine (no side chain)
        view.addStyle({"resi": [selected], "atom": "CA"},
                      {"sphere": {"color": select_color, "radius": 0.7}})
        view.addLabel(selected_label or str(selected),
                      {"backgroundColor": select_color, "fontColor": "white", "fontSize": 12,
                       "backgroundOpacity": 0.85},
                      {"resi": selected, "atom": "CA"})
    view.setHoverable({}, True, *_hover_js(hover_text or {}))
    view.setHoverDuration(100)
    view.zoomTo()
    page = _use_local_3dmol(view._make_html()) + _VIEW_JS % {
        "key": json.dumps(view_key), "viewer": json.dumps(f"viewer_{view.uniqueid}")}
    return (f'<iframe srcdoc="{html.escape(page)}" '
            f'style="width:{width + 10}px;height:{height + 10}px;border:0"></iframe>')
