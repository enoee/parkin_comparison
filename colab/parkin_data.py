"""Data loading and derived quantities for the Parkin MAVE comparison.

Only numpy + stdlib (no pandas in the environment).

Data flow
---------
load_data(csv)                  -> dict of column arrays (one entry per variant)
residue_table(data)             -> per-residue annotations (wt, SS, domain, zf)
domain_segments(domain)         -> [(start, end, name), ...]
ss_segments(ss, ...)            -> [(start, end, 'helix'|'strand'), ...]
zinc_sites_from_column(zf)      -> [tuple of ligand residue numbers, ...]
residual_data(data, residual)   -> everything plots need for one residual: per-variant
                                   residuals and scores, the (residue x aa) matrix, the
                                   wild-type mask and per-residue medians
"""

import csv
import numpy as np

AA_ORDER = list("CDEKRHNQAGSTVMLIFYWP") + ["*"]   # heatmap rows
MIN_VARIANTS = 5    # per-residue medians need at least this many measured variants

# dropdown label -> (kind, method code); the residual column is f"res_{kind}_{method}"
RESIDUALS = {
    "Activity: Deming": ("act", "dem"),
    "Activity: PC spline": ("act", "pcs"),
    "Activity: PC LOWESS": ("act", "pcl"),
    "Abundance: Deming": ("abun", "dem"),
    "Abundance: PC spline": ("abun", "pcs"),
    "Abundance: PC LOWESS": ("abun", "pcl"),
}
DEFAULT_RESIDUAL = "Activity: PC spline"
SCORE_COLUMNS = {"act": ("activity", "greg_act"),   # kind -> (MAVE 1 = x-axis, MAVE 2)
                 "abun": ("abun", "greg_abun")}
KIND_LABEL = {"act": "activity", "abun": "abundance"}
# names shown in plots/hover text for the score columns (what each assay measures)
COLUMN_LABELS = {"activity": "mitophagy", "greg_act": "pSer65Ub"}


def column_label(col):
    return COLUMN_LABELS.get(col, col)


def residual_column(residual):
    """Dropdown label (key of RESIDUALS) -> csv column, e.g. 'res_act_pcs'."""
    kind, method = RESIDUALS[residual]
    return f"res_{kind}_{method}"


# --------------------------------------------------------------------------- loading
def _to_float(values):
    out = np.full(len(values), np.nan)
    for i, v in enumerate(values):
        try:
            out[i] = float(v)
        except (TypeError, ValueError):
            pass
    return out


def load_data(path):
    """Read the comparison csv into a dict of numpy arrays (variants only, WT row dropped).

    Columns that parse as numbers become float arrays (NA -> NaN); others stay strings.
    """
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh) if r["resi"] not in ("", "NA")]
    data = {}
    for col in rows[0]:
        raw = [r[col] for r in rows]
        num = _to_float(raw)
        non_missing = [v for v in raw if v not in ("", "NA")]
        is_num = non_missing and np.isfinite(num).sum() == len(non_missing)
        data[col] = num if is_num else np.array(raw, dtype=object)
    data["resi"] = data["resi"].astype(int)
    return data


def residue_table(data):
    """Per-residue annotations (index 0 = residue 1).

    Keys: wt, SS, domain, zf, resi. The domain annotation is read from the csv column
    'domain' or, if absent, 'pfam'.
    """
    n = int(data["resi"].max())
    source = {"wt": "wt", "SS": "SS", "zf": "zf",
              "domain": "domain" if "domain" in data else "pfam"}
    tab = {k: np.array([""] * n, dtype=object) for k in source}
    for k, col in source.items():
        for r, v in zip(data["resi"], data[col]):
            if tab[k][r - 1] == "":
                tab[k][r - 1] = "" if v in ("NA", None) else v
    tab["resi"] = np.arange(1, n + 1)
    return tab


# ----------------------------------------------------------------------- annotation
def _runs(labels):
    """Run-length encode a per-residue label array -> [(start, end, label)] (1-based)."""
    runs, start = [], 0
    for i in range(1, len(labels) + 1):
        if i == len(labels) or labels[i] != labels[start]:
            runs.append((start + 1, i, labels[start]))
            start = i
    return runs


def domain_segments(domain):
    return [(s, e, d) for s, e, d in _runs(list(domain)) if d not in ("", "NA")]


def ss_segments(ss, helix_codes="HGI", strand_codes="E", min_helix=3, min_strand=2):
    """Collapse DSSP codes into helix / strand segments.

    helix_codes, strand_codes: which DSSP letters count as helix / strand
    min_helix, min_strand: minimum segment length to be drawn
    """
    cls = ["helix" if c in helix_codes else "strand" if c in strand_codes else "coil" for c in ss]
    min_len = {"helix": min_helix, "strand": min_strand}
    return [(s, e, c) for s, e, c in _runs(cls) if c != "coil" and e - s + 1 >= min_len[c]]


def zinc_sites_from_column(zf):
    """Group residues by zinc-finger label (index 0 = residue 1),
    e.g. zf2 -> (238, 241, 260, 263)."""
    sites = {}
    for r, z in enumerate(zf, start=1):
        if z not in ("", "NA"):
            sites.setdefault(z, []).append(r)
    return [tuple(sorted(v)) for _, v in sorted(sites.items(), key=lambda kv: min(kv[1]))]


# ------------------------------------------------------------------------ residuals
def variant_matrix(data, column, n_res=None, aa_order=AA_ORDER):
    """Score matrix of shape (n_res, len(aa_order)) and a boolean wild-type mask."""
    n_res = int(data["resi"].max()) if n_res is None else n_res
    col_idx = {a: j for j, a in enumerate(aa_order)}
    mat = np.full((n_res, len(aa_order)), np.nan)
    wt = np.zeros_like(mat, dtype=bool)
    for r, w, m, v in zip(data["resi"], data["wt"], data["mut_aa"], data[column]):
        if m in col_idx:
            mat[r - 1, col_idx[m]] = v
        if w in col_idx:
            wt[r - 1, col_idx[w]] = True
    return mat, wt


def per_residue_median(data, column, min_n=MIN_VARIANTS, n_res=None):
    """Median over missense variants per position; NaN when fewer than `min_n`."""
    n_res = int(data["resi"].max()) if n_res is None else n_res
    vals = [[] for _ in range(n_res)]
    for r, w, m, v in zip(data["resi"], data["wt"], data["mut_aa"], data[column]):
        if np.isfinite(v) and m not in ("*", "=", w, ""):
            vals[r - 1].append(v)
    return np.array([np.median(v) if len(v) >= min_n else np.nan for v in vals])


def residual_data(data, residual, min_n=MIN_VARIANTS):
    """Everything the figures need for one residual (key of RESIDUALS), computed once.

    Returns a dict with: residual, kind, column, res (per variant), x / y (MAVE 1 / 2
    scores per variant), xcol / ycol, mat + wt (from variant_matrix), med (per residue),
    min_n.
    """
    kind, _ = RESIDUALS[residual]
    column = residual_column(residual)
    xcol, ycol = SCORE_COLUMNS[kind]
    mat, wt = variant_matrix(data, column)
    return dict(residual=residual, kind=kind, column=column, res=data[column],
                x=data[xcol], y=data[ycol], xcol=xcol, ycol=ycol, mat=mat, wt=wt,
                med=per_residue_median(data, column, min_n), min_n=min_n)
