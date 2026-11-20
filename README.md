# Parkin: comparison of two MAVE experiments

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/<github-user>/<repo>/blob/main/parkin_mave_comparison_colab.ipynb)

Interactive comparison of two sets of multiplexed assays of variant effects (MAVE) on the
E3 ubiquitin ligase Parkin (PRKN, UniProt O60260): activity (mitophagy vs. pSer65Ub) and
abundance. The relation between the experiments is fitted with Deming regression and principal
curves (spline and LOWESS smoothing), and the per-variant residuals are shown as

- an A4-width heatmap along the sequence with domain, secondary-structure and zinc-site tracks
  and per-residue median residuals (saved as PDF);
- an interactive explorer: a scatter of the two experiments linked to the heatmap by hover and
  click, and the AlphaFold structure coloured by domain or median residual.

## Run on Google Colab

Click the badge above. Run the setup cell first: it installs `ipympl` and `py3Dmol`, enables
Colab's widget manager and downloads only the files the notebook needs from this repository
(the four `parkin_*.py` modules, `comparison_res.csv` and the AlphaFold model). If the
interactive figures show only as text, use *Runtime → Restart session* and run all cells
again.

## Run locally

Requires Python with numpy, matplotlib, ipywidgets, ipympl, py3Dmol and JupyterLab. Open
`parkin_mave_comparison.ipynb` and run all cells. The structure viewer loads 3Dmol.js from the
web; to use it offline, save `3Dmol-min.js` (version 2.5.4) as `js/3Dmol-min.js` in the
project folder and it is used instead. Quick checks of the code:

```bash
python tests/test_explorer.py
```

## Files

| File | Content |
|---|---|
| `parkin_mave_comparison.ipynb` | notebook for local use |
| `parkin_mave_comparison_colab.ipynb` | the same notebook with a Colab setup section |
| `parkin_data.py` | data loading, annotation and per-residual quantities |
| `parkin_plots.py` | colours and matplotlib figures (heatmap and tracks) |
| `parkin_interactive.py` | linked scatter / heatmap / structure explorer |
| `parkin_structure.py` | py3Dmol structure view |
| `comparison_res.csv` | variant scores, residuals and per-residue annotation |
| `AF-O60260-F1-model_v2.pdb` | AlphaFold model of Parkin |
| `LICENSE` | licence terms for this repository (beta, all rights reserved) |
| `tests/test_explorer.py` | headless checks of data, figures and explorer |

## Licence and credits

This is a beta version shared for viewing and testing: all rights reserved, see `LICENSE`.
Please do not reuse the code, data or figures without permission.

Third-party material keeps its own licence:
- **3Dmol.js** (structure viewer, via py3Dmol; loaded at run time, not part of this
  repository), BSD-3-Clause. Rego N & Koes D (2015) 3Dmol.js: molecular visualization with WebGL.
  *Bioinformatics* 31:1322–1324. https://doi.org/10.1093/bioinformatics/btu829
- **AlphaFold model** of Parkin (`AF-O60260-F1-model_v2.pdb`, AlphaFold Protein Structure
  Database, https://alphafold.ebi.ac.uk/entry/O60260), CC BY 4.0
  (https://creativecommons.org/licenses/by/4.0/).
  Jumper J *et al.* (2021) Highly accurate protein structure prediction with AlphaFold.
  *Nature* 596:583–589. https://doi.org/10.1038/s41586-021-03819-2
  Varadi M *et al.* (2022) AlphaFold Protein Structure Database. *Nucleic Acids Res*
  50:D439–D444. https://doi.org/10.1093/nar/gkab1061
