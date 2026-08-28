# Ritchie Lab SpatialBench Spleen Portal

**Live:** https://ritchie-spleen-portal.vercel.app (access-gated)

Interactive spatial-transcriptomics portal for **GEO [GSE254652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254652)**
(Du, **Ritchie** et al., *Genome Biology* 2025 — *Benchmarking spatial transcriptomics
technologies with the multi-sample SpatialBenchVisium dataset*), built in the same
**Explorer** chrome as the Aguado aortic-valve and May mouse-prostate portals.

The link opens on **Live** — a cinematic walkthrough of the analysis with a streaming
console. The **Portal** segment is the full interactive portal described here, one click
away. The chrome is an indigo gradient header with `Live / Portal / Insights / Report`
mode segments, a **Workspace** section rail on the left, a coordinated **multi-pane
Vitessce grid** in the centre, and the analysis tabs on the right. Fully static: no
backend, no API keys, no analytics.

## The dataset

**23 10x Visium sections, 43,747 in-tissue spots**, of spleen from mice infected with
*Plasmodium berghei* ANKA. The series is a protocol benchmark that also carries a real
genotype experiment, and the portal leads with the biology.

| Group | Sections | Spots | Median UMI | Role |
|---|---|---|---|---|
| **KO** — *Tbx21*^fl/fl^ *Cd23*^cre^ | 3 | 6,787 | 715 – 10,084 | Genotype contrast |
| **CTRL** — Cre-negative littermate | 3 | 5,839 | 469 – 10,490 | Matched comparator |
| **WT** — C57BL/6 | 17 | 31,121 | 3,594 – 36,142 | Protocol comparison |

*Cd23*-cre deletes *Tbx21* (T-bet) specifically in B cells. T-bet loss in B cells
disrupts germinal-centre polarization and antibody maturation, which is what makes this
a biologically grounded benchmark rather than a purely technical one.

Four preparation protocols are represented, and they are not interchangeable:

| Protocol | Placement | Sections | Feature space |
|---|---|---|---|
| Fresh-frozen OCT | Manual | 12 | 32,285 genes (poly-A, whole transcriptome) |
| FFPE | Manual | 7 | 19,465 genes (probe-based) |
| Fresh-frozen OCT | CytAssist | 2 | 19,465 genes (probe-based) |
| FFPE | CytAssist | 2 | 19,465 genes (probe-based) |

Per-section analysis runs in each section's native feature space. The integrated atlas
runs on the **intersection (19,465 genes)** — mixing the two un-intersected would make
protocol look like biology.

## Three things to read before trusting a number

**1. A spot is not a cell.** A 55 µm Visium spot pools roughly 1–10 cells, so a
"compartment" is the dominant programme of a spot, not a cell-type call. There is **no
deconvolution** here: GSE254652 describes a matched scRNA-seq reference, but only raw
sequencing was deposited (no count matrix), so compartments are marker-scored.

**2. Two sections are barely sequenced.** `OCT_pilot_CTRL_172` retains 47% of its spots
at a median 469 UMI, and `OCT_pilot_KO_167` 90% at 715 — against 6,429–10,490 in the
dedicated KO-vs-CTRL experiment. They are flagged `low depth` in the UI and excluded
from the primary contrast. This is consistent with the paper's own finding that the
fresh-frozen pilot performed worst on QC.

**3. The genotype contrast is n = 2 vs 2, and nothing survives FDR.** The gene ranking
in the KO vs CTRL tab is a ranking, not a significance claim. What makes it trustworthy
is the concordance panel: for the genes the Ritchie lab called significant in their own
limma-voom analysis, this independent pipeline reproduces the direction at **91–100%**
(Spearman ρ +0.79 to +0.94).

That concordance check runs on the **same counts and the same sections** the lab used,
so it is not independent replication of the biology — a shared artefact would be
reproduced, not caught. What differs is everything downstream: their limma-voom model
with RCTD/iSC-MEB cell-type labels versus this pipeline's pseudobulk Welch test with
marker-scored Leiden compartments. High agreement therefore says the published effects
are robust to analytical choices, and that this portal's numbers are trustworthy against
the paper of record.

## What the analysis found

**`Tbx21` falls by log2FC −1.31 in whole-spot pseudobulk** (p = 0.057, n = 2 v 2). The
lab's sorted germinal-centre B cells give −4.53 for the same comparison. Both are right:
*Cd23*-cre deletes *Tbx21* in B cells only, and a 55 µm spot also holds T and NK cells
that keep expressing it, so the whole-spot signal is necessarily diluted. The portal
states this rather than implying a clean knockout readout.

**The isotype switch is visible.** `Ighg2c` is the top CTRL-enriched transcript and
`Ighg1` / `Ighe` are among the top KO-enriched ones. T-bet directs class switching to
IgG2c, so losing it shifts the switched repertoire — the canonical T-bet phenotype,
recovered from spatial pseudobulk.

**Germinal-centre zonation resolves on Visium.** The integrated atlas splits the 2,575
germinal-centre spots into **1,109 dark zone** (Cxcr4-high, proliferative) and **1,466
light zone** (Cd83/Cd86-high). Zonation is a gradient *within* a single Leiden cluster,
so it is assigned per spot, not per cluster — a per-cluster argmax erases the split
entirely. A continuous `GC_zonation_axis` ships alongside the discrete call.

**The niches recover splenic architecture.** k-means over neighbourhood composition
returns red pulp, a B-follicle niche that is the GC-high one, a plasma-cell/T-zone
niche, a red-pulp-macrophage niche and a granulocyte niche — textbook white pulp / red
pulp organisation, derived without supervision.

## Contents

```
app.html                 self-contained portal (Plotly + marked from CDN)
vitessce.html            standalone coordinated multi-view (optional fallback)
data/
  manifest.json          {samples:[...], meta:{sid,gsm,protocol,placement,genotype,animal,
                          experiment,n_cells,low_depth}, pathways, signatures, paper_genes}
  precomputed.json       per-section counts, retention, marker genes, signature + pathway stats
  <sid>.json             {modality, spatial:[[x,y]], umap:[[x,y]],
                          obs:{compartment, leiden, niche, *_score, <pathway keys>}}
  <sid>.expr.json/.bin   quantized (uint8, gene-major) log-norm expression, 900 genes
  tissue/<sid>.jpg       cropped + display-enhanced H&E background, index.json = placement
  genotype.json          KO-vs-CTRL section means, pseudobulk DE, Tbx21 callout,
                         concordance vs the authors' published limma-voom, sensitivity
  atlas*.json/.bin       Harmony-integrated 43,489-spot atlas, compartments, pathways,
                         LR axes, Moran's I
  nhood.json             squidpy neighborhood-enrichment z-scores per section
  niche_summary.json     spatial-niche definitions (k-means over neighborhood composition)
  sq_stats.json          precomputed squidpy statistics per section (the Spatial tab)
  graph/<sid>.json       6-NN spatial graph as index pairs, in portal spot order
  h5ad/<sid>.h5ad        analysis-ready AnnData per section
  vitessce/              per-section AnnData-Zarr + configs + index.json
```

## Reproducing

```bash
PY=/home/ubuntu/analysis-env/bin/python3          # scanpy 1.12, squidpy 1.8.1, harmonypy 2.0.0

$PY pipeline/build_anndata.py       # GEO tar -> per-section h5ad + combined (intersected)
$PY pipeline/process.py             # QC -> Leiden -> compartments -> signatures -> portal data
$PY pipeline/tissue_export.py       # cropped, enhanced H&E backgrounds
$PY pipeline/genotype_export.py     # KO vs CTRL + concordance vs the published DE
$PY pipeline/integrate.py           # Harmony atlas, GC zonation, LR axes, Moran's I
$PY pipeline/spatial_analytics.py   # neighborhood enrichment + spatial niches
$PY pipeline/squidpy_stats.py       # centrality, co-occurrence, Ripley, sepal, histology
$PY pipeline/vitessce_export.py     # coordinated multi-view stores
```

`pipeline/patch_app.py` and `pipeline/patch_ported.py` retarget the shell and the ported
exporters from the May build. Every edit is an exact-string replacement asserted to apply
an exact number of times, so a drifted source fails loudly instead of half-patching.

## Sources

- GEO [GSE254652](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254652) — the Visium series (contact: `mritchie@wehi.edu.au`)
- [Du et al., *Genome Biology* 2025](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03543-4) — the SpatialBenchVisium paper
- [github.com/mritchielab/SpatialBenchVisium](https://github.com/mritchielab/SpatialBenchVisium) — the lab's analysis code and the published DE / GO tables used in the concordance panel
- [Zenodo 10.5281/zenodo.12683146](https://zenodo.org/doi/10.5281/zenodo.12683146) — the lab's `targets` project archive
- [Solano et al., bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.04.29.721531v1) + BioStudies [S-BSST2361](https://www.ebi.ac.uk/biostudies/studies/S-BSST2361) — the newer multi-platform SpatialBench (Visium HD / Xenium / MERSCOPE), not used here
