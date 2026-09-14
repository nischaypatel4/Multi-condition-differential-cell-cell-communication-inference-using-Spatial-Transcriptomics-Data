# Multi-Condition Differential Cell–Cell Communication Inference Using Spatial Transcriptomics Data

**CS690 Computational Genomics — Final Project (IIT Kanpur, Nov 2025)**

------------------------------------------------------------------------

## Overview

Ligand–receptor (LR) signalling drives tumour progression, immune modulation, and tissue remodelling, but most cell–cell communication (CCC) tools (e.g. CellChat, CellPhoneDB) ignore spatial context and treat any cell as able to signal to any other. This project builds a framework that infers **differential** LR communication between matched normal and cancer spatial transcriptomics samples, while explicitly modelling **spatial proximity** and **condition-specific latent structure**.

We take spot-level LR interaction matrices produced by [COMMOT](https://github.com/zcang/COMMOT) and:

1.  Test each LR–cell-type combination for significant differences between conditions using a **Mann–Whitney U test**.
2.  Repurpose [**scDisInFact**](https://github.com/ZhangLabGT/scDisInFact) — originally built for disentangling gene-expression variation — to operate on LR interaction vectors instead of gene counts, in both a **non-spatial** and a **spatially-informed** variant.
3.  Generate normal-like counterfactual LR profiles for cancer edges and quantify a differential score (ΔLR = LR<sub>cancer</sub> − LR<sub>normal-predicted</sub>).
4.  Group significant LR pairs into six biologically interpretable **signalling programs** (immune, angiogenic, stromal/ECM, growth/survival, EMT/invasion, other) and visualise them as networks, dendrograms, and cell-type activity heatmaps.

### Research Questions

1.  How can differential cell–cell communication be inferred from multiple spatial transcriptomics datasets?
2.  How does cancer reshape ligand–receptor signalling programs, and what does spatial context reveal that non-spatial modelling misses?

------------------------------------------------------------------------

## Repository Contents

| File | Description |
|----|----|
| `Method_Implementation.ipynb` | End-to-end pipeline: COMMOT LR matrix extraction for both conditions, edge harmonisation, Mann–Whitney differential testing, LIANA+ exploration, and the non-spatial / spatial scDisInFact adaptations. |
| `Report.pdf` | Full write-up (abstract, methods, results, discussion, references, and per-member contributions). |
| `Presentation.pdf` | Slide deck summarising the motivation, methods, and results. |
| `SpatialPlots_Normal_Cancer.pdf` | Spot-level spatial visualisations (sender in blue, receiver in green, interaction strength as red edges) for selected LR pairs in normal vs. cancer tissue, including a negative-control pair (FGF7→FGFR2) that ranked low in the differential analysis. |

------------------------------------------------------------------------

## Data

Two Visium spatial transcriptomics datasets were obtained from [CZ CELLxGENE Discover](https://cellxgene.cziscience.com/):

- **Normal breast** — adult human breast atlas, Kumar et al., *Nature* 2023. [Dataset link](https://cellxgene.cziscience.com/collections/4195ab4c-20bd-4cd3-8b3d-65601277e731)
- **Breast cancer** — human breast cancer atlas, Wu et al., *Nature Genetics* 2021 (Sample `1160920F`). [Dataset link](https://cellxgene.cziscience.com/collections/dea97145-f712-431c-a223-6b5f565f362a)

Both samples used the Visium Spatial Gene Expression V1 platform (\~4,992 spots each, *Homo sapiens*).

### Input Format

COMMOT produces spot-by-spot adjacency matrices per LR pair. These are reshaped into a long-form edge table:

```         
ligand | receptor | sender | receiver | score | sender_celltype | receiver_celltype | sender_x | sender_y | receiver_x | receiver_y
```

Only LR pairs and edges detectable in **both** conditions (after cleaning section-specific spot-ID prefixes) are retained for downstream comparison. Gene symbols are annotated from Ensembl IDs using `mygene`.

------------------------------------------------------------------------

## Methods

### 1. Mann–Whitney U Test (non-parametric baseline)

For every LR pair and sender→receiver cell-type combination (with ≥3 edges per condition), we compare the distributions of COMMOT scores between cancer and normal using a two-sided Mann–Whitney U test, and compute a median-based log₂ fold-change:

```         
log2FC = log2((median_cancer + 1e-6) / (median_normal + 1e-6))
```

Interactions with `|log2FC| > 1` and `p < 0.05` are called significant. This method is fast and assumption-free but does not model condition-specific latent structure or spatial context.

### 2. LIANA+ (evaluated, not adopted)

LIANA+ was considered as an ensemble/consensus CCC framework but was found unsuitable here: it operates on **cell-type-aggregated** expression matrices and infers cell-type→cell-type communication, whereas this project requires **spot-level**, distance-aware inference. Forcing individual spots into LIANA+'s "cell-type" abstraction produces invalid statistics (near-zero group variance) and ignores spatial constraints entirely. LIANA+ code is retained in the notebook for reference/comparison.

### 3. Repurposed scDisInFact

[scDisInFact](https://github.com/ZhangLabGT/scDisInFact) is a generative model that disentangles gene-expression data into a **shared-biology** factor, a **condition-specific** factor, and a **batch** factor, enabling counterfactual prediction. We repurpose it by replacing its `cells × genes` input with a `sender_receiver-edge × LR-pair` matrix, treating each spot-pair as an "observation" and each LR score as a "feature." Two variants were trained:

- **Method 2A (non-spatial):** LR score matrix only.
- **Method 2B (spatial):** LR score matrix augmented with spatial descriptors per edge:
  - Euclidean sender→receiver distance
  - Log-transformed distance
  - k-nearest-neighbour (k=10) hotspot scores for sender and receiver spots
  - Combined edge-level spatial score (mean of sender/receiver hotspot scores)

For each cancer edge, a normal-condition counterfactual LR profile is decoded from its latent representation, and the differential score `ΔLR = LR_observed_cancer − LR_predicted_normal` is computed. LR pairs are then clustered (Ward hierarchical clustering) and mapped to six functional programs — **immune, angiogenic, stromal/ECM, growth/survival, EMT/invasion, other** — curated from literature on pathway-specific LR axes.

------------------------------------------------------------------------

## Key Results

**Mann–Whitney testing** found immune-regulatory signalling markedly *downregulated* in cancer (e.g. `GRN→SORT1`, `IL16→CD4`) alongside strong upregulation of growth/EMT pathways (`PGF→FLT1`, `PDGFC→PDGFRA`, `ANGPTL4→CDH11/SDC3`, `GAS6/PROS1→AXL`), consistent with angiogenesis, EMT, and immune suppression in the tumour microenvironment.

**Non-spatial vs. spatial scDisInFact** agreed on the top cancer-enriched axes (MDK-driven ECM remodelling, `NAMPT→ITGA5+ITGB1`, `GRN→SORT1`, `CXCL12→CXCR4`) but diverged in cell-type routing:

- The **non-spatial** model reported diffuse, tissue-wide signalling (e.g. broad B-cell and epithelial→myeloid crosstalk) and produced looser LR-pair clustering.
- The **spatial** model sharpened cluster boundaries, reassigned several interactions to spatially plausible short-range routes (e.g. epithelial→myeloid, plasmablast→epithelial, smooth-muscle→plasmablast), removed apparent "long-range" artefacts (e.g. epithelial→myeloid and B-cell→epithelial signals that disappear once spatial segregation is modelled), and revealed spatially confined niches such as a smooth-muscle desmoplastic reaction and tertiary-lymphoid-like B-cell restriction.

Consistently cancer-enriched pathways across methods:

- `MDK → SDC1/SDC2/SDC4/LRP1/NCL`
- `ANGPTL4 → syndecans (SDC3/SDC4/CDH11)`
- `CXCL12 → CXCR4`
- `PROS1/GAS6 → AXL`
- `PDGF/VEGF family → PDGFR/FLT1/KDR`

These map onto matrix remodelling, invasion, immune modulation, and angiogenesis — canonical hallmarks of tumour progression. Full figures (networks, dendrograms, per-program heatmaps) are in `Report.pdf` and `Presentation.pdf`; example spot-level visualisations are in `SpatialPlots_Normal_Cancer.pdf`.

------------------------------------------------------------------------

## Requirements

``` bash
pip install scanpy anndata commot mygene liana
pip install torch seaborn networkx scikit-learn pandas numpy scipy matplotlib

# scDisInFact (installed from source)
git clone https://github.com/ZhangLabGT/scDisInFact/
cd scDisInFact
pip install .
```

## Usage

Open `Method_Implementation.ipynb` and run sequentially. The notebook is organised into clearly marked sections:

1.  **COMMOT: generation of spot-level scores** — per-dataset LR matrix extraction and export to edge CSVs (`cancer_edges.csv`, `normal_edges.csv`).
2.  **Differential cell–cell communication analysis** — loading harmonised edges and helper functions (spot-ID cleaning, LR-string parsing, Ensembl→symbol mapping, program-level analysis).
3.  **Method 1: Mann–Whitney test.**
4.  **Method 2: repurposed scDisInFact** — dataset construction, non-spatial (2A) and spatial (2B) training, counterfactual prediction, and program-level aggregation.
5.  **LIANA+ attempt** — exploratory cell-type-level analysis kept for comparison/reference.

Note: paths in the notebook point to Kaggle input directories (`/kaggle/input/...`); update these to your local data locations before rerunning.

------------------------------------------------------------------------

## Conclusion

Integrating spatial proximity with condition-aware latent factor modelling (via a repurposed scDisInFact) yields a more mechanistically grounded, higher-resolution view of tumour communication than non-spatial or purely statistical approaches. The framework generalises to other spatial transcriptomics datasets and lays groundwork for future extensions to perturbation experiments and multimodal spatial data.

## References

Key methods and datasets referenced in this work: CellChat, CellPhoneDB, COMMOT, scDisInFact, LIANA+, and the CZ CELLxGENE Discover breast tissue atlases (Kumar et al. 2023; Wu et al. 2021). Full citations are provided in `Report.pdf`.
