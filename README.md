# Analysis code: splice-event and gene-abundance layers in chronic allergic airway inflammation

Code accompanying the manuscript **"Splice-Event and Gene-Abundance Layers Are Independent and
Differ in Human Asthma Relevance"** (submitted to *Allergy, Asthma & Immunology Research*).

The repository contains the analysis scripts for

1. the gene-level and splice-event-level layers of a chronic ovalbumin mouse model,
2. the splicing-factor / RNA-binding-protein layer, transcription-factor layer and novel-transcript
   annotation check,
3. cross-species mapping to two human asthma cohorts (GSE152004, GSE158752),
4. signature scoring, receiver-operating-characteristic analysis, DeLong comparison of correlated
   areas under the curve, bootstrap confidence intervals, and type-2 stratification.

## Data

| Dataset | Accession | Use |
|---|---|---|
| Mouse RNA sequencing (this study) | processed data: https://doi.org/10.5281/zenodo.23504525 | discovery cohort |
| GSE295506 | mouse ovalbumin cohort | directional concordance |
| GSE152004 | human nasal airway epithelium, 441 asthma / 254 control | concordance, discrimination, type-2 strata |
| GSE158752 | human bronchial epithelial cells, 17 healthy / 50 asthma | replication, severity strata |

The processed tables needed by the scripts (count matrix, FPKM, DESeq2 output, rMATS outputs,
annotation tables) are deposited in the GEO series and in this repository under `data/`.

## Requirements

R 4.6.1 or later, with: `dplyr`, `tidyr`, `readr`, `purrr`, `tibble`, `DESeq2`, `pROC`,
`babelgene`, `clusterProfiler`, `fgsea`, `msigdbr`, `ggplot2`, `patchwork`. Package versions used
for the submitted results are recorded in `docs/session_info.txt`.

## How to run

Scripts are numbered and were run from the directory that contains `lib_theme.R` and `lib_data.R`.

```r
source("scripts/lib_theme.R")     # plotting theme, sample names, group labels
source("scripts/lib_data.R")      # loads the processed tables and derives all objects
```

| Script | What it does |
|---|---|
| `scripts/30_human_analysis.R` | cross-species mapping, concordance, severity, ROC, drug-target overlap; caches `human_results.rds` |
| `scripts/51_aair_extra.R` | AUC per cohort, bootstrap of the AUC difference, type-2 strata |
| `scripts/52_aair_delong.R` | fixed-orientation ROC, DeLong test, bootstrap interval, strata DeLong |
| `scripts/53_fig_t2.R` | Supplementary Figure S11 (ROC of both signatures within type-2 strata) |

`lib_data.R` expects a `DATA` variable pointing at the directory holding the processed tables; set it
at the top of the script or in your R profile:

```r
DATA <- "/path/to/processed_tables"
```

Output of the submitted run is in `docs/aair_delong_report.txt`.

## Notes and caveats

* Numbers in this repository reproduce those reported in the manuscript; the archived numbers in
  `docs/` are the exact console output of the submitted version.
* Two microarray cohorts used in an earlier version of this work (GSE67472, GSE43696) are **not**
  analysed here. Those series deposit probe-level matrices, and probe identifiers cannot be mapped
  to gene symbols without an additional platform annotation step, so they cannot support
  ortholog-based concordance as implemented. The human arm of the manuscript therefore rests on the
  two RNA-seq cohorts above.
* The permutation test for gene-set overlap uses hypergeometric draws from the 23,874 assembled
  genes; the reported two-sided P (0.302) equals twice the lower tail probability (0.299
  analytically).

## License

MIT (see `LICENSE`). If you use this code, please cite the article and the archived release
(`CITATION.cff`).
