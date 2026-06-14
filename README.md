# Co-mutation and Synthetic Lethality Analysis

Pan-cancer co-mutation pattern discovery and synthetic lethal interaction nomination using the MSK-IMPACT 50K cohort (53,654 tumors).

## Repository Contents

| File | Description |
|---|---|
| `02_co_mutation_analysis.R` | Pairwise Fisher's exact test, FDR correction, mutual exclusivity detection |
| `03_survival_analysis.R` | Four-group Cox regression for co-mutated gene pairs |
| `04_figures_network.R` | Co-mutation network visualization and publication figures |
| `05_tcga_validation.R` | TCGA pan-cancer replication of co-mutation patterns |

## Key Features

- Panel-corrected co-mutation analysis (IMPACT468+505 only)
- 187 significant co-mutation pairs identified (FDR<0.05)
- Synthetic lethal candidate nomination (OR<0.75, pathway convergence, drug targetability)
- Classical validation pairs confirmed (KRAS-EGFR, KRAS-BRAF)
- Network visualization via igraph

## Data Sources

- MSK-IMPACT 50K (Bandlamudi et al., *Cancer Cell* 2026)
- TCGA Pan-Cancer (via cBioPortal/GDC)

## Author

Wenjie Cai (蔡文杰), Department of Radiation Oncology, First Hospital of Quanzhou Affiliated to Fujian Medical University
