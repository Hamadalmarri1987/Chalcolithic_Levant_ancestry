# Chalcolithic Levant Ancestry: Mesopotamian Genetic Contribution in Peqi'in

This repository contains the analysis code for our study on the genetic ancestry of Chalcolithic Levant populations, with a focus on the Peqi'in population. We investigate the genetic contributions from neighboring regions including Mesopotamia, Anatolia, and the Zagros using ancient DNA data.

## Citation

**Manuscript Status:** in final stages

If you use this code or findings in your research, please cite:

(will write citation here after publication / if doi is available after publishing it in the bioarchives)

**Note:** Citation will be updated with DOI and publication details upon acceptance.

---

## Overview

This project uses multiple population genetic methods to model the ancestry of Israel Chalcolithic (Israel_C) populations:

- **Principal Component Analysis (PCA)** to visualize genetic relationships
- **qpWave** to test for cladality between ancient populations
- **qpAdm** to model admixture proportions from source populations
- **Rare allele sharing** to detect fine-scale genetic affinities
- **ADMIXTURE** for unsupervised clustering analysis

The analyses use the Allen Ancient DNA Resource (AADR) v62.0 dataset with both 1240k and Human Origins SNP panels.

---

## Repository Structure

```
.
├── PCA_aDNA.R                      # Principal component analysis
├── RareAllele_Admixture.R          # Rare allele sharing + ADMIXTURE visualization
├── qpAdm_complete_analysis.R       # Comprehensive 3-way and 4-way qpAdm models
├── QpWave/                         # Regional qpWave analyses
│   ├── qpwave_iran.R
│   ├── qpwave_levant.R
│   ├── qpwave_anitolia.R
│   ├── Iran_iraq_qpwave.R
│   ├── iraq_anitolia_qpwave.R
│   ├── qpwave_iran_chalcolithic.R
│   └── qpave_iraq.R
└── figure4/                        # Publication figure generation scripts
    ├── figure4a.py
    ├── figure4b.py
    ├── figure4c.py
    ├── figure4d.py
    └── graph_merge/
```

---

## Requirements

### Software Dependencies

- **R** (≥ 4.3.0)
  - `admixtools` (v2.0.10 or later)
  - `tidyverse`
  - `ggplot2`, `ggrepel`, `ggtext`
  - `RColorBrewer`
  - `data.table`
  - `BEDMatrix`
  - `pheatmap`
  - `writexl`, `readxl`

- **Python** (≥ 3.8)
  - `pandas`
  - `numpy`
  - `matplotlib`
  - `openpyxl` (for Excel file handling)

- **EIGENSOFT** (for `smartpca` and `convertf`)
- **PLINK** (v1.9 or v2.0)

### Data Requirements

Download the AADR v62.0 dataset from [Reich Lab](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW):

- `v62.0_1240k_public.{geno,snp,ind,anno}`
- `v62.0_HO_public.{geno,snp,ind,anno}`

Expected directory structure:
```
AADR_v62/
├── v62.0_1240k_public.geno
├── v62.0_1240k_public.snp
├── v62.0_1240k_public.ind
├── v62.0_1240k_public.anno
└── (similar files for HO dataset)
```

## Installation

### R Package Installation

```r
# Install CRAN packages
install.packages(c("tidyverse", "ggplot2", "ggrepel", "ggtext", 
                   "RColorBrewer", "data.table", "BEDMatrix", 
                   "pheatmap", "writexl", "readxl", "magrittr"))

# Install admixtools from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("uqrmaie1/admixtools")
```

### Python Package Installation

```bash
pip install pandas numpy matplotlib openpyxl
```

### External Tools

- **EIGENSOFT:** Follow installation instructions at [github.com/DReichLab/EIG](https://github.com/DReichLab/EIG)
- **PLINK:** Download from [www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)

Ensure `smartpca`, `convertf`, and `plink` are in your system PATH.

---

## Usage

### Step 1: Configure Paths

Before running any scripts, update the data paths in each file to match your local setup:

**In R scripts:**
```r
# Change this line to your AADR directory
aadr_path <- "path/to/your/AADR_v62/"
results_dir <- "path/to/results_aDNA/"
```

**In Python scripts:**
```python
# Update Excel file paths
XLSX_PATH = r"path/to/your/results/S1_qpAdm.xlsx"
```

### Step 2: Run qpWave Analyses (Optional)

These scripts test which populations can be merged based on cladality:

```r
# Test Iran populations
Rscript QpWave/qpwave_iran.R

# Test other regions
Rscript QpWave/qpwave_levant.R
Rscript QpWave/qpwave_anitolia.R
Rscript QpWave/Iran_iraq_qpwave.R
```

**Output:** Console output showing p-values for rank tests. Use these results to inform population merging in qpAdm.

**Runtime:** ~30-60 minutes per script depending on SNP overlap.

### Step 3: Principal Component Analysis

```r
Rscript PCA_aDNA.R
```

**What it does:**
- Computes PCs on modern populations
- Projects ancient Near Eastern samples onto PC space
- Generates publication-quality PCA plots

**Outputs:**
- `results_aDNA/pca_analysis/Figure_PCA_HO_projected.png`
- `results_aDNA/pca_analysis/HO_pca.evec` (eigenvectors)
- `results_aDNA/pca_analysis/HO_pca.eval` (eigenvalues)

**Runtime:** ~1-2 hours (requires running `smartpca`)

### Step 4: Rare Allele Sharing Analysis

```r
Rscript RareAllele_Admixture.R
```

**What it does:**
- Filters SNPs to rare variants (MAF 0.01-0.05)
- Computes pairwise rare allele sharing
- Visualizes ADMIXTURE cross-validation results
- Generates population-level admixture barplots

**Outputs:**
- `results_aDNA/rare_allele_sharing_results.xlsx`
- `results_aDNA/out_population_rare_sharing_heatmap.pdf`
- `results_aDNA/admixture_barplots_by_population2.pdf`
- `results_aDNA/Admixture_CV.pdf`

**Runtime:** ~2-4 hours (rare allele computation is intensive)

### Step 5: qpAdm Modeling

```r
Rscript qpAdm_complete_analysis.R
```

**What it does:**
- Tests 3-way models (Levant + Anatolia + Iran/Mesopotamia)
- Tests 4-way models (Levant + Anatolia + Iran + Mesopotamia)
- Uses merged populations based on qpWave results
- Tests multiple outgroup configurations

**Outputs:**
- `results_aDNA/qpadm_3way_results.xlsx`
- `results_aDNA/qpadm_4way_results.xlsx`
- `results_aDNA/qpadm_priority_models.xlsx`

**Runtime:** 4-8 hours for full analysis (hundreds of models tested)

**Note:** The script includes both simplified and comprehensive model sets. For initial exploration, use the priority models section.

### Step 6: Generate Publication Figures

```bash
# Generate Figure 4 panels
python figure4/figure4a.py
python figure4/figure4b.py
python figure4/figure4c.py
python figure4/figure4d.py
```

**Requirements:** These scripts expect Excel files from qpAdm analysis in a specific format. Update the `XLSX_PATH` variable in each script.

---

## Key Populations

### Target Population
- **Israel_C.AG** (Chalcolithic Levant, Peqi'in Cave)

### Source Populations (by region)

**Levant:**
- Jordan_PPNB.AG
- Israel_PPNB.AG
- Jordan_Baja_Late_PPNB.AG

**Anatolia:**
- Turkey_Marmara_Barcin_N.SG
- Turkey_Central_Catalhoyuk_N.SG
- Turkey_Central_Boncuklu_PPN.WGC.SG

**Zagros (Iran):**
- Iran_GanjDareh_N.AG
- Iran_HajjiFiruz_N.AG
- Iran_TepeAbdulHosein_N.SG
- Iran_Wezmeh_N.SG
- Iran_SehGabi_C.AG (Chalcolithic)

**Mesopotamia:**
- Iraq_Shanidar.AG
- Iraq_PPNA.AG
- Turkey_Southeast_Cayonu_PPN.SG
- Turkey_Southeast_Mardin_PPN.AG

### Outgroup Populations

Multiple outgroup sets are tested (see `qpAdm_complete_analysis.R` for details):
- Ancient West Eurasians (Kostenki14, Ust'Ishim, MA1)
- East Asians (Han, Papuan)
- Native Americans (Karitiana)
- Africans (Mbuti, Ethiopia_4500BP)
- Regional hunter-gatherers (Natufian, Pinarbasi, Kotias)

## Data Availability

- **AADR v62.0:** [Reich Lab Downloads](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data)

---

## Acknowledgments

This work uses data from the Allen Ancient DNA Resource (AADR) curated by the Reich Lab at Harvard Medical School. We thank the admixtools developers for their excellent software package.

---
## Data Availibility Statement
The datasets analysed during the current study are publicly available in the Harvard Dataverse repository at https://doi.org/10.7910/DVN/FFIDCW.
No new primary data were generated. All scripts and analytical workflows used for data processing and statistical analyses are available in the GitHub repository: https://github.com/Hamadalmarri1987/Chalcolithic_Levant_ancestry

---

## License

This code is released under the MIT License. See LICENSE file for details.

---

## Contact

For questions about the analysis or to report issues:
- Open an issue on this repository
- Contact: [corresponding author email - add upon publication]

---

**Last updated:** February 2026
