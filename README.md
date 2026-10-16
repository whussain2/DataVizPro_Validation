<p align="center">
  <img src="assets/banner.svg" alt="DataVizPro: validation and reproducibility" width="100%">
</p>

<p align="center">
  <a href="https://app.plantnura.com"><img alt="Web app" src="https://img.shields.io/badge/app-app.plantnura.com-7B3FE4?style=for-the-badge"></a>
  <a href="https://plantnura.com/tools/datavizpro.html"><img alt="Product page" src="https://img.shields.io/badge/product-DataVizPro-2E8B3E?style=for-the-badge"></a>
  <img alt="Paper" src="https://img.shields.io/badge/paper-Comput.%20Electron.%20Agric.%20(submitted)-667EEA?style=for-the-badge">
</p>
<p align="center">
  <img alt="DataVizPro" src="https://img.shields.io/badge/DataVizPro-v1.1.1-7B3FE4">
  <img alt="R" src="https://img.shields.io/badge/R-%E2%89%A5%204.5-2F7BD9?logo=r&logoColor=white">
  <img alt="Python" src="https://img.shields.io/badge/Python-3.10%2B-2F7BD9?logo=python&logoColor=white">
  <img alt="Licence" src="https://img.shields.io/badge/licence-MIT-2E8B3E">
  <img alt="Reproducible" src="https://img.shields.io/badge/results-fully%20reproducible-D63F8C">
</p>

<h3 align="center">DataVizPro: an integrated no-code cloud platform for publication-ready visualisation in crop science, from field trials to genomics</h3>
<p align="center"><b>Waseem Hussain</b> · PlantNura Consulting &amp; Solutions, Hyderabad, India · <a href="https://orcid.org/0000-0002-6861-0193">ORCID</a></p>

---

This repository is the **companion to the DataVizPro article**. It contains the code, data and outputs behind every computed result in the paper: the agreement of DataVizPro with independent reference software (**Table 4**) and the three case studies on public barley and rice data (**Sections 3.5–3.7**, **Figs 5–7**).

The analyses are run with the **DataVizPro computation engine itself**. The functions that the app's modules call are copied verbatim from the application (v1.1.1) into [`dvp_engine/`](dvp_engine). Anyone can therefore regenerate the published numbers without an account, then check them in the web app at [app.plantnura.com](https://app.plantnura.com).

> **DataVizPro in one line.** 31 modules and 150+ plot types in seven families (distribution, comparative, relationship, multivariate, heatmap with field-spatial analysis, genomic visualisation, special plots). Data preparation, genotype QC, analysis and publication-grade figures, all in the browser, with no code and no installation.

## Contents

- [Highlights](#highlights)
- [Validation against reference software (Table 4)](#validation-against-reference-software-table-4)
- [Case studies](#case-studies)
- [Quick start](#quick-start)
- [Repository structure](#repository-structure)
- [From the paper to the code](#from-the-paper-to-the-code)
- [The DataVizPro engine extract](#the-datavizpro-engine-extract)
- [Data](#data)
- [Software environment](#software-environment)
- [Citation](#citation)
- [Licence and contact](#licence-and-contact)

## Highlights

| | |
|:--|:--|
| **8 / 8** | analyses agree with independent implementations, with Pearson *r* ≥ 0.992 |
| **3 × 10⁻¹⁴** | largest difference between the DataVizPro relationship matrix and rrBLUP |
| **38.2 %** | of barley yield variance captured by the P-spline field trend (ΔAIC = 408) |
| **13 / 27** | top-10 % barley genotypes that stay in the top 10 % after spatial adjustment |
| **32,644** | rice SNPs × 413 accessions retained after DataVizPro genotype QC |
| ***GS3*** | seed-length association peak recovered on chromosome 3 (16.73 Mb) |

## Validation against reference software (Table 4)

Each DataVizPro output was compared with an independent implementation on identical inputs, after identical quality control (call rate ≥ 0.8, MAF ≥ 0.05).

| Analysis | DataVizPro function | Reference | Dataset | Pearson *r* | Max. abs. difference |
|:--|:--|:--|:--|--:|--:|
| LD *r*² | `ld_compute()` | squared genotype correlation (PLINK `--r2` definition) | Rice 44K, chr 1 (500 SNPs; 124,748 pairs) | 1.000000 | < 10⁻¹⁵ |
| LD *r*² | `ld_compute()` | snpStats 1.62 (EM haplotype *r*²) | same | 0.99987 | 0.28 ¹ |
| LD *D*′ | `ld_compute(…, "dprime")` | snpStats 1.62 (EM *D*′) | same | 1.000000 | 1.6 × 10⁻¹⁰ |
| LD decay (half-decay distance) | `ldd_make_entry()` | independent least-squares Hill–Weir fit | Rice 44K, 2 ancestry groups | — | 0.0001 kb |
| Ancestry coefficients, EM engine | `ps_run_parametric_admixture()` | fastmixture 1.3 (ADMIXTURE model) | Rice 44K, *K* = 5 | 0.99958 | 0.094 |
| Ancestry coefficients, sNMF engine | `ps_run_snmf()` | fastmixture 1.3 (ADMIXTURE model) | Rice 44K, *K* = 5 | 0.99813 | 0.21 |
| Genomic relationship matrix | `gh_kinship()` | VanRaden formula (base R); rrBLUP 4.6 | Rice 44K (85,078 pairs) | 1.000000 | 3.2 × 10⁻¹⁴ |
| PCA scores | `compute_pca()` | `prcomp` (base R) | Rice phenotypes (359 × 18) | 1.000000 | 0.009 ² |
| Regression coefficients | `fit_regression_model()` | statsmodels 0.15 (OLS, logistic) | Rice phenotypes | 1.000000 | 7.6 × 10⁻¹² |
| Spatially adjusted genotype values | `fs_fit_trend()` | SpATS 1.0 (PSANOVA) | Barley trial (272 genotypes) | 0.99178 | 0.14 ³ |

<sub>¹ DataVizPro reports composite (genotype-correlation) *r*², which is the PLINK definition. snpStats estimates haplotype *r*² by EM, and the two estimators legitimately differ only for pairs involving heterozygous or missing calls. ² FactoMineR standardises with 1/*n*, `prcomp` with 1/(*n* − 1). ³ The spline bases differ: mgcv uses a 10 × 10 tensor-product P-spline, SpATS PSANOVA uses 17 × 8 segments. Replicates are random in DataVizPro and fixed in SpATS. Spearman rank correlation = 0.989.</sub>

<p align="center"><img src="results/figures/fig_validation.png" alt="DataVizPro versus reference implementations" width="100%"></p>

The full table, with *n* for each comparison, is in [`results/tables/table4_agreement.csv`](results/tables/table4_agreement.csv).

## Case studies

### 1 · Spatial analysis of a barley field trial

The data are the row–column trial of Durban et al. (2003): 272 genotypes × 2 replicates on a 16 × 34 grid. The model has genotype fixed, replicate random and a tensor-product P-spline surface.

<p align="center"><img src="results/figures/fig_barley_spatial.png" alt="Barley field: raw yield, fitted trend and adjusted yield" width="100%"></p>

| Result | Value |
|:--|--:|
| Spatial variance share | 38.2 % |
| ΔAIC (no-spatial − spatial model) | 408.4 |
| Moran's *I* of residuals: design only → design + spatial | 0.255 (*P* < 10⁻¹⁶) → 0.037 (*P* = 0.19) |
| LISA plots significant at FDR 0.05: before → after adjustment | 1 → 0 |
| Spherical semivariogram: nugget / partial sill / range | 0.046 / 0.021 / ≈ 10 plots |
| Spearman correlation, raw vs adjusted genotype means | 0.863 |
| Top-10 % genotypes retained after adjustment | 13 of 27 |

Outputs are in [`case1_barley_summary.csv`](results/tables/case1_barley_summary.csv), [`case1_barley_genotypes.csv`](results/tables/case1_barley_genotypes.csv) and [`case1_barley_variograms.csv`](results/tables/case1_barley_variograms.csv).

### 2 · Genomic characterisation of the rice 44K diversity panel

The data are the panel of Zhao et al. (2011): 413 accessions and 36,901 SNPs, of which 32,644 SNPs remain after QC.

<p align="center"><img src="results/figures/fig_ancestry_K5.png" alt="Ancestry proportions at K = 5" width="100%"></p>

- **Population structure.** Genotypic PC1 and PC2 explain 36.5 % and 10.4 % of the variance. sNMF cross-entropy falls steeply up to *K* = 4 with no clear minimum. At *K* = 5, 361 accessions (87 %) have a dominant ancestry coefficient ≥ 0.7. They form five groups (G1–G5: 14, 102, 86, 58 and 101 accessions), and 52 accessions are admixed.
- **Kinship.** The mean VanRaden relationship is 1.14 within groups and −0.36 between groups.
- **LD decay.** The fitted Hill–Weir curve falls to half its initial value at 636 kb in G2 and 921 kb in G5.
- **GWAS** (rrBLUP mixed model, kinship + 3 PCs). For **seed length**, 36 SNPs on 9 chromosomes pass the Bonferroni threshold. One of the strongest signals lies at ***GS3*** (chromosome 3, 16.73 Mb). For **flowering time** (Arkansas), one SNP on chromosome 6 (9.65 Mb) passes, about 0.3 Mb from ***Hd1***. Inflation factors are λ = 0.93 and 0.85.

<p align="center">
  <img src="results/figures/fig_manhattan_seed_length.png" alt="Seed-length Manhattan plot" width="62%">
  <img src="results/figures/fig_ld_decay.png" alt="LD decay" width="36%">
</p>

### 3 · Multi-trait phenotypic relationships

The data are 18 agronomic, grain and quality traits of the same panel; 359 accessions have complete records.

- **Correlations.** 109 of 153 trait pairs are significant after Benjamini–Hochberg adjustment. The strongest are seed length vs brown-rice seed length (*r* = 0.977, 95 % bootstrap CI 0.972–0.982), seed volume vs surface area (0.909) and seeds vs florets per panicle (0.899).
- **PCA.** PC1 (25.3 %) reflects panicle size and phenology, and PC2 (19.5 %) reflects grain size and panicle number. Ancestry groups explain 32 % and 43 % of their variance.
- **Clustering.** The traits form five clusters (cophenetic *r* = 0.82). The smaller of two accession clusters is mostly ancestry group G2.

## Quick start

```bash
git clone https://github.com/whussain2/DataVizPro_Validation.git
cd DataVizPro_Validation

# 1. R packages (R >= 4.5)
Rscript install_packages.R

# 2. Python references (fastmixture and statsmodels)
pip install -r python/requirements.txt

# 3. Everything: Table 4, case studies, figures (about 30 min on 2 cores)
Rscript run_all.R
```

Each step can also be run on its own from the repository root, for example `Rscript scripts/01_barley_spatial.R`. The scripts call the Python references automatically when `python3` and `fastmixture` are on the PATH. Otherwise, run `python python/references.py` after `02_rice_qc.R` and `09_phenotypes.R`.

## Repository structure

```text
DataVizPro_Validation/
├── run_all.R                  # runs scripts 01-10 in order
├── install_packages.R         # CRAN + Bioconductor dependencies
├── dvp_engine/                # DataVizPro v1.1.1 computation functions (verbatim extract)
│   ├── load_engine.R
│   ├── fieldspatial.R         # field-spatial P-spline model, Moran's I, LISA, variograms
│   ├── ps_io.R · ps_inference.R   # genotype QC, EM / sNMF ancestry, PCA
│   ├── ld_heatmap.R · ld_decay.R  # pairwise LD (r2, D'), LD decay models
│   ├── genomics_heatmap.R     # genomic relationship matrices
│   ├── pca.R · pca_helpers.R  # PCA (FactoMineR)
│   ├── correlation.R · regression.R
├── scripts/
│   ├── 00_setup.R             # paths, engine loader, helpers
│   ├── 01_barley_spatial.R    # case study 1 + SpATS
│   ├── 02_rice_qc.R           # genotype QC
│   ├── 03_ancestry.R          # EM + sNMF vs fastmixture
│   ├── 04_kinship.R           # VanRaden vs formula and rrBLUP
│   ├── 05_genotype_pca.R
│   ├── 06_ld_heatmap.R        # r2 / D' vs snpStats
│   ├── 07_ld_decay.R          # Hill-Weir vs independent fit
│   ├── 08_gwas.R              # rrBLUP GWAS (seed length, flowering time)
│   ├── 09_phenotypes.R        # case study 3 + PCA + regression checks
│   └── 10_figures.R           # figures in results/figures
├── python/
│   ├── references.py          # fastmixture and statsmodels references
│   └── requirements.txt
├── data/                      # public datasets (see data/README.md)
│   ├── barley/Durban2003Barley.csv
│   └── rice/Geno44k.txt · sativas413.map · sativas413.fam · rice_phenotypes.csv
├── results/
│   ├── tables/                # every number reported in the paper, as CSV
│   ├── figures/               # figures shown above
│   └── sessionInfo.txt
├── assets/                    # banner and logos
├── CITATION.cff
└── LICENSE
```

## From the paper to the code

| Manuscript | Script | Output |
|:--|:--|:--|
| Table 4, SpATS row · Section 3.5 · Fig. 5 | `01_barley_spatial.R` | `case1_barley_*.csv`, `fig_barley_spatial.png` |
| Section 3.6: QC | `02_rice_qc.R` | `case2_rice_qc.csv` |
| Table 4, ancestry rows · Section 3.6 · Fig. 6D | `03_ancestry.R` | `case2_rice_ancestry_K5.csv`, `case2_rice_snmf_cross_entropy.csv`, `fig_ancestry_K5.png` |
| Table 4, relationship matrix · Fig. 6B | `04_kinship.R` | `case2_rice_kinship.csv` |
| Section 3.6: PCA · Fig. 6C | `05_genotype_pca.R` | `case2_rice_genotype_pca.csv` |
| Table 4, LD rows | `06_ld_heatmap.R` | `table4_agreement.csv` |
| Table 4, LD decay · Fig. 6E | `07_ld_decay.R` | `case2_rice_ld_decay.csv`, `fig_ld_decay.png` |
| Section 3.6: GWAS · Fig. 6F | `08_gwas.R` | `case2_rice_gwas_summary.csv`, `case2_rice_gwas_loci.csv`, `fig_manhattan_seed_length.png` |
| Table 4, PCA and regression · Section 3.7 · Fig. 7 | `09_phenotypes.R` | `case3_*.csv` |

## The DataVizPro engine extract

[`dvp_engine/`](dvp_engine) holds the computation functions behind the DataVizPro modules used in the paper, **copied verbatim** from the application source (v1.1.1, commit `af537d9`). The scripts call exactly the code that runs when a user clicks *Run* in the app.

- **What is included:** the analysis functions and their internal helpers: 70 functions, about 2,700 lines.
- **What is not included:** the user interface, the settings-driven figure engine, plan and account handling, and anything else from the web application. DataVizPro itself is available at [app.plantnura.com](https://app.plantnura.com), and reviewers and editors can get temporary full access on request.

## Data

All datasets are public and included in [`data/`](data) with their sources and licences ([`data/README.md`](data/README.md)):

- **Barley:** Durban et al. (2003), *J. Agric. Biol. Environ. Stat.* 8, 48–66, via the R package agridat.
- **Rice 44K:** Zhao et al. (2011), *Nat. Commun.* 2, 467 ([ricediversity.org](http://www.ricediversity.org/data/)). The phenotypes are the version redistributed in the R package enhancer.

## Software environment

| Component | Version used |
|:--|:--|
| R | 4.6.1 (application: 4.5.3) |
| mgcv · spdep · gstat | 1.9-1 · 1.4-2 · 2.1-6 |
| SpATS · LEA · snpStats | 1.0-20 · 3.24.0 · 1.62.0 |
| rrBLUP · FactoMineR · clue | 4.6.3 · 2.18 · 0.3-68 |
| fastmixture · statsmodels (Python) | 1.3.0 · 0.15.0 |

The full session information is in [`results/sessionInfo.txt`](results/sessionInfo.txt). The results are deterministic: random steps use fixed seeds (EM and sNMF: 42; fastmixture: 1–5; bootstrap: 1).

## Citation

If you use DataVizPro, these scripts or the results, please cite:

```bibtex
@article{hussain2026datavizpro,
  author  = {Hussain, Waseem},
  title   = {{DataVizPro}: an integrated no-code cloud platform for publication-ready
             visualisation in crop science, from field trials to genomics},
  journal = {Computers and Electronics in Agriculture},
  year    = {2026},
  note    = {Submitted}
}

@software{hussain2026datavizpro_app,
  author    = {Hussain, Waseem},
  title     = {{DataViz Pro}: publication-ready data visualization for plant science},
  version   = {1.1.1},
  year      = {2026},
  publisher = {PlantNura Consulting \& Solutions LLP},
  url       = {https://plantnura.com/tools/datavizpro.html}
}
```

GitHub's *Cite this repository* button uses [`CITATION.cff`](CITATION.cff).

## Licence and contact

The scripts, the engine extract and the documentation are released under the [MIT licence](LICENSE). The datasets keep the terms of their original sources. The DataVizPro web application, its interface and its trademarks are not covered by this licence.

Questions and reproducibility issues: please open an [issue](https://github.com/whussain2/DataVizPro_Validation/issues) or write to **waseemhussain@plantnura.com**.

<p align="center">
  <a href="https://plantnura.com"><img src="assets/plantnura_logo.png" alt="PlantNura" height="56"></a>
</p>
<p align="center"><sub>© 2026 PlantNura Consulting &amp; Solutions LLP · Hyderabad, India</sub></p>
