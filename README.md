# B12-level-and-Parkinson-s-disease

In this project, we investigated the role of common and rare variants in genes related to B12 metabolism and the causal relationships between B12 levels and PD risk, age at onset, and motor/cognitive progression.

## Pathway-Specific Polygenic Risk Score

## Rare variants analysis

## Mendelian randomization (MR.r)

We applied a two-sample MR framework to examine whether B12 serum levels have a causal relationship with Parkinson's disease (PD) risk, age at onset and progression.

Main steps:

1. For the construction of IVs, we used GWAS significant (P-value < 5e-8) SNPs that are associated with B12 levels (exposure). Сlumping was performed with default settings; 
2. For the outcomes, we used PD GWAS full summary statistics on PD risk (with UK Biobank cohort; _PD_yesUKBB_), PD age at onset (_PD_AAO_) and several phenotypes related to PD progression: cognitive impairment (baseline analysis; _CI_base_), Mini-Mental State Examination scores (_MMSE_), Part III Unified Parkinson's Disease Rating Scale scores (_UPDRS3_) and Hoehn and Yahr Scale scores (_HY_);
3. MR analysis was performed in for loops, with one exposure and all outcomes;
4. We performed harmonization of datasets. Steiger filtering was performed to exclude SNPs that explain more variance in the outcome than in the exposure;
5. We calculated F-statistics and R2 for exposure;
6. The MR report was generated. The report included the results for an Inverse variance weighted (IVW, primary method), MR Egger, Weighted median, Simple mode and Weighted mode methods;
7. Heterogeneity was tested using Cochran’s Q test in the IVW and MR-Egger methods.
8. We performed MR-PRESSO test to detect horizontal pleiotropy and detect possible outliers;
9. To illustrate the results, we generated scatter plots, forest plots, funnel plots, and leave-one-out sensitivity analysis (LOO) plots.

The `MR.r` script contains all the code we used for the MR analysis. [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/) v0.6.4 and [MR-PRESSO](https://github.com/rondolab/MR-PRESSO) v1.0 R packages were utilized. The script was implemented using Rstudio.

## Genetic Correlation

The `LDSC.sh` script provides the complete workflow to calculate genetic correlation between exposure and outcome traits. This script automates the process of preparing summary statistics, munging the data for LDSC compatibility, and running the LDSC analysis to evaluate the genetic relationships between Vitamin B12 levels and Parkinson's Disease (PD) risk or age at onset.

**Key Features:**
- **Data Preparation:** Processes GWAS summary statistics to extract essential columns required for LDSC.
- **Munge Summary Statistics:** Formats the data using `munge_sumstats.py` to ensure compatibility with LDSC.
- **Run LDSC Analysis:** Executes `ldsc.py` to calculate genetic correlations, outputting results for further interpretation.
