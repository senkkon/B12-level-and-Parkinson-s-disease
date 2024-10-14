library(TwoSampleMR)
require(ggplot2)
library(dplyr)
library(MRPRESSO)
library(data.table)
library(httr)



### Loading exposure and outcome datasets.

# B12 levels, Paerkinson's disease (PD) risk and age of onset.
# Links to public datasets:
# B12 - http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90012001-GCST90013000/GCST90012772/
# PD_AAO - https://www.pdgenetics.org/resources (rsid column was generated separately)

B12 <- fread("Vi.B12.fastGWA")
PD_yesUKBB <- fread("META_no23_yesUKBB.txt")
PD_AAO <- fread("IPDGC_AAO_GWAS_sumstats_april_2018_rsid.txt")
setnames(PD_AAO, c("CHR", "POS", "rsid", "A1", "A2", "Freq1", "FreqSE", "MinFreq", "MaxFreq",
                            "Effect","StdErr","P-value","Direction", "HetISq","HetChiSq","HetDf","HetPVal"))

# PD progression phenotypes and allele reference.
# Datasets are available at https://pdgenetics.shinyapps.io/pdprogmetagwasbrowser/

CI_base_not_merged_with_ref <- fread("base_DEMENTIA.txt")
HY_not_merged_with_ref <- fread("cont_HY.txt")
MMSE_not_merged_with_ref <- fread("cont_MMSE.txt")
UPDRS3_not_merged_with_ref <- fread("cont_UPDRS3_scaled.txt")
reference <- fread("reference.txt")

CI_base <- merge(CI_base_not_merged_with_ref, reference, by="SNP")
HY <- merge(HY_not_merged_with_ref, reference, by="SNP")
MMSE <- merge(MMSE_not_merged_with_ref, reference, by="SNP")
UPDRS3 <- merge(UPDRS3_not_merged_with_ref, reference, by="SNP")



### Preparing datasets for analysis: renaming column names to the default names used in the TwoSampleMR package.

B12_renamed_for_MR <- B12 %>%
    rename(chr = chromosome, pos = base_pair_location, SNP = variant_id, effect_allele = A1, other_allele = A2, eaf = AF1,
           beta = BETA, pval = p_value, se = SE, samplesize = N)

PD_yesUKBB_renamed_for_MR <- PD_yesUKBB %>%
    rename(effect_allele = A1, other_allele = A2, eaf = freq, beta = b, pval = p, ncase = N_cases, ncontrol = N_controls) %>%
    mutate(chr = sub("chr", "", sub(":.*", "", chrpos)), pos = as.numeric(sub(".*:", "", chrpos)))

PD_AAO_renamed_for_MR <- PD_AAO %>%
    rename(chr = CHR, pos = POS, SNP = rsid, effect_allele = A1, other_allele = A2, eaf = Freq1, beta = Effect, pval = 'P-value',
           se = StdErr) %>%
    mutate(ncase = 16502,  ncontrol = 17996, samplesize = 34498)

CI_base_renamed_for_MR <- CI_base %>%
    rename(chrpos = SNP, SNP = RSID, chr = CHR, pos = START, effect_allele = ALT, other_allele = REF, eaf = MAF, beta = BETA,
           pval = P, se = SE, samplesize = N, gene = NearGENE)

HY_renamed_for_MR <- HY %>%
    rename(chrpos = SNP, SNP = RSID, chr = CHR, pos = START, effect_allele = ALT, other_allele = REF, eaf = MAF, beta = BETA,
           pval = P, se = SE, samplesize = N, gene = NearGENE)

MMSE_renamed_for_MR <- MMSE %>%
    rename(chrpos = SNP, SNP = RSID, chr = CHR, pos = START, effect_allele = ALT, other_allele = REF, eaf = MAF, beta = BETA,
           pval = P, se = SE, samplesize = N, gene = NearGENE)

UPDRS3_renamed_for_MR <- UPDRS3 %>%
    rename(chrpos = SNP, SNP = RSID, chr = CHR, pos = START, effect_allele = ALT, other_allele = REF, eaf = MAF, beta = BETA,
           pval = P, se = SE, samplesize = N, gene = NearGENE)



### Clumping

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')  # clumping's html options required for function to work

exp_data <- B12_renamed_for_MR %>% filter(B12_renamed_for_MR$pval < 5e-8)
exp_data <- format_data(data.frame(exp_data), type = "exposure", snp_col = "SNP", min_pval = 1e-200)
B12_clumped <- clump_data(exp_data, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1)



### Function for performing MR-PRESSO analysis and formatting its results

run_mr_presso <- function(sig) {
  tryCatch({
    presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", 
                        SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sig, 
                        NbDistribution = 1000, SignifThreshold = 0.05)

    capture.output(print(presso), file = paste(folder, "presso.txt", sep = "/"))
    
    # Extract values from presso results (with null handling)
    presso_p_dist <- ifelse(is.null(presso$`MR-PRESSO results`$`Distortion Test`$Pvalue), NA, presso$`MR-PRESSO results`$`Distortion Test`$Pvalue)
    presso_p_glob <- ifelse(is.null(presso$`MR-PRESSO results`$`Global Test`$Pvalue), NA, presso$`MR-PRESSO results`$`Global Test`$Pvalue)      
    
      presso_table <- data.frame(
          id.exposure = sig$id.exposure[1],
          id.outcome = sig$id.outcome[1],
          presso_p_dist = presso_p_dist,
          presso_p_glob = presso_p_glob
    )
    write.table(presso_table, file=file.path(paste(path_part, folder, sep=''), "presso_table.txt"), 
                row.names=FALSE, quote=FALSE, sep = "\t")
    
  },
  error = function(e) {
    capture.output(print(e), file = paste(folder, "presso_failed.txt", sep = "/")) 
    presso_table <- data.frame(
      id.exposure = sig$id.exposure[1],
      id.outcome = sig$id.outcome[1],
      presso_p_dist = NA,
      presso_p_glob = NA
    )
    write.table(presso_table, file=file.path(paste(path_part, folder, sep=''), "presso_table.txt"), 
                row.names=FALSE, quote=FALSE, sep = "\t")

  })
}



### MR analysis using TwoSampleMR and MR-PRESSO packages.

dir.create("MR_results")
setwd("MR_results")

exp <-list(B12_clumped)
out <- list(PD_yesUKBB_renamed_for_MR, CI_base_renamed_for_MR, HY_renamed_for_MR, MMSE_renamed_for_MR,
            UPDRS3_renamed_for_MR, PD_AAO_renamed_for_MR)

exp_names <- c('B12')
out_names <- c('PD_yesUKBB', 'CI_base', 'HY', 'MMSE', 'UPDRS3', 'PD_AAO')

# Setting the disease prevalence to calculate R-squared.
# References: PD - https://pubmed.ncbi.nlm.nih.gov/33121002/, CI - https://pubmed.ncbi.nlm.nih.gov/33121002/
prevalance_exp <- c(NA)
prevalance_out <- c(0.02802, 0.12, NA, NA, NA, NA)

for(i in 1:length(exp)){
   for(j in 1:length(out)){
   
      # Separate folder for each exposure-outcome pair. Harmonization of datasets
      folder <- paste(exp_names[i], " VS ", out_names[j], sep="")
      print(folder)
      dir.create(folder)
      exp_data <- exp[[i]]
      out_data <- format_data(data.frame(out[[j]]), type = "outcome", snps = exp_data$SNP)
      harm_data <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
      harm_data <-subset(harm_data, harm_data$eaf.exposure!="NA")

      # R-squared calculation for exposure and outcome. get_r_from_lor() was used for binary traits,
      # get_r_from_pn() for continuous traits. Since ncase and ncontrol are unknown for CI,
      # in this case we compute R using get_r_from_lor()
      harm_data$r.exposure <- get_r_from_pn(harm_data$pval.exposure, harm_data$samplesize.exposure)
      harm_data$units.exposure <-"SD units"
      if (j %in% c(1)){
          harm_data$prevalence.outcome <- prevalance_out[j]
          harm_data$units.outcome <-"log odds"
          harm_data$r.outcome <- get_r_from_lor(
              lor = harm_data$beta.outcome,
              af = harm_data$eaf.outcome,
              ncase = harm_data$ncase.outcome,
              ncontrol = harm_data$ncontrol.outcome,
              prevalence = harm_data$prevalence.outcome,
              model = "logit",
              correction = FALSE)
          }
      else {
          harm_data$r.outcome <- get_r_from_pn(harm_data$pval.outcome, harm_data$samplesize.outcome)
          harm_data$units.outcome <-"SD units" }

      # Steiger filtering
      steiger <- steiger_filtering(harm_data)
      sig <-subset(steiger, steiger$steiger_dir==TRUE)
      write.table(steiger, file=file.path(folder, "steiger.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
      run_mr_presso(sig)

      # MR analysis. Saving MR results (with calculated F-statistics) in a table
      mr_res <- mr(sig)
      mr_res_with_OR <- generate_odds_ratios(mr_res)
      n <- mean(sig$samplesize.exposure)
      k <- nrow(subset(sig, sig$ambiguous == FALSE))
      R2 <- mean(sig$rsq.exposure)
      F <- (R2*(n-1-k))/((1-R2)*k)
      capture.output(print(R2), file = paste(folder, "r2.txt", sep = "/"))
      capture.output(print(F), file = paste(folder, "f.txt", sep = "/"))
      mr_res_with_OR$R2 <- R2
      mr_res_with_OR$F <- F
      write.table(mr_res_with_OR, file=file.path(folder, "mr_res_with_OR.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
	  
      # Creating a report with MR main results
      mr_report(sig, study = folder, output_path = folder)

      # Saving sensitivity analysis results in a separate tables
      pleiotropy <- mr_pleiotropy_test(sig)
      heterogeneity <- mr_heterogeneity(sig)
      steiger_direction <- directionality_test(sig)
      write.table(pleiotropy, file=file.path(folder, "pleiotropy.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
      write.table(heterogeneity, file=file.path(folder, "heterogeneity.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
      write.table(steiger_direction, file=file.path(folder, "steiger_direction.txt"), row.names=FALSE, quote=FALSE, sep = "\t")

	  # MR-PRESSO analysis
      run_mr_presso(sig)

      # Saving plots with different sizing
      p1 <- mr_scatter_plot(mr_res, sig)
      ggsave(p1[[1]], file = paste("scatter_plot[", folder, "].png", sep=""), path=folder, width=7, height=7)
      mr_res_single <- mr_singlesnp(sig)
      p4 <- mr_funnel_plot(mr_res_single)
      ggsave(p4[[1]], file = paste("funnel_plot[", folder, "].png", sep=""), path=folder, width=7, height=7)
      p2 <- mr_forest_plot(mr_res_single)
      ggsave(p2[[1]], file = paste("forest_plot[", folder, "].png", sep=""), path=folder, width=7, height=12)
      mr_res_loo <- mr_leaveoneout(sig)
      p3 <- mr_leaveoneout_plot(mr_res_loo)
      ggsave(p3[[1]], file = paste("leaveoneout_plot[", folder, "].png", sep=""), path=folder, width=7, height=12)

      gc()
  }
}



### Collecting all tabulated MR results into one table. Archiving all results into a zip-file.

all_data <- list()
folder_names <- list.dirs(full.names = FALSE, recursive = FALSE)

for(i in 1:length(folder_names)){
    folder <- folder_names[[i]]
    mr_res_with_OR <- read.table(paste0(folder, "/mr_res_with_OR.txt"), sep="\t", header = TRUE)
    mr_res_with_OR$MR_name <- folder
    pleiotropy <- read.table(paste0(folder, "/pleiotropy.txt"), header = TRUE, sep = "\t")
    heterogeneity <- read.table(paste0(folder, "/heterogeneity.txt"), header = TRUE, sep = "\t")
    steiger_direction <- read.table(paste0(folder, "/steiger_direction.txt"), header = TRUE, sep = "\t")
    presso_table <- read.table(paste0(folder, "/presso_table.txt"), header = TRUE, sep = "\t")
    mr_res_with_OR <- mr_res_with_OR %>%
       full_join(pleiotropy, by = c("id.exposure", "id.outcome")) %>%
       full_join(heterogeneity, by = c("id.exposure", "id.outcome", "method")) %>%
       full_join(steiger_direction, by = c("id.exposure", "id.outcome")) %>%
       full_join(presso_table, by = c("id.exposure", "id.outcome"))
    all_data[[i]] <- mr_res_with_OR
}

MR_results_table <- do.call(rbind, all_data)
write.table(MR_results_table, "MR_results_table.tsv", sep="\t", row.names = FALSE)

setwd("..")
zip(zipfile = "MR_results.zip", 'MR_results')