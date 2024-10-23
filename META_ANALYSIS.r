# load package
library(metaviz)
library(tidyverse)
library(glue)
library(meta)

# load covariates and phenos
df = read_tsv("path/to/covariates_phenos")

#load PRS from .best file and convert to zscore.
cohorts <- c("APDGC","IPDGC","NIND","NGRC","PPMI","McGill","UKB") # "McGill","UKB"
PRS = data.frame()
for (cohort in cohorts){
  print(cohort)
  temp = read_delim(glue("{cohort}.best"),delim=" ")
  temp = mutate(PRS_z = scale(PRS))
  temp["cohort"] = cohort
  PRS = rbind(PRS,temp)
}


# merge with covariates and pheno
final = merge(df,PRS,by = c("IID","cohort"))
# Logistic regression to estimate the effect of PRS on PD risk

formula = "pheno ~ PRS_z + age + Sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"
risk_results = data.frame()
for (cohort_name in cohorts){
  data = final %>% filter(cohort == cohort_name1)
  print(dim(data))
  model = glm(formula, data = data,family = "binomial")
  summary_model <- summary(model)
  coef = summary_model$coefficients[2, 1]
  se = summary_model$coefficients[2, 2]
  t = summary_model$coefficients[2, 3]
  p = summary_model$coefficients[2, 4]
  risk_results <- rbind(risk_results, data.frame(cohort = cohort_name, t_value = t, p_value = p, BETA = coef, SE = se,Ncase = Ncase, Ncontrol=Ncontrol))
}
# convert beta to OR
risk_results = risk_results %>% mutate(OR = exp(BETA), OR_SE = OR * SE)
write.table(risk_results,"PD_B12_risk.csv",sep=",",row.names=F)

# META analysis
meta_result_risk <- metagen(TE = risk_results$BETA, # log-transformed ORs
                       seTE = risk_results$SE, # standard errors of ORs
                       studlab = risk_results$cohort, # study labels
                       method.tau = "DL", # DerSimonian-Laird (random effects)
                       sm = "OR") # Odds Ratio effect size


# Display summary of meta-analysis and make plots
summary(meta_result_risk)



viz_forest(x = risk_results[1:nrow(risk_results), c("OR", "OR_SE")], study_labels = risk_results[1:nrow(risk_results), c("cohort")],method = "DL",text_size=7,
           summary_label = "RANDOM EFFECT", xlab = "OR",variant = "classic",annotate_CI = T)



