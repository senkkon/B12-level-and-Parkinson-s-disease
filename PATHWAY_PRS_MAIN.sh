#!/bin/bash


script_dir=.
gene_list=B12.gene.list # txt file with gene names only
output_dir=../B12
name=B12


covariate=path/to/covariate
pheno=path/to/phenos
cohort=cohort_name
target=/path/to/QCed/genotyping/data
bash $script_dir/run_pathway_prs.sh $gene_list $output_dir $cohort $name $covariate $pheno $target 



