#!/bin/bash



# this script calculates pathway PRS based on genes


gene_list=$1
output_dir=$2
cohort=$3 # make sure you've QCed your genotyping data
name=$4 # name of the pathway
covariate=$5 # first column need to be IID, but I dont think we will use the results from PRSice
pheno=$6 # first column need to be IID, but I dont think we will use the results from PRSice
target=$7 # you should make sure bim files contain rsid only instead of positions or rsid:alleles

output_folder=$output_dir/$name 
mkdir -p $output_folder


script_dir=.

if [[ -z $gene_list ]]; then echo "ERROR: gene_list (1st arg) not specified"; exit 42; fi
if [[ -z $output_dir ]]; then echo "ERROR: out directory (2nd arg) not specified"; exit 42; fi
if [[ -z $cohort ]]; then echo "ERROR: cohort name (3rd arg) not specified"; exit 42; fi
if [[ -z $name ]]; then echo "ERROR: name (4th arg) not specified"; exit 42; fi
if [[ ! -f $covariate ]]; then echo "ERROR: covariate $covariate does not EXIST"; exit 42; fi
if [[ ! -f $pheno ]]; then echo "ERROR: pheno $pheno does not EXIST"; exit 42; fi
if [[ -z $target ]]; then echo "ERROR: target (7th arg) not specified"; exit 42; fi

if [ ! -f bed_files/$name.GRCh38.bed ];then 
module load StdEnv/2020 python/3.8.10 scipy-stack/2020a
# most cohorts in our folders are based on GRCh38, UKB is on GRCh37
echo "getting positions of genes"
mkdir -p bed_files/
python ${script_dir}/map_genes_GRCh37.py -i $gene_list -o bed_files/ -p $name
python ${script_dir}/map_genes_GRCh38.py -i $gene_list -o bed_files/ -p $name
fi 



resource="-c 5 --mem=10g -t 2:30:0"
bed=bed_files/$name.GRCh38.bed


# bash ${script_dir}/pathway_prs_by_bed.sh $covariate $pheno $target $bed $output_folder $cohort
command="bash ${script_dir}/pathway_prs_by_bed.sh $covariate $pheno $target $bed $output_folder $cohort"
sbatch $resource --wrap "$command" --out $output_folder/$cohort.out --job-name=${name}_${cohort}

