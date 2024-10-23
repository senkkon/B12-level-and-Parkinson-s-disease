## map gene names to their coordinates (GRCh38)

#step1
print("loading packages")
import argparse
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser()
parser.add_argument("-o","--output")
parser.add_argument("-i","--input")
parser.add_argument("-p","--pathname")

args = parser.parse_args()

path = args.input
with open(path,"r") as f:
    genes = f.read().strip()
genes = genes.split("\n")
pathname = args.pathname

if "GBA" in genes:
    print("identify GBA gene, renaming it to GBA1")
    for i in range(len(genes)):
        if genes[i] == "GBA":
            genes[i] = "GBA1"
    print(genes)
# gtf file can be downloaded from https://www.gencodegenes.org/human/

## GRCh38
ref=pd.read_csv("~/runs/go_lab/gencode/gencode.v44.annotation.gtf",sep = "\t",skiprows=5,header=None)
## GRCh37 (uncomment this for GRCh37)
#ref=pd.read_csv("~/runs/go_lab/gencode/gencode.v40lift37.annotation.gtf",sep = "\t",skiprows=5,header=None)

ref.columns = ["seqname","source","feature","start","end","score","strand","frame","attribute"]
ref = ref[ref.feature == "gene"]#only gene is the one we are looking at
##retrieve coordinates for genes
def get_gene_name(s):
    s = s.split(";")
    temp = list(map(lambda x: x.startswith(" gene_name"), s))
    s = pd.Series(s)
    gene_name = s[temp].iloc[0].split()[1].strip('"')
    return gene_name

ref["gene_name"] = list(map(lambda x: get_gene_name(x),ref.attribute))
ref = ref[["seqname","start","end","gene_name"]]
ref.columns = ["chr","start","end","gene"]


##get all genes and their coordinates
all = genes
all_df = ref[ref.gene.isin(set(all))]

### Problem solved
## check mismatched genes 49/50 matched. there is an empty string (GBA should be GBA1)
set(all) ^ set(list(all_df.gene))
# for chrX. rename it to chr23
print("WARNING! RENAMING chrX to chr23")
all_df.loc[all_df.chr=="chrX","chr"] = "chr23"

#make bed files
p = pathname
out = args.output
print("the results bed file")
print(all_df)
if all_df.shape[0] == len(genes):
    print("we get the correct bed file with input genes")
else:
    print(f"please check the gene names, there are {len(genes)} genes in the input, but we got {all_df.shape[0]} genes in the output ")
    exit()
all_df.to_csv(f"{out}/{p}.GRCh38.bed",sep = "\t", index=False, header=None)