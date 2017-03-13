from chunkypipes.components import *
import os
import pandas
import matplotlib.pyplot as plt
import matplotlib.collections as collections

# need map file as well so convert bed, bim, fam, to ped/map
#./../TOOLS/plink --bfile CCPM-Ex_validation_02-20-17-related-removed --recode --out CCPM-Ex_validation_02-20-17-related-removed
# filter sample and snps by call rate
# ./../TOOLS/plink --file CCPM-Ex_validation_02-20-17-related-removed --geno 0.05 --mind 0.05 --make-bed --out CCPM-Ex_validation_02-20-17-related-removed-snp-and-sample-callRate-filter-applied

# LD-pruning via plink

# window size, variant step size, r^2 threshold
#./../TOOLS/plink --bfile CCPM-Ex_validation_02-20-17-related-removed-snp-and-sample-callRate-filter-applied --indep-pairwise 50 5 0.5 --out LDpruned

#het from plink to determine heterozygosity rates and inbreeding coefficeints (LD prune first!)