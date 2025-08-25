##############################
# Sourcing GAPIT#
##############################

setwd("WBDC_geno")
source("WBDC_geno/gapit_functions.txt")

################
# Loading SNP Data #
################

myG <- read.delim("genotypes.hmp.txt", head = FALSE)

################
# Loading pheno Data #
################

myY <- read.table("phenotype_RGB.txt", head = TRUE)

###########################################################
# Running GAPIT #
############################################################
print("Starting GAPIT")

myGAPIT <- GAPIT( Y=myY[,c(1,2)],
                  G=myG,
                  PCA.total=2,   
                  model = c("MLM","FarmCPU", "Blink")
)

