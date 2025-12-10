# Simulate SNP effect sizes for 2 phenotypes that share (p) fraction of variants with a given genetic correlation (shared_corr)

library('MASS')
library('scales')

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 6) {stop("not enough Arguments received")}


h1= as.numeric(args[1])
h2= as.numeric(args[2])
phe_regLoc= args[3]
phe_extraLoc= args[4]
seed = as.numeric(args[5])
outputName= args[6]


# Debug Varrs
# h1=0.45
# h2=0.05
# phe_regLoc= "C:/softwares/Cluster/0shaPRS/0ukbb/1_regular.pheno"
# phe_extraLoc= "C:/softwares/Cluster/0shaPRS/0ukbb/1_extra.pheno"
# outputName = "C:/softwares/Cluster/0shaPRS/0ukbb/pheFinal"
# seed=1


set.seed(seed)
phe_reg = read.table(phe_regLoc, header =F)
phe_extra = read.table(phe_extraLoc, header =F)
N=nrow(phe_reg)



# Phenotype simulation
# scale each of y1, y2, y3 to have variance 1
#generate final phenotype Y = sqrt(h1) * y1 + sqrt(h2) * y2 + sqrt(1-h1-h2) * y3
#(note var(Y) = h1 + h2 + 1-h1-h2 = 1)
#where h1 + h2 = total heritability
#h2/(h1+h2) is the proportion of that heritability you wish to come from the extreme snps

noise = rnorm(N)
pheA = sqrt(h1) * phe_reg$V3 + sqrt(h2) * phe_extra$V3 + sqrt(1-h1-h2) * noise

noise = rnorm(N)
pheB = sqrt(h1) * phe_reg$V4 + sqrt(h2) * phe_extra$V4 + sqrt(1-h1-h2) * noise


pheFinal = phe_reg
pheFinal$V3=pheA
pheFinal$V4=pheB

# sanity check if var(Y) == 1
var(pheFinal$V3)
var(pheFinal$V4)

# write out effect sizes on rows (for LDAK)
filen=paste0(outputName,".pheno")
write.table(pheFinal, filen, sep = " ", row.names = F, col.names = F, quote = FALSE)
print(paste("written final pheno to",filen))

