##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 6) {stop("not enough Arguments received")}

rho = as.numeric(args[1])
subphenoLoc = args[2]
subpheno_otherLoc = args[3]
blendFactorLoc = args[4]
outputLoc = args[5]
counter = args[6]

# load exteral functions
source('/nfs/users/nfs_m/mk23/scripts/shaPRS.R')
library(qvalue)


# Debug Vars
#rho=0.2852094
# subphenoLoc="C:/0Datasets/shaPRS/corrDiag/"
# subpheno_otherLoc="C:/0Datasets/shaPRS/corrDiag/"
# blendFactorLoc ="C:/0Datasets/shaPRS/corrDiag/GWAS_f1_lFDR_meta_SNP_lFDR"
# outputLoc ="C:/0Datasets/shaPRS/corrDiag/corrDiag"
#counter=1


##############################################################
#                         MAIN SCRIPT                        #
##############################################################

# 1. load phenos
subpheno= read.table(subphenoLoc, header = T)

# load other subpheno to perform meta analysis
subpheno_other= read.table(subpheno_otherLoc, header = T)


# INVERSE VARIANCE FIXED EFFECT META ANALYSIS:
CombinedPheno = inverse_metaAnalaysis(subpheno, subpheno_other)


blendingFactors= read.table(blendFactorLoc, header = T)



subpheno_CombinedPheno = merge(subpheno,CombinedPheno,by.x = "SNP",by.y = "SNP")
subpheno_CombinedPheno_blending = merge(subpheno_CombinedPheno,blendingFactors, by.x = "SNP", by.y = "SNP")

# 2. Align PheB/B alleles
misalignedAlleleIndices = which( as.character(subpheno_CombinedPheno_blending$A1.x) != as.character(subpheno_CombinedPheno_blending$A1.y) ) # compare as character, as if we have non-SNPs with different alleles factors will break
subpheno_CombinedPheno_blending$b.y[misalignedAlleleIndices] = -subpheno_CombinedPheno_blending$b.y[misalignedAlleleIndices] # flip effects
if(length(misalignedAlleleIndices) > 0) message(paste0(length(misalignedAlleleIndices)), " misaligned allele(s) effects were flipped" )

# sanitize each input of NAs
subpheno_CombinedPheno_blending <- na.omit(subpheno_CombinedPheno_blending) # remove any NAs of SNPs


cor_orig = cor(subpheno_CombinedPheno_blending$b.x,subpheno_CombinedPheno_blending$b.y)
cor_new = cor(subpheno_CombinedPheno_blending$b.x/subpheno_CombinedPheno_blending$se.x,subpheno_CombinedPheno_blending$b.y/subpheno_CombinedPheno_blending$se.y)


sigmaBar_orig = sqrt(  (1-subpheno_CombinedPheno_blending$lFDR)^2 * subpheno_CombinedPheno_blending$se.x^2 +  subpheno_CombinedPheno_blending$se.y^2 * subpheno_CombinedPheno_blending$lFDR^2 + 2*subpheno_CombinedPheno_blending$lFDR*(1-subpheno_CombinedPheno_blending$lFDR)*cor(subpheno_CombinedPheno_blending$b.x,subpheno_CombinedPheno_blending$b.y)* subpheno_CombinedPheno_blending$se.x*subpheno_CombinedPheno_blending$se.y )



#(1-lFDR)^2 * se.x^2 +  se.y^2 * lFDR^2 + 2*lFDR*(1-lFDR)* cor(b.x/se.x, b.y/se.y) * se.x*se.y 

#(1-lFDR)^2 * se.x^2 +  se.y^2 * lFDR^2 + 2*lFDR*(1-lFDR)* cor(b.x,b.y) * se.x*se.y 


# X-axis
CovAB_orig = with(subpheno_CombinedPheno_blending, cor(b.x,b.y) * se.x*se.y )

# Theoretical against this on Y axis:
CovAB_theory = with(subpheno_CombinedPheno_blending, (1 + rho *se.x/se.y) / (1/se.x^2 + 1/se.y^2 ))


# another Y-axis
CovAB_Empirical = with(subpheno_CombinedPheno_blending, cor(b.x/se.x, b.y/se.y) * se.x*se.y )



# 5. write blended stats to disk

write.table(cbind(cor_orig,cor_new), outputLoc, sep = "\t", row.names = F, col.names = F, quote = FALSE,append=TRUE)
print(paste("written diagnostics to",outputLoc))





filename=paste0(outputLoc,"_CovAB_CovAB_theory.png",counter)
png(filename, width=640 , height=640, res =128);
plot(CovAB_orig,CovAB_theory)
abline(a = 0, b = 1)
dev.off()

filename=paste0(outputLoc,"_CovAB_CovAB_Empirical.png",counter)
png(filename, width=640 , height=640, res =128);
plot(CovAB_orig,CovAB_Empirical)
abline(a = 0, b = 1)
dev.off()

filename=paste0(outputLoc,"_CovAB_Empirical_CovAB_theory.png",counter)
png(filename, width=640 , height=640, res =128);
plot(CovAB_Empirical,CovAB_theory)
abline(a = 0, b = 1)
dev.off()

cor(CovAB_orig,CovAB_Empirical) # 1
cor(CovAB_orig,CovAB_theory) # 0.9996641
cor(CovAB_Empirical,CovAB_theory) # 0.9996641

filen=paste0(outputLoc,"_CovAB",counter)
write.table(cbind(CovAB_orig,CovAB_theory, CovAB_Empirical), filen, sep = "\t", row.names = F, col.names = F, quote = FALSE)
print(paste("written diagnostics to",filen))


