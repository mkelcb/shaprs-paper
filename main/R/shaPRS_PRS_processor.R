##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 4) {stop("not enough Arguments received")}

subphenoLoc = args[1]
CombinedPhenoLoc = args[2]
blendFactorLoc = args[3]
outputLoc = args[4]


thresholds = c()
if (length(args) > 5) {
  thresholdedSNPsLoc = args[5]
  for (i in 6:length(args)) { # loop through the rest of the arguments, where each is supposed to be a threshold
    thresholds = c(thresholds,as.numeric(args[i]) )
  }
}
print(paste0("thresholds that will be processed for composite shaPRS:") )
print(thresholds )


# load exteral functions
source('/nfs/users/nfs_m/mk23/scripts/shaPRS.R')
library(qvalue)


# Debug Vars
# subphenoLoc="C:/softwares/Cluster/0shaPRS/0ibd/pheA_PLINK_PRS"
# CombinedPhenoLoc="C:/softwares/Cluster/0shaPRS/0ibd/pheCombined_PLINK_PRS"
# blendFactorLoc ="C:/softwares/Cluster/0shaPRS/0ibd/blendingFactors_PLINK"
# outputLoc ="C:/softwares/Cluster/0shaPRS/0ibd/testOut"
# thresholds=c(0.95,0.96,0.99)


# subphenoLoc="C:/softwares/Cluster/0shaPRS/0ibd/buga/GWAS_f1ldpred2_PRS_CD"
# CombinedPhenoLoc="C:/softwares/Cluster/0shaPRS/0ibd/buga/GWAS_f1ldpred2_PRS"
# blendFactorLoc ="C:/softwares/Cluster/0shaPRS/0ibd/buga/GWAS_f1_lFDR_meta_SNP_lFDR"
# outputLoc ="C:/softwares/Cluster/0shaPRS/0ibd/buga/testOut"
# thresholdedSNPsLoc="C:/softwares/Cluster/0shaPRS/0ibd/buga/GWAS_f1_lFDR_meta_SNPs"
# thresholds=c(0.96,0.99)


# subphenoLoc="C:/softwares/Cluster/0shaPRS/0ibd/buga/phenoA_sumstats"
# CombinedPhenoLoc="C:/softwares/Cluster/0shaPRS/0ibd/buga/Combined_sumstats"
# blendFactorLoc ="C:/softwares/Cluster/0shaPRS/0ibd/buga/myOutput_SNP_lFDR"

# test data
# subphenoLoc="C:/softwares/Cluster/0shaPRS/0ibd/uga/phenoA_sumstats_PLINK"
# CombinedPhenoLoc="C:/softwares/Cluster/0shaPRS/0ibd/uga/Combined_sumstats_PLINK"
#results$hardThresholds
#inputDataLoc = "C:/softwares/Cluster/0shaPRS/0ibd/uga/shapersToydata.txt" 
#hard_threshold_results = results$hardThresholds
#subpheno_PLINK = subpheno[,c(3,4,7)]
#CombinedPheno_PLINK = CombinedPheno[,c(3,4,7)]
#write.table(subpheno_PLINK,subphenoLoc, col.names = F, row.names = F, quote= F)
#write.table(CombinedPheno_PLINK,CombinedPhenoLoc, col.names = F, row.names = F, quote= F)
#subpheno = subpheno_PLINK
#CombinedPheno =CombinedPheno_PLINK

#write.table(results$lFDRTable,blendFactorLoc, col.names = T, row.names = F, quote= F)
#blendingFactColnames = colnames(blendingFactors)
#colnames(blendingFactors) = blendingFactColnames

#blendingFactors= read.table(blendFactorLoc, header = F)
#blendingFactors$V3 = c( rnorm() )



##############################################################
#                         MAIN SCRIPT                        #
##############################################################

# 1. load phenos
subpheno= read.table(subphenoLoc, header = F)

# load other subpheno to perform meta analysis
CombinedPheno= read.table(CombinedPhenoLoc, header = F)


blendingFactors= read.table(blendFactorLoc, header = T)

# 4. create new data frame to store the new summary stats
blendedSumstats = shaPRS_blend_PLINK(subpheno, CombinedPheno, blendingFactors)


# 5. write blended stats to disk
write.table(blendedSumstats, outputLoc, sep = "\t", row.names = F, col.names = F, quote = FALSE)
print(paste("written blended sumstats to",outputLoc))



hard_threshold_results= list()
# 6. process composite ones, if there were any
if (length(args) > 4) {
  # load list of heterogeneous SNPs for each threshold
  
  for (i in 1:length(thresholds)) { 
    heterogenousSNPs= read.table(paste0(thresholdedSNPsLoc,"_",thresholds[i]), header = F)
    hard_threshold_results[[i]] = heterogenousSNPs$V1
    print(paste0("loaded ", length(hard_threshold_results[[i]]), " heterogenous SNPs for threshold ",thresholds[i] ))
  }
  
  allCompositeResults = shaPRS_composite_PLINK(subpheno, CombinedPheno, hard_threshold_results)
  

  
  for (i in 1:length(thresholds)) { 
    currentThreshold = thresholds[i]
    filename = paste(outputLoc, "_comp_",currentThreshold , sep="")
    write.table( allCompositeResults[i], filename, sep = "\t", row.names = F, col.names = F, quote = FALSE)
    print(paste("for threshold : ",currentThreshold,", written PLINK PRS for ", nrow( allCompositeResults[[i]]), " SNPs to:", filename, sep=""))
  }
}
