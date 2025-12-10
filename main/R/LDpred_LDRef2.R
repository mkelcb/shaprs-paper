# pre compiling functions: https://stackoverflow.com/questions/23817341/faster-i-j-matrix-cell-fill

##############################


# Debug Vars
subphenoLoc="C:/0Datasets/shaPRS/0shaPRS/0ukbb/debug/1_pheA_sumstats"
subpheno_otherLoc="C:/0Datasets/shaPRS/0shaPRS/0ukbb/debug/1_pheB_sumstats"
blendFactorLoc ="C:/0Datasets/shaPRS/0shaPRS/0ukbb/debug/1_pheBlend_lFDR_meta_SNP_lFDR"
outputLoc ="C:/0Datasets/shaPRS/0shaPRS/0ukbb/debug/1_pheBlend_sumstats_meta"

##############################################################
#                         MAIN SCRIPT                        #
##############################################################

# 1. load phenos
subpheno= read.table(subphenoLoc, header = T)

# load other subpheno to perform meta analysis
subpheno_other= read.table(subpheno_otherLoc, header = T)

blendingFactors= read.table(blendFactorLoc, header = T)

# Here execute 'shaPRS_blend_overlap = function(subpheno, subpheno_other, blendingFactors, rho = 0) {', until we get
# just after line: '# remove any NAs of SNPs' subpheno_otherPheno_blending
head(subpheno_otherPheno_blending)
neededCols = c(1,7,8,16,17,20)
dataSubset = subpheno_otherPheno_blending[1:5,neededCols]
set.seed(1)
dataSubset$lFDR = runif(nrow(dataSubset))

# create fake correlations for 5 SNPs for pop 1
pop1 = rbind(               c(1,0.5,0.2,0,0.3) )
pop1 = rbind(pop1, c(0.5,1,1,1,1) )
pop1 = rbind(pop1, c(0.2,1,1,1,0) )
pop1 = rbind(pop1, c(0,1,1,1,-0.7) )
pop1 = rbind(pop1, c(0.3,1,0,-0.7,1) )
pop1LDmatrix = new("dsCMatrix",  Matrix(pop1, sparse = TRUE) )

# same for pop 2
pop2 = rbind(      c(1  ,0.2,0.8,0   ,0.6 ) )
pop2 = rbind(pop2, c(0.2,1  ,0.1,0.2 ,0.9 ) )
pop2 = rbind(pop2, c(0.8,0.1,1  ,0.5 ,0   ) )
pop2 = rbind(pop2, c(0  ,0.2,0.5,1   ,-0.4) )
pop2 = rbind(pop2, c(0.6,0.9,0  ,-0.4,1   ) )
pop2LDmatrix = new("dsCMatrix",  Matrix(pop2, sparse = TRUE) )
#################################################################

remove(subphenoLoc)

remove(subpheno_otherLoc)
remove(blendFactorLoc)
remove(outputLoc)
remove(subpheno)


remove(subpheno_other)
remove(blendingFactors)
remove(neededCols)
remove(pop1)
remove(pop2)
remove(subpheno_otherPheno)
remove(subpheno_otherPheno_blending)
remove(misalignedAlleleIndices)



############################################

















