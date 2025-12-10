# Simulate SNP effect sizes for 2 phenotypes that share (p) fraction of variants with a given genetic correlation (shared_corr)

library('MASS')
library('scales')

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 7) {stop("not enough Arguments received")}


num_SNPs= as.numeric(args[1])
p= as.numeric(args[2])
shared_corr= as.numeric(args[3])
SNPsLoc= args[4]

seed = as.numeric(args[5])
numHeterogSNPs= as.numeric(args[6])
outputName = args[7]

# Debug Varrs
# num_SNPs = 1000
# p= 0.5
# shared_corr = 1.0
# SNPsLoc= "C:/softwares/Cluster/0shaPRS/0ukbb/1"
# outputName = "C:/softwares/Cluster/0shaPRS/0ukbb/pheOut"
# seed=1
# numHeterogSNPs=5

set.seed(seed)
SNPs = read.table(SNPsLoc, header = F)
effectSize_var = 1.0 # variance is hardcoded to be 1, as these get scaled by LDAK later
num_shared = round(num_SNPs * p)
num_nonShared= num_SNPs - num_shared 

covMatrix=matrix(c(1, shared_corr  , shared_corr , 1), nrow=2)
data = mvrnorm(n=num_shared, mu=c(0, 0), Sigma=covMatrix, empirical=TRUE)
phenoA_shared_Beta = data[, 1] *sqrt(effectSize_var) # scale the variances after generating the correlated data, as otherwise mvrnorm would give an error of 'Sigma' is not positive definite'
phenoB_shared_Beta = data[, 2] *sqrt(effectSize_var) 


# confirm that correlation is what we wanted
corshared = cor(phenoA_shared_Beta, phenoB_shared_Beta)

phenoA_nonShared_Beta = rnorm(num_nonShared, 0,sd=sqrt(effectSize_var) )
phenoB_nonShared_Beta = rnorm(num_nonShared, 0,sd=sqrt(effectSize_var) )

corNonshared = cor(phenoA_nonShared_Beta, phenoB_nonShared_Beta) # no correlation

phenoA_Beta = c(phenoA_shared_Beta,phenoA_nonShared_Beta)
phenoB_Beta = c(phenoB_shared_Beta,phenoB_nonShared_Beta)
SNP_effects = cbind(phenoA_Beta,phenoB_Beta) # create final data table that has both the shared and nonshared SNP effects
overalCorr =  cor(phenoA_Beta, phenoB_Beta) # 0.3354613



message("num phenoA_nonShared_Beta: ",length(phenoA_nonShared_Beta), " | num phenoA_shared_Beta: ", length(phenoA_shared_Beta))

message("expected overall rG (p * shared_corr): ",round(p * shared_corr,2), " | actual rG: ", round(overalCorr,2))





filen =paste(outputName,"phenoA_B",".png", sep="" )
png(filen, width=500  , height=500);
# plot what such a correlation looks like
plot(phenoA_Beta,phenoB_Beta, ylab=paste("phenoB Beta (var: ",round(var(phenoB_Beta),2),")",sep=""), xlab=paste("phenoA Beta (var: ",round(var(phenoA_Beta),2),")",sep=""), main = paste("Pheno A and Pheno B SNP effects (rG: ",round(cor(phenoA_Beta, phenoB_Beta),2) ,")",sep="" ), type="n" )
points(phenoA_nonShared_Beta,phenoB_nonShared_Beta, col="black")
points(phenoA_shared_Beta,phenoB_shared_Beta, col = alpha("red", 0.5))

legend("topleft", legend=c("non-shared", "shared"),col=c("black","red"), text.col= c("black","red"), lty=1, lwd=3)

dev.off()

print(paste("written graphic  to",filen))


# Now split the 5 heterogeneous SNPs off from the 
regularSNP_effects_indices = 1:(nrow(SNP_effects)-numHeterogSNPs)
heterogSNP_effects_indices = (nrow(SNP_effects)-numHeterogSNPs +1):nrow(SNP_effects)
regularSNP_effects = SNP_effects[regularSNP_effects_indices,]
heterogSNP_effects = SNP_effects[heterogSNP_effects_indices,]

#create a dataframe that has the SNPs and their effects aligned
simDF = cbind.data.frame(SNPs, SNP_effects)
colnames(simDF) = c("SNP", "pheA_Beta", "pheB_Beta")

filen=paste0(outputName,"_SNPs_regular")
write.table(simDF[regularSNP_effects_indices,], filen, sep = " ", row.names = F, col.names = T, quote = FALSE)

filen=paste0(outputName,"_SNPs_extra")
write.table(simDF[heterogSNP_effects_indices,], filen, sep = " ", row.names = F, col.names = T, quote = FALSE)

# write out the actual SNP Ids, LDAK expect these as separate files, with the rsids on the rows, 
# two rows, for each pheno
LDAK_SNPlist =cbind.data.frame(simDF$SNP,simDF$SNP)  

filen=paste0(outputName,"_SNPs_regular_causals")
write.table( t(LDAK_SNPlist[regularSNP_effects_indices,]) , filen, sep = "\t", row.names = F, col.names = F, quote = FALSE)
filen=paste0(outputName,"_SNPs_extra_causals")
write.table(  t(LDAK_SNPlist[heterogSNP_effects_indices,]), filen, sep = "\t", row.names = F, col.names = F, quote = FALSE)


# write out effect sizes on rows (for LDAK)
filen=paste0(outputName,"_regular")
write.table(t(SNP_effects[regularSNP_effects_indices,]), filen, sep = " ", row.names = F, col.names = F, quote = FALSE)
print(paste("written regular effect sizes to",filen))

filen=paste0(outputName,"_extra")
write.table(t(SNP_effects[heterogSNP_effects_indices,]), filen, sep = " ", row.names = F, col.names = F, quote = FALSE)
print(paste("written extra effect sizes to",filen))

