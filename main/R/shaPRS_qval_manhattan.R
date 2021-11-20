if (!require("qqman")) install.packages("qqman")
library("qqman")
library(scales)
##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 3) {stop("not enough Arguments received")}

inputDataLoc = args[1]
blendFactorsLoc = args[2]
outputLoc = args[3]
plotTitle= args[4]
plotTitle = gsub("_", " ", plotTitle) #remove underscores
B12Loc= args[5]

# DebugVars
 


#  inputData= read.table("C:/0Datasets/shaPRS/ibd/Qvals_new_diagDat", header = T)

##############################################################
#                         MAIN SCRIPT                        #
##############################################################

# 1. load data
inputData= read.table(inputDataLoc, header = T)
blendFactors= read.table(blendFactorsLoc, header = T)
B12= read.table(B12Loc, header = T)
inputData_blendFactors = merge(inputData,blendFactors, by ="SNP") # merge to make sure they are aligned
B12_blendFactors = merge(B12,blendFactors, by ="SNP") # merge to make sure they are aligned

# merge all 3, the two phenos raw, the lFDR, and the shaPRS coefs
inputData_blendFactors_shaPRS = merge(inputData_blendFactors,B12, by ="SNP") # merge to make sure they are aligned


### FIND OUT WHY WE DONT HAVE THE SAME DATA IN IBD

# remove non numeric data
inputData_blendFactors_shaPRS <- inputData_blendFactors_shaPRS[!is.na(as.numeric(as.character(inputData_blendFactors_shaPRS$Beta_B))),]
inputData_blendFactors_shaPRS <- inputData_blendFactors_shaPRS[!is.na(as.numeric(as.character(inputData_blendFactors_shaPRS$Beta_A))),]
inputData_blendFactors_shaPRS <- inputData_blendFactors_shaPRS[!is.na(as.numeric(as.character(inputData_blendFactors_shaPRS$b))),]


# now actually cast them to numeric
inputData_blendFactors_shaPRS$Beta_B = as.numeric(as.character(inputData_blendFactors_shaPRS$Beta_B ))
inputData_blendFactors_shaPRS$Beta_A = as.numeric(as.character(inputData_blendFactors_shaPRS$Beta_A ))
inputData_blendFactors_shaPRS$b = as.numeric(as.character(inputData_blendFactors_shaPRS$b ))



# 2. Align PheB/B alleles
misalignedAlleleIndices = which( as.character(inputData_blendFactors_shaPRS$A1.x) != as.character(inputData_blendFactors_shaPRS$A1.y) ) # compare as character, as if we have non-SNPs with different alleles factors will break
inputData_blendFactors_shaPRS$Beta_A[misalignedAlleleIndices] = -inputData_blendFactors_shaPRS$Beta_A[misalignedAlleleIndices] # flip effects
if(length(misalignedAlleleIndices) > 0) message(paste0(length(misalignedAlleleIndices)), " misaligned allele(s) effects were reversed" )

cor(inputData_blendFactors_shaPRS$Beta_A,inputData_blendFactors_shaPRS$Beta_B)
cor(inputData_blendFactors_shaPRS$Beta_A,inputData_blendFactors_shaPRS$b)
cor(inputData_blendFactors_shaPRS$Beta_B,inputData_blendFactors_shaPRS$b)

# checking if the shaPRS A1/A2 are also aligned, it seems not, but this seems like an artefact, as the correlations are correct, so I am not flipping them
#misalignedAlleleIndices = which( as.character(inputData_blendFactors_shaPRS$A1.x) != as.character(inputData_blendFactors_shaPRS$A1) ) # compare as character, as if we have non-SNPs with different alleles factors will break
#if(length(misalignedAlleleIndices) > 0) message(paste0(length(misalignedAlleleIndices)), " misaligned allele(s) effects were reversed" )




diagData = subset(inputData_blendFactors_shaPRS, select=c("SNP", "lFDR", "Beta_A" , "SE_A"  ,  "Beta_B"  ,"SE_B", "b"   ,    "se" ))
names(diagData)[names(diagData) == 'b'] <- 'Beta_ShaPRS'
names(diagData)[names(diagData) == 'se'] <- 'SE_ShaPRS'


cor(diagData$Beta_A,diagData$Beta_B)
cor(diagData$Beta_A,diagData$Beta_ShaPRS)
cor(diagData$Beta_B,diagData$Beta_ShaPRS)

write.table(diagData,paste0(outputLoc,"_diagDat"), col.names = T, row.names = F, quote = F)


# MANHATTAN PLOT
inputData_blendFactors$P=inputData_blendFactors$Qval # add the adjusted Q vals as 'P', as that is col I would be plotting next
lfdr_2_blending = 1-inputData_blendFactors$lFDR
base_colour1 = 0.20
base_colour2 = 0.20
manhattanBaseColours = c(rgb(0,base_colour2,0,1),rgb(0,0,base_colour1,1) )



allIndices = inputData_blendFactors$BP

# need to offset the SNP BP indices, by the previous number of indices in all previous chromosomes
inputData_blendFactors$Offsets = 0
for (i in 1:21) { # we always set the next offset, so we dont loop til last Chrom
  message(i)
  CHR_SNPs = inputData_blendFactors[inputData_blendFactors$CHR == i,]
  maxBPCurrentChrom = max(CHR_SNPs$BP)
  currentOffset = CHR_SNPs$Offsets[1]
  nextOffset = currentOffset + maxBPCurrentChrom
  inputData_blendFactors[inputData_blendFactors$CHR == (i+1),9] = nextOffset
  
}
allIndices = allIndices +inputData_blendFactors$Offsets


# NAIndices = which(is.na(CHR_SNPs$BP) )
# NAIndices[1]
# 
# CHR_SNPs[214,]
# 
# 
# CHR_SNPs[205:215,]
# SNP CHR      BP    Beta_A      SE_A A_effectAllele  Beta_B    SE_B rs10001793   
# 
# 4 4341573 0.0684807 0.0371074              C 0.00728 0.00831
# B_effectAllele lFDR     Qval        P Offsets
# 214              T    1 3.969311 3.969311       0
# 
# 
# inputData_blendFactors[which(inputData_blendFactors$SNP== "rs10001793"),]
# 
# SNP CHR      BP    Beta_A      SE_A A_effectAllele  Beta_B    SE_B 214 rs10001793   4 4341573 0.0684807 0.0371074              C 0.00728 0.00831 B_effectAllele lFDR     Qval        P Offsets 214              T    1 3.969311 3.969311       0
# 

# get distribution
filename = paste(outputLoc, "_lFDR_hist.png" , sep="")
png(filename, width=640 , height=640, res =128);

hist(B12_blendFactors$lFDR, probability = F, col ="red", xlab = "lFDR", main ="", cex.axis=1.5, cex.lab=1.5)
dev.off()

# get distribution
# was getting 'tol' must be strictly positive and finite, error, this because IQR(x) as 0, so I capped tol
# https://girlincomputerscience.blogspot.com/2016/05/error-when-using-smoothspline.html
filename = paste(outputLoc, "_lFDR_B12.png" , sep="")
png(filename, width=640 , height=640, res =128);
plot(B12_blendFactors$lFDR, abs(B12_blendFactors$b), cex.axis=1.5, cex.lab=1.5, col =alpha("red", 0.3), xlab = "lFDR", ylab = "abs(SNP coef)", main ="")
smoothingSpline = smooth.spline(B12_blendFactors$lFDR, abs(B12_blendFactors$b), spar=1, tol=max(1e-6, 1e-6 * IQR(B12_blendFactors$lFDR) ) )
lines(smoothingSpline, col='blue', lwd=5)
dev.off()
#plot(B12_blendFactors$lFDR[1:5000], B12_blendFactors$b[1:5000], col =alpha("red", 0.3), xlab = "lFDR", ylab = "SNP coef", main ="")
#smoothingSpline = smooth.spline(B12_blendFactors$lFDR, abs(B12_blendFactors$b), spar=1)



#smoothingSpline = smooth.spline(B12_blendFactors$lFDR, abs(B12_blendFactors$b), tol=0.1)

#plot(B12_blendFactors$lFDR[1:10000], abs(B12_blendFactors$b[1:10000]), cex.axis=1.5, cex.lab=1.5, col =alpha("red", 0.3), xlab = "lFDR", ylab = "abs(SNP coef)", main ="")
#smoothingSpline = smooth.spline(B12_blendFactors$lFDR[1:10000], abs(B12_blendFactors$b[1:10000]), spar=1)
#lines(smoothingSpline, col='blue', lwd=5)


filename = paste(outputLoc, "_QAdjusted_manhattan.png" , sep="")
png(filename, width=1280 , height=640, res =128);
par(mar = c(5.1, 4.5, 4.1, 5.1))
plotTitle= ""
manhattan(inputData_blendFactors, col = manhattanBaseColours, cex.lab = 1.8, cex.axis = 1.5, cex.names = 1.2, cex.main = 2, main = plotTitle, suggestiveline=F, genomewideline=F, logp = F, ylab="Q adjusted")
par(new=TRUE)              # key: ask for new plot without erasing old  , mar = c(5.1, 4.1, 8, 5.1)
plot(allIndices,lfdr_2_blending, ylim= c(0,1),  xlab="", ylab="", yaxs = "i", type="n"  , axes=FALSE) # , type="n"  , xaxs = "i" , xaxt='n', ann=FALSE


axis(side=4, at = pretty(range(lfdr_2_blending)) , cex.axis = 1.5)
mtext("1-lFDR", side = 4, line = 3, cex.lab = 1.8, cex.axis = 1.5, cex = 1.5)
dev.off()

#inputData_blendFactors_old = inputData_blendFactors
#inputData_blendFactors = inputData_blendFactors[1:1000,]
#lfdr_2_blending = lfdr_2_blending[1:1000]
#length(lfdr_2_blending)
#allIndices = allIndices[1:1000]
