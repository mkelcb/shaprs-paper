# plots LD pattern
#if (!require("infotheo")) install.packages("infotheo")
#library("infotheo")

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 3) {stop("not enough Arguments received")} 
 
gwasLoc = args[1]
outputLoc = args[2]
outTitle = args[3]


# Debug wars
# outputLoc = "C:/softwares/Cluster/GIANT/miniPRS/debugfolder/"
#gwasLoc = "C:/softwares/Cluster/GIANT/miniPRS/debugfolder/all_SNPS_chiSQ_p"
# outTitle = "GWAS1"





#plotName = paste(outTitle,": D distribution of significant associations", sep="")


# 1. load phenos
gwasRes= read.table(gwasLoc, header = T) 
gwasRes = na.omit(gwasRes)


# Find variants
D_sd = sd(gwasRes$D)
outliers_indices = which( abs(gwasRes$D) > 2*D_sd )
outliers = gwasRes[outliers_indices,]

outliers_indices_medium = which( abs(gwasRes$D) > 1*D_sd & abs(gwasRes$D) < 2*D_sd )
outliers_medium = gwasRes[outliers_indices_medium,]


#filename = paste(outputLoc,plotName, "_LD_test.png" , sep="")
#png(filename, width=1280 , height=640, res =128);
filename = paste(outputLoc,outTitle, "_LD_test.pdf" , sep="")
pdf(filename, width=6.4 , height=6.4);
plot(gwasRes$D, xlab="SNPs", ylab="D", main=paste(outTitle, ": Significant SNPs LD support (sd(D):",round(D_sd,3), ")",sep="") )
points(outliers_indices,outliers$D,col=c("red"), pch =19)

legend("bottomright", legend=c("SNPs with D >2SD"),col=c("red"), text.col= c("red"), lty=1, lwd=3)


text(outliers_indices,outliers$D, labels=outliers$RSID, col="red", cex= 1.0, pos=2)





points(outliers_indices_medium,outliers_medium$D,col=c("green"), pch =19)

#model = lm( -log10(gwasRes$Pval)~ gwasRes$r.2)
#summary(model)


dev.off()

#gwasRes$V1[250:255]


filename = paste(outputLoc,outTitle, "_LD_test_hist2.pdf" , sep="")
pdf(filename, width=6.4 , height=6.4);
#hist(gwasRes$V2, xlab="D", main="Distribution of SNPs LD support", col ="red" , breaks="Scott")
hist(gwasRes$D, xlab="D", main=paste(outTitle, ": Distribution of SNPs LD support (sd(D):",round(D_sd,3), ")",sep=""), xlim=c(-5,5), col ="red" , breaks="FD")


#?hist

dev.off()

print(paste("saved plot to: ", filename ))



filename = paste(outputLoc,outTitle, "_LD_2D.pdf" , sep="")
pdf(filename, width=6.4 , height=6.4);

plot(gwasRes$D, -log10(gwasRes$p_value), xlab = "D", ylab = "-log10(p)" , main=paste(outTitle, ":Significant SNPs LD support (sd(D): ",round(D_sd,3), ")" , sep="") )
D_threshold = 1
outliers_indices_all = which( abs(gwasRes$D) > D_threshold*D_sd )
outliers_all = gwasRes[outliers_indices_all,]
points(outliers_all$D,-log10(outliers_all$p_value),col=c("red"), pch =19)
legend("topleft", legend=c( paste( "SNPs with abs(D) >",D_threshold,"SD", sep="" )  ),col=c("red"), text.col= c("red"), lty=1, lwd=3)
dev.off()
# medium bad: variants between > 1SD < 2SD


filename = paste(outputLoc,outTitle, "_LD_outliers" , sep="")
write.table(outliers_all$RSID, filename, sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
print(paste("saved outlier variants to: ", filename ))
# find example cases
#negMaxOutlier = outliers_all[which(outliers_all$D < -55),] # rs17035797
#posMaxOutlier  = outliers_all[which(outliers_all$D == max(outliers_all$D) ),] # rs11790250

#leastbadOutlier =  outliers_all[which(outliers_all$D == min( abs(outliers_all$D) ) ),] # rs6669582


# 
# #487409 - 110899
# raw= 1-0.0001
# conv=1-0.00005
# diff = raw - conv
# diff
# (sqrt(diff*diff) / raw )
# # 
# # 101981/179466
# # 1 - 0.0001
# #1018/179466
# 
# raw= 1-0.001
# conv=1-0.0005
# raw= 0.001
# conv=0.0005
# diff = raw - conv
# diff
# (sqrt(diff*diff) / raw )
# 
# 48735/144279
