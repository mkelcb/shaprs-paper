# produces a plot comparing the 2 baselines (subpheno and combined) against the advanced methods shaPRS and SMTPred


# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 5) {stop("not enough Arguments received")} 


plotName = args[1]
rG = args[2]
outputLoc = args[3]

options(error=traceback)

# load input files for each method 
combined= NULL
subpheno= NULL
shaPRS= NULL
SMTPred= NULL
xlabels=vector()
for(i in seq(from=4, to=length(args), by=4)){ # 4 as we also add the 'regular'
  baseLoc=args[i]
  print(paste0("baseLoc is: ", baseLoc))
  current_combined = read.table(paste0(baseLoc,"combined")  ,header=F)
  current_subpheno = read.table(paste0(baseLoc,"subpheno")  ,header=F)
  current_shaPRS = read.table(paste0(baseLoc,"shaPRS_meta")  ,header=F)
  current_SMTPred = read.table(paste0(baseLoc,"SMTPred")  ,header=F)

  # replace NAs with col mean
  current_combined[is.na(current_combined[,1]), 1] <- mean(current_combined[,1], na.rm = TRUE)
  current_subpheno[is.na(current_subpheno[,1]), 1] <- mean(current_subpheno[,1], na.rm = TRUE)
  current_shaPRS[is.na(current_shaPRS[,1]), 1] <- mean(current_shaPRS[,1], na.rm = TRUE)
  current_SMTPred[is.na(current_SMTPred[,1]), 1] <- mean(current_SMTPred[,1], na.rm = TRUE)
    
  if (is.null(combined)) {
    combined= current_combined
    subpheno= current_subpheno
    shaPRS= current_shaPRS
    SMTPred= current_SMTPred
  } else {
    combined= cbind( combined,current_combined )
    subpheno= cbind( subpheno, current_subpheno)
    shaPRS= cbind(shaPRS, current_shaPRS)
    SMTPred= cbind( SMTPred, current_SMTPred)
  }

  
  p_current = round(as.numeric(args[(i+1)]),2)
  shared_corr = round(as.numeric(args[(i+2)]),2)
  print(paste0("p_current: ",p_current, " / shared_corr: ", shared_corr, " | baseLoc: ", baseLoc))
  xlabels = c(xlabels, paste0("p:",p_current,"/r:",shared_corr) )
}

#############################################
# LOCAL DEBUG
# args=vector()
# args =c(args,"#causals_1000_rG_0.1_A50_B50_size_half")
# args = c(args,"0.1")
# args = c(args,"C:/softwares/Cluster/0shaPRS/debug/#causals_1000_rG_0.1_A50_B50_size_half")
# 
# args = c(args,"C:/softwares/Cluster/0shaPRS/debug/1000/10/0.1_1.0/A50_B50/size_half/")
# args = c(args,"0.1")
# args = c(args,"1.0")
# args = c(args,"C:/softwares/Cluster/0shaPRS/debug/1000/10/0.55_0.1818182/A50_B50/size_half/")
# args = c(args,"0.55")
# args = c(args,"0.1818182")
# args = c(args, "C:/softwares/Cluster/0shaPRS/debug/1000/10/1.0_0.1/A50_B50/size_half/")
# args = c(args,"1.0")
# args = c(args,"0.1")
# 



########################################

combinedIndex=1
subphenoIndex=2
shaPRSIndex=3
SMTPredIndex=4
boxNames=c("combined", "subpheno", "shaPRS", "SMTPred")
methodCols = c("#11d7d5", "#166938", "#c5b000", "#8e01c2")



xaxisname=paste0("Genetic heterogeneity (rG: ",rG,")")

# Debug vars
# combinedLoc="C:/softwares/Cluster/GIANT/miniPRS/shaPRS/combined"
# subphenoLoc="C:/softwares/Cluster/GIANT/miniPRS/shaPRS/subpheno"
# shaprsLoc="C:/softwares/Cluster/GIANT/miniPRS/shaPRS/shaPRS"
# plotName = "effect of rG on shaPRS"
# outputLoc = "C:/softwares/Cluster/GIANT/miniPRS/shaPRS/out"


plotName = gsub("_", " ", plotName) #remove underscores

########################################################



x_axis_marks=seq(1:ncol(combined))

combined_means= colMeans(combined, na.rm = T)
subpheno_means= colMeans(subpheno, na.rm = T)
shaprs_means= colMeans(shaPRS, na.rm = T)
smtpred_means= colMeans(SMTPred, na.rm = T)

combined_sd = sapply(combined, sd, na.rm = TRUE)
subpheno_sd = sapply(subpheno, sd, na.rm = TRUE)
shaprs_sd = sapply(shaPRS, sd, na.rm = TRUE)
smtpred_sd = sapply(SMTPred, sd, na.rm = TRUE)

sampleSize_of_mean = nrow(combined)


# plot the 95% confidence interval: SE = SD/sqrt(n) , and then mean +- 1.96xSE
subpheno_se = subpheno_sd / sqrt(sampleSize_of_mean)  # / sqrt(samplesizeNums) # the sample size does not mean the sample size of the NN model, but the number of results
combined_se = combined_sd / sqrt(sampleSize_of_mean)  
shaprs_se = shaprs_sd / sqrt(sampleSize_of_mean)  
smtpred_se = smtpred_sd / sqrt(sampleSize_of_mean)  

subpheno_bound = 1.96 *subpheno_se
combined_bound = 1.96 *combined_se
shaprs_bound = 1.96 *shaprs_se
smtpred_bound = 1.96 *smtpred_se
#print(paste0("smtpred_bound: ", smtpred_bound))
#print(paste0("combined_bound: ", combined_bound))
#print(paste0("shaprs_bound: ", shaprs_bound))
#print(paste0("subpheno_bound: ", subpheno_bound))

# y limit
maxSD=max(subpheno_se,combined_se,shaprs_se,smtpred_se)
overallMax= (max(subpheno_means,combined_means,shaprs_means,smtpred_means) + maxSD ) *1.1

# print(paste("subpheno_means is", subpheno_means))
# print(paste("combined_means is", combined_means))
# print(paste("shaprs_means is", shaprs_means))
# print(paste("overallMax is", overallMax))

# perform paired 2 sample t-tests and get their % differences
for ( i in 1:ncol(combined)) {    # https://www.mathsisfun.com/percentage-difference.html

  
  print(paste("at heterogeneity",xlabels[i] ,":" ))
  percDiff = round( (combined_means[i] - subpheno_means[i]) / ( (combined_means[i] + subpheno_means[i])/2 ) * 100)
  print(paste(boxNames[combinedIndex] ,"v", boxNames[subphenoIndex], "| t-test:", t.test(combined[,i], subpheno[,i], paired = T)$p.value , "/ diff:", percDiff,"%" ))
  
  percDiff = round( (combined_means[i] - shaprs_means[i]) / ( (combined_means[i] + shaprs_means[i])/2 ) * 100)
  print(paste(boxNames[combinedIndex] ,"v", boxNames[shaPRSIndex], "| t-test:", t.test(combined[,i], shaPRS[,i], paired = T)$p.value , "/ diff:", percDiff,"%" ))
  
  percDiff = round( (subpheno_means[i] - shaprs_means[i]) / ( (subpheno_means[i] + shaprs_means[i])/2 ) * 100)
  print(paste(boxNames[subphenoIndex] ,"v", boxNames[shaPRSIndex], "| t-test:", t.test(subpheno[,i], shaPRS[,i], paired = T)$p.value , "/ diff:", percDiff,"%" ))
  
  # shaPRS vs SMTPred
  percDiff = round( (smtpred_means[i] - shaprs_means[i]) / ( (smtpred_means[i] + shaprs_means[i])/2 ) * 100)
  print(paste(boxNames[SMTPredIndex] ,"v", boxNames[shaPRSIndex], "| t-test:", t.test(SMTPred[,i], shaPRS[,i], paired = T)$p.value , "/ diff:", percDiff,"%" ))
  
  print("______________________________")
}







#filen =paste(outputLoc,".pdf", sep="" )
#pdf(filen, width=6.4 , height=6.4);
filen =paste(outputLoc,".png", sep="" )
png(filen, width=680  , height=200*3);
plotName="" # plotname disabled for publication
plot(combined_means~x_axis_marks, xlab=xaxisname, type="o", pch=19,ylab = "",xaxt="n", col = methodCols[combinedIndex],lwd = 4, cex =2.0, ylim=c(0,overallMax), main=plotName, cex.lab = 1.7, cex.axis = 1.5, cex.main = 2) # , cex.names = 1.2
arrows(x_axis_marks,combined_means, x_axis_marks,combined_means+ combined_bound, lwd=2, length=0.05, angle=90, col = methodCols[combinedIndex])
arrows(x_axis_marks,combined_means, x_axis_marks,combined_means- combined_bound, lwd=2, length=0.05, angle=90, col = methodCols[combinedIndex])

par(new=T)
plot(subpheno_means~x_axis_marks, type="o", ann=F, pch=19,xaxt="n",yaxt="n", col =methodCols[subphenoIndex],lwd = 4, cex =2.0, ylim=c(0,overallMax))
arrows(x_axis_marks,subpheno_means, x_axis_marks,subpheno_means+ subpheno_bound, lwd=2, length=0.05, angle=90, col =methodCols[subphenoIndex])
arrows(x_axis_marks,subpheno_means, x_axis_marks,subpheno_means- subpheno_bound, lwd=2, length=0.05, angle=90, col =methodCols[subphenoIndex])

par(new=T)
plot(shaprs_means~x_axis_marks, type="o", ann=F, pch=19,xaxt="n",yaxt="n", col = methodCols[shaPRSIndex],lwd = 4, cex =2.0, ylim=c(0,overallMax))
arrows(x_axis_marks,shaprs_means, x_axis_marks,shaprs_means+ shaprs_bound, lwd=2, length=0.05, angle=90, col = methodCols[shaPRSIndex])
arrows(x_axis_marks,shaprs_means, x_axis_marks,shaprs_means- shaprs_bound, lwd=2, length=0.05, angle=90, col = methodCols[shaPRSIndex])


par(new=T)
plot(smtpred_means~x_axis_marks, type="o", ann=F, pch=19,xaxt="n",yaxt="n", col = methodCols[SMTPredIndex],lwd = 4, cex =2.0, ylim=c(0,overallMax))
arrows(x_axis_marks,smtpred_means, x_axis_marks,smtpred_means+ smtpred_bound, lwd=2, length=0.05, angle=90, col = methodCols[SMTPredIndex])
arrows(x_axis_marks,smtpred_means, x_axis_marks,smtpred_means- smtpred_bound, lwd=2, length=0.05, angle=90, col = methodCols[SMTPredIndex])




legend("bottomright", inset=0, legend=boxNames, horiz=TRUE,
       col= c(methodCols[combinedIndex],methodCols[subphenoIndex],methodCols[shaPRSIndex]), text.col=c(methodCols[combinedIndex],methodCols[subphenoIndex],methodCols[shaPRSIndex], methodCols[SMTPredIndex]), fill=c(methodCols[combinedIndex],methodCols[subphenoIndex],methodCols[shaPRSIndex], methodCols[SMTPredIndex]), cex=1.6, pt.cex = 1,box.lty=0, x.intersp=0.25) # ,inset=-0.01


axis(1,at=x_axis_marks,labels=xlabels, cex.axis = 2,mgp=c(3,1.3,0))
title(ylab="Predictive power (r^2)", line=2.6, cex.lab=1.7)

dev.off()

print(paste("written plot to",filen))

 
 

