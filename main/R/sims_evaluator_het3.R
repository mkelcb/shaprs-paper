# create heatmap for simulation results
library(pheatmap)

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 2) {stop("not enough Arguments received")} 
 

outputLoc= args[1]
allResultsLoc= args[2]
allPredictorsLoc= args[3]
# 
# # Debug vars Orig
#allResultsLoc = "C:/0Datasets/shaPRS/0shaPRS/debug/all_results"
#allPredictorsLoc= "C:/0Datasets/shaPRS/0shaPRS/debug/all_predictors"
#outputLoc =   "C:/0Datasets/shaPRS/0shaPRS/debug/sims_eval"

# UKBB
#allResultsLoc = "C:/0Datasets/shaPRS/0shaPRS/debug/ukbb/all_results"
#allPredictorsLoc= "C:/0Datasets/shaPRS/0shaPRS/debug/ukbb/all_predictors"
#outputLoc =   "C:/0Datasets/shaPRS/0shaPRS/debug/ukbb/sims_eval_ukbb"

allResultsLoc = "C:/0Datasets/shaPRS/0shaPRS/0ukbb/results/all_results"
allPredictorsLoc= "C:/0Datasets/shaPRS/0shaPRS/0ukbb/results/all_predictors"
outputLoc =   "C:/0Datasets/shaPRS/0shaPRS/0ukbb/results/sims_eval_ukbb"

################################
padTo = function(x, digits, padder = " ") { # pads a string to a certain length
  x = as.character(x)
  while (nchar(x) < digits ){
    x = paste0(x,padder)
  }
  return(x)
}  

padTo_dec = function(x, digits) { # pads a string to a certain length
  x = as.character(x)
  while (nchar(x) < digits ){
    if (nchar(x) == 1) {
      x = paste0(x,".")
    } else { x = paste0(x,"0") }
  }
  return(x)
}  
################################

allResults = allResults_orig = read.csv(allResultsLoc, header = T)
allPredictors = read.csv(allPredictorsLoc, header = T)
modelNames = names(allResults)

# go through each method
for (i in 1:length(modelNames)) { 
  method = modelNames[i]
  
  # perform regression analyses for each method
  lmodel =  lm(allResults[,i] ~ ., data = allPredictors)
  sink(paste0(outputLoc,"_",method,"txt"))
  print(summary(lmodel))
  sink()  # returns output to the console
  
}

for (i in 1:length(modelNames)) { 
  method = modelNames[i]
  
  # perform regression analyses for each method
  lmodel =  lm(allResults[,i] ~ .^2, data = allPredictors)
  sink(paste0(outputLoc,"_int_",method,"txt"))
  print(summary(lmodel))
  sink()  # returns output to the console
  
}

# go through each method but only look at the 1000 SNP runs
only1000SNPs = which(allPredictors$X.causals =="1000")
for (i in 1:length(modelNames)) { 
  method = modelNames[i]
  
  # perform regression analyses for each method
  lmodel =  lm(allResults[only1000SNPs,i] ~ p_current + shared_corr, data = allPredictors[only1000SNPs,])
  sink(paste0(outputLoc,"_1000_",method,"txt"))
  print(summary(lmodel))
  sink()  # returns output to the console
  
}

for (i in 1:length(modelNames)) { 
  method = modelNames[i]
  
  # perform regression analyses for each method
  lmodel =  lm(allResults[,i] ~ p_current + shared_corr, data = allPredictors[,])
  sink(paste0(outputLoc,"_all_",method,"txt"))
  print(summary(lmodel))
  sink()  # returns output to the console
  
}
sampleSizeAsN = allPredictors$sample_size
sampleSizeAsN[which(sampleSizeAsN == "_half")] = 8000
sampleSizeAsN[which(sampleSizeAsN == "_full")] = 16000
sampleSizeAsN[which(sampleSizeAsN == "_double")] = 32000
sampleSizeAsN = as.numeric(sampleSizeAsN)

######################################################
# I) Plot all of the data


# overall results
median(allResults$combined) # 0.05962532
median(allResults$subpheno) # 0.07367951
median(allResults$shaPRS) # 0.1139283
median(allResults$SMTPred) # 0.07585284

mean(allResults$combined) #  0.07426966
mean(allResults$subpheno) #  0.09124099
mean(allResults$shaPRS) #   0.1189231
mean(allResults$SMTPred) #   0.09280919

# Combined seems a bit low, check when it is better than subpheno
combined_better_than_subpheno_indices = which(allResults$combined > allResults$subpheno)
combined_better_than_subpheno = allResults[combined_better_than_subpheno_indices,]

SMTPred_better_than_shapRS_indices = which(allResults$SMTPred > allResults$shaPRS)
SMTPred_better_than_shapRS = allResults[SMTPred_better_than_shapRS_indices,]

# create heatmaps
i=1
# create rownames
coln = c()
rGs = c()
ps= c()
cors= c()
splits= c()
sample_sizes= c()
for (i in 1:length(allPredictors[,1])) { 
  print(i)
  
  sample_size = gsub("_", "", allPredictors[i,7]) #remove underscores
  rG= padTo_dec( allPredictors[i,2], 4)
  p= padTo_dec( round(allPredictors[i,3], 2) ,4)
  corre = padTo_dec( round(allPredictors[i,4], 2) ,4)
  
  annotations = paste0("rG:",rG, " p:",p, " cor:",corre, " ",allPredictors[i,5], " split:",allPredictors[i,6], " size:",sample_size)
  
  coln = c(coln,annotations)
  
  splits = c(splits,allPredictors[i,6] )
  ps = c(ps,p )
  cors = c(cors,corre )
  sample_sizes = c(sample_sizes,sample_size )
  rGs = c(rGs,paste0(rG," ") ) # add an extra space for padding, otherwise legend will be cut off as pheatmap is shit
}
rownames(allResults) = coln
# want to sort the data rows by rG
orderByRG = order( as.numeric(rGs))
allResults = allResults[orderByRG,]

splits = splits[orderByRG]
rGs = rGs[orderByRG]
ps = ps[orderByRG]
cors = cors[orderByRG]
sample_sizes = sample_sizes[orderByRG]

rGs_DF <- data.frame( rGs,as.numeric(ps),as.numeric(cors),sample_sizes, splits,row.names=rownames(allResults)) # need to match the rownames to the data for the annotation to work
colnames(rGs_DF) <- c("rG", "p", "cor","N", "split") # this is the header for the annotation


# Visualise ALL, to be used in the Supplementary
pheatmap(allResults, main = "shaPRS - all simulations", filename=paste(outputLoc,".png", sep="" ),annotation_row=rGs_DF, height=25, width=10 , cex=1 ,cluster_rows=F, cluster_cols=F)




# II) create a plot of a subset of rG=0.5, full
# splits = splits[orderByRG]
# rGs = rGs[orderByRG]
# ps = ps[orderByRG]
# cors = cors[orderByRG]
#sample_sizes = sample_sizes[orderByRG]

# function to separate regular/extra results:
filterOutByTerm_all = function(allResults,splits,ps,cors,filterTerm = "rG:0.50") {
  subsetResults = allResults

  ps_subset = ps
  cors_subset = cors
  splits_subset = splits
  indices_kept = c()
  i=1
  for (i in 1:nrow(subsetResults) ) { 
    rowname = rownames(subsetResults[i,])
    
    # check if rowname includes extra/regular, 
    if ( grepl( filterTerm, rowname, fixed = TRUE) ) {
      # if yes, we replace it with nothing, and keep it
      rowname_new = gsub(filterTerm, "",rowname) #remove underscores
      rownames(subsetResults)[rownames(subsetResults) == rowname] <- rowname_new
      indices_kept = c(indices_kept, i)
      
    } # discard it otherwise
  }
  subsetResults = subsetResults[indices_kept,]

  ps_subset = ps_subset[indices_kept]
  cors_subset = cors_subset[indices_kept]
  splits_subset = splits_subset[indices_kept]
  
  results = NULL
  results$subsetResults = subsetResults
  results$ps_subset = ps_subset
  results$cors_subset = cors_subset
  results$splits_subset = splits_subset
  
  return(results)
}


# Filter to keep the main interesting scenarios, rG 0.5, regular, full
results_RG05 = filterOutByTerm_all(allResults,splits,ps,cors,filterTerm = "rG:0.50")
results_regular = filterOutByTerm_all(results_RG05$subsetResults,results_RG05$splits_subset,results_RG05$ps_subset,results_RG05$cors_subset,filterTerm = "regular")
results_full = filterOutByTerm_all(results_regular$subsetResults,results_regular$splits_subset,results_regular$ps_subset,results_regular$cors_subset,filterTerm = "size:full")


subsetResults = results_full$subsetResults

subset_DF <- data.frame( results_full$ps_subset, results_full$cors_subset,results_full$splits_subset,row.names=rownames(subsetResults)) # ,row.names=rownames(subsetResults) # need to match the rownames to the data for the annotation to work
colnames(subset_DF) <- c("p","cor","split") # this is the header for the annotation

#plotName = "shaPRS - rG:0.5, n:full, no extra" # no plotname for final publication
plotName =""
pheatmap(subsetResults, main = plotName , filename=paste(outputLoc,"_subset.png", sep="" ),annotation_row=subset_DF, show_rownames = F, height=5, width=5 , cex=1 ,cluster_rows=F, cluster_cols=F)


#######################################################################################

# Quantify simulation results

# gets percentage difference between 2 vectors, positive results indicate first vector is better, negative result that the second is better
percDifference = function (myData) {
  mean1 = mean(myData[,1])
  mean2 = mean(myData[,2])
  percDiff = round( (mean1 - mean2) / ( (mean1 + mean2)/2 ) * 100)
  return(percDiff)
}


# find the largest difference between subpheno and shaPRS, in the full set
results_full_all = filterOutByTerm_all(allResults,splits,ps,cors,filterTerm = "size:full")$subsetResults
results_half_all = filterOutByTerm_all(allResults,splits,ps,cors,filterTerm = "size:half")$subsetResults
results_double_all = filterOutByTerm_all(allResults,splits,ps,cors,filterTerm = "size:double")$subsetResults


getDifferences = function(current_res, otherMethodIndex = 2) # 1 = combined, 2 = subpheno, 4= SMTPred
{
 # print(paste0("otherMethodIndex: ", otherMethodIndex))
  # get the perc of method relative to shaPRS
  othermethodName = colnames(current_res)[otherMethodIndex]
  otherMethod_AsPercOfShaPRS = current_res[,otherMethodIndex] /current_res$shaPRS

  largestDiff_index = which(otherMethod_AsPercOfShaPRS == min(otherMethod_AsPercOfShaPRS)) # when shaPRS was best
  smallestDiff_index = which(otherMethod_AsPercOfShaPRS == max(otherMethod_AsPercOfShaPRS)) # when shaPRS was worst

  pDiff_best = percDifference( cbind(current_res$shaPRS[largestDiff_index], current_res[largestDiff_index,otherMethodIndex]) ) # -127
  print(paste0("best shaPRS is better by ", pDiff_best, "% than ",othermethodName, " at scenario: " ))
  print(current_res[largestDiff_index,])
  
  
  pDiff_worst = percDifference( cbind(current_res$shaPRS[smallestDiff_index], current_res[smallestDiff_index,otherMethodIndex]) ) # -127
  print(paste0("worst shaPRS is better by ", pDiff_worst, "% than ",othermethodName, " at scenario: " ))
  print(current_res[smallestDiff_index,])
}

# find out the % of scenarios where shaPRS was the best
totalSimScenarios = nrow(allResults)
bestResults = apply(allResults, MARGIN = 1, max) # find max of each row, IE t he best method
shaPRS_bestIndices = which(allResults$shaPRS >= bestResults)
shAprsNumBest = length(shaPRS_bestIndices)

percShaPRsBest = shAprsNumBest / totalSimScenarios

print(paste0("shaPRS was best in ",round(percShaPRsBest * 100),"% of the simulation scenarios"))


getDifferences(allResults,1)
[1] "best shaPRS is better by 154% than combined at scenario: "
combined   subpheno     shaPRS    SMTPred
rG:0.10 p:0.10 cor:1.00 extra split:20:80 size:half 0.004496213 0.02696653 0.03489879 0.02410484
[1] "worst shaPRS is better by 0% than combined at scenario: "
combined    subpheno     shaPRS     SMTPred
rG:0.25 p:0.25 cor:1.00 regular split:20:80 size:half 0.02339498 0.004262317 0.02337319 0.007495933


getDifferences(allResults,2)
[1] "best shaPRS is better by 178% than subpheno at scenario: "
combined    subpheno     shaPRS    SMTPred
rG:0.50 p:0.50 cor:1.00 regular split:20:80 size:half 0.04253765 0.002508403 0.04290523 0.01456814
[1] "worst shaPRS is better by -15% than subpheno at scenario: "
combined   subpheno     shaPRS    SMTPred
rG:0.10 p:1.00 cor:0.10 regular split:50:50 size:half 0.01935849 0.03283975 0.02838844 0.03451448


getDifferences(allResults,4)
[1] "best shaPRS is better by 103% than SMTPred at scenario: "
combined    subpheno     shaPRS     SMTPred
rG:0.25 p:0.25 cor:1.00 regular split:20:80 size:half 0.02339498 0.004262317 0.02337319 0.007495933
[1] "worst shaPRS is better by -19% than SMTPred at scenario: "
combined   subpheno     shaPRS    SMTPred
rG:0.10 p:1.00 cor:0.10 regular split:50:50 size:half 0.01935849 0.03283975 0.02838844 0.03451448




getDifferences(results_full_all,1)
getDifferences(results_full_all,2)
getDifferences(results_full_all,4)

















#########################################################


# function to separate regular/extra results:
filterOutByTerm = function(allResults, rGs,heterogeneity_term = " regular") {
  allResults_heterogeneity = allResults
  rGs_heterogeneity = rGs
  indices_kept = c()
  
  for (i in 1:nrow(allResults_heterogeneity) ) { 
    rowname = rownames(allResults_heterogeneity[i,])
    
    # check if rowname includes extra/regular, 
    if ( grepl( heterogeneity_term, rowname, fixed = TRUE) ) {
      # if yes, we replace it with nothing, and keep it
      rowname_new = gsub(heterogeneity_term, "",rowname) #remove underscores
      rownames(allResults_heterogeneity)[rownames(allResults_heterogeneity) == rowname] <- rowname_new
      indices_kept = c(indices_kept, i)
      
    } # discard it otherwise
  }
  allResults_heterogeneity = allResults_heterogeneity[indices_kept,]
  rGs_heterogeneity = rGs_heterogeneity[indices_kept]
  
  results = NULL
  results$allResults_heterogeneity = allResults_heterogeneity
  results$rGs_heterogeneity = rGs_heterogeneity
  
  return(results)
}

# get the regular ones apart
results = filterOutByTerm(allResults, rGs,heterogeneity_term = " regular")
allResults_heterogeneity_regular = allResults_heterogeneity = results$allResults_heterogeneity
rGs_heterogeneity = results$rGs_heterogeneity
rGs_heterogeneity_DF <- data.frame( rGs_heterogeneity,row.names=rownames(allResults_heterogeneity)) # need to match the rownames to the data for the annotation to work
colnames(rGs_heterogeneity_DF) <- c("rG") # this is the header for the annotation

pheatmap(allResults_heterogeneity, main = "shaPRS - extra 5 large effect SNPs", filename=paste(outputLoc,"_reg.png", sep="" ),annotation_row=rGs_heterogeneity_DF, height=25, width=10 , cex=1 ,cluster_rows=F, cluster_cols=F)


# get the extra ones apart
results = filterOutByTerm(allResults, rGs,heterogeneity_term = " extra")
allResults_heterogeneity_extra = allResults_heterogeneity = results$allResults_heterogeneity
rGs_heterogeneity = results$rGs_heterogeneity
rGs_heterogeneity_DF <- data.frame( rGs_heterogeneity,row.names=rownames(allResults_heterogeneity)) # need to match the rownames to the data for the annotation to work
colnames(rGs_heterogeneity_DF) <- c("rG") # this is the header for the annotation

pheatmap(allResults_heterogeneity, main = "shaPRS - regular", filename=paste(outputLoc,"_extra.png", sep="" ),annotation_row=rGs_heterogeneity_DF, height=25, width=10 , cex=1 ,cluster_rows=F, cluster_cols=F)

View(allResults_heterogeneity)
# see what effect adding the extra 5 large effect SNPs have





median(allResults_heterogeneity_regular$shaPRS)
median(allResults_heterogeneity_extra$shaPRS)

combined_index=1
subpheno_index=2
shaPRS_index=3
SMTPred_index=4

boxNames=c("combined", "subpheno", "shaPRS", "SMTPred")
methodCols = c("#11d7d5", "#166938", "#c5b000", "#8e01c2")
overalXmax = max(allResults_heterogeneity_regular)
overalYmax = max(allResults_heterogeneity_extra)

# find out the difference between the performance of the methods of regular vs extra
legendNames = c()
for (i in 1:length(boxNames) ) { 
  regularrSQ = median(allResults_heterogeneity_regular[,i])
  extrarSQ = median(allResults_heterogeneity_extra[,i])
  averag = (regularrSQ + extrarSQ) /2
  percDiff= ( (extrarSQ - regularrSQ ) / averag ) * 100
  legendNames = c(legendNames, paste0(boxNames[i], " (", round(percDiff), "%)") )
}


plot(allResults_heterogeneity_regular$shaPRS,allResults_heterogeneity_extra$shaPRS, col =methodCols[shaPRS_index], xlab ="regular architecture", ylab ="5 large effect SNPs", main="effect of large effect SNPs", xlim = c(0,overalXmax), ylim =c(0,overalYmax))
par(new=TRUE)              # key: ask for new plot without erasing old  , mar = c(5.1, 4.1, 8, 5.1)
plot(allResults_heterogeneity_regular$combined,allResults_heterogeneity_extra$combined, col =methodCols[combined_index], xlab ="", ylab ="",xlim = c(0,overalXmax), ylim =c(0,overalYmax) , axes=FALSE)

par(new=TRUE) 
plot(allResults_heterogeneity_regular$subpheno,allResults_heterogeneity_extra$subpheno, col =methodCols[subpheno_index], xlab ="", ylab ="",xlim = c(0,overalXmax), ylim =c(0,overalYmax) , axes=FALSE)

par(new=TRUE) 
plot(allResults_heterogeneity_regular$SMTPred,allResults_heterogeneity_extra$SMTPred, col =methodCols[SMTPred_index], xlab ="", ylab ="",xlim = c(0,overalXmax), ylim =c(0,overalYmax) , axes=FALSE)
legend(x="bottomright", legend =legendNames, col=methodCols,pch=c(1,1))
##########################################################

# plot heterogeneity vs performance:
# X: low, medium, high (IE when when p is the lowest  )


high_heterog_indices = c()
medium_heterog_indices = c()
low_heterog_indices = c()
for (i in seq(1,nrow(allPredictors), by=6) ) { 
  high_heterog_indices = c(high_heterog_indices,i, i+1) #high heterogeneity, every 1st and 2nd result
  medium_heterog_indices = c(medium_heterog_indices,i+2, i+3) #medium heterogeneity, every 3rd and 4th result
  low_heterog_indices = c(low_heterog_indices,i+4, i+5) #low heterogeneity, every 5th and 6th result
}

allResults_highHet = allResults_orig[high_heterog_indices,]
allPredictors_highHet = allPredictors[high_heterog_indices,]

allResults_mediumHet = allResults_orig[medium_heterog_indices,]
allPredictors_mediumHet = allPredictors[medium_heterog_indices,]

allResults_lowHet = allResults_orig[low_heterog_indices,]
allPredictors_lowHet = allPredictors[low_heterog_indices,]

# write out time series data
timeSeriesLoc ="C:/0Datasets/shaPRS/0shaPRS/0ukbb/timeseries/"

combined = cbind.data.frame(allResults_lowHet$combined,allResults_mediumHet$combined,allResults_highHet$combined)
subpheno = cbind.data.frame(allResults_lowHet$subpheno,allResults_mediumHet$subpheno,allResults_highHet$subpheno)
shaPRS = cbind.data.frame(allResults_lowHet$shaPRS,allResults_mediumHet$shaPRS,allResults_highHet$shaPRS)
SMTPred = cbind.data.frame(allResults_lowHet$SMTPred,allResults_mediumHet$SMTPred,allResults_highHet$SMTPred)

write.table(combined,paste0(timeSeriesLoc,"combined"), col.names = F, row.names = F, quote = F)
write.table(subpheno,paste0(timeSeriesLoc,"subpheno"), col.names = F, row.names = F, quote = F)
write.table(shaPRS,paste0(timeSeriesLoc,"shaPRS"), col.names = F, row.names = F, quote = F)
write.table(SMTPred,paste0(timeSeriesLoc,"SMTPred"), col.names = F, row.names = F, quote = F)

allResults_orig
allPredictors


######################################################

# II) pulling out individual subsets

half_indices = which(allPredictors$sample_size == "_half")
full_indices = which(allPredictors$sample_size == "_full")
double_indices = which(allPredictors$sample_size == "_double")
titleExtra ="(ALL causals)"

half_indices = which(allPredictors$sample_size == "_half" & allPredictors$X.causals == 1000)
full_indices = which(allPredictors$sample_size == "_full"& allPredictors$X.causals == 1000)
double_indices = which(allPredictors$sample_size == "_double"& allPredictors$X.causals == 1000)
titleExtra ="(1000 causals)"



# Find when SMPTred is best
current_indices = which(allPredictors$sample_size == "_double"& allPredictors$X.causals == 1000 & allPredictors$splitA.splitB == "50:50" )
half_indices = full_indices = double_indices = current_indices
mainOverride ="(UKBB) SMTPred best, sample size double, split 50-50, Causals 1000"
#titleExtra="SMTPred"
titleExtra="UKBB_SMTPred"

# find when shaPRS is best
current_indices = which(allPredictors$sample_size == "_full"& allPredictors$X.causals == 1000 & allPredictors$splitA.splitB == "50:50" )
half_indices = full_indices = double_indices = current_indices
mainOverride ="(UKBB) shaPRS best, sample size full, split 50-50, Causals 1000"
titleExtra="shaPRS"

i=1

plotSimRres = function (half_indices,full_indices, double_indices, titleExtra, mainOverride ="") {
  boxNames=c("combined", "subpheno", "shaPRS", "SMTPred")
  methodCols = c("#11d7d5", "#166938", "#c5b000", "#8e01c2")
  sample_size_labels = c("half", "full", "double")
  sample_size = list(half_indices, full_indices, double_indices)
  
  
  for (i in 1:length(sample_size)) {
    current_indices = sample_size[[i]]
    currentData = cbind.data.frame(allResults$combined[current_indices],allResults$subpheno[current_indices],allResults$shaPRS[current_indices],allResults$SMTPred[current_indices])
    
    resultLabels=c()
    for (j in 1:length(boxNames)) {
    meanMethod = round(mean(currentData[,j]),3)
    resultLabels = c(resultLabels, paste0(boxNames[j],"\n(",meanMethod,")"))
   
    }
   
    colnames(currentData) = resultLabels
    
    filename=paste(outputLoc,"_",titleExtra,"_",sample_size_labels[i],".png", sep="" )
    png(filename, width=600  , height=400);
    par(mgp=c(3,2,0))
    if(mainOverride == "") {mainTitle = paste0("Methods at ",sample_size_labels[i]," sample size " , titleExtra)} else { 
      mainTitle= mainOverride}
    boxplot(currentData, col =methodCols, xlab="method", main =mainTitle , ylab ="r^2" )
    dev.off()
    
    print(paste("written plot to",filename))
  }
}



#######################################################################################


# find when shaPRS is best
current_indices = which(allPredictors$sample_size == "_full"& allPredictors$X.causals == 1000 & allPredictors$splitA.splitB == "50:50" )
half_indices = full_indices = double_indices = current_indices
mainOverride ="shaPRS best, sample size full, split 50-50, Causals 1000"
titleExtra="shaPRS"

# compare shaPRS against itself, in the 'tight' and 'even' rG cases
tight_rG = current_indices[c(1,4,7)]
even_rG = current_indices[c(3,6,9)]
currentData = cbind.data.frame(allResults$shaPRS[tight_rG],allResults$shaPRS[even_rG])
colnames(currentData) = c("shaPRS_tightRG", "shaPRS_evenRG")
filename=paste(outputLoc,"_shaPRS_tight_vs_even_",".png", sep="" )
png(filename, width=600  , height=400);
boxplot(currentData, col =methodCols[1:2], xlab="genetic architecture", main ="shaPRS, when rG is made up of tight or even distribution of SNPs" , ylab ="r^2" )
dev.off()



