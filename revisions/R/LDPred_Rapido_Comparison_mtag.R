args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)

inputDataLoc = args[1]

# DebugVars
#inputDataLoc="C:/0Datasets/shaPRS/ldPredComparison/"

# first col: RapidoPG, 2nd col LDpred2
subpheno="subpheno_orig_LDpred"
combined="combined_orig_LDpred"
shaprs="shaPRS_meta_orig_LDpred"
smtpred="SMTPred_orig_LDpred"
mtag="MTAG_orig_LDpred"

# gets percentage difference between 2 vectors, positive results indicate first vector is better, negative result that the second is better
#myData=combinedResults
percDifference = function (myData) {
  mean1 = mean(myData[,1])
  mean2 = mean(myData[,2])
  percDiff = round( (mean1 - mean2) / ( (mean1 + mean2)/2 ) * 100)
  return(percDiff)
}


subphenoResults = read.table(paste0(inputDataLoc,subpheno), header = F)
combinedResults = read.table(paste0(inputDataLoc,combined), header = F)
shaprsResults = read.table(paste0(inputDataLoc,shaprs), header = F)
smtpredResults = read.table(paste0(inputDataLoc,smtpred), header = F)
mtagResults = read.table(paste0(inputDataLoc,mtag), header = F)


# concat all data to get overall correlations
allData = rbind(subphenoResults,combinedResults,shaprsResults,smtpredResults,mtagResults)

 #  0.8285714
#cor(allData$V1,allData$V2) # 0.8417679
TestRes = t.test(allData$V1,allData$V2, paired = T) # p-value = 0.04034

print(paste0("Spearman corr: ", cor(allData$V1,allData$V2, method = "spearman"),  " / p: ", TestRes$p.value ))

# get absolute difference, this is is best calculated per method, as their baselines will be different
subpheno_meanDiff = percDifference(subphenoResults)
combined_meanDiff = percDifference(combinedResults)
shaprs_meanDiff = percDifference(shaprsResults)
smtpred_meanDiff =  percDifference(smtpredResults)
mtagResults_meanDiff =  percDifference(mtagResults)


all_diff = c(subpheno_meanDiff,combined_meanDiff,shaprs_meanDiff, smtpred_meanDiff, mtagResults_meanDiff)
print(paste0("method means")) # 6
print(paste0("subpheno Rapido/LDpred2 : ", mean(subphenoResults$V1)," / ", mean(subphenoResults$V2) )) #0.101852799"
print(paste0("combined Rapido/LDpred2 : ", mean(combinedResults$V1)," / ", mean(combinedResults$V2) )) #  0.115281442"
print(paste0("shaPRS Rapido/LDpred2 : ", mean(shaprsResults$V1)," / ", mean(shaprsResults$V2) )) # 0.13083603"
print(paste0("SMTPred Rapido/LDpred2 : ", mean(smtpredResults$V1)," / ", mean(smtpredResults$V2) )) # 0.13083603"
print(paste0("MTAG Rapido/LDpred2 : ", mean(mtagResults$V1)," / ", mean(mtagResults$V2) )) # 0.13083603"


print(paste0("overall difference favouring Rapido: ",mean(all_diff),"%")) # 6

