# Performs Step II of PRS-CSx, by finding the best linear combination of 2 supplied population's PRS', then evaluates them on the witheld test set
#options(error=traceback)
# grab command line arguments

library(r2redux)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 6) {stop("not enough Arguments received")} 


PRS1ProfileLoc= args[1] # individual PRS profiles for EUR indis
PRS2ProfileLoc = args[2] # individual PRS profiles for EAS indis
famFileLoc = args[3] # the true phenotypes file (in PLINK 1 format)
testSetListLoc= args[4] # list of indis to be used in the test set
isBin= args[5] == "1" # flag to decide if binary (1) or quantitative pheno (0)
outputLoc = args[6] # where to write results


# Load data:
PRS1Profile = read.table(PRS1ProfileLoc ,header=T)
PRS2Profile = read.table(PRS2ProfileLoc ,header=T)
famFile = read.table(famFileLoc ,header=F)
testSetList= read.table(testSetListLoc ,header=F)

# if the famfile is actually a .phe file, we will have different columns, in the .phe file, V3 has the phenotype, not V6, but the rest of the code expects it to be V6
pheFile=grepl(".phe", famFileLoc, fixed = TRUE)
if(pheFile){
  names(famFile)[names(famFile) == 'V3'] <- 'V6'
}

# standardize to PRS to zscores
PRS1Profile$SCORESUM = scale(PRS1Profile$SCORESUM)
PRS2Profile$SCORESUM = scale(PRS2Profile$SCORESUM)


# merge into common DF
TargetDF = merge(PRS1Profile, PRS2Profile, by="IID")
TargetDF = merge(TargetDF, famFile, by.x="IID", by.y="V2")



# exlclude missing pheno people (-9 in UKBB)
numbefore=nrow(TargetDF)
TargetDF = TargetDF[which(TargetDF$V6 != -9),]
TargetDF$origPheno = TargetDF$V6 # save away the original, unscaled pheno for the AUC calculation to avoid the warning msg
TargetDF$V6 = scale(TargetDF$V6)
print(paste0("exluded number of missing pheno indis :" , numbefore -nrow(TargetDF) ))


# find the indices of the test indis
TESTIndices = match(testSetList$V1, TargetDF$IID) # find the indices in the targetDF that refer to the test set

# drop NAs (indis that could not be found)
TESTIndices= TESTIndices[is.na(TESTIndices) == F]

# subset full UKBB into Validation and Testing sets
ValidDF = TargetDF[-TESTIndices,]
TestdDF = TargetDF[TESTIndices,]

# Fit model on Validation sets
#if(isBin) {
#model = glm(V6 ~ SCORESUM.x + SCORESUM.y, data = ValidDF, family = "binomial")
#} else {
#  model = lm(V6 ~ SCORESUM.x + SCORESUM.y, data = ValidDF)
#}
# model without intercept
model = lm(V6 ~ SCORESUM.x + SCORESUM.y -1, data = ValidDF)

print(summary(model))

# write to disk model params 
sink(paste0(outputLoc,"_model.txt"))
print(summary(model))
sink()  # returns output to the console

# Predict on test set
PRS_final = predict(model, newdata=TestdDF)


####################

# get confidence intervals via r2redux
dat1 =  cbind(PRS_final,TestdDF$V6)
nv= nrow(dat1)
output1=r2_var(dat1,1,nv)


correlation = cor(PRS_final, TestdDF$V6)
correlation_sq = correlation^2
test = cor.test(PRS_final, TestdDF$V6)


testResults = cbind.data.frame(correlation,correlation_sq,test$p.value, output1$r2_based_p, output1$lower_r2, output1$upper_r2	)
colnames(testResults) = c("correlation","correlation^2","corr_p","r2redux_p", "CI_LOW", "CI_HIGH")



write.table(testResults,file=outputLoc, quote = F, row.names = F, sep = "\t",append=F)

print(paste("written: correlation_sq",signif(correlation_sq,3), "(CI_LOW:",signif(output1$lower_r,3), "CI_HIGH:",signif(output1$upper_r2,3),")", "to", outputLoc ))

# also output final PRS in the PLINK .sscore format
#head(TestdDF)
PLINKprofile =cbind.data.frame(TestdDF[,1], TestdDF[,2], TestdDF[,3], TestdDF[,4], TestdDF[,5], PRS_final)
colnames(PLINKprofile) = c("FID","IID","PHENO","CNT","CNT2","SCORESUM")
write.table(PLINKprofile,file=paste0(outputLoc,".sscore"), col.names = T , quote = F, row.names = F, sep = "\t",append=F)




# evaluate model accuracy via AUC (if binary)
if(isBin) {
  library(pROC)

  category= TestdDF$origPheno
  prediction= PRS_final
  #category <- c(1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0)
  #prediction <- rev(seq_along(category))
  #prediction[9:10] <- mean(prediction[9:10])
  roc_obj <- roc(category, prediction)
  myAUC= auc(roc_obj)[1]
  ci_res <- ci(roc_obj)
  filen= paste0(outputLoc, "_AUC")

  AUCTable=cbind(myAUC[1],ci_res[1] ,ci_res[3])
  colnames(AUCTable) = c("AUC", "CI_LOW", "CI_HIGH")
  write(AUCTable,file=filen)
  

  print(paste("written: AUC",signif(myAUC,3), "CI low:",signif(ci_res[1],4), " / high: ",signif(ci_res[3],4) , "to", filen ))
  
}
