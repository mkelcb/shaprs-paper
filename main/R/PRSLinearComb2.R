# Performs Step II of PRS-CSx, by finding the best linear combination of 2 supplied population's PRS', then evaluates them on the witheld test set
#options(error=traceback)
# grab command line arguments
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


# standardize to PRS to zscores
PRS1Profile$SCORESUM = scale(PRS1Profile$SCORESUM)
PRS2Profile$SCORESUM = scale(PRS2Profile$SCORESUM)


# merge into common DF
TargetDF = merge(PRS1Profile, PRS2Profile, by="IID")
TargetDF = merge(TargetDF, famFile, by.x="IID", by.y="V2")


# exlclude missing pheno people (-9 in UKBB)
numbefore=nrow(TargetDF)
TargetDF = TargetDF[which(TargetDF$V6 != -9),]
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

# Elena's r^2 SD instructions:
# I have a data table with columns trait and score. From there I sample 1000 
# times with replacement the score column and look at the correlation squared 
# with the trait. I then calculate the sd and that is the SD estimate I gave you
# in the file above.
#sd <- sd(unlist(lapply(1:N, function(i) cor(dt[sample(1:nrow(dt), replace=T), score], dt[,trait])^2)))
library('data.table')
dt =  cbind.data.frame ( PRS_final,TestdDF$V6)
colnames(dt) = c("trait" ,"score")
setDT(dt)
N=1000
sd <- sd(unlist(lapply(1:N, function(i) cor(dt[sample(1:nrow(dt), replace=T), score], dt[,trait])^2)))


# Evaluate model accuracy as r^2
correlation = cor(PRS_final, TestdDF$V6)
correlation_sq = correlation^2
test = cor.test(PRS_final, TestdDF$V6)

testResults = cbind.data.frame(correlation,correlation_sq,test$p.value, sd)
colnames(testResults) = c("correlation","correlation^2","corr_p","corr^2_SD")
write.table(testResults,file=outputLoc, quote = F, row.names = F, sep = "\t",append=F)

print(paste("written: correlation_sq",signif(correlation_sq,3), "(sd:",signif(sd,3),")", "to", outputLoc ))


# evaluate model accuracy via AUC (if binary)
if(isBin) {
  library(pROC)

  category= TestdDF$V6
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
