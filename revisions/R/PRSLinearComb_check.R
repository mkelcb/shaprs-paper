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


# Predict on test set
PRS_final = predict(model, newdata=TestdDF)




correlation_AFR = cor(PRS_final, TestdDF$SCORESUM.x)
correlation_EUR = cor(PRS_final, TestdDF$SCORESUM.y)
print(paste("correlation_AFR",signif(correlation_AFR,3), " / correlation_EUR",signif(correlation_EUR,3) ))

