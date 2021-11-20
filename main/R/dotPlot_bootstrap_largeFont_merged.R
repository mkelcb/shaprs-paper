# produces dot plots, from a list of files (that have a col of values) (these are the bootstrap samples), adds horizontal line for the main prediction

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 6) {stop("not enough Arguments received")} 



outputLoc= args[1]
outputName = args[2]
plotName = args[3]
ceilingValue = as.double(args[4])

colour1 = args[5]
colour2 = args[6]

inputFiles = vector(length = (length(args) -6) )
counter = 1
for (i in 7:length(args)) { # loop through the rest of the arguments, where each is supposed to be an input file that stores the r^2 values for each fold
  inputFiles[counter] = args[i]
  counter = counter +1
}

#fullPath="C:/softwares/Cluster/GIANT/miniPRS/CD_from_CD"
# returns a box name, by cutting out file path and file extension
getBoxName  = function(fullPath) {
  name = unlist( strsplit(fullPath, "/", fixed = TRUE) )
  name = name[length(name)]  # get rid of the path
  #name = unlist( strsplit(name, ".", fixed = TRUE) )[1] # get rid of the file extension (if any)

  print(name)
  name = paste("predicted: ", gsub("_from_", "\ntrained: ", name), sep="" )
  if ( grepl("Comp_", name, fixed = TRUE) ) { # want to replace the cryptic 'Comp_.099' with 'Composite'
    withoutComp = unlist(strsplit(name, "Comp_"))[1]
    name = paste(withoutComp, "Composite", sep="")
  }
  return(name)
}
# uga = "\ntrained: Comp_0.99"
# uga2 = "trained:Blend"
# gsub("_", " ", plotName)
# 
# grepl("Comp_", uga, fixed = TRUE)



# # Debug vars
# inputFiles = c("C:/softwares/Cluster/GIANT/miniPRS/CD_CD_PRS", "C:/softwares/Cluster/GIANT/miniPRS/CD_IBD_PRS", "C:/softwares/Cluster/GIANT/miniPRS/UC_UC_PRS" ,"C:/softwares/Cluster/GIANT/miniPRS/UC_IBD_PRS")
# outputName = "IBD_subpheno_PRS"
# plotName = "IBD_PRS_performance_against_subphenotypes"
# ceilingValue = -1
# outputLoc = "C:/softwares/Cluster/GIANT/miniPRS/IBD_Subphenos"
# colour="#166938"
  # ASHG green: 166938 # UC
  # ASHG yellow: c5b000  # CD

#
# 
plotName = gsub("_", " ", plotName) #remove underscores

########################################################

## ____________________________________________________________
i=1
# attempt to load in as many input files as they exist
inputDatas = list()

for(i in 1:length(inputFiles)) {
  inputDatas[[i]] =read.table( inputFiles[i]  ,header=F)

}

#allColours = rainbow(length(inputFiles)) # get some nice distinct colours based on the number of elements we want to plot
#allColours[allColours== "#FFFF00FF"] <- "black" # yellow is poor for reading with letters, replace it
allColours1 = rep(colour1, length(inputFiles)/2)
allColours2 = rep(colour2, length(inputFiles)/2)
allColours = c(allColours1,allColours2)

# load ^ process data
allData = NULL
allMeans= list() # = vector(length = length(inputFiles) )
boxNames = vector(length = length(inputFiles) )
for (i in 1:length(inputFiles)) {
  data = as.matrix(inputDatas[[i]])
  numNAs = length( which(is.na(data)) )
  if(numNAs > 0 ) { print("Input has NAs which were replaced by column mean")}
  data[is.na(data)] <- mean(data, na.rm = TRUE) # replace NA's by data mean
  
  allData = cbind(allData, data )
  allMeans[[i]] = mean( data )
  boxNames[i] = getBoxName(inputFiles[i])
  #boxNames[i] = paste(boxNames[i],"\n(u=", round(allMeans[[i]], 3),")", sep="")
}

allData[is.na(allData)] <- 0
minVal = min(allData) * 0.9
maxVal = max(allData) * 1.1
overallMean = mean( unlist(allMeans) )
if (ceilingValue != -1 ) { maxVal = ceilingValue * 1.1 }


# perform paired 2 sample t-tests
allCombs= combn(1:ncol(allData), 2)
for ( i in 1:ncol(allCombs)) {
  index_1 = allCombs[1,i]
  index_2 = allCombs[2,i]
  mean1 = mean(allData[,index_1])
  mean2 = mean(allData[,index_2])
  percDiff = round( (mean1 - mean2) / ( (mean1 + mean2)/2 ) * 100)
  
  print( paste(boxNames[index_1],"v",boxNames[index_2],":", t.test(allData[,index_1],allData[,index_2], paired = T)$p.value, "/ diff:", percDiff, "%") )

  print(paste0(boxNames[index_1],": ",mean1 ))
}


print("JUST THE MEANS")
for ( i in 1:ncol(allData)) {
  mean1 = mean(allData[,i])
  print(paste0(boxNames[i],": ",mean1 ))
}







#dev.off()

#Offset=c(-1,+1,-1,+1)  # for loo cohorts:  c(-1,+1,-1,+1,-1)
#Offset = Offset * 0.25




inputDatas_mat = list()
for (i in 1:length(inputFiles)) {
    inputDatas_mat[[i]] = unlist(inputDatas[[i]]) 
}

plotName="" # plotname disabled for publication
filen =paste(outputLoc,outputName,".pdf", sep="" )
pdf(filen, width=3.5* length(boxNames) , height=6.4);

#par( mfrow = c( 1, 1 ) ,mgp = c(1, 2.5, 0))
borderEnabled=allColours
boxplot(allMeans , outline = FALSE, boxlty = 0, whisklty = 0, staplelty = 0,xaxt="n", border =borderEnabled, names=boxNames, ylab=parse(text="r^2"), main=plotName, ylim=c(minVal, maxVal), cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.2, cex.main = 1.8)
stripchart( inputDatas_mat, vertical = TRUE, pch = 16, col = allColours, add = TRUE, method = "jitter", cex = 2) #
axis(1,at=c(1:length(boxNames)),labels=boxNames, cex.axis = 2,mgp=c(3,3.5,0))
if (ceilingValue != -1 ) {
abline(h = ceilingValue, col="green" , lwd=3)
legend("topleft", legend=c("heritability"),col=c("green"), text.col= c("green"), lty=1, lwd=3)
}
#abline(h = overallMean, col="grey" , lwd=2,lty=3) 

dev.off()




filen =paste(outputLoc,outputName,".png", sep="" )
png(filen, width=258.75 * length(boxNames) , height=200*3);

#par( mfrow = c( 1, 1 ) ,mgp = c(1, 2.5, 0))
borderEnabled=allColours
boxplot(allMeans , outline = FALSE, boxlty = 0, whisklty = 0, staplelty = 0,xaxt="n", border =borderEnabled, names=boxNames, ylab=parse(text="r^2"), main=plotName, ylim=c(minVal, maxVal), 
        cex.lab = 1.8, cex.axis = 2, cex.names = 1.2, cex.main = 2, par(mar = c(4.5, 5, 2, 2)))
stripchart( inputDatas_mat, vertical = TRUE, pch = 16, col = allColours, add = TRUE, method = "jitter", cex = 2) #
axis(1,at=c(1:length(boxNames)),labels=boxNames, cex.axis = 2.15,mgp=c(3,3.5,0) )
if (ceilingValue != -1 ) {
  abline(h = ceilingValue, col="green" , lwd=3)
  legend("topleft", legend=c("heritability"),col=c("green"), text.col= c("green"), lty=1, lwd=3)
}
#abline(h = overallMean, col="grey" , lwd=2,lty=3)

dev.off()

print(paste("written plot to",filen))

        
#boxplot(count ~ spray, data = InsectSprays, col = "lightgray", ylab=paste0("predictive power (",parse(text="r^2"),")")   )

#boxplot(count ~ spray, data = InsectSprays, col = "lightgray", ylab=parse(text="r^2")   )
