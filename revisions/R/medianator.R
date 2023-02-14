# calculates the median of a list of supplied files

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 2) {stop("not enough Arguments received")} 
 

outputLoc= args[1]

# Debug vars
# gwasLoc = "C:/softwares/Cluster/0shaPRS/debug/1000/10/1.0_0.1/A50_B50/size_half/combined"  #"C:/softwares/Cluster/GIANT/miniPRS/PRS1"
# outputLoc =   "C:/softwares/Cluster/GIANT/miniPRS/shaPRS/1_pheA" #"C:/softwares/Cluster/GIANT/miniPRS/PRS_ALL"


# load input files for each method 
medians=""
for(i in seq(from=2, to=length(args))){
  res= read.csv(args[i], header = F) 
  res[is.na(res[,1]), 1] <- mean(res[,1], na.rm = TRUE)
  res = unlist(res)
  if (i ==2) {medians = median(res) }
  else {medians = paste0(medians,",",median(res)) }
}  


write(medians,file=outputLoc) # ,append=TRUE

print(paste("written: medianss to", outputLoc ))

