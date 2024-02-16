# gets  the median of a list of supplied files, assumes that it will have 3 columns

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
lrt_diffs=""
for(i in seq(from=2, to=length(args))){
  res= read.table(args[i], header = F) 
  res[is.na(res[,1]), 1] <- mean(res[,1], na.rm = TRUE) # mean impute missing values
  
  res[is.na(res[,2]), 2] <- mean(res[,2], na.rm = TRUE) # mean impute missing values
  res[is.na(res[,3]), 3] <- mean(res[,3], na.rm = TRUE) # mean impute missing values
  
  # get the median r2diff p value
  r2diff_p = median(res[,1])
  r2diff_p = -log10(r2diff_p)
  
  # get the -log10p lrt difference
  prs1_p = median(res[,2])
  prs2_p = median(res[,3])
  
  prs1_p = -log10(prs1_p)
  prs2_p = -log10(prs2_p)
 prs2_prs1_diff = prs2_p - prs1_p
   
 # res = unlist(res)
  if (i ==2) {
    medians = r2diff_p 
    lrt_diffs = prs2_prs1_diff
  } else {
    medians = paste0(medians,",",r2diff_p) 
    lrt_diffs = paste0(lrt_diffs,",",prs2_prs1_diff) 
    }
}  


write(medians,file=paste0(outputLoc,"/r2diffs") ) # ,append=TRUE
write(lrt_diffs,file=paste0(outputLoc,"/lrts") ) # ,append=TRUE

print(paste("written: median r2diff and lrt_diffs p to", outputLoc ))

