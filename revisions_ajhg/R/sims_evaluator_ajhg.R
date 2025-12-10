# create heatmap for simulation results
library(pheatmap)
library(colorspace)
library(grid)
# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 3) {stop("not enough Arguments received")} 
 

outputLoc= args[1]
allResultsLoc= args[2]
allPredictorsLoc= args[3]
isDiff=F
diff_x = 0.72
if( args[4] == "1") {
  isDiff=  T
  diff_x = 0.78
  }


# # Debug vars Orig
# outputLoc =   "C:/0Datasets/shaPRS/0revisions/results/5000_h2_0.75sims_eval_ukbb_revision"
# allResultsLoc = "C:/0Datasets/shaPRS/0revisions/results/causals5000_h2_0.75_all_results"
# allPredictorsLoc= "C:/0Datasets/shaPRS/0revisions/results/causals5000_h2_0.75_all_predictors"
# timeSeriesLoc ="C:/0Datasets/shaPRS/0revisions/results/5000_h2_0.75timeSeriesLoc"

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
dim(allResults)
allPredictors = read.csv(allPredictorsLoc, header = T)

# original order with shaPRS: "combined,subpheno,shaPRS,SMTPred,MTAG"
# order without shaPRS, as now all values are RELATIVE to shaPRS: "combined,subpheno,SMTPred,MTAG"

allResults <- allResults[, c(2, 1, 3, 4)] # reorder it, make single to be first, meta next, then smtpred, then mtag

modelNames = names(allResults)
modelNames[[1]] = "proximal only"
modelNames[[2]] = "meta-analysis"

names(allResults) = names(allResults_orig) = modelNames

i=4


sampleSizeAsN = as.character(allPredictors$sample_size)
sampleSizeAsN[which(sampleSizeAsN == "_half")] = 8000
sampleSizeAsN[which(sampleSizeAsN == "_full")] = 16000
sampleSizeAsN[which(sampleSizeAsN == "_double")] = 32000
sampleSizeAsN = as.numeric(sampleSizeAsN)

######################################################
# I) Plot all of the data




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

  sample_size = gsub("_", "", allPredictors[i,7]) #remove underscores
  rG= padTo_dec( allPredictors[i,2], 4)
  p= padTo_dec( round(allPredictors[i,3], 2) ,4)
  corre = padTo_dec( round(allPredictors[i,4], 2) ,4)
  
  annotations = paste0("rG:",rG, " p:",p, " cor:",corre, " ",allPredictors[i,5], " split:",allPredictors[i,6], " size:",sample_size)
  
  coln = c(coln,annotations)
  
  splits = c(splits,as.character(allPredictors[i,6] ) )
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
# pheatmap(allResults, main = "shaPRS - all simulations", filename=paste(outputLoc,".png", sep="" ),annotation_row=rGs_DF, height=25, width=10 , cex=1 ,cluster_rows=F, cluster_cols=F, angle_col=45)


filename=paste(outputLoc,".png", sep="" )
png(filename, width=740 , height=640 * 8/5, res =128);
labs.row = rep ("", nrow(allResults))
pheatmap(allResults, labels_row = labs.row, main = "",annotation_row=rGs_DF, height=8, width=5 , cex=1, fontsize_col = 14, fontsize = 12, cellwidth=40 ,cluster_rows=F, cluster_cols=F, angle_col=45)

if(isDiff) { grid.text( expression(bold('-log'[10]*'p diff') )  , x= diff_x, y = 0.98, gp = gpar(fontsize = 16, fontfamily = "sans"))
} else {grid.text( expression(bold('-log'[10]*'p') )  , x= diff_x, y = 0.98, gp = gpar(fontsize = 16, fontfamily = "sans"))}

dev.off()



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
plotName =""



# Filter to keep the main interesting scenarios, rG 0.5, regular, full
results_RG05 = filterOutByTerm_all(allResults,splits,ps,cors,filterTerm = "rG:0.50")
results_regular = filterOutByTerm_all(results_RG05$subsetResults,results_RG05$splits_subset,results_RG05$ps_subset,results_RG05$cors_subset,filterTerm = "regular")
results_full = filterOutByTerm_all(results_regular$subsetResults,results_regular$splits_subset,results_regular$ps_subset,results_regular$cors_subset,filterTerm = "size:full")


subsetResults = results_full$subsetResults

subset_DF <- data.frame( results_full$ps_subset, results_full$cors_subset,results_full$splits_subset,row.names=rownames(subsetResults)) # ,row.names=rownames(subsetResults) # need to match the rownames to the data for the annotation to work
colnames(subset_DF) <- c("p","cor","split") # this is the header for the annotation





###

# Filter to keep the main interesting scenarios, rG 0.5, extra, full
results_RG05 = filterOutByTerm_all(allResults,splits,ps,cors,filterTerm = "rG:0.50")
results_regular = filterOutByTerm_all(results_RG05$subsetResults,results_RG05$splits_subset,results_RG05$ps_subset,results_RG05$cors_subset,filterTerm = "extra")
results_full = filterOutByTerm_all(results_regular$subsetResults,results_regular$splits_subset,results_regular$ps_subset,results_regular$cors_subset,filterTerm = "size:full")


subsetResults_extra = results_full$subsetResults

subset_DF_extra <- data.frame( results_full$ps_subset, results_full$cors_subset,results_full$splits_subset,row.names=rownames(subsetResults_extra)) # ,row.names=rownames(subsetResults_extra) # need to match the rownames to the data for the annotation to work
colnames(subset_DF_extra) <- c("p","cor","split") # this is the header for the annotation

#plotName = "shaPRS - rG:0.5, n:full, no extra" # no plotname for final publication

# make pheatmap on the same colour scale:
#https://stackoverflow.com/questions/40226268/r-pheathmap-how-to-set-the-same-scale-in-two-heatmap-images
Breaks <- seq(min(subsetResults, subsetResults_extra), max(subsetResults, subsetResults_extra), length = 100)

## draw   
#pheatmap(test1, breaks = Breaks, cluster_rows=FALSE, cluster_cols=FALSE)
#pheatmap(test2, breaks = Breaks, cluster_rows=FALSE, cluster_cols=FALSE)


# we actually do want the categories to be shown as 'categorical' (IE explicitly showing that 'p=0.50' etc, but want the colours on a gradient)
# we we keep the dataframe's columns as 'as.character', but generate colors
cl <- c("#CE6BAF", "#797CBA")
HCL = rbind("0%" = cl, "40%" = lighten(cl, 0.4), "80%" = lighten(cl, 0.8))
#swatchplot(HCL)
# overwrite annotation colors: https://stackoverflow.com/questions/48678910/change-colors-in-multiple-annotations
annotation_colors = list(
  cor = c("0.50"=as.character(HCL[3,1]), "0.67"=as.character(HCL[2,1]), "1.00"=as.character(HCL[1,1]) ),
  p   = c("0.50"=as.character(HCL[3,2]), "0.75"=as.character(HCL[2,2]), "1.00"=as.character(HCL[1,2]) ))

# to manipulate fontsize and margins, have to use the cellwidth args: https://stackoverflow.com/questions/44701207/column-labels-cropped-when-using-pheatmap
filename=paste(outputLoc,"_subset.png", sep="" )  
png(filename, width=740 , height=640, res =128);
pheatmap(subsetResults, breaks = Breaks, main = plotName         ,annotation_row=subset_DF      , show_rownames = F, legend= T, angle_col = 45, height=5, width=5 , cex=1, fontsize_col = 14, fontsize = 12 ,cluster_rows=F, cluster_cols=F, cellwidth=45, annotation_colors = annotation_colors, display_numbers = T)
if(isDiff) { grid.text( expression(bold('-log'[10]*'p diff') )  , x= diff_x, y = 0.98, gp = gpar(fontsize = 16, fontfamily = "sans"))
} else {grid.text( expression(bold('-log'[10]*'p') )  , x= diff_x, y = 0.98, gp = gpar(fontsize = 16, fontfamily = "sans"))}
dev.off()


#font family not found in Windows font database
filename=paste(outputLoc,"_subset_extra.png", sep="" )
png(filename, width=740 , height=640, res =128);
pheatmap(subsetResults_extra, breaks = Breaks, main = plotName ,annotation_row=subset_DF_extra, show_rownames = F, angle_col = 45, height=5, width=5 , cex=1, fontsize_col = 14, fontsize = 12 ,cluster_rows=F, cluster_cols=F, cellwidth=45, annotation_colors = annotation_colors, display_numbers = T)
if(isDiff) { grid.text( expression(bold('-log'[10]*'p diff') )  , x= diff_x, y = 0.98, gp = gpar(fontsize = 16, fontfamily = "sans"))
} else {grid.text( expression(bold('-log'[10]*'p') )  , x= diff_x, y = 0.98, gp = gpar(fontsize = 16, fontfamily = "sans")) }
dev.off()

# windowsFonts() # the only way to specify font family, is to look here
# set expression in bold: https://stackoverflow.com/questions/31702480/grid-text-fails-to-display-bold-expression
# add text to arbitrary location via grid: https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grid.text.html
#######################################################################################




######################################################


