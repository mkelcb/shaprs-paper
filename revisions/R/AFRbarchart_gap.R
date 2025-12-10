
# 6th place is white so that it won't show!
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#332288", "#AA4499", "#117733", "#FFFFFF", "#DDCC77")

#myDf = prsResults_r2_sd
# adds a blank dummy result, also moves PRS-CSx to be last
addGap = function(myDf) {
  myDf_new = rbind(myDf, rep(0, ncol(myDf))) ## add dummy 
  myDf_new = myDf_new[-5,]
  myDf_new =rbind(myDf_new, myDf[5,])
  origRownames = rownames(myDf_new)
  origRownames[6] = "" # make 
  origRownames[7] = "PRS-CSx" # make 
  rownames(myDf_new) = origRownames
  return(myDf_new)
}

createPlot = function(prsResults, prsResults_sd_low, prsResults_sd_high, ylabel, filn="r2", showHeight =T, plotWidth=1280, xlimit = c(0,45), minY = 0, cexV = 1.5, showValues =F) {
  
  col <- safe_colorblind_palette[1: nrow(prsResults)]
  
  if(showHeight == F) { # remove height from AUC results, as that is a continuous trait
    index = which(colnames(prsResults) == "height") # we need to check against no match..
    if(length(index) > 0)
    {
      prsResults = subset(prsResults, select = -which(colnames(prsResults) == "height"))
      prsResults_sd_low = subset(prsResults_sd_low, select = -which(colnames(prsResults_sd_low) == "height"))
      prsResults_sd_high = subset(prsResults_sd_high, select = -which(colnames(prsResults_sd_high) == "height"))
    }
  }
  
  
  prsResults.matrix <- as.matrix(prsResults)
  # create matrix of errors
  prsResults.error_low = as.matrix(prsResults_sd_low)
  prsResults.error_high = as.matrix(prsResults_sd_high)
  if(minY == 0 & showHeight == F) {
    minY=0.5 # for AUCs we only care  whats above 0.5
  }
  
  # for the legend we dont want to show the dummy 
  rowWithNoName = which(rownames(prsResults.matrix) == "")
  legendMatrix = prsResults.matrix[-rowWithNoName,]
  legendcols = col[-rowWithNoName]
  borders = rep("black", nrow(prsResults.matrix)) # c("red", "green", "yellow", "blue")
  borders[rowWithNoName] = "#FFFFFF"
  
  filename=paste0(outputLoc,filn, ".png")
  
  png(filename, width=plotWidth , height=480, res =128)
  par(mar=c(2.0,6,0.5,2)+0.1,mgp=c(5,1,0)) # par(mar = c(bottom, left, top, right))
  # create bar plot
  #par(mar=c(0.0,0,0.0,0),mgp=c(0,0,0)) # par(mar = c(bottom, left, top, right))
  
  bp <- barplot(prsResults.matrix, xpd = F,
                beside=T,
                col=col,
                ylim=c(minY,max(prsResults) + max(prsResults) * 0.1),
                xlab="",
                border = borders,
                space = c(0, 2.5), # adds more space between groups 
                # ylab=parse(text=paste(ylabel)),
                #legend.text=row.names(prsResults.matrix),
                #args.legend=list("topright",bty="n", ncol = 1, cex =cexV, y.intersp=1.8),
                las=1,
                xlim = xlimit,
                cex=1.2,
                cex.lab = 1.2,
                cex.axis=1.2)
  title(ylab=parse(text=paste(ylabel)), line=3.7, cex.lab=1.5)
  if(showValues) {
    text(x=bp, y= prsResults.matrix+0.04, cex = 1.5, labels=as.character(prsResults.matrix))
  }
  
  
  
  legend("topright", bty="n", ncol = 1, cex =cexV, y.intersp=1.8,
         legend=row.names(legendMatrix), fill=legendcols)
  
  
  if(showHeight == F) { # dont plot error bars for height runs, IE when it is r^2, as those are too small
    # add error bars
    arrows(bp, prsResults.error_low, bp, prsResults.error_high,
           code = 3, angle = 90, length = 0.15)
  }
  
  # Restore default clipping rect
  #par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  
}

###############################################
# AFR results
# loads and creates barchart for the AFR, and then it merges them

prsResults_r2= read.csv("C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/AFR_res.csv", header = T,row.names = 1)
prsResults_r2_sd= read.csv("C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/AFR_res_sd.csv", header = T,row.names = 1)
outputLoc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/AFR_res_gap"

# add dummy gaps to separate out PRS-CSx
prsResults_r2 = addGap(prsResults_r2)
prsResults_r2_sd = addGap(prsResults_r2_sd)

createPlot(prsResults_r2, prsResults_r2_sd, prsResults_r2_sd, "r^2", filn="r2", T, cexV = 1.2)

######################

# Original trans-ethnic results:
# parse results into 'copy pastable' from Elena's output
resultsLoc= read.table("C:/0Datasets/shaPRS/Elena/Elena_161121/auc_r2.txt", header = T)

#outputLoc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/final/gap"
outputLoc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/gap"

pheno="height"
resultsLoc_pheno = cbind.data.frame(paste0(pheno,"_EUR_EUR_LDpred2_CSxsubset_E"), 0,0,0,0.0976, 6.31e-06,pheno)
df = cbind.data.frame(paste0(pheno,"_EUR_EUR_PRSCS_E"), 0,0,0,0.116, 5.96e-06,pheno)
colnames(df) = colnames(resultsLoc_pheno)
resultsLoc_pheno = rbind(resultsLoc_pheno,df)


df = cbind.data.frame(paste0(pheno,".PRSCSx_E"), 0,0,0,0.123, 1.1e-05,pheno)
colnames(df) = colnames(resultsLoc_pheno)
resultsLoc_pheno = rbind(resultsLoc_pheno,df)

df = cbind.data.frame(paste0(pheno,"_EUR_PRSCSx_E"), 0,0,0,0.121, 5.48e-06,pheno)
colnames(df) = colnames(resultsLoc_pheno)
resultsLoc_pheno = rbind(resultsLoc_pheno,df)

df = cbind.data.frame(paste0(pheno,"_shaPRS_LDpred2_CSxsubset_E"), 0,0,0,0.122, 5.37e-06,pheno)
colnames(df) = colnames(resultsLoc_pheno)
resultsLoc_pheno = rbind(resultsLoc_pheno,df)

df = cbind.data.frame(paste0(pheno,"_shaPRS_EUR_PRSCS_E"), 0,0,0,0.122, 5.71e-06,pheno)
colnames(df) = colnames(resultsLoc_pheno)
resultsLoc_pheno = rbind(resultsLoc_pheno,df)

colnames(resultsLoc_pheno) = colnames(resultsLoc)
resultsLoc = rbind(resultsLoc,resultsLoc_pheno)


###################################


pheno="asthma"
resultsLoc_pheno = cbind.data.frame(paste0(pheno,"_EUR_EUR_LDpred2_CSxsubset_E"), 0.595,0.5919,0.5989 ,0.0115, 5.52e-06,pheno)


df = cbind.data.frame(paste0(pheno,"_EUR_EUR_PRSCS_E"), 0.588,0.5849,0.5918,0.00986, 5.65e-06,pheno)
colnames(df) = colnames(resultsLoc_pheno)
resultsLoc_pheno = rbind(resultsLoc_pheno,df)


df = cbind.data.frame(paste0(pheno,".PRSCSx_E"), 0.603,0.5981,0.6079,0.0134,1.11e-05,pheno)
colnames(df) = colnames(resultsLoc_pheno)
resultsLoc_pheno = rbind(resultsLoc_pheno,df)

df = cbind.data.frame(paste0(pheno,"_EUR_PRSCSx_E"), 0.595,0.5916,0.5985 ,0.0113, 6.35e-06,pheno)
colnames(df) = colnames(resultsLoc_pheno)
resultsLoc_pheno = rbind(resultsLoc_pheno,df)

df = cbind.data.frame(paste0(pheno,"_shaPRS_LDpred2_CSxsubset_E"), 0.604,0.6004,0.6073,0.0136, 5.53e-06,pheno)
colnames(df) = colnames(resultsLoc_pheno)
resultsLoc_pheno = rbind(resultsLoc_pheno,df)

df = cbind.data.frame(paste0(pheno,"_shaPRS_EUR_PRSCS_E"), 0.599,0.5951,0.602,0.0123, 5.55e-06,pheno)
colnames(df) = colnames(resultsLoc_pheno)
resultsLoc_pheno = rbind(resultsLoc_pheno,df)

colnames(resultsLoc_pheno) = colnames(resultsLoc)
resultsLoc = rbind(resultsLoc,resultsLoc_pheno)


# create data frame for relevant vars for r^2
pheno="asthma"


prsResults_r2 = NULL
prsResults_AUC = NULL
prsResults_r2_sd = NULL
prsResults_AUC_sd_high = NULL
prsResults_AUC_sd_low = NULL
traits = c("asthma", "T2D", "BRCA", "CAD", "height")
for (i in 1:length(traits)) {
  pheno = traits[i]
  
  resultsLoc_pheno = resultsLoc[grep(pheno, resultsLoc$model),]
  
  
  # primary Ldpred2
  primary_ldpred2_r2 = resultsLoc_pheno$r2[grep("_EUR_EUR_LDpred2_CSxsubset", resultsLoc_pheno$model)]
  primary_ldpred2_AUC = resultsLoc_pheno$AUC[grep("_EUR_EUR_LDpred2_CSxsubset", resultsLoc_pheno$model)]
  primary_ldpred2_r2_sd = resultsLoc_pheno$sd[grep("_EUR_EUR_LDpred2_CSxsubset", resultsLoc_pheno$model)]
  primary_ldpred2_AUC_sd_high = resultsLoc_pheno$ci_hight[grep("_EUR_EUR_LDpred2_CSxsubset", resultsLoc_pheno$model)]
  primary_ldpred2_AUC_sd_low = resultsLoc_pheno$ci_low[grep("_EUR_EUR_LDpred2_CSxsubset", resultsLoc_pheno$model)]
  
  
  # primary PRS-CS
  primary_prscs_r2 = resultsLoc_pheno$r2[grep("_EUR_EUR_PRSCS_E", resultsLoc_pheno$model)]
  primary_prscs_AUC =resultsLoc_pheno$AUC[grep("_EUR_EUR_PRSCS_E", resultsLoc_pheno$model)]
  primary_prscs_r2_sd = resultsLoc_pheno$sd[grep("_EUR_EUR_PRSCS_E", resultsLoc_pheno$model)]
  primary_prscs_AUC_sd_high = resultsLoc_pheno$ci_hight[grep("_EUR_EUR_PRSCS_E", resultsLoc_pheno$model)]
  primary_prscs_AUC_sd_low = resultsLoc_pheno$ci_low[grep("_EUR_EUR_PRSCS_E", resultsLoc_pheno$model)]
  
  
  # primary+adjunct PR-SCSx
  pooledData_prscsx_r2 = resultsLoc_pheno$r2[grep(".PRSCSx_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_prscsx_AUC = resultsLoc_pheno$AUC[grep(".PRSCSx_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_prscsx_r2_sd = resultsLoc_pheno$sd[grep(".PRSCSx_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_prscsx_AUC_sd_high = resultsLoc_pheno$ci_hight[grep(".PRSCSx_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_prscsx_AUC_sd_low = resultsLoc_pheno$ci_low[grep(".PRSCSx_E", resultsLoc_pheno$model, fixed=T)]
  
  
  
  # primary+adjunct PRS-CSx-unweighted
  pooledData_prscsx_unweighted_r2 = resultsLoc_pheno$r2[grep("_EUR_PRSCSx_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_prscsx_unweighted_AUC =resultsLoc_pheno$AUC[grep("_EUR_PRSCSx_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_prscsx_unweighted_r2_sd = resultsLoc_pheno$sd[grep("_EUR_PRSCSx_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_prscsx_unweighted_AUC_sd_high = resultsLoc_pheno$ci_hight[grep("_EUR_PRSCSx_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_prscsx_unweighted_AUC_sd_low = resultsLoc_pheno$ci_low[grep("_EUR_PRSCSx_E", resultsLoc_pheno$model, fixed=T)]
  
  
  # primary+adjunct Ldpred2 + shaPRS
  pooledData_shaprs_ldpred2_unweighted_r2 = resultsLoc_pheno$r2[grep("_shaPRS_LDpred2_CSxsubset_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_shaprs_ldpred2_unweighted_AUC = resultsLoc_pheno$AUC[grep("_shaPRS_LDpred2_CSxsubset_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_shaprs_ldpred2_unweighted_r2_sd = resultsLoc_pheno$sd[grep("_shaPRS_LDpred2_CSxsubset_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_shaprs_ldpred2_unweighted_AUC_sd_high = resultsLoc_pheno$ci_hight[grep("_shaPRS_LDpred2_CSxsubset_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_shaprs_ldpred2_unweighted_AUC_sd_low = resultsLoc_pheno$ci_low[grep("_shaPRS_LDpred2_CSxsubset_E", resultsLoc_pheno$model, fixed=T)]
  
  
  # primary+adjunct PRS-CS +shaPRS-unweighted
  pooledData_shaprs_prcs_unweighted_r2 = resultsLoc_pheno$r2[grep("_shaPRS_EUR_PRSCS_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_shaprs_prcs_unweighted_AUC = resultsLoc_pheno$AUC[grep("_shaPRS_EUR_PRSCS_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_shaprs_prcs_unweighted_r2_sd = resultsLoc_pheno$sd[grep("_shaPRS_EUR_PRSCS_E", resultsLoc_pheno$model, fixed=T)]
  pooledData_shaprs_prcs_unweighted_AUC_sd_high = resultsLoc_pheno$ci_hight[grep("_shaPRS_EUR_PRSCS_E", resultsLoc_pheno$model, fixed=T)] 
  pooledData_shaprs_prcs_unweighted_AUC_sd_low = resultsLoc_pheno$ci_low[grep("_shaPRS_EUR_PRSCS_E", resultsLoc_pheno$model, fixed=T)] 
  
  
  # cols are traits, 
  # PRS types of the rows
  
  prsResults_r2 = cbind(prsResults_r2,   c(primary_ldpred2_r2, primary_prscs_r2, pooledData_shaprs_ldpred2_unweighted_r2, pooledData_shaprs_prcs_unweighted_r2, pooledData_prscsx_r2, pooledData_prscsx_unweighted_r2 ))
  prsResults_AUC = cbind(prsResults_AUC,   c(primary_ldpred2_AUC, primary_prscs_AUC, pooledData_shaprs_ldpred2_unweighted_AUC, pooledData_shaprs_prcs_unweighted_AUC, pooledData_prscsx_AUC, pooledData_prscsx_unweighted_AUC ))
  
  prsResults_r2_sd = cbind(prsResults_r2_sd,   c(primary_ldpred2_r2_sd, primary_prscs_r2_sd, pooledData_shaprs_ldpred2_unweighted_r2_sd, pooledData_shaprs_prcs_unweighted_r2_sd, pooledData_prscsx_r2_sd, pooledData_prscsx_unweighted_r2_sd ))
  prsResults_AUC_sd_high = cbind(prsResults_AUC_sd_high,   c(primary_ldpred2_AUC_sd_high, primary_prscs_AUC_sd_high, pooledData_shaprs_ldpred2_unweighted_AUC_sd_high, pooledData_shaprs_prcs_unweighted_AUC_sd_high, pooledData_prscsx_AUC_sd_high, pooledData_prscsx_unweighted_AUC_sd_high ))
  prsResults_AUC_sd_low = cbind(prsResults_AUC_sd_low,   c(primary_ldpred2_AUC_sd_low, primary_prscs_AUC_sd_low, pooledData_shaprs_ldpred2_unweighted_AUC_sd_low, pooledData_shaprs_prcs_unweighted_AUC_sd_low, pooledData_prscsx_AUC_sd_low, pooledData_prscsx_unweighted_AUC_sd_low ))
  
}

rownames_trait = c("LDpred2", "PRS-CS", "shaPRS+LDpred2", "shaPRS+PRS-CS", "PRS-CSx", "PRS-CSx-stage1")


rownames(prsResults_r2) = rownames_trait
rownames(prsResults_AUC) = rownames_trait
colnames(prsResults_r2) = traits
colnames(prsResults_AUC) = traits

rownames(prsResults_r2_sd) = rownames_trait
rownames(prsResults_AUC_sd_high) = rownames_trait
rownames(prsResults_AUC_sd_low) = rownames_trait
colnames(prsResults_r2_sd) = traits
colnames(prsResults_AUC_sd_high) = traits
colnames(prsResults_AUC_sd_low) = traits


prsResults_r2_v2 = addGap(prsResults_r2)
# add dummy gaps to separate out PRS-CSx
prsResults_r2 = addGap(prsResults_r2)
prsResults_r2_sd = addGap(prsResults_r2_sd)
prsResults_AUC = addGap(prsResults_AUC)
prsResults_AUC_sd_low = addGap(prsResults_AUC_sd_low)
prsResults_AUC_sd_high = addGap(prsResults_AUC_sd_high)

createPlot(prsResults_r2, prsResults_r2_sd, prsResults_r2_sd, "r^2", filn="r2", T, xlimit = c(0,65), cexV = 1.2)
createPlot(prsResults_AUC,prsResults_AUC_sd_low, prsResults_AUC_sd_high, "AUC", filn="auc", F, xlimit = c(0,52.5), cexV = 1.2)



##########################
# merge them:
# GenerateImage() relies on ImageCombiner.R
library(magick)

outputLoc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/AFR_res_gap"
baseloc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/gap"


# CROSS ANCESTRY BARPLOTS
imagesPerRow= 1 # how many images per row
outName="crossAncestryBarplot_publish_rev_gap.png"


image2 <- image_read(paste0(baseloc,"r2.png"))
image3 <- image_read(paste0(outputLoc,"r2.png"))

images = list(image2,image3) # image1, 
labels = list("a", "b") # , "c"
img = GenerateImage(images, labels, imagesPerRow)
image_write(img, paste0("C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/",outName))



