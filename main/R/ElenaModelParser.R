# parse results into 'copy pastable' from Elena's output

resultsLoc= read.table("C:/0Datasets/shaPRS/Elena/Elena_161121/auc_r2.txt", header = T)

outputLoc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/final/"

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








#######

pheno="height"
pheno="asthma"
pheno="BRCA"
pheno="T2D"
pheno="CAD"

# break out pheno:
resultsLoc_pheno = resultsLoc[grep(pheno, resultsLoc$model),]

if(pheno != "height" & pheno != "asthma")resultsLoc_pheno = resultsLoc_pheno[-grep("EAS", resultsLoc_pheno$model),] # remove EAS

print(pheno)

# add my results in from my bash script for height / asthma


# primary Ldpred2
print(paste0("primary Ldpred2 r2: ", signif(resultsLoc_pheno$r2[grep("_EUR_EUR_LDpred2_CSxsubset", resultsLoc_pheno$model)],3), " / AUC: ", signif(resultsLoc_pheno$AUC[grep("_EUR_EUR_LDpred2_CSxsubset", resultsLoc_pheno$model)],3) ))
#resultsLoc_pheno = resultsLoc_pheno[-grep("_EUR_EUR_LDpred2_CSxsubset", resultsLoc_pheno$model),]

# (resultsLoc_pheno$AUC - resultsLoc_pheno$ci_low) - (resultsLoc_pheno$ci_hight - resultsLoc_pheno$AUC)

# primary PRS-CS
print(paste0("primary PRS-CS r2: ", signif(resultsLoc_pheno$r2[grep("_EUR_EUR_PRSCS_E", resultsLoc_pheno$model)],3), " / AUC: ", signif(resultsLoc_pheno$AUC[grep("_EUR_EUR_PRSCS_E", resultsLoc_pheno$model)],3) ))
#resultsLoc_pheno = resultsLoc_pheno[-grep("_EUR_EUR_PRSCS_E", resultsLoc_pheno$model),] 


# primary+adjunct PR-SCSx
print(paste0("primary+adjunct PRS-CSx r2: ", signif(resultsLoc_pheno$r2[grep(".PRSCSx_E", resultsLoc_pheno$model, fixed=T)],3), " / AUC: ", signif(resultsLoc_pheno$AUC[grep(".PRSCSx_E", resultsLoc_pheno$model, fixed=T)],3) ))
#resultsLoc_pheno = resultsLoc_pheno[-grep(".PRSCSx_E", resultsLoc_pheno$model, fixed=T),] 


# primary+adjunct PRS-CSx-unweighted
print(paste0("primary+adjunct PRS-CS- unweighted r2: ", signif(resultsLoc_pheno$r2[grep("_EUR_PRSCSx_E", resultsLoc_pheno$model, fixed=T)],3), " / AUC: ", signif(resultsLoc_pheno$AUC[grep("_EUR_PRSCSx_E", resultsLoc_pheno$model, fixed=T)],3) ))
#resultsLoc_pheno = resultsLoc_pheno[-grep("_EUR_PRSCSx_E", resultsLoc_pheno$model, fixed=T),] 


# primary+adjunct Ldpred2 + shaPRS
print(paste0("primary+adjunct LDpred2 + ShaPRS r2: ", signif(resultsLoc_pheno$r2[grep("_shaPRS_LDpred2_CSxsubset_E", resultsLoc_pheno$model, fixed=T)],3), " / AUC: ", signif(resultsLoc_pheno$AUC[grep("_shaPRS_LDpred2_CSxsubset_E", resultsLoc_pheno$model, fixed=T)],3) ))
#resultsLoc_pheno = resultsLoc_pheno[-grep("_shaPRS_LDpred2_CSxsubset_E", resultsLoc_pheno$model, fixed=T),] 

# primary+adjunct PRS-CS +shaPRS-unweighted
print(paste0("primary+adjunct PRS-CS + ShaPRS-unweighted r2: ", signif(resultsLoc_pheno$r2[grep("_shaPRS_EUR_PRSCS_E", resultsLoc_pheno$model, fixed=T)],3), " / AUC: ", signif(resultsLoc_pheno$AUC[grep("_shaPRS_EUR_PRSCS_E", resultsLoc_pheno$model, fixed=T)],3) ))
#resultsLoc_pheno = resultsLoc_pheno[-grep("_shaPRS_EUR_PRSCS_E", resultsLoc_pheno$model, fixed=T),] 



####################################
# we no longer care about the 'weighted' ones as those did not work
# primary+adjunct PRS-CS +shaPRS
#print(paste0("primary+adjunct PRS-CS + ShaPRS r2: ", signif(resultsLoc_pheno$r2[grep(".PRSCS_E", resultsLoc_pheno$model, fixed=T)],2), " / AUC: ", signif(resultsLoc_pheno$AUC[grep(".PRSCS_E", resultsLoc_pheno$model, fixed=T)],2) ))
#resultsLoc_pheno = resultsLoc_pheno[-grep(".PRSCS_E", resultsLoc_pheno$model, fixed=T),] 


# primary+adjunct PRS-CSx + ShaPRS
#print(paste0("primary+adjunct PRS-CSx + ShaPRS r2: ", signif(resultsLoc_pheno$r2[grep(".PRSCSx_ShaPRS_E", resultsLoc_pheno$model, fixed=T)],2), " / AUC: ", signif(resultsLoc_pheno$AUC[grep(".PRSCSx_ShaPRS_E", resultsLoc_pheno$model, fixed=T)],2) ))
#resultsLoc_pheno = resultsLoc_pheno[-grep(".PRSCSx_ShaPRS_E", resultsLoc_pheno$model, fixed=T),] 


# one remained: _EUR_PRSCSx_ShaPRS_E
# primary+adjunct PRS-CSx + shaPRS-unweighted, IE the unweighted version of the above, this is unused
#
########################################################################

# create data frame of the relevant vars for AUC
# colorblind palette from: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

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

  prsResults_r2 = cbind(prsResults_r2,   c(primary_ldpred2_r2, primary_prscs_r2, pooledData_prscsx_r2, pooledData_prscsx_unweighted_r2, pooledData_shaprs_ldpred2_unweighted_r2, pooledData_shaprs_prcs_unweighted_r2 ))
  prsResults_AUC = cbind(prsResults_AUC,   c(primary_ldpred2_AUC, primary_prscs_AUC, pooledData_prscsx_AUC, pooledData_prscsx_unweighted_AUC, pooledData_shaprs_ldpred2_unweighted_AUC, pooledData_shaprs_prcs_unweighted_AUC ))
  
  prsResults_r2_sd = cbind(prsResults_r2_sd,   c(primary_ldpred2_r2_sd, primary_prscs_r2_sd, pooledData_prscsx_r2_sd, pooledData_prscsx_unweighted_r2_sd, pooledData_shaprs_ldpred2_unweighted_r2_sd, pooledData_shaprs_prcs_unweighted_r2_sd ))
  prsResults_AUC_sd_high = cbind(prsResults_AUC_sd_high,   c(primary_ldpred2_AUC_sd_high, primary_prscs_AUC_sd_high, pooledData_prscsx_AUC_sd_high, pooledData_prscsx_unweighted_AUC_sd_high, pooledData_shaprs_ldpred2_unweighted_AUC_sd_high, pooledData_shaprs_prcs_unweighted_AUC_sd_high ))
  prsResults_AUC_sd_low = cbind(prsResults_AUC_sd_low,   c(primary_ldpred2_AUC_sd_low, primary_prscs_AUC_sd_low, pooledData_prscsx_AUC_sd_low, pooledData_prscsx_unweighted_AUC_sd_low, pooledData_shaprs_ldpred2_unweighted_AUC_sd_low, pooledData_shaprs_prcs_unweighted_AUC_sd_low ))
  
}

#rownames_trait = c("primary\nLDpred2", "primary\nPRS-CS", "pooled\nPRS-CSx", "pooled PRS-CSx\n(unweighted)", "pooled shaPRS+LDpred2\n(unweighted)", "pooled shaPRS+PRS-CS\n(unweighted)")
rownames_trait = c("LDpred2", "PRS-CS", "PRS-CSx", "PRS-CSx\n(unweighted)", "shaPRS+LDpred2\n(unweighted)", "shaPRS+PRS-CS\n(unweighted)")

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




# barchart with par
createPlot = function(prsResults, prsResults_sd_low, prsResults_sd_high, ylabel, filn="r2", showHeight =T) {
  
col <- safe_colorblind_palette[1: nrow(prsResults)]

if(showHeight == F) { # remove height from AUC results, as that is a continuous trait
  prsResults = subset(prsResults, select = -which(colnames(prsResults) == "height"))
  prsResults_sd_low = subset(prsResults_sd_low, select = -which(colnames(prsResults_sd_low) == "height"))
  prsResults_sd_high = subset(prsResults_sd_high, select = -which(colnames(prsResults_sd_high) == "height"))
  
}


prsResults.matrix <- as.matrix(prsResults)
# create matrix of errors
prsResults.error_low = as.matrix(prsResults_sd_low)
prsResults.error_high = as.matrix(prsResults_sd_high)
minY = 0
if(showHeight == F) {
  minY=0.5 # for AUCs we only care  whats above 0.5
}


filename=paste0(outputLoc,filn, ".png")
png(filename, width=1280 , height=480, res =128)
par(mar=c(2.0,6,0.5,2)+0.1,mgp=c(5,1,0)) # par(mar = c(bottom, left, top, right))
# create bar plot
bp <- barplot(prsResults.matrix, xpd = F,
              beside=T,
              col=col,
              ylim=c(minY,max(prsResults) + max(prsResults) * 0.1),
              xlab="",
             # ylab=parse(text=paste(ylabel)),
              legend.text=row.names(prsResults.matrix),
              args.legend=list("topright",bty="n", ncol = 1, y.intersp=1.8),
              las=1,
              xlim = c(0, 45),
              cex=1.2,
              cex.lab = 1.2,
              cex.axis=1.2)
title(ylab=parse(text=paste(ylabel)), line=3.7, cex.lab=1.5)

if(showHeight == F) { # dont plot error bars for height runs, IE when it is r^2, as those are too small
# add error bars
arrows(bp, prsResults.error_low, bp, prsResults.error_high,
       code = 3, angle = 90, length = 0.15)
}

# Restore default clipping rect
#par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()

}

createPlot(prsResults_r2, prsResults_r2_sd, prsResults_r2_sd, "r^2", filn="r2", T)
createPlot(prsResults_AUC,prsResults_AUC_sd_low, prsResults_AUC_sd_high, "AUC", filn="auc", F)

