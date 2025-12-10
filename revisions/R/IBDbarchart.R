# loads and creates barchart for the AFR, and then it merges them


prsResults_r2= read.csv("C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/IBD_res_r2.csv", header = T,row.names = 1)
prsResults_r2_sd= read.csv("C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/IBD_res_sd.csv", header = T,row.names = 1)


outputLoc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/IBD_res_"

# createPlot() is loaded from ElenaModelParser.R
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#332288", "#AA4499", "#DDCC77", "#117733")

createPlot(prsResults_r2, prsResults_r2_sd, prsResults_r2_sd, "r^2", filn="r2", T, plotWidth=2070/2, xlimit = c(0,18) )


prsResults_AUC= read.csv("C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/IBD_res_AUC.csv", header = T,row.names = 1)
prsResults_AUC_LOW= read.csv("C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/IBD_res_AUC_LOW.csv", header = T,row.names = 1)
prsResults_AUC_HIGH= read.csv("C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/IBD_res_AUC_HIGH.csv", header = T,row.names = 1)


createPlot(prsResults_AUC, prsResults_AUC_LOW, prsResults_AUC_HIGH, "AUC", filn="auc", F, plotWidth=2070/2, xlimit = c(0,18), minY=0.5)

createPlot(prsResults_AUC, prsResults_AUC_LOW, prsResults_AUC_HIGH, "AUC", filn="auc2", F, plotWidth=2070, xlimit = c(1, 13.7), minY=0.5, cexV = 1.5, showValues =T)


###################
# get perc diffs
percDiff = function (mean1, mean2, othermethodName) {
  perc = round( (mean1 - mean2) / ( (mean1 + mean2)/2 ) * 100)
  
  print(paste0("shaPRS is better than ", othermethodName, " by: ", perc,"%" ))
}

# CD
percDiff(prsResults_r2$CD[5], prsResults_r2$CD[1], "proximal (CD)") # 4%
percDiff(prsResults_r2$CD[5], prsResults_r2$CD[2], "meta (CD)") # 12%
percDiff(prsResults_r2$CD[5], prsResults_r2$CD[3], "SMTPred (CD)") # 7% 
percDiff(prsResults_r2$CD[5], prsResults_r2$CD[4], "MTAG (CD)") # 11%
#UC
percDiff(prsResults_r2$UC[5], prsResults_r2$UC[1], "proximal (UC)") # 22%
percDiff(prsResults_r2$UC[5], prsResults_r2$UC[2], "meta (UC)") # 6%
percDiff(prsResults_r2$UC[5], prsResults_r2$UC[3], "SMTPred (UC)") # 10%
percDiff(prsResults_r2$UC[5], prsResults_r2$UC[4], "MTAG (UC)") # 39%


##########################
# merge them:
# GenerateImage() relies on ImageCombiner.R
library(magick)


baseloc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/"

imagesPerRow= 1 # how many images per row
outName="IBD_rev2.png"

image1 <- image_read(paste0(baseloc,"manhattan_lFDR_QAdjusted_manhattan.png"))
image2 <- image_read(paste0(baseloc,"IBD_res_auc2.png"))
images = list(image1, image2)
labels = list("a", "b")
img = GenerateImage(images, labels, imagesPerRow)
image_write(img, paste0(baseloc,outName))


# OLD plot, that had 3 images with the r2 too
# imagesPerRow= 2 # how many images per row
# outName="IBD_partial.png"
# 
# image1 <- image_read(paste0(baseloc,"IBD_res_r2.png"))
# image2 <- image_read(paste0(baseloc,"IBD_res_auc.png"))
# 
# 
# images = list(image1, image2)
# labels = list("b", "c")
# img = GenerateImage(images, labels, imagesPerRow)
# image_write(img, paste0(baseloc,outName))
# 
# 
# imagesPerRow= 1
# image_loaded <- image_read(paste0(baseloc,outName))
# 
# 
# image3 <- image_read(paste0(baseloc,"manhattan_lFDR_QAdjusted_manhattan.png"))
# 
# labels = list("a", "")
# images = list(image3,image_loaded )
# outName="IBD_pub_rev.png"
# img = GenerateImage(images, labels, imagesPerRow)
# image_write(img, paste0(baseloc,outName))
