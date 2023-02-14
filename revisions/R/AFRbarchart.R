# loads and creates barchart for the AFR, and then it merges them

prsResults_r2= read.csv("C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/AFR_res.csv", header = T,row.names = 1)

prsResults_r2_sd= read.csv("C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/AFR_res_sd.csv", header = T,row.names = 1)


outputLoc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/AFR_res"

# createPlot() is loaded from ElenaModelParser.R
createPlot(prsResults_r2, prsResults_r2_sd, prsResults_r2_sd, "r^2", filn="r2", T)


##########################
# merge them:
# GenerateImage() relies on ImageCombiner.R
library(magick)

outputLoc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/AFR_res"
baseloc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/final/"
# CROSS ANCESTRY BARPLOTS
imagesPerRow= 1 # how many images per row
outName="crossAncestryBarplot_publish_rev.png"

#image1 <- image_read(paste0(baseloc,"auc.png"))
image2 <- image_read(paste0(baseloc,"r2.png"))
image3 <- image_read(paste0(outputLoc,"r2.png"))

images = list(image2,image3) # image1, 
labels = list("a", "b") # , "c"
img = GenerateImage(images, labels, imagesPerRow)
image_write(img, paste0("C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/",outName))



