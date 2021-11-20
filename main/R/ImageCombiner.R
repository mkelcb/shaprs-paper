# Composites multiple images into a single image
#
# sources:
# https://cran.r-project.org/web/packages/magick/vignettes/intro.html#Combining
# https://www.ddrive.no/post/miracles-with-magick-and-bunny/

###########################

# install.packages("magick")
library(magick)



# I) figure out the total dimensions of the image:
GenerateImage = function(images, labels, imagesPerRow) {
  
  # to avoid crops/overlaps, need to find out the largest width/height among all images:
  largestHeight = 0
  largestWidth = 0
  for (i in 1:length(images)) {
    imageInfo = image_info(images[[i]])
    if(imageInfo$width > largestWidth ) {largestWidth = imageInfo$width}
    if(imageInfo$height > largestHeight ) {largestHeight = imageInfo$height}
  }
  
 # imageInfo = image_info(images[[1]])

  # Find out the total size of the image
  width = largestWidth *  imagesPerRow
  height = largestHeight * length(images) /imagesPerRow
  
  
  print(paste0("total image widht/height: ",width, "/", height, " (numRows: ",imagesPerRow,")"))
  # II) create blank image of the right size
  img = image_blank(width, height) # , color = "none", pseudo_image = "", defines = NULL
  print(img)
  
  
  # III) composit each image into the total
  i=1
  rowCount=0
  Xloc= 0
  Yloc=0
  for (i in 1:length(images)) {
    imageInfo = image_info(images[[i]])
    
    if(i != 1) {
      if(rowCount == imagesPerRow) {
        print(paste0("next row at i ", i))
        rowCount = 0
        Xloc = 0
        Yloc = Yloc + largestHeight # imageInfo$height
      } else {
        print(paste0("next col col i ", i))
        Xloc = Xloc + largestWidth # imageInfo$width
      }
    }
    print(paste0("adding image at: X/Y: ", Xloc , "/",Yloc))
    location=paste0("+",Xloc,"+",Yloc)
    img = image_composite(img, images[[i]], offset = location , operator = "Over" )
    img = image_annotate(img, labels[[i]], size = 32, location  = location, color = "black", weight = 700)
    
    rowCount = rowCount + 1
  }
  print(img)
  
  return(img)
}


baseloc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/final/"

##############################
# SIMULATION DATA
outName="sims_publish.png"
imagesPerRow= 2 # how many images per row
sims1 <- image_read(paste0(baseloc,"sims_eval_ukbb_subset.png"))
sims2 <- image_read(paste0(baseloc,"sims_eval_ukbb_subset_extra.png"))
images = list(sims1, sims2)
labels = list("a", "b")

img = GenerateImage(images, labels, imagesPerRow)
image_write(img, paste0(baseloc,outName))

##############################


##############################
# SIMULATION 2
outName="sims_publish2.png"
imagesPerRow= 2 # how many images per row
sims1 <- image_read(paste0(baseloc,"rG_0.5_het_regular_A50_B50_sizefull.png"))
sims2 <- image_read(paste0(baseloc,"rG_0.5_het_extra_A50_B50_sizefull.png"))
images = list(sims1, sims2)
labels = list("a", "b")

img = GenerateImage(images, labels, imagesPerRow)
image_write(img, paste0(baseloc,outName))

##############################




##############################
# IBD
outName="IBD_publish.png"
imagesPerRow= 1 # how many images per row
image1 <- image_read(paste0(baseloc,"CD_1_QAdjusted_manhattan.png"))
image2 <- image_read(paste0(baseloc,"CD_UC_PRS_PUB.png"))
images = list(image1, image2)
labels = list("a", "b")

img = GenerateImage(images, labels, imagesPerRow)
image_write(img, paste0(baseloc,outName))

##############################


##############################
# Cross ancestry diag
outName="crossAncestryDiag_publish2.png"
imagesPerRow= 3 # how many images per row
image1 <- image_read(paste0(baseloc,"CrossAncestryDiag/EUR_JAP_asthma_lFDR_hist.png"))
image2 <- image_read(paste0(baseloc,"CrossAncestryDiag/EUR_JAP_height_lFDR_hist.png"))

image3 <- image_read(paste0(baseloc,"CrossAncestryDiag/EUR_JAP_BRCA_lFDR_hist.png"))
image4 <- image_read(paste0(baseloc,"CrossAncestryDiag/EUR_JAP_CAD_lFDR_hist.png"))
image5 <- image_read(paste0(baseloc,"CrossAncestryDiag/EUR_JAP_T2D_lFDR_hist.png"))
image6 <- image_read(paste0(baseloc,"CD_diag/CD_1_lFDR_hist.png"))


images = list(image1, image2,image3, image4,image5, image6)
labels = list("a", "b","c", "d","e", "f")

img = GenerateImage(images, labels, imagesPerRow)
image_write(img, paste0(baseloc,outName))

##############################

# CROSS ANCESTRY BARPLOTS
imagesPerRow= 1 # how many images per row
outName="crossAncestryBarplot_publish.png"

image1 <- image_read(paste0(baseloc,"auc.png"))
image2 <- image_read(paste0(baseloc,"r2.png"))

images = list(image1, image2)
labels = list("a", "b")
img = GenerateImage(images, labels, imagesPerRow)
image_write(img, paste0(baseloc,outName))
