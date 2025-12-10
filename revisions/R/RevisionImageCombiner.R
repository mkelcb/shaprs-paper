
# GenerateImage() relies on ImageCombiner.R

library(magick)
baseloc="C:/Users/mk23/GoogleDrive_Cam/0Publications/shaPRS/Figs/revisions/sims/"

##########################
# SIMS 1 maintext

image1 <- image_read(paste0(baseloc,"maintext/1000_h2_0.5sims_eval_ukbb_revision_subset.png"))
image2 <- image_read(paste0(baseloc,"maintext/1000_h2_0.5sims_eval_ukbb_revision_subset_extra.png"))

# all images have the "r^2" label misasligned, so I crop and paste that area of the image shifting it by 20 pixels to a better aesthetic result
r2_bit <- image_crop(image1, "92x37+483") # this cuts the r2 bit

sims1 =  image_composite(image1,r2_bit, offset ="+503" , operator = "Over" )
sims2 =  image_composite(image2,r2_bit, offset ="+503" , operator = "Over" )

outName="sims_publish.png"
imagesPerRow= 2 # how many images per row
images = list(sims1, sims2)
labels = list("a", "b")

img = GenerateImage(images, labels, imagesPerRow)
image_write(img, paste0(baseloc,"pub/",outName))


###################################
# SIMULATION 2
outName="sims_publish2.png"
imagesPerRow= 2 # how many images per row
sims1 <- image_read(paste0(baseloc,"causals_1000_h2_0.5_rG_0.5_het_regular_A50_B50_sizefull_revision.png"))
sims2 <- image_read(paste0(baseloc,"causals_1000_h2_0.5_rG_0.5_het_extra_A50_B50_sizefull_revision.png"))
images = list(sims1, sims2)
labels = list("a", "b")

img = GenerateImage(images, labels, imagesPerRow)
image_write(img, paste0(baseloc,"pub/",outName))



##################################
# All Sim Data:
image1 <- image_read(paste0(baseloc,"1000_h2_0.5sims_eval_ukbb_revision.png"))

sims1 =  image_composite(image1,r2_bit, offset ="+503" , operator = "Over" )

image_write(sims1, paste0(baseloc,"pub/","sims_eval_ukbb.png"))


##########################
# SIMS aux figure


image1 <- image_read(paste0(baseloc,"1000_h2_0.25sims_eval_ukbb_revision_subset.png"))
image2 <- image_read(paste0(baseloc,"1000_h2_0.75sims_eval_ukbb_revision_subset.png"))
image3 <- image_read(paste0(baseloc,"3000_h2_0.5sims_eval_ukbb_revision_subset.png"))
image4 <- image_read(paste0(baseloc,"5000_h2_0.5sims_eval_ukbb_revision_subset.png"))

sims1 =  image_composite(image1,r2_bit, offset ="+503" , operator = "Over" )
sims2 =  image_composite(image2,r2_bit, offset ="+503" , operator = "Over" )
sims3 =  image_composite(image3,r2_bit, offset ="+503" , operator = "Over" )
sims4 =  image_composite(image4,r2_bit, offset ="+503" , operator = "Over" )

outName="sims_publish_SFig4.png"
imagesPerRow= 2 # how many images per row
images = list(sims1, sims2, sims3, sims4)
labels = list("a", "b", "c", "d")

img = GenerateImage(images, labels, imagesPerRow)
image_write(img, paste0(baseloc,"pub/",outName))
