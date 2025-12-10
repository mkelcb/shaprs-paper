


options(error=traceback)
# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 3) {stop("not enough Arguments received")}


# arguments
PCsLoc= args[1]
existingIndisLoc= args[2]
outFileLoc = args[3]

# PCsLoc= "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/customLDRef/QC.csv"
# existingIndisLoc= "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt"
# outFileLoc = "C:/0Datasets/ancestrypca/UKBB_PCA_indis"


PC_UKBB_Id <- read.table(PCsLoc, header=T)

existingIndis <- read.table(existingIndisLoc, header=F)

PC_UKBB_Id_existing = merge (existingIndis, PC_UKBB_Id, by.x = "V2", by.y = "eid")

PC_UKBB_Id = PC_UKBB_Id_existing[,c(-2,-3)]
colnames(PC_UKBB_Id) = rep(paste0("V",1:ncol(PC_UKBB_Id)))



PC_UKBB = PC_UKBB_Id[-1]




all_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)

# compute square distances from the pop center for each PC
all_sq_dist <- apply(all_centers[-1], 1, function(one_center) {
  rowSums(sweep(PC_UKBB, 2, one_center, '-')^2)
})


thr_sq_dist <- max(dist(all_centers[-1])^2) * 0.002 / 0.16
group <- apply(all_sq_dist, 1, function(x) {
  grp <- NA
  ind <- which.min(x)
  if (isTRUE(x[ind] < thr_sq_dist)) {
    grp <- all_centers$Ancestry[ind]
    # We used a more stringent cutoff for the Ashkenazi group
    if (grp == "Ashkenazi" && x[ind] > 12.5^2) grp <- NA
  }
  grp
})

table(group, exclude = NULL)

# Ashkenazi      Caribbean          China          India           Iran
# 2348           2477           1815           6370           1214
# Italy        Nigeria         Poland United Kingdom           <NA>
# 6490           3937           4127         412209          10858

# exclude NAs, IE mixed people who could not  be assigned to any one group
PC_UKBB_Id_group = cbind.data.frame(PC_UKBB_Id$V1, 1:length(group),group)
colnames(PC_UKBB_Id_group) = c("IID", "index","pop")


AFR_indis= PC_UKBB_Id_group$IID[which(PC_UKBB_Id_group$pop == "Nigeria" | PC_UKBB_Id_group$pop == "Caribbean")]


write.table(AFR_indis, outFileLoc, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
print(paste("written ", length(AFR_indis), " AFR indis to: ", outFileLoc))

