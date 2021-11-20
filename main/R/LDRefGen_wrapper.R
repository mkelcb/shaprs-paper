# Blends the LD reference panel for allchromosomes between 2 studies, relying on shaPRS summary stats and LD ref panels

##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################
options(error=traceback)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 7) {stop("not enough Arguments received")}

Pop1LDRefLoc = args[1]
Pop2LDRefLoc = args[2]
blendFactorLoc = args[3]
sumstatsLoc = args[4]
outputLoc = args[5]
produceBasicLDreftoo = (args[6] == '1')
shaPRSscriptLoc= args[7]

# load exteral functions 
source(shaPRSscriptLoc) # '/nfs/users/nfs_m/mk23/scripts/shaPRS.R'
library("bigsnpr")
library("pryr")


# Debug Vars
# Pop1LDRefLoc= "C:/0Datasets/LdPred2/"
# Pop2LDRefLoc="C:/0Datasets/LdPred2/"
# chromNum=21
# blendFactorLoc ="C:/0Datasets/shaPRS/GWAS_f1_lFDR_meta_SNP_lFDR"
# sumstatsLoc="C:/0Datasets/shaPRS/crossAncestry/EUR_JAP_asthma_SE_meta"
# outputLoc ="C:/0Datasets/shaPRS/0shaPRS/0ukbb/debug/1_pheBlend_sumstats_meta"
# produceBasicLDreftoo=T

##############################################################
#                         MAIN SCRIPT                        #
##############################################################


# 1. load Map data
pop1_map_rds = readRDS(file = paste0(Pop1LDRefLoc,"map.rds") )
pop2_map_rds = readRDS(file = paste0(Pop2LDRefLoc,"map.rds") )

# load raw blending factors and summary stats for the entire genome
blendingFactors= read.table(blendFactorLoc, header = T)
sumsData= read.table(sumstatsLoc, header = T)

dir.create(file.path(paste0(outputLoc,"/") ), showWarnings = FALSE)


if(produceBasicLDreftoo) {
  dir.create(file.path(paste0(outputLoc,"_basic/") ), showWarnings = FALSE)
}

# merge them
sumstatsDataAll =  merge(blendingFactors,sumsData,by.x = "SNP",by.y = "SNP")

# remove non numeric data
sumstatsDataAll <- sumstatsDataAll[!is.na(as.numeric(as.character(sumstatsDataAll$lFDR))),]
sumstatsDataAll <- sumstatsDataAll[!is.na(as.numeric(as.character(sumstatsDataAll$SE_A))),]
sumstatsDataAll <- sumstatsDataAll[!is.na(as.numeric(as.character(sumstatsDataAll$SE_B))),]


# now actually cast them to numeric
sumstatsDataAll$lFDR = as.numeric(as.character(sumstatsDataAll$lFDR ))
sumstatsDataAll$SE_A = as.numeric(as.character(sumstatsDataAll$SE_A ))
sumstatsDataAll$SE_B = as.numeric(as.character(sumstatsDataAll$SE_B ))

# align the summary for phe A and B
sumstatsDataAll = alignStrands(sumstatsDataAll)
#chromNum=21
# go through each chrom
for(chromNum in 1:22){

  
  # load the two chromosomes from each population
  pop1LDmatrix = readRDS(file = paste0(Pop1LDRefLoc,"LD_chr",chromNum,".rds") )
  pop2LDmatrix = readRDS(file = paste0(Pop2LDRefLoc,"LD_chr",chromNum,".rds") )
  

  # 2. grab the RSids from the map for the SNPS on this chrom, each LD mat has a potentiall different subset of SNPs
  pop1_chrom_SNPs = pop1_map_rds[ which(pop1_map_rds$chr == chromNum),] # this is guaranteed to be the same order as the pop1LDmatrix
  pop2_chrom_SNPs = pop2_map_rds[ which(pop2_map_rds$chr == chromNum),] # this is guaranteed to be the same order as the pop2LDmatrix
  pop1_chrom_SNPs$pop1_id = 1:nrow(pop1_chrom_SNPs)
  pop2_chrom_SNPs$pop2_id = 1:nrow(pop2_chrom_SNPs)
  
 
  # intersect the 2 SNP lists so that we only use the ones common to both LD matrices by merging them
  chrom_SNPs_df  <- merge(pop1_chrom_SNPs,pop2_chrom_SNPs, by = "rsid")

  # align the two LD matrices
  chrom_SNPs_df = alignStrands(chrom_SNPs_df, A1.x ="a1.x", A2.x ="a0.x", A1.y ="a1.y", A2.y ="a0.y")
  

  
  # subset sumstats data to the same chrom
  sumstatsData = sumstatsDataAll[which(sumstatsDataAll$CHR == chromNum ),]
  
  if(nrow(sumstatsData) > 0) {
    

  
  # merge sumstats with common LD map data
  sumstatsData  <- merge(chrom_SNPs_df,sumstatsData, by.x="rsid", by.y = "SNP")
  
  # remove duplicates
  sumstatsData = sumstatsData[ !duplicated(sumstatsData$rsid) ,]
  # use the effect alleles for the sumstats data with the effect allele of the LD mat
  # as we are aligning the LD mats against each other, not against the summary stats
  # we only use the lFDR /SE from the sumstats, which are directionless, so those dont need to be aligned
  sumstatsData$A1.x =sumstatsData$a1.x 
  sumstatsData$A1.y =sumstatsData$a1.y
  

  # make sure the sumstats is ordered the same way as the LD matrix: https://stackoverflow.com/questions/17878048/merge-two-data-frames-while-keeping-the-original-row-order
  sumstatsData = sumstatsData[order(sumstatsData$pop1_id), ] # it doesn't matter which matrix to use to order the sumstats as they are the same
  

  # subset the LD matrices to the SNPs we actualy have
  pop1LDmatrix = pop1LDmatrix[sumstatsData$pop1_id,sumstatsData$pop1_id]
  pop2LDmatrix = pop2LDmatrix[sumstatsData$pop2_id,sumstatsData$pop2_id]
  
  # generate the blended LD matrix
  cormat = LDRefBlend(pop1LDmatrix,pop2LDmatrix, sumstatsData)
  
  fileLoc= paste0(outputLoc,"/LD_chr",chromNum,".rds")
  saveRDS(cormat,file = fileLoc)
  print(paste0("written PRS specific LD mat to ",fileLoc ))
  

  # optionally, if the 'basicr' PRS mat was requested too, we do that as well
  if(produceBasicLDreftoo) {
    print("Also producing basic LD matrix")
    cormat = LDRefBlend_basic(pop1LDmatrix,pop2LDmatrix, sumstatsData)
    
    fileLoc= paste0(outputLoc,"_basic/LD_chr",chromNum,".rds")
    saveRDS(cormat,file = fileLoc)
    print(paste0("written basic PRS specific LD mat to ",fileLoc ))
  }
  
  } else {print(paste0("no variants on chrom", chromNum))}
  
  
  # also need to write out the list of SNPs that made it into the final subset, as after all LD matrices are done, we need to create a map.rds too
  write.table(sumstatsData$rsid, paste0(outputLoc,chromNum,"_snps"), sep = "\t", row.names = F, col.names = F, quote = FALSE)
  
  
  # map the final list of SNPs back to the original map file's indices
  map_rds_new = pop1_map_rds[which(pop1_map_rds$chr == chromNum),]
  map_rds_new2 = map_rds_new[which(map_rds_new$rsid %in% sumstatsData$rsid),] # match the first to the second
  
  fileLoc= paste0(outputLoc,"/LD_chr",chromNum,"_map.rds")
  saveRDS(map_rds_new2,file = fileLoc)
  print(paste0("written chr map to ",fileLoc ))

  
  mem_used()
}

# at the end concat all of the map files into a single file and write it to disk
all_map_rds = NULL
for(chromNum in 1:22){
  chr_map_rds = readRDS(file = paste0(outputLoc,"/LD_chr",chromNum,"_map.rds") )
  all_map_rds = rbind(all_map_rds,chr_map_rds)
}

fileLoc= paste0(outputLoc,"/map.rds")
saveRDS(all_map_rds,file = fileLoc)
print(paste0("written overall map to ",fileLoc ))

if(produceBasicLDreftoo) {
  fileLoc= paste0(outputLoc,"_basic/map.rds")
  saveRDS(all_map_rds,file = fileLoc)
}


