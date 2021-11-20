# converts an LDpred2 LDmatrix into a PRS-CS formatted matrix

#  install.packages("BiocManager")
#  BiocManager::install("rhdf5")
  
  # Call the R HDF5 Library
  library("rhdf5")
  library(Matrix)
  
  args = commandArgs(trailingOnly = TRUE)
  print("Arguments received:")
  print(args)
  if(length(args) < 3) {stop("not enough Arguments received")}
  
  PRSCSRefLoc = args[1]
  LDpred2Loc = args[2]
  outputLoc = args[3]
  
  dir.create(outputLoc)
  
# debug vars
# PRSCSRefLoc = "C:/0Datasets/shaPRS/PRSCSref/"
# LDpred2Loc = "C:/0Datasets/shaPRS/PRSCSref/ldpred/origEUR/"
#  LDpred2Loc = "C:/0Datasets/shaPRS/PRSCSref/ldpred/"
# outputLoc = "C:/0Datasets/shaPRS/PRSCSref/Out/"

#


  i=21
  j=1
  # write conversion script:
  numTotal = 0
  numFound = 0
  
  infoFileLoc=paste0(outputLoc,"snpinfo_ukbb_hm3") # write the first line of the output SNP info file
  write.table(rbind(c("CHR","SNP","BP","A1","A2","MAF")), file=infoFileLoc,  sep = "\t", row.names = F, col.names = F, quote = FALSE)
  
  # load the LDpred SNP info file
  LDpred2_Info = readRDS(file = paste0(LDpred2Loc,"map.rds") )
  
  # load the original PRS-CS SNP info file
  PRCSInfo = read.table(paste0(PRSCSRefLoc,"snpinfo_ukbb_hm3"), header = T)
  
  # the PRS-CS LD format is aligned against the minor allele, IE A1 there is always the minor
  # whereas the UKBB a1 is NOT always the minor. To get the correct orientation, we need to match
  # the A1/A2 alleles of the UKBB data

  # find the instances where the AF in the UKBB is not the MAF
  flipped_indices = LDpred2_Info$af_UKBB>=0.5
  flipped_mask =  rep(1, nrow(LDpred2_Info) )
  flipped_mask[flipped_indices] = -1 # set
  
  # create A1/A2
  LDpred2_Info$A1 = LDpred2_Info$a1
  LDpred2_Info$A2 = LDpred2_Info$a0
  LDpred2_Info$MAF = LDpred2_Info$af_UKBB
  

  LDpred2_Info$A1[flipped_indices] = LDpred2_Info$a0[flipped_indices] # swap the A1s with the A2s
  LDpred2_Info$A2[flipped_indices] = LDpred2_Info$a1[flipped_indices] # swap the A2s with the A1s
  LDpred2_Info$MAF[flipped_indices] = 1 - LDpred2_Info$MAF[flipped_indices]
  
  convertDSCToDense = function(LDpredChromLDMatrix, numparts = 3) {
    # convert dsc to dense matrix in a non-crashing way
    firstHalf = floor(ncol(LDpredChromLDMatrix)/numparts)
    start = 0;
    end = 0
    r1 = NULL
    for (i in 1:numparts) {
      start = end
      if (i == numparts) { # last one has to go to the end
        end = ncol(LDpredChromLDMatrix)
      } else {
        end = end + firstHalf
      }
      
      r1_p1 = as.matrix(LDpredChromLDMatrix[,(start+1):end])
      r1 = cbind(r1,r1_p1)
      remove(r1_p1)
    }
    return(r1)
  }
  
  
  for (i in 1:22) { # loop each chrom
    ldpredChromLoc = paste0(LDpred2Loc,"LD_chr",i,".rds")
    if ( file.exists(ldpredChromLoc)) {

      # I) Load LDpred2 ref into memory
      # 1) load the chrom
      LDpredChromLDMatrix = readRDS(file =  ldpredChromLoc)
      LDpredChromLDMatrix = convertDSCToDense(LDpredChromLDMatrix, 10) # if LD mat was sparse this could cause error
      # 2. grab the RSids from the map for the SNPS on this chrom, each LD mat has a potentiall different subset of SNPs
      #LDpred2_Info_chrom = readRDS(file = paste0(LDpred2Loc,"LD_chr",i,"_map.rds") ) # this is guaranteed to be the same order as the LDpredChromLDMatrix
      chrom_indices = which(LDpred2_Info$chr == i) # get the indices to subset the data to this chrom only
      LDpred2_Info_chrom = LDpred2_Info[ chrom_indices,] # this is guaranteed to be the same order as the LDpredChromLDMatrix
      PRCSInfo_chr = PRCSInfo[PRCSInfo$CHR == i,]
      PRCSInfo_chr$PRSCSid = 1:nrow(PRCSInfo_chr)

      LDpred2_Info_chrom$id = 1:nrow(LDpred2_Info_chrom) # save away the original ordering of the SNPs in the LDpred2 corr matrix, (indices refer to within chrom!)
      
      
      # II) Load the PRS-CS data
      hdf5Loc =paste0(PRSCSRefLoc,"ldblk_ukbb_chr",i,".hdf5")
      hdf5Loc_shaPRS = paste0(outputLoc,"ldblk_ukbb_chr",i,".hdf5")
      # copy it:
      file.copy(hdf5Loc, hdf5Loc_shaPRS, overwrite = TRUE)
      
  
      blocks = c()
      h5structure = h5ls(hdf5Loc)
      for(j in 1:length(h5structure$group)) {  # find each block
        if(h5structure$otype[j] == "H5I_GROUP") { # only care about groups
       # print(paste0(h5structure$name[j]))
          blocks = c(blocks,h5structure$name[j])
        }
      }
    j=1
    
      print(paste0("found ", length(blocks), " blocks for chrom",i ))
      for (j in 1:length(blocks)) { # iterate each block
        print(paste0("processing block ", j, " for chrom",i,"..." ))
        # get their rsid SNP list for that block from PRS-CS
        blockData <- h5read(hdf5Loc, name = blocks[j])
        numTotal = numTotal + length(blockData$snplist)
        
        # II) match PRS-CS and Ldpred2 data
        # get a submatrix from the LDpred2 overall chrom matrix that matches the PRS-CS block
        blockData_snplist = cbind(blockData$snplist, 1:length(blockData$snplist)) # create dummy df
 
        
       
        # merge it with the ldpred2 chrom to get their intersection
        SNPsFoundInLDpred2Ref = merge(LDpred2_Info_chrom, blockData_snplist, by.x = "rsid", by.y="V1" )
        SNPsFoundInLDpred2Ref$V2 = as.numeric(SNPsFoundInLDpred2Ref$V2) # V2 is the index of the SNPs within the block for PRS-CS (which is WITHIN the Chrom)

        numFound = numFound+nrow(SNPsFoundInLDpred2Ref)
        # order them back to be ascending
        SNPsFoundInLDpred2Ref <- SNPsFoundInLDpred2Ref[order(SNPsFoundInLDpred2Ref$id),] 
        # also merge the map data with the original PRSCS info
        PRCSInfo_chr_LDpredFound = merge(PRCSInfo_chr,SNPsFoundInLDpred2Ref, by.x="SNP", by.y="rsid")
        
        # extract the relevant corr sub  matrix from LDpred2 crhom data
        LDpredChromLDMatrix_matched = LDpredChromLDMatrix[SNPsFoundInLDpred2Ref$id,SNPsFoundInLDpred2Ref$id]
 

        
        #LDpredChromLDMatrix_matched[1:5, 1:5]
        
        # create a mask of the flipped alleles: this is just the outer product of the above vector
        flipped_mask_matched = flipped_mask[chrom_indices][SNPsFoundInLDpred2Ref$id]
        flippedMat=outer(flipped_mask_matched,flipped_mask_matched,"*")  # create a matrix that can be used to flip correlations via elementwise multiplications
        LDpredChromLDMatrix_matched = LDpredChromLDMatrix_matched * flippedMat # apply the flipping to the LD mat
   
        ## replace the original matrix at the indies where we have new data
        blockDataReplaced = blockData$ldblk
        blockDataReplaced[SNPsFoundInLDpred2Ref$V2,SNPsFoundInLDpred2Ref$V2] = LDpredChromLDMatrix_matched
        
       # cor(blockDataReplaced[3,],blockData$ldblk[3,])
      #  plot(blockDataReplaced[3,],blockData$ldblk[3,])
  
        # III) write new hdf5 file: https://bioconductor.riken.jp/packages/devel/bioc/vignettes/rhdf5/inst/doc/rhdf5.html
        # create new replacement object
        blockData_shaPRS = list()
        blockData_shaPRS = blockData
        blockData_shaPRS$ldblk = blockDataReplaced
        blockData_shaPRS$snplist = blockData$snplist
       
        #######################################################
        # check the alignment of the 2 LD matrices
        #PRSCSorigBlock = blockData$ldbl
        #PRSCSorigBlock = PRSCSorigBlock[SNPsFoundInLDpred2Ref$V2,SNPsFoundInLDpred2Ref$V2]
        #LDpredChromLDMatrix_matched[1:5,1:5]
        #PRSCSorigBlock[1:5,1:5]
        
       # LDpredChromLDMatrix_matched[1:10,5]
       # PRSCSorigBlock[1:10,5]
        
       # cor(LDpredChromLDMatrix_matched[1:20,1], PRSCSorigBlock[1:20,1] )
       # cor(LDpredChromLDMatrix_matched[1,], PRSCSorigBlock[1,]) # -0.2936863 # as many alleles are flipped
       # cor(  abs(LDpredChromLDMatrix_matched[1,]),   abs(PRSCSorigBlock[1,])) # 0.9999304
        #######################################################
        
       
        # first need to delete the old object, as we cannot overwrite the same object with different dimensions
        h5delete(hdf5Loc_shaPRS, name = blocks[j])
        h5write(blockData_shaPRS, hdf5Loc_shaPRS, blocks[j], native  = T)
        # verify the update
        # h5structure_shaPRS = h5ls(hdf5Loc_shaPRS)
        # blockData_new <- h5read(hdf5Loc_shaPRS, name = blocks[j])
        
  
        # append SNP list to "snpinfo_ukbb_hm3"
        # use the A1/A2 and MAFs from the LDPred2 data to extract signature:
        # CHR	SNP	BP	A1	A2	MAF
       # subsetSNPData= SNPsFoundInLDpred2Ref[, c(2,1,3,11,12,13)]
      # head(subsetSNPData)
      #  head(PRCSInfo_chr)
        # overwrite the found indices in the SNP info file
        PRCSInfo_chr[PRCSInfo_chr_LDpredFound$PRSCSid,] = PRCSInfo_chr_LDpredFound[,c(2,1,3,17,18,19)]
        
      #   write.table(subsetSNPData, file=infoFileLoc, append=TRUE,  sep = "\t", row.names = F, col.names = F, quote = FALSE)
        #subsetSNPData$rsid[1:5]
        #blockData_shaPRS$snplist[1:5]
        # compare to the original SNP infos
       # PRCSInfo_chr_ordered = PRCSInfo_chr[SNPsFoundInLDpred2Ref$V2,]
        #PRCSInfo_chr_ordered = PRCSInfo_chr[PRCSInfo_chr_LDpredFound$PRSCSid,]
        #PRCSInfo_chr_ordered <- PRCSInfo_chr_ordered[order(PRCSInfo_chr_ordered$PRSCSid),] 
        
        #cor(subsetSNPData$MAF, PRCSInfo_chr_ordered$MAF) # 0.9998349
        #sum(subsetSNPData$A1 != PRCSInfo_chr_ordered$A1) # 0
        #sum(subsetSNPData$A2 != PRCSInfo_chr_ordered$A2) # 0
      }
      
      write.table(PRCSInfo_chr[,c(1,2,3,4,5,6)], file=infoFileLoc, append=TRUE,  sep = "\t", row.names = F, col.names = F, quote = FALSE)
      
    
      
    }
    
  }
  
  infoText=paste0("In total number of SNPs converted ", numFound, " / ", numTotal , " (",round(numFound/numTotal,2),")")
  cat(infoText, file=paste0(outputLoc,"_conversionInfo.txt"))
  
  print(infoText)
  print(paste0("written results for all chroms to ",outputLoc))
  
# check by position:(from: https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/fourier_ls-chr21.bed)
# this only found the same number of SNPs as by matching via rsid
#  block1SNPs = LDpred2_Info_chrom[LDpred2_Info_chrom$pos >= 9411243 & LDpred2_Info_chrom$pos <= 15950982,]
#  nrow(block1SNPs) # 222
   	 
  

  
  
  #    h5f$"blk_1/snplist" = blockData_shaPRS$snplist # this fails
  #    h5f$"blk_1/ldblk" = blockData_shaPRS$ldblk
  
  
  #  h5f$"blk_1/snplist" =  blockData$snplist # this works as it is the right dimensions
  
  
  
  
  # ?h5write
  # h5f = H5Fopen(hdf5Loc_shaPRS, native = T)
  # h5f
  # # h5f$"blk_1" = blockData_shaPRS
  # h5f$"blk_1" = blockData_shaPRS
  # h5f[eval(blocks[j])] = blockData_shaPRS
  # 
  # 
  # 
  # h5f[blocks[j]]
  # # close the file
  # H5Fclose(h5f)
  
  
  # h5closeAll()
  # blockDataSNPs_orig = cbind(blockData$snplist, 1:length(blockData$snplist) )
  # blockDataSNPs_new = cbind(blockData_new$snplist, 1:length(blockData_new$snplist) )
  # 
  # blockDataSNPs_orig_new_merged = merge(blockDataSNPs_orig,blockDataSNPs_new, by="V1")
  # blockDataSNPs_orig_new_merged$V2.x = as.numeric(blockDataSNPs_orig_new_merged$V2.x)
  # blockDataSNPs_orig_new_merged$V2.y = as.numeric(blockDataSNPs_orig_new_merged$V2.y)
  # 
  # 
  # blockData_orig_ld_matched = blockData$ldblk[blockDataSNPs_orig_new_merged$V2.x,blockDataSNPs_orig_new_merged$V2.x]
  # blockData_new_ld_matched = blockData_new$ldblk[blockDataSNPs_orig_new_merged$V2.y,blockDataSNPs_orig_new_merged$V2.y]
  # 
  # cor(blockData_orig_ld_matched[1,],blockData_new_ld_matched[1,])
  # cor(blockData_orig_ld_matched[2,],blockData_new_ld_matched[2,])
  # 
  # plot(blockData_orig_ld_matched[1,],blockData_new_ld_matched[1,])
  # 
  # SNPsFoundInLDpred2Ref = merge(LDpred2_Info_chrom, blockData_snplist, by.x = "rsid", by.y="V1" )
  # SNPsFoundInLDpred2Ref$V2 = as.numeric(SNPsFoundInLDpred2Ref$V2) # V2 is the index of the SNPs within the block (which is WITHIN the Chrom)
  # 
  # 
  