# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 4) {stop("not enough Arguments received")} 

famLoc = args[1]
outputLoc= args[2]
outputName = args[3]
num_samples = as.integer(args[4])

# Debug wars
# famLoc = "C:/softwares/Cluster/Imputation/Impute4/shapeit3/example/gwas.fam"
# outputLoc = "C:/softwares/Cluster/UKBB/"
# outputName = "myFold_"
# num_samples =5




# 1. load .fam file
famFile = as.matrix( read.table(famLoc) )
indi_indices = seq(1,nrow(famFile))

set.seed(42)

for(i in 1:num_samples) {
  #  grab the relevant samples
  counter = i
  bootstrap_indices = sample(indi_indices, replace = T)
  foldSamples = famFile[bootstrap_indices,]
  
  validIndices = setdiff(indi_indices, bootstrap_indices)
  duplicates = duplicated(bootstrap_indices)
  uniques = !duplicates # this actually includes ALL individuals, even the duplicate ones ,but just for once

  
  foldSamples_uniques = foldSamples[uniques,]

  foldSamples_valid = famFile[validIndices,] # these are NOT on either the unique OR the duplicate list
  # because plink does not allow people with duplicate IDs, we have to separate out the duplicates
  # and the unique ones into 2 lists...
  perc_duplicates =  sum ( duplicated(bootstrap_indices) ) /   length(indi_indices)
  
  
  
  filename = paste(outputLoc,"/", outputName, counter, ".fam" , sep="")
  write.table(foldSamples, filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  print(paste("written bootstrap sample to:", filename , " / duplicate individuals:", round(perc_duplicates * 100, 0), "%"))
  
  

  filename = paste(outputLoc,"/", outputName, counter, "_uniques", ".fam" , sep="")
  write.table(foldSamples_uniques, filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  filename = paste(outputLoc,"/", outputName, counter, "_valid", ".fam" , sep="")
  write.table(foldSamples_valid, filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  # PLINK cannot export files where we have the same individual 2x
  # the resampling process allows the sample to be selected more than twice, so even the duplicates list will have duplicates
  # so we need to export out a list for each 'level' of duplicates, and rename them separately, and merge them back together in a loop
  # give new IDs to the duplicates so they can be merged back into the same file without plink noticing
  
  #duplicateNames = seq(1, sum ( duplicated(bootstrap_indices) )  )
  #maxNum = max(duplicateNames)
  #numDigits = floor( log10(maxNum) ) +1
  #duplicateNames = sprintf(paste("%0",numDigits,"d",sep=""),duplicateNames)
  #duplicateNames = paste("duplicate",duplicateNames , sep="")
  
  #duplicatesExist= TRUE
  duplicate_indices = bootstrap_indices[duplicates]
  dupeCounter = 0
  lastDupe = 0
  
  
  #filename = paste(outputLoc,"/", outputName,"duplicate_indices", counter, sep="")
  #write.table(duplicate_indices, filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  #filename = paste(outputLoc,"/", outputName,"bootstrap_indices", counter, sep="")
  #write.table(bootstrap_indices, filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  #foldSamples = as.matrix( read.table("C:/softwares/Cluster/UKBB/GWAS_f3.fam") )
  #duplicate_indices = unlist( read.table("C:/softwares/Cluster/UKBB/GWAS_fduplicate_indices3") )
  #bootstrap_indices =  unlist( read.table("C:/softwares/Cluster/UKBB/GWAS_fbootstrap_indices3") )
  #dupe_orig = read.table("C:/softwares/Cluster/UKBB/bootstrap/GWAS_f1_dupes_level1.fam") 
  #dupenewName = read.table("C:/softwares/Cluster/UKBB/bootstrap/GWAS_f1_dupes_level1_newname.fam") 
  #duplicate_colname = duplicated(dupe_orig$V1)
  #duplicated_fam = dupe_orig[duplicate_colname,]
  
  while( length(duplicate_indices) > 0 ) {
    print(paste("duplicates found:",length(duplicate_indices) ))
    dupeCounter = dupeCounter +1
    duplicates = duplicated(duplicate_indices)
    unique_duplicates_thisLevel = !duplicates
    unique_duplicates_thisLevel_indices = duplicate_indices[unique_duplicates_thisLevel]
    duplicate_indices = duplicate_indices[duplicates]
    
 
    foldSamples_duplicates = rbind( famFile[unique_duplicates_thisLevel_indices,] ) # rbind, otherwise if we have only 1 index, it will turn this into a column not a row 
    filename = paste(outputLoc,"/", outputName, counter, "_dupes_level",dupeCounter, ".fam" , sep="")
    write.table(foldSamples_duplicates, filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    ###################################
    # foldSamples_duplicates
    # duplicate_colname = duplicated(foldSamples_duplicates[,1])
    # duplicated_fam = foldSamples_duplicates[duplicate_colname,] 
    # length(sum(duplicate_colname) )
    # length(duplicated_fam)
    # duplicate_colname
    # typeof(duplicate_colname)
    ###################################
    duplicateNames = seq( (lastDupe +1), ( (lastDupe ) +  length(unique_duplicates_thisLevel_indices) )  )
    duplicateNames = paste("duplicate",duplicateNames , sep="")
    foldSamples_duplicates[,1] = foldSamples_duplicates[,2] =  duplicateNames
    filename = paste(outputLoc,"/", outputName, counter, "_dupes_level",dupeCounter,"_newname", ".fam" , sep="")
    write.table(foldSamples_duplicates, filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
    print(paste("written",length(duplicateNames)," unique dupes to:",filename ))

    lastDupe = lastDupe +  length(unique_duplicates_thisLevel_indices)
    typeof(foldSamples_duplicates)
    
 
    #length(unique_duplicates_thisLevel_indices)
    #length(duplicate_indices)
    #View(as.matrix(unique_duplicates_thisLevel_indices))
    #View(as.matrix(duplicate_indices))
    
    #length(duplicate_indices)
    
  }
# need to also save how many duplicate levels we had so we can loop these in bash
filename = paste(outputLoc,"/", outputName, counter, "_numDupes", sep="")
write.table(paste("numDupeLevels=",dupeCounter,sep=""), filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
#  duplicateNames = seq( (lastDupe +1), length(unique_duplicates_thisLevel_indices)  )

  
  
  #foldSamples_duplicates[,1] = foldSamples_duplicates[,2] =  duplicateNames

#  filename = paste(outputLoc,"/", outputName, counter, "_dupes_newname", ".fam" , sep="")
#  write.table(foldSamples_duplicates, filename, sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}




