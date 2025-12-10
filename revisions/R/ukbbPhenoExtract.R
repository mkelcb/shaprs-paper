# extracts a phenotype into a .pheno format, from a ukbb tab file, given a target pheno ID
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 3) {stop("not enough Arguments received")} 


ukbb_tab_loc = args[1]
targetPheno= args[2]
outputLoc = args[3]
secondPheno = NULL
if(length(args) > 3) secondPheno = args[4]

# Debug Vars:
#ukbb_tab_loc="C:/softwares/Cluster/UKBB/ukb7725.tab" # "/lustre/scratch115/teams/anderson/ukbiobank/ibs/data/raw/ukb7725.tab"
#targetPheno= "f.34.0.0"  # "f.20127.0.0"
#secondPheno = "f.53.0.0" # 
#outputLoc="C:/softwares/Cluster/UKBB/ukbb_neuroExport.txt"
#bd <- read.table(ukbb_tab_loc, header=TRUE, sep="\t", nrows=100) # /lustre/scratch115/teams/anderson/ukbiobank/ibs/data/raw/ukb7725.tab


bd <- read.table(ukbb_tab_loc, header=TRUE, sep="\t") # /lustre/scratch115/teams/anderson/ukbiobank/ibs/data/raw/ukb7725.tab

# if there is an optional second pheno, we will try to produce a final export as their difference: SecondPheno - firstPheno
if ( is.null(secondPheno) == FALSE ) {
  secondPhenoData = bd[secondPheno]  
  if(length(args) > 4) { # if there are more arguments, that will mean to cast the second pheno as date, and export only the date
    secondPhenoData = as.numeric( format(as.Date(secondPhenoData[,1]),"%Y") )
    }
  exportPheno = secondPhenoData -  bd[targetPheno]
} else {
  exportPheno = bd[targetPheno]
}


# mean( as.numeric(exportPheno[,1]) ) # age: 57.29

exportTable = cbind(bd$f.eid, exportPheno)


write.table(exportTable, outputLoc, sep=" ", row.names = F, col.names = F, quote = FALSE) 
print(paste("Written phenotypes to:",outputLoc   ))
