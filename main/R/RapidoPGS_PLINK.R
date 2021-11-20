#############################
# RapidoPGS script: takse a PLINK GWAS result and produce a RapidoPGS PRS file
#############################

# imports
library("RapidoPGS")
options(error=traceback)
# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 5) {stop("not enough Arguments received")}


# arguments
sumstatsLoc=args[1] # location of the PLINK .assoc or .qassoc file
n_eff =  as.numeric(args[2]) # number of individuals in study
MAFloc=args[3] # PLINK minor allele freq location
outLoc=args[4] # output location
outName=args[5] # outputfile name
LDrefPanel=args[6] # an optional LDref panel

#############################
# LOCAL DEBUG VARS
# sumstatsLoc="/lustre/scratch114/projects/crohns/mk23/missense_eqtl_epistasis/data_analysis/eqtl_missense_merged_all_sanity.assoc"
# sumstatsLoc="/lustre/scratch114/projects/crohns/mk23/0Thesis/shaprs/simPhe2/5000/0.5_1.0/A50_B50/size/4_pheA.qassoc"
# 
# n_eff=17554
# MAFloc="/lustre/scratch114/projects/crohns/mk23/0Thesis/shaprs/raw/_MAF.frq"
# outLoc= "/lustre/scratch114/projects/crohns/mk23/0Thesis/shaprs/scratch/ldpred2/"
# outName= "rapidotest"
# 
#########################################################



# 1. load PLINK .assoc or .qassoc
sumstats <- read.table(sumstatsLoc, header = T)

# determine if binary or continuous phenotypes were supplied by checking if the assoc file has an 'OR' column
if("OR" %in% colnames(sumstats)) {
message("Binary phenotype inferred")
  trait = "cc"
  sumstats$BETA = log(sumstats$OR)
} else { 
  message("Continuous phenotype inferred")
  trait = "quant"
}


# 2. load nad add MAF to sumstats
MAFs <- read.table(MAFloc, header = T)
sumstats_MAF = merge(sumstats,MAFs)

# 3. set names that RapidoPGS expects
setnames(sumstats_MAF, old = 
           c("SNP"  , "CHR", "BP", "A2", "A1" , "BETA", "MAF"     , "SE", "P"), new=
           c("SNPID","CHR" , "BP","REF", "ALT", "BETA", "ALT_FREQ", "SE", "P"))
                                                                                                                                                                      



# 4, run RapidoPGS, single or multi, depending on if an LD ref was supplied
if( is.na(LDrefPanel) == F ) { # exists("LDrefPanel")
  print(paste0("rapidopgs_multi from: ", LDrefPanel))
 # sumstats_MAF_NOSNPID = sumstats_MAF
 # sumstats_MAF_NOSNPID$SNPID <- NULL
 # full_PGS <- rapidopgs_multi(sumstats_MAF_NOSNPID, N =n_eff , trait = trait, LDmatrices = LDrefPanel)
  sumstats_MAF_DF = sumstats_MAF
  setDT(sumstats_MAF_DF)
full_PGS <- rapidopgs_multi(sumstats_MAF_DF, N =n_eff , trait = trait, LDmatrices = LDrefPanel)
# recover the original rsids, as currently rapidopgs_multi returns SNP ids in the format of chr:pos
OrigSNPID <- sumstats_MAF$SNPID[match(  full_PGS$SNPID ,paste(sumstats_MAF_DF$CHR,sumstats_MAF_DF$BP, sep=":"))]

#OrigSNP <- sumstats_MAF[match(  full_PGS$SNPID ,paste(sumstats_MAF_DF$CHR,sumstats_MAF_DF$BP, sep=":")),]
full_PGS$SNPID = OrigSNPID



} else {
  print("rapidopgs_single")
full_PGS <- rapidopgs_single(sumstats_MAF, N =n_eff , trait = trait)
}

# 9. produce the final output in PLINK PRS format that would work for PLINK --score sum  (columns: rsid, effect_allele, beta)
PRSFILE = cbind.data.frame(full_PGS$SNPID, full_PGS$ALT, full_PGS$weight)
write.table(PRSFILE, paste0(outLoc,outName), sep = "\t", row.names = F, col.names = F, quote = FALSE)
print(paste("written ",nrow(PRSFILE)," RapidoPGS adjusted betas to", paste0(outLoc,outName)))




