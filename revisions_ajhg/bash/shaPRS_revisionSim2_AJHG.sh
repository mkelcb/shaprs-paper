####################################################
#
# UKBB Simulation studies: Revision, replaces shaPRS_ukbb_extra.sh
# this adds MTAG
# and additional scenarios 
#- 40:60 for primary:adjunct, at h=0.5, numCausals=1000
#- numCausals, 3000, 5000, at h=0.5 at 50:50
#- h =0.25,0.75, at 50:50 with numcausals 1000
#
# screens:  
# screen -r -D 16461.pts-0.node-11-1-3
#
#########################
# static vars
largePlinkMem=200000
bigMem=400000
largeMem=60000
plinkMem=6000  # how much RAM PLINK would use
largerMem=10000 
homeBase='/nfs/users/nfs_m/mk23/'
plink=$homeBase$'software/plink/plink'
plink2=$homeBase$'software/plink/plink2'
scriptsLoc=$homeBase$'scripts/'
#Rscript='/nfs/team152/mk23/software/R/R-3.6.1/bin/Rscript'
#Rscript='/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript'
Rscript='/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript'
shaPRSscriptLoc="/nfs/users/nfs_m/mk23/scripts/shaPRS.R"
baseLDpredLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/"
flashpca='/nfs/users/nfs_m/mk23/software/flashpca/flashpca_x86-64'
qctool2_new='/nfs/team152/mk23/software/qctool2/qctool_v2.0.6-Ubuntu16.04-x86_64/qctool'
qctoolmem=10000

NCORES_LDPred2=15
ldPred2Mem=125000
LDPred_FileStem='LDPred2.PRS'

i=1
j=1
c=1
z=0
numFolds=5
numBootStraps=20
numChroms=22

HLA_chrom=6
HLA_start=28477797
HLA_end=33448354

baseLoc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/'

 

h2s=( 0.25 0.5 0.75 )
pheno_splitA=( 50 20 40 ) # 1= 50/50, 2=20/80
pheno_splitB=( 50 80 60 ) # 1= 50/50, 2=20/80
h2=0.5
numCausals=( 1000 3000 5000 ) 
 
rGs=( 0.1 0.25 0.5 ) 
samplesizes=( '_half' '' '_double' ) # 4000= half of original, 8000=original, 16000= indis are duplicated (within their original phenos)
heterogeneity=( 'extra' 'regular' ) # extra= where the last 5 SNPs have large effect size explaining 5% variance, and 'regular' where this is not the case
numHeterogSNPs=5
percVarHetSNPsExplain=0.05

replicates=50

#p / shared_corr (shared fraction of SNPs / correlation of effect sizes): min/middle/max for each rG
#	rG=0.1  : 0.1/1    , 0.55/0.1818182    , 1/0.1
#	rG=0.25 : 0.25/1   , 0.625/0.4         , 1/0.25
#	rG=0.5  : 0.5/1    , 0.75/0.6666667    , 1/0.5

rG01_ps=( 0.1 0.55 1.0 ) 
rG01_shared_corrs=( 1.0 0.1818182 0.1 ) 

rG025_ps=( 0.25 0.625 1.0 ) 
rG025_shared_corrs=( 1.0 0.4 0.25 ) 

rG05_ps=( 0.5 0.75 1.0 ) 
rG05_shared_corrs=( 1.0 0.6666667 0.5 ) 

rG01_ps=( "0.1" "0.55" "1.0" ) 
rG01_shared_corrs=( "1.0" "0.1818182" "0.1" ) 

rG025_ps=( "0.25" "0.625" "1.0" ) 
rG025_shared_corrs=( "1.0" "0.4" "0.25" ) 

rG05_ps=( "0.5" "0.75" "1.0" ) 
rG05_shared_corrs=( "1.0" "0.6666667" "0.5" ) 


# LD ref panel

ldak=$homeBase$'software/ldak/ldak5.linux.fast'   # $homeBase$'software/ldak/ldak.4.9' 

ldak51=$homeBase$'software/ldak/ldak5.1.linux.fast'



source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh
conda activate mk23_pytorch


#wget http://dougspeed.com/wp-content/uploads/ldak5.1.linux_.zip




######################################

# fix h2 at 0.5 
baseLoc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/'
shaprsLoc=$baseLoc$'shaprs_ukbb/'
shaprsScratch=$shaprsLoc$'scratch/'
shaprsRaw=$shaprsLoc$'raw/'
shaprsResults=$shaprsLoc$'results/'
shaprsSimPhe=$shaprsLoc$'simPhe/'
testSetLdpredLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/raw/pheA_testSet.bed"
pheA_testSet=$shaprsRaw$'pheA_testSet'
pheB_testSet=$shaprsRaw$'pheB_testSet'
NCORES=1
PRS_FileStem="RAPIDOPGS.PRS"
plinkMem=2000  # used for LDpred2, as a 1Mx1M sparse matrix uses 32GB of RAM
logdir=$shaprsSimPhe'logs/'

# 5% of h2 explained by 5 SNPs
numLargeEffectNonSharedSNPs=5
numLargeEffectNonSharedSNPs_h2=0.05
ukbbhapmap3loc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/pgslite2/height/data/GWAS_ALL'

ukbb32K=$shaprsRaw$'ukbb32K'


mkdir -p $logdir
mkdir -p $shaprsRaw
mkdir -p $shaprsScratch
mkdir -p $shaprsResults
mkdir -p $shaprsResults

##################################
# 0) Get a subset of the UKBB:
######################

# randomise people order and pick the same number of people as we had in IBD
shuf -n 31598 $ukbbhapmap3loc$'.fam' >  $ukbb32K

# subset 
arguments=' --memory '$largeMem$' --bfile '$ukbbhapmap3loc$' -make-bed --out '$ukbb32K$' --keep '$ukbb32K$' --allow-extra-chr --allow-no-sex'
$plink $arguments

##################################
# I) Phenotype simulation:
######################
# 1) find out MAFs
arguments=' --memory '$plinkMem$' --bfile '$ukbb32K$' --freq --out '$shaprsScratch$'gwas3freq --allow-extra-chr --allow-no-sex'
$plink $arguments

awk '{ if ($5 > 0.01) print $2 }' $shaprsScratch$'gwas3freq.frq' > $shaprsRaw$'commonMAF_all'
head $shaprsRaw$'commonMAF_all'
wc -l $shaprsRaw$'commonMAF_all' # 1042687

# intersect it with HapMap3: we already have that
cp $shaprsRaw$'commonMAF_all' $shaprsRaw$'commonMAF'


# produce MAF ref data
arguments=' --memory '$plinkMem$' --bfile '$ukbb32K$' --freq --out '$shaprsRaw$'_MAF --allow-extra-chr --allow-no-sex'
$plink $arguments


####################################
# Organise cohorts: divide into pheA/B

# create PheA/PheB indis 
shuf $ukbb32K$'.fam' > $shaprsScratch$'allIndis'

# split Test set off from the training set (20% test, 80% training), but since we have a genuine 32K people, we only really want 20% of the original 16K
numIndis=$(wc -l < "$shaprsScratch"allIndis)
num_20p=$((numIndis / 10)) # 10, as we want 10% of 32K which is the same as 20% of 16K
awk -v num_20p=$num_20p -v top_file=$shaprsRaw$'/phenoBoth_test' -v bottom_file=$shaprsRaw$'/phenoBoth_train' 'NR <= num_20p { print > top_file; next } {print > bottom_file }' $shaprsScratch$'/allIndis'

# split the test set in two for phenotpes A and B
numIndis=$(wc -l < "$shaprsRaw"phenoBoth_test)
num_half=$((numIndis / 2))
awk -v num_half=$num_half -v top_file=$shaprsRaw$'/phenoA_indis_test' -v bottom_file=$shaprsRaw$'/phenoB_indis_test' 'NR <= num_half { print > top_file; next } {print > bottom_file }' $shaprsRaw$'/phenoBoth_test'
wc -l $shaprsRaw$'/phenoA_indis_test' # 1579


# 1) Training split 20/80 for the double
numIndis=$(wc -l < "$shaprsRaw"phenoBoth_train)
num_20p=$((numIndis / 5))
awk -v num_20p=$num_20p -v top_file=$shaprsRaw$'/phenoA_train20_double' -v bottom_file=$shaprsRaw$'/phenoB_train80_double' 'NR <= num_20p { print > top_file; next } {print > bottom_file }' $shaprsRaw$'/phenoBoth_train'
wc -l $shaprsRaw$'/phenoA_train20_double' # 5687
wc -l $shaprsRaw$'/phenoB_train80_double' # 22752


# 2) Training split 50/50 for the double
numIndis=$(wc -l < "$shaprsRaw"phenoBoth_train)
num_20p=$((numIndis / 2))
awk -v num_20p=$num_20p -v top_file=$shaprsRaw$'/phenoA_train50_double' -v bottom_file=$shaprsRaw$'/phenoB_train50_double' 'NR <= num_20p { print > top_file; next } {print > bottom_file }' $shaprsRaw$'/phenoBoth_train'
wc -l $shaprsRaw$'/phenoA_train50_double'  # 14219
wc -l $shaprsRaw$'/phenoB_train50_double'  # 14220


# 3) Training split 40/60 for the double
numIndis=$(wc -l < "$shaprsRaw"phenoBoth_train)
num_20p=$((numIndis / 10))
num_20p=$((num_20p * 4))
awk -v num_20p=$num_20p -v top_file=$shaprsRaw$'/phenoA_train40_double' -v bottom_file=$shaprsRaw$'/phenoB_train60_double' 'NR <= num_20p { print > top_file; next } {print > bottom_file }' $shaprsRaw$'/phenoBoth_train'
wc -l $shaprsRaw$'/phenoA_train40_double'  # 11372
wc -l $shaprsRaw$'/phenoB_train60_double'  # 17067



wc -l $shaprsRaw$'/phenoA_train20' # 5687
wc -l $shaprsRaw$'/phenoB_train80' # 22752
wc -l $shaprsRaw$'/phenoA_train50'  # 14219
wc -l $shaprsRaw$'/phenoB_train50'  # 14220


# 3) Create Normal and Half sized cohorts:  half the _double set to get _full, then half that to get _half
numIndis=$(wc -l < "$shaprsRaw"phenoA_train20_double)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoA_train20_double' >  $shaprsRaw$'/phenoA_train20'
wc -l $shaprsRaw$'/phenoA_train20' # 2843

numIndis=$(wc -l < "$shaprsRaw"phenoB_train80_double)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoB_train80_double' >  $shaprsRaw$'/phenoB_train80'
wc -l $shaprsRaw$'/phenoB_train80' # 11376

numIndis=$(wc -l < "$shaprsRaw"phenoA_train50_double)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoA_train50_double' >  $shaprsRaw$'/phenoA_train50'
wc -l $shaprsRaw$'/phenoA_train50' # 7109

numIndis=$(wc -l < "$shaprsRaw"phenoB_train50_double)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoB_train50_double' >  $shaprsRaw$'/phenoB_train50'
wc -l $shaprsRaw$'/phenoB_train50' # 7110


numIndis=$(wc -l < "$shaprsRaw"phenoA_train40_double)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoA_train40_double' >  $shaprsRaw$'/phenoA_train40'
wc -l $shaprsRaw$'/phenoA_train40' # 5686

numIndis=$(wc -l < "$shaprsRaw"phenoB_train60_double)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoB_train60_double' >  $shaprsRaw$'/phenoB_train60'
wc -l $shaprsRaw$'/phenoB_train60' # 8533





numIndis=$(wc -l < "$shaprsRaw"phenoA_train20)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoA_train20' >  $shaprsRaw$'/phenoA_train20_half'
wc -l $shaprsRaw$'/phenoA_train20_half'  # 1421

numIndis=$(wc -l < "$shaprsRaw"phenoB_train80)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoB_train80' >  $shaprsRaw$'/phenoB_train80_half'
wc -l $shaprsRaw$'/phenoB_train80_half' # 5688

numIndis=$(wc -l < "$shaprsRaw"phenoA_train50)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoA_train50' >  $shaprsRaw$'/phenoA_train50_half'
wc -l $shaprsRaw$'/phenoA_train50_half' # 3554

numIndis=$(wc -l < "$shaprsRaw"phenoB_train50)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoB_train50' >  $shaprsRaw$'/phenoB_train50_half'
wc -l $shaprsRaw$'/phenoB_train50_half' # 3555


numIndis=$(wc -l < "$shaprsRaw"phenoA_train40)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoA_train40' >  $shaprsRaw$'/phenoA_train40_half'
wc -l $shaprsRaw$'/phenoA_train40_half' # 2843

numIndis=$(wc -l < "$shaprsRaw"phenoB_train60)
num_half=$((numIndis / 2))
awk -v num_half=$num_half 'NR <= num_half { print $0 }'  $shaprsRaw$'/phenoB_train60' >  $shaprsRaw$'/phenoB_train60_half'
wc -l $shaprsRaw$'/phenoB_train60_half' # 4266



# create 'Combined' pheno lists
cat $shaprsRaw$'/phenoA_train20_double' $shaprsRaw$'/phenoB_train80_double' > $shaprsRaw$'/phenoB_train2080_double_all'
cat $shaprsRaw$'/phenoA_train50_double' $shaprsRaw$'/phenoB_train50_double' > $shaprsRaw$'/phenoB_train5050_double_all'
cat $shaprsRaw$'/phenoA_train40_double' $shaprsRaw$'/phenoB_train60_double' > $shaprsRaw$'/phenoB_train4060_double_all'

cat $shaprsRaw$'/phenoA_train20_half' $shaprsRaw$'/phenoB_train80_half' > $shaprsRaw$'/phenoB_train2080_half_all'
cat $shaprsRaw$'/phenoA_train50_half' $shaprsRaw$'/phenoB_train50_half' > $shaprsRaw$'/phenoB_train5050_half_all'
cat $shaprsRaw$'/phenoA_train40_half' $shaprsRaw$'/phenoB_train60_half' > $shaprsRaw$'/phenoB_train4060_half_all'

cat $shaprsRaw$'/phenoA_train20' $shaprsRaw$'/phenoB_train80' > $shaprsRaw$'/phenoB_train2080_all'
cat $shaprsRaw$'/phenoA_train50' $shaprsRaw$'/phenoB_train50' > $shaprsRaw$'/phenoB_train5050_all'
cat $shaprsRaw$'/phenoA_train40' $shaprsRaw$'/phenoB_train60' > $shaprsRaw$'/phenoB_train4060_all'


# Create Test Set genotype file for phe A (LDPred cannot subset indis)
arguments=' --memory '$plinkMem$' --threads 5 --bfile '$ukbb32K$' --allow-no-sex --keep '$shaprsRaw$'phenoA_indis_test --make-bed --out '$pheA_testSet
$plink $arguments

#Create Test Set genotype file for phe B
arguments=' --memory '$plinkMem$' --threads 5 --bfile '$ukbb32K$' --allow-no-sex --keep '$shaprsRaw$'phenoB_indis_test --make-bed --out '$pheB_testSet
$plink $arguments


# create hapmap3 version of only the unique training indis
arguments=' --memory '$plinkMem$' --threads 10 --bfile '$ukbb32K$' --keep '$shaprsRaw$'phenoBoth_train --indiv-sort f '$shaprsRaw$'phenoBoth_train --make-bed --out '$shaprsRaw$'phenoBoth_train_hapmap3 --allow-no-sex'
$plink $arguments

# create tagging file for h2 sumher analyses
arguments=' --thin '$shaprsRaw$'thin --bfile '$shaprsRaw$'phenoBoth_train_hapmap3 --window-prune .98 --window-kb 100 '
$ldak $arguments
awk '{print $1, 1}'  $shaprsRaw$'thin.in' > $shaprsRaw$'weights.thin'
arguments=' --calc-tagging '$shaprsRaw'LDAK --bfile '$shaprsRaw$'phenoBoth_train_hapmap3 --weights '$shaprsRaw$'weights.thin --power -.25 --window-kb 1000'
$ldak $arguments
####################################


##################################
# II) Pick causal variants: reuse same causals for all rG/splits/samplesizes/heterog
##########################

c=1

arraylength=${#numCausals[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
causals_all=${numCausals[$j-1]}

mkdir -p $shaprsRaw$causals_all$'/'
# go through each pheno replicate
for ((currentPhe=1; currentPhe<=50; currentPhe++)); do 
echo $currentPhe

# pick the causals
shuf -n $causals_all $shaprsRaw$'commonMAF' >  $shaprsRaw$causals_all$'/'$currentPhe

# LDAK expects the phenotypes A/B in a format of 
# PHE1 SNP1 SNP2 SNP3...
# PHE2 SNP1 SNP2 SNP3...
paste -s -d ' ' $shaprsRaw$causals_all$'/'$currentPhe > $shaprsScratch$'/'$currentPhe$'_'$causals_all
cat  $shaprsScratch$'/'$currentPhe$'_'$causals_all > $shaprsScratch$'/'$currentPhe$'_'$causals_all$'_PheAB_casuals'
cat  $shaprsScratch$'/'$currentPhe$'_'$causals_all >> $shaprsScratch$'/'$currentPhe$'_'$causals_all$'_PheAB_casuals'

done # loop num replicates

done # loop num causals

#############################
# III) Simulate phenotypes (this is independent of splits and sample sizes)


# Reusable pheno gen function for the non-heterogeneous case
function pheno_gen { 
causals_all=$1
p_current=$2
shared_corr=$3
rGText=$4
h2=$5 # '0.5'
h2Text=$(awk -v h2=$h2 'BEGIN{ print h2*100}') 
current_het='regular' # hardcode heterogeneity modifier

plinkMem=6000  # how much RAM PLINK would use


# hardcoded vars to have fewer function params
plink='/nfs/users/nfs_m/mk23/software/plink/plink'
shaprsRaw='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/raw/'
shaprsSimPhe='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/simPhe/'
shaprsScratch='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/scratch/'
ukbb32K='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/raw/ukbb32K'

workingDir=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het$'/'$h2Text$'/'
mkdir -p $workingDir


# go through each pheno replicate
for ((currentPhe=1; currentPhe<=20; currentPhe++)); do #replicates

# I) simulate SNP effect sizes via R script, given p and shared_corr
# (this is independent of SplitA/B), and also of sample size, but not independent of replicates  # 0.04
arguments='/nfs/users/nfs_m/mk23/scripts/rG_effectSim_SNPs.R '$causals_all$' '$p_current$' '$shared_corr$' '$shaprsRaw$causals_all$'/'$currentPhe$' '$currentPhe$' '$workingDir$currentPhe
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# 4) use LDAK to generate actual phenos, based on the genetic architecture (but independent of split/sample size and how we split them, IE all indis will have all phenos, we will just use smaller subsets of them later)
arguments=' --make-phenos '$workingDir$currentPhe$' --bfile '$ukbb32K$' --ignore-weights YES --power -1 --her '$h2$' --num-phenos 2 --effects '$workingDir$currentPhe$'  --causals '$shaprsScratch$'/'$currentPhe$'_'$causals_all$'_PheAB_casuals'
/nfs/users/nfs_m/mk23/software/ldak/ldak5.linux.fast $arguments
rm -rf $workingDir$currentPhe$'.progress'

#######################
# (optional) Sanity check: verify rG via LDAK as a sanity check 
# GWAS with LDAK
# awk '{print $1" "$2" "$3}' $workingDir$currentPhe$'.pheno' > $workingDir$currentPhe$'_A'
# awk '{print $1" "$2" "$4}' $workingDir$currentPhe$'.pheno' > $workingDir$currentPhe$'_B'
# arguments=' --linear '$workingDir$currentPhe$'_A_GWAS --bfile '$ukbb32K$' --pheno '$workingDir$currentPhe$'_A'
# $ldak $arguments
# arguments=' --linear '$workingDir$currentPhe$'_B_GWAS --bfile '$ukbb32K$' --pheno '$workingDir$currentPhe$'_B'
# $ldak $arguments

# arguments=' --sum-cors '$workingDir$currentPhe$'_rG  --check-sums NO --summary '$workingDir$currentPhe$'_A_GWAS.summaries --summary2 '$workingDir$currentPhe$'_B_GWAS.summaries --tagfile '$shaprsRaw'LDAK.tagging --allow-ambiguous YES'
# $ldak $arguments # Overall correlation (SD): 0.489306 (0.047509)
#######################

done # end of loop num replicates

}
export -f pheno_gen # this makes local functions executable when bsubbed



# Reusable pheno gen function for the extra heterogeneous case
function pheno_gen_extra { 
causals_all=$1
p_current=$2
shared_corr=$3
rGText=$4
h2=$5 # '0.5'
h2Text=$(awk -v h2=$h2 'BEGIN{ print h2*100}') 
current_het='extra' # hardcode heterogeneity modifier

numHeterogSNPs=5
percVarHetSNPsExplain=0.05
plinkMem=6000  # how much RAM PLINK would use


# hardcoded vars to have fewer function params

plink='/nfs/users/nfs_m/mk23/software/plink/plink'
shaprsRaw='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/raw/'
shaprsSimPhe='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/simPhe/'
shaprsScratch='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/scratch/'
ukbb32K='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/raw/ukbb32K'


workingDir=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het$'/'$h2Text$'/'
mkdir -p $workingDir

# go through each pheno replicate
for ((currentPhe=1; currentPhe<=20; currentPhe++)); do #replicates

# I) simulate SNP effect sizes via R script, given p and shared_corr
# (this is independent of SplitA/B), and also of sample size, but not independent of replicates  # 0.04
arguments='/nfs/users/nfs_m/mk23/scripts/rG_effectSim_SNPs_heterog.R '$causals_all$' '$p_current$' '$shared_corr$' '$shaprsRaw$causals_all$'/'$currentPhe$' '$currentPhe$' '$numHeterogSNPs$' '$workingDir$currentPhe
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments



# 4) use LDAK to generate bivariate phenos, but with heritability 1.0 for the '_regular' case (IE the 995 SNPs)
arguments=' --make-phenos '$workingDir$currentPhe$'_regular --bfile '$ukbb32K$' --ignore-weights YES --power -1 --her 1.0 --num-phenos 2 --effects '$workingDir$currentPhe$'_regular --causals '$workingDir$currentPhe$'_SNPs_regular_causals'
/nfs/users/nfs_m/mk23/software/ldak/ldak5.linux.fast $arguments
rm -rf $workingDir$currentPhe$'_regular.progress'

# and for the 5 SNPs
# 4) use LDAK to generate bivariate phenos, but with heritability 1.0 for the '_regular' case (IE the 995 SNPs)
arguments=' --make-phenos '$workingDir$currentPhe$'_extra --bfile '$ukbb32K$' --ignore-weights YES --power -1 --her 1.0 --num-phenos 2 --effects '$workingDir$currentPhe$'_extra --causals '$workingDir$currentPhe$'_SNPs_extra_causals'
/nfs/users/nfs_m/mk23/software/ldak/ldak5.linux.fast $arguments
rm -rf $workingDir$currentPhe$'_extra.progress'


# 5) Scale the 2 h2=1.0 phenotypes, so that the total comes to the desired overal h2=0.5
h2Min5=$(awk -v h2=$h2 'BEGIN{ print h2-0.05}') 
arguments='/nfs/users/nfs_m/mk23/scripts/phenoScaler.R '$h2Min5$' 0.05 '$workingDir$currentPhe$'_regular.pheno '$workingDir$currentPhe$'_extra.pheno '$currentPhe$' '$workingDir$currentPhe
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments



####################### aa
# (optional) Sanity check: verify rG via LDAK as a sanity check 
# GWAS with LDAK
# awk '{print $1" "$2" "$3}' $workingDir$currentPhe$'.pheno' > $workingDir$currentPhe$'_A'
# awk '{print $1" "$2" "$4}' $workingDir$currentPhe$'.pheno' > $workingDir$currentPhe$'_B'


# arguments=' --linear '$workingDir$currentPhe$'_A_GWAS --bfile '$ukbb32K$' --pheno '$workingDir$currentPhe$'_A'
# $ldak $arguments
# arguments=' --linear '$workingDir$currentPhe$'_B_GWAS --bfile '$ukbb32K$' --pheno '$workingDir$currentPhe$'_B'
# $ldak $arguments

# arguments=' --sum-cors '$workingDir$currentPhe$'_rG  --check-sums NO --summary '$workingDir$currentPhe$'_A_GWAS.summaries --summary2 '$workingDir$currentPhe$'_B_GWAS.summaries --tagfile '$shaprsRaw'LDAK.tagging --allow-ambiguous YES'
# $ldak $arguments # Overall correlation (SD): 0.428021 (0.111451)
# # Final   0.4646  0.4646  1.0105  -0.0000 0.000100
# # Final   0.4473  0.4473  1.0210  -0.0000 0.000100

# # the overal rG turned to 0, as when we add 5 SNPs that explain 5% of h2, this will mess up the overall rG
# # 
# # try to reestimate the h2/rG with only the 995 'regular' SNPs

# # remove the 5 extra SNPs from the GWAS results:
# awk 'FNR == NR { file1[ $1 ] = $1; next; } FNR <= NR {  if (FNR == 1) {print $0} else {if ($1 in file1 ==0) {print $0}}  }
# ' $workingDir$currentPhe$'_SNPs_extra' $workingDir$currentPhe$'_A_GWAS.summaries' > $workingDir$currentPhe$'_A_GWAS_regular.summaries'

# awk 'FNR == NR { file1[ $1 ] = $1; next; } FNR <= NR {  if (FNR == 1) {print $0} else {if ($1 in file1 ==0) {print $0}}  }
# ' $workingDir$currentPhe$'_SNPs_extra' $workingDir$currentPhe$'_B_GWAS.summaries' > $workingDir$currentPhe$'_B_GWAS_regular.summaries'


# arguments=' --sum-cors '$workingDir$currentPhe$'_rG_regular  --check-sums NO --summary '$workingDir$currentPhe$'_A_GWAS_regular.summaries --summary2 '$workingDir$currentPhe$'_B_GWAS_regular.summaries --tagfile '$shaprsRaw'LDAK.tagging --allow-ambiguous YES'
# $ldak $arguments # Overall correlation (SD): -0.229354 (2.093698)


#######################

done # end of loop num replicates

}
export -f pheno_gen_extra # this makes local functions executable when bsubbed


h2s=( 0.25 0.75 )
numCausals=( 1000 3000 5000 ) 


arraylength=${#numCausals[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
causals_all=${numCausals[$j-1]}
mkdir -p $shaprsSimPhe$causals_all$'/'

# go through each genetic correlation
arraylength_rG=${#rGs[@]}
for (( k=1; k<${arraylength_rG}+1; k++ )); do
rG=${rGs[$k-1]}
rGText=$(awk -v rG=$rG 'BEGIN{ print rG*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point
mkdir -p $shaprsSimPhe$causals_all$'/'$rGText$'/'

# based on RG, select the heterogeneity architecture
case $rG in
0.1)
  Ps=( "${rG01_ps[@]}" )
  shared_corrs=( "${rG01_shared_corrs[@]}" )
  echo "rG01"
  ;;
0.25)
  Ps=( "${rG025_ps[@]}" )
  shared_corrs=( "${rG025_shared_corrs[@]}" )
  echo "rG025"
  ;;
*)
  Ps=( "${rG05_ps[@]}" )
  shared_corrs=( "${rG05_shared_corrs[@]}" )
  echo "rG05"
  ;;
esac


# go through each heterogeneity architecture
arraylength_Ps=${#Ps[@]}
for (( b=1; b<${arraylength_Ps}+1; b++ )); do
p_current=${Ps[$b-1]}
shared_corr=${shared_corrs[$b-1]}


# go through the 2 extra/regular
arraylength_heterogeneity=${#heterogeneity[@]}
for (( e=1; e<${arraylength_heterogeneity}+1; e++ )); do
current_het=${heterogeneity[$e-1]}


arraylength_h2=${#h2s[@]}
for (( r=1; r<${arraylength_h2}+1; r++ )); do
h2=${h2s[$r-1]}


h2Text=$(awk -v h2=$h2 'BEGIN{ print h2*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point


echo "numCausals: "$causals_all$" / rG: "$rG$" | p: "$p_current$" / shared_corr: "$shared_corr$" | current_het: "$current_het$" | h2: "$h2

mkdir -p $shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het$'/'$h2Text$'/'
if [[ "$current_het" == 'extra' ]]; then
echo 'extra'
bsub -q normal -m "modern_hardware" -G team152 -n1 -J PHESIM${sample_size}_p${p_current}_sc${shared_corr}_${current_het} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $logdir$sample_size$'_p'$p_current$'_sc'$shared_corr$'_het'$current_het$'_h2'$h2$'.out' -e $logdir$sample_size$'_p'$p_current$'_sc'$shared_corr$'_het'$current_het$'_h2'$h2$'.err'  "pheno_gen_extra $causals_all $p_current $shared_corr $rGText $h2"
else
echo 'regular'
bsub -q normal -m "modern_hardware" -G team152 -n1 -J PHESIM${sample_size}_p${p_current}_sc${shared_corr}_${current_het} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $logdir$sample_size$'_p'$p_current$'_sc'$shared_corr$'_het'$current_het$'_h2'$h2$'.out' -e $logdir$sample_size$'_p'$p_current$'_sc'$shared_corr$'_het'$current_het$'_h2'$h2$'.err'  "pheno_gen $causals_all $p_current $shared_corr $rGText $h2"
fi

done # end of h2s

done # end of extra/regular

done # end of heterogeneity architecture loop

done # end of rG loop

done # num causals



# need to move the directory of the h2=0.5 for the 1000 run 


#################################
#        REUSABLE FUNCTIONS     #
#################################
$shaprsRaw$'/phenoA_train20_double'
$shaprsRaw$'/phenoB_train80_double'
$shaprsRaw$'/phenoA_train50_double'
$shaprsRaw$'/phenoB_train50_double'
$shaprsRaw$'/phenoA_train40_double'
$shaprsRaw$'/phenoB_train60_double'

$shaprsRaw$'/phenoA_train20_half'
$shaprsRaw$'/phenoB_train80_half'
$shaprsRaw$'/phenoA_train50_half'
$shaprsRaw$'/phenoB_train50_half'
$shaprsRaw$'/phenoA_train40_half'
$shaprsRaw$'/phenoB_train60_half'

$shaprsRaw$'/phenoA_train20'
$shaprsRaw$'/phenoB_train80'
$shaprsRaw$'/phenoA_train50'
$shaprsRaw$'/phenoB_train50'
$shaprsRaw$'/phenoA_train40'
$shaprsRaw$'/phenoB_train60'

# 'Combined' pheno lists
$shaprsRaw$'/phenoB_train2080_double_all'
$shaprsRaw$'/phenoB_train5050_double_all'
$shaprsRaw$'/phenoB_train4060_double_all'
$shaprsRaw$'/phenoA_train2080_half_all'
$shaprsRaw$'/phenoA_train5050_half_all'
$shaprsRaw$'/phenoA_train4060_half_all'
$shaprsRaw$'/phenoB_train2080_all'
$shaprsRaw$'/phenoB_train5050_all'
$shaprsRaw$'/phenoB_train4060_all'
$ukbb32K


current_het='extra'
pheA_=$shaprsRaw$'/phenoA_train50'
pheB_=$shaprsRaw$'/phenoB_train50'
pheCombined_=$shaprsRaw$'/phenoB_train2080_all'
outputDir=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het$'/A'$splitA$'_B'$splitB$'/size'$sample_size$'/'$currentPhe
phenoAllLoc=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het$'/'$h2Text$'/'$currentPhe$'.pheno'








# Reusable GWAS function
function perform_GWAS { 
shaprsRaw=$1
pheA_=$2
pheB_=$3
phenoAllLoc=$4
outputDir=$5

plinkMem=6000  # how much RAM PLINK would use
plink='/nfs/users/nfs_m/mk23/software/plink/plink'
ukbb32K='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/raw/ukbb32K'

# export out PheA and PheB for the relevant indis
awk 'FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {  if ($2 in file1) {print $1" "$2" "$3}  }
' $pheA_ $phenoAllLoc > $outputDir$'_pheA.phe'

awk 'FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {  if ($2 in file1) {print $1" "$2" "$4}  }
' $pheB_ $phenoAllLoc > $outputDir$'_pheB.phe'

# create composite phenos: from the FID IID PheA PheB, we take only PheA for CohortA indis 
cat $outputDir$'_pheA.phe' $outputDir$'_pheB.phe' > $outputDir$'_pheCombined.phe'
  
# perform GWAS phe A+B
arguments=' --memory '$plinkMem$' --threads 5 --bfile '$ukbb32K$' --keep-allele-order --pheno '$outputDir$'_pheCombined.phe --assoc --allow-no-sex --keep '$outputDir$'_pheCombined.phe --out '$outputDir$'_pheCombined'
$plink $arguments

# perform GWAS phe A
arguments=' --memory '$plinkMem$' --threads 5 --bfile '$ukbb32K$' --keep-allele-order --pheno '$outputDir$'_pheA.phe --assoc --allow-no-sex --keep '$outputDir$'_pheA.phe --out '$outputDir$'_pheA'
$plink $arguments

# perform GWAS phe B
arguments=' --memory '$plinkMem$' --threads 5 --bfile '$ukbb32K$' --keep-allele-order --pheno '$outputDir$'_pheB.phe --assoc --allow-no-sex --keep '$outputDir$'_pheB.phe --out '$outputDir$'_pheB'
$plink $arguments

rm -rf $outputDir$'_pheA.log'
rm -rf $outputDir$'_pheB.log'
rm -rf $outputDir$'_pheCombined.log'


# convert the above association into a summary statistics format suitable for the shaPRS blending script
# PLINK qassoc signature:
# CHR          SNP         BP    NMISS       BETA         SE         R2        T            P
avgNumIndis=$(wc -l < "$pheA_")
awk -v avgNumIndis="$avgNumIndis" ' 
FNR == NR {
	A1[ $2 ] = $5;
	A2[ $2 ] = $6;
	next;}
FNR <= NR {
# PLINK .bim file is A1/A2 (usually ALT/REF)
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { print $1"\t"$3"\t"$2"\t"A1[ $2 ]"\t"A2[ $2 ]"\tX\t"$5"\t"$6"\t"$9"\t"avgNumIndis }
}
' OFS='\t' $ukbb32K$'.bim' $outputDir$'_pheA.qassoc' > $outputDir$'_pheA_sumstats'

avgNumIndis=$(wc -l < "$pheB_")
awk -v avgNumIndis="$avgNumIndis" ' 
FNR == NR {
	A1[ $2 ] = $5;
	A2[ $2 ] = $6;
	next;}
FNR <= NR {
# PLINK .bim file is A1/A2 (usually ALT/REF)
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { print $1"\t"$3"\t"$2"\t"A1[ $2 ]"\t"A2[ $2 ]"\tX\t"$5"\t"$6"\t"$9"\t"avgNumIndis }
}
' OFS='\t' $ukbb32K$'.bim' $outputDir$'_pheB.qassoc' > $outputDir$'_pheB_sumstats'

avgNumIndis=$(wc -l < "$outputDir"_pheCombined.phe)
awk -v avgNumIndis="$avgNumIndis" ' 
FNR == NR {
	A1[ $2 ] = $5;
	A2[ $2 ] = $6;
	next;}
FNR <= NR {
# PLINK .bim file is A1/A2 (usually ALT/REF)
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { print $1"\t"$3"\t"$2"\t"A1[ $2 ]"\t"A2[ $2 ]"\tX\t"$5"\t"$6"\t"$9"\t"avgNumIndis }
}
' OFS='\t' $ukbb32K$'.bim' $outputDir$'_pheCombined.qassoc' > $outputDir$'_pheCombined_sumstats'

}

export -f perform_GWAS # this makes local functions executable when bsubbed


###########################
# Blend runs:
################

# shaprsRaw=$shaprsRaw
# cohortFileStem=$cohortFileStem
# avgNumIndis=$avgNumIndis
# outputDir=$outputDir
# outPRSName=$outPRSName


    


# creates sumstats from two PLINK GWAS performing fixed effect meta analysis, then coordinates it then adjusts them to produce a PRS
function rapidopgs_Blend_PRS_meta { 
shaprsRaw=$1
cohortFileStem=$2
avgNumIndis=$3
outputDir=$4
outPRSName=$5

# Create input file for adjusting script:
awk 'FNR == NR { file1[ $3 ] = $7"\t"$8"\t"$4"\t"$5; next; } 
FNR <= NR { { 
if (FNR == 1) {print "SNP\tCHR\tBP\tBeta_A\tSE_A\tA1.x\tA2.x\tBeta_B\tSE_B\tA1.y\tA2.y" } else {
if ( $3 in file1) { print $3"\t"$1"\t"$2"\t"file1[$3]"\t"$7"\t"$8"\t"$4"\t"$5} }   }
}
' $outputDir$'_pheA_sumstats' $outputDir$'_pheB_sumstats' > $cohortFileStem$'_SE_meta'


# Export out the lFDR values
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_adjust_wrapper.R '$cohortFileStem$'_SE_meta '$cohortFileStem$'_lFDR_meta /nfs/users/nfs_m/mk23/scripts/shaPRS.R'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# R script to blend between the two association stats
arguments='/nfs/users/nfs_m/mk23/scripts/sumstatsBlender_shaPRS_meta_wrapper.R '$outputDir$'_pheA_sumstats '$outputDir$'_pheB_sumstats '$cohortFileStem$'_lFDR_meta_SNP_lFDR 0 '$cohortFileStem$'_sumstats_meta /nfs/users/nfs_m/mk23/scripts/shaPRS.R'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


## recreate the PLINK .qassoc formatted summary stats
awk 'FNR == NR { 
b[ $3 ] = $7; 
se[ $3 ] = $8; 
p[ $3 ] = $9; 
next; } 
FNR <= NR { { 
if (FNR == 1) {print "CHR\tSNP\tBP\tNMISS\tBETA\tSE\tR2\tT\tP" } else {
if ( $2 in b) { print $1"\t"$2"\t"$3"\t"$4"\t"b[ $2 ]"\t"se[ $2 ]"\t"$7"\t"$8"\t"p[ $2 ] }  }   }
}
'  $cohortFileStem$'_sumstats_meta' $outputDir$'_pheA.qassoc' > $cohortFileStem$'_meta.qassoc'

# produce the final PRS
arguments='/nfs/users/nfs_m/mk23/scripts/RapidoPGS_PLINK.R '$cohortFileStem$'_meta.qassoc '$avgNumIndis$' '$shaprsRaw$'_MAF.frq '$outputDir$' '$outPRSName
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments

}
export -f rapidopgs_Blend_PRS_meta # this makes local functions executable when bsubbed





# Reusable SMTPred function
function perform_SMPTpred { 
shaprsRaw=$1
avgNumIndis=$2
outputDir=$3
currentPhe=$4
rG=$5

plink='/nfs/users/nfs_m/mk23/software/plink/plink'
smtpred="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/smtpred/smtpred/"
ldscDir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/smtpred/ldsc/"
w_hm3_snplist="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/ldsc/w_hm3.snplist"
ldscRefpanelLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/ldsc/eur_w_ld_chr/"
pheA_testSet="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/raw/pheA_testSet"





# obtain h2 and rG estimates for Pheno A and B via LDSC
mkdir -p $outputDir
/usr/bin/python2  $smtpred$'ldsc_wrapper.py' \
    --out $outputDir \
    --files $outputDir$'_pheA_sumstats' \
            $outputDir$'_pheB_sumstats' \
    --ldscpath $ldscDir \
    --snplist $w_hm3_snplist \
    --ref_ld $ldscRefpanelLoc \
    --w_ld $ldscRefpanelLoc

# sometimes the above fails, for small n
if [ -s $outputDir$'/ldsc_rgs.txt' ]; then
b=1
else
echo "The file is empty."
#18_pheA_sumstats	18_pheB_sumstats	0.4348	2.4903e-01	0.0808
echo -e $currentPhe$'_pheA_sumstats\t'$currentPhe$'_pheB_sumstats\t'$rG$'\t0.1\t0.005' > $outputDir$'/ldsc_rgs.txt'
#head $outputDir$'/ldsc_rgs.txt'
fi




# f.write('\t'.join([n1, n2, row['rg'], row['se'], row['p']]) + '\n')

# produce reweighted PRS files with rapido PRS
arguments='/nfs/users/nfs_m/mk23/scripts/RapidoPGS_PLINK.R '$outputDir$'_pheA.qassoc '$avgNumIndis$' '$shaprsRaw$'_MAF.frq '$outputDir$'/ '$currentPhe$'_pheA_rapido.PRS'
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments
arguments='/nfs/users/nfs_m/mk23/scripts/RapidoPGS_PLINK.R '$outputDir$'_pheB.qassoc '$avgNumIndis$' '$shaprsRaw$'_MAF.frq '$outputDir$'/ '$currentPhe$'_pheB_rapido.PRS'
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments

# produce PLINK PRS score files (both use the same Test set)
arguments=' --memory 6000 --bfile '$pheA_testSet$' --score '$outputDir$'/'$currentPhe$'_pheA_rapido.PRS sum --out '$outputDir$'/'$currentPhe$'_pheA_sumstats'
$plink $arguments
arguments=' --memory 6000 --bfile '$pheA_testSet$' --score '$outputDir$'/'$currentPhe$'_pheB_rapido.PRS sum --out '$outputDir$'/'$currentPhe$'_pheB_sumstats'
$plink $arguments


# generate the new 
/usr/bin/python2 $smtpred$'smtpred.py' \
  --h2file $outputDir$'/ldsc_h2s.txt' \
  --rgfile $outputDir$'/ldsc_rgs.txt' \
  --nfile $outputDir$'/ldsc_ns.txt' \
  --scorefiles $outputDir$'/'$currentPhe$'_pheA_sumstats.profile' \
               $outputDir$'/'$currentPhe$'_pheB_sumstats.profile' \
  --out $outputDir$'/'
  
# this writes a final per indi PRS score file:
# $outputDir$'/multi_trait.score'
cp $outputDir$'/multi_trait.score' $outputDir$'_SMTPRED.profile'
}
export -f perform_SMPTpred # this makes local functions executable when bsubbed	




function rapidopgs_MTAG_DIAG { 
outputDir=$1
shaprsRaw=$2
avgNumIndis=$3
outPRSName=$4

pheA_sumstats=$outputDir$'_pheA_sumstats'
pheB_sumstats=$outputDir$'_pheB_sumstats'
frqFileA=$shaprsRaw$'_MAF.frq'
frqFileB=$shaprsRaw$'_MAF.frq'

# first run MTAG to produce the MTAG raw sumstats
MTAG_DIAG3 $outputDir $pheA_sumstats $pheB_sumstats $frqFileA $frqFileB "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/revision/diags/mtag_simfail_log"


}
export -f rapidopgs_MTAG_DIAG # this makes local functions executable when bsubbed



# maps my sumstats format to the format needed for MTAG, and performs MTAG, output is  $outputDir$"_mtag" (# --use_beta_se does not work in MTAG for now...)
function RUN_MTAG { 
outputDir=$1
pheA_sumstats=$2
pheB_sumstats=$3
frqFileA=$4
frqFileB=$5
# map my sumstats format:
# 1       2        3      4      5            6          7        8      9       10
#chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N
#10      100146895       rs2296438       T       C       0.6544  0.0024  0.0082  0.7694  964057
# to MTAG
#   1             2        3              4       5        6             7           8           9             
# snpid           chr     bpos            a1      a2      freq           z          pval          n
# make sure to avoid division by zero, otherwise awk may truncate file and thus create incorrect result
awk ' FNR == NR { file1[ $2 ] = $5; next; } 
NR <= NR {
if (FNR == 1) {print "snpid\tchr\tbpos\ta1\ta2\tfreq\tz\tpval\tn" }
else { if( $3 in file1 ) {

is_7_numeric = $7 + 0 == $7
is_8_numeric = $8 + 0 == $8
if(is_7_numeric && is_8_numeric) {
print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"file1[$3]"\t"$7/$8"\t"$9"\t"$10} 
}
}
}' $frqFileA $pheA_sumstats > $pheA_sumstats$'_MTAG'


awk ' FNR == NR { file1[ $2 ] = $5; next; } 
NR <= NR {
if (FNR == 1) {print "snpid\tchr\tbpos\ta1\ta2\tfreq\tz\tpval\tn" }
else { if( $3 in file1 ) {

is_7_numeric = $7 + 0 == $7
is_8_numeric = $8 + 0 == $8
if(is_7_numeric && is_8_numeric) {
print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"file1[$3]"\t"$7/$8"\t"$9"\t"$10} 
}
}
}' $frqFileB $pheB_sumstats > $pheB_sumstats$'_MTAG'


#head $pheA_sumstats$'_MTAG'

/usr/bin/python '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/revision/mtag/mtag/mtag.py' --sumstats $pheA_sumstats$'_MTAG',$pheB_sumstats$'_MTAG'  --out $outputDir --n_min 0.0 --force > $outputDir$'.log'


# convert the MTAG output into my own format 
# MTAG: 
#  1               2       3       4      5       6               7        8          9                     10                       11                     12
# SNP             CHR     BP      A1      A2      Z               N       FRQ     mtag_beta               mtag_se                 mtag_z                  mtag_pval
# rs4040617       1       779322  G       A       0.314132        10874   0.1134  0.0059077854624592615   0.030391777030177466    0.19438762848888813     0.8458723762003009

# map my sumstats format:
# 1       2        3      4      5            6          7        8      9       10
#chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N

awk '{if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN"}
else{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$8"\t"$9"\t"$10"\t"$12"\t"$7} }' $outputDir$"_trait_1.txt" > $outputDir$"_mtag"

}
export -f RUN_MTAG # this makes local functions executable when bsubbed



# creates sumstats from two PLINK GWAS performing fixed effect meta analysis, then coordinates it then adjusts them to produce a PRS
# relies on RUN_MTAG being also loaded in
function rapidopgs_MTAG_PRS { 
outputDir=$1
shaprsRaw=$2
avgNumIndis=$3
outPRSName=$4

pheA_sumstats=$outputDir$'_pheA_sumstats'
pheB_sumstats=$outputDir$'_pheB_sumstats'
frqFileA=$shaprsRaw$'_MAF.frq'
frqFileB=$shaprsRaw$'_MAF.frq'

# first run MTAG to produce the MTAG raw sumstats
RUN_MTAG $outputDir $pheA_sumstats $pheB_sumstats $frqFileA $frqFileB

# this produces $outputDir$'_mtag', with my own sumstats format

## recreate the PLINK .qassoc formatted summary stats
awk 'FNR == NR { 
b[ $3 ] = $7; 
se[ $3 ] = $8; 
p[ $3 ] = $9; 
next; } 
FNR <= NR { { 
if (FNR == 1) {print "CHR\tSNP\tBP\tNMISS\tBETA\tSE\tR2\tT\tP" } else {
if ( $2 in b) { print $1"\t"$2"\t"$3"\t"$4"\t"b[ $2 ]"\t"se[ $2 ]"\t"$7"\t"$8"\t"p[ $2 ] }  }   }
}'  $outputDir$"_mtag" $outputDir$'_pheA.qassoc' > $outputDir$'_mtag.qassoc'

# produce the final PRS
arguments='/nfs/users/nfs_m/mk23/scripts/RapidoPGS_PLINK.R '$outputDir$'_mtag.qassoc '$avgNumIndis$' '$shaprsRaw$'_MAF.frq '$outputDir$' '$outPRSName
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments

}
export -f rapidopgs_MTAG_PRS # this makes local functions executable when bsubbed




############################
# (PART 2) Train Models:
############################

# debug values
rG05_ps=( 0.5 0.75 1.0 ) 
rG05_shared_corrs=( 1.0 0.6666667 0.5 ) 

causals_all=1000
splitA=50
splitB=50
sample_size=''

currentPhe=1
p_current=0.25
shared_corr=1.0
rG=0.5
rGText='50'
current_het='extra'

rG=0.25
rGText='25'
h2=0.5
h2Text="50"



h2s=( 0.25 0.5 0.75 )
numCausals=( 1000 3000 5000 ) 



h2s=( 0.25 0.75 )
numCausals=( 1000 ) 



h2s=( 0.25 0.75 )
numCausals=( 1000 ) 




j=1
s=1
p=1
k=1
b=1
e=1
r=1


h2s=( 0.25 0.5 0.75 )
numCausals=( 3000 5000 ) 


arraylength=${#numCausals[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
causals_all=${numCausals[$j-1]}
mkdir -p $shaprsSimPhe$causals_all$'/'

# go through each phenoA/B split scenario
arraylength_split=${#pheno_splitA[@]}
for (( s=1; s<${arraylength_split}+1; s++ )); do
splitA=${pheno_splitA[$s-1]}
splitB=${pheno_splitB[$s-1]}


# go through each sample size scenario
arraylength_sammplesize=${#samplesizes[@]}
for (( p=1; p<${arraylength_sammplesize}+1; p++ )); do
sample_size=${samplesizes[$p-1]}


# go through each genetic correlation
arraylength_rG=${#rGs[@]}
for (( k=1; k<${arraylength_rG}+1; k++ )); do
rG=${rGs[$k-1]}
rGText=$(awk -v rG=$rG 'BEGIN{ print rG*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point
mkdir -p $shaprsSimPhe$causals_all$'/'$rGText$'/'

# based on RG, select the heterogeneity architecture
case $rG in
0.1)
  Ps=( "${rG01_ps[@]}" )
  shared_corrs=( "${rG01_shared_corrs[@]}" )
  echo "rG01"
  ;;
0.25)
  Ps=( "${rG025_ps[@]}" )
  shared_corrs=( "${rG025_shared_corrs[@]}" )
  echo "rG025"
  ;;
*)
  Ps=( "${rG05_ps[@]}" )
  shared_corrs=( "${rG05_shared_corrs[@]}" )
  echo "rG05"
  ;;
esac


# go through each heterogeneity architecture
arraylength_Ps=${#Ps[@]}
for (( b=1; b<${arraylength_Ps}+1; b++ )); do
p_current=${Ps[$b-1]}
shared_corr=${shared_corrs[$b-1]}
echo "rG: "$rG$" | p: "$p_current$" / shared_corr: "$shared_corr


# go through the 2 extra/regular
arraylength_heterogeneity=${#heterogeneity[@]}
for (( e=1; e<${arraylength_heterogeneity}+1; e++ )); do
current_het=${heterogeneity[$e-1]}


arraylength_h2=${#h2s[@]}
for (( r=1; r<${arraylength_h2}+1; r++ )); do
h2=${h2s[$r-1]}

h2Text=$(awk -v h2=$h2 'BEGIN{ print h2*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point


mkdir -p $shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size
# go through each pheno replicate
for ((currentPhe=1; currentPhe<=20; currentPhe++)); do #replicates



# I) Perform GWAS' on PheA/B and together
pheA_=$shaprsRaw$'/phenoA_train'$splitA$sample_size
pheB_=$shaprsRaw$'/phenoB_train'$splitB$sample_size
pheCombined_=$shaprsRaw$'/phenoB_train'$splitA$splitB$sample_size$'_all'
outputDir=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size$'/'$currentPhe
phenoAllLoc=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het$'/'$h2Text$'/'$currentPhe$'.pheno'


# if we have GWAS results on all 3
if [ -s "$outputDir"_pheCombined.qassoc ] ; then

# II) 
######################
# a) Combined runs:
avgNumIndis=$(wc -l < "$pheCombined_")
cohort='_pheCombined'
cohortFileStem=$outputDir$cohort
outPRSName=$cohort$PRS_FileStem
outfile=$cohortFileStem$PRS_FileStem 
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else
rm -rf $outputDir$'_pheCombined.out'
rm -rf $outputDir$'_pheCombined.err'
arguments='/nfs/users/nfs_m/mk23/scripts/RapidoPGS_PLINK.R '$cohortFileStem$'.qassoc '$avgNumIndis$' '$shaprsRaw$'_MAF.frq '$outputDir$' '$outPRSName
bsub -q normal -m "modern_hardware" -G team152 -n${NCORES} -J RAPID${sample_size}_${splitA}${splitB}_p${p_current}_sc${shared_corr}_${current_het}_${currentPhe} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $outputDir$'_pheCombined.out' -e $outputDir$'_pheCombined.err'  "$Rscript $arguments"


fi

#####################
# b) Sub-pheno runs: pheno A
avgNumIndis=$(wc -l < "$pheA_")
cohort='_pheA'
cohortFileStem=$outputDir$cohort
outPRSName=$cohort$PRS_FileStem
outfile=$cohortFileStem$PRS_FileStem 

if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else
rm -rf $outputDir$'_pheA.out'
rm -rf $outputDir$'_pheA.err'
arguments='/nfs/users/nfs_m/mk23/scripts/RapidoPGS_PLINK.R '$cohortFileStem$'.qassoc '$avgNumIndis$' '$shaprsRaw$'_MAF.frq '$outputDir$' '$outPRSName
bsub -q normal -m "modern_hardware" -G team152 -n${NCORES} -J RAPID${sample_size}_${splitA}_p${p_current}_sc${shared_corr}_${current_het}_${currentPhe} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $outputDir$'_pheA.out' -e $outputDir$'_pheA.err'  "$Rscript $arguments"
fi

#####################
# b) Sub-pheno runs: pheno B
avgNumIndis=$(wc -l < "$pheB_")
cohort='_pheB'
cohortFileStem=$outputDir$cohort
outPRSName=$cohort$PRS_FileStem
outfile=$cohortFileStem$PRS_FileStem 

if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
dumm=1
else

rm -rf $outputDir$'_pheB.out'
rm -rf $outputDir$'_pheB.err'

arguments='/nfs/users/nfs_m/mk23/scripts/RapidoPGS_PLINK.R '$cohortFileStem$'.qassoc '$avgNumIndis$' '$shaprsRaw$'_MAF.frq '$outputDir$' '$outPRSName
bsub -q normal -m "modern_hardware" -G team152 -n${NCORES} -J RAPID${sample_size}_${splitB}_p${p_current}_sc${shared_corr}_${current_het}_${currentPhe} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $outputDir$'_pheB.out' -e $outputDir$'_pheB.err'  "$Rscript $arguments"
fi


########################
# c) shaPRS Blended via Meta analysis:
avgNumIndis=$(wc -l < "$pheCombined_")
cohort='_pheBlend'
cohortFileStem=$outputDir$cohort
outPRSName=$cohort$'_meta_'$PRS_FileStem
outfile=$cohortFileStem$'_meta_'$PRS_FileStem 

if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile$
dumm=1
else
dumm=1
rm -rf $outputDir$'_pheBlend_meta.out'
rm -rf $outputDir$'_pheBlend_meta.err'
bsub -q normal -m "modern_hardware" -G team152 -n${NCORES} -J RAPID_META_BLEND${sample_size}_${splitA}_p${p_current}_sc${shared_corr}_${current_het}_${currentPhe} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem}  -o $outputDir$'_pheBlend_meta.out' -e $outputDir$'_pheBlend_meta.err' "rapidopgs_Blend_PRS_meta $shaprsRaw $cohortFileStem $avgNumIndis $outputDir $outPRSName"
fi





########################
# d) SMTPred:
avgNumIndis=$(wc -l < "$pheCombined_")
outfile=$outputDir$'_SMTPRED.profile'

if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile$
dumm=1
else
dumm=1
rm -rf $outputDir$'SMPTRED.out' 
rm -rf $outputDir$'SMPTRED.err'
bsub -q normal -m "modern_hardware" -G team152 -n${NCORES} -J SMPTRED${sample_size}_${splitA}_p${p_current}_sc${shared_corr}_${current_het}_${currentPhe} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem}  -o $outputDir$'SMPTRED.out' -e $outputDir$'SMPTRED.err' "perform_SMPTpred $shaprsRaw $avgNumIndis $outputDir $currentPhe $rG"
fi


########################
# e) MTAG:
avgNumIndis=$(wc -l < "$pheCombined_")
outPRSName='_mtag_'$PRS_FileStem
outfile=$outputDir$outPRSName
#rapidopgs_MTAG_DIAG $outputDir $shaprsRaw $avgNumIndis $outPRSName
#rm -rf $outfile

if [ -s "$outfile" ] ; then
#echo "MTAG EXISTS: "$outfile$
dumm=1
else
dumm=1
echo "MTAG NOT EXIST: "$outfile$
rm -rf $outputDir$'MTAG.out' 
rm -rf $outputDir$'MTAG.err'
#rapidopgs_MTAG_PRS $outputDir $shaprsRaw $avgNumIndis $outPRSName
bsub -q normal -m "modern_hardware" -G team152 -n${NCORES} -J MTAG${sample_size}_${splitA}_p${p_current}_sc${shared_corr}_${current_het}_${currentPhe} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem}  -o $outputDir$'MTAG.out' -e $outputDir$'MTAG.err' "rapidopgs_MTAG_PRS $outputDir $shaprsRaw $avgNumIndis $outPRSName"
fi



################################################
else # if GWAS has NOT yet been done
echo "NO GWAS YET for "$outputDir
#perform_GWAS $shaprsRaw $pheA_ $pheB_ $phenoAllLoc $outputDir
dumm=1
bsub -q normal -m "modern_hardware" -G team152 -n5 -J GWAS${sample_size}_${splitA}${splitB}_p${p_current}_sc${shared_corr}_${current_het}_${currentPhe} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $outputDir$'_GWAS.out' -e $outputDir$'_GWAS.err'  "perform_GWAS $shaprsRaw $pheA_ $pheB_ $phenoAllLoc $outputDir"
# this creates the following 3 files: $outputDir$'_pheCombined.qassoc' , $outputDir$'_pheA.qassoc' and $outputDir$'_pheB.qassoc'
fi


done # end of num replicates loop

done # end of h2s

done # end of extra/regular

done # end of heterogeneity architecture loop

done # end of rG loop

done # end of sample size loop

done # end of split loop

done # end of num causals loop





##################################################
# 5) LDPred, shaPRS:
####################


$shaprsRaw$'/phenoB_indis_test'
$shaprsRaw$'/phenoA_indis_test'
$shaprsRaw$'/phenoA_indis_train'
$shaprsRaw$'/phenoB_indis_train'
$shaprsRaw$'/phenoB_indis_train_all'

$shaprsSimPhe$causals_all$'/'$rGText$'/heterog/'$currentPhe$'_B.pheno'
$shaprsSimPhe$causals_all$'/'$rGText$'/heterog/'$currentPhe$'_A.pheno'
$shaprsSimPhe$causals_all$'/'$rGText$'/heterog/'$currentPhe$'phen1_2_AB.pheno'





$cohortDataFolder'ibd_all.fam'



###############
# Loop and create PRS for 
# Subpheno A from Subpheno A (subA_subA)
# SubPheno A from Combined
# Subpheno A from Blend 
# (do not care about hard thresholds as we want to illustrate the principle)
# (also do not care about SubphenoB, as these are symmetric)



#####################################
#  (PART 4) PRODUCE FINAL RESULTS:  #  
#####################################
#perform_test_evaluation $outste $PRS_FileStem $pheA_testSet $shaprsSimPhe $pheste $resultsDir

# Reusable GWAS function
function perform_test_evaluation { 
outste=$1
PRS_FileStem=$2
pheA_testSet=$3
shaprsSimPhe=$4
pheste=$5
resultsDir=$6

plink='/nfs/users/nfs_m/mk23/software/plink/plink'
plinkMem=2000  # how much RAM PLINK would use
# go through each pheno replicate
for ((currentPhe=1; currentPhe<=20; currentPhe++)); do #replicates


outputDir=$outste$currentPhe

#####################
# a) Combined
cohort='_pheCombined'
profile=$outputDir$cohort$PRS_FileStem 



# build the per individual PRS
arguments=' --memory '$plinkMem$' --bfile '$pheA_testSet$' --score '$profile$' sum --out '$profile
$plink $arguments

# check if profile does not exist, we put NA in the final results
if [ -s "$profile".profile ] ; then
# get correlation between Test set true phenos and PRS
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $profile$'.profile' FS=" " $pheste$'/'$currentPhe$'.pheno' > $outste$currentPhe$'_combined.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$outste$currentPhe$'_combined.csv '$resultsDir$'/combined 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
else
echo "NA" >> $resultsDir$'/combined'
fi

#####################
# b) Subpheno
cohort='_pheA'
profile=$outputDir$cohort$PRS_FileStem 

# build the per individual PRS
arguments=' --memory '$plinkMem$' --bfile '$pheA_testSet$' --score '$profile$' sum --out '$profile
$plink $arguments

# check if profile does not exist, we put NA in the final results
if [ -s "$profile".profile ] ; then
# get correlation between Test set true phenos and PRS
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $profile$'.profile' FS=" " $pheste$'/'$currentPhe$'.pheno' > $outste$currentPhe$'_subpheno.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$outste$currentPhe$'_subpheno.csv '$resultsDir$'/subpheno 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
else
echo "NA" >> $resultsDir$'/subpheno'
fi


#####################
# c) shaPRS Meta
cohort='_pheBlend'
profile=$outputDir$cohort$'_meta_'$PRS_FileStem 


# build the per individual PRS
arguments=' --memory '$plinkMem$' --bfile '$pheA_testSet$' --score '$profile$' sum --out '$profile
$plink $arguments

# check if profile does not exist, we put NA in the final results
if [ -s "$profile".profile ] ; then
# get correlation between Test set true phenos and PRS
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $profile$'.profile' FS=" " $pheste$'/'$currentPhe$'.pheno' > $outste$currentPhe$'_shaPRS_meta.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$outste$currentPhe$'_shaPRS_meta.csv '$resultsDir$'/shaPRS_meta 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
else
echo "NA" >> $resultsDir$'/shaPRS_meta'
fi

#####################
# d) SMPTpred
profile=$outputDir$'_SMTPRED'

# check if profile does not exist, we put NA in the final results
if [ -s "$profile".profile ] ; then
# get correlation between Test set true phenos and PRS
awk 'FNR == NR { file1[ $1 ] = $3; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $profile$'.profile' FS=" " $pheste$'/'$currentPhe$'.pheno' > $outste$currentPhe$'_SMTPred.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$outste$currentPhe$'_SMTPred.csv '$resultsDir$'/SMTPred 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
else
echo "NA" >> $resultsDir$'/SMTPred'
fi


#####################
# e) MTAG
profile=$outputDir$'_mtag_'$PRS_FileStem 

#/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/simPhe/1000/50/0.5_1.0/extra/A50_B50/size/7_mtag_RAPIDOPGS.PRS.profile

rm -rf $profile$'.profile'

# build the per individual PRS
arguments=' --memory '$plinkMem$' --bfile '$pheA_testSet$' --score '$profile$' sum --out '$profile
$plink $arguments

# check if profile does not exist, we put NA in the final results
if [ -s "$profile".profile ] ; then
# get correlation between Test set true phenos and PRS
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $profile$'.profile' FS=" " $pheste$'/'$currentPhe$'.pheno' > $outste$currentPhe$'_mtag.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$outste$currentPhe$'_mtag.csv '$resultsDir$'/_mtag 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
#echo "profile DID exist for phe"$currentPhe
else
#echo "profile did not exist for phe"$currentPhe
echo "NA" >> $resultsDir$'/_mtag'
fi


done # end of loop num replicates

}

export -f perform_test_evaluation # this makes local functions executable when bsubbed



h2s=( 0.25 0.5 0.75 )
numCausals=( 1000 3000 5000 ) 


h2s=( 0.25 0.5 0.75 )
numCausals=( 3000 5000 )


h2s=( 0.25 0.75 )
numCausals=( 1000 )



#cp -r /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/results/1000/ /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/results/1000_bk/



h2s=( 0.5 )
numCausals=( 1000 )



h2s=( 0.25 0.75 )
numCausals=( 1000 )



h2s=( 0.25 0.5 0.75 )
numCausals=( 3000 5000 )


arraylength=${#numCausals[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
causals_all=${numCausals[$j-1]}
mkdir -p $shaprsSimPhe$causals_all$'/'

# go through each phenoA/B split scenario
arraylength_split=${#pheno_splitA[@]}
for (( s=1; s<${arraylength_split}+1; s++ )); do
splitA=${pheno_splitA[$s-1]}
splitB=${pheno_splitB[$s-1]}


# go through each sample size scenario
arraylength_sammplesize=${#samplesizes[@]}
for (( p=1; p<${arraylength_sammplesize}+1; p++ )); do
sample_size=${samplesizes[$p-1]}

# go through each genetic correlation
arraylength_rG=${#rGs[@]}
for (( k=1; k<${arraylength_rG}+1; k++ )); do
rG=${rGs[$k-1]}
rGText=$(awk -v rG=$rG 'BEGIN{ print rG*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point

# based on RG, select the heterogeneity architecture
case $rG in
0.1)
  Ps=( "${rG01_ps[@]}" )
  shared_corrs=( "${rG01_shared_corrs[@]}" )
  echo "rG01"
  ;;
0.25)
  Ps=( "${rG025_ps[@]}" )
  shared_corrs=( "${rG025_shared_corrs[@]}" )
  echo "rG025"
  ;;
*)
  Ps=( "${rG05_ps[@]}" )
  shared_corrs=( "${rG05_shared_corrs[@]}" )
  echo "rG05"
  ;;
esac


# go through each heterogeneity architecture
arraylength_Ps=${#Ps[@]}
for (( b=1; b<${arraylength_Ps}+1; b++ )); do
p_current=${Ps[$b-1]}
shared_corr=${shared_corrs[$b-1]}
echo "rG: "$rG$" | p: "$p_current$" / shared_corr: "$shared_corr

# go through the 2 extra/regular
arraylength_heterogeneity=${#heterogeneity[@]}
for (( e=1; e<${arraylength_heterogeneity}+1; e++ )); do
current_het=${heterogeneity[$e-1]}

arraylength_h2=${#h2s[@]}
for (( r=1; r<${arraylength_h2}+1; r++ )); do
h2=${h2s[$r-1]}

h2Text=$(awk -v h2=$h2 'BEGIN{ print h2*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point


resultsDir=$shaprsResults$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size
echo $resultsDir
mkdir -p $resultsDir

rm -rf $resultsDir$'/combined'
rm -rf $resultsDir$'/subpheno'
rm -rf $resultsDir$'/shaPRS'
rm -rf $resultsDir$'/shaPRS_meta'
rm -rf $resultsDir$'/SMTPred'
rm -rf $resultsDir$'/_mtag'


outste=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size$'/'
pheste=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text



bsub -q normal -m "modern_hardware" -G team152 -n1 -J TEST${sample_size}_${splitA}${splitB}_p${p_current}_sc${shared_corr}_${current_het}_${currentPhe} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $outste$'_TEST.out' -e $outste$'_TEST.err'  "perform_test_evaluation $outste $PRS_FileStem $pheA_testSet $shaprsSimPhe $pheste $resultsDir"
# perform_test_evaluation $outste $PRS_FileStem $pheA_testSet $shaprsSimPhe $pheste $resultsDir




done # end of h2

done # end of extra/regular

done # end of heterogeneity architecture loop

done # end of rG loop

done # end of sample size loop



done # end of split loop

done # num causals




###########################
# Plot results
################

h2s=( 0.25 0.5 0.75 )
numCausals=( 1000 3000 5000 ) 


causals_all=5000
splitA=20
splitB=80
sample_size='_half'
p_current=1.0
shared_corr=0.5
rG=0.5
rGText='50'
currentPhe=1


 




 

arraylength_heterogeneity=${#heterogeneity[@]}
for (( e=1; e<${arraylength_heterogeneity}+1; e++ )); do
current_het=${heterogeneity[$e-1]}

arraylength_h2=${#h2s[@]}
for (( r=1; r<${arraylength_h2}+1; r++ )); do
h2=${h2s[$r-1]}

h2Text=$(awk -v h2=$h2 'BEGIN{ print h2*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point


arraylength=${#numCausals[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
causals_all=${numCausals[$j-1]}


# go through each phenoA/B split scenario
arraylength_split=${#pheno_splitA[@]}
for (( s=1; s<${arraylength_split}+1; s++ )); do
splitA=${pheno_splitA[$s-1]}
splitB=${pheno_splitB[$s-1]}



# go through each sample size scenario
arraylength_sammplesize=${#samplesizes[@]}
for (( p=1; p<${arraylength_sammplesize}+1; p++ )); do
sample_size=${samplesizes[$p-1]}

# go through each genetic correlation
arraylength_rG=${#rGs[@]}
for (( k=1; k<${arraylength_rG}+1; k++ )); do
rG=${rGs[$k-1]}
rGText=$(awk -v rG=$rG 'BEGIN{ print rG*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point

# based on RG, select the heterogeneity architecture
case $rG in
0.1)
  Ps=( "${rG01_ps[@]}" )
  shared_corrs=( "${rG01_shared_corrs[@]}" )
  echo "rG01"
  ;;
0.25)
  Ps=( "${rG025_ps[@]}" )
  shared_corrs=( "${rG025_shared_corrs[@]}" )
  echo "rG025"
  ;;
*)
  Ps=( "${rG05_ps[@]}" )
  shared_corrs=( "${rG05_shared_corrs[@]}" )
  echo "rG05"
  ;;
esac


results_params=''
results_params_temp=$results_params
# go through each heterogeneity architecture
arraylength_Ps=${#Ps[@]}
for (( b=1; b<${arraylength_Ps}+1; b++ )); do
p_current=${Ps[$b-1]}
shared_corr=${shared_corrs[$b-1]}
echo "rG: "$rG$" | p: "$p_current$" / shared_corr: "$shared_corr


resultsDir=$shaprsResults$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size'/'

# concat the X-axis labels as: "directoryNameOfResults heterogeneityArchicture"  eg: myDirc p:0.5/r:1.0
results_params_temp=$results_params$' '$resultsDir$' '$p_current$' '$shared_corr$' '$current_het # '/'$h2Text
results_params=$results_params_temp

done # end of heterogeneity architecture loop


# use a more human understandable display for the sample size than what we used 
if [[ "$sample_size" == "" ]]; then
sample_size_text="full"
else
sample_size_text=$sample_size
fi
#filename='#causals_'$causals_all$'_rG_'$rG$'_het_'$current_het$'_A'$splitA$'_B'$splitB$'_size'$sample_size_text
#plotName='#causals:_'$causals_all$',_rG:_'$rG$',_het:_'$current_het$',_A:'$splitA$'_/_B:'$splitB$',_size:'$sample_size_text

filename='causals_'$causals_all$'_h2_'$h2$'_rG_'$rG$'_het_'$current_het$'_A'$splitA$'_B'$splitB$'_size'$sample_size_text$'_revision'
plotName='causals_'$causals_all$'_h2_'$h2$'_rG:_'$rG$',_het:_'$current_het$',_A:'$splitA$'_/_B:'$splitB$',_size:'$sample_size_text

echo "---------------------------------------------------"
#echo "FILENAME: "$filename
#echo $results_params
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_sims_smtpred_meta_4.R -1 -1 '$plotName$' '$rG$' '$shaprsResults$filename$''$results_params

echo $arguments
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

echo "---------------------------------------------------"



done # end of rG loop


done # end of sample size loop

done # end of split loop

done # num causals





done # end of h2

done # end of extra/regular


#########################################################

# pull out the two representatives, but with a fixed scale from 0.095 - 0.2

arraylength_h2=${#h2s[@]}
for (( r=1; r<${arraylength_h2}+1; r++ )); do
h2=${h2s[$r-1]}

h2Text=$(awk -v h2=$h2 'BEGIN{ print h2*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point


arraylength=${#numCausals[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
causals_all=${numCausals[$j-1]}

# we have hardcoded limits for the h2=50% and 1000 causals as that go into the paper and needs to have the same scale
if [ $h2Text -eq "50" ] && [ $causals_all -eq "1000" ]; then
limitsText="0.085 0.2"
else
limitsText="-1 -1"
fi

arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_sims_smtpred_meta_4.R '$limitsText$' rG:_0.5,_het:_regular,_A:50_/_B:50,_size:full 0.5 /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/results/h2_'$h2$'_causals_'$causals_all$'rG_0.5_het_regular_A50_B50_sizefull_revision /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/results/'$causals_all$'/50/0.5_1.0/regular/'$h2Text$'/A50_B50/size/ 0.5 1.0 regular /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/results/'$causals_all$'/50/0.75_0.6666667/regular/'$h2Text$'/A50_B50/size/ 0.75 0.6666667 regular /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/results/'$causals_all$'/50/1.0_0.5/regular/'$h2Text$'/A50_B50/size/ 1.0 0.5 regular'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_sims_smtpred_meta_4.R '$limitsText$' rG:_0.5,_het:_extra,_A:50_/_B:50,_size:full 0.5 /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/results/h2_'$h2$'_causals_'$causals_all$'rG_0.5_het_extra_A50_B50_sizefull_revision /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/results/'$causals_all$'/50/0.5_1.0/extra/'$h2Text$'/A50_B50/size/ 0.5 1.0 extra /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/results/'$causals_all$'/50/0.75_0.6666667/extra/'$h2Text$'/A50_B50/size/ 0.75 0.6666667 extra /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/results/'$causals_all$'/50/1.0_0.5/extra/'$h2Text$'/A50_B50/size/ 1.0 0.5 extra'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


done

done 



###################################################


# Visualise all of the data

pheno_splitA=( 50 40 20 ) # 1= 50/50, 2=20/80
pheno_splitB=( 50 60 80 ) # 1= 50/50, 2=20/80
h2s=( 0.25 0.5 0.75 )
numCausals=( 1000 3000 5000 ) 

arraylength=${#numCausals[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
causals_all=${numCausals[$j-1]}
mkdir -p $shaprsSimPhe$causals_all$'/'

arraylength_h2=${#h2s[@]}
for (( r=1; r<${arraylength_h2}+1; r++ )); do
h2=${h2s[$r-1]}

h2Text=$(awk -v h2=$h2 'BEGIN{ print h2*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point

echo -e "combined,subpheno,shaPRS,SMTPred,MTAG" > $shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_results'
echo -e "#causals,rG,p_current,shared_corr,heterogeneity,splitA:splitB,sample_size" > $shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_predictors'



# go through each phenoA/B split scenario
arraylength_split=${#pheno_splitA[@]}
for (( s=1; s<${arraylength_split}+1; s++ )); do
splitA=${pheno_splitA[$s-1]}
splitB=${pheno_splitB[$s-1]}



# go through each sample size scenario
arraylength_sammplesize=${#samplesizes[@]}
for (( p=1; p<${arraylength_sammplesize}+1; p++ )); do
sample_size=${samplesizes[$p-1]}

# go through each genetic correlation
arraylength_rG=${#rGs[@]}
for (( k=1; k<${arraylength_rG}+1; k++ )); do
rG=${rGs[$k-1]}
rGText=$(awk -v rG=$rG 'BEGIN{ print rG*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point

# based on RG, select the heterogeneity architecture
case $rG in
0.1)
  Ps=( "${rG01_ps[@]}" )
  shared_corrs=( "${rG01_shared_corrs[@]}" )
  echo "rG01"
  ;;
0.25)
  Ps=( "${rG025_ps[@]}" )
  shared_corrs=( "${rG025_shared_corrs[@]}" )
  echo "rG025"
  ;;
*)
  Ps=( "${rG05_ps[@]}" )
  shared_corrs=( "${rG05_shared_corrs[@]}" )
  echo "rG05"
  ;;
esac

# go through each heterogeneity architecture
arraylength_Ps=${#Ps[@]}
for (( b=1; b<${arraylength_Ps}+1; b++ )); do
p_current=${Ps[$b-1]}
shared_corr=${shared_corrs[$b-1]}
echo "rG: "$rG$" | p: "$p_current$" / shared_corr: "$shared_corr


# go through the 2 extra/regular
arraylength_heterogeneity=${#heterogeneity[@]}
for (( e=1; e<${arraylength_heterogeneity}+1; e++ )); do
current_het=${heterogeneity[$e-1]}


resultsDir=$shaprsResults$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size
mkdir -p $resultsDir

# R script: create median of each, and append htem to an ongoig file: creates file with signature: combined, subpheno, shaPRS, SMTPred, MTAG
arguments='/nfs/users/nfs_m/mk23/scripts/medianator.R '$resultsDir$'/medians '$resultsDir$'/combined '$resultsDir$'/subpheno '$resultsDir$'/shaPRS_meta '$resultsDir$'/SMTPred '$resultsDir$'/_mtag'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# use a more human understandable display for the sample size than what we used 
if [[ "$sample_size" == "" ]]; then
sample_size_text="_full"
else
sample_size_text=$sample_size
fi

# signature of predictors: #causals, rg, p_current, shared_corr, current_het, splitA:splitB, sample_size
echo -e "$causals_all,$rG,$p_current,$shared_corr,$current_het,$splitA:$splitB,$sample_size_text" > $resultsDir$'/predictors' 

# append both to ongoing files to
cat $resultsDir$'/medians' >> $shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_results'

cat $resultsDir$'/predictors'  >> $shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_predictors'


done # end of split loop

done # end of extra/regular

done # end of heterogeneity architecture loop

done # end of rG loop

done # end of sample size loop


done # end of h2s

done # num causals


########################################################
# 7) HEATMAP

arraylength=${#numCausals[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
causals_all=${numCausals[$j-1]}
mkdir -p $shaprsSimPhe$causals_all$'/'

arraylength_h2=${#h2s[@]}
for (( r=1; r<${arraylength_h2}+1; r++ )); do
h2=${h2s[$r-1]}

h2Text=$(awk -v h2=$h2 'BEGIN{ print h2*100}')  # dividing by 1 to get scale to round it off: https://stackoverflow.com/questions/32797418/bash-calculation-using-bc-and-round-up-floating-point


arguments='/nfs/users/nfs_m/mk23/scripts/sims_evaluator_het6.R '$shaprsResults$causals_all$'_h2_'$h2$'sims_eval_ukbb_revision '$shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_results '$shaprsResults$'/causals'$causals_all$'_h2_'$h2$'_all_predictors'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments



done # end of h2s

done # num causals




###############################################################
# Supplementary analysis:


# I) produce a random subset of LDpred2 generated PRS
# 1) use Full set, 50:50, rG=0.5, regular, so we will have 10 entries per method
# 2) spearman rank correlation between those and rapidoPGS (also calculate the mean difference between rapido and LDPred2 for interest)




# 2) run LDpred2 on pheA,pheB and meta

# 3) run shaPRS and SM

causals_all=1000
splitA=50
h2=0.5
h2Text=50
splitB=50
sample_size=''
p_current=0.5
shared_corr=1.0
rG=0.5
rGText='50'
current_het='regular'

currentPhe=1

currentPhe=9
# 1) Select the test cases to be re-done with LDpred2
mkdir -p $shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size
# go through each pheno replicate
for ((currentPhe=1; currentPhe<=10; currentPhe++)); do #replicates


# I) Main vars
pheA_=$shaprsRaw$'/phenoA_train'$splitA$sample_size
pheB_=$shaprsRaw$'/phenoB_train'$splitB$sample_size
pheCombined_=$shaprsRaw$'/phenoB_train'$splitA$splitB$sample_size$'_all'
outputDir=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size$'/'
phenoAllLoc=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het$'/'$h2Text$'/'$currentPhe$'.pheno'


# II) 
###################### 
# a) Combined runs:
avgNumIndis=$(wc -l < "$pheCombined_")
cohort='_pheCombined'
cohortFileStem=$outputDir$currentPhe$cohort
outPRSName=$currentPhe$cohort$LDPred_FileStem
outfile=$outputDir$outPRSName 
if [ -s "$outfile" ] ; then
#echo "LDPred2 ALREADY EXISTS: "$outfile
b=1
else

arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_v4.R '$baseLDpredLoc$' '$cohortFileStem$'.qassoc '$NCORES_LDPred2$' '$avgNumIndis$' '$outputDir$' '$outPRSName$' '$pheA_testSet$'.bed'
#bsub -q long -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J combined${currentPhe} -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $outfile$'.out' -e $outfile$'.err'  "$Rscript $arguments"

fi

#####################
# b) Sub-pheno runs: pheno A
avgNumIndis=$(wc -l < "$pheA_")
cohort='_pheA'
cohortFileStem=$outputDir$currentPhe$cohort
outPRSName=$currentPhe$cohort$LDPred_FileStem
outfile=$outputDir$outPRSName 

if [ -s "$outfile" ] ; then
#echo "LDPred2 ALREADY EXISTS: "$outfile
b=1
else
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_v4.R '$baseLDpredLoc$' '$cohortFileStem$'.qassoc '$NCORES_LDPred2$' '$avgNumIndis$' '$outputDir$' '$outPRSName$' '$pheA_testSet$'.bed'
#bsub -q long -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J pheA${currentPhe} -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $outfile$'.out' -e $outfile$'.err'  "$Rscript $arguments"
fi

#####################
# c) Sub-pheno runs: pheno B
avgNumIndis=$(wc -l < "$pheB_")
cohort='_pheB'
cohortFileStem=$outputDir$currentPhe$cohort
outPRSName=$currentPhe$cohort$LDPred_FileStem
outfile=$outputDir$outPRSName 

if [ -s "$outfile" ] ; then
#echo "LDPred2 ALREADY EXISTS: "$outfile
b=1
else
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_v4.R '$baseLDpredLoc$' '$cohortFileStem$'.qassoc '$NCORES_LDPred2$' '$avgNumIndis$' '$outputDir$' '$outPRSName$' '$pheB_testSet$'.bed'
#bsub -q long -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J pheB${currentPhe} -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $outfile$'.out' -e $outfile$'.err'  "$Rscript $arguments"
fi




#####################
# d) MTAG LDpred2
#rm -rf $outfile

# get the post MTAG processed summary stats, 
outfile_raw=$outputDir$currentPhe$'_mtag'
#head $outfile_raw


# go through each chrom and run LDPred2
for ((c=1; c<=22; c++)); do


#####################################

avgNumIndis=$(wc -l < "$pheCombined_")
outPRSName=$currentPhe$'_mtag_'$LDPred_FileStem$'_'$c
outfile=$outputDir$outPRSName
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
echo MISSING ${MTAG}_${c}
rm -rf $outputDir$'MTAG_Ldpred2.out' 
rm -rf $outputDir$'MTAG_Ldpred2.err'
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_chrom.R '$baseLDpredLoc$' '$outfile_raw$' '$c$' '$avgNumIndis$' '$outputDir$' '$outPRSName$' '$pheA_testSet$' '$shaPRSscriptLoc$' 0'
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments

# bsub -q long -m "modern_hardware" -G team152 -n1 -J ${MTAG}_${c} -R "span[hosts=1] select[mem>${largerMem}] rusage[mem=${largerMem}]" -M${largerMem} -o $outputDir$'MTAG_Ldpred2_'$c$'.out' -e $outputDir$'MTAG_Ldpred2_'$c$'.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"
fi

done # end of chroms


done # end of num replicates loop

# 909606
# 8 --score: 802622 valid predictors loaded
# 9 --score: 693116 valid predictors loaded


###############################################################################################################

# WAIT UNTIL LDPRED2 JOBS FINISH ON THE CLUSTER


# 3) produce LDpred2 .profile files for each test set: pattern is: $pheA_LDPredPRS'.profile'
#  generate PRS profiles for MTAG
currentPhe=9
for ((currentPhe=1; currentPhe<=10; currentPhe++)); do #replicates

outputDir=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size$'/'

outPRSName=$currentPhe$'_mtag_'$LDPred_FileStem
outfile=$outputDir$outPRSName

# b) LDpred2
cat $outfile$'_'[!a-z] $outfile$'_'[!a-z][!a-z] > $outfile$'_all'

awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $outfile$'_all' > $outfile$'_all_no_dupes'

arguments=' --memory '$plinkMem$' --bfile '$pheA_testSet$' --score '$outfile$'_all_no_dupes sum --out '$outfile
$plink $arguments

done # end of num replicates loop




for ((currentPhe=1; currentPhe<=10; currentPhe++)); do #replicates

# I) Main vars
pheA_=$shaprsRaw$'/phenoA_train'$splitA$sample_size
pheB_=$shaprsRaw$'/phenoB_train'$splitB$sample_size
pheCombined_=$shaprsRaw$'/phenoB_train'$splitA$splitB$sample_size$'_all'
outputDir=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size$'/'
phenoAllLoc=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het$'/'$h2Text$'/'$currentPhe$'.pheno'


# II) 
###################### 
# a) Combined runs:
avgNumIndis=$(wc -l < "$pheCombined_")
cohort='_pheCombined'
cohortFileStem=$outputDir$currentPhe$cohort
outPRSName=$currentPhe$cohort$LDPred_FileStem
outfile=$outputDir$outPRSName 

# produce .profile score for the LDPRED2 betas, this is needed for SMTPred
arguments=' --memory '$plinkMem$' --bfile '$pheA_testSet$' --score '$outfile$' sum --out '$outfile
$plink $arguments


#####################
# b) Sub-pheno runs: pheno A
avgNumIndis=$(wc -l < "$pheA_")
cohort='_pheA'
cohortFileStem=$outputDir$currentPhe$cohort
outPRSName=$currentPhe$cohort$LDPred_FileStem
outfile=$outputDir$outPRSName 

# produce .profile score for the LDPRED2 betas, this is needed for SMTPred
arguments=' --memory '$plinkMem$' --bfile '$pheA_testSet$' --score '$outfile$' sum --out '$outfile
$plink $arguments

#####################
# c) Sub-pheno runs: pheno B
avgNumIndis=$(wc -l < "$pheB_")
cohort='_pheB'
cohortFileStem=$outputDir$currentPhe$cohort
outPRSName=$currentPhe$cohort$LDPred_FileStem
outfile=$outputDir$outPRSName 

# produce .profile score for the LDPRED2 betas, this is needed for SMTPred
#                                          MUST BUILD PROFILE SCORES FOR THE SAME TEST INDIVIDUALS AS PheA, we just use the sumstats trained on PheB
arguments=' --memory '$plinkMem$' --bfile '$pheA_testSet$' --score '$outfile$' sum --out '$outfile
$plink $arguments

done # end of num replicates loop




#############################################################



# 4) 
# Produce SMTPred and shaPRS final files

for ((currentPhe=1; currentPhe<=10; currentPhe++)); do #replicates

# I) Main vars
pheA_=$shaprsRaw$'/phenoA_train'$splitA$sample_size
pheB_=$shaprsRaw$'/phenoB_train'$splitB$sample_size
pheCombined_=$shaprsRaw$'/phenoB_train'$splitA$splitB$sample_size$'_all'
outputDir=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size$'/'
phenoAllLoc=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het$'/'$h2Text$'/'$currentPhe$'.pheno'


#####################
avgNumIndis=$(wc -l < "$pheA_")
cohort='_pheA'
outPRSName=$currentPhe$cohort$LDPred_FileStem
outfile=$outputDir$outPRSName 

# Grab PheB derived LDPred file
cohort_pheB='_pheB'
outPRSName_pheB=$currentPhe$cohort_pheB$LDPred_FileStem
pheB_LDPredPRS=$outputDir$outPRSName_pheB 
SMTPred_LDPred2Outdir=$outputDir$currentPhe$'_LDpred2'

cohortCombined='_pheCombined'
outPRSName_pheCombined=$currentPhe$cohortCombined$LDPred_FileStem
pheCombined_LDPredPRS=$outputDir$outPRSName_pheCombined 
shaPRS_LDPred2Outdir=$outputDir$currentPhe$'_LDpred2'
blendingFactorLoc=$outputDir$currentPhe$'_pheBlend_lFDR_meta_SNP_lFDR'

perform_SMPTpred_LDPred2_sim $SMTPred_LDPred2Outdir $outfile $pheB_LDPredPRS $outputDir$currentPhe$'_pheA_sumstats' $outputDir$currentPhe$'_pheB_sumstats'

shaPRS_LDPRED2_sim $shaPRS_LDPred2Outdir $outfile $pheCombined_LDPredPRS $blendingFactorLoc $pheA_testSet


done # end of num replicates loop



outputFolder=$shaPRS_LDPred2Outdir 
pheA_LDPredPRS=$outfile 
pheCombined_LDPredPRS=$pheCombined_LDPredPRS 
blendingFactorLoc=$blendingFactorLoc
pheA_testSet=$pheA_testSet

# creates sumstats from two PLINK GWAS performing fixed effect meta analysis, then coordinates it then adjusts them to produce a PRS
function shaPRS_LDPRED2_sim { 
outputFolder=$1
pheA_LDPredPRS=$2
pheCombined_LDPredPRS=$3
blendingFactorLoc=$4
pheA_testSet=$5

plink='/nfs/users/nfs_m/mk23/software/plink/plink'
smtpred="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/smtpred/smtpred/"
ldscDir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/smtpred/ldsc/"
w_hm3_snplist="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/ldsc/w_hm3.snplist"
ldscRefpanelLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/ldsc/eur_w_ld_chr/"


# assume that the Q-vals and the lFDR values have already been done, so we can just reuse those and proceed to blending the LDPred2 Betas


# we blend/switch between the post-LDPred2 adjusted betas, to produce the PRS we don't need anything else (eg SE)
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_PRS_processor.R '$pheA_LDPredPRS$' '$pheCombined_LDPredPRS$' '$blendingFactorLoc$' '$outputFolder$'_shaprs'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# create .profile for test set
arguments=' --memory '$plinkMem$' --bfile '$pheA_testSet$' --score '$outputFolder$'_shaprs sum --out '$outputFolder$'_shaprs'
$plink $arguments

}
export -f shaPRS_LDPRED2_sim # this makes local functions executable when bsubbed




outputFolder=$SMTPred_LDPred2Outdir
pheA_LDPredPRS=$outfile
pheB_LDPredPRS=$pheB_LDPredPRS
pheA_sumstats=$outputDir$currentPhe$'_pheA_sumstats'
pheB_sumstats=$outputDir$currentPhe$'_pheB_sumstats'

# Reusable SMTPred function
function perform_SMPTpred_LDPred2_sim { 
outputFolder=$1
pheA_LDPredPRS=$2
pheB_LDPredPRS=$3 # the 'pheB' .profile file is for the SAME indis as PheA, but trained on 'pheB'
pheA_sumstats=$4
pheB_sumstats=$5

plink='/nfs/users/nfs_m/mk23/software/plink/plink'
smtpred="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/smtpred/smtpred/"
ldscDir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/smtpred/ldsc/"
w_hm3_snplist="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/ldsc/w_hm3.snplist"
ldscRefpanelLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/ldsc/eur_w_ld_chr/"


# the filenames for the profile score files must match to the sumstats files
filename_pheA=$(awk '{idx = split(FILENAME, parts, "/"); print parts[idx]; nextfile}' $pheA_sumstats)
filename_pheB=$(awk '{idx = split(FILENAME, parts, "/"); print parts[idx]; nextfile}' $pheB_sumstats)
cp $pheA_LDPredPRS$'.profile' $outputFolder$'/'$filename_pheA$'.profile'
cp $pheB_LDPredPRS$'.profile' $outputFolder$'/'$filename_pheB$'.profile'

# obtain h2 and rG estimates for Pheno A and B via LDSC ( uses the original GWAS assoc stats)
mkdir -p $outputFolder
/usr/bin/python2  $smtpred$'ldsc_wrapper.py' \
    --out $outputFolder \
    --files $pheA_sumstats \
            $pheB_sumstats \
    --ldscpath $ldscDir \
    --snplist $w_hm3_snplist \
    --ref_ld $ldscRefpanelLoc \
    --w_ld $ldscRefpanelLoc

# generate the new .PRS, from the 2 supplied LDPred2
/usr/bin/python2 $smtpred$'smtpred.py' \
  --h2file $outputFolder$'/ldsc_h2s.txt' \
  --rgfile $outputFolder$'/ldsc_rgs.txt' \
  --nfile $outputFolder$'/ldsc_ns.txt' \
  --scorefiles $outputFolder$'/'$filename_pheA$'.profile' \
               $outputFolder$'/'$filename_pheB$'.profile' \
  --out $outputFolder$'/'
  
# this writes a final per indi PRS score file:
# $outputFolder$'/multi_trait.score'
cp $outputFolder$'/multi_trait.score' $outputFolder$'_SMTPRED.profile'

}
export -f perform_SMPTpred_LDPred2_sim # this makes local functions executable when bsubbed	


########################

causals_all=1000
splitA=50
h2=0.5
h2Text=50
splitB=50
sample_size=''
p_current=0.5
shared_corr=1.0
rG=0.5
rGText='50'
current_het='regular'


resultsDir=$shaprsResults$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size
#rm -rf $resultsDir$'/subpheno_LDpred'
#rm -rf $resultsDir$'/combined_LDpred'
#rm -rf $resultsDir$'/shaPRS_meta_LDpred'
#rm -rf $resultsDir$'/SMTPred_LDpred'
rm -rf $resultsDir$'/MTAG_LDpred'

# compute correlations between predicted and observed phenos for LDPred 
currentPhe=1
for ((currentPhe=1; currentPhe<=10; currentPhe++)); do #replicates

# I) Main vars
pheA_=$shaprsRaw$'/phenoA_train'$splitA$sample_size
pheB_=$shaprsRaw$'/phenoB_train'$splitB$sample_size
pheCombined_=$shaprsRaw$'/phenoB_train'$splitA$splitB$sample_size$'_all'
outputDir=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size$'/'
phenoAllLoc=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het$'/'$h2Text$'/'$currentPhe$'.pheno'


# Grab PheB derived LDPred file
cohort_pheB='_pheB'
outPRSName_pheB=$currentPhe$cohort_pheB$LDPred_FileStem
pheB_LDPredPRS=$outputDir$outPRSName_pheB 
SMTPred_LDPred2Outdir=$outputDir$currentPhe$'_LDpred2'

cohortCombined='_pheCombined'
outPRSName_pheCombined=$currentPhe$cohortCombined$LDPred_FileStem
pheCombined_LDPredPRS=$outputDir$outPRSName_pheCombined 
shaPRS_LDPred2Outdir=$outputDir$currentPhe$'_LDpred2'

pheste=$shaprsSimPhe$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text

# The final .profile files for LDpred2
#$outputDir$currentPhe$'_pheALDPred2.PRS.profile'
#$outputDir$currentPhe$'_pheCombinedLDPred2.PRS.profile'
#$shaPRS_LDPred2Outdir$'_shaprs.profile'
#       FID       IID  PHENO    CNT   CNT2 SCORESUM
#  4730970   4730970     -9 1811570 652664 0.340457

#$shaPRS_LDPred2Outdir$'_SMTPRED.profile' # has different signature
# FID     IID     1_pheA_sumstats
# 4730970 4730970 5.328991702762222e-05




#awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $outputDir$currentPhe$'_pheALDPred2.PRS.profile' FS=" " $pheste$'/'$currentPhe$'.pheno' > $outputDir$currentPhe$'_pheALDPred2.csv'
#arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$outputDir$currentPhe$'_pheALDPred2.csv '$resultsDir$'/subpheno_LDpred 0 0 1'
#/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $outputDir$currentPhe$'_pheCombinedLDPred2.PRS.profile' FS=" " $pheste$'/'$currentPhe$'.pheno' > $outputDir$currentPhe$'_pheCombinedLDPred2.csv'
# arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$outputDir$currentPhe$'_pheCombinedLDPred2.csv '$resultsDir$'/combined_LDpred 0 0 1'
# /nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $shaPRS_LDPred2Outdir$'_shaprs.profile' FS=" " $pheste$'/'$currentPhe$'.pheno' > $shaPRS_LDPred2Outdir$'_shaprs.csv'
# arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$shaPRS_LDPred2Outdir$'_shaprs.csv '$resultsDir$'/shaPRS_meta_LDpred 0 0 1'
# /nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# awk 'FNR == NR { file1[ $1 ] = $3; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $shaPRS_LDPred2Outdir$'_SMTPRED.profile' FS=" " $pheste$'/'$currentPhe$'.pheno' > $shaPRS_LDPred2Outdir$'_SMTPRED.csv'
# arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$shaPRS_LDPred2Outdir$'_SMTPRED.csv '$resultsDir$'/SMTPred_LDpred 0 0 1'
# /nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


outPRSName=$currentPhe$'_mtag_'$LDPred_FileStem
outfile=$outputDir$outPRSName
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $outfile$'.profile' FS=" " $pheste$'/'$currentPhe$'.pheno' > $shaPRS_LDPred2Outdir$'_MTAG.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$shaPRS_LDPred2Outdir$'_MTAG.csv '$resultsDir$'/MTAG_LDpred 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

#head $resultsDir$'/MTAG_LDpred'

done # end of num replicates loop



########################
causals_all=1000
splitA=50
h2=0.5
h2Text=50
splitB=50
sample_size=''
p_current=0.5
shared_corr=1.0
rG=0.5
rGText='50'
current_het='regular'

currentPhe=1
resultsDir=$shaprsResults$causals_all$'/'$rGText$'/'$p_current$'_'$shared_corr$'/'$current_het'/'$h2Text$'/A'$splitA$'_B'$splitB$'/size'$sample_size
echo $resultsDir


# compute absolute difference and rank correlation between LDPred2 and RapidoPGS results

# grab the same 10 r^2s for the RapidoPGS ones too
head -n 10 $resultsDir$'/subpheno' > $resultsDir$'/subpheno_comparison'
head -n 10 $resultsDir$'/combined' > $resultsDir$'/combined_comparison'
head -n 10 $resultsDir$'/shaPRS_meta' > $resultsDir$'/shaPRS_meta_comparison'
head -n 10 $resultsDir$'/SMTPred' > $resultsDir$'/SMTPred_comparison'
head -n 10 $resultsDir$'/_mtag' > $resultsDir$'/MTAG_comparison'

# paste them next to each other in the same file: first col, RapidoPG, 2nd col LDpred2
paste $resultsDir$'/subpheno_comparison' $resultsDir$'/subpheno_LDpred' > $resultsDir$'/subpheno_orig_LDpred'
paste $resultsDir$'/combined_comparison' $resultsDir$'/combined_LDpred' > $resultsDir$'/combined_orig_LDpred'
paste $resultsDir$'/shaPRS_meta_comparison' $resultsDir$'/shaPRS_meta_LDpred' > $resultsDir$'/shaPRS_meta_orig_LDpred'
paste $resultsDir$'/SMTPred_comparison' $resultsDir$'/SMTPred_LDpred' > $resultsDir$'/SMTPred_orig_LDpred'
paste $resultsDir$'/MTAG_comparison' $resultsDir$'/MTAG_LDpred' > $resultsDir$'/MTAG_orig_LDpred'

arguments='/nfs/users/nfs_m/mk23/scripts/LDPred_Rapido_Comparison_mtag.R '$resultsDir$'/'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# OLD
[1] "Spearman corr: 0.795121951219512 / p: 0.00385965957843545"
[1] "method means"
[1] "subpheno Rapido/LDpred2 : 0.109129532 / 0.101852799"
[1] "combined Rapido/LDpred2 : 0.12615408 / 0.115281442"
[1] "shaPRS Rapido/LDpred2 : 0.15077524 / 0.13083603"
[1] "SMTPred Rapido/LDpred2 : 0.104389398 / 0.110615476"
[1] "overall difference favouring Rapido: 6%"

# REVISION
[1] "Arguments received:"
[1] "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs_ukbb/results/1000/50/0.5_1.0/regular/50/A50_B50/size/"
[1] "Spearman corr: 0.781224489795918 / p: 0.123914053876823"
[1] "method means"
[1] "subpheno Rapido/LDpred2 : 0.109129532 / 0.101852799"
[1] "combined Rapido/LDpred2 : 0.12615408 / 0.115281442"
[1] "shaPRS Rapido/LDpred2 : 0.15077524 / 0.13083603"
[1] "SMTPred Rapido/LDpred2 : 0.104389398 / 0.110615476"
[1] "MTAG Rapido/LDpred2 : 0.127375489199218 / 0.139275255145665"
[1] "overall difference favouring Rapido: 3%"

# so the overall difference favouring rapido has decreased from 6 to 3%, but the spearman also decreased slightly to 0.781


$resultsDir$'/subpheno_LDpred'
$resultsDir$'/combined_LDpred'
$resultsDir$'/shaPRS_meta_LDpred'
$resultsDir$'/SMTPred_LDpred'
$resultsDir$'/MTAG_LDpred'

#cor(allData$V1,allData$V2, method = "spearman") #  0.8285714
#cor(allData$V1,allData$V2) # 0.8417679
#t.test(allData$V1,allData$V2, paired = T) # p-value = 0.04034
#mean(all_diff) # 4.75

# IE RapidoPGS is 4.8% better


#############################
198-162
# Run R script to generate pulication ready images: RevisionImageCombiner.R

# calculations to find the number of simulations
# new number of simulations total
27*6   # 162
4*9 # the  4 extra scenarios

27*6 + 4*9 # 198

27*6 - 3*3*6 # 108, the original scenarios