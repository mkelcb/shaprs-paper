####################################################
#
# shaPRS using Laura's new version of the GWAS1/2/3 dataset
# where we use the full summary statistics and expand beyond HM3
#
# screen -r -D 65157.pts-0.node-11-1-2

#########################
# static vars

# static vars
largePlinkMem=200000
bigMem=400000
ldPred2Mem=125000
largeMem=60000
plinkMem=6000  # how much RAM PLINK would use
homeBase='/nfs/users/nfs_m/mk23/'
plink=$homeBase$'software/plink/plink'
plink2=$homeBase$'software/plink/plink2'
scriptsLoc=$homeBase$'scripts/'
Rscript='/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript'
flashpca='/nfs/users/nfs_m/mk23/software/flashpca/flashpca_x86-64'
qctool2_new='/nfs/team152/mk23/software/qctool2/qctool_v2.0.6-Ubuntu16.04-x86_64/qctool'
qctoolmem=10000
shaPRSscriptLoc="/nfs/users/nfs_m/mk23/scripts/shaPRS.R"
NCORES_LDPred2=15
ldPred2Mem=125000
NCORES_LDPred2=10
plinkformat='.PHENO1.glm.logistic.hybrid'
baseLDpredLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/"
i=1
j=1
c=1
z=0
numFolds=5
numBootStraps=20
numChroms=22
rho=0.2418 # rho=0.2852094
rho_orig=0.2852094
HLA_chrom=6
HLA_start=28477797
HLA_end=33448354

baseLoc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/'

filenameStem="GWAS"
ibdLoc=$baseLoc$'ibd/'
ibdDataLoc=$ibdLoc$'data/'
ibdCommonDataLoc=$ibdDataLoc$'0_common/'
ibdCommonDataScratch=$ibdCommonDataLoc$'scratch/'
ibd_all_cleaned=$ibdCommonDataLoc$'ibd_all_cleaned' # stores all 5 cohorts merged and cleaned/filtered


covariates_filteredIBD=$ibdDataLoc$"covariates"

ibd_gwas=( 'gwas3_cd' 'gwas3_uc' ) 
revisionDir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/revisions/"
nonoverlaps_gwas=( $revisionDir'controls_UC.fam' $revisionDir'controls_CD.fam' ) 
name_gwas=( 'CD_noOverlap' 'UC_noOverlap' ) 

newResultsLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/newWave/"

newScratchLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/scratch/"

crossAncestryLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/"
traits=( 'cd' 'uc' ) 
N_effs=( '14174' '12812' ) # these were calculated from  CD cases=5400, controls=10308 //  UC cases=4647, controls=10308
hapmap3_b37bim='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/0_common/data/twas/raw/hapmap3_r1_b37_fwd_consensus.qc.poly.recode.bim' # this is the full hapmap3, that has the HLA

bcftools='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/software/bcftools/./bcftools'
prscsrefs="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsrefs/"
LABEL="1kg"

mkdir -p $newScratchLoc
mkdir -p $newResultsLoc


##########################
# 1) Convert Laura's sumstats into my own format / QC 



j=1
# loop 
arraylength=${#ibd_gwas[@]}
for (( j=1; j<${arraylength}+1; j++ )); do # 
trait=${traits[$j-1]}
N_eff=${N_effs[$j-1]}
echo $trait

#head -n 2 $covariates

# copy and decompress whole genome results
cp /lustre/scratch123/hgi/mdt2/projects/ibdgwas/IIBDGC/post_imputation/analysis/regenie/humancoreexome/$trait$'/allchr_humancoreexome_step2_'$trait$'_eur_sex_PCs_firthse.regenie.gz' $newScratchLoc$trait$'.gz'
gunzip $newScratchLoc$trait$'.gz'
head $newScratchLoc$trait

wc -l $newScratchLoc$trait # CD: 14,056,620 , UC: same

# extract variant infos into
# file_with_chr_position
#1          24928  A,G
#1          29126  T,C
awk '{if(FNR > 1) {print $1"\t"$2"\t"$4"\t"$5 } }' $newScratchLoc$trait > $newScratchLoc$trait$'_pos'
head $newScratchLoc$trait$'_pos'


# get the rsids via bcftool
/software/team152/bcftools-1.9/bcftools view -T $newScratchLoc$trait$'_pos' /lustre/scratch123/hgi/projects/ibdgwas/IIBDGC/resources/dbsnp154/GRCh38.dbSNP154_multiallelic_split.vcf.gz  | cut -f 1-5 > $newScratchLoc$trait$'_rsids'
# your_output file
#1          918108            rs77644389     G          A
wc -l $newScratchLoc$trait$'_rsids' # 20,483,071

# create a map between the B38 chr:pos reginie format to rsid, note this also uses the B37 positions from the hapmap3_b37bim
awk  'FNR == NR {  file1[$2 ] = $0; alleles1[$2]=$5":"$6; alleles2[$2]=$6":"$5; next; } 
FNR <= NR {   
lookup =  "chr"$1":"$2":"$4":"$5;
A1A2=$4":"$5
A2A1=$5":"$4
if ($3 in file1 && ( A1A2 == alleles1[$3] || A2A1 == alleles1[$3] || A2A1 == alleles2[$3] || A2A1 == alleles2[$3]  ) ) {print file1[$3]"\t"lookup}} 
' $hapmap3_b37bim $newScratchLoc$trait$'_rsids'  > $newScratchLoc$trait$'_rsids_hm3map'
head $newScratchLoc$trait$'_rsids_hm3map'
wc -l $newScratchLoc$trait$'_rsids_hm3map' # 1,251,510 out of 1,403,851


# subset to Hapmap3 and convert to _sumstats format (also only keep INFO > 0.8
#  1      2          3                    4          5           6             7           8        9           10         11           12
#CHROM    GENPOS     ID                   ALLELE0    ALLELE1     A1FREQ        INFO        TEST     BETA.Y1     SE.Y1      CHISQ.Y1     LOG10P.Y1
#21       10039788   chr21:10039788:A:G   A          G           0.00490155    0.445521    ADD      0.0392166   0.275327   0.0202881    0.052206

# map my sumstats format:
# 1       2        3      4      5            6          7        8      9       10
#chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N

awk  -v N_eff=$N_eff  'FNR == NR {  file1[$7 ] = $1"\t"$4"\t"$2; next; } 
FNR <= NR {   
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN"}

else if ($7 >= 0.8 && $3 in file1 ) {print file1[$3]"\t"$5"\t"$4"\t"$6"\t"$9"\t"$10"\t"10^(-$12)"\t"N_eff}} 
' $newScratchLoc$trait$'_rsids_hm3map' $newScratchLoc$trait  > $newScratchLoc$trait$'_sumstats'
head $newScratchLoc$trait$'_sumstats'
wc -l $newScratchLoc$trait$'_sumstats' # 1174041

# output data
newResultsLocCohort=$newResultsLoc${ibd_gwas[$j-1]}/
#echo $newResultsLocCohort

mkdir -p $newResultsLocCohort
mkdir -p $newResultsLocCohort$'bootstrap/'

produce_manhattan_ss $newScratchLoc$trait$'_sumstats' $newResultsLocCohort$'LAURA_GWAS' ${ibd_gwas[$j-1]} 7.30103
remove_ambiguous_alleles $newScratchLoc$trait$'_sumstats'

# QC
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_QC.R '$baseLDpredLoc$'/ '$newScratchLoc$trait$'_sumstats_noAmbiguousAlleles '$newResultsLocCohort${ibd_gwas[$j-1]}'_keep 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


done

################
# 2) subset the SNPs to only keep the QC passed ones from both UC and CD (to make sure the baselines won't have an advantage due to they having more SNPs

# intersect the 2 keeplists
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( $1 in sums  ) {print $0 } }
' $newResultsLoc$'gwas3_cd/gwas3_cd_keep' $newResultsLoc$'gwas3_uc/gwas3_uc_keep'  > $newResultsLoc$'UC_CD_keep'


awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $newResultsLoc$'UC_CD_keep' $newScratchLoc$'cd_sumstats_noAmbiguousAlleles'  > $newResultsLoc$'gwas3_cd/gwas3_cd_hm3'
head $newResultsLoc$'gwas3_cd/gwas3_cd_hm3'
wc -l $newResultsLoc$'gwas3_cd/gwas3_cd_hm3' # 988770


awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $newResultsLoc$'UC_CD_keep' $newScratchLoc$'uc_sumstats_noAmbiguousAlleles'  > $newResultsLoc$'gwas3_uc/gwas3_uc_hm3'
head $newResultsLoc$'gwas3_uc/gwas3_uc_hm3'
wc -l $newResultsLoc$'gwas3_uc/gwas3_uc_hm3' # 988770



##########################

# get empirical rho
GetEmpiricalRho $newResultsLoc$'newRho' $newResultsLoc$'gwas3_cd/gwas3_cd_hm3' $newResultsLoc$'gwas3_uc/gwas3_uc_hm3'
# this was 0.2418



#############################################

# EXPAND shaPRS BEYOND HM3 TO INCLUDE ALL HETEROGENEOUS SNPs:
#############################################################

# NOD2 in hg19: chr16:50,731,050 - 50,766,987

j=1
i=21
# loop 
arraylength=${#ibd_gwas[@]}
for (( j=1; j<${arraylength}+1; j++ )); do # 
trait=${traits[$j-1]}
N_eff=${N_effs[$j-1]}
echo $trait



for ((i=1; i<=$numChroms; i++)); do

# create a map between the B38 chr:pos reginie format to rsid
awk -v chrom=$i 'FNR <= NR {   
lookup =  "chr"$1":"$2":"$4":"$5;
if((index($1, "##") == 0) && $1 ==chrom)  {print $1"\t"$3"\t0\t"$2"\t"$4"\t"$5"\t"lookup}} 
' $newScratchLoc$trait$'_rsids'  > $newScratchLoc$trait$'_rsids_map_'$i
head $newScratchLoc$trait$'_rsids_map_'$i
wc -l $newScratchLoc$trait$'_rsids_map_'$i # 


# map my sumstats format:
# 1       2        3      4      5            6          7        8      9       10
#chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N

awk  -v N_eff=$N_eff  'FNR == NR {  file1[$7 ] = $1"\t"$4"\t"$2; next; } 
FNR <= NR {   
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN"}

else if ($7 >= 0.8 && $3 in file1 ) {print file1[$3]"\t"$5"\t"$4"\t"$6"\t"$9"\t"$10"\t"10^(-$12)"\t"N_eff}} 
' $newScratchLoc$trait$'_rsids_map_'$i $newScratchLoc$trait  > $newScratchLoc$trait$'_sumstats_'$i
head $newScratchLoc$trait$'_sumstats_'$i
wc -l $newScratchLoc$trait$'_sumstats_'$i # 1174041


done # end of chroms loop


done

cat  $newScratchLoc$'cd_rsids_map_'* > $newScratchLoc$'cd_rsids_map_all'
head $newScratchLoc$'cd_rsids_map_all'
wc -l $newScratchLoc$'cd_rsids_map_all'

################
# 2) subset the SNPs to only keep the QC passed ones from both UC and CD (to make sure the baselines won't have an advantage due to they having more SNPs

for ((i=1; i<=$numChroms; i++)); do

# intersect the 2 keeplists
awk 'FNR == NR { sums[ $3 ] = $3; next; } FNR <= NR { if( $3 in sums  ) {print $0 } }
' $newScratchLoc$'cd_sumstats_'$i $newScratchLoc$'uc_sumstats_'$i  > $newResultsLoc$'UC_CD_keep_'$i


awk 'FNR == NR { sums[ $3 ] = $3; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $newResultsLoc$'UC_CD_keep_'$i $newScratchLoc$'cd_sumstats_'$i  > $newScratchLoc$'cd_sumstats_'$i$'_keep'
head $newScratchLoc$'cd_sumstats_'$i$'_keep'
wc -l $newScratchLoc$'cd_sumstats_'$i$'_keep' # 988770


awk 'FNR == NR { sums[ $3 ] = $3; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $newResultsLoc$'UC_CD_keep_'$i $newScratchLoc$'uc_sumstats_'$i  > $newScratchLoc$'uc_sumstats_'$i$'_keep'
head $newScratchLoc$'uc_sumstats_'$i$'_keep'
wc -l $newScratchLoc$'uc_sumstats_'$i$'_keep' # 988770

done # end of chroms loop


##########################
# TEST RUN TO SEE HM3 COVERAGE OF HETEROGENEOUS SNPS:

mkdir -p $newScratchLoc$'gwas3_cd/'
mkdir -p $newScratchLoc$'gwas3_uc/'
for ((i=1; i<=$numChroms; i++)); do

# 3) run MTAG/SHAPRS
RUN_MTAG $newScratchLoc$'gwas3_cd/CD_MTAG_'$i $newScratchLoc$'cd_sumstats_'$i$'_keep' $newScratchLoc$'uc_sumstats_'$i$'_keep'

RUN_MTAG $newScratchLoc$'gwas3_uc/UC_MTAG_'$i $newScratchLoc$'uc_sumstats_'$i$'_keep' $newScratchLoc$'cd_sumstats_'$i$'_keep'


shaPRS_fixed $newScratchLoc$'gwas3_cd/CD_SHAPRS_'$i $newScratchLoc$'cd_sumstats_'$i$'_keep' $newScratchLoc$'uc_sumstats_'$i$'_keep' $rho 1

shaPRS_fixed $newScratchLoc$'gwas3_uc/UC_SHAPRS_'$i $newScratchLoc$'uc_sumstats_'$i$'_keep' $newScratchLoc$'cd_sumstats_'$i$'_keep' $rho 1

done # end of chroms loop


# concat results
head -n 1 $newScratchLoc$'cd_sumstats_1_keep' > $newScratchLoc$'cd_sumstats_all'
tail -n +2 -q $newScratchLoc$'cd_sumstats_'*'_keep' >> $newScratchLoc$'cd_sumstats_all'
head -n 1 $newScratchLoc$'uc_sumstats_1_keep' > $newScratchLoc$'uc_sumstats_all'
tail -n +2 -q $newScratchLoc$'uc_sumstats_'*'_keep' >> $newScratchLoc$'uc_sumstats_all'


head -n 1 $newScratchLoc$'gwas3_cd/CD_SHAPRS_1_lFDR_meta_SNP_lFDR' > $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all'
tail -n +2 -q $newScratchLoc$'gwas3_cd/CD_SHAPRS_'*'_lFDR_meta_SNP_lFDR' >> $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all'
head -n 1 $newScratchLoc$'gwas3_uc/UC_SHAPRS_1_lFDR_meta_SNP_lFDR' > $newScratchLoc$'gwas3_uc/UC_SHAPRS_lFDR_meta_SNP_lFDR_all'
tail -n +2 -q $newScratchLoc$'gwas3_uc/UC_SHAPRS_'*'_lFDR_meta_SNP_lFDR' >> $newScratchLoc$'gwas3_uc/UC_SHAPRS_lFDR_meta_SNP_lFDR_all'

head -n 1 $newScratchLoc$'gwas3_cd/CD_SHAPRS_1_shaprs' > $newScratchLoc$'gwas3_cd/CD_SHAPRS_shaprs_all'
tail -n +2 -q $newScratchLoc$'gwas3_cd/CD_SHAPRS_'*'_shaprs' >> $newScratchLoc$'gwas3_cd/CD_SHAPRS_shaprs_all'
head -n 1 $newScratchLoc$'gwas3_uc/UC_SHAPRS_1_shaprs' > $newScratchLoc$'gwas3_uc/UC_SHAPRS_shaprs_all'
tail -n +2 -q $newScratchLoc$'gwas3_uc/UC_SHAPRS_'*'_shaprs' >> $newScratchLoc$'gwas3_uc/UC_SHAPRS_shaprs_all'


#head -n 1 $newScratchLoc$'gwas3_cd/CD_MTAG_1_mtag' > $newScratchLoc$'gwas3_cd/CD_MTAG_mtag_all'
#tail -n +2 -q $newScratchLoc$'gwas3_cd/CD_MTAG_'*'_mtag' >> $newScratchLoc$'gwas3_cd/CD_MTAG_mtag_all'
#head -n 1 $newScratchLoc$'gwas3_uc/UC_MTAG_1_mtag' > $newScratchLoc$'gwas3_uc/UC_MTAG_mtag_all'
#tail -n +2 -q $newScratchLoc$'gwas3_uc/UC_MTAG_'*'_mtag' >> $newScratchLoc$'gwas3_uc/UC_MTAG_mtag_all'



# build table of 
# RSID p_CD p_UC

# my sumstats format:
# 1       2        3      4      5            6          7        8      9       10
#chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N
awk '{if (FNR == 1 || $9 <  5 * 10e-8) {print $0}}' $newScratchLoc$'cd_sumstats_all' > $newScratchLoc$'cd_sumstats_all_sig'
awk '{if (FNR == 1 || $9 <  5 * 10e-8) {print $0}}' $newScratchLoc$'uc_sumstats_all' > $newScratchLoc$'uc_sumstats_all_sig'
head $newScratchLoc$'cd_sumstats_all_sig'
wc -l $newScratchLoc$'cd_sumstats_all_sig' # 2398

head $newScratchLoc$'uc_sumstats_all_sig'
wc -l $newScratchLoc$'uc_sumstats_all_sig' # 3491


# subset to NOD2
awk '{if($1 == 16 && $2 >= 50731050 && $2 <= 50766987) {print $0}}' $newScratchLoc$'cd_sumstats_all_sig' > $newScratchLoc$'cd_sumstats_all_sig_NOD2'
awk '{if($1 == 16 && $2 >= 50731050 && $2 <= 50766987) {print $0}}' $newScratchLoc$'uc_sumstats_all_sig' > $newScratchLoc$'uc_sumstats_all_sig_NOD2'
head $newScratchLoc$'cd_sumstats_all_sig_NOD2'
wc -l $newScratchLoc$'cd_sumstats_all_sig_NOD2' # 23

head $newScratchLoc$'uc_sumstats_all_sig_NOD2'
wc -l $newScratchLoc$'uc_sumstats_all_sig_NOD2' # 0

# RSID lFDR<0.5
awk '{if (FNR == 1 || $2 < 0.5) print $0}' $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all' > $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig'
head $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig'
wc -l $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig' # 9101

# Also check for lFDR < 1, IE any other SNP that COULD be heterogeneous
awk '{if (FNR == 1 || $2 < 1) print $0}' $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all' > $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig2'
head $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig2'
wc -l $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig2' # 61,338


# subset to NOD2
awk 'FNR == NR { if($1 == 16 && $2 >= 50731050 && $2 <= 50766987) {sums[ $3 ] = $3;} next; } FNR <= NR { if( FNR == 1 || $1 in sums  ) {print $0 } }
' $newScratchLoc$'cd_sumstats_all_sig' $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig'  > $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_sig_NOD2'
head $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_sig_NOD2'
wc -l $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_sig_NOD2' # 23



#########
#HM3
awk '{if (FNR == 1 || $9 <  5 * 10e-8) {print $0}}' $newResultsLoc$'gwas3_cd/gwas3_cd_hm3' > $newScratchLoc$'cd_sumstats_hm3_sig'
awk '{if (FNR == 1 || $9 <  5 * 10e-8) {print $0}}' $newResultsLoc$'gwas3_uc/gwas3_uc_hm3' > $newScratchLoc$'uc_sumstats_hm3_sig'
head $newScratchLoc$'cd_sumstats_hm3_sig'
wc -l $newScratchLoc$'cd_sumstats_hm3_sig' # 443

head $newScratchLoc$'uc_sumstats_hm3_sig'
wc -l $newScratchLoc$'uc_sumstats_hm3_sig' # 288

# RSID lFDR<0.5
awk '{if (FNR == 1 || $2 < 0.5) print $0}' $newResultsLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR' > $newScratchLoc$'gwas3_cd/CD_SHAPRS_hm3_lFDR_meta_SNP_lFDR_sig'
head $newScratchLoc$'gwas3_cd/CD_SHAPRS_hm3_lFDR_meta_SNP_lFDR_sig'
wc -l $newScratchLoc$'gwas3_cd/CD_SHAPRS_hm3_lFDR_meta_SNP_lFDR_sig' # 652

awk '{if (FNR == 1 || $2 < 1) print $0}' $newResultsLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR' > $newScratchLoc$'gwas3_cd/CD_SHAPRS_hm3_lFDR_meta_SNP_lFDR_sig2'
head $newScratchLoc$'gwas3_cd/CD_SHAPRS_hm3_lFDR_meta_SNP_lFDR_sig2'
wc -l $newScratchLoc$'gwas3_cd/CD_SHAPRS_hm3_lFDR_meta_SNP_lFDR_sig2' # 15,179


# subset to NOD2
awk 'FNR == NR { if($1 == 16 && $2 >= 50731050 && $2 <= 50766987) {sums[ $3 ] = $3;} next; } FNR <= NR { if( FNR == 1 || $1 in sums  ) {print $0 } }
' $newScratchLoc$'cd_sumstats_hm3_sig' $newScratchLoc$'gwas3_cd/CD_SHAPRS_hm3_lFDR_meta_SNP_lFDR_sig'  > $newScratchLoc$'gwas3_cd/CD_SHAPRS_hm3_lFDR_meta_SNP_lFDR_sig_NOD2'
head $newScratchLoc$'gwas3_cd/CD_SHAPRS_hm3_lFDR_meta_SNP_lFDR_sig_NOD2'
wc -l $newScratchLoc$'gwas3_cd/CD_SHAPRS_hm3_lFDR_meta_SNP_lFDR_sig_NOD2' # 11


# OK so I've run through All of the ~14mil SNPs to check how much we loose due to restricting to HM3:
#                 ALL    HM3
# p_CD<5*10e-8    2397   442
# p_UC<5*10e-8    3490   287
# lFDR<0.5        9100   651
# NOD2_lFDR<0.5   22     10  

# So we loose ~93% of heterogeneous in total, and 55% of NOD2, which looks like a lot.
# This of course may be misleading as a lot of these SNPs are in the same locus, so once we factor in LD, we may not loose that much actual effect size as far fewer SNPs could be tagging this well. The only way to be sure is to perform the full analysis.


###########################################################

#subset sumstats to hm3 + all lfdr < 1

# first get a list of SNPs that are NOT in hapmap3
awk  'FNR == NR {  file1[$2 ] = $2; next; }
FNR <= NR {   
if ($1 in file1 == 0 || FNR ==1 ) {print $0}} 
' $hapmap3_b37bim $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig2' > $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig2_nothm3'
head $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig2_nothm3'
wc -l $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig2_nothm3' # 54,085, and yes for 7254

awk  'FNR == NR {  print $2; next; }
FNR <= NR {   
if (FNR  > 1 ) {print $1}} 
' $hapmap3_b37bim $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig2_nothm3' > $newScratchLoc$'hm3_heterog'
head $newScratchLoc$'hm3_heterog'
wc -l $newScratchLoc$'hm3_heterog' # 1,457,935


awk  'FNR == NR {  file1[$1 ] = $1; next; }
FNR <= NR {   
if (FNR  == 1 || $3 in  file1) {print $0}} 
' $newScratchLoc$'hm3_heterog' $newScratchLoc$'cd_sumstats_all' > $newScratchLoc$'cd_sumstats_hm3_het'
head $newScratchLoc$'cd_sumstats_hm3_het'
wc -l $newScratchLoc$'cd_sumstats_hm3_het' # 1,227,343

awk  'FNR == NR {  file1[$1 ] = $1; next; }
FNR <= NR {   
if (FNR  == 1 || $3 in  file1) {print $0}} 
' $newScratchLoc$'hm3_heterog' $newScratchLoc$'uc_sumstats_all' > $newScratchLoc$'uc_sumstats_hm3_het'
head $newScratchLoc$'uc_sumstats_hm3_het'
wc -l $newScratchLoc$'uc_sumstats_hm3_het' #  1,227,342


######################################
# create new test set and LDref panels and new QC panels

# (test set must come first, as if SNP does not exist in target then we may keep SNPs in the LDref step that don't exist in the test set
# whereas we could have kept SNPs that were perhaps less well tagging, but present )

# intersect the final sumstats with what we would have in the test set
awk 'FNR == NR { sums[ $3 ] = $3; next; } FNR <= NR { if( $2 in sums  ) {print $2 } }
' $newScratchLoc$'cd_sumstats_hm3_het' $ibdDataLoc$'gwas1/ibd_all.bim'  > $newScratchLoc$'cd_sumstats_hm3_het_testkeep'
head $newScratchLoc$'cd_sumstats_hm3_het_testkeep'
wc -l $newScratchLoc$'cd_sumstats_hm3_het_testkeep' # 1,110,718


# find out how many of these were the heterogenous ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( $1 in sums  ) {print $1 } }
' $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig2_nothm3' $newScratchLoc$'cd_sumstats_hm3_het_testkeep'  > $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het'
head $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het'
wc -l $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het' # 24,125 # only 24K


awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( $1 in sums  ) {print $1 } }
' $newScratchLoc$'gwas3_cd/CD_SHAPRS_hm3_lFDR_meta_SNP_lFDR_sig2' $newScratchLoc$'cd_sumstats_hm3_het_testkeep'  > $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_old'
head $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_old'
wc -l $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_old' # 13170... so we have 13K anyway???

# what about the <0.5 ones?
# can we find GWAS1/GWAS2 on the new panel?

head $ibdDataLoc$'gwas1/ibd_all.fam' # WTCCC66061_WTCCC66061
tail $ibdDataLoc$'gwas1/ibd_all.fam' # WTCCC189350_WTCCC189350
head $ibdDataLoc$'gwas2/ibd_all.fam' # BLOOD291637
tail $ibdDataLoc$'gwas2/ibd_all.fam' # WTCCCT543557

# need to split IDs for GWAS1
#  $ibdDataLoc$'gwas1/ibd_all.fam'

awk '{
n=split($2,a,/_/);  
print a[1]"\t"a[1]"\t"$3"\t"$4"\t"$5"\t"$6
}' $ibdDataLoc$'gwas1/ibd_all.fam' > $ibdDataLoc$'gwas1/ibd_all_splitId.fam'
head $ibdDataLoc$'gwas1/ibd_all_splitId.fam'


mkdir -p $ibdDataLoc$'gwas1_new/'
mkdir -p $ibdDataLoc$'gwas2_new/'



# map back the rsids to the ids we actually have in the PLINK files
$newScratchLoc$'cd_sumstats_hm3_het'
$newScratchLoc$'cd_rsids_map_all'
awk 'FNR == NR { sums[ $2 ] = $7; next; } FNR <= NR { if( $3 in sums  ) {print sums[$3] } }
' $newScratchLoc$'cd_rsids_map_all' $newScratchLoc$'cd_sumstats_hm3_het'  > $newScratchLoc$'cd_sumstats_hm3_het_pvarid'
head $newScratchLoc$'cd_sumstats_hm3_het_pvarid'
wc -l $newScratchLoc$'cd_sumstats_hm3_het_pvarid'


# GWAS1
# /lustre/scratch123/hgi/mdt2/projects/ibdgwas/IIBDGC/post_imputation/affymetrix500/phenotype_data/affymetrix500_all_studies_merged_eur_phenotype_ibd
# 4642, out of 4115

# GWAS2 
# /lustre/scratch123/hgi/mdt2/projects/ibdgwas/IIBDGC/post_imputation/affymetrix6/phenotype_data/affymetrix6_all_studies_merged_eur_phenotype_ibd
# 11068 out of 4694


i=21
# filter to keep the same set of the fam's we have, and also the SNPs we need
for ((i=1; i<=$numChroms; i++)); do

# GWAS1
pgen="/lustre/scratch123/hgi/mdt2/projects/ibdgwas/IIBDGC/post_imputation/affymetrix500/imputed_data/affymetrix500_chr"$i$"_subset"
arguments=' --pfile '$pgen$' --memory  '$plinkMem$' --threads 16 --make-bed --keep '$ibdDataLoc$'gwas1/ibd_all_splitId.fam --extract '$newScratchLoc$'cd_sumstats_hm3_het_pvarid --out '$ibdDataLoc$'gwas1_new/'$i$' --allow-extra-chr --snps-only --bp-space 1'
$plink2 $arguments

# GWAS2
pgen="/lustre/scratch123/hgi/mdt2/projects/ibdgwas/IIBDGC/post_imputation/affymetrix6/imputed_data/affymetrix6_chr"$i$"_subset"
arguments=' --pfile '$pgen$' --memory  '$plinkMem$' --threads 16 --make-bed --keep '$ibdDataLoc$'gwas2/ibd_all.fam --extract '$newScratchLoc$'cd_sumstats_hm3_het_pvarid --out '$ibdDataLoc$'gwas2_new/'$i$' --allow-extra-chr --snps-only --bp-space 1'
$plink2 $arguments
done


# merge chroms, and also overwrite the old var ids with rsids
plinkFileList=$ibdCommonDataScratch$'mergelist.txt'
plinkFileList2=$ibdCommonDataScratch$'mergelist2.txt'
echo -n > $plinkFileList
for ((i=1; i<=$numChroms; i++)); do
echo $ibdDataLoc$'gwas1_new/'$i >> ${plinkFileList}
echo $ibdDataLoc$'gwas2_new/'$i >> ${plinkFileList2}
done

arguments=' --memory '$plinkMem$' --merge-list '$plinkFileList$' --make-bed --out '$ibdDataLoc$'gwas1_new/ibd_all --allow-extra-chr --allow-no-sex'
$plink $arguments

arguments=' --memory '$plinkMem$' --merge-list '$plinkFileList2$' --make-bed --out '$ibdDataLoc$'gwas2_new/ibd_all --allow-extra-chr --allow-no-sex'
$plink $arguments
# 4673 people, 922331 variants ??? why 922331

# restore rsids
awk 'FNR == NR { sums[ $7 ] = $2; next; } FNR <= NR { if( $2 in sums  ) {print $1"\t"sums[$2]"\t"$3"\t"$4"\t"$5"\t"$6 } else {print $0} }
' $newScratchLoc$'cd_rsids_map_all' $ibdDataLoc$'gwas1_new/ibd_all.bim'  > $ibdDataLoc$'gwas1_new/ibd_all_rsid.bim'
head $ibdDataLoc$'gwas1_new/ibd_all_rsid.bim'
wc -l $ibdDataLoc$'gwas1_new/ibd_all_rsid.bim'
cp $ibdDataLoc$'gwas1_new/ibd_all.bim' $ibdDataLoc$'gwas1_new/ibd_old.bim'
cp $ibdDataLoc$'gwas1_new/ibd_all_rsid.bim' $ibdDataLoc$'gwas1_new/ibd_all.bim'

awk 'FNR == NR { sums[ $7 ] = $2; next; } FNR <= NR { if( $2 in sums  ) {print $1"\t"sums[$2]"\t"$3"\t"$4"\t"$5"\t"$6 } else {print $0} }
' $newScratchLoc$'cd_rsids_map_all' $ibdDataLoc$'gwas2_new/ibd_all.bim'  > $ibdDataLoc$'gwas2_new/ibd_all_rsid.bim'
head $ibdDataLoc$'gwas2_new/ibd_all_rsid.bim'
wc -l $ibdDataLoc$'gwas2_new/ibd_all_rsid.bim'
cp $ibdDataLoc$'gwas2_new/ibd_all.bim' $ibdDataLoc$'gwas2_new/ibd_old.bim'
cp $ibdDataLoc$'gwas2_new/ibd_all_rsid.bim' $ibdDataLoc$'gwas2_new/ibd_all.bim'

# the above positions are still in B38
# need to map these back to B37


##################
# Diagnostics

# check how many of the lFDR < 1 ones we found
awk 'FNR == NR { sums[ $3 ] = $3; next; } FNR <= NR { if( $2 in sums  ) {print $2 } }
' $newScratchLoc$'cd_sumstats_hm3_het' $ibdDataLoc$'gwas1_new/ibd_all.bim'  > $newScratchLoc$'cd_sumstats_hm3_het_testkeep_new'
head $newScratchLoc$'cd_sumstats_hm3_het_testkeep_new'
wc -l $newScratchLoc$'cd_sumstats_hm3_het_testkeep_new' # 933127

awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( $1 in sums  ) {print $1 } }
' $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig2_nothm3' $newScratchLoc$'cd_sumstats_hm3_het_testkeep_new'  > $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_new'
head $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_new'
wc -l $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_new' # 41,380, OK, so we could find 2x as many

# check lFDR < 0.5 ones we found
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( $1 in sums  ) {print $1 } }
' $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig' $newScratchLoc$'cd_sumstats_hm3_het_testkeep_new'  > $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_new_05'
head $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_new_05'
wc -l $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_new_05' # 6709 of the lFDR < 0.5 ones

# check how many of these would we have found in the original GWAS1 datasets
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( $2 in sums  ) {print $2 } }
' $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig' $ibdDataLoc$'gwas1/ibd_all.bim'  > $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_orig'
head $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_orig'
wc -l $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_orig' # 1,844 (so around 3-4x as many in the new one)

# how many in the original new wave run?
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( $3 in sums  ) {print $3 } }
' $newScratchLoc$'gwas3_cd/CD_SHAPRS_lFDR_meta_SNP_lFDR_all_sig' $newResultsLoc$'gwas3_cd/CD_SHAPRS_shaprs'  > $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_neww'
head $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_neww'
wc -l $newScratchLoc$'cd_sumstats_hm3_het_testkeep_het_neww' # 974, ok so that is quite a bit less


############# Diagnostics end

###########################################################################
# I. Liftover
###########################################################################

mkdir -p $prscsrefs$'liftover/'
cd $prscsrefs$'liftover/'
believeScratchLoc=$prscsrefs$'liftoverDir/scratch/'
mkdir -p $believeScratchLoc


# liftover
# https://www.biostars.org/p/252938/
# download database and script
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget https://raw.githubusercontent.com/Shicheng-Guo/GscPythonUtility/master/liftOverPlink.py
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/ibdqc.pl

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz


head -n 2 $ibdDataLoc$'gwas1_new/ibd_all.fam' > $believeScratchLoc$'indi_keeplist_2'

# create a dummy plink file with 2 indis, as we will be transcoding to text plink files
# rebuild plink file to avoid chromsome-miss-order problem
arguments=' --memory 43000 --bfile '$ibdDataLoc$'gwas1_new/ibd_all   --keep-allele-order --make-bed --out '$believeScratchLoc$'liftOverDummy.sort --allow-extra-chr --allow-no-sex'
$plink $arguments

# space to tab to generate bed files for liftOver from hg18 to hg19
$plink --bfile $believeScratchLoc$'liftOverDummy.sort' --recode tab --out $believeScratchLoc$'liftOverDummy.sort.tab'

# apply liftOverPlink.py to update hg18 to hg19 or hg38
# Only works in BIRC10, Not HPC, Caused by Python version
#mkdir liftOver
chmod 0777 liftOver
python2 liftOverPlink.py -m $believeScratchLoc$'liftOverDummy.sort.tab.map' -p $believeScratchLoc$'liftOverDummy.sort.tab.ped' -o $believeScratchLoc$'liftOverDummy.hg19' -c hg38ToHg19.over.chain.gz -e ./liftOver

# convert back to binary bim
arguments=' --memory 43000 --file '$believeScratchLoc$'liftOverDummy.hg19  --keep-allele-order --make-bed --out '$believeScratchLoc$'liftedOver --allow-extra-chr --allow-no-sex'
$plink $arguments
head $believeScratchLoc$'liftedOver.bim' 
wc -l $believeScratchLoc$'liftedOver.bim' # 932792  out of 933127


# because we found fewer SNPs than total, we need to subset the PLINK files, otherwise they would be misaligned when we apply the liftover/
awk '{print $2}' $believeScratchLoc$'liftedOver.bim'  > $believeScratchLoc$'liftedOver_extract' 

arguments=' --memory 43000 --bfile '$ibdDataLoc$'gwas1_new/ibd_all  --keep-allele-order --extract '$believeScratchLoc$'liftedOver_extract --make-bed --out '$ibdDataLoc$'gwas1_new/ibd_all_hg19 --allow-extra-chr --allow-no-sex'
$plink $arguments
arguments=' --memory 43000 --bfile '$ibdDataLoc$'gwas2_new/ibd_all  --keep-allele-order --extract '$believeScratchLoc$'liftedOver_extract --make-bed --out '$ibdDataLoc$'gwas2_new/ibd_all_hg19 --allow-extra-chr --allow-no-sex'
$plink $arguments


# overwrite all pos with the hg19 ones
cp $ibdDataLoc$'gwas1_new/ibd_all_hg19.bim' $ibdDataLoc$'gwas1_new/ibd_all_hg19.b38'
cp $ibdDataLoc$'gwas2_new/ibd_all_hg19.bim' $ibdDataLoc$'gwas2_new/ibd_all_hg19.b38'
cp $newScratchLoc$'cd_sumstats_hm3_het' $newScratchLoc$'cd_sumstats_hm3_het.b38'
cp $newScratchLoc$'uc_sumstats_hm3_het' $newScratchLoc$'uc_sumstats_hm3_het.b38'




# .bim files CD / UC
awk 'FNR == NR { file1[ $2 ] = $4;  next; }
FNR <= NR {  if( $2 in file1 ) {print $1"\t"$2"\t"$3"\t"file1[$2]"\t"$5"\t"$6} }
' $believeScratchLoc$'liftedOver.bim' $ibdDataLoc$'gwas1_new/ibd_all_hg19.b38' > $ibdDataLoc$'gwas1_new/ibd_all_hg19.bim'
head $ibdDataLoc$'gwas1_new/ibd_all_hg19.bim'

awk 'FNR == NR { file1[ $2 ] = $4;  next; }
FNR <= NR {  if( $2 in file1 ) {print $1"\t"$2"\t"$3"\t"file1[$2]"\t"$5"\t"$6} }
' $believeScratchLoc$'liftedOver.bim' $ibdDataLoc$'gwas2_new/ibd_all_hg19.b38' > $ibdDataLoc$'gwas2_new/ibd_all_hg19.bim'
head $ibdDataLoc$'gwas2_new/ibd_all_hg19.bim'

# sumstats CD / UC
awk 'FNR == NR { file1[ $2 ] = $4;  next; }
FNR <= NR {  if (FNR == 1) {print $0 }else if( $3 in file1 ) {print $1"\t"file1[$3]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10} }
' $believeScratchLoc$'liftedOver.bim' $newScratchLoc$'cd_sumstats_hm3_het.b38' > $newScratchLoc$'cd_sumstats_hm3_het'
head $newScratchLoc$'cd_sumstats_hm3_het'
wc -l $newScratchLoc$'cd_sumstats_hm3_het'

awk 'FNR == NR { file1[ $2 ] = $4;  next; }
FNR <= NR {  if (FNR == 1) {print $0 }else if( $3 in file1 ) {print $1"\t"file1[$3]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10} }
' $believeScratchLoc$'liftedOver.bim' $newScratchLoc$'uc_sumstats_hm3_het.b38' > $newScratchLoc$'uc_sumstats_hm3_het'
head $newScratchLoc$'uc_sumstats_hm3_het'
wc -l $newScratchLoc$'uc_sumstats_hm3_het'


# need to create fam files that have the binary phenos 
# (for gwas1, we need to split again the IDs to be able to match them)
awk 'FNR == NR { n=split($2,a,/_/); file1[ a[1] ] = $6;  next; }
FNR <= NR {  {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"file1[$2]} }
' $ibdCommonDataScratch$'gwas1/all_found.fam'  $ibdDataLoc$'gwas1_new/ibd_all_hg19.fam' > $ibdDataLoc$'gwas1_new/all_found.fam'

awk 'FNR == NR {  file1[ $2 ] = $6;  next; }
FNR <= NR {  {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"file1[$2]} }
' $ibdCommonDataScratch$'gwas2/all_found.fam'  $ibdDataLoc$'gwas2_new/ibd_all_hg19.fam' > $ibdDataLoc$'gwas2_new/all_found.fam'




##################################

# Generate the PRS-CS LDref panel from the custom SNPs

# the 2 scripts from here: calc_ld.sh write_ldblk.py: https://github.com/getian107/PRScs/issues/20
# the main instructions from here: https://github.com/getian107/PRScs/issues/13
# the above is a bit vague on what files will be actually needed, here is a worked example for Finnish people:
# https://github.com/FINNGEN/CS-PRS-pipeline/blob/master/scripts/panel/finngen.sh
# https://github.com/FINNGEN/CS-PRS-pipeline/blob/master/scripts/panel/write_ldblk.py

# 1. get LD block boundaries: https://www.biorxiv.org/content/10.1101/020255v1.full.pdf, 1,703 LD blocks

cd $prscsrefs
wget https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/EUR/fourier_ls-all.bed
wget https://raw.githubusercontent.com/FINNGEN/CS-PRS-pipeline/master/scripts/panel/write_ldblk.py
# 2. create snpinfo_ukbb_hm3 file, that splits the test set SNPs into the ldblocks
head $prscsrefs$'fourier_ls-all.bed'


rm -rf $prscsrefs$'ibd/'
mkdir -p $prscsrefs$'ibd/'
mkdir -p $prscsrefs$'ibd/ldblk/'

blk=0
counter=0
rm -rf $prscsrefs$'ibd/blk_size'
rm -rf $prscsrefs$'ibd/blk_chr'
echo -e "CHR\tSNP\tBP\tA1\tA2\tMAF" > $prscsrefs$'ibd/ldblk/snpinfo_'$LABEL$'_hm3'

while read p; do # go through each block
p=$(echo $p | sed -e 's/\r//g')  # the last column contains a newline character, which would fuck up awk string comparison and create unreadable files etc
IFS=' ' read -r -a arrIN <<< "$p"  # https://stackoverflow.com/questions/10586153/split-string-into-an-array-in-bash
chr=${arrIN[0]} 
start=${arrIN[1]}
stop=${arrIN[2]}
chrnum=$(awk -v chr=$chr 'BEGIN { n=split(chr,a,/chr/);   print a[2] }')
counter=$((counter+1))

if [ $counter != "1" ]; then  # skip header

blk=$((blk+1))
echo "processing block"$blk$" "$chr$" start: "$start$" stop: "$stop

# create SNP list that falls within the block
awk -v chr=$chrnum -v start=$start -v stop=$stop  '{if ($1 == chr && $4 >= start && $4 < stop) {print $2}}' $ibdDataLoc$'gwas1_new/ibd_all_hg19.bim' > $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.snplistOrig



if [ -s $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.snplistOrig ]; then
#echo "The file is NOT empty."

#awk -v chr=$chrnum -v start=$start -v stop=$stop  '{if ($1 == chr && $4 >= start && $4 <= stop) {print $2}}' $ibdDataLoc$'gwas1_new/ibd_all_hg19.bim' > $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.snplist
#$plink  --bfile $ibdDataLoc$'gwas1_new/ibd_all_hg19'  --keep-allele-order  --extract $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.snplist  --r square  --freq  --out $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL} --memory 2000
# instead of the above, we extract based on from-to    --chr ${chrnum} --from-bp ${start} --to-bp ${stop} 
# we filter on maf too, as otherwise if we get 0 MAF SNPs those will create NANs in the LD matrix
# Also filter on max maf 0.475. This is because if MAF is near 0.5, then we may get get identical allele counts, which then results in a covariance of 0, and correlation needs a division by cov, which would then produce nans...
# see: https://groups.google.com/u/2/g/plink2-users/c/iLu-8Ji_wic?pli=1   and  https://www.biostars.org/p/409611/
# $plink  --bfile $ibdDataLoc$'gwas1_new/ibd_all_hg19'  --keep-allele-order  --write-snplist --ld rs713041 rs732310  --out $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}_rs713041 --memory 2000
# another problem is that even SNPs which have valid MAF ranges, could have QC issues, so end up getting error " rs16936983 is monomorphic across all valid observations" and hence nans.
# $plink  --bfile $ibdDataLoc$'gwas1_new/ibd_all_hg19'  --keep-allele-order  --write-snplist --ld rs16936983 rs12256543  --out $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}_rs713041 --memory 2000
$plink  --bfile $ibdDataLoc$'gwas1_new/ibd_all_hg19'  --keep-allele-order  --write-snplist --maf 0.01 --max-maf 0.475 --r square  --extract $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.snplistOrig --freq  --out $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}2 --memory 2000

# the only way around this is to post-filter the .ld matrix with R and exclude all nan observations
arguments='/nfs/users/nfs_m/mk23/scripts/ldQC.R '$prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}'2.ld '$prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}'2.snplist '$prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# here I could add an additional check if the .ld matrix file still exists, because maybe  the above would have got rid of all remaining SNPs
if [ -s $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.snplist ]; then
echo "The file is NOT empty."



# --freq IGNORES the --maf threshold, and outputs the frq info for ALL snps, even those that fail the criteria
# therefore we would get incorrect snpinfo_ file...
# so we have to filter the frq file first to snplist, as those are the only ones that qualify
awk 'FNR == NR { file1[ $1 ] = $1;  next; }
FNR <= NR {  if( FNR == 1 || $2 in file1 ) {print $0} }
' $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.snplist $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}2.frq > $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.frq2

# add to snp info file
awk 'FNR == NR { file1[ $2 ] = $5;  next; }
FNR <= NR {  if( $2 in file1 ) {print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"file1[$2]} }
' $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.frq2 $ibdDataLoc$'gwas1_new/ibd_all_hg19.bim' >> $prscsrefs$'ibd/ldblk/snpinfo_'$LABEL$'_hm3'

wc -l < $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.snplist >> $prscsrefs$'ibd/blk_size'
else  # else file not empty after post ld calc QC
echo "The file is empty."
# if file is empty we want to add a 0 into the blk size file, otherwise we would not have an entry, which would then misalign it with the  blk_chr
echo -e 0 >> $prscsrefs$'ibd/blk_size'
cp $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}2.snplist $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.snplist
fi # file not empty after post ld calc QC

else # else file empty after initial SNP filtering
echo "The file is empty."
# if file is empty we want to add a 0 into the blk size file, otherwise we would not have an entry, which would then misalign it with the  blk_chr
echo -e 0 >> $prscsrefs$'ibd/blk_size'
cp $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.snplistOrig $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.snplist
fi # file not empty after initial SNP filtering

# also write out the chromosome number and the block size into a file, as that will be used by the next step
echo $chrnum >> $prscsrefs$'ibd/blk_chr'

#rm -rf $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.nosex
#rm -rf $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.log
#rm -rf $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.frq2
#rm -rf $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}2.frq
#rm -rf $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}2.ld
#rm -rf $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}2.snplist
#rm -rf $prscsrefs$'ibd/ldblk/'ldblk_${blk}_${LABEL}.frq2
fi # not header
done <$prscsrefs$'fourier_ls-all.bed'

# we now have to convert it to hdf5  # (NOTE: this script does not use the '$LABEL' argument, as it just defaults to "1kg", which is what I have used)
python3 $prscsrefs$'write_ldblk.py' $prscsrefs$'ibd/' $LABEL

# move the SNP info file to where the rest of the h5 files are
cp $prscsrefs$'ibd/ldblk/snpinfo_'$LABEL$'_hm3' $prscsrefs$'ibd/ldblk_'$LABEL'_chr/snpinfo_'$LABEL'_hm3'

# rename  folder as prs-cs does NOT accept folders with names like "ldblk_1kg_chr" it has to END with "ldblk_1kg" (otherwise get error: UnboundLocalError: local variable 'ref_dict' referenced before assignment )
mv $prscsrefs$'ibd/ldblk_'$LABEL'_chr' $prscsrefs$'ibd/ldblk_'$LABEL



# 3) run MTAG/SHAPRS
RUN_MTAG $newScratchLoc$'gwas3_cd/CD_MTAG_het' $newScratchLoc$'cd_sumstats_hm3_het' $newScratchLoc$'uc_sumstats_hm3_het' 1
RUN_MTAG $newScratchLoc$'gwas3_uc/UC_MTAG_het' $newScratchLoc$'uc_sumstats_hm3_het' $newScratchLoc$'cd_sumstats_hm3_het' 1

shaPRS_fixed $newScratchLoc$'gwas3_cd/CD_SHAPRS_het' $newScratchLoc$'cd_sumstats_hm3_het' $newScratchLoc$'uc_sumstats_hm3_het' $rho 1
shaPRS_fixed $newScratchLoc$'gwas3_uc/UC_SHAPRS_het' $newScratchLoc$'uc_sumstats_hm3_het' $newScratchLoc$'cd_sumstats_hm3_het' $rho 1


# 5) PRS-CS (wait until cluster job finishes)
PRSCS_allChroms_custom $newScratchLoc$'gwas3_cd/CD_MTAG_het_mtag' "gwas3_cd_mtag_het" $newScratchLoc $crossAncestryLoc $prscsrefs$'ibd/ldblk_'$LABEL
PRSCS_allChroms_custom $newScratchLoc$'gwas3_cd/CD_SHAPRS_het_shaprs' "gwas3_cd_shaprs_het" $newScratchLoc $crossAncestryLoc $prscsrefs$'ibd/ldblk_'$LABEL

PRSCS_allChroms_custom $newScratchLoc$'gwas3_uc/UC_MTAG_het_mtag' "gwas3_uc_mtag_het" $newScratchLoc $crossAncestryLoc $prscsrefs$'ibd/ldblk_'$LABEL
PRSCS_allChroms_custom $newScratchLoc$'gwas3_uc/UC_SHAPRS_het_shaprs' "gwas3_uc_shaprs_het" $newScratchLoc $crossAncestryLoc $prscsrefs$'ibd/ldblk_'$LABEL


# get the phenotypes for these
# 6) Build PRS
Concat_BuildPRS_custom 'gwas3_cd_mtag_het' $crossAncestryLoc $ibdDataLoc$'gwas1_new/ibd_all_hg19' $ibdDataLoc$'gwas1_new/all_found.fam' $newResultsLoc
Concat_BuildPRS_custom 'gwas3_cd_shaprs_het' $crossAncestryLoc $ibdDataLoc$'gwas1_new/ibd_all_hg19' $ibdDataLoc$'gwas1_new/all_found.fam' $newResultsLoc

Concat_BuildPRS_custom 'gwas3_uc_mtag_het' $crossAncestryLoc $ibdDataLoc$'gwas2_new/ibd_all_hg19' $ibdDataLoc$'gwas2_new/all_found.fam' $newResultsLoc
Concat_BuildPRS_custom 'gwas3_uc_shaprs_het' $crossAncestryLoc $ibdDataLoc$'gwas2_new/ibd_all_hg19' $ibdDataLoc$'gwas2_new/all_found.fam' $newResultsLoc

## 7) Evaluate
Evaluate_PRS $newResultsLoc$'PRSProfiles/gwas3_cd_mtag_het_PRSCS/gwas3_cd_mtag_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_cd_mtag_het' '0' # correlation_sq  0.0963 (sd: 0.000333 ) /  AUC  AUC 0.694 CI low: 0.6757  / high:  0.7116
Evaluate_PRS $newResultsLoc$'PRSProfiles/gwas3_cd_shaprs_het_PRSCS/gwas3_cd_shaprs_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_cd_shaprs_het' '0' # correlation_sq 0.107 (sd: 0.000399 ) /  AUC 0.706 CI low: 0.6886  / high:  0.7239

Evaluate_PRS $newResultsLoc$'PRSProfiles/gwas3_uc_mtag_het_PRSCS/gwas3_uc_mtag_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_uc_mtag_het' '0' # correlation_sq 0.0436 (sd: 0.000291 ) /  AUC 0.622 CI low: 0.6053  / high:  0.6379
Evaluate_PRS $newResultsLoc$'PRSProfiles/gwas3_uc_shaprs_het_PRSCS/gwas3_uc_shaprs_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_uc_shaprs_het' '0' # correlation_sq 0.0646 (sd: 0.000308 ) /   AUC 0.649 CI low: 0.6332  / high:  0.6651

$newScratchLoc$'gwas3_cd/CD_SHAPRS_het_shaprs'

##########################################
# II) run SMTPred, meta and baselines too:
###########################

# perform PRS-CS
PRSCS_allChroms_custom $newScratchLoc$'gwas3_cd/CD_SHAPRS_het_shaprsMerged_combined' "gwas3_cd_meta_het" $newScratchLoc $crossAncestryLoc $prscsrefs$'ibd/ldblk_'$LABEL
PRSCS_allChroms_custom $newScratchLoc$'gwas3_uc/UC_SHAPRS_het_shaprsMerged_combined' "gwas3_uc_meta_het" $newScratchLoc $crossAncestryLoc $prscsrefs$'ibd/ldblk_'$LABEL

PRSCS_allChroms_custom $newScratchLoc$'cd_sumstats_hm3_het' "gwas3_cd_primary_het" $newScratchLoc $crossAncestryLoc $prscsrefs$'ibd/ldblk_'$LABEL
PRSCS_allChroms_custom $newScratchLoc$'uc_sumstats_hm3_het' "gwas3_uc_primary_het" $newScratchLoc $crossAncestryLoc $prscsrefs$'ibd/ldblk_'$LABEL

# wait until jobs finish on cluster...
Concat_BuildPRS_custom 'gwas3_cd_meta_het' $crossAncestryLoc $ibdDataLoc$'gwas1_new/ibd_all_hg19' $ibdDataLoc$'gwas1_new/all_found.fam' $newResultsLoc
Concat_BuildPRS_custom 'gwas3_uc_meta_het' $crossAncestryLoc $ibdDataLoc$'gwas2_new/ibd_all_hg19' $ibdDataLoc$'gwas2_new/all_found.fam' $newResultsLoc

Concat_BuildPRS_custom 'gwas3_cd_primary_het' $crossAncestryLoc $ibdDataLoc$'gwas1_new/ibd_all_hg19' $ibdDataLoc$'gwas1_new/all_found.fam' $newResultsLoc
Concat_BuildPRS_custom 'gwas3_uc_primary_het' $crossAncestryLoc $ibdDataLoc$'gwas2_new/ibd_all_hg19' $ibdDataLoc$'gwas2_new/all_found.fam' $newResultsLoc

Evaluate_PRS $newResultsLoc$'PRSProfiles/gwas3_cd_meta_het_PRSCS/gwas3_cd_meta_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_cd_meta_het' '0' # correlation_sq 0.095 (sd: 0.000351 ) /   AUC 0.693 CI low: 0.6756  / high:  0.7112
Evaluate_PRS $newResultsLoc$'PRSProfiles/gwas3_uc_meta_het_PRSCS/gwas3_uc_meta_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_uc_meta_het' '0' # correlation_sq 0.0612 (sd: 0.000325 ) /    AUC 0.646 CI low: 0.63  / high:  0.662

Evaluate_PRS $newResultsLoc$'PRSProfiles/gwas3_cd_primary_het_PRSCS/gwas3_cd_primary_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_cd_primary_het' '0' # correlation_sq 0.103 (sd: 0.000344 ) /    AUC 0.701 CI low: 0.6834  / high:  0.7187
Evaluate_PRS $newResultsLoc$'PRSProfiles/gwas3_uc_primary_het_PRSCS/gwas3_uc_primary_het_PRSCS_all.sscore' $newResultsLoc$'gwas3_uc_primary_het' '0' # correlation_sq 0.052 (sd: 0.00028 ) /     AUC 0.632 CI low: 0.6155  / high:  0.6478


###############
# SMTPred also requires the 'PheB' PRS be built from the PheA PRS file
$crossAncestryLoc$'gwas3_cd_primary_het_PRSCS_no_dupes'
# here we build a CD PRS, from UC
BuildPRS_generic 'gwas3_cd_primary_het_pheB' $crossAncestryLoc$'gwas3_uc_primary_het_PRSCS_no_dupes' $ibdDataLoc$'gwas1_new/ibd_all_hg19' $ibdDataLoc$'gwas1_new/all_found.fam' $newResultsLoc 1

# here we build a UC PRS, from CD
BuildPRS_generic 'gwas3_uc_primary_het_pheB' $crossAncestryLoc$'gwas3_cd_primary_het_PRSCS_no_dupes' $ibdDataLoc$'gwas2_new/ibd_all_hg19' $ibdDataLoc$'gwas2_new/all_found.fam' $newResultsLoc 1

perform_SMPTpred_PRSCS $newResultsLoc$'PRSProfiles/gwas3_cd_het_SMTPRED/' $newResultsLoc$'PRSProfiles/gwas3_cd_primary_het_PRSCS/gwas3_cd_primary_het_PRSCS_all.sscore' $newResultsLoc$'PRSProfiles/gwas3_cd_primary_het_pheB/gwas3_cd_primary_het_pheB_all.sscore' $newScratchLoc$'cd_sumstats_hm3_het' $newScratchLoc$'uc_sumstats_hm3_het'
perform_SMPTpred_PRSCS $newResultsLoc$'PRSProfiles/gwas3_uc_het_SMTPRED/' $newResultsLoc$'PRSProfiles/gwas3_uc_primary_het_PRSCS/gwas3_uc_primary_het_PRSCS_all.sscore' $newResultsLoc$'PRSProfiles/gwas3_uc_primary_het_pheB/gwas3_uc_primary_het_pheB_all.sscore' $newScratchLoc$'uc_sumstats_hm3_het' $newScratchLoc$'cd_sumstats_hm3_het'

Evaluate_PRS $newResultsLoc$'PRSProfiles/gwas3_cd_het_SMTPRED/_SMTPRED.profile' $newResultsLoc$'gwas3_smtpred_cd_het' '0' #  correlation_sq 0.1 (sd: 0.000318 ) / AUC 0.702 CI low: 0.6844  / high:  0.7196 
Evaluate_PRS $newResultsLoc$'PRSProfiles/gwas3_uc_het_SMTPRED/_SMTPRED.profile' $newResultsLoc$'gwas3_smtpred_uc_het' '0' #   correlation_sq 0.0586 (sd: 0.000323 ) / AUC 0.642 CI low: 0.626  / high:  0.6581

# num SNPs for UC: 856877

############################

#CD
               r2      AUC
CD_primary      0.103   0.701
CD_Meta         0.095   0.693
CD_SMTPred      0.100   0.702
CD_MTAG         0.096   0.694
CD_shaPRS       0.107   0.706

# SD            SD        LOW      HIGH
CD_primary      0.000344  0.6834   0.7187
CD_Meta         0.000351  0.6756   0.7112
CD_SMTPred      0.000318  0.6844   0.7196     
CD_MTAG         0.000333  0.6757   0.7116
CD_shaPRS       0.000399  0.6886   0.7239


#UC
                r2      AUC
UC_primary      0.052   0.632
UC_Meta         0.061   0.646
UC_SMTPred      0.059   0.642
UC_MTAG         0.044   0.622
UC_shaPRS       0.065   0.649

# SD            SD        LOW      HIGH
UC_primary      0.00028   0.6155   0.6478 
UC_Meta         0.000325  0.63     0.662
UC_SMTPred      0.000323  0.626    0.6581    
UC_MTAG         0.000291  0.6053   0.6379
UC_shaPRS       0.000308  0.6332   0.6651



# need to emphasise that MTAG's main limitation is that learns a single weight applied to all SNPs. thus where power allows as dataset sizes grow, shaPRS is expected to outperform MTAG
                 
##############################################


##### Reusable  functions


# converts my sumstats format into the PRS-CS format, postfixing input with '_PRSCS' (specify if binary)
function convertToPRSCS { 
sumstatsLoc=$1
isBinary=$2

# convert my _sumstats format:
# chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N
# to PRS-CS:
# SNP          A1   A2   BETA/OR      P

awk -v isBinary="$isBinary" ' 
FNR <= NR {
if (FNR == 1) {
coefName="BETA"
if(isBinary == "1") {coefName="OR"}
print "SNP\tA1\tA2\t"coefName"\tP" 
}
else { 
COEF=$7
if(isBinary == "1") {COEF=exp($7)}
print $3"\t"$4"\t"$5"\t"COEF"\t"$9 
}
}
' OFS='\t' $sumstatsLoc > $sumstatsLoc$'_PRSCS'

#head $sumstatsLoc$'_PRSCS'

}
export -f convertToPRSCS # this makes local functions executable when bsubbed



# generates manhattan plot from a Susmstats formatted file
function produce_manhattan_ss { 
ssFile=$1 # this is the standard plink .phe.glm.logistic.hybrid file
output_loc=$2
plottitle=$3
sigTreshold=$4

#   1    2     3    4    5       6           7   8   9   10 
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N
# TO
# extract summary stats into format into the following signature: # SNP CHR BP         P    zscore
awk '{ 
if ( FNR == 1) {print "SNP\tCHR\tBP\tP\tzscore"}
else  { print $3"\t"$1"\t"$2"\t"$9"\t"$7/$8 }
 }' $ssFile > $output_loc$'_gwasResults'


# call Rscript to produce plot
arguments='/nfs/users/nfs_m/mk23/scripts/make_manhattan.R '$output_loc$'_gwasResults '$output_loc$' '$plottitle$' '$sigTreshold
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments 

#rm -rf $output_loc$'_gwasResults'
}
export -f produce_manhattan_ss # this makes local functions executable when bsubbed




# Reusable SMTPred function
function perform_SMPTpred_PRSCS {
outputDir=$1
pheA_PRSprofile=$2
pheB_PRSprofile=$3 # the 'pheB' .profile file is for the SAME indis as PheA, but trained on 'pheB'
pheA_sumstats=$4
pheB_sumstats=$5

plink='/nfs/users/nfs_m/mk23/software/plink/plink'
smtpred="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/smtpred/smtpred/"
ldscDir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/smtpred/ldsc/"
w_hm3_snplist="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/ldsc/w_hm3.snplist"
ldscRefpanelLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/ldsc/eur_w_ld_chr/"

# if I supply the SMTPred with the same filename in different folders eg pheA/sumstats and pheB/sumstats, then the idiot will overwrite the two
cp $pheB_sumstats $pheB_sumstats$'_pheB'

mkdir -p $outputDir

# the filenames for the profile score files must match to the sumstats files
filename_pheA=$(awk '{idx = split(FILENAME, parts, "/"); print parts[idx]; nextfile}' $pheA_sumstats)
filename_pheB=$(awk '{idx = split(FILENAME, parts, "/"); print parts[idx]; nextfile}' $pheB_sumstats$'_pheB')

# as we use PLINK2 scoring now, the file format is slightly different to PLINK1, but smtpred expects PLINK1 format:
awk '{if (FNR ==1) {print "                    FID                     IID  PHENO    CNT   CNT2 SCORESUM"} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' OFS=" " $pheA_PRSprofile > $outputDir$'/'$filename_pheA$'.profile'
awk '{if (FNR ==1) {print "                    FID                     IID  PHENO    CNT   CNT2 SCORESUM"} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' OFS=" " $pheB_PRSprofile > $outputDir$'/'$filename_pheB$'.profile'


# obtain h2 and rG estimates for Pheno A and B via LDSC ( uses the original GWAS assoc stats)
/usr/bin/python2  $smtpred$'ldsc_wrapper.py' \
    --out $outputDir \
    --files $pheA_sumstats \
            $pheB_sumstats$'_pheB' \
    --ldscpath $ldscDir \
    --snplist $w_hm3_snplist \
    --ref_ld $ldscRefpanelLoc \
    --w_ld $ldscRefpanelLoc

# generate the new .PRS, from the 2 supplied LDPred2
/usr/bin/python2 $smtpred$'smtpred.py' \
  --h2file $outputDir$'/ldsc_h2s.txt' \
  --rgfile $outputDir$'/ldsc_rgs.txt' \
  --nfile $outputDir$'/ldsc_ns.txt' \
  --scorefiles $outputDir$'/'$filename_pheA$'.profile' \
               $outputDir$'/'$filename_pheB$'.profile' \
  --out $outputDir$'/'
  
# this writes a final per indi PRS score file:
# rest of the pipeline works with PLINK2 profile scores, which are 6 columns, but this is PLINK1 style score file with just 3
awk 'FNR == NR { file1[ $1 ] = $3; next; } FNR <= NR { if (FNR == 1) {print "FID\tIID\tPHENO1\tX\tX\tSCORE1_SUM"} else  { if ( $2 in file1 ) { print $1"\t"$2"\t"file1[$2]"\tX\tX\t"$3  } } }' $pheA_PRSprofile $outputDir$'/multi_trait.score' > $outputDir$'_SMTPRED.profile'

}
export -f perform_SMPTpred_PRSCS # this makes local functions executable when bsubbed	



# outPutLocation=$scratchLocV4$'rhotest'
# pheASumstats=$rawLocV4$"AAAGen_Leicester_UKBB_shaprs"
# pheBSumstats=$rawLocV4$'AAD_RELATED_hm3'
# pheASumstats=$rawLocV4$'AAD_RELATED_hm3'
# calculates empriical overlap between studies: load the rho into session via: source $outPutLocation$'rho'
function GetEmpiricalRho {
outPutLocation=$1
pheASumstats=$2
pheBSumstats=$3

# perform meta analysis
arguments='/nfs/users/nfs_m/mk23/scripts/metaAnalysis.R '$pheASumstats$' '$pheBSumstats$' '$outPutLocation$'_meta /nfs/users/nfs_m/mk23/scripts/shaPRS.R'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# in GWAS3 IBD naive GWAS, find SNPs with p > 0.01 # this is to find the empirical correlation, we find the 'mildly to non associated' SNPS
# otherwise we would get 
awk '{ if ( $9 > 0.01 && NR >1) {print $3} }' $outPutLocation$'_meta' > $outPutLocation$'_unassociated'
wc -l $outPutLocation$'_unassociated'
# 6,933,546

# Get same SNPs in the two subphenos
# create a file with signature RSiD, Beta_CD, SE_CD, Beta_UC, SE_UC
awk 'FNR == NR { file1[ $1 ] = $1; next; } 
FNR <= NR { { if ( $3 in file1) { print $3"\t"$7"\t"$8 } }  }
'  $outPutLocation$'_unassociated' $pheASumstats > $outPutLocation$'A_unassociated'

awk 'FNR == NR { file1[ $1 ] = $0; next; } 
FNR <= NR { { 
if (FNR == 1) {print "RSid\tBeta_CD\tSE_CD\tBeta_UC\tSE_UC" }
if ( $3 in file1) { print file1[$3]"\t"$7"\t"$8 } }  
}
'  $outPutLocation$'A_unassociated' $pheBSumstats > $outPutLocation$'A_B_unassociated'


# run R script to get rho (correlation between studies,  this is a scalar)
# quantifies the linear association between estimated coefficients, for NON-disease SNPs
# rho = cor(Beta_CD/SE_CD, Beta_UC/SE_UC)
arguments='/nfs/users/nfs_m/mk23/scripts/rho.R '$outPutLocation$'A_B_unassociated '$outPutLocation$'rho'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# 
#source $outPutLocation$'rho'
}
export -f GetEmpiricalRho # this makes local functions executable when bsubbed


# allScoreLoc=$newResultsLoc$'PRSProfiles/gwas3_cd_mtag_PRSCS/gwas3_cd_mtag_PRSCS_all.sscore' 
# outResultLoc=$newResultsLoc$'gwas3_cd_mtag' 
# doAppend='0'
# evaluates the r^2 and AUC of a PRS on a test set, and writes results to $evaluates$"_r2_AUC" and $evaluates$"_r2"
function Evaluate_PRS { 
allScoreLoc=$1
outResultLoc=$2
doAppend=$3

# we are working with PLINK2 profile scores, which are 6 columns, but the R scripts expects a 2 col one now
awk '{print $2"\t"$3"\t"$6 }' $allScoreLoc > $allScoreLoc$'_3col'


arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$allScoreLoc$'_3col '$outResultLoc$'_r2 1 1 '$doAppend$' 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments 
}
export -f Evaluate_PRS # this makes local functions executable when bsubbed


# phe_name='gwas3_cd_shaPRSLQmtagHQ'
# PRS_rawLoc=$crossAncestryResults$'gwas3_cd_shaPRSLQmtagHQ'
# bfile=$ibdDataLoc$'gwas1/ibd_all'
# famPheno=$ibdCommonDataScratch$'gwas1/all_found.fam'
# outLoca=$newResultsLoc
# PRSCS='1'
# builds a PRS from PRS-CS or LDPred2 formatted score file
function BuildPRS_generic {
phe_name=$1
PRS_rawLoc=$2
bfile=$3
famPheno=$4
outLoca=$5
PRSCS=$6

# generate PRS per each chrom
mkdir -p $outLoca$'PRSProfiles/'
mkdir -p $outLoca$'PRSProfiles/'$phe_name$'/'

# only submit if it doesn't exist yet
outfile=$outLoca$'PRSProfiles/'$phe_name$'/'$phe_name$'_all.sscore'
rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
b=1



#PRS-CS
if [ $PRSCS == "1" ] ; then
echo "PRS-CS PRS"
scoreText='2 4 6'
else
# LDpred2/Rapido
echo "LDpred2 PRS"
scoreText=''
fi 
arguments=' --memory '$plinkMem$' --bfile '$bfile$' --fam '$famPheno' --threads 16 --score '$PRS_rawLoc$' '$scoreText$' ignore-dup-ids cols=-scoreavgs,+scoresums --out '$outLoca$'PRSProfiles/'$phe_name$'/'$phe_name$'_all'
echo $arguments
$plink2 $arguments

fi # if score exists

}
export -f BuildPRS_generic # this makes local functions executable when bsubbed



   
# i=20
# penoNam='gwas3_cd_mtag'
# scratchLoca=$crossAncestryResults
# bfile=$ibdDataLoc$'gwas1/ibd_all'
# famPheno=$ibdCommonDataScratch$'gwas1/all_found.fam'
# outLoca=$newResultsLoc

function Concat_BuildPRS {
penoNam=$1
scratchLoca=$2
bfile=$3
famPheno=$4
outLoca=$5


# cp $scratchLoca$penoNam$'_PRSCS' $scratchLoca$penoNam$'_PRSCS_old'

#echo $scratchLoca
#ls $scratchLoca$penoNam$'_EUR_pst_eff_a1_b0.5_phiauto_chr'*'.txt'

#  concat all chrom PRS into single files
cat $scratchLoca$penoNam$'_EUR_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $scratchLoca$penoNam$'_PRSCS'
head $scratchLoca$penoNam$'_PRSCS'
wc -l $scratchLoca$penoNam$'_PRSCS'

# remove dupes:
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $scratchLoca$penoNam$'_PRSCS' > $scratchLoca$penoNam$'_PRSCS_no_dupes'

PRS_rawLoc=$scratchLoca$penoNam$'_PRSCS_no_dupes'
phe_name=$penoNam$'_PRSCS'
PRSCS='1'


BuildPRS_generic $phe_name $PRS_rawLoc $bfile $famPheno $outLoca $PRSCS
 
}
export -f Concat_BuildPRS # this makes local functions executable when bsubbed


function Concat_BuildPRS_custom {
penoNam=$1
scratchLoca=$2
bfile=$3
famPheno=$4
outLoca=$5


# cp $scratchLoca$penoNam$'_PRSCS' $scratchLoca$penoNam$'_PRSCS_old'

#echo $scratchLoca
#ls $scratchLoca$penoNam$'_EUR_pst_eff_a1_b0.5_phiauto_chr'*'.txt'

#  concat all chrom PRS into single files
cat $scratchLoca$penoNam$'/_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $scratchLoca$penoNam$'_PRSCS'
head $scratchLoca$penoNam$'_PRSCS'
wc -l $scratchLoca$penoNam$'_PRSCS'

# remove dupes:
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $scratchLoca$penoNam$'_PRSCS' > $scratchLoca$penoNam$'_PRSCS_no_dupes'

PRS_rawLoc=$scratchLoca$penoNam$'_PRSCS_no_dupes'
phe_name=$penoNam$'_PRSCS'
PRSCS='1'


BuildPRS_generic $phe_name $PRS_rawLoc $bfile $famPheno $outLoca $PRSCS
 
}
export -f Concat_BuildPRS_custom # this makes local functions executable when bsubbed




#outputDir=$newResultsLoc$'gwas3_cd/bootstrap_nodupes_amb/CD_MTAG'$i
#pheA_sumstats=$newResultsLoc$'gwas3_cd/bootstrap_nodupes/gwas3_cd_raw'$i
#pheB_sumstats=$newResultsLoc$'gwas3_uc/bootstrap_nodupes/gwas3_uc_raw'$i
# maps my sumstats format to the format needed for MTAG, and performs MTAG, output is  $outputDir$"_mtag" (# --use_beta_se does not work in MTAG for now...)
function RUN_MTAG { 
outputDir=$1
pheA_sumstats=$2
pheB_sumstats=$3
includeAmbigSNPs=$4

if [ -n "$4" ] ; then
includeAmbigSNPs=' --incld_ambig_snps '
else 
includeAmbigSNPs=''
fi
echo "includeAmbigSNPs is: "$includeAmbigSNPs

# map my sumstats format:
# 1       2        3      4      5            6          7        8      9       10
#chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N
#10      100146895       rs2296438       T       C       0.6544  0.0024  0.0082  0.7694  964057
# to MTAG
#   1             2        3              4       5        6             7           8           9             
# snpid           chr     bpos            a1      a2      freq           z          pval          n
# make sure to avoid division by zero, otherwise awk may truncate file and thus create incorrect result
awk '{
if (FNR == 1) {print "snpid\tchr\tbpos\ta1\ta2\tfreq\tz\tpval\tn" }
is_7_numeric = $7 + 0 == $7
is_8_numeric = $8 + 0 == $8
if(is_7_numeric && is_8_numeric) {
print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7/$8"\t"$9"\t"$10} 
}'  $pheA_sumstats > $pheA_sumstats$'_MTAG'


awk '{
if (FNR == 1) {print "snpid\tchr\tbpos\ta1\ta2\tfreq\tz\tpval\tn" }
is_7_numeric = $7 + 0 == $7
is_8_numeric = $8 + 0 == $8
if(is_7_numeric && is_8_numeric) {
print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7/$8"\t"$9"\t"$10} 
}'  $pheB_sumstats > $pheB_sumstats$'_MTAG'

#head $pheA_sumstats$'_MTAG'
/usr/bin/python '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/revision/mtag/mtag/mtag.py' --sumstats $pheA_sumstats$'_MTAG',$pheB_sumstats$'_MTAG'  --out $outputDir --n_min 0.0 --force $includeAmbigSNPs > $outputDir$'.log'



# 2022/05/08/07:29:07 PM Dropped 88502 SNPs due to strand ambiguity, 1035282 SNPs remain in intersection after merging trait1

# convert the MTAG output into my own format 
# MTAG: 
#  1               2       3       4      5       6               7        8          9                     10                       11                     12
# SNP             CHR     BP      A1      A2      Z               N       FRQ     mtag_beta               mtag_se                 mtag_z                  mtag_pval
# rs4040617       1       779322  G       A       0.314132        10874   0.1134  0.0059077854624592615   0.030391777030177466    0.19438762848888813     0.8458723762003009

# map my sumstats format:
# 1       2        3      4      5            6          7        8      9       10
#chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N

awk '{if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN"}
else{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$8"\t"$9"\t"$10"\t"$12"\t"$7} }' $outputDir$"_trait_1.txt" > $outputDir$"_mtagMerged"

# in case there were SNPs that were not present in trait 2 these would get removed by the above, we use trait 1's values
awk 'FNR == NR { file1[ $3 ] = $0; next;  }
FNR <= NR { if ($3 in file1) {print file1[$3] } else {print $0} }
' $outputDir$'_mtagMerged' $pheA_sumstats > $outputDir$'_mtag'

#wc -l $outputDir$'_mtagMerged'
#wc -l $outputDir$'_mtag'

}
export -f RUN_MTAG # this makes local functions executable when bsubbed




#  Performs the shaPRS step that produces a _sumstats formatted file for both the blended and the combined sumstats: output $outPutLocation'_shaprs'
function shaPRS_fixed { 
outPutLocation=$1
pheASumstats=$2
pheBSumstats=$3
rhoSet=$4


if [ -n "$5" ] ; then
includeAmbigSNPs='1'
else 
includeAmbigSNPs='0'
fi
echo "includeAmbigSNPs is: "$includeAmbigSNPs

# Create input file for adjusting script:
#       SNP	CHR	BP	Beta_A	SE_A	Beta_B	SE_B
# rs4040617   1  779322 -0.0017630 0.008608 -0.010990 0.008592
awk 'FNR == NR { file1[ $3 ] = $7"\t"$8"\t"$4"\t"$5; next; } 
FNR <= NR { { 
if (FNR == 1) {print "SNP\tCHR\tBP\tBeta_A\tSE_A\tA1.x\tA2.x\tBeta_B\tSE_B\tA1.y\tA2.y" } else {
if ( $3 in file1) { print $3"\t"$1"\t"$2"\t"file1[$3]"\t"$7"\t"$8"\t"$4"\t"$5} }   }
}
' $pheBSumstats $pheASumstats > $outPutLocation$'_SE_meta'
# this still has the B/A mixed up

# Export out the lFDR values
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_adjust_wrapper.R '$outPutLocation$'_SE_meta '$outPutLocation$'_lFDR_meta '$rhoSet$' /nfs/users/nfs_m/mk23/scripts/shaPRS.R '$includeAmbigSNPs
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments



# R script to blend between the two association stats
arguments='/nfs/users/nfs_m/mk23/scripts/sumstatsBlender_shaPRS_meta_wrapper.R '$pheASumstats$' '$pheBSumstats$' '$outPutLocation$'_lFDR_meta_SNP_lFDR '$rhoSet$' '$outPutLocation$'_shaprsMerged /nfs/users/nfs_m/mk23/scripts/shaPRS.R '$includeAmbigSNPs
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# in case there were SNPs that were not present in trait 2 these would get removed by the above, we use trait 1's values
awk 'FNR == NR { file1[ $3 ] = $0; next;  }
FNR <= NR { if ($3 in file1) {print file1[$3] } else {print $0} }
' $outPutLocation$'_shaprsMerged' $pheASumstats > $outPutLocation$'_shaprs'


mv $outPutLocation$'_shaprsMerged_combined' $outPutLocation$'_shaprsMerged_combined2'

awk 'FNR == NR { file1[ $3 ] = $0; next;  }
FNR <= NR { if ($3 in file1) {print file1[$3] } else {print $0} }
' $outPutLocation$'_shaprsMerged_combined2' $pheASumstats > $outPutLocation$'_shaprsMerged_combined'



# find out the max Q-val of the run
awk 'BEGIN {max = -1; maxSNP="";} FNR <= NR { { if ( FNR > 1 && max < $3) { max = $3; maxSNP=$1 } } } END {print "SNP with max Qval: "maxSNP"\t"max;}' FS=" " $outPutLocation$'_lFDR_meta_SNP_lFDR'

}
export -f shaPRS_fixed # this makes local functions executable when bsubbed 



rawPRSLoc=$newResultsLoc$'gwas3_uc/UC_SHAPRS_shaprs'
pheno_name="gwas3_uc_shaprs"
scratchLoca=$newScratchLoc
resultsDir=$crossAncestryLoc
  

  
    
  
rawPRSLoc=$newScratchLoc$'gwas3_cd/CD_SHAPRS_het_shaprs'
pheno_name="gwas3_cd_shaprs_het"
scratchLoca=$newScratchLoc
resultsDir=$crossAncestryLoc
refLoc=$prscsrefs$'ibd/ldblk_'$LABEL # must 
k=21

function PRSCS_allChroms { 
rawPRSLoc=$1
pheno_name=$2
scratchLoca=$3
resultsDir=$4

convertToPRSCS $rawPRSLoc '1'
mkdir -p $scratchLoca$'logs/'

# get AVERAGE sample size, this was better than minimum
Neur=$(awk '{ if(FNR > 1) {sum += $10; n++ } } END { printf("%0.f", sum / n) ; }' $rawPRSLoc)
echo $Neur


validBim=$rawPRSLoc$'_validBim'
# create a 'fake' .bim file for PRSx from my sumstats file that has all the info
awk '{if (FNR !=1) {print $1"\t"$3"\t0\t"$2"\t"$4"\t"$5} }' $rawPRSLoc > $validBim$'.bim'

k=21
mkdir -p $scratchLoca$pheno_name$'/'
for ((k=1; k<=22; k++)); do # $numChroms


#resultsDir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/"
# 1) PRS-CS
# check if output doesn't already exist
outfile=$resultsDir$pheno_name$'_EUR_pst_eff_a1_b0.5_phiauto_chr'$k$'.txt' # $scratchLoca$pheno_name$'/'$pheno_name$'_EUR_pst_eff_a1_b0.5_phiauto_chr'$k$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
b=1
else
b=1
echo "PRS-CS SUBMITTING: "$pheno_name$' '$k


prscsmem=43000
bsub -q normal -m "modern_hardware" -G team152 -n16 -J ${pheno_name}_${i} -R "span[hosts=1] select[mem>${prscsmem}] rusage[mem=${prscsmem}]" -M${prscsmem} -o $scratchLoca$'logs/'$pheno_name$'_'$k$'.out' -e $scratchLoca$'logs/'$pheno_name$'_'$k$'.err'  "performPRSCS $rawPRSLoc$'_PRSCS' $Neur $validBim $pheno_name $k 16 $resultsDir"


fi

done # end of chrom loop

} 
export -f PRSCS_allChroms # this makes local functions executable when bsubbed



eursums=$rawPRSLoc$'_PRSCS'
Neur=$Neur
validBim=$validBim
mypheno=$pheno_name
chrom=$k
ncores=16
out_dir=$resultsDir

      

# Performs PRS-CS for 1 population (EUR)
function performPRSCS { 
eursums=$1
Neur=$2
validBim=$3
mypheno=$4
chrom=$5
ncores=$6
out_dir=$7

if [ -n "$7" ] ; then
b=1
else 
out_dir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/"
fi
echo "out dir is: "$out_dir


export MKL_NUM_THREADS=$ncores
export NUMEXPR_NUM_THREADS=$ncores
export OMP_NUM_THREADS=$ncores

prscsxdir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsscript/PRScsx/"
prscsrefs="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsrefs/"

python $prscsxdir/PRScsx.py --ref_dir=$prscsrefs --bim_prefix=$validBim --sst_file=$eursums --n_gwas=$Neur --pop=EUR --out_dir=$out_dir --out_name=$mypheno --chrom=$chrom  --seed=42

}
export -f performPRSCS # this makes local functions executable when bsubbed



function PRSCS_allChroms_custom { 
rawPRSLoc=$1
pheno_name=$2
scratchLoca=$3
resultsDir=$4
refDir=$5

convertToPRSCS $rawPRSLoc '1'
mkdir -p $scratchLoca$'logs/'

# get AVERAGE sample size, this was better than minimum
Neur=$(awk '{ if(FNR > 1) {sum += $10; n++ } } END { printf("%0.f", sum / n) ; }' $rawPRSLoc)
echo $Neur


validBim=$rawPRSLoc$'_validBim'
# create a 'fake' .bim file for PRSx from my sumstats file that has all the info
awk '{if (FNR !=1) {print $1"\t"$3"\t0\t"$2"\t"$4"\t"$5} }' $rawPRSLoc > $validBim$'.bim'

k=21
mkdir -p $scratchLoca$pheno_name$'/'
for ((k=1; k<=22; k++)); do # $numChroms


#resultsDir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/"
# 1) PRS-CS
# check if output doesn't already exist
outfile=$resultsDir$pheno_name$'/_pst_eff_a1_b0.5_phiauto_chr'$k$'.txt' # $scratchLoca$pheno_name$'/'$pheno_name$'_EUR_pst_eff_a1_b0.5_phiauto_chr'$k$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
b=1
else
b=1
echo "PRS-CS SUBMITTING: "$pheno_name$' '$k


prscsmem=43000
bsub -q normal -m "modern_hardware" -G team152 -n16 -J ${pheno_name}_${i} -R "span[hosts=1] select[mem>${prscsmem}] rusage[mem=${prscsmem}]" -M${prscsmem} -o $scratchLoca$'logs/'$pheno_name$'_'$k$'.out' -e $scratchLoca$'logs/'$pheno_name$'_'$k$'.err'  "performPRSCS_custom $rawPRSLoc$'_PRSCS' $Neur $validBim $pheno_name $k 16 $resultsDir $refDir"


fi

done # end of chrom loop

} 
export -f PRSCS_allChroms_custom # this makes local functions executable when bsubbed




function performPRSCS_custom { 
sumsFile=$1
NGWAS=$2
validBim=$3
mypheno=$4
chrom=$5
ncores=$6
out_dir=$7
refLoc=$8


export MKL_NUM_THREADS=$ncores
export NUMEXPR_NUM_THREADS=$ncores
export OMP_NUM_THREADS=$ncores

if [ -n "$7" ] ; then
b=1
else 
out_dir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/"
fi
echo "out dir is: "$out_dir

prscssorig="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscss/"
mkdir -p $out_dir$mypheno
python $prscssorig/PRScs/PRScs.py --ref_dir=$refLoc --bim_prefix=$validBim --sst_file=$sumsFile --n_gwas=$NGWAS --out_dir=$out_dir$mypheno/ --chrom=$chrom  --seed=42

}
export -f performPRSCS_custom # this makes local functions executable when bsubbed

chrom=19
mkdir -p $out_dir$mypheno
python $prscssorig/PRScs/PRScs.py --ref_dir=$refLoc --bim_prefix=$validBim --sst_file=$sumsFile --n_gwas=$NGWAS --out_dir=$out_dir$mypheno/ --chrom=$chrom  --seed=42


###########################

shaPRS_new $pheA_sumstats $pheB_sumstats $shaPRSPre$pheno$'_'$i $rho '0'

shaPRS_fixed $newScratchLoc$'gwas3_cd/CD_SHAPRS_het' $newScratchLoc$'cd_sumstats_hm3_het' $newScratchLoc$'uc_sumstats_hm3_het' $rho 1

$newScratchLoc$'gwas3_cd/CD_SHAPRS_het_SE_meta'

# Plot
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_qval_manhattan.R '$newScratchLoc$'gwas3_cd/CD_SHAPRS_het_SE_meta '$newScratchLoc$'gwas3_cd/CD_SHAPRS_het_lFDR_meta_SNP_lFDR '$newResultsLoc$'manhattan_lFDR CD '$newScratchLoc$'gwas3_cd/CD_SHAPRS_het_shaprsMerged_combined'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
# "number of SNPs with lFDR < 0.01: 2794"



$newScratchLoc$'gwas3_cd/CD_SHAPRS_het_sumstats_meta'
$newScratchLoc$'gwas3_cd/CD_SHAPRS_het_sumstats_meta'

arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_qval_manhattan.R '$shaPRSPre$pheno$'_'$i$'_SE_meta '$shaPRSPre$pheno$'_'$i$'_lFDR_meta_SNP_lFDR '$shaPRSPre$pheno$'_'$i$' '$pheno$' '$shaPRSPre$pheno$'_'$i$'_sumstats_meta'

j=5
ibdarray_target_extended=( 'gwas1/' 'gwas2/' 'gwas3/' 'gwas3_cd/' 'gwas3_uc/' ) 
ibdLoc=$baseLoc$'ibd/'
ibdDataLoc=$ibdLoc$'data/'
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
shaPRSPre=$bootStrapLoc$'shaPRSPre/'


head $shaPRSPre$'UC_1_sumstats_meta'
head $shaPRSPre$'UC_1_SE_meta'

$shaPRSPre$pheno$'_'$i$'_SE_meta'
######################

# Sample size calculations:

#Training: (data from Laura)
humancoreexome
EUR
IBD:10888/10308
CD: 5,400 /10,308
UC: 4,647 /10,308

# num overlapping controls: 10308


$newResultsLoc$'gwas3_cd/gwas3_cd_hm3'
$newResultsLoc$'gwas3_uc/gwas3_uc_hm3'
# CD: 14,174
# UC: 12,812


# Test sets:
$ibdDataLoc$'gwas1_new/ibd_all_hg19.fam'
$ibdDataLoc$'gwas1_new/ibd_all_hg19.fam'

#CD 
awk '{count[$6]++} END {for (word in count) print word, count[word]}' $ibdDataLoc$'gwas1_new/all_found.fam'
1 2,896
2 1,181

# UC 
awk '{count[$6]++} END {for (word in count) print word, count[word]}' $ibdDataLoc$'gwas2_new/all_found.fam'
1 2,764
2 1,909

awk '{count[$6]++} END {for (word in count) print word, count[word]}' $ibdCommonDataScratch$'gwas1/all_found.fam'

