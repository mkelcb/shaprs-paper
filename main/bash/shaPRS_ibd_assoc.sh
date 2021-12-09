####################################################
#  IBD Bootstrap generation and GWAS script

####################################################
#bcftools='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/software/bcftools/./bcftools'
ibd_raw_merged='/lustre/scratch115/realdata/mdt0/projects/ibdgwas/new_imputation/merged/all_studies.filtered.sampleids_cleaned_lowercase.vcf.gz'

ibd_raw_cohortbaseLoc='/lustre/scratch115/realdata/mdt0/projects/ibdgwas/new_imputation/'


ibdLoc=$baseLoc$'ibd/'
ibdDataLoc=$ibdLoc$'data/'
ibdCommonDataLoc=$ibdDataLoc$'0_common/'
ibdCommonDataScratch=$ibdCommonDataLoc$'scratch/'
ibd_all_cleaned=$ibdCommonDataLoc$'ibd_all_cleaned' # stores all 5 cohorts merged and cleaned/filtered
ibd_raw_cohortbaseLoc='/lustre/scratch115/realdata/mdt0/projects/ibdgwas/new_imputation/'

covariates_filteredIBD=$ibdDataLoc$"covariates"
ibdResultsLoc=$ibdLoc$'results/'

rawIBDLoc=$ibdDataLoc$'raw/'
# phenosIBD_all=$ibdDataLoc$'phenos_all'
# ibd_all=$ibdDataLoc$'ibd_all'

# the 4 'core' cohorts, where UC/CD are considered together as IBD
ibdarray=( 'GWAS1M25.vcfs/' 'GWAS2M25.vcfs/' 'GWAS3M25.vcfs/' ) #  
ibdarray_target=( 'gwas1/' 'gwas2/' 'gwas3/' ) # 'ichip/'

# the 4 subpheno arrays
ibdarray_target_subpheno=( 'gwas3_cd/' 'gwas3_uc/' ) # 
ibdarray_subpeno_parent=( 'gwas3/' 'gwas3/' ) # where the subpheno parent cohorts are located

# where the original team152 UC/CD phenos are 
ibdarray_subpheno_origsamples=( $ibdCommonDataScratch$'GWAS3_orig.cd' $ibdCommonDataScratch$'GWAS_orig.uc' ) # 


# extended version of the above arrays that have the 2 subphenos
ibdarray_extended=( 'GWAS1M25.vcfs/' 'GWAS2M25.vcfs/' 'GWAS3M25.vcfs/' 'GWAS3M25.vcfs/'  'GWAS3M25.vcfs/' ) #  'ichipM25.vcfs/'
ibdarray_target_extended=( 'gwas1/' 'gwas2/' 'gwas3/' 'gwas3_cd/' 'gwas3_uc/' ) # 'ichip/'

ibdarray_subphenos_extended=( 'CD' 'UC' 'IBD' 'CD' 'UC' ) # index matched subpheno lookups for when plotting
subpheno_h2=( -1 -1 0.2541 0.3525 0.2228 ) # index matched LDSC h2 estimates
#qthresholds=( 0.25 0.5 0.75 0.8 0.85 0.9 0.95 ) # the lFDR thresholds used to export SNPs
qthresholds=( 0.96 0.97 0.98 0.99 1.0 ) # the lFDR thresholds used to export SNPs
Composite_PRS_workdir=${ibdCommonDataScratch}Composite_PRS_workdir/
GWAS3_IBD_assoc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[2]}PLINK.qassoc
GWAS3_CD_assoc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[3]}PLINK.qassoc
GWAS3_UC_assoc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[4]}PLINK.qassoc


# CD : 0.3525 (0.0804)
# UC : 0.2228 (0.0737)
# IBD: 0.2541 (0.0656)


########################################################################

# convert merged dataset into PLINK1
arguments=' --memory '$plinkMem$' --vcf '$ibd_raw_merged$' --make-bed --out '$ibdCommonDataScratch$'_temp --allow-no-sex '
$plink $arguments
#16154914 variants and 50612 people pass filters and QC.
# Total genotyping rate is 0.640815.
# the merged file does not have INFO tag, as that is specific for each cohort


# have a file for HLA exclusion range
echo -e $HLA_chrom$" "$HLA_start$" "$HLA_end$" HLA" > $ibdCommonDataScratch$'HLA_range'



arraylength=${#ibdarray[@]}
for (( j=1; j<${arraylength}+1; j++ )); do

mkdir -p ${ibdCommonDataScratch}${ibdarray_target[$j-1]}
mkdir -p ${ibdCommonDataScratch}${ibdarray_target[$j-1]}logs/
echo ${ibdarray_target[$j-1]}

#rm -rf ${ibdCommonDataScratch}${ibdarray_target[$j-1]}_allSNPs
for ((i=1; i<=$numChroms; i++)); do

# FARM3!!  loop through all datasets, export out VCFs again, but only with INFO > 0.8
####commands="/software/team152/bcftools/bcftools1.6/bcftools/./bcftools view -i 'MIN(INFO>0.8)' ${ibd_raw_cohortbaseLoc}${ibdarray[$j-1]}${i}.filtered_dec.vcf.gz  -Oz -o ${ibdCommonDataScratch}${ibdarray_target[$j-1]}${i}_info8.vcf.gz"
####bsub -q long -G team152 -n3 -J ${ibdarray_target[$j-1]}${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o ${ibdCommonDataScratch}${ibdarray_target[$j-1]}logs/${i}.out -e ${ibdCommonDataScratch}${ibdarray_target[$j-1]}logs/${i}.err  "$commands"

done # end of chroms loop
done # end of cohorts loop



# convert these to PLINK1
#########################
arraylength=${#ibdarray[@]}
for (( j=1; j<${arraylength}+1; j++ )); do

plinkFileList=$ibdCommonDataScratch$'mergelist.txt'
rm -f $plinkFileList
#rm -rf ${ibdCommonDataScratch}${ibdarray_target[$j-1]}_allSNPs
for ((i=1; i<=$numChroms; i++)); do

rangeExcl=''
if [ $i == $HLA_chrom ] ; then  
rangeExcl=' --exclude range '$ibdCommonDataScratch$'HLA_range'
fi 
# filter merged dataset, apply a geno filter of 0.98, and MAF filter of 0.1% , and the remaining SNP filter
filename=${ibdCommonDataScratch}${ibdarray_target[$j-1]}${i}_info8
arguments=' --memory '$plinkMem$' --vcf '$filename$'.vcf.gz  --geno 0.02 --maf 0.001 --make-bed --snps-only --bp-space 1 '$rangeExcl$' --out '$filename$' --allow-no-sex --double-id'
$plink $arguments

# need to find duplicate SNPs such as: 
# 16	rs34055956	0	65158106	A	G
# 16	rs34055956	0	65158107	T	C
# 8	rs34077128	0	144312508	A	G
# 8	rs34077128	0	144312509	T	C
# 4	rs34235198	0	159398228	A	G
# 4	rs34235198	0	159398229	C	T

# find duplicate IDs, AND filter by the baspair position ,so that we keep one of them
# only keep uniques:
awk ' { a[$2]++; if( a[$2] > 1){ print $0} }' $filename$'.bim' > $filename$'_duplicatedSNPs'
arguments=' --memory '$plinkMem$' --bfile '$filename$' --exclude '$filename$'_duplicatedSNPs --make-bed --out '$filename$'_nodupSNPs --allow-no-sex --double-id'
$plink $arguments

rm -rf $filename$'.bed'
rm -rf $filename$'.bim'
rm -rf $filename$'.fam'

echo $filename$'_nodupSNPs' >> ${plinkFileList}

done # end of chrom loop

# merge all .bims, enforce same order as chrom 1
#plinkFileList=$ibdCommonDataScratch$'mergelist.txt'
filename=${ibdCommonDataScratch}${ibdarray_target[$j-1]}
arguments=' --memory '$plinkMem$' --merge-list '$plinkFileList$' --make-bed --out '$filename$'all --allow-extra-chr --allow-no-sex --indiv-sort f '$filename$'1_info8_nodupSNPs.fam'
$plink $arguments
done # end of cohorts loop



# ADDITIONAL QC STEPS:
#####################
arraylength=${#ibdarray[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
cohort_folder=${ibdarray_target[$j-1]}
cohort_name=${cohort_folder::-1} # need to remove the last character so that gwas1/ -> gwas1
filename=${ibdCommonDataScratch}${ibdarray_target[$j-1]}



# 1) Compare HWE/MAF before/after converting to PLINK( if it is within +- 5%, otherwise throw out )
# this is to ensure that the lossy conversion does not change the data too much #

# the data converted QC-d so far may have a different subset of indis/SNPs to the raw
awk '{  {print $2} }' $filename$'all.fam' > $filename$'allIndis'
awk '{  {printf $2" "} }' $filename$'all.bim' > $filename$'all_snps'
rm -rf  $filename$'ALL_HWE_MAF_fail'
for ((i=1; i<=$numChroms; i++)); do


rawFile=${ibd_raw_cohortbaseLoc}${ibdarray[$j-1]}${i}.filtered_dec.vcf.gz
arguments=' -g '$rawFile$' -filetype vcf -incl-rsids '$filename$'all_snps -incl-samples '$filename$'allIndis -snp-stats -osnp '$filename$i$'_rawStats'
$qctool2_new $arguments

arguments=' --memory '$plinkMem$' --bfile '$filename$'all --freq --hardy --out '$filename$i$'_convertedStats --chr '$i$' --allow-extra-chr --allow-no-sex'
$plink $arguments

# cache RSID -> MAF ( skip the 11 header lines ( col 14 = MAF)
awk 'FNR == NR {  if(FNR>11) { file1[ $2 ] = $14; } next; } FNR <= NR { { if ( $2 in file1 ) { print $2"\t"file1[$2]"\t"$5 } } }' $filename$i$'_rawStats' $filename$i$'_convertedStats.frq' > $filename$i$'_MAF_diag'

# cache RSID -> HWE ( skip the 11 header lines ( col 8 = HWE)
awk 'FNR == NR {  if(FNR>11) { file1[ $2 ] = $8; } next; } FNR <= NR { { if ( $2 in file1  ) { print $2"\t"file1[$2]"\t"$9 } } }' $filename$i$'_rawStats' $filename$i$'_convertedStats.hwe' > $filename$i$'_HWE_diag'

# check if raw/converted MAF/HWE are within +- 0.95, and if not add it to exclude list
awk -v log_10=2.302585 '{ 
# for p-values we use their -log10 
raw=-log($2) / log_10;
conv=-log($3) / log_10;
if (raw == 0) { raw = 0.0001; conv = conv +0.0001 }
diff = raw - conv
if (sqrt(diff*diff) / raw > 0.05) {print $1} }' $filename$i$'_HWE_diag' > $filename$i$'_HWE_fail'

awk  '{ 
raw=$2
conv=$3
if (raw == 0) { raw = 0.0001; conv = conv +0.0001 }
diff = raw - conv
if (sqrt(diff*diff) / raw > 0.05) {print $1} }' $filename$i$'_MAF_diag' > $filename$i$'_MAF_fail'


# concat exclude list for each chrom
cat $filename$i$'_HWE_fail' $filename$i$'_MAF_fail' > $filename$i$'_HWE_MAF_fail'
cat $filename$i$'_HWE_MAF_fail' >>  $filename$'ALL_HWE_MAF_fail'

done # end of chroms loo


# find the unique 
awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $filename$'ALL_HWE_MAF_fail' > $filename$'ALL_HWE_MAF_fail_uniqes'
numConvFailed=$(wc -l < "$filename"ALL_HWE_MAF_fail_uniqes)
echo -e $numConvFailed > $filename$'NUM_HWE_MAF_fails'
echo "Number of SNPs failed thresholding conversion to PLINK1: "$numConvFailed

#GWAS1 12014
#GWAS2 15736
#GWAS3 12539
# (12014 + 15736 + 12539 )/3 = 13429.67, 13,430

# remove SNPs that failed the conversion test
arguments=' --memory '$plinkMem$' --bfile '$filename$'all --exclude '$filename$'ALL_HWE_MAF_fail_uniqes --make-bed --out '$filename$'all_convQC --allow-extra-chr --allow-no-sex'
$plink $arguments

#####################


# the names of the phenotype files, these are of the format: $ibdCommonDataScratch$cohort_name$'_found' eg: $ibdCommonDataScratch$'newwave_found'
awk ' FNR == NR {
{ 
n=split($2,a,/_/);  # the old pheno files have retarded ID_ID nomenclature
len=int(length(a));
IID = tolower(a[len])
sex[IID] = $5;
phe[IID] = $6;

 next;
 }
} FNR <= NR { {  
n=split($2,a,/_/);  # the old pheno files have retarded ID_ID nomenclature
len=int(length(a));
IID = tolower(a[len])
if ( IID in phe ) { print $1" "$2" "$3" "$4" "sex[IID]" "phe[IID] }}  }' $ibdCommonDataScratch$cohort_name$'_found' $filename$'all_convQC.fam' > $filename$'all_convQC_found.keep'



wc -l $filename$'all_convQC.fam'  # 5880
wc -l $filename$'all_convQC_found.keep' # 
#head $filename$'all_convQC_found.keep'
#head $filename$'all_convQC.fam'

# match phenos and SEX for all_convQC cohorts  (only keep those that have both sex and the pheno)
arguments=' --memory '$plinkMem$' --bfile '$filename$'all_convQC --keep '$filename$'all_convQC_found.keep --indiv-sort f '$filename$'all_convQC_found.keep --make-bed --allow-no-sex --out '$filename$'all_found'
$plink $arguments

# overwrite the phenos with the ones that have sex/pheno
cp $filename$'all_found.fam' $filename$'all_found.fam.old'
cp $filename$'all_found.keep' $filename$'all_found.fam'

done # end of cohort loop

###########
# Branch off point for subphenotypes

# Loop for the 4 subphenos, 
arraylength=${#ibdarray_target_subpheno[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
cohort_folder=${ibdarray_target_subpheno[$j-1]}
cohort_name=${cohort_folder::-1} # need to remove the last character so that gwas1/ -> gwas1
filename=${ibdCommonDataScratch}${ibdarray_target_subpheno[$j-1]}



# create their directories
mkdir -p $filename  # the scratch
mkdir -p $ibdDataLoc$cohort_folder  # their main data
mkdir -p $ibdDataLoc$cohort_folder  # their main data
mkdir -p $ibdDataLoc$cohort_folder$'scratch/'  # their main data scratch


echo $cohort_folder


# create the 'all found' versions by subsetting their 'parent' cohorts to the subphenotypes
parentData=${ibdCommonDataScratch}${ibdarray_subpeno_parent[$j-1]}all_found

#                          keep the original 5 cols, that includes SEX from the original fam
awk 'FNR == NR {	file1[ $2 ] = $1"\t"$2"\t"$3"\t"$4"\t"$5;	next;}FNR <= NR {
 if ( $2 in file1 ) {  print file1[$2]"\t"$5+1  }  }
' OFS='\t' $parentData$'.fam' ${ibdarray_subpheno_origsamples[$j-1]} > $filename$'all_found.fam'
# print the 5th col +1 (as the .sample files are 0/1, where as in PLINk we need 1/2 for control/case


arguments=' --memory '$plinkMem$' --bfile '$parentData$' --keep '$filename$'all_found.fam --indiv-sort f '$filename$'all_found.fam --make-bed --allow-no-sex --out '$filename$'all_found'
$plink $arguments


done




###########
arraylength=${#ibdarray_extended[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
cohort_folder=${ibdarray_target_extended[$j-1]}
cohort_name=${cohort_folder::-1} # need to remove the last character so that gwas1/ -> gwas1
filename=${ibdCommonDataScratch}${ibdarray_target_extended[$j-1]}


###########
# filter out SNPs with p<10e-5 HWE in controls
awk ' { { if ( $6 == 1 ) { print $0 } } }' OFS='\t' $filename$'all_found.fam' > $filename$'all_found.controls'  # keeplist for controls
arguments=' --memory '$plinkMem$' --bfile '$filename$'all_found --hardy --out '$filename$'all_found_HWE --allow-no-sex --keep '$filename$'all_found.controls' # perform HWE test in controls only
$plink $arguments

# create a keeplist of variants, where the p>10e-5 
awk ' { { if ( $3 == "ALL" && $9 < 0.00001 ) { print $2 } } }' $filename$'all_found_HWE.hwe' > $filename$'all_found_HWE_fail'
wc -l $filename$'all_found_HWE_fail' # 10e-7: 19908,  10e-5:32550


# filter out SNPs with p<10e-7 HWE overall:
arguments=' --memory '$plinkMem$' --bfile '$filename$'all_found --hardy --out '$filename$'all_found_HWE_all --allow-no-sex' # perform HWE test on everyone
$plink $arguments

awk ' { { if ( $3 == "ALL" && $9 < 0.0000001 ) { print $2 } } }' $filename$'all_found_HWE_all.hwe' > $filename$'all_found_HWE_fail_all'
wc -l $filename$'all_found_HWE_fail_all' # 30304

cat $filename$'all_found_HWE_fail' $filename$'all_found_HWE_fail_all' > $filename$'all_found_HWE_fail_caseControls'



arguments=' --memory '$plinkMem$' --bfile '$filename$'all_found --exclude '$filename$'all_found_HWE_fail_caseControls --make-bed --allow-no-sex --out '$filename$'all_found_HWE'
$plink $arguments

# remove all previous temporary files
# rm -rf $filename$'all_found_HWE.hwe'
# rm -rf $filename$'all_found.bed'
# rm -rf $filename$'all_found.bim'
# rm -rf $filename$'all_found.fam'
# rm -rf $filename$'all.bed'
# rm -rf $filename$'all.bim'
# rm -rf $filename$'all.fam'


###########################################

# PCA on the TOTAL ( using the SNPs selected for this in the UKBB)


# perform PCA on the above in FLashPCA
arguments=' --bfile '$filename$'all_found_HWE --extract '$ibdCommonDataScratch$'SNP_PCAsubset --indiv-sort 0 --allow-no-sex --make-bed --out '$filename$'_PCA'
$plink $arguments
# 77897 variants remaining.
# UKBB: 147,604
# GWAS1 73844
# GWAS2 98097 
# GWAS3 78814
# ( 73844 + 98097 + 78814) / 3 = 83585
wc -l $filename$'_PCA.bim'

#awk '{ print $1" "$2" "$3" "$4" "$5" -9"}' $ibdCommonDataScratch$'_INFOfiltered_pheno_PCA.fam' > $ibdCommonDataScratch$'_INFOfiltered_pheno_PCA.fam2'

arguments=' -n7 --verbose -d 20 --bfile '$filename$'_PCA --outpc '$filename$'_PCA_pc  --outload '$filename$'_PCA_loadings --outmeansd '$filename$'_PCA_meansd'
$flashpca $arguments


#############################################
# build covariate design matrix: SEX, BATCH, PC1-20

# export out the phenos
awk ' { print $1 "\t" $2 "\t" $6 }' $filename$'all_found_HWE.fam' > $filename$'all_found_HWE.pheno'
awk ' { print $3-1 }' $filename$'all_found_HWE.pheno' > $filename$'pheno' # binary phenos need to be 0-1

# export out SEX ( add a header)
awk 'BEGIN{print "SEX"}; { print $5 };' $filename$'all_found_HWE.fam' > $filename$'all_found_HWE.sex'

# cut the FID/IID bits from the PCA file
cut -f1,2 --complement $filename$'_PCA_pc' > $filename$'_PCA_pc_NOID'

# paste together the overall covars
paste $filename$'all_found_HWE.sex' $filename$'_PCA_pc_NOID' > $filename$'_covariates'


done # end of cohort loop





#############################################
# build covariate design matrix: SEX, BATCH, PC1-20
arraylength=${#ibdarray_extended[@]}
for (( j=1; j<${arraylength}+1; j++ )); do

filename=${ibdCommonDataScratch}${ibdarray_target_extended[$j-1]}
# pheno regress on ALL to create pheno residuals (this has to be for all
awk ' { print $3-1 }' $filename$'all_found_HWE.pheno' > $filename$'phenotypes'
#arguments='/nfs/users/nfs_m/mk23/scripts/phenoRegress_binary.R '$filename$'phenotypes '$filename$' residuals '$filename$'_covariates 1 1 1 0 0'
#$Rscript  $arguments # using pdf output now, which works on dgx server


# Check if we get the same result via Backward selelction@ YES
#                                                                        1                    2               3                      4           5 6 7 8 9
arguments='/nfs/users/nfs_m/mk23/scripts/phenoRegress_backward.R '$filename$'phenotypes '$filename$' residuals_backwards '$filename$'_covariates 1 1 1 0 0'
$Rscript $arguments 
head $filename$'residuals_significantCovariates.txt'
head $filename$'residuals_backwards_significantCovariates.txt'


# copy the regression results over 
regresultDir=$ibdResultsLoc$'regression/'${ibdarray_target_extended[$j-1]}
mkdir -p $regresultDir
cp $filename$'original.pdf' $regresultDir$'original.pdf'
cp $filename$'adjusted.pdf' $regresultDir$'adjusted.pdf'
cp $filename$'residuals_significantCovariates.txt' $regresultDir$'residuals_significantCovariates.txt'
cp $filename$'residuals_modelSummary_REDUCED' $regresultDir$'residuals_modelSummary'
cp $filename$'residuals_modelSummary_inner_REDUCED' $regresultDir$'residuals_modelSummary_inner'


# produce tables of effect p-val / Beta
# get the p values
awk 'FNR == NR { file1[ $1 ] = $1; next; } FNR <= NR { { if ( $1 in file1 ) 
{
covname=$1

if ($6 == "<" ) { print covname"\t"$6$7  }
else {print covname"\t"$6}
} } }
' OFS='\t' $regresultDir$'residuals_significantCovariates.txt' $regresultDir$'residuals_modelSummary' > $regresultDir$'regressiontable_temp'

awk 'FNR == NR { file1[ $1 ] = $1"\t"$2; next; } FNR <= NR { 
{ if ( $1 in file1 ) 
{ 
print file1[ $1 ]"\t"$2
} } }
' OFS='\t' $regresultDir$'regressiontable_temp' $regresultDir$'residuals_modelSummary_inner' > $regresultDir$'regressiontable_temp2'


# also want to include the covariates that do not have an overall effect (like the centre) so we just copy them over as 'NA'
awk 'FNR == NR { file1[ $1 ] = $0; next; } FNR <= NR { { 
if ( $1 in file1 ) { print file1[ $1 ]}
else {print $0"\tNA"} 
}}
' OFS='\t' $regresultDir$'regressiontable_temp2' $regresultDir$'regressiontable_temp' > $regresultDir$'regressiontable_temp3'


echo -e "COV\tP-value\tEffect_size" > $regresultDir$'regressiontable_header'
cat $regresultDir$'regressiontable_header' $regresultDir$'regressiontable_temp3' > $regresultDir$'regressiontable'

# cleanup
rm -rf $regresultDir$'regressiontable_temp'
rm -rf $regresultDir$'regressiontable_temp2'
rm -rf $regresultDir$'regressiontable_temp3'
rm -rf $regresultDir$'regressiontable_header'

# produce keeplist for the publicids
awk '{print $1"\t"$2}' $filename$'all_found_HWE.fam' > $filename$'publicIDkeeplist'

cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
mkdir -p $cohortDataFolder

# move the final 'all' files into the data folder
cp $filename$'all_found_HWE.bed' $cohortDataFolder$'ibd_all.bed'
cp $filename$'all_found_HWE.bim' $cohortDataFolder$'ibd_all.bim'

# Generate final phenotype file for this cohort to use: 
paste $filename$'publicIDkeeplist' $filename$'residuals' > $cohortDataFolder$'ibd_all.pheno'

# add in the phenotype residuals
cp $filename$'all_found_HWE.fam' $filename$'all_found_HWE.fam_origpheno'
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' $filename$'all_found_HWE.fam' > $filename$'all_found_HWE.fam_nopheno'
paste $filename$'all_found_HWE.fam_nopheno' $filename$'residuals' > $cohortDataFolder$'ibd_all.fam'


done

#####################################


# Move the 'ibd_all's back into the scratch folder, as they are NOT the final genotype files,
#  as I need to subset them so that they have the 'lowest common denominator' panel of SNPs acorss them
arraylength=${#ibdarray_extended[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
backupLoc=${ibdCommonDataScratch}${ibdarray_target_extended[$j-1]}

mv $cohortDataFolder$'ibd_all.bed' $backupLoc$'ibd_all.bed'
mv $cohortDataFolder$'ibd_all.bim' $backupLoc$'ibd_all.bim'
mv $cohortDataFolder$'ibd_all.fam' $backupLoc$'ibd_all.fam'

echo $cohortDataFolder
echo $backupLoc
done


# find the 'lowest common denominator for all 3+2 datasets
arraylength=${#ibdarray_extended[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
cohort_folder=${ibdarray_target_extended[$j-1]}
cohort_name=${cohort_folder::-1} # need to remove the last character so that gwas1/ -> gwas1
backupLoc=${ibdCommonDataScratch}${ibdarray_target_extended[$j-1]}


awk '{  {print $2}  }' $backupLoc$'ibd_all.bim' > $backupLoc$'_allSNPs'
#wc -l $backupLoc$'_allSNPs'
echo $j $cohort_name

# for all but the first GWAS' we intersect this cohorts SNPs with the previous ones
if [ $j -gt 1 ] ; then
#echo "later GWASes"


currentGWAS_SNPs=$backupLoc$'_allSNPs'

# find out how  to always pull in the previous lowest common denom SNPs rather than the previous cohorts SNPs
#cohort_folder_prev=${ibdarray_target_extended[$j-2]}
#cohort_name_prev=${cohort_folder_prev::-1} # need to remove the last character so that gwas1/ -> gwas1

if [ $j == 2 ] ; then
backupLoc_prev=${ibdCommonDataScratch}${ibdarray_target_extended[$j-2]}
prevGWAS_SNPs=$backupLoc_prev$'_allSNPs'
lastCommonSNPs=$prevGWAS_SNPs
fi # end of if it is actually the 2nd cohort

#echo "lastCommonSNPs "$lastCommonSNPs
currentCommonSNPs=$lastCommonSNPs$'_'$cohort_name
comm -12  <(sort $lastCommonSNPs) <(sort $currentGWAS_SNPs) >  $currentCommonSNPs
lastCommonSNPs=$currentCommonSNPs
wc -l $currentCommonSNPs

fi # end of if greater than 2nd Cohort
done
# 7,046,463 SNPS in tota

# extract the final common SNPs 
arraylength=${#ibdarray_extended[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
backupLoc=${ibdCommonDataScratch}${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
arguments=' --memory '$plinkMem$' --bfile '$backupLoc$'ibd_all  --extract '$currentCommonSNPs$' --allow-no-sex --make-bed --out '$cohortDataFolder$'ibd_all'
$plink $arguments
done

#####################################

# perform assoc in each cohort, on the 100% of the dataset. We don't withhold anything, this is just used for the QC steps and the manhattan plots
arraylength=${#ibdarray_extended[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
mkdir -p ${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
arguments=' --memory '$plinkMem$' --bfile '$cohortDataFolder$'ibd_all --assoc --allow-no-sex --out '$outputloc$'PLINK'
$plink $arguments
done


# Additional QC step : to remove 'noisy' associations, IE SNPs that do not have enough LD friends supporting the association
#################

arraylength=${#ibdarray_extended[@]}
for (( j=1; j<${arraylength}+1; j++ )); do

cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohort_folder=${ibdarray_target_extended[$j-1]}
cohort_name=${cohort_folder::-1} # need to remove the last character so that gwas1/ -> gwas1
diagOutFolder=$outputloc$'diag/'

scratchLoc=${ibdCommonDataScratch}${ibdarray_target_extended[$j-1]}
mkdir -p $diagOutFolder

rm -rf $diagOutFolder$cohort_name$'_ALL_LD_outliers_nofriends'
# need to split the data and results into choromosomes, as it makes no sense to produce an LD based exclusion criteria that is genome wide
for ((i=1; i<=$numChroms; i++)); do
# create a temporary assoc file 
awk -v chrom="$i" ' { if ($1 ==chrom || NR == 1 ){ print $0 } } ' $outputloc$'PLINK.qassoc' > $scratchLoc$'PLINK'$i$'.qassoc'
falsePositiveLDTest $cohortDataFolder$'ibd_all' $diagOutFolder $scratchLoc$'PLINK'$i$'.qassoc' $cohort_name$'_'$i
cat $diagOutFolder$cohort_name$'_'$i'_LD_outliers_nofriends' >> $diagOutFolder$cohort_name$'_ALL_LD_outliers_nofriends'
done # end of chrom loop

done # end of cohort loop
#########################
# Issue with intra-chromosome D tests: what if there is only 1 SNP, there is no SD to measure it against to find out if it is outlier
arraylength=${#ibdarray_extended[@]}

for (( j=1; j<${arraylength}+1; j++ )); do
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortScratch=${ibdCommonDataScratch}${ibdarray_target_extended[$j-1]}
cohort_folder=${ibdarray_target_extended[$j-1]}
cohort_name=${cohort_folder::-1} # need to remove the last character so that gwas1/ -> gwas1
diagOutFolder=$outputloc$'diag/'


genomeWide_D=$diagOutFolder$'genomeWide_D'
rm -rf $genomeWide_D
echo -e "RSID\tp_value\tD" > $genomeWide_D

# concat all chrom Ds and nofriends into complete lists, and rerun the Rscripts on them
for ((i=1; i<=$numChroms; i++)); do
title=$cohort_name$'_'$i


if [ -s "$diagOutFolder$title"_allSigSNPs ] ; then
while read p; do # go through each SNP
p=$(echo $p | sed -e 's/\r//g')  # the last column contains a newline character, which would fuck up awk string comparison and create unreadable files etc
IFS=' ' read -r -a arrIN <<< "$p"  # https://stackoverflow.com/questions/10586153/split-string-into-an-array-in-bash
targetVar=${arrIN[0]} 
targetVar_P=${arrIN[1]}


# the Rscript writes a file that has the chisq p val in the format of "D=pval"
source $diagOutFolder$title$targetVar$'.chisqp'



# add this to a running total 
echo -e $targetVar$"\t"$targetVar_P$"\t"$D >> $genomeWide_D


done <$diagOutFolder$title$'_allSigSNPs'
fi # significant SNPs exist for chrom



done # end of chrom loop



# run Rscript to produce a genome-wide outliers etc
arguments='/nfs/users/nfs_m/mk23/scripts/LD_distr.R '$genomeWide_D$' '$diagOutFolder$' '$cohort_name$'_genomeWide'
$Rscript $arguments 

# need to add the SNPs with no friends, otherwise they won't get excluded:
# actually not as these will be on the per chromosome list too, as that does not depend on overall D or SD(D)
# cat $diagOutFolder$cohort_name'_genomeWide_LD_outliers' $diagOutFolder$cohort_name$'_ALL_LD_outliers_nofriends' > $diagOutFolder$cohort_name'_genomeWide_LD_outliers_nofriends'


# cat this list and the per-Chrom exclude list, and get the uniques
cat $diagOutFolder$cohort_name$'_ALL_LD_outliers_nofriends' $diagOutFolder$cohort_name'_genomeWide_LD_outliers' > $diagOutFolder$cohort_name'_genomeWide_LD_outliers_nofriends_COMBINED_dupes'
awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $diagOutFolder$cohort_name'_genomeWide_LD_outliers_nofriends_COMBINED_dupes' > $diagOutFolder$cohort_name'_genomeWide_LD_outliers_nofriends_COMBINED'

wc -l $diagOutFolder$cohort_name$'_ALL_LD_outliers_nofriends'
wc -l $diagOutFolder$cohort_name'_genomeWide_LD_outliers_nofriends_COMBINED'

done # end of cohorts loop





#########################



# create --exclude variants list that failed the LD-friend test, in ANY of the 5 cohorts (so that the PRS will be compatible
arraylength=${#ibdarray_extended[@]}
rm -rf ${ibdCommonDataScratch}SNP_ALL_LD_outliers_nofriends_dupes
for (( j=1; j<${arraylength}+1; j++ )); do
cohort_folder=${ibdarray_target_extended[$j-1]}
cohort_name=${cohort_folder::-1} # need to remove the last character so that gwas1/ -> gwas1
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
diagOutFolder=$outputloc$'diag/'
# create a UNION of all SNPs that failed across the 5 cohorts
wc -l $diagOutFolder$cohort_name$'_genomeWide_LD_outliers_nofriends_COMBINED'
cat $diagOutFolder$cohort_name$'_genomeWide_LD_outliers_nofriends_COMBINED' >> ${ibdCommonDataScratch}SNP_ALL_LD_outliers_nofriends_dupes
done

wc -l ${ibdCommonDataScratch}SNP_ALL_LD_outliers_nofriends_dupes 
awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' ${ibdCommonDataScratch}SNP_ALL_LD_outliers_nofriends_dupes > ${ibdCommonDataScratch}SNP_ALL_LD_outliers_nofriends
wc -l ${ibdCommonDataScratch}SNP_ALL_LD_outliers_nofriends # 748 in total


### Remove the above variants from both the results as well as from the plink files 
arraylength=${#ibdarray_extended[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohort_folder=${ibdarray_target_extended[$j-1]}
cohort_name=${cohort_folder::-1} # need to remove the last character so that gwas1/ -> gwas1
diagOutFolder=$outputloc$'diag/'
cohortScratch=${ibdCommonDataScratch}${ibdarray_target_extended[$j-1]}

# remove the above from ALL of the results as well as from the ibd_all's to produce a uniformly QCd dataset
awk ' FNR == NR {file1[ $1 ] = $1; next;} FNR <= NR { {  if ( $2 in file1 == 0 || $1 == "CHR" ) { print $0 } } }' ${ibdCommonDataScratch}SNP_ALL_LD_outliers_nofriends $outputloc$'PLINK.qassoc' > $outputloc$'PLINK.qassoc_filtered'

# replace the original association file with the new clean one
mv $outputloc$'PLINK.qassoc' $outputloc$'PLINK.qassoc.orig'
mv $outputloc$'PLINK.qassoc_filtered' $outputloc$'PLINK.qassoc'


# REMOVE QC_FAIL_SNPs, as the $cohortDataFolder$'ibd_all' will be used to produce PWAS/TWAS and other derived datasets, so it cannot have any bad SNPs in it

mv $cohortDataFolder$'ibd_all.bim' $cohortScratch$'ibd_all2.bim'
mv $cohortDataFolder$'ibd_all.fam' $cohortScratch$'ibd_all2.fam'
mv $cohortDataFolder$'ibd_all.bed' $cohortScratch$'ibd_all2.bed'
arguments=' --memory '$plinkMem$' --bfile '$cohortScratch$'ibd_all2 --exclude '$ibdCommonDataScratch$'SNP_ALL_LD_outliers_nofriends --allow-no-sex --make-bed --out '$cohortDataFolder$'ibd_all'
$plink $arguments

done # end of cohort loop


#################

## MANHATTAN FOR IBD: FARM5
arraylength=${#ibdarray_extended[@]}
for (( j=1; j<${arraylength}+1; j++ )); do
mkdir -p ${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}

cohort_folder=${ibdarray_target_extended[$j-1]}
cohort_name=${cohort_folder::-1} # need to remove the last character so that gwas1/ -> gwas1

#produceManhattanPlot_PLINK $outputloc$'PLINK.qassoc' $outputloc$'PLINK_'$cohort_name 'GWAS_'$cohort_name$' 7.30103'

pheno_current=${ibdarray_subphenos_extended[$j-1]}
title=$cohort_name$'_'$pheno_current
echo "MANHATTAN FOR: "$title

bsub -q yesterday -G team152 -n1 -J ${cohort_name} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o ${outputloc}Manh.out -e ${outputloc}Manh.err  "/nfs/users/nfs_m/mk23/scripts/dgx/produceManhattanPlot_PLINK.sh $outputloc$'PLINK.qassoc' $outputloc $title$' 7.30103'"
done


#################################################################################
### Bootstraps

# Need different Bootstraps for GWAS3 IBD/CD/UC (as CD/UC have different subsets of indis
# we also need to actually generate and keep PLINK files, due to how PLINK doesn't allow --keep with the same indis features multiple time...
arraylength=${#ibdarray_extended[@]}
for (( j=4; j<${arraylength}+1; j++ )); do

cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
rm -rf $bootStrapLoc
mkdir -p $bootStrapLoc

arguments=$scriptsLoc$'bootStrapGen_plink.R  '$cohortDataFolder'ibd_all.fam '$bootStrapLoc$' '$filenameStem$'_f '$numBootStraps
$Rscript $arguments




# 1) intersect the final QCd genotype data with the  HAPmap3 b37 panel of SNPs ($hapmap3_b37bim) -> existing_hapmap3
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortScratch=${ibdCommonDataScratch}${ibdarray_target_extended[$j-1]}
awk 'NR>1 {print $1"\t"$2"\t"$3}' $outputloc$'PLINK.qassoc' > $cohortScratch$'ALL_SNPs'

awk 'FNR == NR {
coordId=$1"_"$4
file1[ coordId ] = coordId; next; } 
FNR <= NR { { coordId_current=$1"_"$3; if ( coordId_current in file1) { print $2 } }  }
'  $hapmap3_b37bim $cohortScratch$'ALL_SNPs' > $cohortScratch$'existing_hapmap3_ibd'
wc -l $cohortScratch$'existing_hapmap3_ibd' # 1124113 of the 1,403,851




# generate the Bootstraps
for ((i=1; i<=numBootStraps; i++)); do

foldFam=$filenameStem$'_f'$i
arguments=' --memory '$plinkMem$' --threads 5 --bfile '$cohortDataFolder$'ibd_all --keep '$bootStrapLoc$foldFam$'_uniques.fam --make-bed --out '$bootStrapLoc$foldFam$'_uniques --allow-no-sex --extract '$cohortScratch$'existing_hapmap3_ibd'
$plink $arguments

# find out how many level of duplicates there were  
source $bootStrapLoc$foldFam$'_numDupes'

plinkFileList=$bootStrapLoc'/mergelist'$i$'.txt'
rm -f $plinkFileList

for ((c=1; c<=numDupeLevels; c++)); do
cp $bootStrapLoc$foldFam$'_dupes_level'$c$'.fam' $bootStrapLoc$foldFam$'_dupes_level'$c$'.fam.orig' # need to create a backup copy BEFORE we make PLINK dataset below, as that will overwrite these .fams with different order, and the mapping between duplicate123 -> UKBBID would be lost
arguments=' --memory '$plinkMem$' --threads 50 --bfile '$cohortDataFolder$'ibd_all --keep '$bootStrapLoc$foldFam$'_dupes_level'$c$'.fam  --indiv-sort f '$bootStrapLoc$foldFam$'_dupes_level'$c$'.fam --make-bed --out '$bootStrapLoc$foldFam$'_dupes_level'$c$' --allow-no-sex --extract '$cohortScratch$'existing_hapmap3_ibd'
$plink $arguments

# give the duplicate individuals a different name, so plink wont notice when we merge them... we carefully preserved the order of the .fam and new_name.fam so that simply overwriting the .fam file, will NOT permute the phenotypes
cp $bootStrapLoc$foldFam$'_dupes_level'$c$'_newname.fam' $bootStrapLoc$foldFam$'_dupes_level'$c$'.fam'
echo $bootStrapLoc$foldFam$'_dupes_level'$c >> ${plinkFileList}
done # end of dupoeLevels



echo $bootStrapLoc$foldFam$'_uniques' >> ${plinkFileList}
arguments=' --memory '$plinkMem$' --merge-list '$plinkFileList$' --make-bed --out '$bootStrapLoc$foldFam$' --allow-extra-chr --allow-no-sex'
$plink $arguments


#arguments=' --memory '$plinkMem$' --bfile '$bootStrapLoc$foldFam$'_uniques  --threads 50 --bmerge '$bootStrapLoc$foldFam$'_dupes --make-bed --out '$bootStrapLoc$foldFam$' --allow-no-sex'
#$plink $arguments


# remove the 'dupes' and unique .beds/bims to save space
for ((c=1; c<=numDupeLevels; c++)); do
rm -f $bootStrapLoc$foldFam$'_dupes_level'$c$'.bed'
rm -f $bootStrapLoc$foldFam$'_dupes_level'$c$'.bim'
done # end of dupe levels
rm -f $bootStrapLoc$foldFam$'_uniques.bed'
rm -f $bootStrapLoc$foldFam$'_uniques.bim'



# create the validation set for the Bootstrap too (the .fam file for which R created)
foldFam=$filenameStem$'_f'$i
arguments=' --memory '$plinkMem$' --bfile '$cohortDataFolder$'ibd_all --keep '$bootStrapLoc$foldFam$'_valid.fam --allow-no-sex --make-bed --out '$bootStrapLoc$foldFam$'_valid --allow-extra-chr --extract '$cohortScratch$'existing_hapmap3_ibd'
$plink $arguments

# extract phenos 
awk '{print $1"\t"$2"\t"$6}' $bootStrapLoc$foldFam$'.fam' > $bootStrapLoc$foldFam$'.pheno'
awk '{print $1"\t"$2"\t"$6}' $bootStrapLoc$foldFam$'_valid.fam' > $bootStrapLoc$foldFam$'_valid.pheno'

done # end of Bootstrap loop


done # end of cohort loop


