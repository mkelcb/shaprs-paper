####################################################
# dgx-server screens:


# screen -r -D 14485.pts-0.node-11-1-1

# screen -r -D 5099.pts-0.node-11-1-2  # QC test 2


# Switching IBD from 'post' processing to pre-processing

################################

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

gencode_map='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/0_common/data/gencode_map'
ukbb_imputed='/lustre/scratch115/realdata/mdt3/projects/ukbiobank/FullRelease/Imputed/EGAD00010001474/'

commonLoc=$baseLoc$'0_common/'
commonDataLoc=$commonLoc$'data/'
commonScratchLoc=$commonDataLoc$'scratch/'
SNP_conv_HWE_FAIL=$commonDataLoc$'SNP_conv_HWE_FAIL' # the ukbb SNPs that failed conversion or HWE p < 10e-7 filters
ukbbQCkeeplist=$commonDataLoc$'ukbb_hq_keeplist.txt' # Headerless file, col1: TeamID (for imputed data) col2: PublicID (for the genotypedata): 1587700 5172879 // 376510 indis
ukbbQCkeeplist_publicDs=$ukbbQCkeeplist$'_pubIDs'



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

qthresholds_string=$(echo ${qthresholds[*]// /|})  # convert an array to string: https://stackoverflow.com/questions/13470413/converting-a-bash-array-into-a-delimited-string
num_thresholds=${#qthresholds[@]}
Composite_PRS_workdir=${ibdCommonDataScratch}Composite_PRS_workdir/
GWAS3_IBD_assoc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[2]}PLINK.qassoc
GWAS3_CD_assoc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[3]}PLINK.qassoc
GWAS3_UC_assoc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[4]}PLINK.qassoc

filenameStem="GWAS"
ldrefhapmap3="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/ldrefhapmap3.snps"

outLdpredName="ldpred2_PRS"
baseLDpredLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/"
NCORES=15
GWAS1Folder=$ibdDataLoc${ibdarray_target_extended[0]}
GWAS2Folder=$ibdDataLoc${ibdarray_target_extended[1]}

rho=0.2852094 # calculated via the approximation formula


diagFolder='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/scratch/'
diagOutLoc=$diagFolder$'corrDiag'

results_GWAS2UC_UC_PRS=$ibdResultsLoc$'PRS2/GWAS2UC_UC_PRS/'
results_GWAS1CD_CD_PRS=$ibdResultsLoc$'PRS2/GWAS1CD_CD_PRS/'
results_GWAS2UC_IBD_PRS=$ibdResultsLoc$'PRS2/GWAS2UC_IBD_PRS/'
results_GWAS1CD_IBD_PRS=$ibdResultsLoc$'PRS2/GWAS1CD_IBD_PRS/'
results_GWAS2UC_Blend_PRS=$ibdResultsLoc$'PRS2/GWAS2UC_Blend_PRS/'
results_GWAS1CD_Blend_PRS=$ibdResultsLoc$'PRS2/GWAS1CD_Blend_PRS/'
results_GWAS2UC_Comp_PRS=$ibdResultsLoc$'PRS2/GWAS2UC_Comp_PRS/'
results_GWAS1CD_Comp_PRS=$ibdResultsLoc$'PRS2/GWAS1CD_Comp_PRS/'
results_GWAS2UC_SMTPred_PRS=$ibdResultsLoc$'PRS2/GWAS2UC_SMTPred_PRS/'
results_GWAS1CD_SMTPred_PRS=$ibdResultsLoc$'PRS2/GWAS1CD_SMTPred_PRS/'
N_effs=( 0 0 16852 10874 10782 ) 


source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh
conda activate mk23_pytorch

mkdir -p $diagFolder



############################################################################
# I) extract the UC/ CD test sets
# need to make sure that the Test Set PLINK files have the same number of predictors in the same order as the .qassoc files
arraylength=${#ibdarray_extended[@]}
for (( j=3; j<${arraylength}+1; j++ )); do # GWAS3 IBD/CD/UC
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
LDPredOutLoc=$bootStrapLoc$'LDpred/'


mkdir -p $LDPredOutLoc

for ((i=1; i<=numBootStraps; i++)); do # go through all bootstrap
foldFam=$filenameStem$'_f'$i


if [ $j == 3 ] || [ $j == 4 ] ; then
echo "Combined or CD"
# CD - GWAS1
testSet=$bootStrapLoc$foldFam$'TEST_hapmap3'
arguments=' --memory '$plinkMem$' --threads 1 --bfile '$GWAS1Folder$'ibd_all --extract '$bootStrapLoc$foldFam$'.bim  --make-bed --out '$testSet
$plink $arguments
fi 



if [ $j == 3 ] || [ $j == 5 ] ; then
echo "Combined or UC"
# UC - GWAS2
testSet=$bootStrapLoc$foldFam$'TEST_hapmap3_UC'
arguments=' --memory '$plinkMem$' --threads 1 --bfile '$GWAS2Folder$'ibd_all --extract '$bootStrapLoc$foldFam$'.bim  --make-bed --out '$testSet
$plink $arguments
fi 

# convert the above association into a summary statistics format suitable for the shaPRS/smtpred
# PLINK qassoc signature:
# CHR          SNP         BP    NMISS       BETA         SE         R2        T            P
avgNumIndis=$(wc -l < "$bootStrapLoc$foldFam".fam)
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
' OFS='\t' $testSet$'.bim' $bootStrapLoc$foldFam$'PLINK.qassoc' > $bootStrapLoc$foldFam$'_phe_sumstats'


done
done



foldFam=$filenameStem$'_f2'
$bootStrapLoc$foldFam$'PLINK.qassoc'
head $bootStrapLoc$foldFam$'PLINK.qassoc'

############################################################################

# arguments=' --memory '$plinkMem$' --threads 5 --bfile '$bootStrapLoc$foldFam$' --assoc --out '$bootStrapLoc$foldFam$'PLINK'
# $plink $arguments
j=5
i=1

# 2) Adjust the SNP betas in LDPred for the 3 baselines: UC, CD, IBD
arraylength=${#ibdarray_extended[@]}
for (( j=3; j<${arraylength}+1; j++ )); do # GWAS3 IBD/CD/UC
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
LDPredOutLoc=$bootStrapLoc$'LDpred/'
binPheno='0' # we actually ran the association as quantiative phenotypes
avgNumIndis=${N_effs[$j-1]}
mkdir -p $LDPredOutLoc


for ((i=1; i<=numBootStraps; i++)); do # go through all bootstraps to create association files from the training set
foldFam=$filenameStem$'_f'$i

# convert the above association into my summary statistics format
# PLINK qassoc signature:
# CHR          SNP         BP    NMISS       BETA         SE         R2        T            P
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
' OFS='\t' $testSet$'.bim' $bootStrapLoc$foldFam$'PLINK.qassoc' > $bootStrapLoc$foldFam$'_phe_sumstats'


# go through each chrom
for ((c=1; c<=22; c++)); do



#####################################
# LDPred2 Auto does not really use the test set, the 'test set' I used was for CD (it is only used to eliminate extreme PRS, without looking at the phenotypes)

arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_chrom.R '$baseLDpredLoc$' '$bootStrapLoc$foldFam$'_phe_sumstats '$c$' '$avgNumIndis$' '$LDPredOutLoc$' '$foldFam$outLdpredName$'_'$c$' '$bootStrapLoc$foldFam$'TEST_hapmap3 '$shaPRSscriptLoc$' '$binPheno



pheA_LDPredPRS=$LDPredOutLoc$foldFam$outLdpredName$'_'$c

# LDPRED BASELINES:
# if it is GWAS3_IBD, then both GWAS1_CD / GWAS2_UC
if [ $j == 3 ] ; then


if [ -s "$pheA_LDPredPRS" ] ; then
#echo "ALREADY EXISTS: "$pheA_LDPredPRS
b=1
else

echo "Combined PRS FOR GWAS1/GWAS2 from GWAS3 IBD "$i$"_chrom"$c
bsub -q normal -m "modern_hardware" -G team152 -n1 -J BASE_IBD${i}_${c} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $LDPredOutLoc$'BASE_IBD_'$i$'_'$c$'.out' -e $LDPredOutLoc$'BASE_IBD_'$i$'_'$c$'.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"
fi

fi


# if it is GWAS3_CD then just GWAS1_CD
if [ $j == 4 ] ; then
if [ -s "$pheA_LDPredPRS" ] ; then
#echo "ALREADY EXISTS: "$pheA_LDPredPRS
b=1
else

echo "PRS FOR GWAS1 from GWAS3 CD "$i$"_chrom"$c
bsub -q normal -m "modern_hardware" -G team152 -n1 -J BASE_CD_CD${i}_${c} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $LDPredOutLoc$'BASE_CD_CD'$i$'_'$c$'.out' -e $LDPredOutLoc$'BASE_CD_CD'$i$'_'$c$'.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"
fi

fi


# if it is GWAS3_UC then just GWAS2_UC
if [ $j == 5 ] ; then
if [ -s "$pheA_LDPredPRS" ] ; then
#echo "ALREADY EXISTS: "$pheA_LDPredPRS
b=1
else

echo "PRS FOR GWAS2 from GWAS3 UC "$i$"_chrom"$c
bsub -q normal -m "modern_hardware" -G team152 -n1 -J BASE_UC_UC${i}_${c} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $LDPredOutLoc$'BASE_UC_UC'$i$'_'$c$'.out' -e $LDPredOutLoc$'BASE_UC_UC'$i$'_'$c$'.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"
fi

fi

done # end of chroms

done # end of bootstraps
done # end of cohort loop

###################################################################


###############################################################################################################
# WAIT UNTIL LDPRED2 JOBS FINISH ON THE CLUSTER



j=4
i=1
# Aggregate LDpred2 results for Baselines
for (( j=4; j<${arraylength}+1; j++ )); do # GWAS3 IBD/CD/UC
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
LDPredOutLoc=$bootStrapLoc$'LDpred/'


mkdir -p $LDPredOutLoc

for ((i=1; i<=numBootStraps; i++)); do # numBootStraps go through all bootstraps
foldFam=$filenameStem$'_f'$i

pheA_LDPredPRS=$LDPredOutLoc$foldFam$outLdpredName

# 1) concat all chrom PRS into single files
cat $LDPredOutLoc$foldFam$outLdpredName$'_'* > $pheA_LDPredPRS
head $pheA_LDPredPRS
wc -l $pheA_LDPredPRS # 958305


# depending on which pheno it was, we pick the other's ldpred adjusted sumstats file, and pick a different test set
if [ $j == 4 ] ; then
cohortDataFolder_pheB=$ibdDataLoc${ibdarray_target_extended[4]}
testSet=$bootStrapLoc$foldFam$'TEST_hapmap3'
fi

if [ $j == 5 ] ; then
cohortDataFolder_pheB=$ibdDataLoc${ibdarray_target_extended[3]}
testSet=$bootStrapLoc$foldFam$'TEST_hapmap3_UC'
fi

# produce .profile score for the LDPRED2 betas, this is needed for SMTPred
# we need to produce 2 .profile score for each pheno, one trained on the pheno itself, and another on pheB, but still produced for the SAME target Test

# trained PheA -> predicted PheA
arguments=' --memory '$plinkMem$' --bfile '$testSet$' --score '$pheA_LDPredPRS$' sum --out '$pheA_LDPredPRS
$plink $arguments


# trained PheB -> predicted PheA
bootStrapLoc_pheB=$cohortDataFolder_pheB$'bootstrap/'
LDPredOutLoc_pheB=$bootStrapLoc_pheB$'LDpred/'
pheB_LDPredPRS=$LDPredOutLoc_pheB$foldFam$outLdpredName
arguments=' --memory '$plinkMem$' --bfile '$testSet$' --score '$pheB_LDPredPRS$' sum --out '$pheA_LDPredPRS$'_pheB'
$plink $arguments

done # end of bootstraps
done # end of cohort loop



# Baseline IBD:
j=3
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
LDPredOutLoc=$bootStrapLoc$'LDpred/'


mkdir -p $LDPredOutLoc

for ((i=1; i<=numBootStraps; i++)); do # numBootStraps go through all bootstraps
foldFam=$filenameStem$'_f'$i

pheA_LDPredPRS=$LDPredOutLoc$foldFam$outLdpredName

# 1) concat all chrom PRS into single files
#cat $LDPredOutLoc$foldFam$outLdpredName$'_'* > $pheA_LDPredPRS
cat $LDPredOutLoc$foldFam$outLdpredName$'_'[!a-z] $LDPredOutLoc$foldFam$outLdpredName$'_'[!a-z][!a-z] > $pheA_LDPredPRS
head $pheA_LDPredPRS
wc -l $pheA_LDPredPRS # 958305


# CD
arguments=' --memory '$plinkMem$' --bfile '$bootStrapLoc$foldFam$'TEST_hapmap3 --score '$pheA_LDPredPRS$' sum --out '$pheA_LDPredPRS$'_CD_IBD'
$plink $arguments


# UC
arguments=' --memory '$plinkMem$' --bfile '$bootStrapLoc$foldFam$'TEST_hapmap3_UC --score '$pheA_LDPredPRS$' sum --out '$pheA_LDPredPRS$'_UC_IBD'
$plink $arguments

done # end of bootstraps


###############################################################################################################



j=4
i=1
#rm -rf $diagOutLoc

# 4) SMTPRED
arraylength=${#ibdarray_extended[@]}
for (( j=4; j<${arraylength}+1; j++ )); do # GWAS3 IBD/CD/UC
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
LDPredOutLoc=$bootStrapLoc$'LDpred/'

mkdir -p $LDPredOutLoc

for ((i=1; i<=numBootStraps; i++)); do # go through all bootstraps to create association files from the training set
foldFam=$filenameStem$'_f'$i
cohortDataFolder_pheCombined=$ibdDataLoc${ibdarray_target_extended[2]}
bootStrapLoc_pheCombined=$cohortDataFolder_pheCombined$'bootstrap/'
LDPredOutLoc_pheCombined=$bootStrapLoc_pheCombined$'LDpred/'
pheCombined_LDPredPRS=$LDPredOutLoc_pheCombined$foldFam$outLdpredName

pheA_LDPredPRS=$LDPredOutLoc$foldFam$outLdpredName
pheA_sumstats=$bootStrapLoc$foldFam$'_phe_sumstats'


#####################################a

# if it is GWAS3_CD then just GWAS1_CD
if [ $j == 4 ] ; then
echo "PRS FOR GWAS1 from GWAS3 CD"
#PheB, is UC, if PheA was CD
cohortDataFolder_pheB=$ibdDataLoc${ibdarray_target_extended[4]}
bootStrapLoc_pheB=$cohortDataFolder_pheB$'bootstrap/'
LDPredOutLoc_pheB=$bootStrapLoc_pheB$'LDpred/'
pheB_sumstats=$bootStrapLoc_pheB$foldFam$'_phe_sumstats'

bsub -q normal -m "modern_hardware" -G team152 -n5 -J SMTPRED_CD${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $LDPredOutLoc$'SMTPRED_CD'$i$'.out' -e $LDPredOutLoc$'SMTPRED_CD'$i$'.err'  "perform_SMPTpred_LDPred2 $bootStrapLoc$foldFam $pheA_LDPredPRS $pheA_LDPredPRS$'_pheB' $pheA_sumstats $pheB_sumstats"

#perform_SMPTpred_LDPred2 $bootStrapLoc$foldFam $pheA_LDPredPRS $pheA_LDPredPRS$'_pheB' $pheA_sumstats $pheB_sumstats

fi


# if it is GWAS3_UC then just GWAS2_UC
if [ $j == 5 ] ; then
echo "PRS FOR GWAS2 from GWAS3 UC"
#PheB, is UC, if PheA was CD
cohortDataFolder_pheB=$ibdDataLoc${ibdarray_target_extended[3]}
bootStrapLoc_pheB=$cohortDataFolder_pheB$'bootstrap/'
LDPredOutLoc_pheB=$bootStrapLoc_pheB$'LDpred/'

pheB_sumstats=$bootStrapLoc_pheB$foldFam$'_phe_sumstats'

bsub -q normal -m "modern_hardware" -G team152 -n5 -J SMTPRED_UC${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $LDPredOutLoc$'SMTPRED_UC'$i$'.out' -e $LDPredOutLoc$'SMTPRED_UC'$i$'.err'  "perform_SMPTpred_LDPred2 $bootStrapLoc$foldFam $pheA_LDPredPRS $pheA_LDPredPRS$'_pheB' $pheA_sumstats $pheB_sumstats"
fi

done # end of bootstraps
done # end of cohort loop


###############################################################################################################



# 4) SHAPRS
arraylength=${#ibdarray_extended[@]}
for (( j=4; j<${arraylength}+1; j++ )); do # GWAS3 CD/UC
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
shaPRSPre=$bootStrapLoc$'shaPRSPre/'
binPheno='0' # we actually ran the association as quantiative phenotypes
#avgNumIndis=${N_effs[$j-1]}

mkdir -p $shaPRSPre

for ((i=1; i<=numBootStraps; i++)); do # go through all bootstraps to create association files from the training set
foldFam=$filenameStem$'_f'$i
pheA_LDPredPRS=$shaPRSPre$foldFam$outLdpredName
pheA_sumstats=$bootStrapLoc$foldFam$'_phe_sumstats'


#####################################


# I) Apply shaPRS pre-processing steps, all chroms at once
# if it is GWAS3_CD then just GWAS1_CD
if [ $j == 4 ] ; then
#echo "PRS FOR GWAS1 from GWAS3 CD"
#PheB, is UC, if PheA was CD
pheno='CD'
cohortDataFolder_pheB=$ibdDataLoc${ibdarray_target_extended[4]}
bootStrapLoc_pheB=$cohortDataFolder_pheB$'bootstrap/'
pheB_sumstats=$bootStrapLoc_pheB$foldFam$'_phe_sumstats'

#bsub -q normal -m "modern_hardware" -G team152 -n1 -J SHAPRS_CD${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $shaPRSPre$'SHAPRS_CD'$i$'.out' -e $shaPRSPre$'SHAPRS_CD'$i$'.err'  "shaPRS_LDPRED2 $rho $bootStrapLoc$foldFam $pheA_LDPredPRS $pheCombined_LDPredPRS $pheA_sumstats $pheB_sumstats"


fi



# if it is GWAS3_UC then just GWAS2_UC
if [ $j == 5 ] ; then
#echo "PRS FOR GWAS2 from GWAS3 UC"
#PheB, is UC, if PheA was CD
pheno='UC'
cohortDataFolder_pheB=$ibdDataLoc${ibdarray_target_extended[3]}
bootStrapLoc_pheB=$cohortDataFolder_pheB$'bootstrap/'
pheB_sumstats=$bootStrapLoc_pheB$foldFam$'_phe_sumstats'

#bsub -q normal -m "modern_hardware" -G team152 -n1 -J SHAPRS_UC${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $shaPRSPre$'SHAPRS_UC'$i$'.out' -e $shaPRSPre$'SHAPRS_UC'$i$'.err'  "shaPRS_LDPRED2 $rho $bootStrapLoc$foldFam $pheA_LDPredPRS $pheCombined_LDPredPRS $pheA_sumstats $pheB_sumstats"

fi




# 1) Apply the preprocessing shaPRS step ( already done)
shaPRS_new $pheA_sumstats $pheB_sumstats $shaPRSPre$pheno$'_'$i $rho '0'



# visualise diagnostics
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_qval_manhattan.R '$shaPRSPre$pheno$'_'$i$'_SE_meta '$shaPRSPre$pheno$'_'$i$'_lFDR_meta_SNP_lFDR '$shaPRSPre$pheno$'_'$i$' '$pheno$' '$shaPRSPre$pheno$'_'$i$'_sumstats_meta'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# 2) Perform QC
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_QC.R '$baseLDpredLoc$'/ '$shaPRSPre$pheno$'_'$i$'_sumstats_meta '$shaPRSPre$pheno$'_'$i$'_sumstats_meta_QC_keep '$binPheno
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# subset the SNPs to be the same panel for all 3 sumstats 
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $shaPRSPre$pheno$'_'$i$'_sumstats_meta_QC_keep' $shaPRSPre$pheno$'_'$i$'_sumstats_meta'  > $shaPRSPre$pheno$'_'$i$'_sumstats_meta_QC'


# go through each chrom
for ((c=1; c<=22; c++)); do


#####################################
# LDPred2 Auto does not really use the test set, the 'test set' I used was for CD (it is only used to eliminate extreme PRS, without looking at the phenotypes)
pheA_LDPredPRS=$LDPredOutLoc$foldFam$outLdpredName$'_'$c

# III) Build PRS via LDpred2
outfile=$shaPRSPre$pheno$'_'$i$'_'$c
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
echo MISSING ${pheno}_${i}_chrom${c}
avgNumIndis=$(tail -n 1 "$shaPRSPre$pheno"_"$i"_sumstats_meta_QC | awk '{ print $10}')
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_chrom.R '$baseLDpredLoc$' '$shaPRSPre$pheno$'_'$i$'_sumstats_meta_QC '$c$' '$avgNumIndis$' '$shaPRSPre$' '$pheno$'_'$i$'_'$c$' '$bootStrapLoc$foldFam$'TEST_hapmap3 '$shaPRSscriptLoc$' '$binPheno

#$Rscript $arguments

bsub -q normal -m "modern_hardware" -G team152 -n1 -J ${pheno}_${i}_${c} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $shaPRSPre$pheno$'_'$i$'_'$c$'.out' -e $shaPRSPre$pheno$'_'$i$'_'$c$'.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"
fi

done # end of chroms

done # end of bootstraps
done # end of cohort loop


###############################################################


# aggregate per chrom shaPRS+LDpred2 results:
for (( j=4; j<${arraylength}+1; j++ )); do # GWAS3 IBD/CD/UC
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
LDPredOutLoc=$bootStrapLoc$'LDpred/'
shaPRSPre=$bootStrapLoc$'shaPRSPre/'


for ((i=1; i<=numBootStraps; i++)); do # numBootStraps go through all bootstraps
foldFam=$filenameStem$'_f'$i

pheA_LDPredPRS=$LDPredOutLoc$foldFam$outLdpredName$'_shaprs'

# 1) concat all chrom PRS into single files
cat $LDPredOutLoc$foldFam$outLdpredName$'_'[!a-z] $LDPredOutLoc$foldFam$outLdpredName$'_'[!a-z][!a-z] > $pheA_LDPredPRS
head $pheA_LDPredPRS
wc -l $pheA_LDPredPRS # 958305


# depending on which pheno it was, we pick the other's ldpred adjusted sumstats file, and pick a different test set
if [ $j == 4 ] ; then
arguments=' --memory '$plinkMem$' --bfile '$bootStrapLoc$foldFam$'TEST_hapmap3 --score '$pheA_LDPredPRS$' sum --out '$pheA_LDPredPRS$'_shaprs_CD'
$plink $arguments

fi

if [ $j == 5 ] ; then
arguments=' --memory '$plinkMem$' --bfile '$bootStrapLoc$foldFam$'TEST_hapmap3_UC --score '$pheA_LDPredPRS$' sum --out '$pheA_LDPredPRS$'_shaprs_UC'
$plink $arguments

fi

done # end of bootstraps
done # end of cohort loop


###############################################################################

# Evaluate ShaPRS test sets and recreate final results...

# port code that starts from below: # shaPRS BLEND: # I accidentally included a header


# 2 subpheno results
results_GWAS2UC_UC_PRS=$ibdResultsLoc$'PRS3/GWAS2UC_UC_PRS/'
rm -rf $results_GWAS2UC_UC_PRS
mkdir -p $results_GWAS2UC_UC_PRS

results_GWAS1CD_CD_PRS=$ibdResultsLoc$'PRS3/GWAS1CD_CD_PRS/'
rm -rf $results_GWAS1CD_CD_PRS
mkdir -p $results_GWAS1CD_CD_PRS

# 2 combined results
results_GWAS2UC_IBD_PRS=$ibdResultsLoc$'PRS3/GWAS2UC_IBD_PRS/'
rm -rf $results_GWAS2UC_IBD_PRS
mkdir -p $results_GWAS2UC_IBD_PRS

results_GWAS1CD_IBD_PRS=$ibdResultsLoc$'PRS3/GWAS1CD_IBD_PRS/'
rm -rf $results_GWAS1CD_IBD_PRS
mkdir -p $results_GWAS1CD_IBD_PRS


# 2 ShaPRS Blend
results_GWAS2UC_Blend_PRS=$ibdResultsLoc$'PRS3/GWAS2UC_Blend_PRS/'
rm -rf $results_GWAS2UC_Blend_PRS
mkdir -p $results_GWAS2UC_Blend_PRS

results_GWAS1CD_Blend_PRS=$ibdResultsLoc$'PRS3/GWAS1CD_Blend_PRS/'
rm -rf $results_GWAS1CD_Blend_PRS
mkdir -p $results_GWAS1CD_Blend_PRS


# 2 SMTPRed
results_GWAS2UC_SMTPred_PRS=$ibdResultsLoc$'PRS3/GWAS2UC_SMTPred_PRS/'
rm -rf $results_GWAS2UC_SMTPred_PRS
mkdir -p $results_GWAS2UC_SMTPred_PRS

results_GWAS1CD_SMTPred_PRS=$ibdResultsLoc$'PRS3/GWAS1CD_SMTPred_PRS/'
rm -rf $results_GWAS1CD_SMTPred_PRS
mkdir -p $results_GWAS1CD_SMTPred_PRS





j=3
i=1
j=4


# 6) Produce Test results:
arraylength=${#ibdarray_extended[@]}
for (( j=3; j<${arraylength}+1; j++ )); do # GWAS3 IBD/CD/UC
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
LDPredOutLoc=$bootStrapLoc$'LDpred/'



for ((i=1; i<=numBootStraps; i++)); do # go through all bootstraps
foldFam=$filenameStem$'_f'$i

pheA_LDPredPRS=$LDPredOutLoc$foldFam$outLdpredName
pheA_sumstats=$bootStrapLoc$foldFam$'_phe_sumstats'

cohortDataFolder_pheCombined=$ibdDataLoc${ibdarray_target_extended[2]}
bootStrapLoc_pheCombined=$cohortDataFolder_pheCombined$'bootstrap/'
LDPredOutLoc_pheCombined=$bootStrapLoc_pheCombined$'LDpred/'
pheCombined_LDPredPRS=$LDPredOutLoc_pheCombined$foldFam$outLdpredName


#####################################aa
if [ $j == 3 ] ; then
echo "PRS FOR GWAS1/2 from Combined"

# Combined -> GWAS1/CD 
echo "CD Meta"
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'_CD_IBD.profile' FS=" " $bootStrapLoc$foldFam$'TEST_hapmap3.fam' > $pheA_LDPredPRS$'_CD_IBD.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_CD_IBD.csv '$results_GWAS1CD_IBD_PRS$'/_CD_IBD 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# Combined -> GWAS2/UC
echo "UC Meta"
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'_UC_IBD.profile' FS=" " $bootStrapLoc$foldFam$'TEST_hapmap3_UC.fam' > $pheA_LDPredPRS$'_UC_IBD.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_UC_IBD.csv '$results_GWAS2UC_IBD_PRS$'/_UC_IBD 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

fi

#####################################

# GWA
if [ $j == 4 ] ; then
# load up best threshold for composite

echo "CD SHAPRS BLEND"
testSet=$bootStrapLoc$foldFam$'TEST_hapmap3'

echo "CD ShaPRS"
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'_shaprs_CD.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_shaprs_CD.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_shaprs_CD.csv '$results_GWAS1CD_Blend_PRS$'/_shaprs_CD 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

echo "CD VANILLA"
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_CD.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_CD.csv '$results_GWAS1CD_CD_PRS$'/_CD 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

echo "CD SMTPRED"
awk 'FNR == NR { file1[ $1 ] = $3; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $bootStrapLoc$foldFam$'_SMTPRED.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_SMTPred.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_SMTPred.csv '$results_GWAS1CD_SMTPred_PRS$'/SMTPred 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
fi

# if it is GWAS3_UC then just GWAS2_UC
if [ $j == 5 ] ; then
# load up best threshold for composite

echo "UC SHAPRS BLEND"
testSet=$bootStrapLoc$foldFam$'TEST_hapmap3_UC'

echo "UC ShaPRS"
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'_shaprs_UC.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_shaprs_UC.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_shaprs_UC.csv '$results_GWAS2UC_Blend_PRS$'/_shaprs_UC 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

echo "UC VANILLA"
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_UC.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_UC.csv '$results_GWAS2UC_UC_PRS$'/_UC 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

echo "UC SMTPRED"
awk 'FNR == NR { file1[ $1 ] = $3; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $bootStrapLoc$foldFam$'_SMTPRED.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_SMTPred.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_SMTPred.csv '$results_GWAS2UC_SMTPred_PRS$'/SMTPred 0 0 1'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

fi

done # end of bootstraps
done # end of cohort loop



head /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/gwas3_uc/bootstrap/LDpred/GWAS_f20ldpred2_PRS_UC.csv

head /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/gwas3_uc/bootstrap/LDpred/GWAS_f19ldpred2_PRS_UC.csv


#########################
# Plots:
# need to rename the filenames as those will be used in the plot
mv $results_GWAS1CD_CD_PRS$'/_CD' $results_GWAS1CD_CD_PRS$'/CD_from_CD'
mv $results_GWAS1CD_IBD_PRS$'/_CD_IBD' $results_GWAS1CD_IBD_PRS$'/CD_from_IBD'
mv $results_GWAS1CD_Blend_PRS$'/_shaprs_CD' $results_GWAS1CD_Blend_PRS$'/CD_from_Blend'
mv $results_GWAS1CD_SMTPred_PRS$'/SMTPred' $results_GWAS1CD_SMTPred_PRS$'/CD_from_SMTPred'

mv $results_GWAS2UC_UC_PRS$'/_UC' $results_GWAS2UC_UC_PRS$'/UC_from_UC'
mv $results_GWAS2UC_IBD_PRS$'/_UC_IBD' $results_GWAS2UC_IBD_PRS$'/UC_from_IBD'
mv $results_GWAS2UC_Blend_PRS$'/_shaprs_UC' $results_GWAS2UC_Blend_PRS$'/UC_from_Blend'
mv $results_GWAS2UC_SMTPred_PRS$'/SMTPred' $results_GWAS2UC_SMTPred_PRS$'/UC_from_SMTPred'

cp $results_GWAS1CD_IBD_PRS$'/CD_from_IBD' $results_GWAS1CD_IBD_PRS$'/CD_from_Meta'
cp $results_GWAS1CD_Blend_PRS$'/CD_from_Blend' $results_GWAS1CD_Blend_PRS$'/CD_from_shaPRS'

cp $results_GWAS2UC_IBD_PRS$'/UC_from_IBD' $results_GWAS2UC_IBD_PRS$'/UC_from_Meta'
cp $results_GWAS2UC_Blend_PRS$'/UC_from_Blend' $results_GWAS2UC_Blend_PRS$'/UC_from_shaPRS'

 
# Single merged plot
arguments='/nfs/users/nfs_m/mk23/scripts/dotPlot_bootstrap_largeFont_merged.R '$ibdResultsLoc$'PRS3/ CD_UC_PRS_PUB CD_UC_prediction_accuracy -1 #c5b000 #166938 '$results_GWAS1CD_CD_PRS$'/CD_from_CD '$results_GWAS1CD_IBD_PRS$'/CD_from_Meta '$results_GWAS1CD_Blend_PRS$'/CD_from_shaPRS '$results_GWAS1CD_SMTPred_PRS$'/CD_from_SMTPred '$results_GWAS2UC_UC_PRS$'/UC_from_UC '$results_GWAS2UC_IBD_PRS$'/UC_from_Meta '$results_GWAS2UC_Blend_PRS$'/UC_from_shaPRS '$results_GWAS2UC_SMTPred_PRS$'/UC_from_SMTPred'
$Rscript $arguments

[1] "predicted: CD\ntrained: CD v predicted: CD\ntrained: Meta : 0.0955640413983717 / diff: -9 %"
[1] "predicted: CD\ntrained: CD: 0.0238835677157861"
[1] "predicted: CD\ntrained: CD v predicted: CD\ntrained: shaPRS : 9.14458482722662e-06 / diff: -23 %" X
[1] "predicted: CD\ntrained: CD: 0.0238835677157861"
[1] "predicted: CD\ntrained: CD v predicted: CD\ntrained: SMTPred : 0.00170369203956435 / diff: -5 %"
[1] "predicted: CD\ntrained: CD: 0.0238835677157861"
[1] "predicted: CD\ntrained: CD v predicted: UC\ntrained: UC : 5.08113530786309e-11 / diff: 68 %"
[1] "predicted: CD\ntrained: CD: 0.0238835677157861"
[1] "predicted: CD\ntrained: CD v predicted: UC\ntrained: Meta : 3.32216123763691e-08 / diff: 47 %"
[1] "predicted: CD\ntrained: CD: 0.0238835677157861"
[1] "predicted: CD\ntrained: CD v predicted: UC\ntrained: shaPRS : 1.35488175856169e-07 / diff: 40 %"
[1] "predicted: CD\ntrained: CD: 0.0238835677157861"
[1] "predicted: CD\ntrained: CD v predicted: UC\ntrained: SMTPred : 1.84109513183679e-10 / diff: 56 %"
[1] "predicted: CD\ntrained: CD: 0.0238835677157861"
[1] "predicted: CD\ntrained: Meta v predicted: CD\ntrained: shaPRS : 8.38188297233817e-08 / diff: -14 %" X
[1] "predicted: CD\ntrained: Meta: 0.0262347720766502"
[1] "predicted: CD\ntrained: Meta v predicted: CD\ntrained: SMTPred : 0.420942567877758 / diff: 4 %"
[1] "predicted: CD\ntrained: Meta: 0.0262347720766502"
[1] "predicted: CD\ntrained: Meta v predicted: UC\ntrained: UC : 1.52852715630786e-10 / diff: 76 %"
[1] "predicted: CD\ntrained: Meta: 0.0262347720766502"
[1] "predicted: CD\ntrained: Meta v predicted: UC\ntrained: Meta : 3.5961661243796e-09 / diff: 56 %"
[1] "predicted: CD\ntrained: Meta: 0.0262347720766502"
[1] "predicted: CD\ntrained: Meta v predicted: UC\ntrained: shaPRS : 1.89350043501853e-08 / diff: 49 %"
[1] "predicted: CD\ntrained: Meta: 0.0262347720766502"
[1] "predicted: CD\ntrained: Meta v predicted: UC\ntrained: SMTPred : 2.28933380403699e-10 / diff: 65 %"
[1] "predicted: CD\ntrained: Meta: 0.0262347720766502"
[1] "predicted: CD\ntrained: shaPRS v predicted: CD\ntrained: SMTPred : 0.000247451909977452 / diff: 18 %" X
[1] "predicted: CD\ntrained: shaPRS: 0.0300693790975887"
[1] "predicted: CD\ntrained: shaPRS v predicted: UC\ntrained: UC : 8.40442359741098e-14 / diff: 88 %"
[1] "predicted: CD\ntrained: shaPRS: 0.0300693790975887"
[1] "predicted: CD\ntrained: shaPRS v predicted: UC\ntrained: Meta : 5.18381044690828e-12 / diff: 68 %"
[1] "predicted: CD\ntrained: shaPRS: 0.0300693790975887"
[1] "predicted: CD\ntrained: shaPRS v predicted: UC\ntrained: shaPRS : 1.53642004890883e-11 / diff: 62 %"
[1] "predicted: CD\ntrained: shaPRS: 0.0300693790975887"
[1] "predicted: CD\ntrained: shaPRS v predicted: UC\ntrained: SMTPred : 6.08155385629614e-14 / diff: 77 %"
[1] "predicted: CD\ntrained: shaPRS: 0.0300693790975887"
[1] "predicted: CD\ntrained: SMTPred v predicted: UC\ntrained: UC : 1.10915641702109e-11 / diff: 73 %"
[1] "predicted: CD\ntrained: SMTPred: 0.0250802925311907"
[1] "predicted: CD\ntrained: SMTPred v predicted: UC\ntrained: Meta : 4.64577965253539e-09 / diff: 52 %"
[1] "predicted: CD\ntrained: SMTPred: 0.0250802925311907"
[1] "predicted: CD\ntrained: SMTPred v predicted: UC\ntrained: shaPRS : 1.23743713795488e-08 / diff: 45 %"
[1] "predicted: CD\ntrained: SMTPred: 0.0250802925311907"
[1] "predicted: CD\ntrained: SMTPred v predicted: UC\ntrained: SMTPred : 7.23831809505536e-11 / diff: 61 %"
[1] "predicted: CD\ntrained: SMTPred: 0.0250802925311907"

[1] "predicted: UC\ntrained: UC v predicted: UC\ntrained: Meta : 0.00231872679487593 / diff: -23 %"
[1] "predicted: UC\ntrained: UC: 0.0117192054306032"
[1] "predicted: UC\ntrained: UC v predicted: UC\ntrained: shaPRS : 3.72103299945217e-06 / diff: -30 %" X
[1] "predicted: UC\ntrained: UC: 0.0117192054306032"
[1] "predicted: UC\ntrained: UC v predicted: UC\ntrained: SMTPred : 0.000122361417174365 / diff: -13 %"
[1] "predicted: UC\ntrained: UC: 0.0117192054306032"
[1] "predicted: UC\ntrained: Meta v predicted: UC\ntrained: shaPRS : 0.00701273992769753 / diff: -7 %" X
[1] "predicted: UC\ntrained: Meta: 0.0148044492562028"
[1] "predicted: UC\ntrained: Meta v predicted: UC\ntrained: SMTPred : 0.059897566030372 / diff: 10 %"
[1] "predicted: UC\ntrained: Meta: 0.0148044492562028"
[1] "predicted: UC\ntrained: shaPRS v predicted: UC\ntrained: SMTPred : 0.000347721282708467 / diff: 17 %" X


#CD
[1] "predicted: CD\ntrained: CD v predicted: CD\ntrained: shaPRS : 9.14458482722662e-06 / diff: -23 %" X - cited in apper
[1] "predicted: CD\ntrained: Meta v predicted: CD\ntrained: shaPRS : 8.38188297233817e-08 / diff: -14 %" X
[1] "predicted: CD\ntrained: shaPRS v predicted: CD\ntrained: SMTPred : 0.000247451909977452 / diff: 18 %" X - cited in apper


# UC
[1] "predicted: UC\ntrained: UC v predicted: UC\ntrained: shaPRS : 3.72103299945217e-06 / diff: -30 %" X - cited in apper
[1] "predicted: UC\ntrained: Meta v predicted: UC\ntrained: shaPRS : 0.00701273992769753 / diff: -7 %" X
[1] "predicted: UC\ntrained: shaPRS v predicted: UC\ntrained: SMTPred : 0.000347721282708467 / diff: 17 %" X - cited in apper


[1] "predicted: CD\ntrained: CD: 0.0238835677157861"
[1] "predicted: CD\ntrained: Meta: 0.0262347720766502"
[1] "predicted: CD\ntrained: shaPRS: 0.0300693790975887"
[1] "predicted: CD\ntrained: SMTPred: 0.0250802925311907"

[1] "predicted: UC\ntrained: UC: 0.0117192054306032"
[1] "predicted: UC\ntrained: Meta: 0.0148044492562028"
[1] "predicted: UC\ntrained: shaPRS: 0.015859817080825"
[1] "predicted: UC\ntrained: SMTPred: 0.0133628818931277"



###################################
# REUSABLE FUNCTIONS

# Reusable SMTPred function
function perform_SMPTpred_LDPred2 {
outputDir=$1
pheA_LDPredPRS=$2
pheB_LDPredPRS=$3
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
cp $pheA_LDPredPRS$'.profile' $outputDir$'/'$filename_pheA$'.profile'
cp $pheB_LDPredPRS$'.profile' $outputDir$'/'$filename_pheB$'.profile'



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
# $outputDir$'/multi_trait.score'
cp $outputDir$'/multi_trait.score' $outputDir$'_SMTPRED.profile'
}
export -f perform_SMPTpred_LDPred2 # this makes local functions executable when bsubbed	


#$python2
#import pip
#pip.main(['install', "pandas", "--no-cache-dir"])
#python2 -m pip install bitarray
#python2 -m pip install scipy

# creates sumstats from two PLINK GWAS performing fixed effect meta analysis, then coordinates it then adjusts them to produce a PRS
function shaPRS_LDPRED2 { 
rho=$1
outputDir=$2
pheA_LDPredPRS=$3
pheCombined_LDPredPRS=$4
pheA_sumstats=$5
pheB_sumstats=$6


# To calculate the Q-vals and their lFDR, we use the raw summary stats
# Create input file for adjusting script:
#       SNP	CHR	BP	Beta_A	SE_A	Beta_B	SE_B
# rs4040617   1  779322 -0.0017630 0.008608 -0.010990 0.008592
#awk 'FNR == NR { file1[ $3 ] = $7"\t"$8"\t"$4; next; } 
#FNR <= NR { { 
#if (FNR == 1) {print "SNP\tCHR\tBP\tBeta_A\tSE_A\tA_effectAllele\tBeta_B\tSE_B\tB_effectAllele" } else {
#if ( $3 in file1) { print $3"\t"$1"\t"$2"\t"file1[$3]"\t"$7"\t"$8"\t"$4} }   }
#}
#'  $pheA_sumstats $pheB_sumstats > $outputDir$'_SE_meta'


awk 'FNR == NR { file1[ $3 ] = $7"\t"$8"\t"$4"\t"$5; next; } 
FNR <= NR { { 
if (FNR == 1) {print "SNP\tCHR\tBP\tBeta_A\tSE_A\tA1.x\tA2.x\tBeta_B\tSE_B\tA1.y\tA2.y" } else {
if ( $3 in file1) { print $3"\t"$1"\t"$2"\t"file1[$3]"\t"$7"\t"$8"\t"$4"\t"$5} }   }
}
' $pheA_sumstats $pheB_sumstats > $outputDir$'_SE_meta'



qthresholds=( 0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.43 0.39 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6 ) # the lFDR thresholds used to export SNPs
qthresholds_string=$(echo ${qthresholds[*]// /|})  # convert an array to string: https://stackoverflow.com/questions/13470413/converting-a-bash-array-into-a-delimited-string

# Export out the lFDR values, and the thresholds
# arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_adjust_wrapper.R '$outputDir$'_SE_meta '$rho$' '$outputDir$'_lFDR_meta '$qthresholds_string
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_adjust_wrapper.R '$outputDir$'_SE_meta '$outputDir$'_lFDR_meta /nfs/users/nfs_m/mk23/scripts/shaPRS.R'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# this outputs following:
# GWAS_f1_lFDR_meta_SNP_lFDR
# and one for each threshold like:
# GWAS_f1_lFDR_meta_SNPs_0.99

# we blend/switch between the post-LDPred2 adjusted betas, to produce the PRS we don't need anything else (eg SE)
# arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_PRS_processor.R '$pheA_LDPredPRS$' '$pheCombined_LDPredPRS$' '$outputDir$'_lFDR_meta_SNP_lFDR '$outputDir$'_shaprs '$outputDir$'_lFDR_meta_SNPs '$qthresholds_string
arguments='/nfs/users/nfs_m/mk23/scripts/sumstatsBlender_shaPRS_meta_wrapper.R '$pheA_sumstats$' '$pheB_sumstats$' '$outputDir$'_lFDR_meta_SNP_lFDR '$rho$' '$outputDir$'_shaprs /nfs/users/nfs_m/mk23/scripts/shaPRS.R'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# this produces 1 PRS for the blend, and 5 outputs for each hard threshold postfixed like $outputDir$'_shaprs_0.96' etc

}
export -f shaPRS_LDPRED2 # this makes local functions executable when bsubbed



# performs diagnostics that compares if cor(B1,B2) == cor(B1/SE1,B2/SE2) ?
function shaPRS_CorrDiag { 
rho=$1
outputDir=$2
pheA_sumstats=$3
pheB_sumstats=$4
diagOutLoc=$5
i=$6

# we blend/switch between the post-LDPred2 adjusted betas, to produce the PRS we don't need anything else (eg SE)
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_CorrDiag.R '$rho$' '$pheA_sumstats$' '$pheB_sumstats$' '$outputDir$'_lFDR_meta_SNP_lFDR '$diagOutLoc$' '$i
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

}
export -f shaPRS_CorrDiag # this makes local functions executable when bsubbed



# find the effective sample size for CD/UC/IBD
#GetN_eff_Fam /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/scratch/gwas3_cd/all_found.fam
#numCases 3810 / numcontrols: 9492
#10874


#GetN_eff_Fam /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/scratch/gwas3_uc/all_found.fam
#numCases 3765 / numcontrols: 9492
#10782

# For IBD
#cases=3810+3765 
#controls= 9492

#N_eff = 4 / (1 / cases + 1 / controls) 
#round(N_eff) # 16852

#fam="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/scratch/gwas3_uc/all_found.fam"
# calculates N_eff based on a fam file
function GetN_eff_Fam {
fam=$1 # stores the cases/controls counts

# 3) filter for QC failed in UKBB
awk '
FNR <= NR { 


if($6 == 1) controls++
else if ($6 == 2) cases++

 } END { 
 print "numCases "cases" / numcontrols: "controls
 # calculate N_eff
N_eff = 4 / (1 / cases + 1 / controls) 
print int(N_eff) 
} ' $fam
}


