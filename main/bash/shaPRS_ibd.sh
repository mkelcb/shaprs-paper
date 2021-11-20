####################################################
# dgx-server screens:

# screen -r -D 64096.pts-0.node-11-1-3

# screen -r -D 53278.pts-0.node-11-1-3

# Scripts for the IBD ShaPRS analyses: that builds on the QC from chapter2, 
# but adds aligning the alleles, SMTPred, and LDPred2 output

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

twasCommonLoc=$commonDataLoc$'twas/'
twasRawLoc=$twasCommonLoc$'raw/'
hapmap3_b37bim=$twasRawLoc$'hapmap3_r1_b37_fwd_consensus.qc.poly.recode.bim' # this is from here: https://www.biostars.org/p/306616/


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
#qthresholds=( 0.96 0.97 0.98 0.99 1.0 ) # the lFDR thresholds used to export SNPs
#qthresholds=( 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.0 )
#qthresholds=( 0.80 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.90 0.95 0.99 )

#qthresholds=( 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.80 0.90 0.95 0.99 )
#qthresholds=( 0.3 0.4 0.5 0.6 )
#qthresholds=( 0 0.1 0.2 )

#qthresholds=( 0.5 0.6 0.70 0.80 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.0 )
qthresholds=( 0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.43 0.39 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6 )

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



source /software/hgi/installs/anaconda3/etc/profile.d/conda.sh
conda activate mk23_pytorch

mkdir -p $diagFolder



####################################################

# I) Baseline LDPred2
###################################


############################################################################
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

# 2) Adjust the SNP betas in LDPred
arraylength=${#ibdarray_extended[@]}
for (( j=3; j<${arraylength}+1; j++ )); do # GWAS3 IBD/CD/UC
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
LDPredOutLoc=$bootStrapLoc$'LDpred/'


mkdir -p $LDPredOutLoc


for ((i=1; i<=numBootStraps; i++)); do # go through all bootstraps to create association files from the training set
foldFam=$filenameStem$'_f'$i

# create PLINK assocs: this is already done and does not need repeating: $bootStrapLoc$foldFam$'PLINK.qassoc'
# arguments=' --memory '$plinkMem$' --threads 5 --bfile '$bootStrapLoc$foldFam$' --assoc --out '$bootStrapLoc$foldFam$'PLINK'
# $plink $arguments


#####################################
# LDPred2 Auto does not really use the test set, the 'test set' I used was for CD (it is only used to eliminate extreme PRS, without looking at the phenotypes)
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_v2.R '$baseLDpredLoc$' '$bootStrapLoc$foldFam$'PLINK.qassoc '$NCORES$' '$avgNumIndis$' '$LDPredOutLoc$' '$foldFam$outLdpredName$' '$bootStrapLoc$foldFam$'TEST_hapmap3.bed'


# LDPRED BASELINES:
# if it is GWAS3_IBD, then both GWAS1_CD / GWAS2_UC
if [ $j == 3 ] ; then
echo "Combined PRS FOR GWAS1/GWAS2 from GWAS3 IBD"

# Want a combined LDpred2, as otherwise I would have to re-run LDpred2 6x times (1x for the blended and 5x for each composite, if I went the meta analysis way)
# so I will blend/swap between the post-LDPred2 processed sumstats, but for that I need the subphenos as well as the Combined estimates at hand

bsub -q long -m "modern_hardware" -G team152 -n${NCORES} -J BASE_IBD_CD -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $LDPredOutLoc$'BASE_IBD_CD.out' -e $LDPredOutLoc$'BASE_IBD_CD.err'  "$Rscript $arguments"
fi


# if it is GWAS3_CD then just GWAS1_CD
if [ $j == 4 ] ; then
echo "PRS FOR GWAS1 from GWAS3 CD"
bsub -q long -m "modern_hardware" -G team152 -n${NCORES} -J BASE_CD_CD${i} -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $LDPredOutLoc$'BASE_CD_CD'$i$'.out' -e $LDPredOutLoc$'BASE_CD_CD'$i$'.err'  "$Rscript $arguments"
fi


# if it is GWAS3_UC then just GWAS2_UC
if [ $j == 5 ] ; then
echo "PRS FOR GWAS2 from GWAS3 UC"
bsub -q long -m "modern_hardware" -G team152 -n${NCORES} -J BASE_UC_UC${i} -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $LDPredOutLoc$'BASE_UC_UC'$i$'.out' -e $LDPredOutLoc$'BASE_UC_UC'$i$'.err'  "$Rscript $arguments"
fi

done # end of bootstraps''
done # end of cohort loop




###############################################################################################################
# WAIT UNTIL LDPRED2 JOBS FINISH ON THE CLUSTER
# 3) produce LDpred2 .profile files for each test set: pattern is: $pheA_LDPredPRS'.profile'
j=5
i=1

for (( j=4; j<${arraylength}+1; j++ )); do # GWAS3 IBD/CD/UC
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
LDPredOutLoc=$bootStrapLoc$'LDpred/'


mkdir -p $LDPredOutLoc

for ((i=1; i<=numBootStraps; i++)); do # numBootStraps go through all bootstraps
foldFam=$filenameStem$'_f'$i

pheA_LDPredPRS=$LDPredOutLoc$foldFam$outLdpredName


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



###############################################################################################################
j=4
i=1
rm -rf $diagOutLoc

# 4) SMTPRED + SHAPRS
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


#####################################

# if it is GWAS3_CD then just GWAS1_CD
if [ $j == 4 ] ; then
echo "PRS FOR GWAS1 from GWAS3 CD"
#PheB, is UC, if PheA was CD
cohortDataFolder_pheB=$ibdDataLoc${ibdarray_target_extended[4]}
bootStrapLoc_pheB=$cohortDataFolder_pheB$'bootstrap/'
LDPredOutLoc_pheB=$bootStrapLoc_pheB$'LDpred/'
pheB_sumstats=$bootStrapLoc_pheB$foldFam$'_phe_sumstats'

#bsub -q normal -m "modern_hardware" -G team152 -n5 -J SHAPRS_CD${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $LDPredOutLoc$'SHAPRS_CD'$i$'.out' -e $LDPredOutLoc$'SHAPRS_CD'$i$'.err'  "shaPRS_LDPRED2 $rho $bootStrapLoc$foldFam $pheA_LDPredPRS $pheCombined_LDPredPRS $pheA_sumstats $pheB_sumstats"
#bsub -q normal -m "modern_hardware" -G team152 -n5 -J SMTPRED_CD${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $LDPredOutLoc$'SMTPRED_CD'$i$'.out' -e $LDPredOutLoc$'SMTPRED_CD'$i$'.err'  "perform_SMPTpred_LDPred2 $bootStrapLoc$foldFam $pheA_LDPredPRS $pheA_LDPredPRS$'_pheB' $pheA_sumstats $pheB_sumstats"


shaPRS_CorrDiag $rho $bootStrapLoc$foldFam $pheA_sumstats $pheB_sumstats $diagOutLoc $i

fi


# if it is GWAS3_UC then just GWAS2_UC
if [ $j == 5 ] ; then
echo "PRS FOR GWAS2 from GWAS3 UC"
#PheB, is UC, if PheA was CD
cohortDataFolder_pheB=$ibdDataLoc${ibdarray_target_extended[3]}
bootStrapLoc_pheB=$cohortDataFolder_pheB$'bootstrap/'
LDPredOutLoc_pheB=$bootStrapLoc_pheB$'LDpred/'

pheB_sumstats=$bootStrapLoc_pheB$foldFam$'_phe_sumstats'

#bsub -q normal -m "modern_hardware" -G team152 -n5 -J SHAPRS_UC${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $LDPredOutLoc$'SHAPRS_UC'$i$'.out' -e $LDPredOutLoc$'SHAPRS_UC'$i$'.err'  "shaPRS_LDPRED2 $rho $bootStrapLoc$foldFam $pheA_LDPredPRS $pheCombined_LDPredPRS $pheA_sumstats $pheB_sumstats"
#bsub -q normal -m "modern_hardware" -G team152 -n5 -J SMTPRED_UC${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $LDPredOutLoc$'SMTPRED_UC'$i$'.out' -e $LDPredOutLoc$'SMTPRED_UC'$i$'.err'  "perform_SMPTpred_LDPred2 $bootStrapLoc$foldFam $pheA_LDPredPRS $pheA_LDPredPRS$'_pheB' $pheA_sumstats $pheB_sumstats"
fi

done # end of bootstraps
done # end of cohort loop


###############################################################################################################


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



#####################################

j=4
i=1
z=3

j=5


# 5) find best performing PRS for composite shaPRS via the validation sets
arraylength=${#ibdarray_extended[@]}
for (( j=4; j<${arraylength}+1; j++ )); do # GWAS3 IBD/CD/UC
outputloc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[$j-1]}
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
LDPredOutLoc=$bootStrapLoc$'LDpred/'

mkdir -p $LDPredOutLoc


# delete any previous results
for (( z=0; z<${num_thresholds}; z++ )); do 
currentThreshold=${qthresholds[$z]}
if [ $currentThreshold == "1.0" ] ; then
currentThreshold=1
fi
if [ $currentThreshold == "0.90" ] ; then
currentThreshold=0.9
fi

if [ $currentThreshold == "0.80" ] ; then
currentThreshold=0.8
fi

if [ $currentThreshold == "0.70" ] ; then
currentThreshold=0.7
fi
rm -rf $bootStrapLoc$'/shaPRS_composite_'$currentThreshold
done # end of qthresholds loop


for ((i=1; i<=numBootStraps; i++)); do # go through all bootstraps to create association files from the training set
foldFam=$filenameStem$'_f'$i

pheA_LDPredPRS=$LDPredOutLoc$foldFam$outLdpredName
pheA_sumstats=$bootStrapLoc$foldFam$'_phe_sumstats'

cohortDataFolder_pheCombined=$ibdDataLoc${ibdarray_target_extended[2]}
bootStrapLoc_pheCombined=$cohortDataFolder_pheCombined$'bootstrap/'
LDPredOutLoc_pheCombined=$bootStrapLoc_pheCombined$'LDpred/'
pheCombined_LDPredPRS=$LDPredOutLoc_pheCombined$foldFam$outLdpredName


for (( z=0; z<${num_thresholds}; z++ )); do # go through all Q-val thresholds
currentThreshold=${qthresholds[$z]}


# problem, when writing the results for 1.0 and 0.90, R rounded it to '1', cutting the '.0' bit off, so the filename only has '1'
if [ $currentThreshold == "1.0" ] ; then
currentThreshold=1
fi
if [ $currentThreshold == "0.90" ] ; then
currentThreshold=0.9
fi

if [ $currentThreshold == "0.80" ] ; then
currentThreshold=0.8
fi

if [ $currentThreshold == "0.70" ] ; then
currentThreshold=0.7
fi


# grab the composite PRS file
currentThreshold_PRS=$bootStrapLoc$foldFam$'_shaprs_comp_'$currentThreshold
head $currentThreshold_PRS

# PROBLEM: the composite score had an r^2 too high between validation and bootstrap, that increased with lower thresholds
# subpheno-subpheno, not affected, only shaPRS, this indicated some sort of overfitting
# IDEA1: incorrect bootstrap sampling, that would include validation indis in the bootstrap: NO, this wasn't the case 
# IDEA2: the combined SNP estimates come from DIFFERENT RANDOM sampled bootstraps, those that may include the subpheno validation set indis too
# IE blending towards the combined means more overfitting
# to prevent this, need to cross filter the subpheno validation sets, to exclude individuals that were in the same bootstrap for the combined sample
awk 'FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {  
 if( $2 in file1 ); else  {print $0}
}
' $bootStrapLoc_pheCombined$foldFam$'.fam' $bootStrapLoc$foldFam$'_valid.fam' > $bootStrapLoc$foldFam$'_valid_noCombined.fam'



# build a .profile score for each on the VALIDATION set (both UC and CD)
arguments=' --memory '$plinkMem$' --bfile '$bootStrapLoc$foldFam$'_valid --keep '$bootStrapLoc$foldFam$'_valid_noCombined.fam --score '$currentThreshold_PRS$' sum --out '$currentThreshold_PRS
$plink $arguments

# get correlation between VALID set phenos and PRS
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $currentThreshold_PRS$'.profile' FS="\t" $bootStrapLoc$foldFam$'_valid.pheno' > $bootStrapLoc$foldFam$'_valid.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$bootStrapLoc$foldFam$'_valid.csv '$bootStrapLoc$'/shaPRS_composite_'$currentThreshold
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments




done # end of qthresholds loop

done # end of bootstraps
done # end of cohort loop




# pick the threshold with the highest r^2
j=4
i=1

arraylength=${#ibdarray_extended[@]}
for (( j=4; j<${arraylength}+1; j++ )); do # GWAS3 IBD/CD/UC
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[$j-1]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'

# find best threshold with the highest avg performance: output to $bootStrapLoc$'/shaPRS_composite_best'
arguments='/nfs/users/nfs_m/mk23/scripts/bestThreshold.R '$bootStrapLoc$'/shaPRS_composite_ '$qthresholds_string
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

done # end of cohort loop


# CD
#  "best threshold : 0.56, had highest mean 0.036861398 written to:/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/gwas3_cd/bootstrap//shaPRS_composite_best"

# UC 
#  "best threshold : 0.49, had highest mean 0.025556341 written to:/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/gwas3_uc/bootstrap//shaPRS_composite_best"



####################################################################


# shaPRS BLEND: # I accidentally included a header
 $bootStrapLoc$foldFam$'_shaprs'
# subpheno_CombinedPheno_blending.V1      subpheno_CombinedPheno_blending.V2.x    blendedBeta
# rs4040617       G       -0.000351439680459075


# shaPRS Composite:
bestQThreshold=$(head $bootStrapLoc$'/shaPRS_composite_best')
$bootStrapLoc$foldFam$'_shaprs_comp_'$bestQThreshold
# 4_96394332      G       -5.72575919590387e-06




# 2 subpheno results
results_GWAS2UC_UC_PRS=$ibdResultsLoc$'PRS2/GWAS2UC_UC_PRS/'
rm -rf $results_GWAS2UC_UC_PRS
mkdir -p $results_GWAS2UC_UC_PRS

results_GWAS1CD_CD_PRS=$ibdResultsLoc$'PRS2/GWAS1CD_CD_PRS/'
rm -rf $results_GWAS1CD_CD_PRS
mkdir -p $results_GWAS1CD_CD_PRS

# 2 combined results
results_GWAS2UC_IBD_PRS=$ibdResultsLoc$'PRS2/GWAS2UC_IBD_PRS/'
rm -rf $results_GWAS2UC_IBD_PRS
mkdir -p $results_GWAS2UC_IBD_PRS

results_GWAS1CD_IBD_PRS=$ibdResultsLoc$'PRS2/GWAS1CD_IBD_PRS/'
rm -rf $results_GWAS1CD_IBD_PRS
mkdir -p $results_GWAS1CD_IBD_PRS


# 2 ShaPRS Blend
results_GWAS2UC_Blend_PRS=$ibdResultsLoc$'PRS2/GWAS2UC_Blend_PRS/'
rm -rf $results_GWAS2UC_Blend_PRS
mkdir -p $results_GWAS2UC_Blend_PRS

results_GWAS1CD_Blend_PRS=$ibdResultsLoc$'PRS2/GWAS1CD_Blend_PRS/'
rm -rf $results_GWAS1CD_Blend_PRS
mkdir -p $results_GWAS1CD_Blend_PRS


# 2 ShaPRS Comp
results_GWAS2UC_Comp_PRS=$ibdResultsLoc$'PRS2/GWAS2UC_Comp_PRS/'
rm -rf $results_GWAS2UC_Comp_PRS
mkdir -p $results_GWAS2UC_Comp_PRS

results_GWAS1CD_Comp_PRS=$ibdResultsLoc$'PRS2/GWAS1CD_Comp_PRS/'
rm -rf $results_GWAS1CD_Comp_PRS
mkdir -p $results_GWAS1CD_Comp_PRS


# 2 SMTPRed
results_GWAS2UC_SMTPred_PRS=$ibdResultsLoc$'PRS2/GWAS2UC_SMTPred_PRS/'
rm -rf $results_GWAS2UC_SMTPred_PRS
mkdir -p $results_GWAS2UC_SMTPred_PRS

results_GWAS1CD_SMTPred_PRS=$ibdResultsLoc$'PRS2/GWAS1CD_SMTPred_PRS/'
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


pheA_LDPredPRS=$LDPredOutLoc$foldFam$outLdpredName
pheA_sumstats=$bootStrapLoc$foldFam$'_phe_sumstats'

cohortDataFolder_pheCombined=$ibdDataLoc${ibdarray_target_extended[2]}
bootStrapLoc_pheCombined=$cohortDataFolder_pheCombined$'bootstrap/'
LDPredOutLoc_pheCombined=$bootStrapLoc_pheCombined$'LDpred/'
pheCombined_LDPredPRS=$LDPredOutLoc_pheCombined$foldFam$outLdpredName


#####################################aa
if [ $j == 3 ] ; then
echo "PRS FOR GWAS1/2 from Combined"
cohortDataFolder_CD=$ibdDataLoc${ibdarray_target_extended[3]}
bootStrapLoc_CD=$cohortDataFolder_CD$'bootstrap/'
testSet=$bootStrapLoc$foldFam$'TEST_hapmap3'

# Combined -> GWAS1/CD 
arguments=' --memory '$plinkMem$' --bfile '$testSet$' --score '$pheA_LDPredPRS$' sum --out '$pheA_LDPredPRS$'_CD_IBD'
$plink $arguments
# find correlation 
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'_CD_IBD.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_CD_IBD.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_CD_IBD.csv '$results_GWAS1CD_IBD_PRS$'/_CD_IBD'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


cohortDataFolder_UC=$ibdDataLoc${ibdarray_target_extended[4]}
bootStrapLoc_UC=$cohortDataFolder_UC$'bootstrap/'
testSet=$bootStrapLoc$foldFam$'TEST_hapmap3_UC'
# Combined -> GWAS2/UC
arguments=' --memory '$plinkMem$' --bfile '$testSet$' --score '$pheA_LDPredPRS$' sum --out '$pheA_LDPredPRS$'_UC_IBD'
$plink $arguments
# find correlation 
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'_UC_IBD.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_UC_IBD.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_UC_IBD.csv '$results_GWAS2UC_IBD_PRS$'/_UC_IBD'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
fi

#####################################

# GWA
if [ $j == 4 ] ; then
# load up best threshold for composite
bestQThreshold=$(head $bootStrapLoc$'/shaPRS_composite_best')

echo "CD SHAPRS BLEND"
testSet=$bootStrapLoc$foldFam$'TEST_hapmap3'
arguments=' --memory '$plinkMem$' --bfile '$testSet$' --score '$bootStrapLoc$foldFam$'_shaprs sum --out '$pheA_LDPredPRS$'_shaprs_CD'
$plink $arguments
# find correlation 
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'_shaprs_CD.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_shaprs_CD.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_shaprs_CD.csv '$results_GWAS1CD_Blend_PRS$'/_shaprs_CD'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

echo "CD SHAPRS COMPOSITE"
arguments=' --memory '$plinkMem$' --bfile '$testSet$' --score '$bootStrapLoc$foldFam$'_shaprs_comp_'$bestQThreshold$' sum --out '$pheA_LDPredPRS$'_shaprs_comp_CD'
$plink $arguments
# find correlation 
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'_shaprs_comp_CD.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_shaprs_comp_CD.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_shaprs_comp_CD.csv '$results_GWAS1CD_Comp_PRS$'/_shaprs_comp_CD'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

echo "CD VANILLA"
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_CD.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_CD.csv '$results_GWAS1CD_CD_PRS$'/_CD'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

echo "CD SMTPRED"
awk 'FNR == NR { file1[ $1 ] = $3; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $bootStrapLoc$foldFam$'_SMTPRED.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_SMTPred.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_SMTPred.csv '$results_GWAS1CD_SMTPred_PRS$'/SMTPred'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
fi

# if it is GWAS3_UC then just GWAS2_UC
if [ $j == 5 ] ; then
# load up best threshold for composite
bestQThreshold=$(head $bootStrapLoc$'/shaPRS_composite_best')

echo "UC SHAPRS BLEND"
testSet=$bootStrapLoc$foldFam$'TEST_hapmap3_UC'
arguments=' --memory '$plinkMem$' --bfile '$testSet$' --score '$bootStrapLoc$foldFam$'_shaprs sum --out '$pheA_LDPredPRS$'_shaprs_UC'
$plink $arguments
# find correlation 
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'_shaprs_UC.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_shaprs_UC.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_shaprs_UC.csv '$results_GWAS2UC_Blend_PRS$'/_shaprs_UC'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

echo "UC SHAPRS COMPOSITE"
arguments=' --memory '$plinkMem$' --bfile '$testSet$' --score '$bootStrapLoc$foldFam$'_shaprs_comp_'$bestQThreshold$' sum --out '$pheA_LDPredPRS$'_shaprs_comp_UC'
$plink $arguments
# find correlation 
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'_shaprs_comp_UC.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_shaprs_comp_UC.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_shaprs_comp_UC.csv '$results_GWAS2UC_Comp_PRS$'/_shaprs_comp_UC'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

echo "UC VANILLA"
awk 'FNR == NR { file1[ $2 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $pheA_LDPredPRS$'.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_UC.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_UC.csv '$results_GWAS2UC_UC_PRS$'/_UC'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

echo "UC SMTPRED"
awk 'FNR == NR { file1[ $1 ] = $3; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $bootStrapLoc$foldFam$'_SMTPRED.profile' FS=" " $testSet$'.fam' > $pheA_LDPredPRS$'_SMTPred.csv'
arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$pheA_LDPredPRS$'_SMTPred.csv '$results_GWAS2UC_SMTPred_PRS$'/SMTPred'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

fi

done # end of bootstraps
done # end of cohort loop



########################################

# need to rename the filenames as those will be used in the plot
mv $results_GWAS1CD_CD_PRS$'/_CD' $results_GWAS1CD_CD_PRS$'/CD_from_CD'
mv $results_GWAS1CD_IBD_PRS$'/_CD_IBD' $results_GWAS1CD_IBD_PRS$'/CD_from_IBD'
mv $results_GWAS1CD_Blend_PRS$'/_shaprs_CD' $results_GWAS1CD_Blend_PRS$'/CD_from_Blend'
mv $results_GWAS1CD_Comp_PRS$'/_shaprs_comp_CD' $results_GWAS1CD_Comp_PRS$'/CD_from_Comp_0.56'
mv $results_GWAS1CD_SMTPred_PRS$'/SMTPred' $results_GWAS1CD_SMTPred_PRS$'/CD_SMTPred'

mv $results_GWAS2UC_UC_PRS$'/_UC' $results_GWAS2UC_UC_PRS$'/UC_from_UC'
mv $results_GWAS2UC_IBD_PRS$'/_UC_IBD' $results_GWAS2UC_IBD_PRS$'/UC_from_IBD'
mv $results_GWAS2UC_Blend_PRS$'/_shaprs_UC' $results_GWAS2UC_Blend_PRS$'/UC_from_Blend'
mv $results_GWAS2UC_Comp_PRS$'/_shaprs_comp_UC' $results_GWAS2UC_Comp_PRS$'/UC_from_Comp_0.49'
mv $results_GWAS2UC_SMTPred_PRS$'/SMTPred' $results_GWAS2UC_SMTPred_PRS$'/UC_SMTPred'


mv $results_GWAS1CD_SMTPred_PRS$'/CD_SMTPred' $results_GWAS1CD_SMTPred_PRS$'/CD_from_SMTPred'
mv $results_GWAS2UC_SMTPred_PRS$'/UC_SMTPred' $results_GWAS2UC_SMTPred_PRS$'/UC_from_SMTPred'


#Main plots
# produce plot CD, only include best hard threshold
arguments='/nfs/users/nfs_m/mk23/scripts/dotPlot_bootstrap_largeFont.R '$ibdResultsLoc$'PRS2/ CD_PRS_PUB CD_prediction_accuracy -1 #c5b000 '$results_GWAS1CD_CD_PRS$'/CD_from_CD '$results_GWAS1CD_IBD_PRS$'/CD_from_IBD '$results_GWAS1CD_Blend_PRS$'/CD_from_Blend '$results_GWAS1CD_Comp_PRS$'/CD_from_Comp_0.56 '$results_GWAS1CD_SMTPred_PRS$'/CD_from_SMTPred'
$Rscript $arguments

arguments='/nfs/users/nfs_m/mk23/scripts/dotPlot_bootstrap_largeFont.R '$ibdResultsLoc$'PRS2/ UC_PRS_PUB UC_prediction_accuracy -1 #166938 '$results_GWAS2UC_UC_PRS$'/UC_from_UC '$results_GWAS2UC_IBD_PRS$'/UC_from_IBD '$results_GWAS2UC_Blend_PRS$'/UC_from_Blend '$results_GWAS2UC_Comp_PRS$'/UC_from_Comp_0.49 '$results_GWAS2UC_SMTPred_PRS$'/UC_from_SMTPred'
$Rscript $arguments



# plots that only include Blend a
cp $results_GWAS1CD_IBD_PRS$'/CD_from_IBD' $results_GWAS1CD_IBD_PRS$'/CD_from_Meta'
cp $results_GWAS1CD_Blend_PRS$'/CD_from_Blend' $results_GWAS1CD_Blend_PRS$'/CD_from_shaPRS'

cp $results_GWAS2UC_IBD_PRS$'/UC_from_IBD' $results_GWAS2UC_IBD_PRS$'/UC_from_Meta'
cp $results_GWAS2UC_Blend_PRS$'/UC_from_Blend' $results_GWAS2UC_Blend_PRS$'/UC_from_shaPRS'

arguments='/nfs/users/nfs_m/mk23/scripts/dotPlot_bootstrap_largeFont.R '$ibdResultsLoc$'PRS2/ CD_PRS_PUB CD_prediction_accuracy -1 #c5b000 '$results_GWAS1CD_CD_PRS$'/CD_from_CD '$results_GWAS1CD_IBD_PRS$'/CD_from_Meta '$results_GWAS1CD_Blend_PRS$'/CD_from_shaPRS '$results_GWAS1CD_SMTPred_PRS$'/CD_from_SMTPred'
$Rscript $arguments

arguments='/nfs/users/nfs_m/mk23/scripts/dotPlot_bootstrap_largeFont.R '$ibdResultsLoc$'PRS2/ UC_PRS_PUB UC_prediction_accuracy -1 #166938 '$results_GWAS2UC_UC_PRS$'/UC_from_UC '$results_GWAS2UC_IBD_PRS$'/UC_from_Meta '$results_GWAS2UC_Blend_PRS$'/UC_from_shaPRS '$results_GWAS2UC_SMTPred_PRS$'/UC_from_SMTPred'
$Rscript $arguments



# Single merged plot
arguments='/nfs/users/nfs_m/mk23/scripts/dotPlot_bootstrap_largeFont_merged.R '$ibdResultsLoc$'PRS2/ CD_UC_PRS_PUB CD_UC_prediction_accuracy -1 #c5b000 #166938 '$results_GWAS1CD_CD_PRS$'/CD_from_CD '$results_GWAS1CD_IBD_PRS$'/CD_from_Meta '$results_GWAS1CD_Blend_PRS$'/CD_from_shaPRS '$results_GWAS1CD_SMTPred_PRS$'/CD_from_SMTPred '$results_GWAS2UC_UC_PRS$'/UC_from_UC '$results_GWAS2UC_IBD_PRS$'/UC_from_Meta '$results_GWAS2UC_Blend_PRS$'/UC_from_shaPRS '$results_GWAS2UC_SMTPred_PRS$'/UC_from_SMTPred'
$Rscript $arguments


# Q-value manhattan plot of the 1st CD bootstrapa
i=1
cohortDataFolder=$ibdDataLoc${ibdarray_target_extended[3]}
bootStrapLoc=$cohortDataFolder$'bootstrap/'
foldFam=$filenameStem$'_f'$i
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_qval_manhattan.R '$bootStrapLoc$foldFam$'_SE_meta '$bootStrapLoc$foldFam$'_lFDR_meta_SNP_lFDR '$ibdResultsLoc$'PRS2/Qvals'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# because we did shaPRS as a post-processing step for IBD, we dont have '$bootStrapLoc$foldFam$'_sumstats_meta'
# we used the post-ldpred2 adjusted betas in shaPRS, which means we don't have Standard errors, only Betas
# but we dont need SEs for the diagnostic plots, only the Betas
# so create dummy 'shaPFRS' sumstats file that uses the correct post shaPRS/ldPred adjusted betas


# post shaprs/ldpred adjusted betas
$bootStrapLoc$foldFam$'_shaprs'
$bootStrapLoc$foldFam$'_SE_meta'

awk 'FNR == NR { file1[ $1 ] = $3; next; } 
FNR <= NR {  
if (FNR == 1) {
$7=$7
$8="X"
print $0
}
else {
if ($3 in file1) {
$7=file1[$3]
print $0}
} }' $bootStrapLoc$foldFam$'_shaprs' $bootStrapLoc$foldFam$'_phe_sumstats' > $bootStrapLoc$foldFam$'_sumstats_meta'
head $bootStrapLoc$foldFam$'_sumstats_meta'
wc -l $bootStrapLoc$foldFam$'_sumstats_meta'


# same as above, but with additional diagnostic plots
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_qval_manhattan.R '$bootStrapLoc$foldFam$'_SE_meta '$bootStrapLoc$foldFam$'_lFDR_meta_SNP_lFDR '$ibdResultsLoc$'PRS2/Qvals_new CD '$bootStrapLoc$foldFam$'_sumstats_meta'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# lFDR an Qvalue
head $bootStrapLoc$foldFam$'_lFDR_meta_SNP_lFDR'

head $bootStrapLoc$foldFam$'_SE_meta' # the input data for pheno A and pheno B
head $bootStrapLoc$foldFam$'_phe_sumstats' # this is phe A, NOT the shaPRS!!! 1124114

# need to find phe

head $bootStrapLoc$foldFam$'_sumstats_meta'

$bootStrapLoc$foldFam$'_shaprs'
$bootStrapLoc$foldFam$'_lFDR_meta_SNPs'

/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/results/PRS2/Qvals_new_diagDat




#CD
[1] "predicted: CD\ntrained: CD v predicted: CD\ntrained: IBD : 0.560412835784633 / diff: -3 %"
[1] "predicted: CD\ntrained: CD v predicted: CD\ntrained: Blend : 9.14459459998199e-06 / diff: -23 %"  X
[1] "predicted: CD\ntrained: CD v predicted: CD\ntrained: Composite : 0.00258593069746124 / diff: -15 %"
[1] "predicted: CD\ntrained: CD v predicted: CD\ntrained: SMTPred : 0.000176704385866349 / diff: -6 %"
[1] "predicted: CD\ntrained: IBD v predicted: CD\ntrained: Blend : 2.86974090660627e-11 / diff: -20 %"
[1] "predicted: CD\ntrained: IBD v predicted: CD\ntrained: Composite : 2.28070548191642e-06 / diff: -12 %"
[1] "predicted: CD\ntrained: IBD v predicted: CD\ntrained: SMTPred : 0.625282568301732 / diff: -3 %"
[1] "predicted: CD\ntrained: Blend v predicted: CD\ntrained: Composite : 2.53410473770138e-12 / diff: 8 %"
[1] "predicted: CD\ntrained: Blend v predicted: CD\ntrained: SMTPred : 0.000382641877873765 / diff: 17 %"
[1] "predicted: CD\ntrained: Composite v predicted: CD\ntrained: SMTPred : 0.0584622025758863 / diff: 9 %"

[1] "predicted: CD\ntrained: Composite: 0.027652717"
[1] "predicted: CD\ntrained: CD: 0.0238835675"    XX
[1] "predicted: CD\ntrained: IBD: 0.0246327425"   XX
[1] "predicted: CD\ntrained: Blend: 0.030069379" XX
[1] "predicted: CD\ntrained: Composite: 0.027652717"
[1] "predicted: CD\ntrained: SMTPred: 0.0252843225"

0.024

#UC
[1] "predicted: UC\ntrained: UC v predicted: UC\ntrained: IBD : 0.0383744711456854 / diff: -15 %"
[1] "predicted: UC\ntrained: UC v predicted: UC\ntrained: Blend : 3.72100838505253e-06 / diff: -30 %"
[1] "predicted: UC\ntrained: UC v predicted: UC\ntrained: Composite : 0.000199358926241561 / diff: -25 %"
[1] "predicted: UC\ntrained: UC v predicted: UC\ntrained: SMTPred : 1.68749271161051e-05 / diff: -14 %"
[1] "predicted: UC\ntrained: IBD v predicted: UC\ntrained: Blend : 1.53270813547926e-06 / diff: -15 %"
[1] "predicted: UC\ntrained: IBD v predicted: UC\ntrained: Composite : 5.57058677652205e-05 / diff: -9 %"
[1] "predicted: UC\ntrained: IBD v predicted: UC\ntrained: SMTPred : 0.864579204037575 / diff: 1 %"
[1] "predicted: UC\ntrained: Blend v predicted: UC\ntrained: Composite : 9.27700586698735e-07 / diff: 6 %"
[1] "predicted: UC\ntrained: Blend v predicted: UC\ntrained: SMTPred : 0.000597488548923104 / diff: 16 %"
[1] "predicted: UC\ntrained: Composite v predicted: UC\ntrained: SMTPred : 0.0270014149717807 / diff: 10 %"

[1] "predicted: UC\ntrained: Composite: 0.014997798"
[1] "predicted: UC\ntrained: UC: 0.01171920475"   XX
[1] "predicted: UC\ntrained: IBD: 0.0136385461"    XX
[1] "predicted: UC\ntrained: Blend: 0.0158598175"  XX
[1] "predicted: UC\ntrained: Composite: 0.014997798"
[1] "predicted: UC\ntrained: SMTPred: 0.013515063"


3.72100838505253e-06

# Improvement:
#-15 %"


#############################










