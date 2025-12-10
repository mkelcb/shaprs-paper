####################################################
#
# Follow-up shaPRS Cross-Ancestry analysis that included PRS-CSx
#
#########################
# screen -r -D  28652.pts-0.node-11-2-3 # IBD
# screen -r -D  57179.pts-1.node-11-1-1 # cross ancestry



#########################
prscssorig=$crossAncestryRaw$'prscss/'


prscsrefs=$crossAncestryRaw$'prscsrefs/'
prscsscript=$crossAncestryRaw$'prscsscript/'
crossAncestryScratch=$asthmaCrossAncestry$'scratch/'

prscsxdir=$prscsscript$'PRScsx/'
shaPRS="_shaPRS"
NCORES_LDPred2=10
mkdir -p $prscsrefs
mkdir -p $prscsscript
mkdir -p $prscssorig


# download ref panels for EUR and EAS (need to add 'dl=1' to the file page to make wget work: https://superuser.com/questions/470664/how-to-download-dropbox-files-using-wget-command
cd $prscsrefs
wget https://www.dropbox.com/s/fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz?dl=1
wget https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=1
wget https://www.dropbox.com/s/oyn5trwtuei27qj/snpinfo_mult_ukbb_hm3?dl=1

# decompress them
mv snpinfo_mult_ukbb_hm3?dl=1 snpinfo_mult_ukbb_hm3
mv ldblk_ukbb_eur.tar.gz?dl=1 ldblk_ukbb_eur.tar.gz
mv ldblk_ukbb_eas.tar.gz?dl=1 ldblk_ukbb_eas.tar.gz

tar -xvf ldblk_ukbb_eas.tar.gz
tar -xvf ldblk_ukbb_eur.tar.gz

# get the prsCSx repository:
cd $prscsscript

git clone https://github.com/getian107/PRScsx.git



# run test to confirm that installation works
python $prscsxdir/PRScsx.py --help
python $prscsxdir/PRScsx.py --ref_dir=$prscsrefs --bim_prefix=$prscsxdir/test_data/test --sst_file=$prscsxdir/test_data/EUR_sumstats.txt,$prscsxdir/test_data/EAS_sumstats.txt --n_gwas=200000,100000 --pop=EUR,EAS --chrom=22 --phi=1e-2 --out_dir=$crossAncestryScratch --out_name=test


cd $prscssorig
git clone https://github.com/getian107/PRScs.git

python $prscssorig/PRScs/PRScs.py --help
################

# 0) Create Valid/Test divide in the UKBB
#$crossAncestrySumStats$'all_hm3.fam' 
#$crossAncestrySumStats$'all_height_hm3.fam' # these have the same number of indis in the same order

numItems=$(wc -l < "$crossAncestrySumStats"all_hm3.fam)
mySeed=42
arguments='/nfs/users/nfs_m/mk23/scripts/randomNums.R '$numItems$' '$mySeed$foldFam$' '$crossAncestrySumStats$'randomNums'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# create a a TEST set indis, the first 50% of them
awk -v numItems="$numItems" '{if (FNR < numItems/2) print $0}' $crossAncestrySumStats$'randomNums' > $crossAncestrySumStats$'randomNums_TEST'
head $crossAncestrySumStats$'randomNums_TEST'
wc -l $crossAncestrySumStats$'randomNums_TEST'


$crossAncestrySumStats$'randomNums'

awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $2 } }
' $crossAncestrySumStats$'randomNums_TEST' $crossAncestrySumStats$'all_hm3.fam' > $crossAncestrySumStats$'all_TEST' 
head $crossAncestrySumStats$'all_TEST' 
wc -l $crossAncestrySumStats$'all_TEST' # 125,630 , 125K test indis still


##############
# I) Convert sumstats into PRSCS format + build the initial EUR/EAS PRS


PRS_array=( 'EUR_JAP_asthma' 'EUR_JAP_height' 'EUR_JAP_BRCA' 'EUR_JAP_CAD' 'EUR_JAP_T2D' )
EUR_array=( 'EUR_asthma_hm3' 'EUR_height_hm3' 'EUR_BRCA_hm3' 'EUR_CAD_hm3' 'EUR_T2D_hm3' )
JPT_array=( 'JP_asthma_hm3' 'JP_height_hm3' 'JP_BRCA_hm3' 'JP_CAD_hm3' 'JP_T2D_hm3' )
EUR_CUSTOMLDREFARRAY=( $asthmaEURRef $heightEURRef $brcaEURRef $cadEURRef $t2dEURRef )
JAP_CUSTOMLDREFARRAY=( $asthmaEASRef $heightEASRef $brcaEASRef $cadEASRef $t2dEASRef )



PRS_binary=( '1' '0' '1' '1' '1' )
arraylength=${#PRS_array[@]}

j=1
i=12

for (( j=1; j<${arraylength}+1; j++ )); do
pheno=${PRS_array[$j-1]} 
binPheno=${PRS_binary[$j-1]} 
EURPheno=${EUR_array[$j-1]} 
JPTPheno=${JPT_array[$j-1]} 


# 1) subset full panel to the LDpred2 POST-QC set that I used for shaprs
# subset the SNPs to be the same panel for all 3 sumstats 
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums) {print $0 } }
' $crossAncestryResults$pheno$'_sumstats_meta_QC_keep' $eursumstatsLoc$EURPheno  > $eursumstatsLoc$EURPheno$'_QC'

awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $crossAncestryResults$pheno$'_sumstats_meta_QC_keep' $japsumstatsLoc$JPTPheno  > $japsumstatsLoc$JPTPheno$'_QC'



# 2) Conver them to PRS-CS format
convertToPRSCS $eursumstatsLoc$EURPheno$'_QC' $binPheno
convertToPRSCS $japsumstatsLoc$JPTPheno$'_QC' $binPheno

# 3) run PRS-CSx
Neur=$(tail -n 1 $eursumstatsLoc$EURPheno | awk '{ print $10}')
Neas=$(tail -n 1 $japsumstatsLoc$JPTPheno | awk '{ print $10}')

validBim=$eursumstatsLoc$EURPheno$'_validBim'
# create a 'fake' .bim file for PRSx from my sumstats file that has all the info
awk '{if (FNR !=1) {print $1"\t"$3"\t0\t"$2"\t"$4"\t"$5} }' $eursumstatsLoc$EURPheno$'_QC' > $validBim$'.bim'
mkdir -p $crossAncestryResults$'logs/'
for ((i=1; i<=22; i++)); do # $numChroms


# 1) PRS-CSx
# performPRSCSx $eursumstatsLoc$EURPheno$'_QC_PRSCS' $japsumstatsLoc$JPTPheno$'_QC_PRSCS' $Neur $Neas $eursumstatsLoc$EURPheno$'_validBim' $pheno $i 7

# check if output doesn't already exist
outfile=$crossAncestryResults$pheno$'_EUR_pst_eff_a1_b0.5_phiauto_chr'$i$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
echo "SUBMITTING: "$pheno$' '$i
bsub -q normal -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J ${pheno}_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $crossAncestryResults$'logs/'$pheno$'_'$i$'.out' -e $crossAncestryResults$'logs/'$pheno$'_'$i'.err'  "performPRSCSx $eursumstatsLoc$EURPheno$'_QC_PRSCS' $japsumstatsLoc$JPTPheno$'_QC_PRSCS' $Neur $Neas $eursumstatsLoc$EURPheno$'_validBim' $pheno $i $NCORES_LDPred2"

fi


done # end of chroms

done



# WAIT HERE UNTIL ALL FINISHED ON CLUSTER
###############

##############
# II) Validate & Evaluate


PRS_array=( 'EUR_JAP_asthma' 'EUR_JAP_height' 'EUR_JAP_BRCA' 'EUR_JAP_CAD' 'EUR_JAP_T2D' )
EUR_array=( 'EUR_asthma_hm3' 'EUR_height_hm3' 'EUR_BRCA_hm3' 'EUR_CAD_hm3' 'EUR_T2D_hm3' )
JPT_array=( 'JP_asthma_hm3' 'JP_height_hm3' 'JP_BRCA_hm3' 'JP_CAD_hm3' 'JP_T2D_hm3' )


famFiles=( $crossAncestrySumStats$'all_hm3.fam' $crossAncestrySumStats$'all_height_hm3.fam' )
PRS_binary=( '1' '0' '1' '1' '1' )
arraylength=${#PRS_array[@]}

j=1
for (( j=1; j<${arraylength}+1; j++ )); do
pheno=${PRS_array[$j-1]} 
binPheno=${PRS_binary[$j-1]} 
EURPheno=${EUR_array[$j-1]} 
JPTPheno=${JPT_array[$j-1]} 

# 1) concat all chrom PRS into single files
cat $crossAncestryResults$pheno$'_EUR_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $crossAncestryResults$pheno$'_EUR_PRSCSx'
cat $crossAncestryResults$pheno$'_EAS_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $crossAncestryResults$pheno$'_EAS_PRSCSx'

# remove dupes:
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $crossAncestryResults$pheno$'_EUR_PRSCSx' > $crossAncestryResults$pheno$'_EUR_PRSCSx_no_dupes'
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $crossAncestryResults$pheno$'_EAS_PRSCSx' > $crossAncestryResults$pheno$'_EAS_PRSCSx_no_dupes'



#### Only for the first 2 traits that I have at Sanger: (the others are sent to Elena)

if [ "$j" -lt 3 ]; then
famFile=${famFiles[$j-1]} 
# 2) build PRS profile for both EUR/EAS on the full UKBB

# PRSCSx produces a PRS file of format:
#  1       2            3       4   5         6
# 22	rs9605903	17054720	C	T	4.954646e-04
# PLINK wants to know the variant ID, A1 and coef: these are 2 4 6

#EUR 
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$crossAncestrySumStats$'all_hm3 --fam '$famFile$' --score '$crossAncestryResults$pheno$'_EUR_PRSCSx_no_dupes 2 4 6 sum --out '$crossAncestryResults$pheno$'_EUR_PRSCSx.PRS'
$plink $arguments

#EAS 
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$crossAncestrySumStats$'all_hm3 --fam '$famFile$' --score '$crossAncestryResults$pheno$'_EAS_PRSCSx_no_dupes 2 4 6 sum --out '$crossAncestryResults$pheno$'_EAS_PRSCSx.PRS'
$plink $arguments



# 3) Perform Stage 2: find the best linear combination of the EUR/EAS + evaluate
arguments='/nfs/users/nfs_m/mk23/scripts/PRSLinearComb2.R '$crossAncestryResults$pheno$'_EUR_PRSCSx.PRS.profile '$crossAncestryResults$pheno$'_EAS_PRSCSx.PRS.profile '$famFile$' '$crossAncestrySumStats$'all_TEST '$binPheno$' '$crossAncestryResults$pheno$'_PRSCSx_result'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# PRS-CSx-weighted height:

# PRS-CSx-weighted asthma:
#  correlation_sq 0.0134 / AUC 0.603 
fi

done






###############################################################

# because PRS-CSx does its own QC routine on SNPs, we loose a few K more SNPs if we start from the LDpred2 panel
# to ensure we work off the exact same subset, we need to subset the LDpred2 SNPs with the PRS-CSx, and re-run shaPRS:




PRS_array=( 'EUR_JAP_asthma' 'EUR_JAP_height' 'EUR_JAP_BRCA' 'EUR_JAP_CAD' 'EUR_JAP_T2D' )
EUR_array=( 'EUR_asthma_hm3' 'EUR_height_hm3' 'EUR_BRCA_hm3' 'EUR_CAD_hm3' 'EUR_T2D_hm3' )
JPT_array=( 'JP_asthma_hm3' 'JP_height_hm3' 'JP_BRCA_hm3' 'JP_CAD_hm3' 'JP_T2D_hm3' )
EUR_CUSTOMLDREFARRAY=( $asthmaEURRef $heightEURRef $brcaEURRef $cadEURRef $t2dEURRef )
JAP_CUSTOMLDREFARRAY=( $asthmaEASRef $heightEASRef $brcaEASRef $cadEASRef $t2dEASRef )
PRS_JAP_array=( 'JAP_EUR_asthma' 'JAP_EUR_height' 'JAP_EUR_BRCA' 'JAP_EUR_CAD' 'JAP_EUR_T2D' )
 

famFiles=( $crossAncestrySumStats$'all_hm3.fam' $crossAncestrySumStats$'all_height_hm3.fam' )
PRS_binary=( '1' '0' '1' '1' '1' )
arraylength=${#PRS_array[@]}



j=4
i=21
for (( j=1; j<${arraylength}+1; j++ )); do
pheno=${PRS_array[$j-1]} 
binPheno=${PRS_binary[$j-1]} 
EURPheno=${EUR_array[$j-1]} 
JPTPheno=${JPT_array[$j-1]} 
EUR_PRSCSRef=${EUR_CUSTOMLDREFARRAY[$j-1]} 
JAP_PRSCSRef=${JAP_CUSTOMLDREFARRAY[$j-1]} 

pheno_JAP=${PRS_JAP_array[$j-1]} 

Neur=$(tail -n 1 $eursumstatsLoc$EURPheno | awk '{ print $10}')
Neas=$(tail -n 1 $japsumstatsLoc$JPTPheno | awk '{ print $10}')
validBim=$eursumstatsLoc$EURPheno$'_validBim'



# the final PRS for PRS-CSx is a composite of 2 PRS, which had a different set of SNPs
# we are comparing ourselves most directly to the EUR PRS-CSx, so use that as a subset
# wc -l $crossAncestryResults$pheno$'_EAS_PRSCSx_no_dupes'

awk 'FNR == NR { sums[ $2 ] = $2; next; } FNR <= NR { if( FNR == 1 || $3 in sums) {print $0 } }
' $crossAncestryResults$pheno$'_EUR_PRSCSx_no_dupes' $eursumstatsLoc$EURPheno$'_QC' > $eursumstatsLoc$EURPheno$'_QC_CSxsubset'

awk 'FNR == NR { sums[ $2 ] = $2; next; } FNR <= NR { if( FNR == 1 || $3 in sums) {print $0 } }
' $crossAncestryResults$pheno$'_EUR_PRSCSx_no_dupes' $crossAncestryResults$pheno'_sumstats_meta_QC' > $crossAncestryResults$pheno'_sumstats_meta_QC_CSxsubset'


# Do QC for the JAP ShaPRS sumstats too to bring it in line with the rest
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $crossAncestryResults$pheno$'_sumstats_meta_QC_keep' $crossAncestryResults$pheno_JAP$'_sumstats_meta'  > $crossAncestryResults$pheno_JAP$'_sumstats_meta_QC'

awk 'FNR == NR { sums[ $2 ] = $2; next; } FNR <= NR { if( FNR == 1 || $3 in sums) {print $0 } }
' $crossAncestryResults$pheno$'_EUR_PRSCSx_no_dupes' $crossAncestryResults$pheno_JAP'_sumstats_meta_QC' > $crossAncestryResults$pheno_JAP'_sumstats_meta_QC_CSxsubset'

# 2) Convert them to PRS-CS format
convertToPRSCS $eursumstatsLoc$EURPheno$'_QC_CSxsubset' $binPheno
convertToPRSCS $crossAncestryResults$pheno'_sumstats_meta_QC_CSxsubset' $binPheno
convertToPRSCS $crossAncestryResults$pheno_JAP'_sumstats_meta_QC_CSxsubset' $binPheno





for ((i=1; i<=22; i++)); do # $numChroms

# EUR+LDpred2
outfile=$crossAncestryResults$'EUR_EUR_LD_'$pheno$'_CSxsubset_'$i
#mv $outfile $outfile$'_empirical' # backup originals
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
#  EUR  (EUR LD panel)
echo "LDpred SUBMITTING: "$pheno$'_EUR crhom'$i
avgNumIndis=$(tail -n 1 $eursumstatsLoc$EURPheno | awk '{ print $10}')
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_chrom.R /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/ '$eursumstatsLoc$EURPheno$'_QC_CSxsubset '$i$' '$avgNumIndis$' '$crossAncestryResults$' EUR_EUR_LD_'$pheno$'_CSxsubset_'$i$' '$crossAncestrySumStats$'all_hm3 '$shaPRSscriptLoc$' '$binPheno
#bsub -q normal -m "modern_hardware" -G team152 -n1 -J EUR_EUR_LD_${pheno}_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $crossAncestryResults$'EUR_EUR_LD_'$pheno$'_'$i$'.out' -e $crossAncestryResults$'EUR_EUR_LD_'$pheno$'_'$i$'.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"


fi

# ShaPRS+ LDpred2
outfile=$crossAncestryResults$'shaPRS_LD_'$pheno$'_CSxsubset_'$i
#mv $outfile $outfile$'_empirical' # backup originals
#rm -rf $outfile
if [ -s "$outfile" ] ; then
b=1
#echo "ALREADY EXISTS: "$outfile
else
# EUR-shaPRS shaPRS LD Panel
echo "ShaPRS SUBMITTING: "$pheno$'_EUR crhom'$i
avgNumIndis=$(tail -n 1 "$crossAncestryResults$pheno"_sumstats_meta | awk '{ print $10}')
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_chrom.R '$crossAncestryResults$pheno$'/ '$crossAncestryResults$pheno'_sumstats_meta_QC_CSxsubset '$i$' '$avgNumIndis$' '$crossAncestryResults$' shaPRS_LD_'$pheno$'_CSxsubset_'$i$' '$crossAncestrySumStats$'all_hm3 '$shaPRSscriptLoc$' '$binPheno
#bsub -q normal -m "modern_hardware" -G team152 -n1 -J shaPRS_LD_${pheno}_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $crossAncestryResults$'shaPRS_LD_'$pheno$'_'$i$'.out' -e $crossAncestryResults$'shaPRS_LD_'$pheno$'_'$i$'.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"
fi



# EUR PRS-CS
# check if output doesn't already exist
outfile=$crossAncestryResults$pheno$'_EUR_EUR_pst_eff_a1_b0.5_phiauto_chr'$i$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
b=1
#echo "PRS-CS SUBMITTING: "$pheno$'_EUR '$i
#bsub -q normal -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J ${pheno}_EUR_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $crossAncestryResults$'logs/'$pheno$'_EUR_'$i$'.out' -e $crossAncestryResults$'logs/'$pheno$'_EUR_'$i'.err'  "performPRSCS $eursumstatsLoc$EURPheno$'_QC_CSxsubset_PRSCS' $Neur $eursumstatsLoc$EURPheno$'_validBim' $pheno$'_EUR' $i $NCORES_LDPred2"

fi




# ShaPRS-PRS-CS (EUR)
# check if output doesn't already exist
outfile=$crossAncestryResults$pheno$shaPRS$'_EUR/_pst_eff_a1_b0.5_phiauto_chr'$i$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
b=1
echo "PRS-CS SUBMITTING: EUR "$pheno$shaPRS$' '$i
# custom LD panel
#bsub -q normal -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J ${pheno}${shaPRS}_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $crossAncestryResults$'logs/'$pheno$shaPRS$'_'$i$'.out' -e $crossAncestryResults$'logs/'$pheno$shaPRS$'_'$i'.err'  "performPRSCS_custom $crossAncestryResults$pheno'_sumstats_meta_QC_CSxsubset_PRSCS' $Neur $eursumstatsLoc$EURPheno$'_validBim' $pheno$shaPRS$'_EUR' $i $NCORES_LDPred2 $prscsrefs'ldblk_ukbb_'$EUR_PRSCSRef"

# original LD panel
#bsub -q normal -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J ${pheno}${shaPRS}_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $crossAncestryResults$'logs/'$pheno$shaPRS$'_'$i$'.out' -e $crossAncestryResults$'logs/'$pheno$shaPRS$'_'$i'.err'  "performPRSCS_custom $crossAncestryResults$pheno'_sumstats_meta_QC_CSxsubset_PRSCS' $Neur $eursumstatsLoc$EURPheno$'_validBim' $pheno$shaPRS$'_EUR' $i $NCORES_LDPred2 $prscsrefs'ldblk_ukbb_eur'"

fi


# ShaPRS-PRS-CS (EAS)
# check if output doesn't already exist
outfile=$crossAncestryResults$pheno$shaPRS$'_EAS/_pst_eff_a1_b0.5_phiauto_chr'$i$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
b=1
echo "PRS-CS SUBMITTING: EAS "$pheno$shaPRS$' '$i
# custom LD panel
#bsub -q normal -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J ${pheno}${shaPRS}_EAS_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $crossAncestryResults$'logs/'$pheno$shaPRS$'_EAS_'$i$'.out' -e $crossAncestryResults$'logs/'$pheno$shaPRS$'_EAS_'$i'.err'  "performPRSCS_custom $crossAncestryResults$pheno_JAP'_sumstats_meta_QC_CSxsubset_PRSCS' $Neas $eursumstatsLoc$EURPheno$'_validBim' $pheno$shaPRS$'_EAS' $i $NCORES_LDPred2 $prscsrefs'ldblk_ukbb_'$JAP_PRSCSRef"

# original LD panel
#bsub -q normal -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J ${pheno}${shaPRS}_EAS_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $crossAncestryResults$'logs/'$pheno$shaPRS$'_EAS_'$i$'.out' -e $crossAncestryResults$'logs/'$pheno$shaPRS$'_EAS_'$i'.err'  "performPRSCS_custom $crossAncestryResults$pheno_JAP'_sumstats_meta_QC_CSxsubset_PRSCS' $Neas $eursumstatsLoc$EURPheno$'_validBim' $pheno$shaPRS$'_EAS' $i $NCORES_LDPred2 $prscsrefs'ldblk_ukbb_eas'"


fi


# 4) PRS-CSx of ShaPRS EUR+EAS
# check if output doesn't already exist
outfile=$crossAncestryResults$pheno$'_PRSCSx_EUR_pst_eff_a1_b0.5_phiauto_chr'$i$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
echo "SUBMITTING: ShaPRS PRS-CSx "$pheno$' '$i
#bsub -q normal -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J ${pheno}_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $crossAncestryResults$'logs/'$pheno$'_'$i$'.out' -e $crossAncestryResults$'logs/'$pheno$'_'$i'.err'  "performPRSCSx $crossAncestryResults$pheno'_sumstats_meta_QC_CSxsubset_PRSCS' $crossAncestryResults$pheno_JAP'_sumstats_meta_QC_CSxsubset_PRSCS' $Neur $Neas $eursumstatsLoc$EURPheno$'_validBim' $pheno$'_PRSCSx' $i $NCORES_LDPred2"

fi

done # end of chroms

done # end of traits



# !!! HERE, need to write PRS-CSx-ShaPRS evaluation script (and also to test CAD)


echo $crossAncestryResults$pheno$'_sumstats_meta'
echo $crossAncestryResults$pheno_JAP$'_sumstats_meta'





# EVALUATE SHAPRS
PRS_array=( 'EUR_JAP_asthma' 'EUR_JAP_height' 'EUR_JAP_BRCA' 'EUR_JAP_CAD' 'EUR_JAP_T2D' )
EUR_array=( 'EUR_asthma_hm3' 'EUR_height_hm3' 'EUR_BRCA_hm3' 'EUR_CAD_hm3' 'EUR_T2D_hm3' )
JPT_array=( 'JP_asthma_hm3' 'JP_height_hm3' 'JP_BRCA_hm3' 'JP_CAD_hm3' 'JP_T2D_hm3' )


famFiles=( $crossAncestrySumStats$'all_hm3.fam' $crossAncestrySumStats$'all_height_hm3.fam' )
PRS_binary=( '1' '0' '1' '1' '1' )
arraylength=${#PRS_array[@]}



for (( j=1; j<${arraylength}+1; j++ )); do
pheno=${PRS_array[$j-1]} 
binPheno=${PRS_binary[$j-1]} 
EURPheno=${EUR_array[$j-1]} 
JPTPheno=${JPT_array[$j-1]} 
JAP_PRSCSRef=${JAP_CUSTOMLDREFARRAY[$j-1]} 

# 1) concat all chrom PRS into single files

# a) PRS-CS
cat $crossAncestryResults$pheno$'_EUR_EUR_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $crossAncestryResults$pheno$'_EUR_EUR_PRSCS'
cat $crossAncestryResults$pheno$shaPRS$'_EUR/_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS'
cat $crossAncestryResults$pheno$shaPRS$'_EAS/_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS'

# remove dupes:
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $crossAncestryResults$pheno$'_EUR_EUR_PRSCS' > $crossAncestryResults$pheno$'_EUR_EUR_PRSCS_no_dupes'
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS' > $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS_no_dupes'
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS' > $crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS_no_dupes'

# b) LDpred2
cat $crossAncestryResults$'EUR_EUR_LD_'$pheno$'_CSxsubset_'[!a-z] $crossAncestryResults$'EUR_EUR_LD_'$pheno$'_CSxsubset_'[!a-z][!a-z] > $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset'
cat $crossAncestryResults$'shaPRS_LD_'$pheno$'_CSxsubset_'[!a-z] $crossAncestryResults$'shaPRS_LD_'$pheno$'_CSxsubset_'[!a-z][!a-z] > $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset'

awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset' > $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset_no_dupes'
awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset' > $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset_no_dupes'


# c) PRS-CSx ShaPRS:
cat $crossAncestryResults$pheno$'_PRSCSx_EUR_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS'
cat $crossAncestryResults$pheno$'_PRSCSx_EAS_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS'

# remove dupes:
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS' > $crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS_no_dupes'
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS' > $crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS_no_dupes'






#### Only for the first 2 traits that I have at Sanger: (the others are sent to Elena)

if [ "$j" -lt 3 ]; then
famFile=${famFiles[$j-1]} 
# 2) build PRS profile for both EUR/EAS on the full UKBB

# a) PRS-CS
# PRSCSx produces a PRS file of format:
#  1       2            3       4   5         6
# 22	rs9605903	17054720	C	T	4.954646e-04
# PLINK wants to know the variant ID, A1 and coef: these are 2 4 6

#EUR EUR
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$crossAncestrySumStats$'all_hm3 --fam '$famFile$' --score '$crossAncestryResults$pheno$'_EUR_EUR_PRSCS_no_dupes 2 4 6 sum --out '$crossAncestryResults$pheno$'_EUR_EUR_PRSCS.PRS'
$plink $arguments

#SHAPRS EUR
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$crossAncestrySumStats$'all_hm3 --fam '$famFile$' --score '$crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS_no_dupes 2 4 6 sum --out '$crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS.PRS'
$plink $arguments

#SHAPRS EAS
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$crossAncestrySumStats$'all_hm3 --fam '$famFile$' --score '$crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS_no_dupes 2 4 6 sum --out '$crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS.PRS'
$plink $arguments


# evaluate the unweighted PRS-CS ones   
evaluatePRS_NOVALID $famFile $pheno$'_EUR_EUR_PRSCS.PRS' $binPheno

# asthma: correlation_sq 0.0099 / AUC 0.5884
# height: correlation_sq 0.1161

evaluatePRS_NOVALID $famFile $pheno$shaPRS$'_EUR_PRSCS.PRS' $binPheno
# asthma: 
# UNI LD  : correlation_sq 0.0123 / AUC 0.5986     / with EAS panel: correlation_sq 0.0125 /  AUC 0.5993
# EUR LD  : correlation_sq 0.0121 / AUC 0.5982


# height: 
# UNI LD  : correlation_sq 0.1216
# EUR LD  : correlation_sq 0.1226  


#SHAPRS PRS BLEND
# 3) Perform Stage 2: find the best linear combination of the EUR/EAS + evaluate
arguments='/nfs/users/nfs_m/mk23/scripts/PRSLinearComb2.R '$crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS.PRS.profile '$crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS.PRS.profile '$famFile$' '$crossAncestrySumStats$'all_TEST '$binPheno$' '$crossAncestryResults$pheno$shaPRS$'_EUR_EAS_PRSCS_result'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# asthma:
# UNI LD  : correlation_sq 0.0126 / AUC 0.6001
# EUR LD  : correlation_sq 0.0123 / AUC 0.5992


# height: 
# UNI LD  : correlation_sq 0.122
# EUR LD  : correlation_sq 0.1232


# b) LDpred2
# EUR EUR baseline
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$crossAncestrySumStats$'all_hm3 --fam '$famFile$' --score '$crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset_no_dupes sum --out '$crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset.PRS'
$plink $arguments

# shaPRS
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$crossAncestrySumStats$'all_hm3 --fam '$famFile$' --score '$crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset_no_dupes sum --out '$crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS'
$plink $arguments

# evaluate the LDpred2 PRS too
evaluatePRS_NOVALID $famFile $pheno$'_EUR_EUR_LDpred2_CSxsubset.PRS' $binPheno
# asthma: correlation_sq 0.0115 /  AUC 0.5954
# Height: correlation_sq 0.0976 

evaluatePRS_NOVALID $famFile $pheno$'_shaPRS_LDpred2_CSxsubset.PRS' $binPheno
# asthma: correlation_sq 0.0136 / AUC 0.6038
# Height: correlation_sq 0.1218


#####################
# PRS-CSx ShaPRS
#EUR 
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$crossAncestrySumStats$'all_hm3 --fam '$famFile$' --score '$crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS_no_dupes 2 4 6 sum --out '$crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS.PRS'
$plink $arguments

#EAS 
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$crossAncestrySumStats$'all_hm3 --fam '$famFile$' --score '$crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS_no_dupes 2 4 6 sum --out '$crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS.PRS'
$plink $arguments

# 3) Perform Stage 2: find the best linear combination of the EUR/EAS + evaluate
arguments='/nfs/users/nfs_m/mk23/scripts/PRSLinearComb2.R '$crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS.PRS.profile '$crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS.PRS.profile '$famFile$' '$crossAncestrySumStats$'all_TEST '$binPheno$' '$crossAncestryResults$pheno$'_PRSCSx_result_ShaPRS'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# asthma: 
# correlation_sq 0.012 /  AUC 0.5975

# height:
#  correlation_sq 0.1245 

#####################

fi

done



###############################################################a
$crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset_no_dupes'
# convert PRS for BRCA, CAD and T2D into Elena's format 

# LDpred2
pheno="EUR_JAP_CAD"
convertPRS $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset'
convertPRS $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset'
# compress them 
zip $crossAncestryResults$pheno'_LDpred2.zip' $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset_E' $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset_E'


pheno="EUR_JAP_T2D"
convertPRS $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset'
convertPRS $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset'
# compress them 
zip $crossAncestryResults$pheno'_LDpred2.zip' $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset_E' $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset_E'


pheno="EUR_JAP_BRCA"
convertPRS $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset'
convertPRS $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset'
# compress them 
zip $crossAncestryResults$pheno'_LDpred2.zip' $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset_E' $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset_E'


# PRS-CSx:
pheno="EUR_JAP_CAD"
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EUR_PRSCSx'
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EAS_PRSCSx'
# compress them 
zip $crossAncestryResults$pheno'_PRSCSx.zip' $crossAncestryResults$pheno$'_EUR_PRSCSx_E' $crossAncestryResults$pheno$'_EAS_PRSCSx_E'

pheno="EUR_JAP_T2D"
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EUR_PRSCSx'
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EAS_PRSCSx'
# compress them 
zip $crossAncestryResults$pheno'_PRSCSx.zip' $crossAncestryResults$pheno$'_EUR_PRSCSx_E' $crossAncestryResults$pheno$'_EAS_PRSCSx_E'

pheno="EUR_JAP_BRCA"
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EUR_PRSCSx'
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EAS_PRSCSx'
# compress them 
zip $crossAncestryResults$pheno'_PRSCSx.zip' $crossAncestryResults$pheno$'_EUR_PRSCSx_E' $crossAncestryResults$pheno$'_EAS_PRSCSx_E'


# PRS-CS:
pheno="EUR_JAP_CAD"
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EUR_EUR_PRSCS'
convertPRS_PRSCSx $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS'
convertPRS_PRSCSx $crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS'
# compress them 
zip $crossAncestryResults$pheno'_PRSCS.zip' $crossAncestryResults$pheno$'_EUR_EUR_PRSCS_E' $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS_E' $crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS_E'



pheno="EUR_JAP_T2D"
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EUR_EUR_PRSCS'
convertPRS_PRSCSx $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS'
convertPRS_PRSCSx $crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS'
# compress them 
zip $crossAncestryResults$pheno'_PRSCS.zip' $crossAncestryResults$pheno$'_EUR_EUR_PRSCS_E' $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS_E' $crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS_E'

pheno="EUR_JAP_BRCA"
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EUR_EUR_PRSCS'
convertPRS_PRSCSx $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS'
convertPRS_PRSCSx $crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS'
# compress them 
zip $crossAncestryResults$pheno'_PRSCS.zip' $crossAncestryResults$pheno$'_EUR_EUR_PRSCS_E' $crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS_E' $crossAncestryResults$pheno$shaPRS$'_EAS_PRSCS_E'




# PRS-CSx ShaPRS:
pheno="EUR_JAP_CAD"
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS'
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS'
# compress them 
zip $crossAncestryResults$pheno'_PRSCSx_ShaPRS.zip' $crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS_E' $crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS_E'

pheno="EUR_JAP_T2D"
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS'
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS'
# compress them 
zip $crossAncestryResults$pheno'_PRSCSx_ShaPRS.zip' $crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS_E' $crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS_E'

pheno="EUR_JAP_BRCA"
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS'
convertPRS_PRSCSx $crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS'
# compress them 
zip $crossAncestryResults$pheno'_PRSCSx_ShaPRS.zip' $crossAncestryResults$pheno$'_EUR_PRSCSx_ShaPRS_E' $crossAncestryResults$pheno$'_EAS_PRSCSx_ShaPRS_E'



# Zip up all final PRS files for all traits

# PRS-CS
$crossAncestryResults$pheno$'_EUR_EUR_PRSCS' 
$crossAncestryResults$pheno$shaPRS$'_EUR_PRSCS' 

# PRS-CSx
$crossAncestryResults$pheno$'_EUR_PRSCSx'
$crossAncestryResults$pheno$'_EAS_PRSCSx'

#LDpred
$crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset'
$crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset'


$crossAncestryResults$pheno$'_EUR_EUR_PRSCS' $crossAncestryResults$pheno$'_shaPRS_EUR_PRSCS' $crossAncestryResults$pheno$'_EUR_PRSCSx' $crossAncestryResults$pheno$'_EAS_PRSCSx' $crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset' $crossAncestryResults$pheno$'_shaPRS_LDpred2_CSxsubset'

zip -j $crossAncestryResults$'ALL_PRS.zip' 

'EUR_JAP_CAD'
'EUR_JAP_T2D'
'EUR_JAP_BRCA'
'EUR_JAP_asthma'
'EUR_JAP_height'

$crossAncestryResults$'EUR_JAP_CAD'$'_EUR_EUR_PRSCS' $crossAncestryResults$'EUR_JAP_CAD'$shaPRS$'_EUR_PRSCS' $crossAncestryResults$'EUR_JAP_CAD'$'_EUR_PRSCSx' $crossAncestryResults$'EUR_JAP_CAD'$'_EAS_PRSCSx' $crossAncestryResults$'EUR_JAP_CAD'$'_EUR_EUR_LDpred2_CSxsubset' $crossAncestryResults$'EUR_JAP_CAD'$'_shaPRS_LDpred2_CSxsubset'



#######################
# MODEL EVALUATION

# PRS-CSx-unweighted:

#height
pheno='EUR_JAP_height'
famFile=$crossAncestrySumStats$'all_height_hm3.fam'

# LDpred2 EUR baseline
evaluatePRS_NOVALID $famFile $pheno$'_EUR_EUR_LDpred2_CSxsubset.PRS' 0
#  correlation_sq 0.0976 (sd: 6.31e-06 )


# PRS-CS baseline
evaluatePRS_NOVALID $famFile $pheno$'_EUR_EUR_PRSCS.PRS' 0
# correlation_sq 0.116 (sd: 5.96e-06 )

# PRS-CSx
arguments='/nfs/users/nfs_m/mk23/scripts/PRSLinearComb2.R '$crossAncestryResults$pheno$'_EUR_PRSCSx.PRS.profile '$crossAncestryResults$pheno$'_EAS_PRSCSx.PRS.profile '$famFile$' '$crossAncestrySumStats$'all_TEST 0 '$crossAncestryResults$pheno$'_PRSCSx_result'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
# correlation_sq 0.123 (sd: 1.1e-05 )

# PRS-CSx-unweighted
evaluatePRS_NOVALID $famFile $pheno$'_EUR_PRSCSx.PRS' '0'  
#correlation_sq 0.121 (sd: 5.48e-06 )


# LDpred2 shaPRS unweighted
evaluatePRS_NOVALID $famFile $pheno$'_shaPRS_LDpred2_CSxsubset.PRS' 0
#  correlation_sq 0.122 (sd: 5.37e-06 )


# PRS-CS- shaPRS unweighted
evaluatePRS_NOVALID $famFile $pheno$shaPRS$'_EUR_PRSCS.PRS' 0
#  correlation_sq 0.122 (sd: 5.71e-06 )





###############################
# asthma:
pheno='EUR_JAP_asthma'
famFile=$crossAncestrySumStats$'all_hm3.fam'

# LDpred2 EUR baseline
evaluatePRS_NOVALID $famFile $pheno$'_EUR_EUR_LDpred2_CSxsubset.PRS' 1
# correlation_sq 0.0115 (sd: 5.52e-06 )   AUC 0.595 CI low: 0.5919  / high:  0.5989

# PRS-CS baseline
evaluatePRS_NOVALID $famFile $pheno$'_EUR_EUR_PRSCS.PRS' 1
#  correlation_sq 0.00986 (sd: 5.65e-06 ) AUC 0.588 CI low: 0.5849  / high:  0.5918


# PRS-CSx
arguments='/nfs/users/nfs_m/mk23/scripts/PRSLinearComb2.R '$crossAncestryResults$pheno$'_EUR_PRSCSx.PRS.profile '$crossAncestryResults$pheno$'_EAS_PRSCSx.PRS.profile '$famFile$' '$crossAncestrySumStats$'all_TEST 1 '$crossAncestryResults$pheno$'_PRSCSx_result'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments
# correlation_sq 0.0134 (sd: 1.11e-05 )  AUC 0.603 CI low: 0.5981  / high:  0.6079

# PRS-CSx-unweighted
evaluatePRS_NOVALID $famFile $pheno$'_EUR_PRSCSx.PRS' '1'  
# correlation_sq 0.0113 (sd: 6.35e-06 )  AUC 0.595 CI low: 0.5916  / high:  0.5985



# LDpred2 shaPRS unweighted
evaluatePRS_NOVALID $famFile $pheno$'_shaPRS_LDpred2_CSxsubset.PRS' 1
#  correlation_sq 0.0136 (sd: 5.53e-06 )  AUC 0.604 CI low: 0.6004  / high:  0.6073


# PRS-CS- shaPRS unweighted
evaluatePRS_NOVALID $famFile $pheno$shaPRS$'_EUR_PRSCS.PRS' 1
#  correlation_sq 0.0123 (sd: 5.55e-06 )   AUC 0.599 CI low: 0.5951  / high:  0.602







             PRS-CSx      PRS-CSx(no VALID)             shaPRS
 asthma:      
       r^2   0.0134        0.0113                       0.0138
       AUC   0.603         0.595                        0.6044

             PRS-CSx      PRS-CSx(no VALID)             shaPRS
 height:   
       r^2   0.1234        0.1207                       0.1146

#Asthma PRS-CSx:
#JP  --score: 651179 valid predictors loaded
#EUR --score: 652984 valid predictors loaded.



#height PRS-CSx:
#JP: --score: 597391 valid predictors loaded.
#EUR --score: 597463  valid predictors loaded.

# RESULTS:



# asthma EUR: --score: 652984 valid predictors loaded.
# correlation_sq 0.0099 
#  AUC 0.5884

# asthma shaPRS: --score: 652984 valid predictors loaded.
#  correlation_sq 0.0121 
# : AUC 0.5982

# height EUR: --score: 597463 valid predictors loade
# correlation_sq 0.1161

# height ShaPRS: --score: 597463 valid predictors 
# correlation_sq 0.1226



             PRS-CSx      PRS-CSx(no VALID)             shaPRS+PRS-CS			shaPRS+LDPred2
 asthma:      
       r^2   0.0134        0.0113                       0.0121 					0.0138
       AUC   0.603         0.595                        0.5982					0.6044

             PRS-CSx      PRS-CSx(no VALID)             shaPRS+PRS-CS			shaPRS+LDPred2
 height:   
       r^2   0.1234        0.1207                       0.1226					0.1146
	   
	   
	   
.............PRS-CSx      PRS-CSx(no VALID)   shaPRS+PRS-CS			shaPRS+LDPred2
 asthma:      
       r^2   0.0134        0.0113             0.0121 					0.0138
       AUC   0.603         0.595              0.5982					0.6044
 height:   
       r^2   0.1234        0.1207             0.1226				    0.1146
	   
 BRCA:
       r^2   0.0084		   0.0068			  0.0073					0.0095
       AUC   0.60	       0.59				  0.60						0.61
	   
 T2D:
       r^2   0.0125		   0.0114			  0.0125					0.0102
       AUC   0.66	  	   0.65 			  0.66						0.64
	   
 CAD:
       r^2   0.0179		   0.0179			  0.0166					0.0171
       AUC   0.67		   0.67				  0.66						0.67




#################################################
# Functions:

# generates an output in Elena's format, with a postfix of '_E' added to the input
function convertPRS_PRSCSx {
PRSFILE=$1

# remove dupes
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $PRSFILE > $PRSFILE$'_no_dupes'

# PRSCSx signature
#  1          2              3    4       5            6
# 10      rs10904494      113934  C       A       1.963357e-06



# Elena's signature: CHR POS SNP A1 A2 WEIGHT
awk '
FNR <= NR { 
{ 
if (FNR == 1) {print "CHR\tPOS\tSNP\tA1\tA2\tWEIGHT"}

print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6


 } }' $PRSFILE$'_no_dupes'  > $PRSFILE$'_E'

#head $PRSFILE$'_E'

}










# evaluates ready made profile scores
function evaluatePRS_NOVALID {
famFile=$1
PRSName=$2
calcAUC3=$3
# call PLINK's score to build PRS for each



# find out accuracy 
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $crossAncestryResults$PRSName$'.profile' FS=" " $famFile  > $outste$crossAncestryResults$PRSName$'.csv'


arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$outste$crossAncestryResults$PRSName$'.csv '$crossAncestryResults$PRSName$'_resultNovalid '$calcAUC3
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments 
}




# need to create a 2 column keep file for plink
paste $crossAncestrySumStats$'all_TEST' $crossAncestrySumStats$'all_TEST'  > $crossAncestrySumStats$'all_TEST_PLINK' 


# Asthma
pheno='EUR_JAP_asthma'
evaluatePRS_TEST $crossAncestrySumStats$'all_hm3.fam' 'shaPRS_LD_'$pheno '1' $crossAncestrySumStats$'all_TEST_PLINK'
#  correlation_sq 0.0138 AUC 0.6044

# height
pheno='EUR_JAP_height'
evaluatePRS_TEST $crossAncestrySumStats$'all_height_hm3.fam' 'shaPRS_LD_'$pheno '0' $crossAncestrySumStats$'all_TEST_PLINK' 


famFile=$crossAncestrySumStats$'all_hm3.fam'
testSet=$crossAncestrySumStats$'all_TEST_PLINK' 
PRSName='shaPRS_LD_'$pheno
calcAUC3='1'


# same as evaluatePRS, but restricts the test set to the same 50% of indis that we used for PRSCSx
function evaluatePRS_TEST {
famFile=$1
PRSName=$2
calcAUC3=$3
testSet=$4
# call PLINK's score to build PRS for each

awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $crossAncestryResults$PRSName > $crossAncestryResults$PRSName$'_no_dupes'

arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$crossAncestrySumStats$'all_hm3 --fam '$famFile$' --score '$crossAncestryResults$PRSName$'_no_dupes sum --keep '$testSet$' --out '$crossAncestryResults$PRSName
$plink $arguments

# find out accuracy 
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $crossAncestryResults$PRSName$'.profile' FS=" " $famFile  > $outste$crossAncestryResults$PRSName$'.csv'
#awk '{count[$6]++} END {for (word in count) print word, count[word]}' $famFile
#awk '{count[$2]++} END {for (word in count) print word, count[word]}' FS="," $outste$crossAncestryResults$PRSName$'.csv'

arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$outste$crossAncestryResults$PRSName$'.csv '$crossAncestryResults$PRSName$'_result '$calcAUC3
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments 
}


# Performs PRS-CS for a custom population
function performPRSCS_custom { 
sumsFile=$1
NGWAS=$2
validBim=$3
mypheno=$4
chrom=$5
ncores=$6
refLoc=$7

export MKL_NUM_THREADS=$ncores
export NUMEXPR_NUM_THREADS=$ncores
export OMP_NUM_THREADS=$ncores



prscssorig="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscss/"
crossAncestryResults="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/"
mkdir -p $crossAncestryResults$mypheno
python $prscssorig/PRScs/PRScs.py --ref_dir=$refLoc --bim_prefix=$validBim --sst_file=$sumsFile --n_gwas=$NGWAS --out_dir=$crossAncestryResults$mypheno/ --chrom=$chrom  --seed=42


}
export -f performPRSCS_custom # this makes local functions executable when bsubbed





# Performs PRS-CS for 1 population (EUR)
function performPRSCS { 
eursums=$1
Neur=$2
validBim=$3
mypheno=$4
chrom=$5
ncores=$6


export MKL_NUM_THREADS=$ncores
export NUMEXPR_NUM_THREADS=$ncores
export OMP_NUM_THREADS=$ncores

prscsxdir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsscript/PRScsx/"
prscsrefs="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsrefs/"
crossAncestryResults="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/"
python $prscsxdir/PRScsx.py --ref_dir=$prscsrefs --bim_prefix=$validBim --sst_file=$eursums --n_gwas=$Neur --pop=EUR --out_dir=$crossAncestryResults --out_name=$mypheno --chrom=$chrom  --seed=42

}
export -f performPRSCS # this makes local functions executable when bsubbed











#eursums=$eursumstatsLoc$EURPheno$'_QC_PRSCS'
#eassums=$japsumstatsLoc$JPTPheno$'_QC_PRSCS'
# Performs PRS-CSx for 2 pops on a given chrom
function performPRSCSx { 
eursums=$1
eassums=$2
Neur=$3
Neas=$4
validBim=$5
mypheno=$6
chrom=$7
ncores=$8


export MKL_NUM_THREADS=$ncores
export NUMEXPR_NUM_THREADS=$ncores
export OMP_NUM_THREADS=$ncores

prscsxdir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsscript/PRScsx/"
prscsrefs="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsrefs/"
crossAncestryResults="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/"
python $prscsxdir/PRScsx.py --ref_dir=$prscsrefs --bim_prefix=$validBim --sst_file=$eursums,$eassums --n_gwas=$Neur,$Neas --pop=EUR,EAS --out_dir=$crossAncestryResults --out_name=$mypheno --chrom=$chrom  --seed=42

}
export -f performPRSCSx # this makes local functions executable when bsubbed






#sumstatsLoc=$eursumstatsLoc$EURPheno$'_QC'
#isBinary=$binPheno

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




wc -l $crossAncestryResults$pheno$'_sumstats_meta_QC_keep'





# --chrom=CHROM 




