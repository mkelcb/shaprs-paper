
# screen -r -D 61726.pts-0.node-11-1-1

################

# EXTRACTING AND COMPRESSING THE SUPPLEMENTARY DATA:
asthaBaseLoc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/'

asthmaCrossAncestry=$asthaBaseLoc$'crossAncestry/'
crossAncestrySumStats=$asthmaCrossAncestry$'sumstats/'
crossAncestryResults=$asthmaCrossAncestry$'results/'

crossAncestryResults='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/DELETE_LATER/asthma/crossAncestry/results/'




# EUR JAP
generateCrossAncstryDiagFile $crossAncestryResults$'EUR_JAP_asthma'
generateCrossAncstryDiagFile $crossAncestryResults$'EUR_JAP_height'
generateCrossAncstryDiagFile $crossAncestryResults$'EUR_JAP_BRCA'
generateCrossAncstryDiagFile $crossAncestryResults$'EUR_JAP_CAD'
generateCrossAncstryDiagFile $crossAncestryResults$'EUR_JAP_T2D'

# rename to be more friendly
mv $crossAncestryResults$'EUR_JAP_asthma_diagDat_relab123456' $crossAncestryResults$'EUR_JAP_asthma_diagdata' 
mv $crossAncestryResults$'EUR_JAP_height_diagDat_relab123456' $crossAncestryResults$'EUR_JAP_height_diagdata' 
mv $crossAncestryResults$'EUR_JAP_asthma_diagDat_relab123456' $crossAncestryResults$'EUR_JAP_asthma_diagdata' 
mv $crossAncestryResults$'EUR_JAP_BRCA_diagDat_relab123456' $crossAncestryResults$'EUR_JAP_BRCA_diagdata' 
mv $crossAncestryResults$'EUR_JAP_CAD_diagDat_relab123456' $crossAncestryResults$'EUR_JAP_CAD_diagdata' 
mv $crossAncestryResults$'EUR_JAP_T2D_diagDat_relab123456' $crossAncestryResults$'EUR_JAP_T2D_diagdata'


# AFR - EUR
generateCrossAncstryDiagFile_AFR $revisionsResults$'AFR_EUR_BMI'
generateCrossAncstryDiagFile_AFR $revisionsResults$'AFR_EUR_height'
generateCrossAncstryDiagFile_AFR $revisionsResults$'AFR_EUR_LDL'


# rename to be more friendly
mv $revisionsResults$'AFR_EUR_BMI_diagDat_relab123456' $revisionsResults$'AFR_EUR_BMI_diagdata' 
mv $revisionsResults$'AFR_EUR_height_diagDat_relab123456' $revisionsResults$'AFR_EUR_height_diagdata' 
mv $revisionsResults$'AFR_EUR_LDL_diagDat_relab123456' $revisionsResults$'AFR_EUR_LDL_diagdata'

# copy original data into safe place
cp $crossAncestryResults$'EUR_JAP_asthma_diagdata'  $revisionsResults$'EUR_JAP_asthma_diagdata' 
cp $crossAncestryResults$'EUR_JAP_height_diagdata'  $revisionsResults$'EUR_JAP_height_diagdata' 
cp $crossAncestryResults$'EUR_JAP_asthma_diagdata'  $revisionsResults$'EUR_JAP_asthma_diagdata' 
cp $crossAncestryResults$'EUR_JAP_BRCA_diagdata'  $revisionsResults$'EUR_JAP_BRCA_diagdata' 
cp $crossAncestryResults$'EUR_JAP_CAD_diagdata'  $revisionsResults$'EUR_JAP_CAD_diagdata' 
cp $crossAncestryResults$'EUR_JAP_T2D_diagdata' $revisionsResults$'EUR_JAP_T2D_diagdata'


cp $crossAncestryResults$'crossAncestryLegend.txt' $revisionsResults$'crossAncestryLegend.txt'
cp $crossAncestryResults$'crossAncestryLegend.txt' $revisionsResults$'crossAncestryLegendAFR.txt'

# edit legend for AFR ones
# nano $revisionsResults$'crossAncestryLegendAFR.txt'

################################################
#  produce the IBD diag data ( executed on the IBD script variables)
crossAncestryLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/DELETE_LATER/asthma/crossAncestry/results/"

# SNP lFDR Beta_UC SE_UC Beta_CD SE_CD Beta_ShaPRS_CD SE_ShaPRS_CD Beta_ShaPRS_UC SE_ShaPRS_UC Meta_LDPred2 CD_LDpred2  UC_LDpred2  CD_shaPRS_LDpred2  UC_shaPRS_LDpred2
# SNP lFDR Beta_A  SE_A  Beta_B  SE_B  Beta_ShaPRS    SE_ShaPRS

# screen -r -D  22289.pts-3.node-11-1-1

$newScratchLoc$'gwas3_cd/CD_SHAPRS_het'


$newResultsLoc$'manhattan_lFDR_diagDat'


# rename this to more friendly
awk '{if (FNR == 1) {print "SNP lFDR Beta_UC SE_UC Beta_CD SE_CD Beta_ShaPRS_CD SE_ShaPRS_CD"}  else {print $0} }' $newResultsLoc$'manhattan_lFDR_diagDat' > $newResultsLoc$'manhattan_lFDR_diagDat2'
head $newResultsLoc$'manhattan_lFDR_diagDat2'

# add UC results for shaPRS
awk 'FNR == NR { file1[ $3 ] = $7" "$8; next; } 
FNR <= NR { if( FNR == 1 ) {print $0" Beta_ShaPRS_UC SE_ShaPRS_UC" } 
else if( $1 in file1) {print $0" "file1[ $1 ] } }
' OFS=" " FS="\t" $newScratchLoc$'gwas3_uc/UC_SHAPRS_het_shaprs' FS=" " $newResultsLoc$'manhattan_lFDR_diagDat2' > $newResultsLoc$'manhattan_lFDR_diagDat3'
head $newResultsLoc$'manhattan_lFDR_diagDat3'

# add CD results for MTAG
awk 'FNR == NR { file1[ $3 ] = $7" "$8; next; } 
FNR <= NR { if( FNR == 1 ) {print $0" Beta_MTAG_CD SE_MTAG_CD" } 
else if( $1 in file1) {print $0" "file1[ $1 ] } }
' OFS=" " FS="\t" $newScratchLoc$'gwas3_cd/CD_MTAG_het_mtag' FS=" " $newResultsLoc$'manhattan_lFDR_diagDat3' > $newResultsLoc$'manhattan_lFDR_diagDat4'
head $newResultsLoc$'manhattan_lFDR_diagDat4'

# add UC results for MTAG
awk 'FNR == NR { file1[ $3 ] = $7" "$8; next; } 
FNR <= NR { if( FNR == 1 ) {print $0" Beta_MTAG_UC SE_MTAG_UC" } 
else if( $1 in file1) {print $0" "file1[ $1 ] } }
' OFS=" " FS="\t" $newScratchLoc$'gwas3_uc/UC_MTAG_het_mtag' FS=" " $newResultsLoc$'manhattan_lFDR_diagDat4' > $newResultsLoc$'manhattan_lFDR_diagDat5'
head $newResultsLoc$'manhattan_lFDR_diagDat5'



# Add 'Meta PRS-CS'
AddPRSCSToDiag $newResultsLoc$'manhattan_lFDR_diagDat5' $crossAncestryLoc$'gwas3_cd_meta_het_PRSCS_no_dupes' 'Meta_PRSCS' '1'

# Add 'CD_PRSCS'
AddPRSCSToDiag $newResultsLoc$'manhattan_lFDR_diagDat51' $crossAncestryLoc$'gwas3_cd_primary_het_PRSCS_no_dupes' 'CD_PRSCS' '1'

# Add 'UC_PRSCS'
AddPRSCSToDiag $newResultsLoc$'manhattan_lFDR_diagDat511' $crossAncestryLoc$'gwas3_uc_primary_het_PRSCS_no_dupes' 'UC_PRSCS' '1'





# Add 'CD_shaPRS_PRSCS'
AddPRSCSToDiag $newResultsLoc$'manhattan_lFDR_diagDat5111' $crossAncestryLoc$'gwas3_cd_shaprs_het_PRSCS_no_dupes' 'CD_shaPRS_PRSCS' '1'


# Add 'UC_shaPRS_PRSCS'
AddPRSCSToDiag $newResultsLoc$'manhattan_lFDR_diagDat51111' $crossAncestryLoc$'gwas3_uc_mtag_het_PRSCS_no_dupes' 'UC_shaPRS_PRSCS' '1'


# Add 'CD_MTAG_PRSCS'
AddPRSCSToDiag $newResultsLoc$'manhattan_lFDR_diagDat511111' $crossAncestryLoc$'gwas3_cd_mtag_het_PRSCS_no_dupes' 'CD_MTAG_PRSCS' '1'

# Add 'UC_MTAG_PRSCS'
AddPRSCSToDiag $newResultsLoc$'manhattan_lFDR_diagDat5111111' $crossAncestryLoc$'gwas3_uc_shaprs_het_PRSCS_no_dupes' 'UC_MTAG_PRSCS' '1'

mv $newResultsLoc$'manhattan_lFDR_diagDat51111111'

grep "rs1000000" $crossAncestryLoc$'gwas3_cd_meta_het_PRSCS_no_dupes'

# rename to be more friendly
mv $newResultsLoc$'manhattan_lFDR_diagDat51111111' $newResultsLoc$'CD_UC_diagdata'
##############################




# zip it up (the legend files were separately uploaded)

# EUR - EAS (this is already OK)
zip -j $revisionsResults$'_shaPRS_asthma_suppData.zip' $revisionsResults$'crossAncestryLegend.txt' $revisionsResults$'EUR_JAP_asthma_diagdata'

zip -j $revisionsResults$'_shaPRS_height_suppData.zip' $revisionsResults$'crossAncestryLegend.txt' $revisionsResults$'EUR_JAP_height_diagdata'

zip -j $revisionsResults$'_shaPRS_BRCA_suppData.zip' $revisionsResults$'crossAncestryLegend.txt' $revisionsResults$'EUR_JAP_BRCA_diagdata'

zip -j $revisionsResults$'_shaPRS_CAD_suppData.zip' $revisionsResults$'crossAncestryLegend.txt' $revisionsResults$'EUR_JAP_CAD_diagdata'

zip -j $revisionsResults$'_shaPRS_T2D_suppData.zip' $revisionsResults$'crossAncestryLegend.txt' $revisionsResults$'EUR_JAP_T2D_diagdata'

# AFR - EUR
zip -j $revisionsResults$'_shaPRS_AFR_BMI_suppData.zip' $revisionsResults$'crossAncestryLegendAFR.txt' $revisionsResults$'AFR_EUR_BMI_diagdata'

zip -j $revisionsResults$'_shaPRS_AFR_height_suppData.zip' $revisionsResults$'crossAncestryLegendAFR.txt' $revisionsResults$'AFR_EUR_height_diagdata'

zip -j $revisionsResults$'_shaPRS_AFR_LDL_suppData.zip' $revisionsResults$'crossAncestryLegendAFR.txt' $revisionsResults$'AFR_EUR_LDL_diagdata'

# IBD
zip -j '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/revision/results/shaPRS_IBD_suppData.zip' '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/newWave/CD_UC_diagdata' '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/revision/results/IBDlegend.txt'

# zip all of it up into a single .zip 
zip -j $revisionsResults$'All_shaPRS_suppData.zip' '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/newWave/CD_UC_diagdata' '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/revision/results/IBDlegend.txt' $revisionsResults$'crossAncestryLegend.txt' $revisionsResults$'crossAncestryLegendAFR.txt' $revisionsResults$'EUR_JAP_asthma_diagdata' $revisionsResults$'EUR_JAP_height_diagdata' $revisionsResults$'EUR_JAP_asthma_diagdata' $revisionsResults$'EUR_JAP_BRCA_diagdata' $revisionsResults$'EUR_JAP_CAD_diagdata' $revisionsResults$'EUR_JAP_T2D_diagdata' $revisionsResults$'_shaPRS_AFR_BMI_suppData.zip' $revisionsResults$'_shaPRS_AFR_height_suppData.zip' $revisionsResults$'_shaPRS_AFR_LDL_suppData.zip'


###########################################


filest=$revisionsResults$'AFR_EUR_BMI'
function generateCrossAncstryDiagFile_AFR {
filest=$1

# convert change the column labels to be explicit:
awk '{
if (FNR == 1) {print "SNP lFDR Beta_EUR SE_EUR Beta_AFR SE_AFR Beta_ShaPRS SE_ShaPRS" }
else { print $0}
}' $filest$'_diagDat' > $filest$'_diagDat_relab'
head $filest$'_diagDat_relab'

# add each 
AddPRSCSToDiag $filest$'_diagDat_relab' $filest$'_AFR_AFR_PRSCS' 'PRS-CS' '1'

AddPRSCSToDiag $filest$'_diagDat_relab1' $filest$'_shaPRS_AFR_PRSCS' 'shaPRS+PRS-CS' '2'

AddPRSCSToDiag $filest$'_diagDat_relab12' $filest$'_AFR_PRSCSx' 'AFR_PRSCSx-stage1' '3'

AddPRSCSToDiag $filest$'_diagDat_relab123' $filest$'_EUR_PRSCSx' 'EUR_PRSCSx-stage1' '4'

AddLDpredToDiag $filest$'_diagDat_relab1234' $filest$'_AFR_AFR_LDpred2_CSxsubset' 'LDpred2' '5'

AddLDpredToDiag $filest$'_diagDat_relab12345' $filest$'_shaPRS_LDpred2_CSxsubset' 'shaPRS+LDpred2' '6'

head $filest$'_diagDat_relab123456'
}


filest=$crossAncestryResults$pheno
function generateCrossAncstryDiagFile {
filest=$1

# convert change the column labels to be explicit:
awk '{
if (FNR == 1) {print "SNP lFDR Beta_EAS SE_EAS Beta_EUR SE_EUR Beta_ShaPRS SE_ShaPRS" }
else { print $0}
}' $filest$'_diagDat' > $filest$'_diagDat_relab'
head $filest$'_diagDat_relab'

# add each 
AddPRSCSToDiag $filest$'_diagDat_relab' $filest$'_EUR_EUR_PRSCS' 'PRS-CS' '1'

AddPRSCSToDiag $filest$'_diagDat_relab1' $filest$'_shaPRS_EUR_PRSCS' 'shaPRS+PRS-CS' '2'

AddPRSCSToDiag $filest$'_diagDat_relab12' $filest$'_EUR_PRSCSx' 'EUR_PRSCSx-stage1' '3'

AddPRSCSToDiag $filest$'_diagDat_relab123' $filest$'_EAS_PRSCSx' 'EAS_PRSCSx-stage1' '4'

AddLDpredToDiag $filest$'_diagDat_relab1234' $filest$'_EUR_EUR_LDpred2_CSxsubset' 'LDpred2' '5'

AddLDpredToDiag $filest$'_diagDat_relab12345' $filest$'_shaPRS_LDpred2_CSxsubset' 'shaPRS+LDpred2' '6'

head $filest$'_diagDat_relab123456'
}


#diagDatSoFar=$crossAncestryResults$pheno$'_diagDat_relab' 
#PRSFile=$crossAncestryResults$pheno$'_EUR_EUR_LDpred2_CSxsubset' 
#colname='_EUR_EUR_LDpred2_CSxsubset'
#outname="1"
function AddLDpredToDiag {
diagDatSoFar=$1
PRSFile=$2 # assumes that weight is 3rd column, and SNP is 1st col
colname=$3
outname=$4
awk -v colname="$colname" 'FNR == NR { file1[ $1 ] = $3; next; } 
FNR <= NR { 
if( FNR == 1 ) {print $0" "colname } 
else if( $1 in file1) {print $0" "file1[ $1 ] } 
else {print $0" 0" } 
}
' OFS=" "  FS="\t" $PRSFile FS=" " $diagDatSoFar > $diagDatSoFar$outname

head $diagDatSoFar$outname
wc -l $diagDatSoFar$outname
wc -l $diagDatSoFar

}


function AddPRSCSToDiag {
diagDatSoFar=$1
PRSFile=$2 # assumes that weight is 6th column, and SNP is 2nd col
colname=$3
outname=$4
awk -v colname="$colname" 'FNR == NR { file1[ $2 ] = $6; next; } 
FNR <= NR { 
if( FNR == 1 ) {print $0" "colname } 
else if( $1 in file1) {print $0" "file1[ $1 ] } 
else {print $0" 0" } 
}
' OFS=" "  FS="\t" $PRSFile FS=" " $diagDatSoFar > $diagDatSoFar$outname

head $diagDatSoFar$outname
wc -l $diagDatSoFar$outname
wc -l $diagDatSoFar

}

#################################################




