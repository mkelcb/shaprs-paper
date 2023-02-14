####################################################
#
# shaPRS revisions AFR

# screen -r -D 38802.pts-0.node-11-1-1  # generating the missing lFDR SNP infos


#########################
# static vars
largePlinkMem=200000
bigMem=400000
largeMem=60000
plinkMem=6000  # how much RAM PLINK would use
largerMem=12000 
homeBase='/nfs/users/nfs_m/mk23/'
plink=$homeBase$'software/plink/plink'
plink2=$homeBase$'software/plink/plink2'
scriptsLoc=$homeBase$'scripts/'
#Rscript='/nfs/team152/mk23/software/R/R-3.6.1/bin/Rscript'
Rscript='/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript'
flashpca='/nfs/users/nfs_m/mk23/software/flashpca/flashpca_x86-64'
qctool2_new='/nfs/team152/mk23/software/qctool2/qctool_v2.0.6-Ubuntu16.04-x86_64/qctool'
qctoolmem=10000
shaPRS="_shaPRS"
python2='python2'
ldsc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/0_software/ldsc/ldsc.py'
UKBBQCFile='/lustre/scratch115/realdata/mdt3/projects/ukbiobank/FullRelease/Imputed/001/ukb_sqc_v2.txt'
asthmaRawPhenos='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/data/raw/pheno_raw.txt'
revisionDir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/revisions/"
hapmap3_b37bim='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/0_common/data/twas/raw/hapmap3_r1_b37_fwd_consensus.qc.poly.recode.bim' # this is the full hapmap3, that has the HLA
baseLDpredLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/"
AFRRef="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/AFRRef/ld-ref/"
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




# Thesis location
baseLoc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/'
commonScratchLoc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/0_common/data/scratch/'
commonDataLoc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/0_common/data/'
ukbbQCkeeplist_publicDs=$commonDataLoc$'ukbb_hq_keeplist.txt_pubIDs'
thesisAsthmaDataLoc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/data/'



shaprsBaseloc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/'
revisions=$shaprsBaseloc$'revision/'
mtagscript=$revisions$'mtag/'
revisionsScratch=$revisions$'scratch/'
revisionsRaw=$revisions$'raw/'
revisionsResults=$revisions$'results/'

asthaBaseLoc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/'

asthmaCrossAncestry=$asthaBaseLoc$'crossAncestry/'
crossAncestrySumStats=$asthmaCrossAncestry$'sumstats/'
crossAncestryResults=$asthmaCrossAncestry$'results/'
crossAncestryRaw=$asthmaCrossAncestry$'raw/'
prscssorig=$crossAncestryRaw$'prscss/'

shaPRSscriptLoc="/nfs/users/nfs_m/mk23/scripts/shaPRS.R"
prscsrefs=$crossAncestryRaw$'prscsrefs/'
prscsscript=$crossAncestryRaw$'prscsscript/'
crossAncestryScratch=$asthmaCrossAncestry$'scratch/'

prscsxdir=$prscsscript$'PRScsx/'


mkdir -p $mtagscript
mkdir -p $revisionsRaw
mkdir -p $revisionsScratch
mkdir -p $revisionsResults


##############################

# PAGE:

# https://www.ebi.ac.uk/gwas/studies/GCST008048
# https://www.ebi.ac.uk/gwas/studies/GCST008053

# screen -r -D 30405.pts-0.node-11-1-1

cd $revisionsRaw
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008048/WojcikG_PMID_T2D_bmi.gz


gunzip WojcikG_PMID_T2D_bmi.gz

head WojcikG_PMID_T2D_bmi
# Chr     Position_hg19   SNP     Other-allele    Effect-allele   Effect-allele-frequency Sample-size     Effect-allele-frequency-cases       Sample-size-cases       Beta    SE      P-val   INFO-score      rsid
#1       10616   rs376342519:10616:CCGCCGTTGCAAAGGCGCGCCG:C      CCGCCGTTGCAAAGGCGCGCCG  C       0.9913283       45725       0.9891043       14042   0.01864849      0.1011172       0.8536804       0.604   rs376342519
#1       10642   1:10642:G:A     G       A       0.005889623     45725   0.007325986     14042   0.2478162       0.1482746   0.09465601      0.441   rs558604819

wc -l WojcikG_PMID_T2D_bmi # 23997820

wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008053/WojcikG_PMID_invn_rheight_alls.gz


gunzip WojcikG_PMID_invn_rheight_alls.gz

head WojcikG_PMID_invn_rheight_alls
wc -l WojcikG_PMID_invn_rheight_alls # 34656551
# Chr     Position_hg19   SNP     Other-allele    Effect-allele   Effect-allele-frequency Sample-size     Beta    SE P-val    INFO-score      rsid
# 1       10539   1:10539:C:A     C       A       0.004370543     49781   0.02978985      0.06364248      0.6397265  0.46     rs537182016



wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST002001-GCST003000/GCST002783/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz

gunzip All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz
# These studies do NOT make available just AFR indi association stats, so they cannot be used

############################################

# download the 3 AFR studies:

# 1. Height https://www.ebi.ac.uk/gwas/studies/GCST009055
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009055/heightannotated.txt.gz

gunzip heightannotated.txt.gz

head -n 2 heightannotated.txt
wc -l heightannotated.txt # 29,270,441

# in their paper in Table 1 they say "p_assoc, p value from RE2", IE they consider "RE2" to be the 'best" meta-analysis so I will be using that
# "mmc6.xlsx" says that
#pval_RE	random-effects meta-analysis
#beta_RE	beta coefficient for random effects meta-analysis
#pval_RE2	p value for Han-Eskin random effects meta-analysis

#  1                     2               3          4            5             6             7           8             9            10          11             12          13        14         15           16         17             18    19        20           21             22             23           24           25             26         27    28          29          30        31         32         33                                                                                                                                                                                     
# snpid             beta_uganda     se_uganda    beta_DCC      se_DCC      beta_DDS       se_DDS     beta_AADM      se_AADM      studyno     pval_fe        beta_fe      se_fe    pval_re     beta_re      se_re     pval_re2       Isquare  Q      pval_Q      tau_square     pval_uganda     pval_DCC     pval_DDS    pval_AADM     af_uganda    af_DCC af_DDS     af_AADM    no_uganda   no_DCC     no_DDS     no_AADM
# 10:100000012:G:A 8.343710e-03 5.639773e-02 -1.733857e-02 9.494213e-02 -1.676534e-02 1.158643e-01 -0.0283 0.06371 4 0.783054 -0.0100797 0.0366079 0.783054 -0.0100797 0.0366079 0.832427 0.00000 0.197677 0.977963 0.00000 0.882387 0.855094 0.884949 0.656898 0.972764 0.960781 0.965522 0.97619 6347 1478 1114 5187

# this SNP is meant to have p val of: 2 x 10-20
grep "rs6787063" heightannotated.txt # not present

grep "3:16594958" heightannotated.txt # not present, this is Build38

grep "3:16636465" heightannotated.txt # this is Build37
# 3:16636465:C:G 1.548675e-02 1.895958e-02 -4.297945e-01 4.226201e-02 8.538201e-03 4.945418e-02 0.02878 0.02055 4 0.107468 -0.0205784 0.0127842 0.266057 -0.0918249 0.0825624 2.17637E-20 97.1012 103.491 2.75850E-22 0.0260419 0.414026 2.70574E-24 0.862928 0.161368 0.409756 0.303451 0.325446 0.367746 6347 1478 1114 5187

grep "16594958" heightannotated.txt

# TO:
#  1     2     3     4    5         6        7    8   9  10
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N



awk 'FNR == NR { file1[ $1":"$4":"$5":"$6 ] = $2; file2[ $1":"$4":"$6":"$5 ] = $2; next; } 
{ 
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
if ( $1 in file1 || $1 in file2 ) {
if($1 in file1) {rsid=file1[$1]}
else {rsid=file2[$1]}
split($1,variant,":")  # split 10:100000012:G:A  -> 10, 100000012, G, A
N=$30+$31+$32+$33
AF=$26 * ($30/N) + $27 * ($31/N) + $28 * ($32/N)  + $29 * ($33/N) 
print variant[1]"\t"variant[2]"\t"rsid"\t"variant[3]"\t"variant[4]"\t"AF"\t"$15"\t"$16"\t"$17"\t"N}}   ' OFS='\t' FS="\t" $hapmap3_b37bim FS=" " heightannotated.txt > $revisionsScratch$'height_AFR_hm3'
head $revisionsScratch$'height_AFR_hm3'
wc -l $revisionsScratch$'height_AFR_hm3' # 1,361,687 # so we have all the HM3 vars!

########
# BMI: https://www.ebi.ac.uk/gwas/studies/GCST009057
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009057/bmiannotated.txt.gz

gunzip bmiannotated.txt.gz

head -n 2 bmiannotated.txt
wc -l bmiannotated.txt # 29,280,511
# snpid beta_uganda se_uganda beta_DCC se_DCC beta_DDS se_DDS beta_AADM se_AADM studyno pval_fe beta_fe se_fe pval_re beta_re se_re pval_re2 Isquare Q pval_Q tau_square pval_uganda pval_DCC pval_DDS pval_AADM af_uganda af_DCC af_DDS af_AADM no_uganda no_DCC no_DDS no_AADM
#10:100000012:G:A 6.955139e-02 5.742565e-02 -1.767886e-01 9.512398e-02 -5.907292e-03 1.166917e-01 -0.4927 0.3718 4 0.912445 -0.00494496 0.0449727 0.555139 -0.0485512 0.0822794 0.458041 55.0057 6.66751 0.0832854 0.0134475 0.225836 0.0630970 0.959626 0.185113 0.972764 0.960781 0.965522 0.97619 6197 1478 1114 5187


awk 'FNR == NR { file1[ $1":"$4":"$5":"$6 ] = $2; file2[ $1":"$4":"$6":"$5 ] = $2; next; } 
{ 
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
if ( $1 in file1 || $1 in file2 ) {
if($1 in file1) {rsid=file1[$1]}
else {rsid=file2[$1]}
split($1,variant,":")  # split 10:100000012:G:A  -> 10, 100000012, G, A
N=$30+$31+$32+$33
AF=$26 * ($30/N) + $27 * ($31/N) + $28 * ($32/N)  + $29 * ($33/N) 
print variant[1]"\t"variant[2]"\t"rsid"\t"variant[3]"\t"variant[4]"\t"AF"\t"$15"\t"$16"\t"$17"\t"N}}   ' OFS='\t' FS="\t" $hapmap3_b37bim FS=" " bmiannotated.txt > $revisionsScratch$'bmi_AFR_hm3'
head $revisionsScratch$'bmi_AFR_hm3'
wc -l $revisionsScratch$'bmi_AFR_hm3' # 1361687

##################
# LDL: https://www.ebi.ac.uk/gwas/studies/GCST009043

wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009043/lowdlipoannotated.txt.gz



gunzip lowdlipoannotated.txt.gz

head -n 2 lowdlipoannotated.txt
wc -l lowdlipoannotated.txt # 29322937
#snpid beta_uganda se_uganda beta_DCC se_DCC beta_DDS se_DDS beta_AADM se_AADM studyno pval_fe beta_fe se_fe pval_re beta_re se_re pval_re2 Isquare Q pval_Q tau_square pval_uganda pval_DCC pval_DDS pval_AADM af_uganda af_DCC af_DDS af_AADM no_uganda no_DCC no_DDS no_AADM
#10:100000012:G:A 2.958460e-02 5.586749e-02 4.300130e-02 9.580216e-02 1.209279e-02 1.168602e-01 0.05586 0.07271 4 0.330083 0.0370312 0.0380218 0.330083 0.0370312 0.0380218 0.388228 0.00000 0.134250 0.987432 0.00000 0.596424 0.653536 0.917581 0.442334 0.972764 0.960781 0.965522 0.97619 6407 1475 1118 4086


awk 'FNR == NR { file1[ $1":"$4":"$5":"$6 ] = $2; file2[ $1":"$4":"$6":"$5 ] = $2; next; } 
{ 
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
if ( $1 in file1 || $1 in file2 ) {
if($1 in file1) {rsid=file1[$1]}
else {rsid=file2[$1]}
split($1,variant,":")  # split 10:100000012:G:A  -> 10, 100000012, G, A
N=$30+$31+$32+$33
AF=$26 * ($30/N) + $27 * ($31/N) + $28 * ($32/N)  + $29 * ($33/N) 
print variant[1]"\t"variant[2]"\t"rsid"\t"variant[3]"\t"variant[4]"\t"AF"\t"$15"\t"$16"\t"$17"\t"N}}   ' OFS='\t' FS="\t" $hapmap3_b37bim FS=" " lowdlipoannotated.txt > $revisionsScratch$'LDL_AFR_hm3'
head $revisionsScratch$'LDL_AFR_hm3'
wc -l $revisionsScratch$'LDL_AFR_hm3' # 1362368



########################





##################################################################
# BLOOD LIPIDS: Willer CJ et al. Discovery and refinement of loci associated with lipid levels
# http://csg.sph.umich.edu/willer/public/lipids2013/
# The tables below summarize the results of the meta-analysis described in Willer et al (2013). 
# Each table includes ten columns: the marker names in build hg18 and hg19, the marker names in rsid format, the two alleles of the SNP, 
# the effect size associated with the SNP, the corresponding standard error, the number of individuals evaluated for the SNP, 
# the combined p-value for the SNP and the frequency of the first allele (the trait increasing allele) in the 1000 Genomes European sample. 
# The files include both genotyped and imputed SNPs and are provided in a compressed text file to conserve space.

wget http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_LDL.txt.gz
gunzip jointGwasMc_LDL.txt.gz
head jointGwasMc_LDL.txt
wc -l jointGwasMc_LDL.txt # 2,437,752
# SNP_hg18        SNP_hg19        rsid    A1      A2      beta    se      N       P-value Freq.A1.1000G.EUR
# chr10:10000135  chr10:9960129   rs4747841       a       g       0.0037  0.0052  89138.00        0.7158  0.4908


#  EUR LDL :
awk 'FNR == NR { file1[ $2 ] = $2;  next; } 
{
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else if($3 in file1) { 
split($2,chrPos,":")  # split chr10:9960129 -> chr10 9960129
split(chrPos[1],chr,"chr") # split chr10  -> 10
print chr[2]"\t"chrPos[2]"\t"$3"\t"toupper($4)"\t"toupper($5)"\t"$10"\t"$6"\t"$7"\t"$9"\t"$8}
} ' OFS='\t' FS="\t" $hapmap3_b37bim jointGwasMc_LDL.txt > $revisionsScratch$'LDL_EUR_hm3'
head $revisionsScratch$'LDL_EUR_hm3'
wc -l $revisionsScratch$'LDL_EUR_hm3' # 1072701




# EUR Height: already done this for the EUR-JP cross ancestry stuff
# https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz
cp '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/sumstats/eur/EUR_height_hm3' $revisionsScratch$'height_EUR_hm3'
$revisionsScratch$'height_EUR_hm3'

# EUR BMI: GWAS Anthropometric 2015 BMI Summary Statistics
wget https://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz

gunzip SNP_gwas_mc_merge_nogc.tbl.uniq.gz
head All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq
wc -l head All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq # 2555511

awk 'FNR == NR { file1[ $2 ] = $1"\t"$4;  next; }
{
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else if($1 in file1) { 
print file1[$1]"\t"$0} 
} ' OFS='\t' FS="\t" $hapmap3_b37bim FS="\t"  All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq > $revisionsScratch$'bmi_EUR_hm3'

head $revisionsScratch$'bmi_EUR_hm3'
wc -l $revisionsScratch$'bmi_EUR_hm3' # 1130678




########################


# get the BMI/height phenos for the AFR indis
'/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/height/data/raw/pheno_raw.txt'
'/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/bmi/data/raw/pheno_raw.txt' # 502,536 indis

grep "5172879" '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/bmi/data/raw/pheno_raw.txt'
# not found
grep "1587700" '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/bmi/data/raw/pheno_raw.txt'
# ok this is this format



baseLoc_UKBB="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/rarePGS/pcaUKBB"
dataLoc_UKBB=$baseLoc_UKBB$'data/scratch/'
resultsLoc_UKBB=$baseLoc_UKBB$'data/'
scratchLoc_UKBB=$baseLoc_UKBB$'scratch/'


head $scratchLoc_UKBB$'UKBBIDs_PC16'
wc -l $scratchLoc_UKBB$'UKBBIDs_PC16' # 
 
# those excluded due to relatedness or sex discordant or withdrawn consent
# /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbbQCexcludelistt.txt
# 40614
# 2445446	1410461

# this can be used to map between
# Team17670 PublicID
# 2791700   5172879  ...
# /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt_oldids
# wc -l 451837
grep "2791700" /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt_oldids
# 5172879 2791700

# map between
grep "1587700" /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt
# 1587700 -> 5172879

grep "5172879" "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/customLDRef/QC_v2.csv"

# this has the '2791700' IDs and the 16 PCs that we could use to identify AFR indis:
# "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/customLDRef/QC_v2.csv"

# this has the '5172879' IDs and the 16 PCs that we could use to identify AFR indis:
grep "5172879" "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/customLDRef/QC.csv"



# phenotypes are in "1587700" format, but QC file is in "5172879" format.
# we need to map the phenos to "5172879"


awk 'FNR == NR { file1[ $1 ] = $2;  next; }
{
if($1 in file1)  {
print file1[$1]"\t"$2 
} }' OFS='\t' FS=" " /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt FS=" " '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/bmi/data/raw/pheno_raw.txt' > $revisionsScratch$'UKBB_BMI_phenos'
head $revisionsScratch$'UKBB_BMI_phenos'
wc -l $revisionsScratch$'UKBB_BMI_phenos'


awk 'FNR == NR { file1[ $1 ] = $2;  next; }
{
if($1 in file1)  {
print file1[$1]"\t"$2 
} }' OFS='\t' FS=" " /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt FS=" " '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/height/data/raw/pheno_raw.txt' > $revisionsScratch$'UKBB_height_phenos'
head $revisionsScratch$'UKBB_height_phenos'
wc -l $revisionsScratch$'UKBB_height_phenos'


#/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt
awk 'FNR == NR { file1[ $1 ] = $2;  next; }
{
if($1 in file1)  {
print file1[$1]"\t"$2 
} }' OFS='\t' FS=" " /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt FS=" " '/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/height/data/raw/pheno_raw.txt' > $revisionsScratch$'UKBB_height_phenos'
head $revisionsScratch$'UKBB_height_phenos'

###################################
# extract LDL phenos (30780)
ukbbphenos='/lustre/scratch123/hgi/mdt1/projects/ukbb_team152/45669/phenotype_data/release_20220704/raw_phenotype/ukb52241/ukb52241.tab' # stores pheno data
arguments='/nfs/users/nfs_m/mk23/scripts/ukbbPhenoExtract.R '$ukbbphenos$' f.30780.0.0 '$revisionsScratch$'LDL_raw.txt'
$Rscript $arguments


# remove withdrawn samples (no need, as the new raw pheno came from a release that already excluded them
awk 'FNR == NR { file1[ $1 ] = $1;  next; }
{if($1 in file1 == 0)  {
print $0 } }' OFS='\t' /lustre/scratch123/hgi/mdt1/projects/ukbb_team152/45669/withdrawn_samples_20220222/w45669_20220222.csv $revisionsScratch$'LDL_raw.txt' > $revisionsScratch$'LDL_raw2.txt'
head $revisionsScratch$'LDL_raw2.txt'
wc -l $revisionsScratch$'LDL_raw2.txt'
wc -l $revisionsScratch$'LDL_raw.txt'

- ukb45669_cal_chr5_v2_s488282.fam: file used to match application 45669 eid with genotyped files
- ukb45669_imp_chr5_v3_s487314.sample: file used to match application 45669 eid with imputed files



# /lustre/scratch123/hgi/mdt1/projects/ukbb_team152/45669/phenotype_data/release_20220704/raw_phenotype/ukb52241 
awk 'FNR == NR { file1[ $1 ] = $2;  next; }
{
if($1 in file1)  {
print file1[$1]"\t"$2 
} }' OFS='\t' FS=" " /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt FS=" " $revisionsScratch$'LDL_raw.txt' > $revisionsScratch$'UKBB_LDL_phenos'
head $revisionsScratch$'UKBB_LDL_phenos'
wc -l $revisionsScratch$'UKBB_LDL_phenos' # 451660


##############################


# identify AFR indis via Prive's UKBB pop centers
arguments='/nfs/users/nfs_m/mk23/scripts/UKBB_findAFR.R /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/customLDRef/QC.csv /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt '$revisionsScratch$'UKBB_AFR_indis'
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments
# "written  6414  AFR indis to: $revisionsScratch$'UKBB_AFR_indis'


################

# get AFR indis phe files

awk 'FNR == NR { file1[ $1 ] = $1;  next; }
{
if($1 in file1)  {
print $1"\t"$1"\t"$2 
} }' OFS='\t' $revisionsScratch$'UKBB_AFR_indis' $revisionsScratch$'UKBB_BMI_phenos' > $revisionsRaw$'UKBB_BMI_AFR.phe'
head $revisionsRaw$'UKBB_BMI_AFR.phe'
wc -l $revisionsRaw$'UKBB_BMI_AFR.phe' # 6411


awk 'FNR == NR { file1[ $1 ] = $1;  next; }
{
if($1 in file1)  {
print $1"\t"$1"\t"$2 
} }' OFS='\t' $revisionsScratch$'UKBB_AFR_indis' $revisionsScratch$'UKBB_height_phenos' > $revisionsRaw$'UKBB_height_AFR.phe'
head $revisionsRaw$'UKBB_height_AFR.phe'
wc -l $revisionsRaw$'UKBB_height_AFR.phe' # 6411
# 4147832

# create an AFR indis keeplist for the TeamID17670
awk 'FNR == NR { file1[ $1 ] = $2;  next; }
{
if($1 in file1)  {
print file1[$1] 
} }' OFS='\t' FS=" " /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt_oldids FS=" " $revisionsScratch$'UKBB_AFR_indis' > $revisionsScratch$'UKBB_AFR_indis_TeamID17670'
head $revisionsScratch$'UKBB_AFR_indis_TeamID17670'
wc -l $revisionsScratch$'UKBB_AFR_indis_TeamID17670' # 6414


# prepare LDL pheno file
awk 'FNR == NR { file1[ $1 ] = $1;  next; }
{
if($1 in file1)  {
print $1"\t"$1"\t"$2 
} }' OFS='\t' $revisionsScratch$'UKBB_AFR_indis' $revisionsScratch$'UKBB_LDL_phenos' > $revisionsRaw$'UKBB_LDL_AFR.phe'
head $revisionsRaw$'UKBB_LDL_AFR.phe'
wc -l $revisionsRaw$'UKBB_LDL_AFR.phe' # 6409


######

# Generate the LDpred2 LDref panel for the AFR population: use this keeplist $revisionsScratch$'UKBB_AFR_indis_TeamID17670'


/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/R
tempDataLoca="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/customLDRef/data/scractch/tmp-data"
refLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/AFRRef/"



map_hapmap3 <- bigreadr::fread2("/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/customLDRef/data/raw/ldblk_1kg_eur/snpinfo_1kg_hm3")

#ind_sub_csv = readRDS(file = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/customLDRef/data/ind_sub_csv.rds" ) # these are the indices
#str(ind_sub_csv)

library(doParallel)
cl <- makeCluster(22)
parallel::clusterExport(cl, "map_hapmap3")
list_snp_id <- parLapply(cl, 1:22, function(chr) {
  mfi <- paste0("/lustre/scratch123/hgi/mdt1/projects/ukbiobank_genotypes/FullRelease/Imputed/EGAD00010001474/ukb_mfi_chr", chr, "_v3.txt")
  infos_chr <- bigreadr::fread2(mfi, showProgress = FALSE)
  joined <- dplyr::inner_join(cbind(chr = chr, infos_chr), map_hapmap3[1:2],
                              by = c("chr" = "CHR", "V2" = "SNP"))
  with(joined[!vctrs::vec_duplicate_detect(joined$V2), ],
       paste(chr, V3, V4, V5, sep = "_"))
})
stopCluster(cl)



sum(lengths(list_snp_id))  # 1,117,493


sample <- bigreadr::fread2("/lustre/scratch123/hgi/mdt1/projects/ukbb_team152/17670/ibs/data/raw/ukb17670_imp_chr1_v3_s487395.sample")
str(sample)
sample <- sample[-1, ]

csv <- "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/revision/scratch/UKBB_AFR_indis_TeamID17670"

df0 <- bigreadr::fread2(csv)
# head(df0)
#       V1
#1 1451051

# match the desired AFR individuals to the ones we have in the samples
ind.indiv <- match(df0$V1, sample$ID_2)
sub <- which(!is.na(ind.indiv) ) # get the indices of the indices which are OK
head(sub)
head(ind.indiv)

sub2 = sub # for compatibility
head(sample$ID_2[ind.indiv[sub2]]  ) # this maps back to the correct actual sample IDs

refLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/AFRRef/"
NCORES <- 32
#install.packages("RSQLite")
#install.packages("dbplyr")
library(dbplyr)
library(RSQLite)
system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles   = glue::glue("/lustre/scratch123/hgi/mdt1/projects/ukbiobank_genotypes/FullRelease/Imputed/EGAD00010001474/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = paste0(refLoc,"/UKBB_imp_HM3"),
    ind_row     = ind.indiv[sub2],
    ncores      = NCORES
  )
) # 43 min with 23 cores




library(bigsnpr)
ukb <- snp_attach(paste0(refLoc,"UKBB_imp_HM3.rds") )
G <- ukb$genotypes
CHR <- as.numeric(ukb$map$chromosome)
POS <- ukb$map$physical.pos
dim(G) # 362,320 x 1,117,493
file.size(G$backingfile) / 1024^3  # 377 GB



# Write bed file
ukb$map <- dplyr::mutate(ukb$map, chromosome = as.integer(chromosome),
                         marked.ID = rsid, genetic.dist = 0)
ukb$fam <- snp_fake(n = nrow(G), m = 1)$fam
ukb$fam$sample.ID <- sample$ID_2[ind.indiv[sub2]]



RNGversion("3.5.1"); set.seed(1)
ind.val <- sort(sample(nrow(G), nrow(G)/5))
# for simulations
ind.gwas <- sort(sample(setdiff(rows_along(G), ind.val), (nrow(G)/5)*4))
ind.test <- setdiff(rows_along(G), c(ind.gwas, ind.val))
save(ind.val, ind.gwas, ind.test, file = paste0(refLoc,"ind_gwas_val_test.RData"))
# for real applications
ind.test <- setdiff(rows_along(G), ind.val)
save(ind.val, ind.test, file = paste0(refLoc,"ind_val_test.RData"))

snp_writeBed(ukb, bedfile = paste0(refLoc,"UKBB_imp_HM3_val.bed"), ind.row = ind.val)


###############################
# https://github.com/privefl/paper-ldpred2/blob/master/code/provide-ld-ref.R
library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)
#install.packages("runonce")
library(runonce)

#### Filter data from allele frequency errors ####

ukb <- snp_attach(paste0(refLoc,"UKBB_imp_HM3.rds"))
G <- ukb$genotypes
map <- transmute(ukb$map,
                 chr = as.integer(chromosome), pos = physical.pos,
                 a0 = allele1, a1 = allele2, rsid)  # reversed somehow..



NCORES <- 32
af <- runonce::save_run(
  big_colstats(G, ncores = NCORES)$sum / (2 * nrow(G)),
  file = paste0(refLoc,"af.rds")
)

#######################
# download_1000G gives me a timeout error after 60 secs so I needed to manually download this from bash:
#tempDataLoca="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/customLDRef/data/scractch/tmp-data"
#cd $tempDataLoca
#wget https://ndownloader.figshare.com/files/17838962
#mv 17838962 1000G_phase3_common_norel.zip
#######################


obj.bed <- bed(download_1000G(tempDataLoca, overwrite = F, delete_zip = F))



map2 <- dplyr::transmute(obj.bed$map, chr = chromosome, pos = physical.pos,
                         a1 = allele1, a0 = allele2, beta = 1)
info_snp <- snp_match(map2, map)
# 1,664,852 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,058,123 variants have been matched; 0 were flipped and 0 were reversed.


fam2 <- bigreadr::fread2(sub_bed(obj.bed$bedfile, ".fam2"))
ind_eur <- which(fam2$`Super Population` == "AFR")
af2 <- bed_MAF(obj.bed, ind.row = ind_eur, ind.col = info_snp$`_NUM_ID_.ss`,
               ncores = NCORES)
numInPop = length(ind_eur)


af_UKBB <- af[info_snp$`_NUM_ID_`]
af_1000G <- af2$af

N <- nrow(G)
af_avg <- (af_UKBB * N + af_1000G * numInPop) / (N + numInPop)
chi2 <- (af_UKBB - af_1000G)^2 / (af_avg * (1 - af_avg) * (1 / N + 1 / numInPop))
hist(pval <- pchisq(chi2, df = 1, lower.tail = FALSE), breaks = 0:50 / 50)

is_bad <- pval < 1e-5 | af_1000G < 0.01 | af_1000G > 0.99 | af_UKBB < 0.005 | af_UKBB > 0.995
qplot(af_UKBB, af_1000G, color = is_bad) +
  theme_bigstatsr() +
  geom_abline(color = "red") +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  labs(x = "Allele frequency in UKBB", y = "Allele frequency in 1000G",
       color = "Removed?")
	  
ggsave(paste0(refLoc,"af-outliers.png"), width = 9, height = 7)




ind_keep <- info_snp$`_NUM_ID_`[which(!is_bad)]
map_keep <- map[ind_keep, ]
map_keep$af_UKBB <- af_UKBB[which(!is_bad)]
vctrs::vec_duplicate_any(map_keep[, 1:2])  # FALSE

bigassertr::assert_dir(paste0(refLoc,"ld-ref")   )



#### Compute LD matrices ####

CHR <- map_keep$chr
POS2 <- snp_asGeneticPos(CHR, map_keep$pos, dir = tempDataLoca, ncores = NCORES)


# install.packages("future.batchtools")
# install.packages("parallelly")
# install.packages("rstudioapi")


# library(parallelly)
# library(future.batchtools)
# NCORES <- 15
# plan(batchtools_slurm(resources = list(
  # t = "12:00:00", c = NCORES, mem = "125g",
  # name = basename(rstudioapi::getSourceEditorContext()$path))))

lapply(1:22, function(chr) {

  ind.chr <- which(CHR == chr)
print(chr)
  corr <- snp_cor(
    G, ind.col = ind_keep[ind.chr], infos.pos = POS2[ind.chr],
    size = 3 / 1000, ncores = NCORES
  )

  saveRDS(corr, file = paste0(paste0(refLoc,"ld-ref/LD_chr"), chr, ".rds"), version = 2)
})


sum(sapply(list.files(paste0(refLoc,"ld-ref"), full.names = TRUE), file.size)) / 1024^3
# 8.5 GB


# Compute LD scores
map_keep$ld <- do.call('c', lapply(1:22, function(chr) {
  cat(chr, ".. ", sep = "")
  corr_chr <- readRDS(paste0(refLoc,"ld-ref/LD_chr", chr, ".rds"))
  Matrix::colSums(corr_chr^2)
}))

# add positions in different builds
options(timeout = 300)
liftOver <- runonce::download_file(
  "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver", tempDataLoca)
bigsnpr:::make_executable(liftOver)
map_keep$pos_hg17 <- snp_modifyBuild(map_keep, liftOver, from = "hg19", to = "hg17")$pos
map_keep$pos_hg18 <- snp_modifyBuild(map_keep, liftOver, from = "hg19", to = "hg18")$pos
map_keep$pos_hg38 <- snp_modifyBuild(map_keep, liftOver, from = "hg19", to = "hg38")$pos

saveRDS(map_keep, paste0(refLoc,"ld-ref/map.rds"), version = 2)
quit()

########################################
# get the genotypes in PLINK1 format for the above on the HM3 panel AFR indis

awk '{print $1"\t"$1}' $revisionsScratch$'UKBB_AFR_indis_TeamID17670' > $revisionsScratch$'UKBB_AFR_indis_keep'



awk '{print $2}' $hapmap3_b37bim > $hapmap3_b37bim$'_SNPs'
head $hapmap3_b37bim$'_SNPs'
wc -l $hapmap3_b37bim$'_SNPs' # 1403851

plinkFileList=$revisionsRaw$'mergelist3.txt'
rm -f $plinkFileList

for ((i=1; i<=$numChroms; i++)); do  # $numChroms

# export UKBB into PLINK1 format, using PLINK2 (as Bgen 1.2 only works with PLINK2), missingness 2%, 10% call threshold, only SNPs kept with unique coordinates, must add MAF filter again, as the MAF stats that come with the UKBB were computed on the 500K total not the 365K EUR
arguments=' --bgen /lustre/scratch123/hgi/mdt1/projects/ukbiobank_genotypes/FullRelease/Imputed/EGAD00010001474/ukb_imp_chr'$i$'_v3.bgen ref-first --sample /lustre/scratch123/hgi/mdt1/projects/ukbb_team152/17670/ibs/data/raw/ukb17670_imp_chr1_v3_s487395.sample  --maf 0.001 --hard-call-threshold 0.1 --geno 0.1 --extract '$hapmap3_b37bim$'_SNPs --keep '$revisionsScratch$'UKBB_AFR_indis_keep --make-bed --out '$revisionsRaw$'qc_'$i$'_AFR_hm3 --allow-extra-chr --allow-no-sex --snps-only --bp-space 1'
$plink2 $arguments

# create merge list
echo $revisionsRaw$'qc_'$i$'_AFR_hm3' >> ${plinkFileList}

done # end of chrom loop

# merge 
arguments=' --memory '$plinkMem$' --merge-list '$plinkFileList$' --make-bed --out '$revisionsRaw$'all_AFR_hm3 --allow-extra-chr --allow-no-sex'
$plink $arguments
head $revisionsRaw$'all_AFR_hm3.bim'
wc -l $revisionsRaw$'all_AFR_hm3.bim' # 1,102,616



# the above uses the publicID, we need to  rewrite the .fam file to use the same IDs that we will have for the .phe files
phenoFile=$revisionsRaw$'all_AFR_hm3.fam'
#phenoFile=$asthmaCrossAncestryScratch$'qc_'$i$'_AFR_hm3.fam'

cp $phenoFile $phenoFile$'.old'

awk 'FNR == NR { file1[ $2 ] = $1;  next; }
{print file1[$1]"\t"file1[$1]"\t"$3"\t"$4"\t"$5"\t"$6
}' OFS='\t' FS=" " /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data_ukbb/ukbb_hq_keeplist.txt_oldids $phenoFile$'.old' > $phenoFile
head $phenoFile
wc -l $phenoFile


# 0) Create Valid/Test divide in the UKBB
#$crossAncestrySumStats$'all_hm3.fam' 
#$crossAncestrySumStats$'all_height_hm3.fam' # these have the same number of indis in the same order

numItems=$(wc -l < "$revisionsRaw"all_AFR_hm3.fam)
mySeed=42
arguments='/nfs/users/nfs_m/mk23/scripts/randomNums.R '$numItems$' '$mySeed$' '$revisionsRaw$'randomNums'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# create a a TEST set indis, the first 50% of them
awk -v numItems="$numItems" '{if (FNR < numItems/2) print $0}' $revisionsRaw$'randomNums' > $revisionsRaw$'randomNums_TEST'
head $revisionsRaw$'randomNums_TEST'
wc -l $revisionsRaw$'randomNums_TEST'



awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $2 } }
' $revisionsRaw$'randomNums_TEST' $revisionsRaw$'all_AFR_hm3.fam' > $revisionsRaw$'all_AFR_hm3_TEST' 
head $revisionsRaw$'all_AFR_hm3_TEST' 
wc -l $revisionsRaw$'all_AFR_hm3_TEST' # 3205


######################
# 4) QC: remove amibguous alleles and low quality SNPs

#########
# Height:
remove_ambiguous_alleles $revisionsScratch$'height_AFR_hm3'
wc -l $revisionsScratch$'height_AFR_hm3_noAmbiguousAlleles' # 1255061

remove_ambiguous_alleles $revisionsScratch$'height_EUR_hm3'
wc -l $revisionsScratch$'height_EUR_hm3_noAmbiguousAlleles' # 1035914


arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_QC.R '$AFRRef$'/ '$revisionsScratch$'height_AFR_hm3_noAmbiguousAlleles '$revisionsScratch$'height_AFR_hm3_keep 0'
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments
#out of 996836 QC excluded 37967, and kept 958869 SNPs 


arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_QC.R '$baseLDpredLoc$'/ '$revisionsScratch$'height_EUR_hm3_noAmbiguousAlleles '$revisionsScratch$'height_EUR_hm3_keep 0'
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments
# out of 963563 QC excluded 541, and kept 963022 SNPs

awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $revisionsScratch$'height_AFR_hm3_keep' $revisionsScratch$'height_EUR_hm3_keep' > $revisionsScratch$'height_keep'


# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $revisionsScratch$'height_keep' $revisionsScratch$'height_AFR_hm3_noAmbiguousAlleles'  > $revisionsRaw$'height_AFR_HM3'
head $revisionsRaw$'height_AFR_HM3'
wc -l $revisionsRaw$'height_AFR_HM3' # 875599


# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $revisionsScratch$'height_keep' $revisionsScratch$'height_EUR_hm3_noAmbiguousAlleles'  > $revisionsRaw$'height_EUR_HM3'
head $revisionsRaw$'height_EUR_HM3'
wc -l $revisionsRaw$'height_EUR_HM3' # 875599



##############
# BMI:
remove_ambiguous_alleles $revisionsScratch$'bmi_AFR_hm3'
wc -l $revisionsScratch$'bmi_AFR_hm3_noAmbiguousAlleles' # 1255060

remove_ambiguous_alleles $revisionsScratch$'bmi_EUR_hm3'
wc -l $revisionsScratch$'bmi_EUR_hm3_noAmbiguousAlleles' # 1043906


arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_QC.R '$AFRRef$'/ '$revisionsScratch$'bmi_AFR_hm3_noAmbiguousAlleles '$revisionsScratch$'bmi_AFR_hm3_keep 0'
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments
#  out of 996879 QC excluded 96272, and kept 900607 SNPs

arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_QC.R '$baseLDpredLoc$'/ '$revisionsScratch$'bmi_EUR_hm3_noAmbiguousAlleles '$revisionsScratch$'bmi_EUR_hm3_keep 0'
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments
# "out of 966019 QC excluded 550, and kept 965469 SNPs

awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $revisionsScratch$'bmi_AFR_hm3_keep' $revisionsScratch$'bmi_EUR_hm3_keep' > $revisionsScratch$'bmi_keep'


# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $revisionsScratch$'bmi_keep' $revisionsScratch$'bmi_AFR_hm3_noAmbiguousAlleles'  > $revisionsRaw$'bmi_AFR_HM3'
head $revisionsRaw$'bmi_AFR_HM3'
wc -l $revisionsRaw$'bmi_AFR_HM3' # 824641

# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $revisionsScratch$'bmi_keep' $revisionsScratch$'bmi_EUR_hm3_noAmbiguousAlleles'  > $revisionsRaw$'bmi_EUR_HM3'
head $revisionsRaw$'bmi_EUR_HM3'
wc -l $revisionsRaw$'bmi_EUR_HM3' # 824643




##############
# LDL:
remove_ambiguous_alleles $revisionsScratch$'LDL_AFR_hm3'
wc -l $revisionsScratch$'LDL_AFR_hm3_noAmbiguousAlleles' # 1255687

remove_ambiguous_alleles $revisionsScratch$'LDL_EUR_hm3'
wc -l $revisionsScratch$'LDL_EUR_hm3_noAmbiguousAlleles' #989592


arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_QC.R '$AFRRef$'/ '$revisionsScratch$'LDL_AFR_hm3_noAmbiguousAlleles '$revisionsScratch$'LDL_AFR_hm3_keep 0'
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments
# out of 997444 QC excluded 31848, and kept 965596 SNPs,
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_QC.R '$baseLDpredLoc$'/ '$revisionsScratch$'LDL_EUR_hm3_noAmbiguousAlleles '$revisionsScratch$'LDL_EUR_hm3_keep 0'
/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments
# out of 927822 QC excluded 370, and kept 927452 SNPs


awk 'NR==FNR { lines[$0]=1; next } $0 in lines' $revisionsScratch$'LDL_AFR_hm3_keep' $revisionsScratch$'LDL_EUR_hm3_keep' > $revisionsScratch$'LDL_keep'


# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $revisionsScratch$'LDL_keep' $revisionsScratch$'LDL_AFR_hm3_noAmbiguousAlleles'  > $revisionsRaw$'LDL_AFR_HM3'
head $revisionsRaw$'LDL_AFR_HM3'
wc -l $revisionsRaw$'LDL_AFR_HM3' # 849698


# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $revisionsScratch$'LDL_keep' $revisionsScratch$'LDL_EUR_hm3_noAmbiguousAlleles'  > $revisionsRaw$'LDL_EUR_HM3'
head $revisionsRaw$'LDL_EUR_HM3'
wc -l $revisionsRaw$'LDL_EUR_HM3' # 849698



##########################

# ShaPRS
$revisionsRaw$'bmi_EUR_HM3'
$revisionsRaw$'bmi_AFR_HM3'

$revisionsRaw$'LDL_AFR_HM3'
$revisionsRaw$'LDL_EUR_HM3'

$revisionsRaw$'height_EUR_HM3'
$revisionsRaw$'height_AFR_HM3'

#phenos
$revisionsRaw$'UKBB_BMI_AFR.phe'
$revisionsRaw$'UKBB_height_AFR.phe'
$revisionsRaw$'UKBB_LDL_AFR.phe'


# test set geno
$revisionsRaw$'all_AFR_hm3'

# LD ref panels
$baseLDpredLoc # EUR
$AFRRef # AFR


revisionsResults=$revisions$'results/'
mkdir -p $revisionsResults


pheno='AFR_EUR_height'
shaPRS_new $revisionsRaw$'height_AFR_HM3' $revisionsRaw$'height_EUR_HM3' $revisionsResults$pheno '0' $AFRRef $baseLDpredLoc

arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_qval_manhattan.R '$revisionsResults$pheno$'_SE_meta '$revisionsResults$pheno$'_lFDR_meta_SNP_lFDR '$revisionsResults$pheno$' '$pheno$' '$revisionsResults$pheno$'_sumstats_meta'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


pheno='AFR_EUR_BMI'
shaPRS_new $revisionsRaw$'bmi_AFR_HM3' $revisionsRaw$'bmi_EUR_HM3' $revisionsResults$pheno '0' $AFRRef $baseLDpredLoc

arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_qval_manhattan.R '$revisionsResults$pheno$'_SE_meta '$revisionsResults$pheno$'_lFDR_meta_SNP_lFDR '$revisionsResults$pheno$' '$pheno$' '$revisionsResults$pheno$'_sumstats_meta'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


pheno='AFR_EUR_LDL'
shaPRS_new $revisionsRaw$'LDL_AFR_HM3' $revisionsRaw$'LDL_EUR_HM3' $revisionsResults$pheno '0' $AFRRef $baseLDpredLoc

arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_qval_manhattan.R '$revisionsResults$pheno$'_SE_meta '$revisionsResults$pheno$'_lFDR_meta_SNP_lFDR '$revisionsResults$pheno$' '$pheno$' '$revisionsResults$pheno$'_sumstats_meta'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# to evaluate PRS-CSx fairly, we need to provide the same AFR reference panel that we used for shaPRS

# convert LDpredformat to PRS-CS format:
cd $prscsrefs

wget https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=1

wget https://www.dropbox.com/s/dtccsidwlb6pbtv/ldblk_ukbb_afr.tar.gz?dl=1
# decompress them
mv ldblk_ukbb_afr.tar.gz?dl=1 ldblk_ukbb_afr.tar.gz
tar -xvf ldblk_ukbb_afr.tar.gz

pheno='PRSCS_AFR'
largePlinkMem=200000
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_afr/ '$AFRRef$'/ '$prscsrefs$'ldblk_ukbb_customAFR/'
bsub -q yesterday -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

# as PRS-CSx uses the hardcoded population labels, and it expects the LDref corresponding to the "AFR" pop to be named 'ldblk_ukbb_afr'
# we have to rename the original, and substitute it with the new ones
mv $prscsrefs$'ldblk_ukbb_afr/' $prscsrefs$'ldblk_ukbb_afr_orig/'

mv $prscsrefs$'ldblk_ukbb_customAFR/' $prscsrefs$'ldblk_ukbb_afr/' 




####################################




##############
# I) Convert sumstats into PRSCS format + build the initial EUR/AFR PRS

# ShaPRS
$revisionsRaw$'bmi_EUR_HM3'
$revisionsRaw$'bmi_AFR_HM3'

$revisionsRaw$'LDL_AFR_HM3'
$revisionsRaw$'LDL_EUR_HM3'

$revisionsRaw$'height_EUR_HM3'
$revisionsRaw$'height_AFR_HM3'


$revisionsRaw$'AFR_EUR_BMI'
$revisionsRaw$'AFR_EUR_LDL'
$revisionsRaw$'AFR_EUR_height'




#phenos
$revisionsRaw$'UKBB_BMI_AFR.phe'
$revisionsRaw$'UKBB_height_AFR.phe'
$revisionsRaw$'UKBB_LDL_AFR.phe'


# test set geno
$revisionsRaw$'all_AFR_hm3'

# LD ref panels
$baseLDpredLoc # EUR
$AFRRef # AFR

$prscsrefs$'ldblk_ukbb_afr/'




#PRS_array=( 'AFR_EUR_BMI' 'AFR_EUR_height' )
PRS_array=( 'AFR_EUR_BMI' 'AFR_EUR_height' 'AFR_EUR_LDL' ) 
EUR_array=( 'bmi_EUR_HM3' 'height_EUR_HM3' 'LDL_EUR_HM3' )
AFR_array=( 'bmi_AFR_HM3' 'height_AFR_HM3' 'LDL_AFR_HM3' )


arraylength=${#PRS_array[@]}

j=3
i=21

for (( j=1; j<${arraylength}+1; j++ )); do
pheno=${PRS_array[$j-1]} 
binPheno='0'
EURPheno=${EUR_array[$j-1]} 
AFRPheno=${AFR_array[$j-1]} 


# 2) Conver them to PRS-CS format
convertToPRSCS $revisionsRaw$EURPheno $binPheno
convertToPRSCS $revisionsRaw$AFRPheno $binPheno



# 3) run PRS-CSx
Neur=$(tail -n 1 $revisionsRaw$EURPheno | awk '{ print sprintf("%0.f",$10)}') # apply rounding to get integer
Nafr=$(tail -n 1 $revisionsRaw$AFRPheno | awk '{ print sprintf("%0.f",$10)}') # apply rounding to get integer


validBim=$revisionsRaw$AFRPheno$'_validBim'
# create a 'fake' .bim file for PRSx from my sumstats file that has all the info
awk '{if (FNR !=1) {print $1"\t"$3"\t0\t"$2"\t"$4"\t"$5} }' $revisionsRaw$AFRPheno > $validBim$'.bim'


mkdir -p $revisionsResults$'logs/'
for ((i=1; i<=22; i++)); do # $numChroms


# 1) PRS-CSx

# check if output doesn't already exist
outfile=$revisionsResults$pheno$'_AFR_pst_eff_a1_b0.5_phiauto_chr'$i$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
b=1
else
echo "SUBMITTING: "$pheno$' '$i
bsub -q normal -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J ${pheno}_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $revisionsResults$'logs/'$pheno$'_'$i$'.out' -e $revisionsResults$'logs/'$pheno$'_'$i'.err'  "performPRSCSx_AFR $revisionsRaw$EURPheno$'_PRSCS' $revisionsRaw$AFRPheno$'_PRSCS' $Neur $Nafr $revisionsRaw$EURPheno$'_validBim' $pheno $i $NCORES_LDPred2 $revisionsResults"

fi

done # end of chroms

done



###############################################################



# II) Validate & Evaluate the PRS-CSx baselines

#PRS_array=( 'AFR_EUR_BMI' 'AFR_EUR_height' )
PRS_array=( 'AFR_EUR_BMI' 'AFR_EUR_height' 'AFR_EUR_LDL' ) 
EUR_array=( 'bmi_EUR_HM3' 'height_EUR_HM3' 'LDL_EUR_HM3' )
AFR_array=( 'bmi_AFR_HM3' 'height_AFR_HM3' 'LDL_AFR_HM3' )

phenoFiles=( $revisionsRaw$'UKBB_BMI_AFR.phe' $revisionsRaw$'UKBB_height_AFR.phe' $revisionsRaw$'UKBB_LDL_AFR.phe' )


arraylength=${#PRS_array[@]}

j=3
for (( j=1; j<${arraylength}+1; j++ )); do
pheno=${PRS_array[$j-1]} 
binPheno='0'
EURPheno=${EUR_array[$j-1]} 
AFRPheno=${AFR_array[$j-1]} 

# 1) concat all chrom PRS into single files
cat $revisionsResults$pheno$'_EUR_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $revisionsResults$pheno$'_EUR_PRSCSx'
cat $revisionsResults$pheno$'_AFR_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $revisionsResults$pheno$'_AFR_PRSCSx'

# remove dupes:
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $revisionsResults$pheno$'_EUR_PRSCSx' > $revisionsResults$pheno$'_EUR_PRSCSx_no_dupes'
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $revisionsResults$pheno$'_AFR_PRSCSx' > $revisionsResults$pheno$'_AFR_PRSCSx_no_dupes'




phenoFile=${phenoFiles[$j-1]} 
# 2) build PRS profile for both EUR/AFR on the full UKBB

# PRSCSx produces a PRS file of format:
#  1       2            3       4   5         6
# 22	rs9605903	17054720	C	T	4.954646e-04
# PLINK wants to know the variant ID, A1 and coef: these are 2 4 6

#EUR 
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$revisionsRaw$'all_AFR_hm3 --pheno '$phenoFile$' --score '$revisionsResults$pheno$'_EUR_PRSCSx_no_dupes 2 4 6 sum --out '$revisionsResults$pheno$'_EUR_PRSCSx.PRS'
$plink $arguments

#AFR 
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$revisionsRaw$'all_AFR_hm3 --pheno '$phenoFile$' --score '$revisionsResults$pheno$'_AFR_PRSCSx_no_dupes 2 4 6 sum --out '$revisionsResults$pheno$'_AFR_PRSCSx.PRS'
$plink $arguments



# 3) Perform Stage 2: find the best linear combination of the AFR/EUR + evaluate
arguments='/nfs/users/nfs_m/mk23/scripts/PRSLinearComb2.R '$revisionsResults$pheno$'_AFR_PRSCSx.PRS.profile '$revisionsResults$pheno$'_EUR_PRSCSx.PRS.profile '$phenoFile$' '$revisionsRaw$'all_AFR_hm3_TEST '$binPheno$' '$revisionsResults$pheno$'_PRSCSx_result'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


done


# BMI:correlation_sq 0.021 (sd: 0.000438 )

# Height: correlation_sq 0.0282 (sd: 0.000484 )

# LDL: correlation_sq 0.0563 (sd: 0.000461 ) 


###############################################################

# because PRS-CSx does its own QC routine on SNPs, we loose a few K more SNPs if we start from the LDpred2 panel
# to ensure we work off the exact same subset, we need to subset the LDpred2 SNPs with the PRS-CSx, and re-run shaPRS:


#PRS_array=( 'AFR_EUR_BMI' 'AFR_EUR_height' )
PRS_array=( 'AFR_EUR_BMI' 'AFR_EUR_height' 'AFR_EUR_LDL' )
EUR_array=( 'bmi_EUR_HM3' 'height_EUR_HM3' 'LDL_EUR_HM3' )
AFR_array=( 'bmi_AFR_HM3' 'height_AFR_HM3' 'LDL_AFR_HM3' )
phenoFiles=( $revisionsRaw$'UKBB_BMI_AFR.phe' $revisionsRaw$'UKBB_height_AFR.phe' $revisionsRaw$'UKBB_LDL_AFR.phe' )

arraylength=${#PRS_array[@]}

j=3
i=2
plinkMem=100000
for (( j=1; j<${arraylength}+1; j++ )); do
pheno=${PRS_array[$j-1]} 
binPheno='0'
EURPheno=${EUR_array[$j-1]} 
AFRPheno=${AFR_array[$j-1]} 

Neur=$(tail -n 1 $revisionsRaw$EURPheno | awk '{ print sprintf("%0.f",$10)}') # apply rounding to get integer
Nafr=$(tail -n 1 $revisionsRaw$AFRPheno | awk '{ print sprintf("%0.f",$10)}') # apply rounding to get integer
validBim=$revisionsRaw$AFRPheno$'_validBim'

# the final PRS for PRS-CSx is a composite of 2 PRS, which had a different set of SNPs
# we are comparing ourselves most directly to the AFR PRS-CSx, so use that as a subset
# wc -l $revisionsResults$pheno$'_AFR_PRSCSx_no_dupes'

awk 'FNR == NR { sums[ $2 ] = $2; next; } FNR <= NR { if( FNR == 1 || $3 in sums) {print $0 } }
' $revisionsResults$pheno$'_AFR_PRSCSx_no_dupes' $revisionsRaw$AFRPheno > $revisionsRaw$AFRPheno$'_QC_CSxsubset'

awk 'FNR == NR { sums[ $2 ] = $2; next; } FNR <= NR { if( FNR == 1 || $3 in sums) {print $0 } }
' $revisionsResults$pheno$'_AFR_PRSCSx_no_dupes' $revisionsResults$pheno'_sumstats_meta' > $revisionsResults$pheno'_sumstats_meta_QC_CSxsubset'


# 2) Convert them to PRS-CS format 
convertToPRSCS $revisionsRaw$AFRPheno$'_QC_CSxsubset' $binPheno
convertToPRSCS $revisionsResults$pheno'_sumstats_meta_QC_CSxsubset' $binPheno



for ((i=1; i<=22; i++)); do #

# LDpred2: baseline primary via LDpred2
outfile=$revisionsResults$'AFR_AFR_LD_'$pheno$'_CSxsubset_'$i
#mv $outfile $outfile$'_empirical' # backup originals
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "LDpred2 ALREADY EXISTS: "$outfile
b=1
else
#  AFR  (AFR LD panel)
echo "LDpred SUBMITTING: "$pheno$'_AFR crhom'$i
avgNumIndis=$(tail -n 1 $revisionsRaw$AFRPheno | awk '{ print sprintf("%0.f",$10)}') # apply rounding to get integer
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_chrom.R '$AFRRef$' '$revisionsRaw$AFRPheno$'_QC_CSxsubset '$i$' '$avgNumIndis$' '$revisionsResults$' AFR_AFR_LD_'$pheno$'_CSxsubset_'$i$' '$revisionsRaw$'all_AFR_hm3 '$shaPRSscriptLoc$' '$binPheno

#rm -rf $revisionsResults$'AFR_AFR_LD_'$pheno$'_'$i$'.out'
#rm -rf $revisionsResults$'AFR_AFR_LD_'$pheno$'_'$i$'.err'
#bsub -q long -m "modern_hardware" -G team152 -n4 -J AFR_AFR_LD_${pheno}_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $revisionsResults$'AFR_AFR_LD_'$pheno$'_'$i$'.out' -e $revisionsResults$'AFR_AFR_LD_'$pheno$'_'$i$'.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

#/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments


fi

# ShaPRS+ LDpred2: shaPRS using the custom LD panel via LDpred2
outfile=$revisionsResults$'shaPRS_LD_'$pheno$'_CSxsubset_'$i
#mv $outfile $outfile$'_empirical' # backup originals
#rm -rf $outfile
if [ -s "$outfile" ] ; then
b=1
#echo "ShaPRS+ LDpred2 ALREADY EXISTS: "$outfile
else
# AFR-shaPRS shaPRS LD Panel
echo "ShaPRS+ LDpred2 SUBMITTING: "$pheno$'_AFR crhom'$i
avgNumIndis=$(tail -n 1 "$revisionsResults$pheno"_sumstats_meta | awk '{ print sprintf("%0.f",$10)}') # apply rounding to get integer
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_chrom.R '$revisionsResults$pheno$'/ '$revisionsResults$pheno'_sumstats_meta_QC_CSxsubset '$i$' '$avgNumIndis$' '$revisionsResults$' shaPRS_LD_'$pheno$'_CSxsubset_'$i$' '$revisionsRaw$'all_AFR_hm3 '$shaPRSscriptLoc$' '$binPheno

#rm -rf $revisionsResults$'shaPRS_LD_'$pheno$'_'$i$'.out'
#rm -rf $revisionsResults$'shaPRS_LD_'$pheno$'_'$i$'.err'
#bsub -q long -m "modern_hardware" -G team152 -n4 -J shaPRS_LD_${pheno}_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $revisionsResults$'shaPRS_LD_'$pheno$'_'$i$'.out' -e $revisionsResults$'shaPRS_LD_'$pheno$'_'$i$'.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"
#/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments
fi



# PRS-CS: baseline PRS-CS, with the built in AFR LDpanel
# check if output doesn't already exist
outfile=$revisionsResults$pheno$'_AFR_AFR_pst_eff_a1_b0.5_phiauto_chr'$i$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "PRS-CS ALREADY EXISTS: "$outfile
b=1
else
b=1
echo "PRS-CS SUBMITTING: "$pheno$'_AFR '$i
#rm -rf $revisionsResults$'logs/'$pheno$'_AFR_'$i$'.out'
#rm -rf $revisionsResults$'logs/'$pheno$'_AFR_'$i'.err'
#bsub -q normal -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J ${pheno}_AFR_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $revisionsResults$'logs/'$pheno$'_AFR_'$i$'.out' -e $revisionsResults$'logs/'$pheno$'_AFR_'$i'.err'  "performPRSCS_AFR $revisionsRaw$AFRPheno$'_QC_CSxsubset_PRSCS' $Nafr $revisionsRaw$AFRPheno$'_validBim' $pheno$'_AFR' $i $NCORES_LDPred2 $revisionsResults"
#performPRSCS_AFR $revisionsRaw$AFRPheno$'_QC_CSxsubset_PRSCS' $Nafr $revisionsRaw$AFRPheno$'_validBim' $pheno$'_AFR' $i $NCORES_LDPred2 $revisionsResults

fi


# ShaPRS-PRS-CS: PRS-CS  (IE not CSx) applied onto the post shaPRS data
# check if output doesn't already exist
outfile=$revisionsResults$pheno$shaPRS$'_AFR/_pst_eff_a1_b0.5_phiauto_chr'$i$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ShaPRS-PRS-CS ALREADY EXISTS: "$outfile
b=1
else
b=1
echo "ShaPRS-PRS-CS SUBMITTING: AFR "$pheno$shaPRS$' '$i
# custom LD panel
#rm -rf $revisionsResults$'logs/'$pheno$shaPRS$'_'$i$'.out' 
#rm -rf $revisionsResults$'logs/'$pheno$shaPRS$'_'$i'.err'
#bsub -q normal -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J ${pheno}${shaPRS}_${i} -R "span[hosts=1] select[mem>${plinkMem}] rusage[mem=${plinkMem}]" -M${plinkMem} -o $revisionsResults$'logs/'$pheno$shaPRS$'_'$i$'.out' -e $revisionsResults$'logs/'$pheno$shaPRS$'_'$i'.err'  "performPRSCS_custom_AFR $revisionsResults$pheno'_sumstats_meta_QC_CSxsubset_PRSCS' $Nafr $revisionsRaw$AFRPheno$'_validBim' $pheno$shaPRS$'_AFR' $i $NCORES_LDPred2 $prscsrefs'ldblk_ukbb_afr' $revisionsResults"
#performPRSCS_custom_AFR $revisionsResults$pheno'_sumstats_meta_QC_CSxsubset_PRSCS' $Nafr $revisionsRaw$AFRPheno$'_validBim' $pheno$shaPRS$'_AFR' $i $NCORES_LDPred2 $prscsrefs'ldblk_ukbb_afr' $revisionsResults

fi

done # end of chroms

done # end of traits


#######################################################


# EVALUATE SHAPRS
#PRS_array=( 'AFR_EUR_BMI' 'AFR_EUR_height' )
PRS_array=( 'AFR_EUR_BMI' 'AFR_EUR_height' 'AFR_EUR_LDL' )
EUR_array=( 'bmi_EUR_HM3' 'height_EUR_HM3' 'LDL_EUR_HM3' )
AFR_array=( 'bmi_AFR_HM3' 'height_AFR_HM3' 'LDL_AFR_HM3' )


phenoFiles=( $revisionsRaw$'UKBB_BMI_AFR.phe' $revisionsRaw$'UKBB_height_AFR.phe' $revisionsRaw$'UKBB_LDL_AFR.phe' )

arraylength=${#PRS_array[@]}




# PRS-CS: baseline PRS-CS, with the built in AFR LDpanel
# ShaPRS-PRS-CS: PRS-CS  (IE not CSx) applied onto the post shaPRS data

j=1
for (( j=1; j<${arraylength}+1; j++ )); do
pheno=${PRS_array[$j-1]} 
binPheno='0'

AFRPheno=${AFR_array[$j-1]} 


# 1) concat all chrom PRS into single files

# a) PRS-CS
cat $revisionsResults$pheno$'_AFR_AFR_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $revisionsResults$pheno$'_AFR_AFR_PRSCS'
cat $revisionsResults$pheno$shaPRS$'_AFR/_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $revisionsResults$pheno$shaPRS$'_AFR_PRSCS'

# remove dupes:
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $revisionsResults$pheno$'_AFR_AFR_PRSCS' > $revisionsResults$pheno$'_AFR_AFR_PRSCS_no_dupes'
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $revisionsResults$pheno$shaPRS$'_AFR_PRSCS' > $revisionsResults$pheno$shaPRS$'_AFR_PRSCS_no_dupes'
# b) LDpred2
cat $revisionsResults$'AFR_AFR_LD_'$pheno$'_CSxsubset_'[!a-z] $revisionsResults$'AFR_AFR_LD_'$pheno$'_CSxsubset_'[!a-z][!a-z] > $revisionsResults$pheno$'_AFR_AFR_LDpred2_CSxsubset'
cat $revisionsResults$'shaPRS_LD_'$pheno$'_CSxsubset_'[!a-z] $revisionsResults$'shaPRS_LD_'$pheno$'_CSxsubset_'[!a-z][!a-z] > $revisionsResults$pheno$'_shaPRS_LDpred2_CSxsubset'

awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $revisionsResults$pheno$'_AFR_AFR_LDpred2_CSxsubset' > $revisionsResults$pheno$'_AFR_AFR_LDpred2_CSxsubset_no_dupes'
awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $revisionsResults$pheno$'_shaPRS_LDpred2_CSxsubset' > $revisionsResults$pheno$'_shaPRS_LDpred2_CSxsubset_no_dupes'


phenoFile=${phenoFiles[$j-1]} 
# 2) build PRS profile for both EUR/AFR on the full UKBB

# a) PRS-CS
# PRSCSx produces a PRS file of format:
#  1       2            3       4   5         6
# 22	rs9605903	17054720	C	T	4.954646e-04
# PLINK wants to know the variant ID, A1 and coef: these are 2 4 6

#PRS-CS
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$revisionsRaw$'all_AFR_hm3 --pheno '$phenoFile$' --score '$revisionsResults$pheno$'_AFR_AFR_PRSCS_no_dupes 2 4 6 sum --out '$revisionsResults$pheno$'_AFR_AFR_PRSCS.PRS'
$plink $arguments

evaluatePRS_NOVALID_phe $phenoFile $pheno$'_AFR_AFR_PRSCS.PRS' $binPheno
# AFR_EUR_BMI
# --score: 638552 valid predictors loaded. correlation_sq 0.0037 (sd: 0.000222 )
# AFR_EUR_height
# --score: 680312 valid predictors loaded. correlation_sq 0.00387 (sd: 0.00024 )
# AFR_EUR_LDL
# --score: 660233 valid predictors loaded. correlation_sq 0.00816 (sd: 0.000223 )

#SHAPRS PRS-CS
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$revisionsRaw$'all_AFR_hm3 --pheno '$phenoFile$' --score '$revisionsResults$pheno$shaPRS$'_AFR_PRSCS_no_dupes 2 4 6 sum --out '$revisionsResults$pheno$shaPRS$'_AFR_PRSCS.PRS'
$plink $arguments

evaluatePRS_NOVALID_phe $phenoFile $pheno$shaPRS$'_AFR_PRSCS.PRS' $binPheno

# AFR_EUR_BMI
# --score: 629146 valid predictors loaded. correlation_sq 0.0249 (sd: 0.000214 )
# AFR_EUR_height
# --score: 680311 valid predictors loaded. correlation_sq 0.019 (sd: 0.000223 )
# AFR_EUR_LDL
# --score: 660011 valid predictors loaded. correlation_sq 0.0187 (sd: 0.000233 )

# LDpred2: baseline primary via LDpred2
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$revisionsRaw$'all_AFR_hm3 --pheno '$phenoFile$' --score '$revisionsResults$pheno$'_AFR_AFR_LDpred2_CSxsubset_no_dupes sum --out '$revisionsResults$pheno$'_AFR_AFR_LDpred2_CSxsubset.PRS'
$plink $arguments

evaluatePRS_NOVALID_phe $phenoFile $pheno$'_AFR_AFR_LDpred2_CSxsubset.PRS' $binPheno

# AFR_EUR_BMI
# --score: 638552 valid predictors loaded  correlation_sq 0.00254 (sd: 0.000216 )
# AFR_EUR_height
# --score: 680312 valid predictors loaded  correlation_sq 0.00349 (sd: 0.000213 )
# AFR_EUR_LDL
# --score: 660233 valid predictors loaded. correlation_sq 0.0559 (sd: 0.000222 )

# ShaPRS+ LDpred2: shaPRS using the custom LD panel via LDpred2
arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$revisionsRaw$'all_AFR_hm3 --pheno '$phenoFile$' --score '$revisionsResults$pheno$'_shaPRS_LDpred2_CSxsubset_no_dupes sum --out '$revisionsResults$pheno$'_shaPRS_LDpred2_CSxsubset.PRS'
$plink $arguments

evaluatePRS_NOVALID_phe $phenoFile $pheno$'_shaPRS_LDpred2_CSxsubset.PRS' $binPheno
# AFR_EUR_BMI
# --score: 629146 valid predictors loaded. correlation_sq 0.0196 (sd: 0.000243 ) 
# AFR_EUR_height
# --score: 680311 valid predictors loaded. correlation_sq 0.0126 (sd: 0.000205 )
# AFR_EUR_LDL
# --score: 660011 valid predictors loaded. correlation_sq 0.0208 (sd: 0.000214 )


# PRS-CSx-unweighted
evaluatePRS_NOVALID_phe $phenoFile $pheno$'_AFR_PRSCSx.PRS' $binPheno 
# AFR_EUR_BMI
#  correlation_sq 0.0053 (sd: 0.000227 )
# AFR_EUR_height
#  correlation_sq 0.00911 (sd: 0.000221 )
# AFR_EUR_LDL
# correlation_sq 0.0296 (sd: 0.000219 )

# evaluatePRS_NOVALID_phe $phenoFile $pheno$'_EUR_PRSCSx.PRS' $binPheno 

done

# PRS-CSx:
# BMI:correlation_sq 0.021 (sd: 0.000438 )
# Height: correlation_sq 0.0282 (sd: 0.000484 )
# LDL: correlation_sq 0.0563 (sd: 0.000461 )  


# RESULTS:
             PRS-CSx      PRS-CSx-stage1       shaPRS+PRS-CS			shaPRS+LDPred2
height       0.0282       0.00911              0.00387                 0.0126
BMI          0.021        0.0053               0.0249                  0.0196
LDL          0.0563       0.0296               0.0187                  0.0208


             LDpred2     PRS-CS                PRS-CSx      PRS-CSx-stage1       shaPRS+PRS-CS			shaPRS+LDPred2
height       0.00349     0.00387               0.0282       0.00911              0.00387                 0.0126
BMI          0.00254     0.0037                0.021        0.0053               0.0249                  0.0196
LDL          0.0559      0.00816               0.0563       0.0296               0.0187                  0.0208





# LDPred2
# AFR_EUR_BMI
# --score: 638552 valid predictors loaded  correlation_sq 0.00254 (sd: 0.000216 )
# AFR_EUR_height
# --score: 680312 valid predictors loaded  correlation_sq 0.00349 (sd: 0.000213 )
# AFR_EUR_LDL
# --score: 660233 valid predictors loaded. correlation_sq 0.0559 (sd: 0.000222 )


# PRS-CS
# AFR_EUR_BMI
# --score: 638552 valid predictors loaded. correlation_sq 0.0037 (sd: 0.000222 )
# AFR_EUR_height
# --score: 680312 valid predictors loaded. correlation_sq 0.00387 (sd: 0.00024 )
# AFR_EUR_LDL
# --score: 660233 valid predictors loaded. correlation_sq 0.00816 (sd: 0.000223 )

########################################################
# FUNCTIONS:


# 3) exclude ambiguous alleles:  A/T or C/G
function remove_ambiguous_alleles { 
sumstatsLoc=$1
#  1     2     3    4     5
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N
awk ' 
FNR <= NR {
if (FNR == 1) {print $0 }
else { 
if( toupper($4) == "A" && toupper($5) == "T" || toupper($5) == "A" && toupper($4) == "T" || toupper($4) == "G" && toupper($5) == "C" || toupper($5) == "G" && toupper($4) == "C") {}
else {print $0}
}
}
' $sumstatsLoc > $sumstatsLoc$'_noAmbiguousAlleles'

head $sumstatsLoc$'_noAmbiguousAlleles'
wc -l $sumstatsLoc$'_noAmbiguousAlleles'
wc -l $sumstatsLoc
}
export -f remove_ambiguous_alleles # this makes local functions executable when bsubbed





# evaluates ready made profile scores
function evaluatePRS_NOVALID_phe {
phenoFile=$1
PRSName=$2
calcAUC3=$3
# call PLINK's score to build PRS for each



# find out accuracy 
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$3","file1[$2]  } } }' OFS=',' FS=" " $revisionsResults$PRSName$'.profile' FS="\t" $phenoFile  > $outste$revisionsResults$PRSName$'.csv'


arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$outste$revisionsResults$PRSName$'.csv '$revisionsResults$PRSName$'_resultNovalid '$calcAUC3
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments 
}





# Performs PRS-CSx for 2 pops on a given chrom
function performPRSCSx_AFR { 
eursums=$1
afrsums=$2
Neur=$3
Nafr=$4
validBim=$5
mypheno=$6
chrom=$7
ncores=$8
revisionsResults=$9

export MKL_NUM_THREADS=$ncores
export NUMEXPR_NUM_THREADS=$ncores
export OMP_NUM_THREADS=$ncores

prscsxdir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsscript/PRScsx/"
prscsrefs="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsrefs/"

python $prscsxdir/PRScsx.py --ref_dir=$prscsrefs --bim_prefix=$validBim --sst_file=$eursums,$afrsums --n_gwas=$Neur,$Nafr --pop=EUR,AFR --out_dir=$revisionsResults --out_name=$mypheno --chrom=$chrom  --seed=42

}
export -f performPRSCSx_AFR # this makes local functions executable when bsubbed



# Performs PRS-CS for a custom population
function performPRSCS_custom_AFR { 
sumsFile=$1
NGWAS=$2
validBim=$3
mypheno=$4
chrom=$5
ncores=$6
refLoc=$7
revisionsResults=$8

export MKL_NUM_THREADS=$ncores
export NUMEXPR_NUM_THREADS=$ncores
export OMP_NUM_THREADS=$ncores



prscssorig="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscss/"

mkdir -p $revisionsResults$mypheno
python $prscssorig/PRScs/PRScs.py --ref_dir=$refLoc --bim_prefix=$validBim --sst_file=$sumsFile --n_gwas=$NGWAS --out_dir=$revisionsResults$mypheno/ --chrom=$chrom  --seed=42


}
export -f performPRSCS_custom_AFR # this makes local functions executable when bsubbed





# Performs PRS-CS for 1 population (AFR)
function performPRSCS_AFR { 
sumsfile=$1
Nn=$2
validBim=$3
mypheno=$4
chrom=$5
ncores=$6
revisionsResults=$7

export MKL_NUM_THREADS=$ncores
export NUMEXPR_NUM_THREADS=$ncores
export OMP_NUM_THREADS=$ncores

prscsxdir="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsscript/PRScsx/"
prscsrefs="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsrefs/"

python $prscsxdir/PRScsx.py --ref_dir=$prscsrefs --bim_prefix=$validBim --sst_file=$sumsfile --n_gwas=$Nn --pop=AFR --out_dir=$revisionsResults --out_name=$mypheno --chrom=$chrom  --seed=42

}
export -f performPRSCS_AFR # this makes local functions executable when bsubbed




