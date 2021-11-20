####################################################
# dgx-server screens:

# screen -r -D  28652.pts-0.node-11-2-3 # IBD
# screen -r -D  57179.pts-1.node-11-1-1 # cross ancestry



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
Rscript='/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript'
flashpca='/nfs/users/nfs_m/mk23/software/flashpca/flashpca_x86-64'
qctool2_new='/nfs/team152/mk23/software/qctool2/qctool_v2.0.6-Ubuntu16.04-x86_64/qctool'
qctoolmem=10000

python2='python2'
ldsc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/0_software/ldsc/ldsc.py'
UKBBQCFile='/lustre/scratch115/realdata/mdt3/projects/ukbiobank/FullRelease/Imputed/001/ukb_sqc_v2.txt'
asthmaRawPhenos='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/data/raw/pheno_raw.txt'



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




####################################################
# ASTHMA-CROSS ANCESTRY
##########################

# 376510 high quality people
ukbbQCkeeplist='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/0_common/data/ukbb_hq_keeplist.txt' # Headerless file, col1: TeamID (for imputed data) col2: PublicID (for the genotypedata): 1587700 5172879 // 376510 indis
ukbbhapmap3loc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/pgslite2/height/data/GWAS_ALL' # this does NOT have the HLA
hapmap3_b37bim='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/0_common/data/twas/raw/hapmap3_r1_b37_fwd_consensus.qc.poly.recode.bim' # this is the full hapmap3, that has the HLA
asthaBaseLoc='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/'

asthmaCrossAncestry=$asthaBaseLoc$'crossAncestry/'
crossAncestrySumStats=$asthmaCrossAncestry$'sumstats/'
crossAncestryResults=$asthmaCrossAncestry$'results/'
crossAncestryRaw=$asthmaCrossAncestry$'raw/'
ref1KG=$crossAncestryRaw$'/ref1KG/'
shaPRSscriptLoc="/nfs/users/nfs_m/mk23/scripts/shaPRS.R"


asthmaCrossAncestryScratch=$asthaBaseLoc$'scratch/'

eursumstatsLoc=$crossAncestrySumStats$'eur/'
japsumstatsLoc=$crossAncestrySumStats$'jap/'

# Find the 1000 Genomes dataset
thousandgenomesLoc="/lustre/scratch118/core/corebio/g1k/vcf/"
vcf_part1="ALL.chr"
vcf_part2=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
ref1KG_PLINK='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/retreat2019/data/scratch/'
indiPopulations="/lustre/scratch118/core/corebio/g1k/vcf/integrated_call_samples_v3.20130502.ALL.panel"
popALDpanel='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/'
popBLDpanel='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/JPTRef/'

asthmaEURRef="asthmaeur"
asthmaEASRef="asthmaeas"

heightEURRef="heighteur"
heightEASRef="heighteas"

t2dEURRef="t2deur"
t2dEASRef="t2deas"

cadEURRef="cadeur"
cadEASRef="cadeas"

brcaEURRef="brcaeur"
brcaEASRef="brcaeas"




mkdir -p $crossAncestryResults

mkdir -p $asthmaCrossAncestryScratch

mkdir -p $crossAncestrySumStats
mkdir -p $ref1KG
cd $crossAncestrySumStats

mkdir -p $japsumstatsLoc
mkdir -p $eursumstatsLoc

cd $eursumstatsLoc

##################################
# T2D:
mkdir -p $crossAncestryRaw$'T2D/'
cd $crossAncestryRaw$'T2D/'


# EUR:
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004773/METAANALYSIS_DIAGRAM_SE1.txt

head METAANALYSIS_DIAGRAM_SE1.txt
wc -l METAANALYSIS_DIAGRAM_SE1.txt  #12,056,347

# calculate Effective sample size from cases/controls: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html
# cases= 26676
# controls=132532
# N_eff = 4 / (1 / cases + 1 / controls) 
# round(N_eff) # 88825


#   1              2       3       4       5        6          7
#Chr:Position    Allele1 Allele2 Effect  StdErr  P-value TotalSampleSize
#5:85928892      T       C       -0.013  0.026   0.61    158186
# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N


# convert to common format ( needed to get the rsids from the 1KG, B37 and Effect for the Allele 1 and corresponding Standard Error I checked in # DIAGRAM_1000G_GWAS.pdf)
awk ' 
FNR == NR { file1[ $1":"$4 ] = $1"\t"$4"\t"$2; next; } 
FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { if( $1 in file1 ) {print file1[$1]"\t"toupper($2)"\t"toupper($3)"\tX\t"$4"\t"$5"\t"$6"\t88825"} }
}
' OFS='\t' $ref1KG$'EUR_ALL.bim' METAANALYSIS_DIAGRAM_SE1.txt > $eursumstatsLoc$'EUR_T2D_raw'

head $eursumstatsLoc$'EUR_T2D_raw'
wc -l $eursumstatsLoc$'EUR_T2D_raw'


# get rid of ambigous alleles
remove_ambiguous_alleles $eursumstatsLoc$'EUR_T2D_raw'


# hapmap 3 subset
awk ' 
FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {
if (FNR == 1) {print $0 }
else { if( $3 in file1 ) {print $0} }
}
' OFS='\t' $hapmap3_b37bim $eursumstatsLoc$'EUR_T2D_raw_noAmbiguousAlleles' > $eursumstatsLoc$'EUR_T2D_hm3'

head $eursumstatsLoc$'EUR_T2D_hm3'
wc -l $eursumstatsLoc$'EUR_T2D_hm3' # 1161973



##############
# JP: 
wget http://jenger.riken.jp/14/
mv index.html td2.gz
gunzip td2.gz
head td2
wc -l td2  # 12557762

#cases= 36614
#controls=155150
#N_eff = 4 / (1 / cases + 1 / controls) 
#N_eff # 118493



#  1              2       3        4       5       6       7       8       9       10      11      12
#SNP             CHR     POS       REF     ALT     Frq     BETA    SE      P       Dir    HetP     N
#1:725932_G_A    1       725932    A       G       0.0040  -0.0737 0.1394  0.597  -?+-     0.4502  166718
#                                  REF and ALT are the wrong way around...
# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N

# convert to common format 
awk ' 
FNR == NR { file1[ $1":"$4 ] = $1"\t"$4"\t"$2; next; } 
FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { 
split($1,chr,":")  # split 1:725932_G_A  -> 1 725932_G_A 
split(chr[2],pos,"_") # split 725932_G_A   -> 725932, G, A 

if (chr[1]":"pos[1] in file1)
print file1[chr[1]":"pos[1]]"\t"toupper($4)"\t"toupper($5)"\tX\t"$7"\t"$8"\t"$9"\t118493" }
}
' OFS='\t' $ref1KG$'EUR_ALL.bim' td2 > $japsumstatsLoc$'JP_T2D_raw'
head $japsumstatsLoc$'JP_T2D_raw'


# get rid of ambigous alleles
remove_ambiguous_alleles $japsumstatsLoc$'JP_T2D_raw'

#cp $japsumstatsLoc$'JP_T2D_raw' $japsumstatsLoc$'JP_T2D_raw_old'
#cp $japsumstatsLoc$'JP_T2D_raw_noAmbiguousAlleles' $japsumstatsLoc$'JP_T2D_raw_noAmbiguousAlleles_old'
#cp $japsumstatsLoc$'JP_T2D_hm3' $japsumstatsLoc$'JP_T2D_hm3_old'

# hapmap 3 subset
awk ' 
FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {
if (FNR == 1) {print $0 }
else { if( $3 in file1 ) {print $0} }
}
' OFS='\t' $hapmap3_b37bim $japsumstatsLoc$'JP_T2D_raw_noAmbiguousAlleles' > $japsumstatsLoc$'JP_T2D_hm3'

head $japsumstatsLoc$'JP_T2D_hm3'
wc -l $japsumstatsLoc$'JP_T2D_hm3' # 1067912


#################################################################################
# HEIGHT:
mkdir -p $crossAncestryRaw$'height/'
cd $crossAncestryRaw$'height/'



# EUR: a
wget https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz
gunzip GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz

wc -l GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt # 2,550,859
head GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt
#MarkerName      Allele1 Allele2 Freq.Allele1.HapMapCEU  b       SE      p      N
#rs4747841       A       G       0.551   -0.0011 0.0029  0.70    253213
# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N


# convert to common format ( needed to get the chr/pos from the 1KGP project bims
awk ' 
FNR == NR { file1[ $2 ] = $1"\t"$4; next; } 
FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { if( $1 in file1 ) {print file1[$1]"\t"$0} }
}
' OFS='\t' $ref1KG$'EUR_ALL.bim' GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt > $eursumstatsLoc$'EUR_height_raw'

head $eursumstatsLoc$'EUR_height_raw'
wc -l $eursumstatsLoc$'EUR_height_raw'

# get rid of ambigous alleles
remove_ambiguous_alleles $eursumstatsLoc$'EUR_height_raw'

# hapmap 3 subset
awk ' 
FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {
if (FNR == 1) {print $0 }
else { if( $3 in file1 ) {print $0} }
}
' OFS='\t' $hapmap3_b37bim $eursumstatsLoc$'EUR_height_raw_noAmbiguousAlleles' > $eursumstatsLoc$'EUR_height_hm3'

head $eursumstatsLoc$'EUR_height_hm3'
wc -l $eursumstatsLoc$'EUR_height_hm3' # 1035914


#################################################################################
# JP:
wget http://jenger.riken.jp:8080/download/Height_autosomes
mv Height_autosomes Height_autosomes.gz
gunzip Height_autosomes.gz

head Height_autosomes
wc -l Height_autosomes # 5,657,428
#   1         2       3       4        5          6           7        8        9
# chrom      pos     ref     alt     rsids   nearest_genes   pval    beta    sebeta maf
#   1       707522  G       C       rs371890604     AL669831.1      0.47    -0.00550.0077   0.12
# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N


# convert to common format

avgNumIndis=159095  #
awk -v avgNumIndis="$avgNumIndis" ' 
FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { print $1"\t"$2"\t"$5"\t"$4"\t"$3"\tX\t"$8"\t"$9"\t"$7"\t"avgNumIndis }
}
' OFS='\t' Height_autosomes > $japsumstatsLoc$'JP_height_raw'
head $japsumstatsLoc$'JP_height_raw'


# get rid of ambigous alleles
remove_ambiguous_alleles $japsumstatsLoc$'JP_height_raw'


# hapmap 3 subset
awk ' 
FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {
if (FNR == 1) {print $0 }
else { if( $3 in file1 ) {print $0} }
}
' OFS='\t' $hapmap3_b37bim $japsumstatsLoc$'JP_height_raw_noAmbiguousAlleles' > $japsumstatsLoc$'JP_height_hm3'

head $japsumstatsLoc$'JP_height_hm3'
wc -l $japsumstatsLoc$'JP_height_hm3' # 912106


#################################a
# BRCA:
mkdir -p $crossAncestryRaw$'BRCA/'
cd $crossAncestryRaw$'BRCA/'

# EUR:
wget 'ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004988/oncoarray_bcac_public_release_oct17 (1).txt.gz'
gunzip 'oncoarray_bcac_public_release_oct17 (1).txt.gz'

wc -l 'oncoarray_bcac_public_release_oct17 (1).txt' # 11,792,43
head -n 2 'oncoarray_bcac_public_release_oct17 (1).txt'

# want to choose the subset that had (14,910 cases, 17,588 controls), from: oncoarray_bcac_summary_statistics_file_fields.xlsx
#bcac_gwas_all_eaf_controls
#bcac_gwas_all_beta    $24
#bcac_gwas_all_se      $25
#bcac_gwas_all_P1df    $26
 
#   1                                                       2                            3           4                  5             6                  7                                      8                                 9                     10                           11                     12                         13                             14                    15                        16                 17                           18                     19                20                        21                    22                                 23                      24                            25                       26
#var_name                                              phase3_1kg_id                    chr     position_b37           a0             a1      bcac_onco_icogs_gwas_eaf_controls       bcac_onco_icogs_gwas_beta       bcac_onco_icogs_gwas_se   bcac_onco_icogs_gwas_P1df       bcac_onco2_r2   bcac_onco2_eaf_controls         bcac_onco2_beta                   bcac_onco2_se   bcac_onco2_P1df_Wald    bcac_onco2_P1df_LRT       bcac_icogs2_r2         bcac_icogs2_eaf_controls   bcac_icogs2_beta      bcac_icogs2_se     bcac_icogs2_P1df_Wald   bcac_icogs2_P1df_LRT        bcac_gwas_all_eaf_controls        bcac_gwas_all_beta           bcac_gwas_all_se        bcac_gwas_all_P1df      bcac_onco_icogs_gwas_erpos_eaf_controls bcac_onco_icogs_gwas_erpos_beta   bcac_onco_icogs_gwas_erpos_se   bcac_onco_icogs_gwas_erpos_P1df bcac_onco2_erpos_r2     bcac_onco2_erpos_eaf_controls   bcac_onco2_erpos_beta     bcac_onco2_erpos_se     bcac_onco2_erpos_P1df_Wald      bcac_onco2_erpos_P1df_LRT       bcac_icogs2_erpos_r2    bcac_icogs2_erpos_eaf_controls    bcac_icogs2_erpos_beta  bcac_icogs2_erpos_se    bcac_icogs2_erpos_P1df_Wald     bcac_icogs2_erpos_P1df_LRT      bcac_gwas_erpos_eaf_controls    bcac_gwas_erpos_beta      bcac_gwas_erpos_se      bcac_gwas_erpos_P1df    bcac_onco_icogs_gwas_erneg_eaf_controls bcac_onco_icogs_gwas_erneg_beta bcac_onco_icogs_gwas_erneg_se     bcac_onco_icogs_gwas_erneg_P1df bcac_onco2_erneg_r2     bcac_onco2_erneg_eaf_controls   bcac_onco2_erneg_beta   bcac_onco2_erneg_se       bcac_onco2_erneg_P1df_Wald      bcac_onco2_erneg_P1df_LRT       bcac_icogs2_erneg_r2    bcac_icogs2_erneg_eaf_controls  bcac_icogs2_erneg_beta  bcac_icogs2_erneg_se      bcac_icogs2_erneg_P1df_Wald     bcac_icogs2_erneg_P1df_LRT      bcac_gwas_erneg_eaf_controls    bcac_gwas_erneg_beta    bcac_gwas_erneg_se        bcac_gwas_erneg_P1df
#1_10616_CCGCCGTTGCAAAGGCGCGCCG_C        rs376342519:10616:CCGCCGTTGCAAAGGCGCGCCG:C      1          10616     CCGCCGTTGCAAAGGCGCGCCG  C                0.9925                         0.0001                                    0.0578              0.9987                        0.837659               0.9925                      0.00009                       0.05784                 0.9987                   0.9987             0.26246                     0.994742              0.01915              0.1312                   0.88393              0.88394                              NULL                      NULL                       NULL                      NULL      0.9925  -0.0437 0.0645  0.4983  0.841136        0.992544        -0.04371        0.06455 0.49827 0.49848 0.26206 0.994756        -0.01889        0.15469   0.9028  0.90283 NULL    NULL    NULL    NULL    0.9925  0.1702  0.1101  0.1221  0.841772        0.992544        0.17021 0.1101  0.1221  0.11615 0.261928  0.994743        0.02675 0.25431 0.91624 0.9161  NULL    NULL    NULL    NULL



#   1                  2          3           4          5        6                    7                                    8                              9                         10                        11                12                       13                             14                      15                    16                 17                                 18                            19                    20                         21                           22                     23                      24
#var_name        phase3_1kg_id   chr     position_b37    a0      a1      bcac_onco_icogs_gwas_eaf_controls       bcac_onco_icogs_gwas_beta       bcac_onco_icogs_gwas_se bcac_onco_icogs_gwas_P1df       bcac_onco2_r2   bcac_onco2_eaf_controls bcac_onco2_betabcac_onco2_se    bcac_onco2_P1df_Wald    bcac_onco2_P1df_LRT     bcac_icogs2_r2  bcac_icogs2_eaf_controls        bcac_icogs2_betabcac_icogs2_se  bcac_icogs2_P1df_Wald   bcac_icogs2_P1df_LRT    bcac_gwas_all_eaf_controls      bcac_gwas_all_beta      bcac_gwas_all_se        bcac_gwas_all_P1df      bcac_onco_icogs_gwas_erpos_eaf_controls bcac_onco_icogs_gwas_erpos_beta bcac_onco_icogs_gwas_erpos_se   bcac_onco_icogs_gwas_erpos_P1df bcac_onco2_erpos_r2     bcac_onco2_erpos_eaf_controls   bcac_onco2_erpos_beta  bcac_onco2_erpos_se      bcac_onco2_erpos_P1df_Wald      bcac_onco2_erpos_P1df_LRT       bcac_icogs2_erpos_r2    bcac_icogs2_erpos_eaf_controls  bcac_icogs2_erpos_beta  bcac_icogs2_erpos_se    bcac_icogs2_erpos_P1df_Wald     bcac_icogs2_erpos_P1df_LRT     bcac_gwas_erpos_eaf_controls     bcac_gwas_erpos_beta    bcac_gwas_erpos_se      bcac_gwas_erpos_P1df    bcac_onco_icogs_gwas_erneg_eaf_controls bcac_onco_icogs_gwas_erneg_beta bcac_onco_icogs_gwas_erneg_se   bcac_onco_icogs_gwas_erneg_P1df bcac_onco2_erneg_r2     bcac_onco2_erneg_eaf_controls   bcac_onco2_erneg_beta   bcac_onco2_erneg_se     bcac_onco2_erneg_P1df_Wald      bcac_onco2_erneg_P1df_LRT       bcac_icogs2_erneg_r2    bcac_icogs2_erneg_eaf_controls  bcac_icogs2_erneg_beta  bcac_icogs2_erneg_se   bcac_icogs2_erneg_P1df_Wald      bcac_icogs2_erneg_P1df_LRT      bcac_gwas_erneg_eaf_controls    bcac_gwas_erneg_beta    bcac_gwas_erneg_se      bcac_gwas_erneg_P1df
#1_10616_CCGCCGTTGCAAAGGCGCGCCG_C        rs376342519:10616:CCGCCGTTGCAAAGGCGCGCCG:C      1       10616   CCGCCGTTGCAAAGGCGCGCCG C0.9925  0.0001  0.0578  0.9987  0.837659        0.9925  0.00009 0.05784 0.9987  0.9987  0.26246 0.994742        0.01915 0.1312 0.88393  0.88394 NULL    NULL    NULL    NULL    0.9925  -0.0437 0.0645  0.4983  0.841136        0.992544        -0.04371       0.06455  0.49827 0.49848 0.26206 0.994756        -0.01889        0.15469 0.9028  0.90283 NULL    NULL    NULL    NULL    0.9925 0.1702   0.1101  0.1221  0.841772        0.992544        0.17021 0.1101  0.1221  0.11615 0.261928        0.994743        0.026750.25431  0.91624 0.9161  NULL    NULL    NULL    NULL
# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N

# convert to common format
#cases= 14910
#controls=17588
#N_eff = 4 / (1 / cases + 1 / controls) 
#N_eff # 32277.32
avgNumIndis=32277 # 32498  #14910  + 17588
awk -v avgNumIndis="$avgNumIndis" ' 
FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { 
split($2,rsid,":")
print $3"\t"$4"\t"rsid[1]"\t"$6"\t"$5"\tX\t"$24"\t"$25"\t"$26"\t"avgNumIndis } }
' OFS='\t' 'oncoarray_bcac_public_release_oct17 (1).txt' > $eursumstatsLoc$'EUR_BRCA_raw'
head $eursumstatsLoc$'EUR_BRCA_raw'


# get rid of ambigous alleles
remove_ambiguous_alleles $eursumstatsLoc$'EUR_BRCA_raw'


# hapmap 3 subset
awk ' 
FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {
if (FNR == 1) {print $0 }
else { if( $3 in file1 ) {print $0} }
}
' OFS='\t' $hapmap3_b37bim $eursumstatsLoc$'EUR_BRCA_raw_noAmbiguousAlleles' > $eursumstatsLoc$'EUR_BRCA_hm3'


head $eursumstatsLoc$'EUR_BRCA_hm3'
wc -l $eursumstatsLoc$'EUR_BRCA_hm3' # 1094221



grep "$(printf '12\t1')" $eursumstatsLoc$'EUR_BRCA_raw_noAmbiguousAlleles'

# it seems that Chrom12 had QC issues here as none of the RSids were present for chrom12. These were excluded to be on the safe side.
grep "$(printf '12\t1')" $eursumstatsLoc$'EUR_BRCA_hm3'



################################
# JP:
wget http://134.160.84.25:8080/download/Breast_cancer

mv Breast_cancer Breast_cancer.gz
gunzip Breast_cancer.gz

head Breast_cancer
wc -l Breast_cancer # 6,045,460

#   1     2        3        4        5          6          7       8         9
# chrom   pos     ref     alt     rsids   nearest_genes   pval    beta    sebeta maf
#1       751343  T       A       rs28544273      AL669831.1      0.049   -0.062 0.031    0.14


# convert to common format
#cases= 5552
#controls=89731
#N_eff = 4 / (1 / cases + 1 / controls) 
#N_eff # 20914

avgNumIndis=20914  # 5,552 + 89,731
awk -v avgNumIndis="$avgNumIndis" ' 
FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { print $1"\t"$2"\t"$5"\t"$4"\t"$3"\tX\t"$8"\t"$9"\t"$7"\t"avgNumIndis }
}
' OFS='\t' Breast_cancer > $japsumstatsLoc$'JP_BRCA_raw'
head $japsumstatsLoc$'JP_BRCA_raw'


# get rid of ambigous alleles
remove_ambiguous_alleles $japsumstatsLoc$'JP_BRCA_raw'


# hapmap 3 subset
awk ' 
FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {
if (FNR == 1) {print $0 }
else { if( $3 in file1 ) {print $0} }
}
' OFS='\t' $hapmap3_b37bim $japsumstatsLoc$'JP_BRCA_raw_noAmbiguousAlleles' > $japsumstatsLoc$'JP_BRCA_hm3'

head $japsumstatsLoc$'JP_BRCA_hm3'
wc -l $japsumstatsLoc$'JP_BRCA_hm3' # 925322



#################################
# CAD:
mkdir -p $crossAncestryRaw$'CAD/'
cd $crossAncestryRaw$'CAD/'


#EUR: this uses the Interim UKBB, IE those will have to be excluded
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004787/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz
gunzip UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz

head -n 100 UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt
wc -l UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt # 9,026,568

# the main paper says of the sample size for some of the meta analysis
# 71,602 cases and 260,875 controls (for exome markers 53,135 and 215,611
# for the discovery samples it says cases= 10801 and controls=137914
# but the released data has reported sample sizes that range from 150K to almost 300K
# to get a sensible N_eff, I average the case:controls ratio from the information and 'impute' 
#71602/(260875 + 71602) # 0.2153593
#53135/(215611+53135) #  0.1977146
#10801/(137914 + 10801) # 0.07262885  # but this is 150K, almost none of the actual data is below 200K
# I settled to 20% cases
# find out the sample sizes 
#awk '{count[$11]++} END {for (word in count) print word, count[word]}' UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt

#     1              2            3        4          5                6                        7                  8       9       10               11
#Markername      snptestid       chr     bp_hg19 effect_allele   noneffect_allele        effect_allele_freq      logOR   se_gc  p-value_gc       n_samples       exome   info_ukbb
#1:569406_G_A    rs561255355      1       569406          G              A                    0.99858           0.05191 0.27358 0.849496        154959            yes     0.41
# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N

#cases= 10801
#controls=137914
#N_eff = 4 / (1 / cases + 1 / controls) 
#round(N_eff) # 40066

# convert to common formata
awk ' 
FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { 
cases= $11 * 0.2
controls=$11 * 0.8
N_eff = int( 4 / (1 / cases + 1 / controls) )

print $3"\t"$4"\t"$2"\t"$5"\t"$6"\tX\t"$8"\t"$9"\t"$10"\t" N_eff

}
}
' OFS='\t' UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt > $eursumstatsLoc$'EUR_CAD_raw'
head $eursumstatsLoc$'EUR_CAD_raw'


# get rid of ambigous alleles
remove_ambiguous_alleles $eursumstatsLoc$'EUR_CAD_raw'


# hapmap 3 subset
awk ' 
FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {
if (FNR == 1) {print $0 }
else { if( $3 in file1 ) {print $0} }
}
' OFS='\t' $hapmap3_b37bim $eursumstatsLoc$'EUR_CAD_raw_noAmbiguousAlleles' > $eursumstatsLoc$'EUR_CAD_hm3'
head $eursumstatsLoc$'EUR_CAD_hm3'
wc -l $eursumstatsLoc$'EUR_CAD_hm3' # 1127191


##############################################
# JP: 
wget http://134.160.84.25:8080/download/Coronary_artery_disease

mv Coronary_artery_disease Coronary_artery_disease.gz
gunzip Coronary_artery_disease.gz

head Coronary_artery_disease
wc -l Coronary_artery_disease # 6046444

#chrom   pos     ref     alt     rsids   nearest_genes   pval    beta    sebeta maf
#1       751343  T       A       rs28544273      AL669831.1      0.75    -0.00480.015    0.14
# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N


# convert to common format
#cases= 29319
#controls=183134
#N_eff = 4 / (1 / cases + 1 / controls) 
#round(N_eff) # 
avgNumIndis=101092  # 29,319	+ 183,134 
awk -v avgNumIndis="$avgNumIndis" ' 
FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else {  print $1"\t"$2"\t"$5"\t"$4"\t"$3"\tX\t"$8"\t"$9"\t"$7"\t"avgNumIndis  }
}
' OFS='\t' Coronary_artery_disease > $japsumstatsLoc$'JP_CAD_raw'
head $japsumstatsLoc$'JP_CAD_raw'


# get rid of ambigous alleles
remove_ambiguous_alleles $japsumstatsLoc$'JP_CAD_raw'


# hapmap 3 subset
awk ' 
FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {
if (FNR == 1) {print $0 }
else { if( $3 in file1 ) {print $0} }
}
' OFS='\t' $hapmap3_b37bim $japsumstatsLoc$'JP_CAD_raw_noAmbiguousAlleles' > $japsumstatsLoc$'JP_CAD_hm3'
head $japsumstatsLoc$'JP_CAD_hm3'
wc -l $japsumstatsLoc$'JP_CAD_hm3' # 925414


#################################
# Asthma:
mkdir -p $crossAncestryRaw$'asthma/'
cd $crossAncestryRaw$'asthma/'

# EUR:
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006862/TAGC_meta-analyses_results_for_asthma_risk.zip

wc -l $crossAncestryRaw$'asthma/TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv' # 2,001,281
head -n 2 $crossAncestryRaw$'asthma/TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv' # 2,001,281

# 1        2        3                   4                       5                         6                      7                      
#chr     rsid    position        reference_allele        alternate_allele        Multiancestry_beta_fix  Multiancestry_se_fix   Multiancestry_pval_fix   Multiancestry_beta_rand Multiancestry_se_rand   Multiancestry_pval_rand Multiancestry_HetQtest  Multiancestry_df_HetQtest       Multiancestry_pval_HetQtest     European_ancestry_beta_fix      European_ancestry_se_fix        European_ancestry_pval_fix      European_ancestry_beta_rand     European_ancestry_se_rand       European_ancestry_pval_rand     European_ancestry_HetQtest      European_ancestry_df_HetQtest   European_ancestry_pval_HetQtest
#1       rs3094315       752566  G       A       -.00830888      .02015612       .6801735        -.00830888      .02015612      .6801735 36.575  51      .9361117        -.0152675       .02376994       .5206766        -.0152675       .02376994       .520676630.743  42      .900529
# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N

# using the European subset estimates: Description_TAGC file_meta-analyses results.xlsx
#European_ancestry_beta_fix  $15
#European_ancestry_se_fix    $16
#European_ancestry_pval_fix  $17

# convert to common format
#cases= 19954
#controls=107715
#N_eff = 4 / (1 / cases + 1 / controls) 
#round(N_eff) # 67341

avgNumIndis=67341  # 19,954 European ancestry cases, 107,715 European ancestry controls
awk -v avgNumIndis="$avgNumIndis" ' 
FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else {  print $1"\t"$3"\t"$2"\t"$5"\t"$4"\tX\t"$15"\t"$16"\t"$17"\t"avgNumIndis  }
}
' OFS='\t' $crossAncestryRaw$'asthma/TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv' > $eursumstatsLoc$'EUR_asthma_raw'
head $eursumstatsLoc$'EUR_asthma_raw'


# get rid of ambigous alleles
remove_ambiguous_alleles $eursumstatsLoc$'EUR_asthma_raw'


# hapmap 3 subset
awk ' 
FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {
if (FNR == 1) {print $0 }
else { if( $3 in file1 ) {print $0} }
}
' OFS='\t' $hapmap3_b37bim $eursumstatsLoc$'EUR_asthma_raw_noAmbiguousAlleles' > $eursumstatsLoc$'EUR_asthma_hm3'
head $eursumstatsLoc$'EUR_asthma_hm3'
wc -l $eursumstatsLoc$'EUR_asthma_hm3' # 995413







#########################
# JP:
wget http://134.160.84.25:8080/download/Asthma

mv Asthma Asthma.gz
gunzip Asthma.gz

head Asthma
wc -l Asthma # 6143488


#chrom   pos     ref     alt     rsids   nearest_genes   pval    beta    sebeta maf
#1       751343  T       A       rs28544273      AL669831.1      0.25    -0.029 0.026    0.14

# convert to common format
#cases= 8216
#controls=201592
#N_eff = 4 / (1 / cases + 1 / controls) 
#round(N_eff) # 31577

avgNumIndis=31577  # 8216+201592 # 8,216 [M=3,932, F=4,284]	201,592 [M=103,089, F=98,503]
awk -v avgNumIndis="$avgNumIndis" ' 
FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else {  print $1"\t"$2"\t"$5"\t"$4"\t"$3"\tX\t"$8"\t"$9"\t"$7"\t"avgNumIndis  }
}
' OFS='\t' Asthma > $japsumstatsLoc$'JP_asthma_raw'
head $japsumstatsLoc$'JP_asthma_raw'




# get rid of ambigous alleles
remove_ambiguous_alleles $japsumstatsLoc$'JP_asthma_raw'


# hapmap 3 subset
awk ' 
FNR == NR { file1[ $2 ] = $2; next; } 
FNR <= NR {
if (FNR == 1) {print $0 }
else { if( $3 in file1 ) {print $0} }
}
' OFS='\t' $hapmap3_b37bim $japsumstatsLoc$'JP_asthma_raw_noAmbiguousAlleles' > $japsumstatsLoc$'JP_asthma_hm3'
head $japsumstatsLoc$'JP_asthma_hm3'
wc -l $japsumstatsLoc$'JP_asthma_hm3' # 944065


###################################



$eursumstatsLoc$'EUR_T2D_hm3'
$japsumstatsLoc$'JP_T2D_hm3'

$eursumstatsLoc$'EUR_height_hm3'
$japsumstatsLoc$'JP_height_hm3'

$eursumstatsLoc$'EUR_BRCA_hm3'
$japsumstatsLoc$'JP_BRCA_hm3'

$eursumstatsLoc$'EUR_CAD_hm3' 
$japsumstatsLoc$'JP_CAD_hm3'

$eursumstatsLoc$'EUR_asthma_hm3'
$japsumstatsLoc$'JP_asthma_hm3'




# Ia)  run shaPRS to produce the coefs and the lFDRs - EUR_JAP
pheno='EUR_JAP_asthma'
shaPRS_new $eursumstatsLoc$'EUR_asthma_hm3' $japsumstatsLoc$'JP_asthma_hm3' $crossAncestryResults$pheno '0' $popALDpanel $popBLDpanel

pheno='EUR_JAP_height'
shaPRS_new $eursumstatsLoc$'EUR_height_hm3' $japsumstatsLoc$'JP_height_hm3' $crossAncestryResults$pheno '0' $popALDpanel $popBLDpanel

pheno='EUR_JAP_T2D'
shaPRS_new $eursumstatsLoc$'EUR_T2D_hm3' $japsumstatsLoc$'JP_T2D_hm3' $crossAncestryResults$pheno '0' $popALDpanel $popBLDpanel

pheno='EUR_JAP_BRCA'
shaPRS_new $eursumstatsLoc$'EUR_BRCA_hm3' $japsumstatsLoc$'JP_BRCA_hm3' $crossAncestryResults$pheno '0' $popALDpanel $popBLDpanel

pheno='EUR_JAP_CAD'
shaPRS_new $eursumstatsLoc$'EUR_CAD_hm3' $japsumstatsLoc$'JP_CAD_hm3' $crossAncestryResults$pheno '0' $popALDpanel $popBLDpanel


# (wait for above to finish on cluster)
# convert LDpredformat to PRS-CS format:
pheno='EUR_JAP_asthma'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eur/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$asthmaEURRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='EUR_JAP_height'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eur/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$heightEURRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='EUR_JAP_T2D'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eur/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$t2dEURRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='EUR_JAP_BRCA'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eur/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$brcaEURRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='EUR_JAP_CAD'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eur/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$cadEURRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"


###############
## pheno='EUR_JAP_height'
## arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eur/ /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/ '$prscsrefs$'ldblk_ukbb_'$heightEURRef$'orig/'
## bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

#PRSCSRefLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsrefs/ldblk_ukbb_eur/"
#LDpred2Loc = "/lustre/scratch123/hgi/Lmdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/"
#outputLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/raw/prscsrefs/ldblk_ukbb_heighteur_orig/"

## pheno='EUR_JAP_height'
## arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v5.R '$prscsrefs$'ldblk_ukbb_eur/ /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/ '$prscsrefs$'ldblk_ukbb_'$heightEURRef$'orig2/'
## bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"


###############


# Ib)  run shaPRS to produce the coefs and the lFDRs - JAP_EUR
# (IE running the other way around)
pheno='JAP_EUR_asthma'
shaPRS_new $japsumstatsLoc$'JP_asthma_hm3' $eursumstatsLoc$'EUR_asthma_hm3' $crossAncestryResults$pheno '0' $popBLDpanel $popALDpanel

pheno='JAP_EUR_height'
shaPRS_new $japsumstatsLoc$'JP_height_hm3' $eursumstatsLoc$'EUR_height_hm3' $crossAncestryResults$pheno '0' $popBLDpanel $popALDpanel

pheno='JAP_EUR_T2D'
shaPRS_new $japsumstatsLoc$'JP_T2D_hm3' $eursumstatsLoc$'EUR_T2D_hm3' $crossAncestryResults$pheno '0' $popBLDpanel $popALDpanel


pheno='JAP_EUR_BRCA'
shaPRS_new $japsumstatsLoc$'JP_BRCA_hm3' $eursumstatsLoc$'EUR_BRCA_hm3' $crossAncestryResults$pheno '0' $popBLDpanel $popALDpanel


pheno='JAP_EUR_CAD'
shaPRS_new $japsumstatsLoc$'JP_CAD_hm3' $eursumstatsLoc$'EUR_CAD_hm3' $crossAncestryResults$pheno '0' $popBLDpanel $popALDpanel



# (wait for above to finish on cluster)
# convert LDpredformat to PRS-CS format:
# NOTE: we are using the PRS-CS EUR format (ldblk_ukbb_eur) as a template, as we want the EUR LD-blocks  (and not the ldblk_ukbb_eas)
pheno='JAP_EUR_asthma'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eur/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$asthmaEASRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='JAP_EUR_height'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eur/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$heightEASRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='JAP_EUR_T2D'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eur/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$t2dEASRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='JAP_EUR_BRCA'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eur/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$brcaEASRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='JAP_EUR_CAD'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eur/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$cadEASRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

#####
## Try it with EAS PRS-CS LD panel

pheno='JAP_EUR_asthma'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eas/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$asthmaEASRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='JAP_EUR_height'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eas/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$heightEASRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='JAP_EUR_T2D'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eas/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$t2dEASRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='JAP_EUR_BRCA'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eas/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$brcaEASRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

pheno='JAP_EUR_CAD'
arguments='/nfs/users/nfs_m/mk23/scripts/ConvertLDpred2ToPRSCS_v4.R '$prscsrefs$'ldblk_ukbb_eas/ '$crossAncestryResults$pheno$'/ '$prscsrefs$'ldblk_ukbb_'$cadEASRef$'/'
bsub -q normal -m "modern_hardware" -G team152 -n7 -J ${pheno}_LDCONV -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $crossAncestryResults$pheno$'_LDCONV.out' -e $crossAncestryResults$pheno$'_LDCONV.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"



head $crossAncestryResults$'EUR_JAP_height_lFDR_meta_SNP_lFDR'
head $crossAncestryResults$'EUR_JAP_height_sumstats_meta'

head $crossAncestryResults$'EUR_JAP_height_lFDR_meta_SNP_lFDR_old'
head $crossAncestryResults$'EUR_JAP_height_sumstats_meta_old'



############################################################################




j=3


PRS_array=( 'EUR_JAP_asthma' 'EUR_JAP_height' 'EUR_JAP_BRCA' 'EUR_JAP_CAD' 'EUR_JAP_T2D' )
EUR_array=( 'EUR_asthma_hm3' 'EUR_height_hm3' 'EUR_BRCA_hm3' 'EUR_CAD_hm3' 'EUR_T2D_hm3' )
JPT_array=( 'JP_asthma_hm3' 'JP_height_hm3' 'JP_BRCA_hm3' 'JP_CAD_hm3' 'JP_T2D_hm3' )

PRS_binary=( '1' '0' '1' '1' '1' )
arraylength=${#PRS_array[@]}


for (( j=1; j<${arraylength}+1; j++ )); do
pheno=${PRS_array[$j-1]} 
binPheno=${PRS_binary[$j-1]} 
EURPheno=${EUR_array[$j-1]} 
JPTPheno=${JPT_array[$j-1]} 


# visualise
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_qval_manhattan.R '$crossAncestryResults$pheno$'_SE_meta '$crossAncestryResults$pheno$'_lFDR_meta_SNP_lFDR '$crossAncestryResults$pheno$' '$pheno$' '$crossAncestryResults$pheno$'_sumstats_meta'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

#done


# II) perform QC separately on each
#$eursumstatsLoc$EURPheno
#$japsumstatsLoc$JPTPheno
#$crossAncestryResults$pheno$'_sumstats_meta_combined'
#$crossAncestryResults$pheno$'_sumstats_meta'

arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_QC.R '$crossAncestryResults$pheno$'/ '$crossAncestryResults$pheno$'_sumstats_meta '$crossAncestryResults$pheno$'_sumstats_meta_QC_keep '$binPheno
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments


# subset the SNPs to be the same panel for all 3 sumstats 
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums) {print $0 } }
' $crossAncestryResults$pheno$'_sumstats_meta_QC_keep' $crossAncestryResults$pheno$'_sumstats_meta_combined'  > $crossAncestryResults$pheno$'_sumstats_meta_combined_QC'

awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $crossAncestryResults$pheno$'_sumstats_meta_QC_keep' $crossAncestryResults$pheno$'_sumstats_meta'  > $crossAncestryResults$pheno$'_sumstats_meta_QC'

awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums ) {print $0 } }
' $crossAncestryResults$pheno$'_sumstats_meta_QC_keep' $eursumstatsLoc$EURPheno  > $eursumstatsLoc$EURPheno$'_QC'

# use the surviving meta/shaPRS SNP list to subset the japanese-only sumstats list so that it works off the same final list of SNPs as the combined/shaPRS
awk 'FNR == NR { sums[ $3 ] = $3; next; } FNR <= NR { if( FNR == 1 || $3 in sums ) {print $0 } }
' $crossAncestryResults$pheno$'_sumstats_meta_QC' $japsumstatsLoc$JPTPheno  > $japsumstatsLoc$JPTPheno$'_QC'



#######################################################################


# IIIa) Build PRS via LDpred2
outfile=$crossAncestryResults$'EUR_JAP_LD_'$pheno
#mv $outfile $outfile$'_empirical' # backup originals
#rm -rf $outfile
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
else
#  JAP  (EUR LD panel)
avgNumIndis=$(tail -n 1 $japsumstatsLoc$JPTPheno | awk '{ print $10}')
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_v5.R /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/ '$japsumstatsLoc$JPTPheno$'_QC '$NCORES_LDPred2$' '$avgNumIndis$' '$crossAncestryResults$' EUR_JAP_LD_'$pheno$' '$crossAncestrySumStats$'all_hm3 '$shaPRSscriptLoc$' '$binPheno
#$Rscript $arguments

#bsub -q long -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J EUR_JAP_LD_${pheno} -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $crossAncestryResults$'EUR_JAP_LD_'$pheno$'.out' -e $crossAncestryResults$'EUR_JAP_LD_'$pheno$'.err'  "$Rscript $arguments"
fi



# IIIa) Build PRS via LDpred2
outfile=$crossAncestryResults$'EUR_EUR_LD_'$pheno
#mv $outfile $outfile$'_empirical' # backup originals
#rm -rf $outfile
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
else
#  EUR  (EUR LD panel)
avgNumIndis=$(tail -n 1 $eursumstatsLoc$EURPheno | awk '{ print $10}')
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_v5.R /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/ '$eursumstatsLoc$EURPheno$'_QC '$NCORES_LDPred2$' '$avgNumIndis$' '$crossAncestryResults$' EUR_EUR_LD_'$pheno$' '$crossAncestrySumStats$'all_hm3 '$shaPRSscriptLoc$' '$binPheno
#$Rscript $arguments

#bsub -q long -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J EUR_EUR_LD_${pheno} -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $crossAncestryResults$'EUR_EUR_LD_'$pheno$'.out' -e $crossAncestryResults$'EUR_EUR_LD_'$pheno$'.err'  "$Rscript $arguments"
fi


outfile=$crossAncestryResults$'Combined_EUR_LD_'$pheno
#mv $outfile $outfile$'_empirical' # backup originals
#rm -rf $outfile
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
else
#  EUR-Combined (EUR+JP) (EUR LD panel)
avgNumIndis=$(tail -n 1 "$crossAncestryResults$pheno"_sumstats_meta_combined | awk '{ print $10}')
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_v5.R /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/ '$crossAncestryResults$pheno'_sumstats_meta_combined_QC '$NCORES_LDPred2$' '$avgNumIndis$' '$crossAncestryResults$' Combined_EUR_LD_'$pheno$' '$crossAncestrySumStats$'all_hm3 '$shaPRSscriptLoc$' '$binPheno
#bsub -q long -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J Combined_EUR_LD_${pheno} -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $crossAncestryResults$'Combined_EUR_LD_'$pheno$'.out' -e $crossAncestryResults$'Combined_EUR_LD_'$pheno$'.err'  "$Rscript $arguments"
fi

outfile=$crossAncestryResults$'shaPRS_LD_'$pheno
#mv $outfile $outfile$'_empirical' # backup originals
#rm -rf $outfile
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
else
# EUR-shaPRS shaPRS LD Panel
avgNumIndis=$(tail -n 1 "$crossAncestryResults$pheno"_sumstats_meta | awk '{ print $10}')
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_v5.R '$crossAncestryResults$pheno$'/ '$crossAncestryResults$pheno'_sumstats_meta_QC '$NCORES_LDPred2$' '$avgNumIndis$' '$crossAncestryResults$' shaPRS_LD_'$pheno$' '$crossAncestrySumStats$'all_hm3 '$shaPRSscriptLoc$' '$binPheno
#bsub -q long -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J shaPRS_LD_${pheno} -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $crossAncestryResults$'shaPRS_LD_'$pheno$'.out' -e $crossAncestryResults$'shaPRS_LD_'$pheno$'.err'  "$Rscript $arguments"
fi


outfile=$crossAncestryResults$'shaPRS_EUR_LD_'$pheno
#mv $outfile $outfile$'_empirical' # backup originals
#rm -rf $outfile
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
else
# EUR-shaPRS EUR LD Panel
avgNumIndis=$(tail -n 1 "$crossAncestryResults$pheno"_sumstats_meta | awk '{ print $10}')
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_v5.R /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/ '$crossAncestryResults$pheno'_sumstats_meta_QC '$NCORES_LDPred2$' '$avgNumIndis$' '$crossAncestryResults$' shaPRS_EUR_LD_'$pheno$' '$crossAncestrySumStats$'all_hm3 '$shaPRSscriptLoc$' '$binPheno
#bsub -q long -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J shaPRS_EUR_LD_${pheno} -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $crossAncestryResults$'shaPRS_EUR_LD_'$pheno$'.out' -e $crossAncestryResults$'shaPRS_EUR_LD_'$pheno$'.err'  "$Rscript $arguments"
fi

outfile=$crossAncestryResults$'Combined_LD_'$pheno
#mv $outfile $outfile$'_empirical' # backup originals
#rm -rf $outfile
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
else
#   EUR-Combined shaPRS LD Panel
avgNumIndis=$(tail -n 1 "$crossAncestryResults$pheno"_sumstats_meta_combined | awk '{ print $10}')
arguments='/nfs/users/nfs_m/mk23/scripts/LDpred2_auto_v5.R '$crossAncestryResults$pheno$'/ '$crossAncestryResults$pheno'_sumstats_meta_combined_QC '$NCORES_LDPred2$' '$avgNumIndis$' '$crossAncestryResults$' Combined_LD_'$pheno$' '$crossAncestrySumStats$'all_hm3 '$shaPRSscriptLoc$' '$binPheno
#bsub -q long -m "modern_hardware" -G team152 -n${NCORES_LDPred2} -J Combined_LD_${pheno} -R "span[hosts=1] select[mem>${ldPred2Mem}] rusage[mem=${ldPred2Mem}]" -M${ldPred2Mem} -o $crossAncestryResults$'Combined_LD_'$pheno$'.out' -e $crossAncestryResults$'Combined_LD_'$pheno$'.err'  "$Rscript $arguments"
fi



done



# !!! HERE

############################################################################

# evaluate Asthma / Height on the test set




# Asthma
pheno='EUR_JAP_asthma'
evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'EUR_JAP_LD_'$pheno '1'    # Theoretical: correlation_sq 0.0015 / AUC 0.5351 653000 valid predictors loaded        # Empirical: correlation_sq 0.0015 / AUC 0.5351 654276 valid predictors loaded
evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'EUR_EUR_LD_'$pheno '1'      # Theoretical:  correlation_sq 0.0113 / AUC 0.5947  # 652999 valid predictors loaded. # Empirical: correlation_sq 0.0113 / AUC 0.5946  # 654275 valid predictors loaded.
evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'Combined_EUR_LD_'$pheno '1' # Theoretical: correlation_sq 0.0112 / AUC 0.5938 # 652999 valid predictors loaded.   # Empirical: correlation_sq 0.0117 / AUC 0.5959 # 654275 valid predictors loaded.
evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'Combined_LD_'$pheno '1'     # Theoretical: correlation_sq  0.0135 / AUC 0.6033 # 652999 valid predictors loaded.  # Empirical: correlation_sq 0.0135 / AUC 0.6029 # 654275 valid predictors loaded.
evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'shaPRS_EUR_LD_'$pheno '1'   # Theoretical: correlation_sq 0.0129 / AUC 0.601 # 652999 valid predictors loaded.    # Empirical: correlation_sq 0.0131 / AUC 0.6016 # 654275 valid predictors loaded.
evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'shaPRS_LD_'$pheno '1'       # Theoretical: correlation_sq  0.0136 / AUC 0.6035 # 652999 valid predictors loaded.  # Empirical: correlation_sq 0.0136 / AUC 0.6035 # 654275 valid predictors loaded.


awk '{count[$6]++} END {for (word in count) print word, count[word]}' $crossAncestrySumStats$'all_hm3.fam'
# 1 222649
# 2 28576


# height
pheno='EUR_JAP_height'
evaluatePRS $crossAncestrySumStats$'all_height_hm3.fam' 'EUR_JAP_LD_'$pheno '0'      # Theoretical: correlation_sq 0.0413 # 597805 valid predictors loaded. # Empirical: correlation_sq 0.0415 # 609645 valid predictors loaded.
evaluatePRS $crossAncestrySumStats$'all_height_hm3.fam' 'EUR_EUR_LD_'$pheno '0'      # Theoretical: correlation_sq 0.105  # 597805 valid predictors loaded. # Empirical: correlation_sq 0.1048 # 609645 valid predictors loaded.
evaluatePRS $crossAncestrySumStats$'all_height_hm3.fam' 'Combined_EUR_LD_'$pheno '0' # Theoretical: correlation_sq 0.107 # 597805 valid predictors loaded   # Empirical: correlation_sq 0.1071 # 609645 valid predictors loaded
evaluatePRS $crossAncestrySumStats$'all_height_hm3.fam' 'Combined_LD_'$pheno '0'     # Theoretical: correlation_sq 0.1086 # 597805 valid predictors loaded. # Empirical: correlation_sq 0.1087 # 609645 valid predictors loaded.
evaluatePRS $crossAncestrySumStats$'all_height_hm3.fam' 'shaPRS_EUR_LD_'$pheno '0'   # Theoretical: correlation_sq 0.115 # 597805 valid predictors loaded.  # Empirical: correlation_sq 0.1152 # 609645 valid predictors loaded.
evaluatePRS $crossAncestrySumStats$'all_height_hm3.fam' 'shaPRS_LD_'$pheno '0'       # Theoretical: correlation_sq 0.1147 # 597805 valid predictors loaded. # Empirical: correlation_sq 0.1148 # 609645 valid predictors loaded


pheno='EUR_JAP_T2D' # EUR_JAP_BRCA  EUR_JAP_CAD EUR_JAP_T2D
wc -l $crossAncestryResults$'EUR_JAP_LD_'$pheno



# Elena's Theoretical results for BRCA/CAD/T2D
# "model" "AUC" "r2"
# "Combined_EUR_LD_EUR_JAP_BRCA_E" 0.606819493237408 0.00953641682342133
# "Combined_EUR_LD_EUR_JAP_CAD_E" 0.668158569004927 0.0176812362973048
# "Combined_EUR_LD_EUR_JAP_T2D_E" 0.635826169039535 0.00920102295950981
# "Combined_LD_EUR_JAP_BRCA_E" 0.607344628946614 0.00961922520781916
# "Combined_LD_EUR_JAP_CAD_E" 0.656116385738032 0.014845980290117
# "Combined_LD_EUR_JAP_T2D_E" 0.655301330347591 0.0120658681370776

# "EUR_EUR_LD_EUR_JAP_BRCA_E" 0.603056788176649 0.00885546221458986
# "EUR_EUR_LD_EUR_JAP_CAD_E" 0.631309141041285 0.0103268715102956
# "EUR_EUR_LD_EUR_JAP_T2D_E" 0.62533866409085 0.00780336300593319
# "EUR_JAP_LD_EUR_JAP_BRCA_E" 0.55067731752388 0.00213075418284815
# "EUR_JAP_LD_EUR_JAP_CAD_E" 0.568755367460607 0.00267550882281642
# "EUR_JAP_LD_EUR_JAP_T2D_E" 0.576982419165456 0.00291110430794624

# "shaPRS_EUR_LD_EUR_JAP_BRCA_E" 0.609018112432376 0.00987836291371829
# "shaPRS_EUR_LD_EUR_JAP_CAD_E" 0.679012693889687 0.0203131524696627
# "shaPRS_EUR_LD_EUR_JAP_T2D_E" 0.638447899961932 0.00956662078211792
# "shaPRS_LD_EUR_JAP_BRCA_E" 0.609401260361629 0.0099414210489418
# "shaPRS_LD_EUR_JAP_CAD_E" 0.669828139766618 0.0178648174887799
# "shaPRS_LD_EUR_JAP_T2D_E" 0.640787244647595 0.00988095961051051

# "theor_Combined_LD_theor_EUR_JAP_BRCA_E" 0.606779829565598 0.00952516633921873
# "theor_Combined_LD_theor_EUR_JAP_CAD_E" 0.655996300293071 0.0148242366546092
# "theor_Combined_LD_theor_EUR_JAP_T2D_E" 0.655258246625216 0.0120649878521466
# "theor_Combined_theor_EUR_LD_EUR_JAP_CAD_E" 0.668048874994293 0.0176581174617735
# "theor_Combined_theor_EUR_LD_EUR_JAP_T2D_E" 0.635544317965913 0.00916088551097992

# "theor_EUR_EUR_LD_EUR_JAP_BRCA_E" 0.602732396497127 0.00880516857921991
# "theor_EUR_EUR_LD_EUR_JAP_CAD_E" 0.631198432091009 0.0103106406626691
# "theor_EUR_EUR_LD_EUR_JAP_T2D_E" 0.624631771520051 0.00770983808498473
# "theor_EUR_JAP_LD_EUR_JAP_BRCA_E" 0.54658519169381 0.00180469442922215
# "theor_EUR_JAP_LD_EUR_JAP_CAD_E" 0.568774643665759 0.00267643229227253
# "theor_EUR_JAP_LD_EUR_JAP_T2D_E" 0.578254370623175 0.00301015468203018

# "theor_shaPRS_LD_theor_EUR_JAP_BRCA_E" 0.609405344809099 0.00995921998548026
# "theor_shaPRS_LD_theor_EUR_JAP_CAD_E" 0.668923321526022 0.0176544438530993
# "theor_shaPRS_LD_theor_EUR_JAP_T2D_E" 0.644651501488707 0.0104439194571882
# "theor_shaPRS_theor_EUR_LD_EUR_JAP_BRCA_E" 0.608802529227551 0.00985886506845825
# "theor_shaPRS_theor_EUR_LD_EUR_JAP_CAD_E" 0.680052616072609 0.0205056971691544
# "theor_shaPRS_theor_EUR_LD_EUR_JAP_T2D_E" 0.638190457273135 0.00953196890638384

# CAD
pheno='EUR_JAP_CAD'
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'EUR_JAP_LD_'$pheno '1'
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'EUR_EUR_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'Combined_EUR_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'shaPRS_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'Combined_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'shaPRS_EUR_LD_'$pheno '1' # 

# convert to Elena's format 
convertPRS $crossAncestryResults$'EUR_JAP_LD_'$pheno
convertPRS $crossAncestryResults$'EUR_EUR_LD_'$pheno
convertPRS $crossAncestryResults$'Combined_EUR_LD_'$pheno
convertPRS $crossAncestryResults$'shaPRS_LD_'$pheno
convertPRS $crossAncestryResults$'Combined_LD_'$pheno
convertPRS $crossAncestryResults$'shaPRS_EUR_LD_'$pheno

#rm -rf $crossAncestryResults$pheno'_compressed.zip'
#zip $crossAncestryResults$pheno'_compressed.zip' $crossAncestryResults$'EUR_JAP_LD_'$pheno$'_E' $crossAncestryResults$'EUR_EUR_LD_'$pheno$'_E' $crossAncestryResults$'Combined_EUR_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_LD_'$pheno$'_E' $crossAncestryResults$'Combined_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_EUR_LD_'$pheno$'_E'
zip $crossAncestryResults$pheno'_compressed_theor.zip' $crossAncestryResults$'EUR_JAP_LD_'$pheno$'_E' $crossAncestryResults$'EUR_EUR_LD_'$pheno$'_E' $crossAncestryResults$'Combined_EUR_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_LD_'$pheno$'_E' $crossAncestryResults$'Combined_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_EUR_LD_'$pheno$'_E'







# BRCA
pheno='EUR_JAP_BRCA'
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'EUR_JAP_LD_'$pheno '1'
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'EUR_EUR_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'Combined_EUR_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'shaPRS_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'Combined_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'shaPRS_EUR_LD_'$pheno '1' # 

# convert to Elena's format 
convertPRS $crossAncestryResults$'EUR_JAP_LD_'$pheno
convertPRS $crossAncestryResults$'EUR_EUR_LD_'$pheno
convertPRS $crossAncestryResults$'Combined_EUR_LD_'$pheno
convertPRS $crossAncestryResults$'shaPRS_LD_'$pheno
convertPRS $crossAncestryResults$'Combined_LD_'$pheno
convertPRS $crossAncestryResults$'shaPRS_EUR_LD_'$pheno

#rm -rf $crossAncestryResults$pheno'_compressed.zip'
#zip $crossAncestryResults$pheno'_compressed.zip' $crossAncestryResults$'EUR_JAP_LD_'$pheno$'_E' $crossAncestryResults$'EUR_EUR_LD_'$pheno$'_E' $crossAncestryResults$'Combined_EUR_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_LD_'$pheno$'_E' $crossAncestryResults$'Combined_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_EUR_LD_'$pheno$'_E'
zip $crossAncestryResults$pheno'_compressed_theor.zip' $crossAncestryResults$'EUR_JAP_LD_'$pheno$'_E' $crossAncestryResults$'EUR_EUR_LD_'$pheno$'_E' $crossAncestryResults$'Combined_EUR_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_LD_'$pheno$'_E' $crossAncestryResults$'Combined_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_EUR_LD_'$pheno$'_E'


# T2D
pheno='EUR_JAP_T2D'
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'EUR_JAP_LD_'$pheno '1'
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'EUR_EUR_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'Combined_EUR_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'shaPRS_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'Combined_LD_'$pheno '1' # 
#evaluatePRS $crossAncestrySumStats$'all_hm3.fam' 'shaPRS_EUR_LD_'$pheno '1' # 

# convert to Elena's format 
convertPRS $crossAncestryResults$'EUR_JAP_LD_'$pheno
convertPRS $crossAncestryResults$'EUR_EUR_LD_'$pheno
convertPRS $crossAncestryResults$'Combined_EUR_LD_'$pheno
convertPRS $crossAncestryResults$'shaPRS_LD_'$pheno
convertPRS $crossAncestryResults$'Combined_LD_'$pheno
convertPRS $crossAncestryResults$'shaPRS_EUR_LD_'$pheno

#rm -rf $crossAncestryResults$pheno'_compressed.zip'
#zip $crossAncestryResults$pheno'_compressed.zip' $crossAncestryResults$'EUR_JAP_LD_'$pheno$'_E' $crossAncestryResults$'EUR_EUR_LD_'$pheno$'_E' $crossAncestryResults$'Combined_EUR_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_LD_'$pheno$'_E' $crossAncestryResults$'Combined_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_EUR_LD_'$pheno$'_E'
zip $crossAncestryResults$pheno'_compressed_theor.zip' $crossAncestryResults$'EUR_JAP_LD_'$pheno$'_E' $crossAncestryResults$'EUR_EUR_LD_'$pheno$'_E' $crossAncestryResults$'Combined_EUR_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_LD_'$pheno$'_E' $crossAncestryResults$'Combined_LD_'$pheno$'_E' $crossAncestryResults$'shaPRS_EUR_LD_'$pheno$'_E'






#############################################################################

famFile=$crossAncestrySumStats$'all_hm3.fam'

PRSName='EUR_EUR_LD_EUR_JAP_asthma'
calcAUC3='1'

PRSName="RapidoPGS_EUR_JAP_asthma"
PRSName="RapidoPGS_EUR_JAP_asthma_combined"


PRSName="RapidoPGS_EUR_JAP_asthma_LD_EUR"
PRSName="RapidoPGS_EUR_JAP_asthma_LD"
PRSName="RapidoPGS_EUR_JAP_asthma_combined_LD_EUR"


$crossAncestryResults$'shaPRS_LD_'$pheno


$crossAncestryResults$PRSName$'.profile'


function evaluatePRS {
famFile=$1
PRSName=$2
calcAUC3=$3
# call PLINK's score to build PRS for each

awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $crossAncestryResults$PRSName > $crossAncestryResults$PRSName$'_no_dupes'

arguments=' --memory '$plinkMem$' --allow-no-sex --bfile '$crossAncestrySumStats$'all_hm3 --fam '$famFile$' --score '$crossAncestryResults$PRSName$'_no_dupes sum --out '$crossAncestryResults$PRSName
$plink $arguments

# find out accuracy 
awk 'FNR == NR { file1[ $1 ] = $6; next; } FNR <= NR { { if ( $2 in file1 ) { print $2","$6","file1[$2]  } } }' OFS=',' FS=" " $crossAncestryResults$PRSName$'.profile' FS=" " $famFile  > $outste$crossAncestryResults$PRSName$'.csv'
#awk '{count[$6]++} END {for (word in count) print word, count[word]}' $famFile
#awk '{count[$2]++} END {for (word in count) print word, count[word]}' FS="," $outste$crossAncestryResults$PRSName$'.csv'

arguments='/nfs/users/nfs_m/mk23/scripts/correlator.R '$outste$crossAncestryResults$PRSName$'.csv '$crossAncestryResults$PRSName$'_result '$calcAUC3
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments 
}


# generates an output in Elena's format, with a postfix of '_E' added to the input
function convertPRS {
PRSFILE=$1

# remove dupes
awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $PRSFILE > $PRSFILE$'_no_dupes'

# bim file signature
# 1         2            3           4    5       6
#CHR       RSID          X         POS   A1      A2 
#1       rs3131969       0       754182  A       G


# PLINK PRS signature: 
# rs3131969       G       0.000306261375221792

# create file with signature: CHR POS SNP A1 A2 WEIGHT
awk 'FNR == NR { A1[ $1 ] = $2; beta[ $1 ] = $3; next; } 
FNR <= NR { 
{ 
if (FNR == 1) {print "CHR\tPOS\tSNP\tA1\tA2\tWEIGHT"}
if ( $2 in A1 ) { 
if ( $5 == A1[$2] ) { effectAllele= $5; refAllele= $6}
else  { effectAllele= $6; refAllele= $5}
print $1"\t"$4"\t"$2"\t"effectAllele"\t"refAllele"\t"beta[$2]


} } }' $PRSFILE$'_no_dupes' $crossAncestrySumStats$'all_hm3.bim' > $PRSFILE$'_E'

#head $PRSFILE$'_E'

}




PRSFILE=$crossAncestryResults$'EUR_JAP_LD_'$pheno
 


'plink --allow-no-sex --bfile '<UKBB_PLINK_FILE>'  --score '<PRS FILE>' --out '<OUT LOCATION>
$plink $arguments




############################################################################


function convertSumstatsToPLINK1 {
sumData=$1
isBinary=$2
# convert my _sumstats format:
# chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N
# to PLINK:
# CHR     SNP     BP      NMISS   BETA   SE   R2   T   P
#                          X                  X    X                (these are not needed)
# recreate the PLINK .assoc formatted summary stats,
# as pheno is binary, but the summary data was given as Log(OR), and the input for both rapido/LDpred2 is expected to be OR (and SE of LOR)
# I need to exp the LOR to be OR 

if [[ "$isBinary" == '1' ]]; then
echo 'binary'
awk 'FNR <= NR { if (FNR == 1) {print "CHR\tSNP\tBP\tNMISS\tOR\tSE\tR2\tT\tP" } else { print $1"\t"$3"\t"$2"\tX\t"exp($7)"\t"$8"\tX\tX\t"$9 }  }
' $sumData > $sumData$'_plinkformat_temp'
else
echo 'continuous'
awk 'FNR <= NR { if (FNR == 1) {print "CHR\tSNP\tBP\tNMISS\tBETA\tSE\tR2\tT\tP" } else { print $1"\t"$3"\t"$2"\tX\t"$7"\t"$8"\tX\tX\t"$9 }  }
' $sumData > $sumData$'_plinkformat_temp'
fi

# also remove any duplicates
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $sumData$'_plinkformat_temp' > $sumData$'_plinkformat'


rm -rf $sumData$'_plinkformat_temp'
}






pheno='EUR_JAP_BRCA'

pheASumstats=$eursumstatsLoc$'EUR_BRCA_hm3' 
pheBSumstats=$japsumstatsLoc$'JP_BRCA_hm3' 
outPutLocation=$crossAncestryResults$pheno 
rho='0'



pheno='EUR_JAP_T2D'
#shaPRS_new $eursumstatsLoc$'EUR_BRCA_hm3' $japsumstatsLoc$'JP_BRCA_hm3' $crossAncestryResults$pheno '0' $popALDpanel $popBLDpanel


subphenoLoc =  "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/sumstats/eur/EUR_asthma_hm3"
subpheno_otherLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/sumstats/jap/JP_asthma_hm3"
blendFactorLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_asthma_lFDR_meta_SNP_lFDR"
rho =0
outputLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_asthma_sumstats_meta"


pheno='EUR_JAP_BRCA'
shaPRS_new $eursumstatsLoc$'EUR_BRCA_hm3' $japsumstatsLoc$'JP_BRCA_hm3' $crossAncestryResults$pheno '0' $popALDpanel $popBLDpanel




outPutLocation=$crossAncestryResults$'EUR_JAP_T2D'
pheASumstats=$eursumstatsLoc$'EUR_T2D_hm3'
pheBSumstats=$japsumstatsLoc$'JP_T2D_hm3'
rho='0'

outPutLocation=$crossAncestryResults$'EUR_JAP_asthma'
pheASumstats=$eursumstatsLoc$'EUR_asthma_hm3'
pheBSumstats=$japsumstatsLoc$'JP_asthma_hm3'
rho='0'

# Performs the shaPRS step that produces a _sumstats formatted file for both the blended and the combined sumstats, together with a PRS specific LD matrix under $outPutLocation
function shaPRS_new { 
pheASumstats=$1
pheBSumstats=$2
outPutLocation=$3
rho=$4
popALDpanel=$5
popBLDpanel=$6

# Create input file for adjusting script:


awk 'FNR == NR { file1[ $3 ] = $7"\t"$8"\t"$4"\t"$5; next; } 
FNR <= NR { { 
if (FNR == 1) {print "SNP\tCHR\tBP\tBeta_A\tSE_A\tA1.x\tA2.x\tBeta_B\tSE_B\tA1.y\tA2.y" } else {
if ( $3 in file1) { print $3"\t"$1"\t"$2"\t"file1[$3]"\t"$7"\t"$8"\t"$4"\t"$5} }   }
}
' $pheBSumstats $pheASumstats > $outPutLocation$'_SE_meta'



# Export out the lFDR values
arguments='/nfs/users/nfs_m/mk23/scripts/shaPRS_adjust_wrapper.R '$outPutLocation$'_SE_meta '$outPutLocation$'_lFDR_meta /nfs/users/nfs_m/mk23/scripts/shaPRS.R'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

# R script to blend between the two association stats
arguments='/nfs/users/nfs_m/mk23/scripts/sumstatsBlender_shaPRS_meta_wrapper.R '$pheASumstats$' '$pheBSumstats$' '$outPutLocation$'_lFDR_meta_SNP_lFDR '$rho$' '$outPutLocation$'_sumstats_meta /nfs/users/nfs_m/mk23/scripts/shaPRS.R'
/nfs/team152/mk23/software/R_X_farm/R-3.6.1/bin/Rscript $arguments

#head -n 11 $outPutLocation$'_SE_meta' > $asthmaCrossAncestryScratch$'_badData'

if [ "$popALDpanel" = "0" ]; then
echo "NO LD PANEL PROVIDED"
else

# create PRS specific LD matrix
arguments='/nfs/users/nfs_m/mk23/scripts/LDRefGen_wrapper.R '$popALDpanel$' '$popBLDpanel$' '$outPutLocation$'_lFDR_meta_SNP_lFDR '$outPutLocation$'_SE_meta '$outPutLocation$' 0 /nfs/users/nfs_m/mk23/scripts/shaPRS.R'
#/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments
#bsub -q normal -m "modern_hardware" -G team152 -n7 -J LDrefgen -R "span[hosts=1] select[mem>${largePlinkMem}] rusage[mem=${largePlinkMem}]" -M${largePlinkMem} -o $outPutLocation$'.out' -e $outPutLocation$'.err'  "/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/Rscript $arguments"

fi

}
export -f shaPRS_new # this makes local functions executable when bsubbed



Pop1LDRefLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/"
Pop2LDRefLoc ="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/JPTRef/"
blendFactorLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_BRCA_lFDR_meta_SNP_lFDR"
sumstatsLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_BRCA_SE_meta"
outputLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_BRCA"
produceBasicLDreftoo = T
shaPRSscriptLoc="/nfs/users/nfs_m/mk23/scripts/shaPRS.R"


#################################################


sumstatsLoc=$eursumstatsLoc$'EUR_asthma_hm3'
MAFLoc=$crossAncestrySumStats$'EUR_MAF_hm3.frq'

# adds MAF to a sumstats file from a .frq file
function addMAF { 
sumstatsLoc=$1
MAFLoc=$2

#MAF file:
#  1       2       3       4       5          6 
# CHR     SNP     A1      A2      MAF     NCHROBS

# Sumstats file:
#  1       2       3       4       5          6           7       8       9       10
# chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N

awk 'FNR == NR { MAF[ $2 ] = $5; A1[ $2 ] = $3; next; }
FNR <= NR { 
if (FNR == 1) {print $0} else {
if( $3 in MAF ) { 
if( A1[$3] == $4 ) { AF = MAF[$3]}
else {AF = 1-MAF[$3]}
$6 = AF
print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10
}
}} 
' $MAFLoc $sumstatsLoc > $sumstatsLoc$'_MAF'
#head $sumstatsLoc$'_MAF'

}
export -f addMAF # this makes local functions executable when bsubbed






# find out the Genetic correlation between the EUR and JAP summary data for the same traits
ldscDir='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/0_software/ldsc_new/'
ldscDirCondaEnv='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/0_software/ldsc_new/ldscDirCondaEnv/'
lsdc=$ldscDir$'ldsc/ldsc.py'
munge_sumstats=$ldscDir$'ldsc/munge_sumstats.py'
eur_w_ld_chr='/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/ibd/data/0_common/ldsc/eur_w_ld_chr/'

####################################################
# needed a fresh install of LDSC as the old one was buggy
cd $ldscDir


git clone https://github.com/bulik/ldsc.git
cd ldsc
# need to activate the python2 miniconda to install it
# https://docs.hpc.cam.ac.uk/hpc/software-tools/python.html#using-anaconda-python
#module load miniconda/2
conda env  create --prefix $ldscDirCondaEnv --file environment.yml
#  --a1 A1               Name of A1 column , this is described as the 
# "A1: Allele 1, interpreted as ref allele for signed sumstat.", the 'ref' here means that LDSC thinks its the effect allele 

# needed to install pandas older version otherwise would get: pandas ImportError "cannot import name AbstractMethodError"
pip install pandas==0.19.2
###################################

# the LDSC process requires MAF in the summary data, so we add them here
addMAF $eursumstatsLoc$'EUR_T2D_hm3' $crossAncestrySumStats$'EUR_MAF_hm3.frq'
addMAF $eursumstatsLoc$'EUR_height_hm3' $crossAncestrySumStats$'EUR_MAF_hm3.frq'
addMAF $eursumstatsLoc$'EUR_BRCA_hm3' $crossAncestrySumStats$'EUR_MAF_hm3.frq'
addMAF $eursumstatsLoc$'EUR_CAD_hm3' $crossAncestrySumStats$'EUR_MAF_hm3.frq'
addMAF $eursumstatsLoc$'EUR_asthma_hm3' $crossAncestrySumStats$'EUR_MAF_hm3.frq'


addMAF $japsumstatsLoc$'JP_T2D_hm3' $crossAncestrySumStats$'JP_EAF.frq'
addMAF $japsumstatsLoc$'JP_height_hm3' $crossAncestrySumStats$'JP_EAF.frq'
addMAF $japsumstatsLoc$'JP_BRCA_hm3' $crossAncestrySumStats$'JP_EAF.frq'
addMAF $japsumstatsLoc$'JP_CAD_hm3' $crossAncestrySumStats$'JP_EAF.frq'
addMAF $japsumstatsLoc$'JP_asthma_hm3' $crossAncestrySumStats$'JP_EAF.frq'


source activate $ldscDirCondaEnv

#cp $japsumstatsLoc$'JP_T2D_hm3_MAF' $japsumstatsLoc$'JP_T2D_hm3_MAF_old'

###################
# Step 1: merge sumstats

$munge_sumstats --sumstats $eursumstatsLoc$'EUR_T2D_hm3_MAF' --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $eur_w_ld_chr$'w_hm3.snplist'  --out $eursumstatsLoc$'EUR_T2D_hm3.sumstats'
$munge_sumstats --sumstats $japsumstatsLoc$'JP_T2D_hm3_MAF' --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $eur_w_ld_chr$'w_hm3.snplist'  --out $japsumstatsLoc$'JP_T2D_hm3.sumstats'

$munge_sumstats --sumstats $eursumstatsLoc$'EUR_height_hm3_MAF' --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $eur_w_ld_chr$'w_hm3.snplist'  --out $eursumstatsLoc$'EUR_height_hm3.sumstats'
$munge_sumstats --sumstats $japsumstatsLoc$'JP_height_hm3_MAF' --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $eur_w_ld_chr$'w_hm3.snplist'  --out $japsumstatsLoc$'JP_height_hm3.sumstats'

$munge_sumstats --sumstats $eursumstatsLoc$'EUR_BRCA_hm3_MAF' --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $eur_w_ld_chr$'w_hm3.snplist'  --out $eursumstatsLoc$'EUR_BRCA_hm3.sumstats'
$munge_sumstats --sumstats $japsumstatsLoc$'JP_BRCA_hm3_MAF' --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $eur_w_ld_chr$'w_hm3.snplist'  --out $japsumstatsLoc$'JP_BRCA_hm3.sumstats'

$munge_sumstats --sumstats $eursumstatsLoc$'EUR_CAD_hm3_MAF' --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $eur_w_ld_chr$'w_hm3.snplist'  --out $eursumstatsLoc$'EUR_CAD_hm3.sumstats'
$munge_sumstats --sumstats $japsumstatsLoc$'JP_CAD_hm3_MAF' --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $eur_w_ld_chr$'w_hm3.snplist'  --out $japsumstatsLoc$'JP_CAD_hm3.sumstats'

$munge_sumstats --sumstats $eursumstatsLoc$'EUR_asthma_hm3_MAF' --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $eur_w_ld_chr$'w_hm3.snplist'  --out $eursumstatsLoc$'EUR_asthma_hm3.sumstats'
$munge_sumstats --sumstats $japsumstatsLoc$'JP_asthma_hm3_MAF' --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $eur_w_ld_chr$'w_hm3.snplist'  --out $japsumstatsLoc$'JP_asthma_hm3.sumstats'

# Step 2: get the rGs

$lsdc --rg $eursumstatsLoc$'EUR_T2D_hm3'.sumstats.sumstats.gz,$japsumstatsLoc$'JP_T2D_hm3'.sumstats.sumstats.gz --ref-ld-chr $eur_w_ld_chr --w-ld-chr $eur_w_ld_chr --out $crossAncestryResults$'/rG_T2D'

#Genetic Correlation: 0.7009 (0.0438)
#Z-score: 15.9962
#P: 1.3579e-57



$lsdc --rg $eursumstatsLoc$'EUR_height_hm3'.sumstats.sumstats.gz,$japsumstatsLoc$'JP_height_hm3'.sumstats.sumstats.gz --ref-ld-chr $eur_w_ld_chr --w-ld-chr $eur_w_ld_chr --out $crossAncestryResults$'/rG_height'
#Genetic Correlation: 0.682 (0.0184)
#Z-score: 37.0571
#P: 1.3811e-300

$lsdc --rg $eursumstatsLoc$'EUR_BRCA_hm3'.sumstats.sumstats.gz,$japsumstatsLoc$'JP_BRCA_hm3'.sumstats.sumstats.gz --ref-ld-chr $eur_w_ld_chr --w-ld-chr $eur_w_ld_chr --out $crossAncestryResults$'/rG_BRCA'
#Genetic Correlation: 0.6634 (0.1512)
#Z-score: 4.3885
#P: 1.1411e-05


$lsdc --rg $eursumstatsLoc$'EUR_CAD_hm3'.sumstats.sumstats.gz,$japsumstatsLoc$'JP_CAD_hm3'.sumstats.sumstats.gz --ref-ld-chr $eur_w_ld_chr --w-ld-chr $eur_w_ld_chr --out $crossAncestryResults$'/rG_CAD'
#Genetic Correlation: 0.7144 (0.0414)
#Z-score: 17.2734
#P: 7.4582e-67


$lsdc --rg $eursumstatsLoc$'EUR_asthma_hm3'.sumstats.sumstats.gz,$japsumstatsLoc$'JP_asthma_hm3'.sumstats.sumstats.gz --ref-ld-chr $eur_w_ld_chr --w-ld-chr $eur_w_ld_chr --out $crossAncestryResults$'/rG_asthma'
#Genetic Correlation: 0.53 (0.0993)
#Z-score: 5.339
#P: 9.3472e-08







#################################################













[1] "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/"                                               
[2]
                                                                                                          
                                
                                                                                 
                             
                                                                         
 

baseLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/"  
sumstatsLoc= "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_height_sumstats_meta_combined"
NCORES= 15                                                                                      
n_eff=244864                                                                                   
outLoc= "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/"                     
outName="Combined_EUR_LD_EUR_JAP_height"                                                                          
testSetLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/sumstats/all"                 
shaPRSscriptLoc= "/nfs/users/nfs_m/mk23/scripts/shaPRS.R"                                                            
binary_outcome=T       


baseLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_BRCA/"        
sumstatsLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_BRCA_sumstats_meta_combined"
NCORES= 15                                                                                      
n_eff=127781                                                                                   
outLoc= "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/"                     
outName="Combined_LD_EUR_JAP_BRCA"                                                                          
testSetLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/sumstats/all"                 
shaPRSscriptLoc= "/nfs/users/nfs_m/mk23/scripts/shaPRS.R"                                                            
binary_outcome=T                                                        




subphenoLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/sumstats/eur/EUR_asthma_hm3"
subpheno_otherLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/sumstats/jap/JP_asthma_hm3"
blendFactorLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_asthma_lFDR_meta_SNP_lFDR"
rho = 0
outputLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_asthma_sumstats_meta"
shaPRSscriptLoc= "/nfs/users/nfs_m/mk23/scripts/shaPRS.R"   



	
	
Pop1LDRefLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/"       
Pop2LDRefLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/JPTRef/"        
blendFactorLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_asthma_lFDR_meta_SNP_lFDR"
sumstatsLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_asthma_SE_meta"
outputLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_asthma2"
produceBasicLDreftoo = T                                                                     
shaPRSscriptLoc="/nfs/users/nfs_m/mk23/scripts/shaPRS.R"

typeof(inputData$SE_A)

as.numeric(inputData$SE_A[1])
nrow(inputData)


non_nomerics <- which(is.na(as.numeric(as.character(inputData$SE_A))))
nonNumericInputDatas = inputData[non_nomerics,]

inputDataLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_BRCA_SE_meta"
inputDataLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/scratch/_badData"
outputLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/results/EUR_JAP_BRCA_lFDR_meta"
shaPRSscriptLoc = "/nfs/users/nfs_m/mk23/scripts/shaPRS.R"


inputData_subset=inputData[1:5,]

inputData_subset$SE_A * inputData_subset$SE_B

as.numeric(inputData_subset$SE_A) * as.numeric(inputData_subset$SE_B)

nonNumericA_incides = which(!is.numeric(inputData$SE_A) )

nonNumericB_incides = which(!is.numeric(inputData$SE_B) )


nonNumericinputData = inputData[which(!is.numeric(inputData$SE_A) ),]

inputData[which(!is.numeric(inputData$SE_B) ),]


SNP        CHR     BP    Beta_A    SE_A          A_effectAllele Beta_B   SE_B
rs1048488   1 7   60912  0.11737  0.8249            T           0.061    0.031

nonNumericinputData$SE_A

as.numeric(nonNumericinputData$SE_A) * as.numeric(nonNumericinputData$SE_B)

inputData_subset$SE_A = as.numeric(as.character(inputData_subset$SE_A))

as.numeric(as.character(nonNumericinputData$SE_A)) * as.numeric(as.character(nonNumericinputData$SE_B))

Pop1LDRefLoc = args[1]
Pop2LDRefLoc = args[2]
blendFactorLoc = args[3]
sumstatsLoc = args[4]
outputLoc = args[5]
produceBasicLDreftoo = (args[6] == '1')
shaPRSscriptLoc= args[7]


Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 102

##############################################################
# 1) Generate JP RDMS LD panel for hapmap3
###########################
# install R 4, as Rapido now requires R > 4.0
cd /nfs/team152/mk23/software/R_X_farm
wget https://cran.r-project.org/src/base/R-4/R-4.1.0.tar.gz

tar zxvf R-4.1.0.tar.gz
cd R-4.1.0/
./configure --prefix=/nfs/team152/mk23/software/R_X_farm/R-4.1.0/
make
###########################


/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/R

cd /lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/JPTRef/

/nfs/team152/mk23/software/R_X_farm/R-4.1.0/bin/R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")
#install.packages("RapidoPGS")
#install.packages("remotes")a

library(remotes)
install_github('GRealesM/RapidoPGS')


library("RapidoPGS")
reference="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/1KG/"

create_1000G_fixed(directory = reference, remove.related=TRUE, qc.maf = 0.00001, qc.hwe=1e-50, qc.geno=1, autosomes.only=TRUE)

create_1000G(directory = reference, remove.related=TRUE, qc.maf = 0.00001, qc.hwe=1e-50, qc.geno=1, autosomes.only=TRUE)


library("utils")
install.packages("R.utils")
library("R.utils")
importFrom(utils,download.file)

@importFrom

@importFrom utils download.file setTxtProgressBar txtProgressBar

utils::download.file

library("bigsnpr")


		  
create_1000G_fixed <- function(directory = "ref-data", remove.related=TRUE, qc.maf = 0.01, qc.hwe=1e-10, qc.geno=0, autosomes.only=TRUE){
  
  # Remove annoying timeout limit in download.file
  timeout <- getOption('timeout')
  options(timeout=10000)
  
  # Max chromosomes to be considered, if autosomes.only is set to FALSE it will
  # be modified below
  max.chr=22
  
  dir.create(directory)
  plink <- download_plink(directory) #  "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/1KG/plink" # 
  message("Downloading 1000 Genomes Phase III files in VCF format.")
  message("These are (very) big files, so they'll take a while to download. Time for a cuppa.")
  for(chr in c(1:max.chr)){
    message("Downloading chromosome ",chr,"...")
    download.file(paste0("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"),
                  destfile = paste0(directory,"/chr",chr,".vcf.gz"), mode = "wb")
    message("Done!")
  }

  message("Transforming files into Plink format...")
  for(chr in 1:max.chr){
    system(paste0(plink, " --vcf ",directory, "/chr", chr, ".vcf.gz --make-bed --out ", directory, "/chr", chr), ignore.stdout = FALSE, ignore.stderr = FALSE)
  }
  # We don't need our hard-earned vcf files anymore, so we can delete them
  unlink(paste0(directory,"/*vcf.gz"))
  message("Done!")
  message("Starting QC procedure...")
  ped <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20200731.ALL.ped")
  pop <- bigreadr::fread2("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv")
  
  for(chr in 1:max.chr){
    chrno <- ifelse(chr <23, chr, ifelse(chr==23, "X", "Y"))
    message("Processing Chr ", chrno, "...")
    bed <- snp_plinkQC(plink, prefix.in = paste0(directory,"/chr",chr),
                       prefix.out = paste0(directory,"/chr",chr,"_QC"),
                       geno = qc.geno, maf = qc.maf, hwe = qc.hwe)
    rds <- snp_readBed(bed)
    snp <- snp_attach(paste0(directory,"/chr",chr,"_QC.rds"))  
    G <- snp$genotypes
    fam <- snp$fam
    fam <- dplyr::left_join(fam, ped, by = c("sample.ID" = "Individual ID"))
    fam <- dplyr::left_join(fam, pop, by = c("Population" = "Population Code"))
    
    snp$fam$family.ID <- paste(fam$`Super Population`, fam$Population, sep = "_")
    snp$fam$paternal.ID <- snp$fam$maternal.ID <- 0L
    if(remove.related){
      unrelated_tags <- c("unrel","unrels","father","mother")
      ind_norel <- which(fam$Relationship %in% unrelated_tags & fam$Siblings == 0  & fam$`Second Order` == 0 & fam$`Third Order` == 0)
      
      maf <- snp_MAF(G, ind.row = ind_norel)
      
      bed <- snp_writeBed(snp, tempfile(fileext = ".bed"),
                          ind.row = ind_norel, ind.col = which(maf > 0.05))
    } else{
      
      maf <- snp_MAF(G)
      
      bed <- snp_writeBed(snp, paste0(directory, "/chr",chr,"_tmp.bed"),
                          ind.col = which(maf > qc.maf))
    }
    rds <- snp_readBed(bed)
    snp <- snp_attach(rds)
    unlink(paste0(directory, "/chr", chr,"\\.*"))
    snp_writeBed(snp, paste0(directory, "/chr",chr,".bed"))
    message("Done!")
  }
  unlink(paste0(directory, "/*_QC.*"))
  unlink(paste0(directory, "/*_tmp.*"))
  unlink(paste0(directory, "/plink"))
  message("Done! Now you can find your reference panel ready at ", directory,"/.")
  on.exit(options(timeout=timeout))
}




###########################################a


library("bigsnpr")
library("R.utils")

reference="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/1KG/"
mapLocation="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/map.rds"
outputLoc="/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/JPTRef2/"
chrs=21
ncores=16


#chr_map_rds_orig = readRDS(file = paste0("/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/JPTRef/","LD_chr",21,"_map.rds") )

# 1. load Map data
map_rds = readRDS(file = mapLocation )


for(chrs in 1:22){
# 2. grab the RSids from the map, that pertain to the given chrom:
chrom_SNPs = paste(map_rds$chr[ which(map_rds$chr == chrs)], map_rds$pos[ which(map_rds$chr == chrs)], sep=":" )

# load Bed file
 if(!file.exists(paste0(reference,"chr",chrs,".rds"))){
        snp_readBed(paste0(reference,"chr",chrs,".bed"))
      }
	  
obj.bigSNP <- snp_attach(paste(reference,"chr",chrs, ".rds", sep=""))


 # find the indices of the SNPs in BigSNP that are also on the LD map
SNP_chrPos = paste(obj.bigSNP$map$chromosome,obj.bigSNP$map$physical.pos, sep=":")
snp.idx <- which(SNP_chrPos %in% chrom_SNPs) 

# find duplicates in the surviving list
SNPIDs = SNP_chrPos[snp.idx] # subset orig bigSNP list to the list we have found
# find the indices of the duplicates
dupe_indices= which(duplicated(SNPIDs) )
if( length(dupe_indices) > 0 ) {
snp.idx = snp.idx[-dupe_indices]
}


# In this case, we'll focus only on individuals of EAS ancestry
euridx  <- grep("EAS", obj.bigSNP$fam$family.ID)
# Get aliases for useful slots
G   <- obj.bigSNP$genotypes  #  2299 109774

# need to create genetic positions from the bp positions, as otherwise the ' size = 3 / 1000'  won't work
POS2 <- snp_asGeneticPos( obj.bigSNP$map$chromosome, infos.pos = obj.bigSNP$map$physical.pos, dir = outputLoc, ncores = ncores)
# 109774
# https://figshare.com/articles/dataset/European_LD_reference/13034123?file=25503857 use 3kb for size as per LDpred 2
LDmat <- snp_cor(G, ind.col = snp.idx, ind.row= euridx, ncores = ncores, size = 3 / 1000, infos.pos = POS2[snp.idx] )


fileLoc= paste0(outputLoc,"LD_chr",chrs,".rds")
saveRDS(LDmat,file = fileLoc)
print(paste0("written LD mat to ",fileLoc ))


# map the final list of SNPs back to the original map file's indices
# so that the map_rds file's indices match the indices in the LDmat 
# (obj_bigSNP_pos is always a SUBSET of map_rds_new_pos)




map_rds_new = map_rds[which(map_rds$chr == chrs),]

matchingIndices = match(paste(obj.bigSNP$map$chromosome[snp.idx],obj.bigSNP$map$physical.pos[snp.idx], sep=":"), paste(map_rds_new$chr,map_rds_new$pos, sep=":"))

map_rds_subset = map_rds_new[matchingIndices,]  # match the first to the second

# overwrite alleles (take PLINK's A1 as a1 here )
map_rds_subset$a1 = obj.bigSNP$map$allele1[snp.idx]
map_rds_subset$a0 = obj.bigSNP$map$allele2[snp.idx]

# overwrite the MAF
map_rds_subset$af_UKBB = snp_MAF(G, ind.col = snp.idx, ind.row= euridx, ncores = ncores)
head(obj.bigSNP$map$allele1[snp.idx])

# overwrite current LD
map_rds_subset$ld = Matrix::colSums(LDmat^2) 

# overwrite alleles

fileLoc= paste0(outputLoc,"LD_chr",chrs,"_map.rds")
saveRDS(map_rds_subset,file = fileLoc)
print(paste0("written chr map to ",fileLoc ))
}


# concat the map from each chroms 
all_map_rds = NULL
for(chrs in 1:22){
chr_map_rds = readRDS(file = paste0(outputLoc,"LD_chr",chrs,"_map.rds") )
all_map_rds = rbind(all_map_rds,chr_map_rds)
}

fileLoc= paste0(outputLoc,"map.rds")
saveRDS(all_map_rds,file = fileLoc)
print(paste0("written overall map to ",fileLoc ))




# 4. execute on console generating chr21 PRS LDRefGen ( asthma)
# 5. loop and submit all 22 chroms as jobs
# 6. script to concat/merge a map.rds from the 22 blended LDref
# 7. execute RapidoMulti

# 2) implement the final version of the LD blending scripts/ a
# this should output an LD ref panel for each run, which will be used 
# also prepare baselines:
# original (single with no LDref)a
# LD_EUR(1-lfDR) + LD_JP(LFDR)


# 3) re-run Asthma with hapmap3 with RapidoPGS-multi a


echo $crossAncestrySumStats$'old/_lFDR_meta_SNP_lFDR_hapmap3_WITH_HLA'
echo $crossAncestrySumStats$'old/_SE_meta_hapmap3_WITH_HLA'





library("bigsnpr")

Pop1LDRefLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/"
Pop2LDRefLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/shaprs/raw/JPTRef/"
chromNum = 21
blendFactorLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/sumstats/old/_lFDR_meta_SNP_lFDR_hapmap3_WITH_HLA"
sumstatsLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/sumstats/old/_SE_meta_hapmap3_WITH_HLA"
outputLoc = "/lustre/scratch123/hgi/mdt1/projects/crohns/mk23/0Thesis/asthma/crossAncestry/sumstats/old/asthmaOut"
produceBasicLDreftoo = T
shaPRSscriptLoc = '/nfs/users/nfs_m/mk23/scripts/shaPRS.R'

chromNum=21







