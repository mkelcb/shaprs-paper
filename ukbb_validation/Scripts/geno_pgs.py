import pandas as pd
import os
import re
import numpy as np
import sys
import pathlib


def snps_qc(pgs_f, ukbb_snps, info, out, chrom_col="CHR", pos_col= "BP", ref_col="REF", alt_col="ALT", weight_col="weight"):
  """ QC SNPs based on info score and matching alleles
Parameters:
  # pgs_f (char): file with pgs score
  # chrom_col (char):column name for chromosome
  # pos_col (char):column name for position built
  # ref_col (char):column name for reference allele
  # alt_col (char):column name for alternative allele
  # weight_col (char): column name for weight for alternative allele
  # ukbb_snps (char):full name for files with ukbb snp info per chromosome, each file lists the minor allele frequency and info score for each of the markers in the imputed data, calculated using QCTOOL. The order of markers in these files is not guaranteed to be the same as the BGEN files.
  Cols:
  Alternate_id
  RS_id
  Position
  Allele1
  Allele2
  Minor Allele
  MAF
  Info score
  # info (int): numeric with info score threshold
  # out (char):directory to save output files
# Returns:
  saves 2 files: first with failed SNPs (if any) and second one pgs file with SNPs to use in UKBB.

"""
  df = pd.read_csv(pgs_f, sep = '\s+')
   
  # make sure df has no NA in weight column
  na_w = df[weight_col].isnull().values.any()
  if na_w:
    sys.exit("NA values in weight column, correct and start again")
   
  regex = re.compile(".*chr([0-9]+)_.*")
  keep=None
  ex=None
  for f in ukbb_snps:
    
    ukbb = pd.read_csv(f, sep = '\t', usecols = list(range(5)) + [7],
              names = ["id", "rsid", "pos", "ref_bb", "alt_bb", "info_score"])
    chromosome = int(regex.findall(os.path.basename(f))[0])
    ukbb['chrom'] = chromosome

    # merge by position, then do QC
    df2 = pd.merge(df.loc[df[chrom_col] == chromosome], ukbb, how= 'left', left_on = [pos_col], right_on = ['pos'])
     
    # Get failed SNPs
    failed = df2[(df2['info_score'].isnull() ) | (df2['info_score'] < info )]

    df2 = df2.loc[df2['info_score'] >= info]

    ## look for multiallelic snps based on duplicated positions (ukbb pos) and mark all occurrances as True (keep=False)

    df2['multiallelic']= df2.duplicated(['pos'], keep=False)

    failed = failed.append(df2[df2['multiallelic'] == True], ignore_index=True)

    failed.fillna({'multiallelic' : False}, inplace=True)

    # remove multiallelic from df2

    df2 = df2.loc[df2['multiallelic'] == False]
     
    ## look for allelic mismatches:

    df2['flip'] = np.where( (df2[ref_col] == df2['ref_bb']) & (df2[alt_col] == df2['alt_bb']) , False, True)

    # make sure when flip=True that the alleles in pgs_f and ukbb are the same, exclude otherwise

    failed = failed.append(df2[~(df2['flip'] == True &(df2[ref_col] == df2['alt_bb']) &(df2[alt_col] == df2['ref_bb']))])

    df2 = df2[df2['flip'] == True & (df2[ref_col] == df2['alt_bb']) & (df2[alt_col] == df2['ref_bb'])]

    if keep is None:
      keep = df2.copy()
    else:
      keep = keep.append(df2)
       
    if failed.shape[0] != 0:
      if ex is None:
        ex = failed.copy()
      else:
        ex = ex.append(failed)

  # sort keep and ex by chrom and position before saving
   
  keep.sort_values(by=[chrom_col, pos_col] , inplace=True)
  failed.sort_values(by=[chrom_col, pos_col] , inplace=True)

  keep.to_csv(out + "/pgs_ukbb_" + os.path.basename(pgs_f), index=False, sep=" ")
  failed.to_csv(out + "/failed_" + os.path.basename(pgs_f), index=False, sep= " ") 
     

def pgs_compute(bgen_open, chrom, df, chr_col="CHR", pos_col = "BP", weight_col = "weight", threads=1):
  """ QC SNPs based on info score and matching alleles
  Parameters:
  # bgen_open (open_bgen object): an open_bgen object to facilitate data extraction from bgen file
  # chrom (char) chromosme number for bgen file
  # df (pandas data frame): data frame with pgs score
  # chr_col (char):column name for chromosome in df
  # pos_col (char):column name for position in df
  # weight_col (char): column name for weigth for the alternative allele
  # threads (int): number of threads 
      """
   
  # Split by "flip" column
  df1 = df.loc[df['flip'] == False]
  df2 = df.loc[df['flip'] == True]

  # reset indices
  df1.reset_index(drop=True, inplace=True)
  df2.reset_index(drop=True, inplace=True)

  pos = bgen_open.positions
   
  # Call pgs_comp_help for df1 and df2 (map)

  score_1 = 0
  score_2 = 0

  if df1.shape[0]:
     
    score_1 = pgs_comp_help(bgen_open, df1, chrom, pos, pos_col, weight_col, flip=False, threads=threads)

  if df2.shape[0]:

    score_2 = pgs_comp_help(bgen_open, df2, chrom, pos, pos_col, weight_col, flip=True, threads=threads)
     
  return(score_1 + score_2)



def pgs_comp_help(bgen_open, df1, chrom, pos, pos_col, weight_col, flip=True, threads=1):

  # I need to get the indices of the ids that correspond to the ids in df. Searching is much faster by position (binary serach) rather than linear by ids. (I first get the indices for matching positions and then I check which of the indices corresponds to the ids I want to extract). Because I have already excluded multiallelic SNPs each position I am retriving is unique and I have already excluded SNPs if alleles dont match. Therefore I dont really need to check weather ids match for SNPs on the same position. This makes the code faster.

  # get the indices for matching positions (convert pandas column into np.array first)
  idx_pos = np.searchsorted(pos, df1.loc[:, pos_col].values)


  sys.stderr.write("Extracting variants from chromosome %s\n" % chrom)

  probs = bgen_open.read(idx_pos, dtype = np.float32, num_threads= threads)

  w=df1.loc[:,weight_col].values
   
  if flip:

    res = np.matmul(probs[:, :, 1] , w) + np.matmul(probs[:, :, 0] , 2*w)
     
  else :

    res = np.matmul(probs[:, :, 1] , w ) + np.matmul(probs[:, :, 2] , 2*w)
     
  return(res)


def sum_subscores(score_files, out_file, samples = None):
  """Given files with np.arrays with the score values for a subset variants across samples, add them to calculate the total score per sample. If samples is not None, add the sample id and save output.
  Parameters
  score_files (list): list with path to files to pre-compute scores (txt files)
  out_file (string): full name to output file
  samples (string): full name to sample file for bgen file, defaults to None.
  """
  score = None

  for f in score_files:
    a = pd.read_csv(f)['score']
     
    if score is None:
      score = a.copy()
    else:
      score += a


  # convert score to df

  df = score.to_frame()
  if samples is not None:
    sam = pd.read_csv(samples, sep = ' ', skiprows=[1], usecols=[0])
    df['f.eid'] = sam['ID_1'] 


  df.to_csv(out_file, sep=" ", index = False)
