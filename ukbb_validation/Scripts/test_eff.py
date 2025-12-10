# Improving efficiency for processing UKBB genotype data
# Given a pgs file compute a pgs score in ukbb (all samples) for a range of lines. The end line is given by the breaks wildcard and the starting line will cover at most the range end_line -chunk, see below. Then to reduce memeory usage a for loop split the range 'b' times. The scores are added up sequentially and saved as an np.array at the end.


# system imports
import os
import pandas as pd
import sys
import pathlib
import numpy as np
import psutil


# import local functions
sys.path.insert(0, '/home/ev250/myPypacks/pyrunPGS/pyrunPGS')
import geno_pgs

# import importlib
# importlib.reload(geno_pgs)

# get Snakefile inputs:
bgen = snakemake.input['bgen']
pgs_f = snakemake.input['pgs']
meta_dir = snakemake.params['meta_dir']
chrom = int(snakemake.params['chrom'])
cols = snakemake.params['cols']
end_line=int(snakemake.wildcards['breaks'])
chunk = int(snakemake.params['chunk'])
b= int(snakemake.params['sub_chunk'])
out=snakemake.output['out']
threads=int(snakemake.threads)


# bgen = "/home/ev250/rds/rds-mrc-bsu/biobank/EGAD00010001474/ukb_imp_chr3_v3.bgen"
# pgs_f =  "/home/ev250/rds/rds-cew54-wallace-share/Projects/rapidoPGS/revision/model_QC/pgs_ukbb_T2D_Scott_28566273_1_qcfilt_RapidoPGSsingle_prior02.full.model"
# chrom = 20
# cols=["CHR", "BP",  "weight"]
# meta_dir =  "/home/ev250/rds/rds-mrc-bsu/ev250/rapidoPGS/revision/meta_data_bgen_py"
# chunk = 10**5
# end_line = 142097
# b=2*10**3
# threads = 32

# pgs_f = '/home/ev250/rds/rds-mrc-bsu/ev250/martinPRS/model_QC/pgs_ukbb_EUR_JAP_LD_EUR_JAP_BRCA_E'
# end_line = 48678
# chrom=3
# cols=["CHR", "POS",  "WEIGHT"]

# set environmental variable before calling bgen_reader
# set dir to save metadata
os.environ["BGEN_READER_CACHE_HOME"] = meta_dir



###############################
## Deal with bgen file
###############################

# get dir with bgen file
bgen_dir = os.path.dirname(bgen)

# get basename of bgen file
bgen_name = os.path.basename(bgen)

# change working directory to dir with bgen files (recommended by Carl from  bgen_reader-py) 
os.chdir(bgen_dir)


###############################
## Deal with pgs
###############################

chr_col=cols[0]
pos_col=cols[1]
weight_col=cols[2]

df = pd.read_csv(pgs_f, sep=" ")
df = df.loc[df[chr_col] == chrom]
# reset index
df.reset_index(drop=True, inplace=True)


###############################
## Process data
###############################

from bgen_reader._metafile import infer_metafile_filepath
metadata_file = infer_metafile_filepath(pathlib.Path(bgen_name))
print(metadata_file)

from bgen_reader import open_bgen
bgen_open = open_bgen(bgen_name, verbose=False)

print(bgen_name)
print(bgen_open.ids[:5])
print(bgen_open.read(index=0, max_combinations=None )[0])

start_idx = end_line - chunk if end_line % chunk == 0 else  end_line - (end_line % chunk)
end_idx = end_line -1

# Split variants in groups of b variants and compute scores sequentially 

start_job = [x for x in range(start_idx, end_idx,b)]
end_job = [x-1 for x in range(start_idx + b, end_idx,b)] + [end_idx]


print(start_job)
print(end_job)

# call functions to compute score


score = None

for i in range(len(end_job)):
    st=start_job[i]
    end=end_job[i]
    temp = geno_pgs.pgs_compute(bgen_open=bgen_open, chrom=chrom, df=df[st: end], chr_col=chr_col, pos_col = pos_col, weight_col = weight_col, threads=threads)
    print("iteration number  = " + str(i))
    if score is None:
        score = temp
    else:
        score += temp
   
if score is not None:

    score_df = pd.DataFrame(score, columns=['score'])

    score_df.to_csv(out, sep=" ", index=False)
        
 
del bgen_open

  
