import sys

sys.path.insert(0, '/home/ev250/myPypacks/pyrunPGS/pyrunPGS')

import geno_pgs


# extract inputs

pgs = snakemake.input['pgs']
bbsnps = snakemake.input['bbsnps']
info = int(snakemake.params['info'])
cols = snakemake.params['cols']
out_dir = snakemake.params['out_dir']

geno_pgs.snps_qc(pgs_f=pgs,
                 chrom_col=cols[0],
                 pos_col =cols[1],
                 ref_col= cols[2],
                 alt_col=cols[3],
                 weight_col=cols[4],
                 ukbb_snps=bbsnps,
                 info=info,
                 out=out_dir)
