# removes nans from an ld matrix

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
#print("Arguments received:")
#print(args)
if(length(args) < 3) {stop("not enough Arguments received")} 
 
ldmatloc = args[1]
snplistloc= args[2]
outLoc = args[3]

# DEBUG vars
# ldmatloc="C:/0Datasets/shaPRS/0revisions/chris/ldblk_1125_1kg.ld"
# snplistloc="C:/0Datasets/shaPRS/0revisions/chris/ldblk_1125_1kg.snplist"
# outLoc="C:/0Datasets/shaPRS/0revisions/chris/ldblk_1125_1kgfixed"

ldmat=read.table(ldmatloc, header  = F)
snplist= read.table(snplistloc, header  = F)
colnames(ldmat) = rownames(ldmat)  = unlist(snplist)


ldmat_noNAs= ldmat[complete.cases(ldmat),complete.cases(ldmat) ]
if(nrow(ldmat_noNAs) < nrow(ldmat)) {
  print(paste0("eliminated ", (nrow(ldmat) - nrow(ldmat_noNAs) ), " nan SNPs" ))
}

write.table(ldmat_noNAs, paste0(outLoc,".ld"), quote = F,col.names = F,row.names = F, sep="\t")
write.table(colnames(ldmat_noNAs), paste0(outLoc,".snplist"), quote = F,col.names = F,row.names = F, sep="\t")

