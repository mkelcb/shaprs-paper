library(data.table)

## subset samples as in ld2pred paper https://github.com/privefl/simus-PRS/blob/master/paper3-SCT/code_real/01-SCT-BRCA.R


## read snakemake input and output files:

score.f <- snakemake@input[['score']]
rel_ind.f <- snakemake@input[['rel_ind']]
pheno.f <- snakemake@input[['eid']]
out.f <- snakemake@output[['score']]

score <- fread(score.f)

## get f.eid , sex, batch  and ethniticy (codes "22001-0.0", "22000-0.0" and "21000-0.0")

dt0 <- fread(pheno.f, header=T, select=c("eid", "22001-0.0", "22000-0.0" , "21000-0.0"),
             col.names = c("f.eid", "sex", "batch", "ethnicity"))

## add related individual info to dt0

rel_ind <- dt0$f.eid %in% fread(rel_ind.f)$ID2

dt0[ , is_rel2 := rel_ind]

## get cancer disease fields
df_cancer0 <- fread(pheno.f, select = paste0("2453-", 0:2, ".0"))
df_cancer1 <- fread(pheno.f, select = c(paste0("20001-0.", 0:5),
                                     paste0("20001-1.", 0:5),
                                     paste0("20001-2.", 0:5)))
df_cancer2 <- fread(pheno.f, select = c(paste0("40001-", 0:2, ".0"),
                                     paste0("40002-0.", 0:13),
                                     paste0("40002-1.", 0:13),
                                     paste0("40002-2.", 0:13),
                                     paste0("40006-", 0:31, ".0"),
                                     paste0("41202-0.", 0:379),
                                     paste0("41204-0.", 0:434)))
## select BRCA
ind_BRCA <- sort(unique(unlist(c(
  lapply(df_cancer1,  function(x) which(x == 1002)),
  lapply(df_cancer2, function(x) which(substr(x, 1, 3) %in% c("C50", "D05")))
))))



## select cases and controls
y <- rep(NA, nrow(dt0))

y[rowMeans(df_cancer0 == 0, na.rm = TRUE) == 1] <- 0 # no cancer

y[ind_BRCA] <- 1

## add phenotype to dt0
dt0[, BRCA:=y]


## merge scores with dt0

score <- merge(dt0, score, by="f.eid")

## select scores for un-related, white british, females and exclude pilot study


score <- score[ethnicity ==1001 & !is_rel2 & !is.na(BRCA) & sex==0 & batch > 0, ]


write.table(score, out.f, row.names=F)
