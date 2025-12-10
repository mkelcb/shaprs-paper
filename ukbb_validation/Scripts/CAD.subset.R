library(data.table)

## subset samples as in ld2pred paper https://github.com/privefl/simus-PRS/blob/master/paper3-SCT/code_real/07-SCT-CAD.R


## read snakemake input and output files:

score.f <- snakemake@input[['score']]
rel_ind.f <- snakemake@input[['rel_ind']]
pheno.f <- snakemake@input[['eid']]
out.f <- snakemake@output[['score']]

score <- fread(score.f)

## get f.eid batch and ethnicity (codes "22000-0.0" and "21000-0.0")

dt0 <- fread(pheno.f, header=T, select=c("eid", "22000-0.0", "21000-0.0"),
             col.names = c("f.eid", "batch", "ethnicity"))


## add related individual info to dt0

rel_ind <- dt0$f.eid %in% fread(rel_ind.f)$ID2

dt0[ , is_rel2 := rel_ind]

## get heart disease fields
df_heart <- fread(pheno.f, select = c(paste0("6150-0.", 0:3),
                                   paste0("6150-1.", 0:3),
                                   paste0("6150-2.", 0:3)))
df_illness <- fread(pheno.f, select = c(paste0("20002-0.", 0:28),
                                     paste0("20002-1.", 0:28),
                                     paste0("20002-2.", 0:28)))
df_ICD10 <- fread(pheno.f, select = c(paste0("40001-", 0:2, ".0"),
                                   paste0("40002-0.", 0:13),
                                   paste0("40002-1.", 0:13),
                                   paste0("40002-2.", 0:13),
                                   paste0("41202-0.", 0:379),
                                   paste0("41204-0.", 0:434)))
ind_CAD <- sort(unique(unlist(c(
  lapply(df_heart,   function(x) which(x == 1)),
  lapply(df_illness, function(x) which(x == 1075)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) %in% paste0("I", 21:23))),
  lapply(df_ICD10,   function(x) which(x == "I252"))
))))

ind_heart <- sort(unique(unlist(c(
  lapply(df_heart,   function(x) which(x %in% 1:3)),
  lapply(df_illness, function(x) which(x %in% 1074:1080)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "I"))
))))

y <- rep(0, nrow(dt0)); y[ind_heart] <- NA; y[ind_CAD] <- 1

## add phenotype to dt0
dt0[, CAD:=y]


## merge scores with dt0

score <- merge(dt0, score, by="f.eid")

## select scores for un-related, white british excluding other heart diseases


score <- score[ethnicity ==1001 & !is_rel2 & !is.na(CAD) & batch > 0, ]


write.table(score, out.f, row.names=F)
