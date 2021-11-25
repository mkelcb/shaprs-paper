library(data.table)
library(ggplot2)
library(pROC)


## read snakemake input and output files:
score_f <- snakemake@input[['score']]
snp.ex <- snakemake@input[['excluded']]
snp.kept <- snakemake@input[['kept']]
out_dir <- snakemake@params[['out_dir']]

## score_f <- list.files('/home/ev250/rds/rds-mrc-bsu/ev250/martinPRS/ukbb_score/ind_subset/', full.names=T)
## snp.ex <- list.files('/home/ev250/rds/rds-mrc-bsu/ev250/martinPRS/model_QC',"failed",  full.names=T)

model <- unlist(lapply(score_f, function(i) gsub(".score.txt", "", basename(i))))

################################
## Perform SNP QC per model ##
################################

failed <- lapply(model, function(i) fread(grep(paste("failed", i, sep = "_"), snp.ex, value=T)))
good <- lapply(model, function(i) fread(grep(paste("pgs_ukbb", i, sep="_"), snp.kept, value=T)))

names(failed) <- names(good) <- model

## Count failed SNPs per model

count.failed <- rbindlist(lapply(model, function(i) data.table(model=i,
                                                     different_alleles=failed[[i]][multiallelic == FALSE & info_score >=0.3,.N],
                                                     multiallelic_alleles=length(unique(failed[[i]][multiallelic==TRUE, rsid])),
                                                     low_info=failed[[i]][multiallelic==FALSE & info_score<0.3, .N])))


count.failed[, total_failed := different_alleles +  multiallelic_alleles + low_info]

total.good <- lapply(good, nrow)

count.failed[ , total_okay:= unlist(total.good)][ ,per_failed:= round(total_failed*100/(total_failed + total_okay), 2)]



##################
## Benchmark scores ##
##################

scores.l <- lapply(model, function(i) fread(grep(paste0("/", i), score_f, value=T)))
names(scores.l) <- model

trait <- unlist(lapply(model, function(i) {
    v = strsplit(i, "_")[[1]]
    v <- v[length(v) - 1]
    if(v == "T2D") v <- "TD2"
    return(v)
}))

names(trait) <- model

## Number of individuals run in each model and cases/controls if appropriate
ninds <- rbindlist(mapply(function(i,j) {
    if(j=="T2D") j <- "TD2"
    dt <- i[,.N, get(j)]
    setnames(dt, 'get', 'CaseCtl')
    dt[ , total:=nrow(i)]
    ## if(length(unique(i[[j]])) == 2) dt[, cases:=i[get(j) == 1,.N]][,contols:=i[get(j)==0,.N]]
    return(dt)
}, i=scores.l, j= trait, SIMPLIFY=F),fill=T, idcol="Model")


## AUC and r2 per model

auc <- rbindlist(lapply(names(trait), function(t) data.table(model=t, AUC=roc(get(trait[[t]]) ~ score, scores.l[[t]])$auc[1])))


r2 <- rbindlist(lapply(names(trait), function(t) data.table(model=t, r2=cor(scores.l[[t]][[trait[[t]]]], scores.l[[t]][['score']])^2)))


auc_r2 <- merge(auc, r2,  by='model')

## r <- rbindlist(lapply(names(trait), function(t) data.table(model=t, r=cor(scores.l[[t]][[trait[[t]]]], scores.l[[t]][['score']]))))

##################
## save output ##
##################

write.table(auc_r2, paste0(out_dir, "/auc_r2.txt"), row.names=FALSE)

write.table(ninds, paste0(out_dir, "/N_ukbb.txt"), row.names=F)

write.table(count.failed, paste0(out_dir, "/snp_qc_ukbb.txt"), row.names=F)
