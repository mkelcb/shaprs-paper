library(data.table)
library(ggplot2)
library(pROC)


## read snakemake input and output files:
score_x <- snakemake@input[['score_prscsx']]
score_rest <- snakemake@input[['score_other']]
snp.ex <- snakemake@input[['excluded']]
snp.kept <- snakemake@input[['kept']]
out_dir <- snakemake@params[['out_dir']]
traits <- snakemake@params[['trait']]


## score_rest <- list.files('/home/ev250/rds/rds-mrc-bsu/ev250/martinPRS/ukbb_score/test_subset', full.names=T)

## score_x <- list.files('/home/ev250/rds/rds-mrc-bsu/ev250/martinPRS/ukbb_score/prscsx_test', full.names=T)

## snp.ex <- list.files('/home/ev250/rds/rds-mrc-bsu/ev250/martinPRS/model_QC',"failed",  full.names=T)
## traits <- c('BRCA', 'T2D', 'CAD')

names(snp.ex) <-  gsub("failed_", "", basename(snp.ex))
model <-  gsub(".score.txt", "", basename(c(score_rest, score_x)))
model_rest <- gsub(".score.txt", "", basename(c(score_rest)))
model_x <-  paste("prscsx", gsub(".score.txt", "", basename(score_x)), sep="_")

################################
## Perform SNP QC per model ##
################################

failed <- lapply(names(snp.ex), function(i) fread(grep(paste("failed", i, sep = "_"), snp.ex, value=T)))
good <- lapply(names(snp.ex), function(i) fread(grep(paste("pgs_ukbb", i, sep="_"), snp.kept, value=T)))

names(failed) <- names(good) <- names(snp.ex)

## Count failed SNPs per model

failed.dt <- rbindlist(failed, idcol="model")

if(nrow(failed.dt)) {

    count.failed <- rbindlist(lapply(names(snp.ex), function(i) data.table(model=i,
                                                                           different_alleles=failed[[i]][multiallelic == FALSE & info_score >=0.3,.N],
                                                                           multiallelic_alleles=length(unique(failed[[i]][multiallelic==TRUE, rsid])),
                                                                           low_info=failed[[i]][multiallelic==FALSE & info_score<0.3, .N])))


    count.failed[, total_failed := different_alleles +  multiallelic_alleles + low_info]

    total.good <- lapply(good, nrow)

    count.failed[ , total_okay:= unlist(total.good)][ ,per_failed:= round(total_failed*100/(total_failed + total_okay), 2)]

}

##################
## Benchmark scores ##
##################

scores.l <- lapply(model, function(i) fread(grep(paste0("/", i), c(score_rest, score_x), value=T)))
names(scores.l) <- c(model_rest, model_x)

trait <- unlist(lapply(names(scores.l), function(i) {
    v = strsplit(i, "_")[[1]]
    w <- v[(v %in% traits)]
    if(w == "T2D") w <- "TD2"
    return(w)
}))


names(trait) <- c(model_rest, model_x)

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


auc_r2 <- merge(auc, r2,  by='model', sort=F)

auc_r2[, trait:= trait]
setkey(auc_r2, trait, model)

##################
## save output ##
##################

write.table(auc_r2, paste0(out_dir, "/auc_r2.txt"), row.names=FALSE)

write.table(ninds, paste0(out_dir, "/N_ukbb.txt"), row.names=F)

if(exists("count.failed")) {
    write.table(count.failed, paste0(out_dir, "/snp_qc_ukbb.txt"), row.names=F)
}
