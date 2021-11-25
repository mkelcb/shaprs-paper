library(data.table)

##' Compute prscsx in test after model training. In this version the trait is scaled as z-score and the linear model is fit without intercept
##'
##' @param id_f file with column f.eid and trait name for inds to use in test set
##' @param model_files path to models with different ethnicities to train and test model
##' @param trait character with trait to predict
##' @param out file name to save score after model training in test set
##' @export
##' @return saves file with scores
##' compute_prscsx2()

compute_prscsx2 <- function(id_f, model_files, trait, out){

    ## ukbb calls T2D as TD2, adjust
    trait <- ifelse(trait == "T2D", "TD2", trait)

    ids <- fread(id_f)$f.eid
    all_l <- lapply(model_files, fread)

    names(all_l) <- sapply(model_files, basename)

    ##  merge and scale scores
    all <- Reduce(function(a,b) merge(a, b, by=names(all_l[[1]])[names(all_l[[1]]) != "score"], suffixes=names(all_l)), all_l)

    score_cols <- grep('score', names(all), value=T)
    all[, (score_cols) := lapply(score_cols, function(i) scale(get(i)))]

    ## Scale trait as z-score

    ztrait <- paste0("z_", trait)

    all[, eval(ztrait) := scale(get(trait))]

    ## Select train inds and merge for running linear regression

    train <- all[!f.eid %in% ids,]    

    ## Fit model to train set (lm without intercept)    

    form <- as.formula(paste(ztrait, "~", paste(score_cols, collapse = " + "), "-1"))
    
    model = lm(form, data = train)

    print(form)
    
    print(summary(model))

    ## Predict on test set

    test <- all[f.eid %in% ids,]

    test[, score := predict(model, test)]


    write.table(test[, (c(score_cols, ztrait)) := NULL], out, row.names=FALSE)

}

## id_f = '/home/ev250/rds/rds-mrc-bsu/ev250/martinPRS/benchmark/id.test.CAD.txt'

## model_files= c('/home/ev250/rds/rds-mrc-bsu/ev250/martinPRS/ukbb_score/ind_subset/EUR_JAP_CAD_EUR_PRSCSx_E.score.txt', '/home/ev250/rds/rds-mrc-bsu/ev250/martinPRS/ukbb_score/ind_subset/EUR_JAP_CAD_EAS_PRSCSx_E.score.txt')

## trait = "CAD"

compute_prscsx2(snakemake@input[['test_ids']], snakemake@input[['scores']], snakemake@wildcards[['trait']], snakemake@output[['out']] )

