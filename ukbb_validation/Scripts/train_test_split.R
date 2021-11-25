library(data.table)

##' Split ids in test and train
##'
##' @param id_f file with column f.eid and trait name
##' @param split numeric with proportion to split, defaults to 0.5
##' @param out file name to save ids for individuals in test set
##' @export
##' @return saves file with ids
##' split4test()

split4test <- function(id_f, split=0.5, out){

    ids <- fread(id_f)
    test <- ids[sample(1:nrow(ids), as.integer(nrow(ids)*split)), .(f.eid)]
    write.table(test, out, row.names=FALSE)

}

## id_f <- '/home/ev250/rds/rds-mrc-bsu/ev250/martinPRS/ukbb_score/ind_subset/Combined_EUR_LD_EUR_JAP_BRCA_E.score.txt'

## split=0.5

split4test(snakemake@input[['total']], snakemake@params[['split']], snakemake@output[['out']])
