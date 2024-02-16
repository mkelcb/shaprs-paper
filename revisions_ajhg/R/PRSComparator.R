
# Formally evaluate if the r2 of 2 PRS methods are significantly different via r2redux
# also if 2 ROC curves are different (for binary traits)
# https://glassboxmedicine.com/2020/02/04/comparing-aucs-of-machine-learning-models-with-delongs-test/
#######################


# install.packages("r2redux")
library(pROC)
library(r2redux)
library(lmtest)


# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 3) {stop("not enough Arguments received")} 

gwasLoc = args[1]
gwasLoc2 = args[2]
outputLoc= args[3]
calcAUC = args[4] == "1"
hasHeader = is.na(args[5]) == F && args[5] == "1"
doAppend = is.na(args[6]) == F && args[6] == "1"
separator= ","
if ( is.na(args[7]) == F && args[7] == "1") { separator = "\t"}

print(paste0("hasHeader: ", hasHeader, " / doAppend: " , doAppend , " / calcAUC: ", calcAUC, " / separator: ", separator))





# 1. load phenos
gwasRes= read.table(gwasLoc, header = hasHeader, sep=separator) 
gwasRes2= read.table(gwasLoc2, header = hasHeader, sep=separator) 

# need to exlude rows that are NA in EITHER file
#gwasRes= read.table("C:/0LocalHPC/PRS_GXE/dummy/dummyPRS_PRS.sscore", header = T, sep="\t") 
#gwasRes2 = gwasRes
#gwasRes[1,3]  = NA
#gwasRes2[2,3]  = NA
NOT_NA = !is.na(gwasRes[,3]) & !is.na(gwasRes2[,3])
gwasRes = gwasRes[NOT_NA,] 
gwasRes2 = gwasRes2[NOT_NA,]


# also remove -9s as those indicate missing pheno for a .fam file
# gwasRes[1,2]  = -9
# gwasRes2[2,2]  = -9
NOT_MISSING = (gwasRes[,2] != -9) & (gwasRes2[,2] != -9)
gwasRes = gwasRes[NOT_MISSING,] 
gwasRes2 = gwasRes2[NOT_MISSING,]



# create data frame as trait, PRS1, PRS2
alldat = cbind(gwasRes[,2], gwasRes[,3], gwasRes2[,3])

# perform r2redux r2 significance test
nv = nrow(alldat)
r_obj = r_diff(alldat, 1, 2, nv)
#r_obj$r_based_p

# LRT, depending on binary or continous outcomes
if(calcAUC) {
nested <- glm( (gwasRes[,2]-1) ~gwasRes[,3], family = "binomial") # PRS1 nested # PRS1 nested # -1 as for logistic regression it must be 0-1, but PLINK is 1-2
nested_2 <- glm( (gwasRes[,2]-1) ~gwasRes2[,3], family = "binomial") # PRS2 nested
complex <- glm( (gwasRes[,2]-1) ~gwasRes[,3]+gwasRes2[,3], family = "binomial")
} else {
  nested <- lm( gwasRes[,2] ~gwasRes[,3]) 
  nested_2 <- lm(gwasRes[,2]~gwasRes2[,3]) # PRS2 nested
  complex <- lm(gwasRes[,2]~gwasRes[,3]+gwasRes2[,3])
}


# nested with the first
lrtestObj = lrtest(nested, complex)
PRS1_Lrt_p = lrtestObj$`Pr(>Chisq)`[2] # 0.01610404

# nested with the second
lrtestObj = lrtest(nested_2, complex)
PRS2_Lrt_p = lrtestObj$`Pr(>Chisq)`[2] # 0.05035712




testResults = cbind.data.frame(r_obj$r_based_p,PRS1_Lrt_p, PRS2_Lrt_p	)
colnames(testResults) = c("r2redux_p","LRT_p_PRS1","LRT_p_PRS2")

if(doAppend == F) {
  write.table(testResults,file=outputLoc, quote = F, row.names = F, sep = "\t",append=doAppend)
} else { 
  write.table(testResults,file=outputLoc, quote = F,col.names = F, row.names = F, sep = "\t",append=doAppend)
}

print(paste("written: PRS1 vs PRS2 r2redux diff p: ",signif(r_obj$r_based_p,3), " / LRT p (PRS1): ",signif(PRS1_Lrt_p,3), " / LRT p (PRS2): ",signif(PRS2_Lrt_p,3), "to", outputLoc ))



if(calcAUC) {
  category= gwasRes[,2]
  prediction= gwasRes[,3]
  prediction2= gwasRes2[,3]
  roc_obj <- roc(category, prediction)
  roc_obj2 <- roc(category, prediction2)
  
  delong = roc.test(roc_obj,roc_obj2,method=c("delong"))
  
  filen= paste0(outputLoc, "_AUC")
  if(doAppend == F) {
    write.table(delong$p.value,file=filen, quote = F, row.names = F, sep = "\t")
  } else { 
    write.table(delong$p.value,file=filen, quote = F,col.names = F, row.names = F, sep = "\t",append=doAppend)
    
  }
  print(paste("written: delong p-val",signif(delong$p.value,3), "to", filen ))
}







# 
# 
# 
# response<-c(0,0,0,0,0,0,1,1,1,1,1,1,1)
# modela<-c(0.1,0.2,0.05,0.3,0.1,0.6,0.6,0.7,0.8,0.99,0.8,0.67,0.5)
# modelb<-c(0.3,0.6,0.2,0.1,0.1,0.9,0.23,0.7,0.9,0.4,0.77,0.3,0.89)
# roca <- roc(response,modela)
# rocb<-roc(response,modelb)
# 
# delong = roc.test(roca,rocb,method=c("delong"))
# 
# ?roc.test
# print(paste0("delong p-val: ",delong$p.value))
# 
# 
# 
# # need to z-score both Y and .profile scores
# ?r2_var
# 
# View(dat1)
# 
# dat=dat1
# nv=length(dat$V1) # the number of observations, ie the sample size
# v1=c(1)
# output1=r2_var(dat1,v1,nv)
# nv=length(dat2$V1)
# v1=c(1)
# output2=r2_var(dat2,v1,nv)
# 
# r_obj = r_diff(dat1, 1, 2, nv)
# r_obj$r_based_p
# ?r_diff
# $var1
# [1] 0.0009247466
# 
# $var2
# [1] 0.0009238836
# 
# r2varObj = r2_var(dat1, 1, nv)
# r2varObj$lower_r2
# 
# ?r2_var
# 
# r2_var

# 
# # according to the description of their package, the dat1, is P+T PRSs, so they would not be independent, so comparing them via r_diff is OK
# function (dat, v1, v2, nv) 
# {
#   dat = scale(dat) # this z-scores it automatically
#   omat = cor(dat)
#   if (length(v1) == 1 & length(v2) == 1) { # if we only select 1 PRS each
#     ord = c(1, (1 + v1), (1 + v2)) # when se specift PRS1 to be '1', it actually selects the 2nd column, as the 1st col is assumed to be the pheno
#     m1 = lm(dat[, 1] ~ dat[, (1 + v1)])
#     s1 = summary(m1)
#     m2 = lm(dat[, 1] ~ dat[, (1 + v2)])
#     s2 = summary(m2)
#     R2 = s1$r.squared
#     mv2 = 1
#     t100 = (1/(nv) * (1 - R2)^2)
#     var1 = t100
#     R2 = s2$r.squared
#     mv2 = 1
#     t100 = (1/(nv) * (1 - R2)^2)
#     var2 = t100
#     dvr = (s1$r.squared^0.5 - s2$r.squared^0.5)
#     r_aoa = r_olkin1_2(omat[ord, ord], nv)
#     r_chi_dum = dvr^2/r_aoa
#     r_p3 = pchisq(r_chi_dum, 1, lower.tail = F)
#     r_uci = dvr + 1.96 * r_aoa^0.5
#     r_lci = dvr - 1.96 * r_aoa^0.5
#     z = list(r1 = s1$r.squared^0.5, r2 = s2$r.squared^0.5, 
#              var1 = var1, var2 = var2, var_diff = r_aoa, r_based_p = r_p3, 
#              r_based_p_one_tail = r_p3/2, mean_diff = dvr, upper_diff = r_uci, 
#              lower_diff = r_lci)
#     return(z)
#   }
# }

