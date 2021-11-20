blendedSE = sqrt(  (1-subpheno_CombinedPheno_blending$lFDR)^2 * subpheno_CombinedPheno_blending$se.x^2 +  subpheno_CombinedPheno_blending$se.y^2 * subpheno_CombinedPheno_blending$lFDR^2 + 2*subpheno_CombinedPheno_blending$lFDR*(1-subpheno_CombinedPheno_blending$lFDR)*cor(subpheno_CombinedPheno_blending$b.x,subpheno_CombinedPheno_blending$b.y)* subpheno_CombinedPheno_blending$se.x*subpheno_CombinedPheno_blending$se.y )

  
  
blendedSE_new = sqrt(  (1-subpheno_CombinedPheno_blending$lFDR)^2 * subpheno_CombinedPheno_blending$se.x^2 +  subpheno_CombinedPheno_blending$se.y^2 * subpheno_CombinedPheno_blending$lFDR^2 + 2*subpheno_CombinedPheno_blending$lFDR*(1-subpheno_CombinedPheno_blending$lFDR)*cor(subpheno_CombinedPheno_blending$b.x,subpheno_CombinedPheno_blending$b.y)* subpheno_CombinedPheno_blending$se.x*subpheno_CombinedPheno_blending$se.y )



cor_orig = cor(subpheno_CombinedPheno_blending$b.x,subpheno_CombinedPheno_blending$b.y)

cor_new = cor(subpheno_CombinedPheno_blending$b.x/subpheno_CombinedPheno_blending$se.x,subpheno_CombinedPheno_blending$b.y/subpheno_CombinedPheno_blending$se.y)







corrDiagLoc="C:/0Datasets/shaPRS/corrDiag/corrDiag"
corrDiag= read.table(corrDiagLoc, header = F)

t.test(corrDiag$V1,corrDiag$V2, paired = T)

meanDiff = mean (abs( corrDiag$V1-corrDiag$V2 ) )
sdDiff = sd (abs( corrDiag$V1-corrDiag$V2 ) )
plot(corrDiag$V1,corrDiag$V2, xlab = "cor(B1,B2)", ylab = "cor(B1/SE1,B2/SE2)", main = paste0("mean abs diff:",signif(meanDiff, digits=2)," ( " ,signif(sdDiff, digits=2)," )") )


