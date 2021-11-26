## code to produce figure 3

files=list.files(".", pattern="EUR_JAP_.*_diagDat$")

## data is loaded from these files
files=c("EUR_JAP_asthma_diagDat", "EUR_JAP_BRCA_diagDat", "EUR_JAP_CAD_diagDat",
        "EUR_JAP_height_diagDat", "EUR_JAP_T2D_diagDat", "IBD_Qvals_new_diagDat.gz")

## sanity checks: beta should be positively correlated between EUR and JAP, more so for significant snps
library(data.table)
for(f in files) {
  d=fread(f)
  d[,z:=Beta_ShaPRS/SE_ShaPRS][,p:=pnorm(-abs(z))]
  d[,zA:=Beta_A/SE_A][,pA:=pnorm(-abs(zA))]
  d[,zB:=Beta_B/SE_B][,pB:=pnorm(-abs(zB))]
  cat(f,
      with(d, cor(Beta_A,Beta_B)),
      with(d[p<1e-6], cor(Beta_A,Beta_B)),
      with(d[pA<1e-6], cor(Beta_A,Beta_B)),
      with(d[pB<1e-6], cor(Beta_A,Beta_B)),
      "\n", sep="\t")
}

library(ggplot2)
library(cowplot)
library(seaborn)
data=lapply(files,fread)
plots=lapply(seq_along(data), function(i) {
  data[[i]][,pA:=pnorm(-abs(Beta_A/SE_A))*2]
  data[[i]][,pB:=pnorm(-abs(Beta_B/SE_B))*2]
  invisible(ggplot(data[[i]][pmin(pA,pB)<1e-4,], aes(x=Beta_A,y=Beta_B,col=lFDR)) +
            geom_point() + ggtitle(files[i]) +
            geom_smooth(method="lm",se=FALSE) +
            geom_abline(linetype="dashed"))
})
plot_grid(plotlist=plots)

## make figure 3
plotter=function(i.files) {
  message("loading ",files[i.files])
  d=data[[i.files]]
  head(d)

  dim(d)
  d[,z:=Beta_ShaPRS/SE_ShaPRS][,p:=pnorm(-abs(z))]
  with(d, cor(Beta_A,Beta_B))
  with(d[p<1e-6], cor(Beta_A,Beta_B))

  use=d$p < 1e-6 | 1:nrow(d) %in% sample(1:nrow(d),size=nrow(d)/100)
  sum(use)
  duse=d[use,]
  duse[,het:=ifelse(lFDR < 0.5, 0, 1)]
  ## duse[,het:=ifelse(p < 5e-8, max(-log10(duse$p)), 0)]

  theme_set(theme_classic())
  ggplot(duse, aes(x=Beta_A,y=Beta_B,col=-log10(p))) + geom_point()

  ## ggplot(duse, aes(x=lFDR,y=Beta_ShaPRS,col=z)) +
  ##   scale_colour_gradient2(mid="lightyellow",midpoint=0) +
  ##   geom_point() +
  ##   geom_rug(sides="b")
  dusep=duse[Beta_B>0]

  dusep[,pcat:=cut(p,c(0,1e-8,1),include.lowest=TRUE)]

  with(dusep, table(p=p< 1e-4, lfdr=lFDR < .5))

  col=function()
    scale_colour_gradient("",limits=c(0,1),low="#BD4733",high="#195DBA")


  A=ggplot(dusep, aes(x=lFDR,fill=pcat)) +
    geom_histogram(col="black") +
    ## geom_histogram(aes(col=lFDR)) +
    ## col() +
    labs(y="number of SNPs") +
    scale_fill_manual("shaPRS P value",#palette="pastel6",
                      values=c("grey","white"),
                      labels=rev(c(expression("non-significant ">=5 *"x"* 10^-8),
                                   expression("significant "<5 * "x" * 10^-8)))) +
    ## scale_y_log10("count (log scale)",breaks=c(100,1e+5,1e+8),labels=c(100,"100,000","100,000,000")) +
    theme(legend.position=c(0.3,0.8))
  A

  ## A2=ggplot(duse, aes(x=lFDR,y=(Beta_ShaPRS-Beta_B),col=lFDR)) +
  ##   col() +
  ##   geom_point() +
  ##   geom_rug(sides="b")
  B=ggplot(duse[p > 5e-8 & lFDR > .5], aes(x=Beta_B,y=Beta_ShaPRS,col=lFDR)) +
    col() +
    geom_point(alpha=0.5) +
    xlim(min(duse$Beta_B),max(duse$Beta_B)) +
    ylim(min(duse$Beta_ShaPRS),max(duse$Beta_ShaPRS)) +
    ## geom_smooth(fill="lightblue",method="lm") +
    geom_smooth(aes(group=het,col=het),method="lm",se=FALSE) +
    ## geom_smooth(col="black",method="lm",se=FALSE) +
    labs(x="Beta_1") +
    background_grid() +
    geom_abline(linetype="dashed")
  B

  C=ggplot(duse[p > 5e-8 & lFDR > .5], aes(x=SE_B,y=SE_ShaPRS,col=lFDR)) +
    col() +
    geom_point(alpha=0.5) +
    ## geom_smooth(fill="lightblue") +
    geom_smooth(aes(group=het,col=het),method="lm",se=FALSE) +
    ## geom_smooth(col="black",method="lm",se=FALSE) +
    scale_x_continuous(breaks=if(max(duse$SE_B)> 0.09) { c(.05,.1) } else { c(.02, .04, .06) },
                       limits=c(min(duse$SE_B),max(duse$SE_B))) +
    ylim(min(duse$SE_ShaPRS),max(duse$SE_ShaPRS)) +
    labs(x="SE_1") +
    background_grid() +
    geom_abline(linetype="dashed")

  D=ggplot(duse[p <= 5e-8 | lFDR<=.5], aes(x=Beta_B,y=Beta_ShaPRS,col=lFDR)) +
    col() +
    geom_point(alpha=0.5) +
    xlim(min(duse$Beta_B),max(duse$Beta_B)) +
    ylim(min(duse$Beta_ShaPRS),max(duse$Beta_ShaPRS)) +
    geom_smooth(aes(group=het,col=het),method="lm",se=FALSE) +
    ## geom_smooth(fill="lightblue",method="lm") +
    labs(x="Beta_1") +
    background_grid() +
    geom_abline(linetype="dashed")
  D

  E=ggplot(duse[p <= 5e-8 | lFDR<=.5], aes(x=SE_B,y=SE_ShaPRS,col=lFDR)) +
    col() +
    geom_point(alpha=0.5) +
    scale_x_continuous(breaks=if(max(duse$SE_B)> 0.09) { c(.05,.1) } else { c(.02, .04, .06) },
                       limits=c(min(duse$SE_B),max(duse$SE_B))) +
    ylim(min(duse$SE_ShaPRS),max(duse$SE_ShaPRS)) +
    ## geom_smooth(fill="lightblue") +
    geom_smooth(aes(group=het,col=het),method="lm",se=FALSE) +
    background_grid() +
    labs(x="SE_1") +
    geom_abline(linetype="dashed")
  bplot=plot_grid(B+theme(legend.position="none"),
                  C+theme(legend.position="none"),ncol=2,labels=letters[2:3])
  cplot=plot_grid(D+theme(legend.position="none"),
                  E+theme(legend.position="none"),ncol=2,labels=letters[4:5])
  legend= get_legend(
                                        # create some space to the left of the legend
    B + theme(legend.box.margin = margin(0, 0, 0, 12))
  )

  list(A=A,B=B,C=C,D=D,E=E, bplot=bplot,cplot=cplot)
}

plots=lapply(seq_along(files), plotter )

## left=with(plots[[1]],
##      plot_grid(A + ggtitle("Asthma in EUR, adjunct: EAS"),
##                   plot_grid(bplot,NULL,rel_widths=c(3,.5)),
##                   plot_grid(cplot,legend,rel_widths=c(3,.5)),
##                   ncol=1,
##                labels=c("a","",""),rel_heights=c(1,1,1)))

## right=with(plots[[6]],
##      plot_grid(A + ggtitle("UC, adjunct: Crohn's"),
##                   plot_grid(bplot,NULL,rel_widths=c(3,.5)),
##                   plot_grid(cplot,legend,rel_widths=c(3,.5)),
##                   ncol=1,
##                labels=c("a","",""),rel_heights=c(1,1,1))

## plot_grid(left,right)

legend= get_legend(
  ## create some space to the left of the legend
  plots[[1]]$B + theme(legend.position="bottom",legend.box.margin = margin(0, 0, 0, 12))
)
top=plot_grid(plots[[1]]$A + ggtitle("Asthma in EUR, adjunct: EAS"),
              plots[[6]]$A + ggtitle("UC, adjunct: Crohn's"),nrow=1, labels=letters[1:2])
mid=plot_grid(plots[[1]]$B +theme(legend.position="none"),
              plots[[1]]$C +theme(legend.position="none"),
              plots[[6]]$B +theme(legend.position="none"),
              plots[[6]]$C +theme(legend.position="none"),
              nrow=1, labels=letters[c(3,4,7,8)])
bot=plot_grid(plots[[1]]$D +theme(legend.position="none"),
              plots[[1]]$E +theme(legend.position="none"),
              plots[[6]]$D +theme(legend.position="none"),
              plots[[6]]$E +theme(legend.position="none"),
              nrow=1, labels=letters[c(5,6,9,10)])

final=plot_grid(top,mid,bot,plot_grid(legend),nrow=4,rel_heights=c(3,3,3,.5)) +
  draw_label("heterogeneous", color = "black", size = 10, x=.35,y=.03) +
  draw_label("homogeneous", color = "black", size = 10, x=.65,y=.03)
save_plot(paste0("shaprs-fig-4.png"),plot=final,base_height=10,base_width=8)

final

## Legend
"
Two examples of shaPRS analyses. The columns show a cross-ancestry analysis of asthma (left), and a cross-trait analysis of UC, leveraging a GWAS of Crohn's disease as an adjunct dataset. a, b show the distribution of lFDR values, where low lFDR corresponds to higher heterogeneity in estimated effects. The middle row compares the input beta (Beta_1) and standard error (SE_1) to its shaPRS-adjusted output (Beta_shaPRS, SE_shaPRS respectively) for SNPs which are both not significant (shaPRS p < 5x10-8) and have higher homogeneity (lFDR > 0.5).  The bottom row shows the same for the remainder of the SNPs which either have heterogeneous effects (lFDR < 0.5) or are significant (shaPRS p < 5x10-8).  Points are coloured by their lFDR value, with lower values representing greater likelihood of heterogeneous effects.
"

## Narrative
"
ShaPRS adapts its behaviour to the pattern of genetic sharing in the studies (Fig 4).
In a cross ancestry analysis, very few SNPs are detected to have genuinely different effects (ie low lFDR) whilst in the cross-trait analysis.
For the majority of SNPs with high homogeneity (lFDR > 0.5), standard errors are shrunk by shaPRS, whilst coefficients are shrunk towards 0 for non-significant SNPs (shaPRS p > 5x10-8) with higher homogeneity (lFDR < 0.5) but left unchanged otherwise.
In the cross-trait analysis,  the lFDR distribution is shifted to the left. Standard errors remain broadly unchanged because both studies already share controls, but beta coefficients adapt at both non-significant and non-significant SNPs.
"

