## code to produce figure 3
## system("wget https://github.com/mkelcb/shaprs-paper/raw/master/Diag/Qvals_diagDat_CD_UC.zip && unzip Qvals_diagDat_CD_UC.zip")
files=list.files(".", pattern="EUR_JAP_.*_diagDat$")

## data is loaded from these files
files=c("EUR_JAP_asthma_diagDat_relab123456", "EUR_JAP_BRCA_diagDat", "EUR_JAP_CAD_diagDat",
        "EUR_JAP_height_diagDat", "EUR_JAP_T2D_diagDat", "CD_UC_diag")

dold=fread("bak/Qvals_diagDat")
dnew=fread("CD_UC_diag")
par(mfrow=c(1,2))
hist(dold$lFDR,main="Qvals_diagDat")
hist(dnew$lFDR,main="CD_UC_diag")


library(data.table)
reader=function(f) {
  d=fread(f)
  if(f=="CD_UC_diag") {
    d[,Beta_B:=Beta_CD][,SE_B:=SE_CD][,Beta_ShaPRS:=Beta_ShaPRS_CD][,SE_ShaPRS:=SE_ShaPRS_CD]
    d[,Beta_A:=Beta_UC][,SE_A:=SE_UC]
  }
  if(f=="EUR_JAP_asthma_diagDat_relab123456") {
    d[,Beta_B:=Beta_EUR][,SE_B:=SE_EUR]
    d[,Beta_A:=Beta_JP][,SE_A:=SE_JP]
  }
  d[,z:=Beta_ShaPRS/SE_ShaPRS][,p:=pnorm(-abs(z))*2]
  d[,zA:=Beta_A/SE_A][,pA:=pnorm(-abs(zA))]
  d[,zB:=Beta_B/SE_B][,pB:=pnorm(-abs(zB))]
  d[,het:=ifelse(lFDR < 0.5, 0, 1)]
  d[,pcat:=cut(p,c(0,1e-8,1),include.lowest=TRUE)]
  levels(d$pcat) <-rev(c(expression("non-significant ">=5 *"x"* 10^-8),
                                   expression("significant "<5 * "x" * 10^-8)))
  copy(d)
}

library(ggplot2)
library(cowplot)
library(seaborn)
data=lapply(files,reader)

## sanity checks: beta should be positively correlated between EUR and JP, more so for significant snps
for(i in seq_along(data)) {
  d=data[[i]]
  cat(files[i],
      nrow(d),
      with(d, cor(Beta_A,Beta_B)),
      with(d[p<1e-6], cor(Beta_A,Beta_B)),
      with(d[pA<1e-6], cor(Beta_A,Beta_B)),
      with(d[pB<1e-6], cor(Beta_A,Beta_B)),
      "\n", sep="\t")
}

## plots=lapply(seq_along(data), function(i) {
##   data[[i]][,pA:=pnorm(-abs(Beta_A/SE_A))*2]
##   data[[i]][,pB:=pnorm(-abs(Beta_B/SE_B))*2]
##   invisible(ggplot(data[[i]][pmin(pA,pB)<1e-4,], aes(x=Beta_A,y=Beta_B,col=lFDR)) +
##             geom_point() + ggtitle(files[i]) +
##             geom_smooth(method="lm",se=FALSE) +
##             geom_abline(linetype="dashed"))
## })
## plot_grid(plotlist=plots)

## make figure 4
plotter=function(i.files, thr=1e-6) {
  message("loading ",files[i.files])
  d=data[[i.files]]
  head(d)

  dim(d)
  with(d, cor(Beta_A,Beta_B))
  with(d[p<1e-6], cor(Beta_A,Beta_B))
  with(d, table(pcat, het))

  ## subsample
  use=d$p < 5e-8 | d$lFDR < .5 | 1:nrow(d) %in% sample(1:nrow(d),size=nrow(d)/100)
  sum(use)
  duse=d[use,]
  ## duse[,het:=ifelse(p < 5e-8, max(-log10(duse$p)), 0)]

  theme_set(theme_classic())
  ## ggplot(d, aes(x=Beta_A,y=Beta_B,col=-log10(p))) + geom_point()

  col=function()
    scale_fill_manual("ShaPRS P value",#palette="pastel6",
                      values=rev(c("grey80","skyblue")),
                      labels=rev(c(expression("non-significant ">=5 *"x"* 10^-8),
                                   expression("significant "<5 * "x" * 10^-8))))
      col2=function()
        scale_colour_manual("ShaPRS P value",#palette="pastel6",
                      values=rev(c("grey40","deepskyblue")),
                      labels=rev(c(expression("non-significant ">=5 *"x"* 10^-8),
                                   expression("significant "<5 * "x" * 10^-8))))
  col3=function()
    scale_fill_manual("ShaPRS P value",#palette="pastel6",
                      values=rev(c(alpha("grey80",0.2),alpha("skyblue",0.2))),
                      labels=rev(c(expression("non-significant ">=5 *"x"* 10^-8),
                                   expression("significant "<5 * "x" * 10^-8))))
  ## col=function()
  ##   scale_colour_gradient("",limits=c(0,1),low="#BD4733",high="#195DBA")

  A=ggplot(d, aes(x=lFDR,fill=pcat)) +
    geom_histogram(col="black") +
    ## geom_histogram(aes(col=lFDR)) +
    ## col() +
    labs(y="number of SNPs") +
    col() + col2() +
    facet_wrap(pcat ~ ., scales="free",ncol=1) +
    ## scale_y_log10() +
    ## scale_y_log10("count (log scale)",breaks=c(100,1e+5,1e+8),labels=c(100,"100,000","100,000,000")) +
    theme(legend.position="none",#c(0.5,0.3),
          legend.title=element_text(size=rel(0.9)),
          legend.text=element_text(size=rel(0.9)),
          strip.background=element_blank(),
          strip.text=element_blank())
  A

  ## A2=ggplot(duse, aes(x=lFDR,y=(Beta_ShaPRS-Beta_B),col=lFDR)) +
  ##   col() +
  ##   geom_point() +
  ##   geom_rug(sides="b")
  B=ggplot(duse[lFDR > .5], aes(x=Beta_B,y=Beta_ShaPRS,col=pcat)) +
    col3() + col2() +
    geom_point(alpha=0.5) +
    xlim(min(duse$Beta_B),max(duse$Beta_B)) +
    ylim(min(duse$Beta_ShaPRS),max(duse$Beta_ShaPRS)) +
    ## geom_smooth(fill="lightblue",method="lm") +
    geom_smooth(aes(group=pcat,col=pcat),method="lm",se=FALSE) +
    ## geom_smooth(col="black",method="lm",se=FALSE) +
    labs(x="Beta_1") +
    background_grid() +
    geom_abline(linetype="dashed")
  B

  C=ggplot(duse[lFDR > .5], aes(x=SE_B,y=SE_ShaPRS,col=pcat)) +
    col3() + col2() +
    geom_point(alpha=0.2) +
    ## geom_smooth(fill="lightblue") +
    geom_smooth(aes(group=pcat,col=pcat),method="lm",se=FALSE) +
    ## geom_smooth(col="black",method="lm",se=FALSE) +
    scale_x_continuous(breaks=if(max(duse$SE_B)> 0.09) { c(.05,.1) } else { c(.02, .04, .06) },
                       limits=c(min(duse$SE_B),max(duse$SE_B))) +
    ylim(min(duse$SE_ShaPRS),max(duse$SE_ShaPRS)) +
    labs(x="SE_1") +
    background_grid() +
    geom_abline(linetype="dashed")
  C

  D=ggplot(duse[lFDR<=.5], aes(x=Beta_B,y=Beta_ShaPRS,col=pcat)) +
    col3() + col2() +
    geom_point(alpha=0.2) +
    xlim(min(duse$Beta_B),max(duse$Beta_B)) +
    ylim(min(duse$Beta_ShaPRS),max(duse$Beta_ShaPRS)) +
    geom_smooth(aes(group=pcat,col=pcat),method="lm",se=FALSE) +
    ## geom_smooth(fill="lightblue",method="lm") +
    labs(x="Beta_1") +
    background_grid() +
    geom_abline(linetype="dashed")
  D

  E=ggplot(duse[lFDR<=.5], aes(x=SE_B,y=SE_ShaPRS,col=pcat)) +
    col3() + col2() +
    geom_point(alpha=0.1) +
    scale_x_continuous(breaks=if(max(duse$SE_B)> 0.09) { c(.05,.1) } else { c(.02, .04, .06) },
                       limits=c(min(duse$SE_B),max(duse$SE_B))) +
    ylim(min(duse$SE_ShaPRS),max(duse$SE_ShaPRS)) +
    ## geom_smooth(fill="lightblue") +
    geom_smooth(aes(group=pcat,col=pcat),method="lm",se=FALSE) +
    background_grid() +
    labs(x="SE_1") +
    geom_abline(linetype="dashed")
  E
  ## bplot=plot_grid(B+theme(legend.position="none"),
  ##                 C+theme(legend.position="none"),ncol=2,labels=letters[2:3])
  ## cplot=plot_grid(D+theme(legend.position="none"),
  ##                 E+theme(legend.position="none"),ncol=2,labels=letters[4:5])
  ## legend= get_legend(
  ##                                       # create some space to the left of the legend
  ##   B + theme(legend.box.margin = margin(0, 0, 0, 12))
  ## )

  list(A=A,B=B,C=C,D=D,E=E)#, bplot=bplot,cplot=cplot)
}

plots=lapply(seq_along(files), plotter )
lapply(data, function(x) sum(x$lFDR < 0.01))

## check histograms make sense
library(magrittr)
A=lapply(plots, "[[", 1)
for(i in seq_along(files))
  A[[i]]=A[[i]] + ggtitle(files[[i]])
plot_grid(plotlist=A)
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
  plots[[1]]$A + theme(legend.position="bottom",legend.box.margin = margin(0, 0, 0, 12))
)
top=plot_grid(plots[[1]]$A + ggtitle("Asthma in EUR, adjunct: EAS"),
              plots[[6]]$A + ggtitle("Crohn's, adjunct: UC"),nrow=1, labels=letters[1:2])
botleft=plot_grid(plots[[1]]$B +theme(legend.position="none"),
                  plots[[1]]$C +theme(legend.position="none"),
              nrow=1, labels=letters[c(3:4)])
botright=plot_grid(plots[[1]]$D +theme(legend.position="none"),
                  plots[[1]]$E +theme(legend.position="none"),
                  nrow=1, labels=letters[c(5:6)])
titleleft <- ggdraw() + draw_label("Asthma, homogeneous SNPs", x=0.1, hjust=0)
titleright <- ggdraw() + draw_label("Asthma, heterogeneous SNPs", x=0.1, hjust=0)
bot=plot_grid(titleleft,titleright,
              botleft , botright,
              nrow=2,
              rel_heights=c(.1,.9))
## mid=plot_grid(plots[[1]]$B +theme(legend.position="none"),
##               plots[[1]]$C +theme(legend.position="none"),
##               plots[[6]]$B +theme(legend.position="none"),
##               plots[[6]]$C +theme(legend.position="none"),
##               nrow=1, labels=letters[c(3,4,7,8)])
## bot=plot_grid(plots[[1]]$D +theme(legend.position="none"),
##               plots[[1]]$E +theme(legend.position="none"),
##               plots[[6]]$D +theme(legend.position="none"),
##               plots[[6]]$E +theme(legend.position="none"),
##               nrow=1, labels=letters[c(5,6,9,10)])

final=plot_grid(top,bot,plot_grid(legend),nrow=3,rel_heights=c(5,3,.5))
## final=plot_grid(top,mid,bot,plot_grid(legend),nrow=4,rel_heights=c(3,3,3,.5)) +
##   draw_label("lFDR: heterogeneous", color = "black", size = 10, x=.35,y=.03) +
##   draw_label("homogeneous", color = "black", size = 10, x=.65,y=.03)
save_plot(paste0("shaprs-fig-4.png"),plot=final,base_height=8,base_width=8)
final
