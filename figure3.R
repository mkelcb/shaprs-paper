## code to produce figure 3

## files=list.files(".", pattern="EUR_JAP_.*_diagDat$")

## data is loaded from these files
files=c("EUR_JAP_asthma_diagDat", "EUR_JAP_BRCA_diagDat", "EUR_JAP_CAD_diagDat",
        "EUR_JAP_height_diagDat", "EUR_JAP_T2D_diagDat")

## sanity checks: beta should be positively correlated between EUR and JAP, more so for significant snps
library(data.table)
for(f in files) {
  d=fread(f)
  d[,z:=Beta_ShaPRS/SE_ShaPRS][,p:=pnorm(-abs(z))]
  cat(f,
      with(d, cor(Beta_A,Beta_B)),
      with(d[p<1e-6], cor(Beta_A,Beta_B)), "\n", sep="\t")
}

data=lapply(files,fread)
plots=lapply(seq_along(data), function(i) {
  data[[i]][,p:=pnorm(-abs(Beta_B/SE_B))*2]
  invisible(ggplot(data[[i]][p<1e-4,], aes(x=Beta_A,y=Beta_B)) +
            geom_point() + ggtitle(files[i]) +
            geom_smooth(method="lm",se=FALSE) +
            geom_abline(linetype="dashed"))
})
plot_grid(plotlist=plots)

## make figure 3
for(i.files in seq_along(files)) {
  message("loading ",files[i.files])
d=fread(files[i.files])
head(d)

dim(d)
d[,z:=Beta_ShaPRS/SE_ShaPRS][,p:=pnorm(-abs(z))]
with(d, cor(Beta_A,Beta_B))
with(d[p<1e-6], cor(Beta_A,Beta_B))

use=d$p < 1e-6 | 1:nrow(d) %in% sample(1:nrow(d),size=nrow(d)/100)
sum(use)
duse=d[use,]

library(ggplot2)
theme_set(theme_minimal())
ggplot(duse, aes(x=Beta_A,y=Beta_B,col=-log10(p))) + geom_point()

ggplot(duse, aes(x=lFDR,y=Beta_ShaPRS,col=z)) +
  scale_colour_gradient2(mid="lightyellow",midpoint=0) +
  geom_point() +
  geom_rug(sides="b")
dusep=duse[Beta_B>0]

library(cowplot)
library(seaborn)
dusep[,pcat:=cut(p,c(0,1e-8,1e-4,1),include.lowest=TRUE)]

with(dusep, table(p=p< 1e-4, lfdr=lFDR < .5))

col=function()
  scale_colour_gradient(limits=c(0,1),low="#BD4733",high="#195DBA")


A=ggplot(dusep, aes(x=lFDR,fill=pcat)) +
  geom_histogram() +
  labs(y="number of SNPs") +
  scale_fill_seaborn("P value",palette="pastel6",
                     labels=rev(c(expression("">=10^-4),expression("">=10^-8),expression(""<10^-8)))) +
  ## scale_y_log10("count (log scale)",breaks=c(100,1e+5,1e+8),labels=c(100,"100,000","100,000,000")) +
  theme(legend.position=c(0.1,0.8))
## A

A2=ggplot(duse, aes(x=lFDR,y=(Beta_ShaPRS-Beta_B),col=lFDR)) +
  col() +
  geom_point() +
  geom_rug(sides="b")

B=ggplot(duse[lFDR>0.5 & p > 1e-8], aes(x=Beta_B,y=Beta_ShaPRS,col=lFDR)) +
  col() +
  geom_point() +
  geom_smooth(fill="lightblue",method="lm") +
  labs(x="Beta_1") +
  geom_abline(linetype="dashed")

C=ggplot(duse[lFDR>0.5 & p > 1e-8], aes(x=SE_B,y=SE_ShaPRS,col=lFDR)) +
  col() +
  geom_point() +
  ## geom_smooth(fill="lightblue") +
  labs(x="SE_1") +
  geom_abline(linetype="dashed")

D=ggplot(duse[lFDR<0.5 | p < 1e-8], aes(x=Beta_B,y=Beta_ShaPRS,col=lFDR)) +
  col() +
  geom_point() +
  geom_smooth(fill="lightblue",method="lm") +
  labs(x="Beta_1") +
  geom_abline(linetype="dashed")

E=ggplot(duse[lFDR<0.5 | p < 1e-8], aes(x=SE_B,y=SE_ShaPRS,col=lFDR)) +
  col() +
  geom_point() +
  ## geom_smooth(fill="lightblue") +
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

final=plot_grid(A,
          plot_grid(bplot,NULL,rel_widths=c(3,.5)),
          plot_grid(cplot,legend,rel_widths=c(3,.5)),
          ncol=1,
          labels=c("a","",""),rel_heights=c(1,1,1))
  save_plot(paste0("shaprs-fig-3-dataset-",i.files,".png"),final,base_height=8,base_width=6)
}

system("cp shaprs-fig-3-dataset-1.png shaprs-fig-3.png") # just use first for illustration
