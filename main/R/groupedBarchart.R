# Grouped barchart with error bars
#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible

San_Diego <- c(65,20,74)
Rosarito <- c(34,35,23)
La_Paz <- c(21,71,28)
Mating_strategy <- c("Ultradominant","Dominant","Sneaker")
col <- c("darkorange1","skyblue3","gold2")        

lizards <- data.frame(row.names=Mating_strategy, San_Diego, 
                      Rosarito, La_Paz)

lizards.matrix <- as.matrix(lizards)
# create matrix of errors
lizards.error = matrix(c(10, 5, 12, 10, 8, 6, 12, 28, 3), ncol = 3)


# create bar plot
bp <- barplot(lizards.matrix,
              beside=T,
              col=col,
              ylim=c(0,100),
              xlab="Site",ylab="Frequency",
              legend.text=row.names(lizards.matrix),
              args.legend=list(x="top",bty="n"),
              las=1,
              cex.axis=1.2)



# add error bars
arrows(bp, lizards.matrix - lizards.error, bp, lizards.matrix + lizards.error,
       code = 3, angle = 90, length = 0.15)

