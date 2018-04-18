#
Sequences_lost <- function(Reads_in=NA, Reads_out, Sample_names, out="", abundance=T, rel=F, col=c("#b2df8a", "#fb9a99", "#33a02c",  "#fb9a99"), main=""){

ORIGscipen <- getOption("scipen")
options(scipen=10)

AddText <- !is.na(Reads_in[1])

if(is.na(Reads_in[1])){
Reads_in <- Reads_out
}

delta <- Reads_in-Reads_out

if(rel==T){
Reads_out <- Reads_out/Reads_in*100
delta <- delta/Reads_in*100
}

data <- as.matrix(rbind(Reads_out, delta))


if(out!=""){
pdf(out, height=c(length(Reads_in)+3)/4, width=9)
}

par(mar=c(3,8,2,1))
x <- barplot(data, col=col, border=0, horiz=T, xlim=c(0,max(data[1,]+data[2,])*1.1), ylim=c(0.042*length(Reads_in), 1.156*length(Reads_in)))
axis(2, x, labels= Sample_names, las=2, tick=F)
mtext(main, side=3, adj=0, line=0.5, cex=1.2, font=2)
#title(main=main, adj=0, line=0, font=2)

if(abundance==T){
text(max(data[1,]+data[2,])/100, x, round(data[1,], 2), adj=0, col=col[3])
if(AddText){
text(max(data[1,]+data[2,])/100+ c(data[1,]+data[2,]), x, round(data[2,], 2), adj=0, col= col[4])
}
}


if(out!=""){
dev.off()
}
options(scipen= ORIGscipen)
}
