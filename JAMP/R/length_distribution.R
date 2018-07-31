# does not currently NOT suppots wrapped fasta files
Length_distribution <- function(sequFile=NA, out="", fastq=NA, col=c("#ffffb2", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a"), maxL=600){

ORIGscipen <- getOption("scipen")
options(scipen=10)

# comand to cound sequence length in each line
cmd <- paste("'{ print length($0) }' \"", sequFile, "\"", sep="")


A <- system2("awk", cmd, stdout=T, stderr=T) # count length
A <- as.numeric(A)

# detect empty file
if(length(A)==0){
message(paste("WARNING: The following file is empty and does NOT contain any sequenced. Thus it is not plotted: ", sequFile, sep=""))
} else {


# autodetect fastq
if(is.na(fastq)){
fastq <- sum(A[1:30][seq(3, 30, 4)])==7
}

# select reads only!
if(fastq){
B <- A[seq(2, length(A), 4)]
} else {
B <- A[seq(2, length(A), 2)]
}

NumberofReads <- length(B)
# count length dirstribution
temptab <- cbind(table(B))

data <- data.frame("length" =as.numeric(row.names(temptab)), "reads"=temptab, rel= round(temptab/sum(temptab)*100, 6), "col"=col[1], stringsAsFactors=F)


data$col[data$rel>0.01] <- col[2]
data$col[data$rel>0.1] <- col[3]
data$col[data$rel>1] <- col[4]
data$col[data$rel>10] <- col[5]


if(out!=""){
pdf(out, height=3.7, width=c(maxL/70+3))
}

par(mar=c(4,4,2,1))
plot(NULL, xlim=c(maxL*0.03, maxL*0.97), ylim=c(1, 99), xlab="Sequence length [bp]", ylab="proportions [%]", yaxt="n")

axis(2, las=2)
mtext(paste(sub(".*/", "", sequFile), ", Nuber of reads=", NumberofReads, sep=""), side=3, adj=0, line=0.5, cex=1.2, font=2)

step <- 1
for (i in data$length){
rect(i-0.5, 0, i+0.5, 100, col=data$col[step], border=F)
rect(i-0.5, 0, i+0.5, data$rel[step], col="black", border=F)
step <- 1 + step
}

meep <- data[order(data$rel, decreasing=T),]
text(meep$length[1], 95, meep$length[1], cex=1, col="Darkblue")

#legend

stop <- maxL*1.1

rect(0, 96, stop/100, 100, col=col[5])
rect(0, 91, stop/100, 95, col=col[4])
rect(0, 86, stop/100, 90, col=col[3])
rect(0, 81, stop/100, 85, col=col[2])
rect(0, 76, stop/100, 80, col=col[1])
rect(0, 71, stop/100, 75, col="white")

text(stop/75, 98, ">10%", adj=0,cex=0.8)
text(stop/75, 93, "10-1%", adj=0,cex=0.8)
text(stop/75, 88, "1-0.1%", adj=0,cex=0.8)
text(stop/75, 83, "0.1-0.01%", adj=0,cex=0.8)
text(stop/75, 78, "<0.01%", adj=0,cex=0.8)
text(stop/75, 73, "0%", adj=0,cex=0.8)



if(out!=""){
dev.off()
}
options(scipen= ORIGscipen)

}
}
