# FastQC

FastQC <- function(files="latest" ){

# get folder path!
if(file.exists("log.txt")){
log <- readLines("log.txt") # load log
step <- which(log=="PROCESSING MODULE:")+1 # detect 

folder <- log[step[length(step)]]
message("Runing module: FastQC")

} else { # if nothing has been run jet
message("Runing module: FastQC")
dir.create("_FastQC", showWarnings=F)
folder <- "_FastQC"
}



if(files=="latest"[1]){
files <- list.files(paste(folder, "/_data", sep=""), full.name=T)
}


dir.create(paste(folder, "/FastQC/_raw_reports/", sep=""), showWarnings=F, recursive=T)

# run fastQC module
cmd <- paste("-o ", folder, "/FastQC/_raw_reports/ ", sep="", paste("\"", files, "\" ", collapse="", sep="") )

A <- system2("fastQC", cmd)

#system2("vsearch", paste("--fastq_eestats ", files[1], sep="", " --output test.txt"))


# run fastqcr
QCpath <- paste(folder, "/FastQC/_raw_reports/", sep="")
QCfiles <- list.files(QCpath, full.names=T, pattern="fastqc.zip")

QClist <- list()
for(i in 1:length(QCfiles)){
QClist[[i]] <- qc_read(QCfiles[i])
}


# plot for each sample
dir.create(paste(folder, "/FastQC/_plots_for_each_sample/", sep=""), showWarnings=F)



for (i in 1:length(QClist)){

temp_file <- paste(folder, "/FastQC/_plots_for_each_sample/", sep="" , sub(".*/(.*).fastq", "\\1.pdf", files[i]))
pdf(temp_file, height=8, width=10, useDingbats=F)

qc_plot(QClist[[i]], "Basic statistics")

#qc_plot(QClist[[i]], "Per base sequence quality")
A <- data.frame(QClist[[i]]$per_base_sequence_quality)
plot(NULL, ylim=c(1, 42), xlim=c(1, nrow(A)), xlab="Possition in Read (bp)", ylab="Phred score", main="Per base sequence quality", xaxt="n")
select <- c(1, 5, seq(10, 60, 10))
axis(1, labels=A$Base[select], at= select)

rect(-10,44, 1000, 28, border=F, col="White")
rect(-10,20, 1000, 28, border=F, col="Gray90")
rect(-10,20, 1000, -10, border=F, col="Gray80")

k <- 1
for (k in 1:nrow(A)){
lines(c(k,k), c(A$X90th.Percentile[k], A$X10th.Percentile[k]), lwd=1, lend=1)
lines(c(k+0.4,k-0.4), c(A$X90th.Percentile[k], A$X90th.Percentile[k]), lwd=1, lend=1)
lines(c(k+0.4,k-0.4), c(A$X10th.Percentile[k], A$X10th.Percentile[k]), lwd=1, lend=1)
rect(k+0.4, A$Median[k], k-0.4, A$Upper.Quartile[k], col="lightblue")
rect(k+0.4, A$Median[k], k-0.4, A$Lower.Quartile[k], col="lightblue")
lines(c(k+0.4,k-0.4), c(A$Median[k], A$Median[k]), lwd=1.5, lend=1, col="Red")

}
lines(1:nrow(A), A$Mean, col="black", lwd=2)


#qc_plot(QClist[[i]], "Per sequence quality scores")

A <- data.frame(QClist[[i]]$per_sequence_quality_scores)
B <- data.frame("Quality"=1:40, "Count"=rep(0, 40))
B$Count[match(A$Quality, B$Quality)] <- A$Count

plot(NULL, xlim=c(1,40), ylim=range(B$Count), xlab="Mean Sequence Quality (Phred Score)", ylab="Counts", main="Per base sequence quality")
for(k in 1:nrow(B)){
if(B$Count[k]>0){
rect(k-0.5, 0, k+0.5, B$Count[k], border=F, col=if(B$Quality[k]<30){"Red"}else{"Green"})
}
}
Q30 <- sum(B$Count[B$Quality>=30])/sum(B$Count)*100
text(30, max(B$Count)*0.8, paste("%<=Q30\n", round(Q30, 2)), adj=0)

#qc_plot(QClist[[i]], "Sequence length distribution")
A <- data.frame(QClist[[i]]$sequence_length_distribution)
if(nrow(A)==1){A <- rbind(c(A$Length-1, 0), A, c(A$Length+1, 0))}
plot(A, typ="l", main="Sequence length distribution", xlab="Sequence length in bp", ylab="Count")

#qc_plot(QClist[[i]], "Per base sequence content")
A <- data.frame(QClist[[i]]$per_base_sequence_content)
plot(NULL, ylim=c(0, 100), xlim=c(1, nrow(A)), xlab="", ylab="Counts", main="Per base sequence content", xaxt="n")
axis(1, labels=A$Base[c(1, 5, seq(10, 60, 10))], at=c(1, 5, seq(10, 60, 10)))
polygon(c(1, 1:nrow(A), nrow(A)), c(0,A$G+A$A+A$T+A$C,0), border=F, col="Blue")
polygon(c(1, 1:nrow(A), nrow(A)), c(0,A$G+A$A+A$T,0), border=F, col="Green")
polygon(c(1, 1:nrow(A), nrow(A)), c(0,A$G+A$A,0), border=F, col="Red")
polygon(c(1, 1:nrow(A), nrow(A)), c(0,A$G,0), border=F, col="Yellow")

#qc_plot(QClist[[i]], "Adapter content")
A <- data.frame(QClist[[i]]$adapter_content)
plot(NULL, ylim=c(1, 100), xlim=c(1, nrow(A)), xlab="Possition in Read (bp)", ylab="% Adapter", main="Adapter content", xaxt="n")
select <- c(1, 5, seq(10, 60, 10))
axis(1, labels=A$Position[select], at= select)

for (k in 2:ncol(A)){
lines(1:nrow(A), A[,k], col="black", lwd=1)
}


#qc_plot(QClist[[i]], "Per base N content")
A <- data.frame(QClist[[i]]$per_base_n_content)
plot(NULL, ylim=c(1, 100), xlim=c(1, nrow(A)), xlab="Possition in Read (bp)", ylab="% N", main="Per base N content", xaxt="n")
select <- c(1, 5, seq(10, 60, 10))
axis(1, labels=A$Base[select], at= select)

lines(1:nrow(A), A[,2], col="Red", lwd=1.5)


qc_plot(QClist[[i]], "Per sequence GC content")
qc_plot(QClist[[i]], "Sequence duplication levels")
#qc_plot(QClist[[i]], "Overrepresented sequences")
#qc_plot(QClist[[i]], "Kmer content")

dev.off()

}




# aggregate




#message(" ")
#message(" Module completed!")

#cat(file="log.txt", paste(Sys.time(), "\n", "*** Module completed!", "", sep="\n"), append=T, sep="\n")


}

