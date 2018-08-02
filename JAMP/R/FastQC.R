# FastQC

FastQC <- function(files="latest", space=10 ){

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



if(files[1]=="latest"){
files <- list.files(paste(folder, "/_data", sep=""), full.name=T)
}


dir.create(paste(folder, "/FastQC/_raw_reports/", sep=""), showWarnings=F, recursive=T)

# run fastQC module
cmd <- paste("-o ", folder, "/FastQC/_raw_reports/ ", sep="", paste("\"", files, "\" ", collapse="", sep="") )

A <- system2("fastQC", cmd)

#system2("vsearch", paste("--fastq_eestats ", files[1], sep="", " --output test.txt"))


# run fastqcr
#QCpath <- paste(folder, "/FastQC/_raw_reports/", sep="")
QCfiles <- paste(folder, "/FastQC/_raw_reports/", sep="", gsub(".*/(.*).fastq.*", "\\1_fastqc.zip", files))

#QCfiles <- list.files(QCpath, full.names=T, pattern="fastqc.zip")

QClist <- list()
for(i in 1:length(QCfiles)){
QClist[[i]] <- qc_read(QCfiles[i])
}


# plot for each sample
dir.create(paste(folder, "/FastQC/_plots_for_each_sample/", sep=""), showWarnings=F,  recursive=T)

# summary table
temp <- rep(NA, length(QClist))
exp <- data.frame("ID"= temp, "SequDepth"= temp, "Length"= temp, "LenghtRel"= temp, "MaxAdapters"= temp, "MaxNRel"= temp, "Diversity"= temp, "GC"=temp, "Q30"=temp, "duplicated"=temp)

# PLOTTING
#
for (i in 1:length(QClist)){

temp_file <- paste(folder, "/FastQC/_plots_for_each_sample/", sep="" , sub(".*/(.*)\\.fastq.*", "\\1.pdf", files[i]))
pdf(temp_file, height=8, width=10, useDingbats=F)

qc_plot(QClist[[i]], "Basic statistics")

#summary stats
exp[i,1] <- sub(".*/", "", files[i]) # file name
exp[i,2] <- QClist[[i]]$basic_statistics[4,2] # sequ depth
exp[i,3] <- QClist[[i]]$basic_statistics[6,2]
exp[i,8] <- QClist[[i]]$basic_statistics[7,2]
exp$duplicated[i] <- 100-QClist[[i]]$total_deduplicated_percentage


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
B <- data.frame("Quality"=1:60, "Count"=rep(0, 60))
B$Count[match(A$Quality, B$Quality)] <- A$Count

plot(NULL, xlim=c(1,41), ylim=range(B$Count), xlab="Mean Sequence Quality (Phred Score)", ylab="Counts", main="Per base sequence quality")
for(k in 1:nrow(B)){
if(B$Count[k]>0){
rect(k-0.5, 0, k+0.5, B$Count[k], border=F, col=if(B$Quality[k]<30){"Red"}else{"Green"})
}
}
Q30 <- sum(B$Count[B$Quality>=30])/sum(B$Count)*100
text(30, max(B$Count)*0.8, paste("%<=Q30\n", round(Q30, 2)), adj=0)

exp$Q30[i] <- round(Q30, 2)

#qc_plot(QClist[[i]], "Sequence length distribution")
A <- data.frame(QClist[[i]]$sequence_length_distribution)
if(nrow(A)==1){A <- rbind(c(A$Length-1, 0), A, c(A$Length+1, 0))}
plot(1:nrow(A), A$Count, xlim=c(1, nrow(A)), ylim=c(0, max(A$Count)), typ="l", main="Sequence length distribution", xlab="Sequence length in bp", ylab="Count")

exp[i,3] <- A$Length[A$Count==max(A$Count)]

exp[i,4] <- round(max(A$Count)/sum(A$Count)*100, 2) # percent of sequences with length

#qc_plot(QClist[[i]], "Per base sequence content")
A <- data.frame(QClist[[i]]$per_base_sequence_content)
plot(NULL, ylim=c(0, 100), xlim=c(1, nrow(A)), xlab="", ylab="Counts", main="Per base sequence content", xaxt="n")
axis(1, labels=A$Base[c(1, 5, seq(10, 60, 10))], at=c(1, 5, seq(10, 60, 10)))
polygon(c(1, 1:nrow(A), nrow(A)), c(0,A$G+A$A+A$T+A$C,0), border=F, col="Blue")
polygon(c(1, 1:nrow(A), nrow(A)), c(0,A$G+A$A+A$T,0), border=F, col="Green")
polygon(c(1, 1:nrow(A), nrow(A)), c(0,A$G+A$A,0), border=F, col="Red")
polygon(c(1, 1:nrow(A), nrow(A)), c(0,A$G,0), border=F, col="Yellow")

SD <- NULL
for (m in 1:nrow(A)){
SD[m] <- sd(abs(A[m,-1]-25))
}

exp[i,7] <- round(mean(SD), 2)



#qc_plot(QClist[[i]], "Adapter content")
A <- data.frame(QClist[[i]]$adapter_content)
plot(NULL, ylim=c(1, 100), xlim=c(1, nrow(A)), xlab="Possition in Read (bp)", ylab="% Adapter", main="Adapter content", xaxt="n")
select <- c(1, 5, seq(10, 60, 10))
axis(1, labels=A$Position[select], at= select)

for (k in 2:ncol(A)){
lines(1:nrow(A), A[,k], col="black", lwd=1)
}

temp <- rowSums(A[,-1])
exp[i,5] <- round(max(temp), 2)


#qc_plot(QClist[[i]], "Per base N content")
A <- data.frame(QClist[[i]]$per_base_n_content)
plot(NULL, ylim=c(1, 100), xlim=c(1, nrow(A)), xlab="Possition in Read (bp)", ylab="% N", main="Per base N content", xaxt="n")
select <- c(1, 5, seq(10, 60, 10))
axis(1, labels=A$Base[select], at= select)

lines(1:nrow(A), A[,2], col="Red", lwd=1.5)

exp[i,6] <- round(max(A$N.Count), 2)

qc_plot(QClist[[i]], "Per sequence GC content")
qc_plot(QClist[[i]], "Sequence duplication levels")
#qc_plot(QClist[[i]], "Overrepresented sequences")
#qc_plot(QClist[[i]], "Kmer content")

dev.off()

}


# aggregate
write.csv(exp, paste(folder, "/FastQC/stats.csv", sep=""))

exp$SequDepth <- as.numeric(exp$SequDepth)
exp$GC <- as.numeric(exp$GC)

space <- 5
maxdepth <- max(exp$SequDepth)

pdf(paste(folder, "/FastQC/plot_tables.pdf", sep=""), width=8, height=12)

g <- 1
while(nrow(exp)>(g*50-49)){

data <- exp[(g*50-49):(g*50),]
data <- data[!is.na(data$ID),]

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(NULL, xlim=c(c(-space-0.1),16), ylim=c(50, 0), xlab="", ylab="", xaxt="n", yaxt="n", axes=F)

text(-space, 1:nrow(data), data$ID, adj=0) # ID

for (i in 1:nrow(data)){
rect(-0.5, i-0.5, 10, i+0.5, col="white", border=F)
rect(0, i-0.5, data$SequDepth[i]/maxdepth*2, i+0.5, border=F, col="Gray") # depth
rect(2, i-0.5, 2+data$LenghtRel[i]/50, i+0.5, border=F, col="Gray") # depth
rect(4, i-0.5, 4+data$Q30[i]/50, i+0.5, border=F, col=if(data$Q30[i]>=80){"lightgreen"} else {"Red"}) # depth
rect(6, i-0.5, 6+data$MaxAdapters[i]/50, i+0.5, border=F, col=if(data$MaxAdapters[i]<5){"lightgreen"} else {"Red"}) # adapters
rect(8, i-0.5, 8+data$MaxNRel[i]/50, i+0.5, border=F, col=if(data$MaxNRel[i]<1){"lightgreen"} else {"Red"}) # depth
rect(10, i-0.5, 10+data$Diversity[i]/50, i+0.5, border=F, col=if(data$Diversity[i]<5){"lightgreen"} else {"Red"}) # depth
rect(12, i-0.5, 12+data$GC[i]/50, i+0.5, border=F, col="Gray") # GC
rect(14, i-0.5, 14+data$duplicated[i]/50, i+0.5, border=F, col="Gray") # GC

}

text(0.1, 1:nrow(data), paste(round(data$SequDepth/1000), "K", sep=""), adj=0) # ID
text(2.1, 1:nrow(data), data$Length, adj=0) # lengt
text(4.1, 1:nrow(data), data$Q30, adj=0) # Q30
text(6.1, 1:nrow(data), data$MaxAdapters, adj=0) # adapters
text(8.1, 1:nrow(data), data$MaxNRel, adj=0) # Ns
text(10.1, 1:nrow(data), data$Diversity, adj=0) # base composition - mean SD from 25%
text(12.1, 1:nrow(data), data$GC, adj=0) # GC
text(14.1, 1:nrow(data), data$duplicated, adj=0) # GC

for(m in c(0,2,4,6,8,10,12, 14)){
rect(m, nrow(data)+0.5,  m+2, 0.5)
}
rect(-space, nrow(data)+0.5,  0, 0.5)

text(-space, 0, "ID", adj=0)
text(0.1, 0, "Depth", adj=0)
text(2.1, 0, "Length", adj=0)
text(4.1, 0, "Q30+", adj=0)
text(6.1, 0, "Adapt", adj=0)
text(8.1, 0, "N", adj=0)
text(10.1, 0, "DivSD", adj=0)
text(12.1, 0, "GC", adj=0)
text(14.1, 0, "Dumpli", adj=0)

g <- g+1
#print(g)
}

dev.off()

message(" ")
message(" Module completed!")

#cat(file="log.txt", paste(Sys.time(), "\n", "*** Module completed!", "", sep="\n"), append=T, sep="\n")


}

