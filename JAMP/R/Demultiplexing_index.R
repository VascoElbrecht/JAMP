# Demultiplexing_indexed
# v0.1

Demultiplexing_index <- function(files="latest", fileR1=NA, fileR2=NA, fileI1=NA, indexTable=NA, software="illumina-utils", exe="iu-demultiplex", revcompI=T, md5=T, OS="autodetect", delete_data=T){

folder <- Core(module="Demultiplexing_index", delete_data=delete_data)
cat(file="log.txt", c("\n", "Version v0.1", "\n"), append=T, sep="\n")

files_to_delete <- NULL

indextab <- read.csv(indexTable, sep="\t", stringsAsFactors=F)

#auto detect files!
if (is.na(fileR1[1])){
source(paste(folder, "/robots.txt", sep=""))
files <- list.files(paste(last_data, "/_data", sep=""), full.names=T)

fileR1 <- files[grep("R1.*_00.\\.fastq\\.*", files)]
fileR2 <- files[grep("R2.*_00.\\.fastq\\.*", files)]
fileI1 <- files[grep("I1.*_00.\\.fastq\\.*", files)]

temp <- paste("Files autoetected!\nfileR1:\n", fileR1, "\nfileR2:\n", fileR2, "\nfileI1\n", fileI1, "\n", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")
}



# check md5
temp <- "Calculating md5 checksums to verify file integrety:"
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

if(md5){

if(OS=="autodetect"){
sys <- Sys.info()[['sysname']]
if(sys=="Darwin"){
sys <- "Mac"
md5_cmd <- "md5"
}
if(sys=="Linux"){
md5_cmd <- "md5sum"
}
} else {
if(OS=="Mac"){md5_cmd <- "md5"} else {md5_cmd <- "md5sum"}
}


# claculate md5 checksums for log
temp <- c(fileR1, fileR2, fileI1)

A <- NULL
for (i in 1:length(temp)){
A[i] <- system2(md5_cmd, paste("\"", temp[i], "\"", sep=""), stdout=T)
}
cat(file="log.txt", A, append=T, sep="\n")
message(paste(A, collapse="\n"))
}




# ungz files if needed!
if(software=="illumina-utils"){
if(sum(grep("\\.gz$", fileR1))>0){ # R1
system2("gunzip", paste(" -k", fileR1[grep("\\.gz$", fileR1)], sep=" "))
fileR1 <- sub("\\.gz$", "", fileR1)
files_to_delete <- c(files_to_delete, fileR1)
}
if(sum(grep("\\.gz$", fileR2))>0){ # R2
system2("gunzip", paste(" -k", fileR2[grep("\\.gz$", fileR2)], sep=" "))
fileR2 <- sub("\\.gz$", "", fileR2)
files_to_delete <- c(files_to_delete, fileR2)
}
if(sum(grep("\\.gz$", fileI1))>0){ # R2
system2("gunzip", paste(" -k", fileI1[grep("\\.gz$", fileI1)], sep=" "))
fileI1 <- sub("\\.gz$", "", fileI1)
files_to_delete <- c(files_to_delete, fileI1)
}
}


# demultiplex

if(software=="illumina-utils"){
#needs to be moved in to bin!
temp <- paste("starting to demultiplex using \"", software, "\" and indexes from the file \"", indexTable, "\"", sep="")
message(temp)
if(md5){cat(file="log.txt", A, append=T, sep="\n")}


system2(exe, paste("-s ", indexTable, " --r1 ", fileR1, " --r2 ", fileR2, " -i ", fileI1, " -o ", folder, "/_data/", if(revcompI){" -x"}, sep=""))


# move stats file and make plots of sequencing depth!
# also add summary stats

file.copy(paste(folder, "/_data/00_DEMULTIPLEXING_REPORT", sep=""), paste(folder, "/_stats/00_DEMULTIPLEXING_REPORT.txt", sep=""))
file.remove(paste(folder, "/_data/00_DEMULTIPLEXING_REPORT", sep=""))

# make squencing depth plot
tab <- read.csv(paste(folder, "/_stats/00_DEMULTIPLEXING_REPORT.txt", sep=""), sep="\t", stringsAsFactors=F)

tab <- tab[order(tab$sample), ]

rel <- tab$num_reads_stored/sum(tab$num_reads_stored)*100
tab <- data.frame(tab, "rel"=rel)


write.csv(tab, paste(folder, "/_stats/stats.csv", sep=""))


# make differnt pages
pdf(paste(folder, "/_stats/plot_sequ_depth.pdf", sep=""), width=8, height=12)
g <- 1
while(nrow(tab)>(g*50-49)){

data <- tab[(g*50-49):(g*50),]
data <- data[!is.na(data$sample),]

par(mar=c(0,0,2,0), oma=c(0,0,0,0))
plot(NULL, xlim=c(-max(rel)*0.48, max(rel)), ylim=c(50, 1), xlab="Number of sequences per sample", ylab="", yaxt="n", xaxt="n", axes=F)
text(-max(rel)*0.5, 1:nrow(data), data$sample, adj=0) # ID
axis(3)

for (i in 1:50){
rect(0, i-0.5, data$rel[i], i+0.5, col="Gray", border=F)
}
text(max(rel)*0.01, 1:nrow(data), data$num_reads_stored, adj=0) # ID

g <- g+1
}

dev.off()

# report number of sequences in files!
ReadsIn <- Count_sequences(fileI1)
temp <- paste("A total of ", ReadsIn, " were processed, of wich ", sum(tab$num_reads_stored), " reads could be assigned to the provided indexes (", round(sum(tab$num_reads_stored)/ReadsIn*100, 2), "%).\nDetected samples: ", nrow(tab), "/", nrow(indextab), " (", round(nrow(tab)/nrow(indextab)*100, 2), "%)", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

if(!nrow(tab)/nrow(indextab)*100==100){
temp <- paste("WARNING: Number of reads in index files does not match the number of samples recovered. Verify that files are correct!")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

}

renamME <- list.files(paste(folder, "/_data/", sep=""), full.names=T)
file.rename(renamME, gsub("-", "_", renamME))

files_to_delete <- c(files_to_delete, gsub("-", "_", renamME))

} # iu-demultiplex

cat(file=paste(folder, "/robots.txt", sep=""), "\n# DELETE_START", files_to_delete, "# DELETE_END", append=T, sep="\n")

message("\nModule completed!\n\n")

cat(file="log.txt", "*** Module completed!\n\n", append=T, sep="\n")

}


