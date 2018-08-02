# Demultiplexing_shifted

Demultiplexing_index <- function(files="latest", fileR1=NA, fileR2=NA, fileI1=NA, indexTable=NA, software="illumina-utils", revcompI=T md5=T, OS="autodetect"){

folder <- Core(module="Demultiplexing_index")
cat(file="log.txt", c("\n", "Version v0.1", "\n"), append=T, sep="\n")


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
}
if(sum(grep("\\.gz$", fileR2))>0){ # R2
system2("gunzip", paste(" -k", fileR2[grep("\\.gz$", fileR2)], sep=" "))
fileR2 <- sub("\\.gz$", "", fileR2)
}
if(sum(grep("\\.gz$", fileI1))>0){ # R2
system2("gunzip", paste(" -k", fileI1[grep("\\.gz$", fileI1)], sep=" "))
fileI1 <- sub("\\.gz$", "", fileI1)
}
}


# demultiplex

if(software=="illumina-utils"){
#needs to be moved in to bin!
temp <- paste("starting to demultiplex using \"", software, "\" and indexes from the file \"", indexTable, "\"", sep="")
message(temp)
cat(file="log.txt", A, append=T, sep="\n")


A <- system2("iu-demultiplex", paste("-s ", indexTable, " --r1 ", fileR1, " --r2 ", fileR1, " -i ", fileI1, " -o ", folder, "/_data/", if(revcompI){" -x"}, sep=""))




cat(file="log.txt", "*** Module completed!\n\n", append=T, sep="\n")

}


