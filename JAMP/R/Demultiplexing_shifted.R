# Demultiplexing_shifted v0.1

Demultiplexing_shifted <- function(file1, file2, tags=NA, combinations=NA, md5=T, OS="autodetect"){

Core(module="Demultiplexing_shifted")
cat(file="../log.txt", c("\n", "Version v0.1", "\n"), append=T, sep="\n")


# get file compression
compression <- gsub(".*\\.(.*)$", "\\1", file1)


# get tagging table
if(tags=="BF_BR"){
barcodes <- read.csv(paste(system.file(package="JAMP"), "/BF_BR.csv", sep=""), stringsAsFactors=F)
} else {barcodes <- read.csv(tags, stringsAsFactors=F)}

combos <- read.csv(combinations, stringsAsFactors=F)

#basic log stats
temp <- paste(Sys.time(), "Starting demultiplexing using:", paste("Barcode table:", tags, sep=""), paste("Searching for ", nrow(combos)/2, " samples as given in table: ", combinations, sep=""), paste("Raw data format: ", paste(unique(compression), collapse=" "), sep=""), paste("Read 1 RAW data: ", paste(file1, collapse=" "), sep=""), paste("Read 2 RAW data: ", paste(file2, collapse=" "), sep=""), "\n", sep="\n")
cat(file="../log.txt", temp, append=T, sep="\n")


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
temp <- c(file1, file2)

A <- NULL
for (i in 1:length(temp)){
A[i] <- system2(md5_cmd, temp[i], stdout=T)
}
cat(file="../log.txt", A, append=T, sep="\n")
}


savesequ <- function(file_name1="N_debres_r1.txt", file_name2="N_debres_r2.txt", rm1=0, rm2=0){
if(length(A)>0){

rm1 <- rm1+1 # remove nucleotides
rm2 <- rm2+1
exp_sequ1 <- data1[sort(c(A*4, A*4-1, A*4-2, A*4-3))]
readL <- nchar(exp_sequ1[seq(2, length(exp_sequ1), 4)])
exp_sequ1[seq(2, length(exp_sequ1), 4)] <- substr(exp_sequ1[seq(2, length(exp_sequ1), 4)], rm1, readL)
exp_sequ1[seq(4, length(exp_sequ1), 4)] <- substr(exp_sequ1[seq(4, length(exp_sequ1), 4)], rm1, readL)

exp_sequ2 <- data2[sort(c(A*4, A*4-1, A*4-2, A*4-3))]
readL <- nchar(exp_sequ2[seq(2, length(exp_sequ2), 4)])
exp_sequ2[seq(2, length(exp_sequ2), 4)] <- substr(exp_sequ2[seq(2, length(exp_sequ2), 4)], rm2, readL)
exp_sequ2[seq(4, length(exp_sequ2), 4)] <- substr(exp_sequ2[seq(4, length(exp_sequ2), 4)], rm2, readL)

cat(exp_sequ1, file= paste("_data/", file_name1, sep=""), append=T, sep="\n")
cat(exp_sequ2, file= paste("_data/", file_name2, sep=""), append=T, sep="\n")}
}
# save function end

i <- 1
chunk <- 100000*4
mytime <- Sys.time()

for (i in 1:length(file1)){

if(compression[i]=="gz"){
con1 <- gzfile(file1[i], "rt")
con2 <- gzfile(file2[i], "rt")
} else if(compression[i]=="bz2"){
con1 <- bzfile(file1[i], "rt")
con2 <- bzfile(file2[i], "rt")
} else if(compression[i]=="fastq"){
con1 <- file(file1[i], "rt")
con2 <- file(file2[i], "rt")
} else {
warning("Raw data has to be compressed in \".gz\" or \".bz2\" or uncompressed in \".fastq\" format! Processing was stopped, please check your raw data format.")
setwd("../")
stop()
}



repeat {
data1 <- readLines(con1, chunk) # read chunck
data2 <- readLines(con2, chunk) # read chunck

data1 <- data1[!is.na(data1)] # remove NA in debres if end of file is reached!
data2 <- data2[!is.na(data2)]

if (!length(data1)) break

sequ1 <- data1[seq(2, length(data1), 4)]
sequ2 <- data2[seq(2, length(data2), 4)]

sequ1 <- substr(sequ1, 1, 5)
sequ2 <- substr(sequ2, 1, 5)

match_sequ1 <- match(sequ1, barcodes$barcode, nomatch=0)
match_sequ2 <- match(sequ2, barcodes$barcode, nomatch=0)

barcodeID <- c("NA", barcodes$ID)

tags <- paste(barcodeID[match_sequ1+1], barcodeID[match_sequ2+1], sep="_")
match_sequ <- match(tags, combos$ID)


# write down number of bases to remove
match_sequ1[match_sequ1==0] <- NA
rm1_master <- barcodes$rm[match_sequ1]

match_sequ2[match_sequ2==0] <- NA
rm2_master <- barcodes$rm[match_sequ2]

for (k in 1:nrow(combos)){ # save all xxx tag combinations based on demultifile (name file)
A <- which(match_sequ==k)
if(length(A)>0){savesequ(combos$File1[k], combos$File2[k], rm1_master[A], rm2_master[A])}
}
A <- which(is.na(match_sequ))
if(length(A)>0){savesequ("N_debris_r1.txt", "N_debris_r2.txt", 0, 0)}

# repeat loop enede
} # looping all sequence files!
close(con1)
close(con2)
}

temp <- paste(Sys.time(), paste("Done in:", message(mytime - Sys.time()), collapse=""), sep="\n")
cat(file="../log.txt", temp, append=T, sep="\n")


message("done in:")
message(mytime - Sys.time())
message(" ")

# count number of reads in each file
abundance <- Count_sequences(list.files("_data", full.names=T))
temp <- data.frame("files_demultiplexed"=list.files("_data"), abundance)

write.csv(temp, "_stats/read_abundance.csv")

temp <- temp[seq(1, length(abundance), 2),]

options("scipen"=100, "digits"=7)
message(paste("A total of ", sum(temp$abundance), "sequences where demultiplexed"), sep="")
message(paste(round(temp$abundance[temp=="N_debris_r1.txt"]/sum(temp$abundance)*100, digits=2)), "% of sequences could not be matched with any of the tagging combinations (e.g. sequencing errors in the tags or PhiX.)", sep="")

temp <- paste(Sys.time(), paste("Number of reads demultiplexed ", sum(temp$abundance), sep=""), paste("Number of not matching reads (N_debris):", temp$abundance[temp=="N_debris_r1.txt"], sep=""), paste("Relative abundance of not matching reads: ", round(c(temp$abundance[temp=="N_debris_r1.txt"]/sum(temp$abundance)*100), digits=2), "\n", "\n", "Module completed!", "\n", "\n", sep=""), sep="\n")
cat(file="../log.txt", temp, append=T, sep="\n")



# end count reads
#options("scipen"=-100, "digits"=7)
setwd("../") # return to base folder
}


?options
