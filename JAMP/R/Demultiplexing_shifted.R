# Demultiplexing_shifted v0.1

Demultiplexing_shifted <- function(file1, file2, tags=NA, combinations=NA){

Core(module="Demultiplexing_shifted")

# get file compression
compression <- gsub(".*\\.(.*)$", "\\1", file1)


# get tagging table
if(tags=="BF_BR"){
barcodes <- read.csv(paste(system.file(package="JAMP"), "/BF_BR.csv", sep=""), stringsAsFactors=F)
} else {barcodes <- read.csv(tags, stringsAsFactors=F)}

combos <- read.csv(combinations, stringsAsFactors=F)


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

for (i in 1:length(file1)){

if(compression=="gz"){
con1 <- gzfile(file1[i], "rt")
con2 <- gzfile(file2[i], "rt")
} else if(compression=="bz2"){
con1 <- bzfile(file1[i], "rt")
con2 <- bzfile(file2[i], "rt")
} else if(compression=="fastq"){
con1 <- file(ngsfile1, "rt")
con2 <- file(ngsfile2, "rt")
} else {
warning("Raw data has to be compressed in \".gz\" or \".bz2\" or uncompressed in \".fastq\" format! Processing was stopped, please check your raw data format.")
setwd("../")
stop()
}

mytime <- Sys.time()
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
rm2_master <- barcodes$rm[match_sequ1]


for (k in 1:nrow(combos)){ # save all xxx tag combinations based on demultifile (name file)
A <- which(match_sequ==k)
if(length(A)>0){savesequ(combos$File1[k], combos$File2[k], rm1_master[A], rm1_master[A])}
}
A <- which(is.na(match_sequ))
if(length(A)>0){savesequ("N_debris_r1.txt", "N_debris_r2.txt", 0, 0)}

# repeat loop enede
} # looping all sequence files!

print("done in:")
print(mytime - Sys.time())
print("")
close(con1)
close(con2)

}


setwd("../") # return to base folder
}
