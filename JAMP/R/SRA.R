# SRA downloader v0.1

SRA <- function(ID=NA, rename=NA, split3=T, exe="fastq-dump", delete_data=T){

folder <- Core(module="SRA", delete_data=delete_data)
cat(file="log.txt", c("\n","Version v0.1", "\n"), append=T, sep="\n")
message(" ")

files_to_delete <- NULL


temp <- paste("Sarting to download " , length(ID), " samples from NCBI SRA!", sep="")
message(temp)
message(" ")
cat(file="log.txt", temp, append=T, sep="\n")


# cmd to DL from SRA
cmd <- paste(if(split3){"--split-3 "}, ID, " -O ", folder,"/_data/", sep="")


tab_exp <- NULL
for (i in 1:length(cmd)){
A <- system2(exe, cmd[i], stdout=T, stderr=T)
cat(c(paste(exe, cmd[i]), "", A, "", ""), file=paste(folder, "/_stats/download_log.txt", sep=""), sep="\n", append=T)

tab_exp[i] <- as.numeric(sub("Written (.*) spots.*", "\\1", A[2]))
temp <- paste(ID[i], ": ", tab_exp[i], if(!is.na(rename)){paste(" -> ", rename[i], sep="")}, sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

}
temp <- paste("Downloading complete!", " Downloaded a toatal of ", sum(tab_exp), " spots from ", length(ID), " IDs.", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")


# rename downloaded files!
if(!is.na(rename[1])){
if(length(rename)==length(ID)){

filelist <- list.files(paste(folder, "/_data", sep=""), full.names=T)

renamed <- filelist
for (i in 1:length(ID)){
renamed <- sub(ID[i], rename[i], renamed)
}

cbind(filelist, renamed)

A <- file.rename(filelist, renamed)
temp <- paste(sum(A), " out of ", length(A), " files renamed! (read 1 and 2)", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

} else {
stop(paste("Number of IDs has the match the number of new names! IDs: ", length(ID), " New Names: ", length(ID), sep=""))
}

# Add files to delete - after renaming!
files_to_delete <- c(files_to_delete, renamed)

} else { # add not renamed files to delete list
files_to_delete <- c(files_to_delete, list.files(paste(folder, "/_data", sep=""), full.names=T))
}


# write stats file!
if(!is.na(rename[1])){
final_names <- sub(".*_data/(.*)_1.fastq", "\\1", renamed[seq(1, length(renamed), 2)])} else {
final_names <- ID
}


tab_exp <- data.frame("SRA_ID"=ID, "spots"= tab_exp, "file_names"=final_names)

write.csv(tab_exp, paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_SRA_download_stats.csv", sep=""))

# make some plots

temp <- read.csv(paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_SRA_download_stats.csv", sep=""), stringsAsFactors=F)

Sequences_lost(temp$spots, temp$spots, temp$file_names, out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_spots_downladed.pdf", sep=""), main=paste(folder, ": Downloaded spots ", sep=""))



message(" ")
message("Module completed!")

if(delete_data){
cat(file=paste(folder, "/robots.txt", sep=""), "\n# DELETE_START", files_to_delete, "# DELETE_END", append=T, sep="\n")
}

cat(file="log.txt", paste(Sys.time(), "*** Module completed!", "", sep="\n"), append=T, sep="\n")
}

