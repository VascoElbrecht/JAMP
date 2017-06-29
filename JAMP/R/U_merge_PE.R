# U_merge_PE v0.1
# maybe add option to split and merge large files automatically?

U_merge_PE <- function(files="latest", file1=NA, file2=NA, fastq_maxdiffs=99, fastq_pctid=90, fastq=T){

Core(module="U_merge_PE")
cat(file="../log.txt", c("\n","Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files[1]=="latest"){
source("robots.txt")
file1 <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T, pattern="_r1.txt")
file2 <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T, pattern="_r2.txt")
} else {
file1 <- sub(".*(_data/.*)", "../\\1", file2)
file2 <- sub(".*(_data/.*)", "../\\1", file2)
}

merge_identical <- sub(".*_data/(.*)_r1.txt", "\\1", file1)==sub(".*_data/(.*)_r2.txt", "\\1",file2)
# merging not identical reads
if(!sum(merge_identical)==length(merge_identical)){
warning("There is a problem with the files you want to merge. Not all fastq files have a matchign pair with identical name. Please check. Package stopped.")
setwd("../")
stop()
}

if(length(grep(".*N_debris_r1.txt", file1))==1){message("N_debris are excluded and not merged.")}

file1 <- file1[!grepl(".*N_debris_r1.txt", file1)] # remove debres from list
file2 <- file2[!grepl(".*N_debris_r2.txt", file2)] # remove debres from list


message(paste("Starting to PE merge ", length(file1), " samples.", sep=""))
message(" ")

# new file names

new_names <- sub(".*(_data/.*)", "\\1", file1)
if(fastq){new_names <- sub("r1.txt", "PE.fastq", new_names)} else {new_names <- sub("r1.txt", "PE.fasta", new_names)}


dir.create("_stats/merge_stats")
log_names <- sub("_data", "_stats/merge_stats", new_names)
log_names <- sub("_PE.fast[aq]", "_PE_log.txt", log_names)

cmd <- paste(" -fastq_mergepairs \"", file1, "\" -reverse \"", file2,  "\" ", if(fastq){"-fastqout"} else {"-fastaout"}, " \"", new_names, "\"", " -report ", log_names, " -fastq_maxdiffs ", fastq_maxdiffs , " -fastq_pctid ", fastq_pctid ," -fastq_trunctail 0", sep="")

tab_exp <- NULL
for (i in 1:length(cmd)){
system2("usearch", cmd[i], stdout=F, stderr=F)
temp <- readLines(log_names[i])

# table export
merged <- as.numeric(sub(".*merged \\((.*)..", "\\1",temp[6]))
median_length <- as.numeric(sub("(.*)Median", "\\1",temp[11]))
temp_count <- Count_sequences(new_names[i], fastq)
short_name <- sub("_data/(.*)_PE.fast.", "\\1", new_names[i])
tab_exp <- rbind(tab_exp, c(short_name, temp_count, merged, median_length))

meep <- paste(short_name, ": ", merged, "% merged - median length: ", median_length, sep="")
message(meep)
cat(file="../log.txt", meep, append=T, sep="\n")
}

cat(file="../log.txt", "\n", append=T, sep="\n")


tab_exp <- data.frame(tab_exp)
names(tab_exp) <- c("Sample", "Sequ_count", "percent_merged", "median_length")

write.csv(tab_exp, "_stats/sequ_length_abund.csv")

# make some plots?

message(" ")
message(" Done with PE merging")

cat(file="../log.txt", paste(Sys.time(), "Done with PE merging", "", sep="\n"), append=T, sep="\n")

setwd("../")
}

