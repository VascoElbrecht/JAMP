# U_merge_PE v0.1
# maybe add option to split and merge large files automatically?

U_merge_PE <- function(files="latest", file1=NA, file2=NA, fastq_maxdiffs=99, fastq_pctid=75, fastq_minovlen=16, fastq=T, LDist=T){

folder <- Core(module="U_merge_PE")
cat(file="log.txt", c("\n","Version v0.2", "\n"), append=T, sep="\n")
message(" ")


if(!is.na(file1[1])&!is.na(file2[1])){
files[1] <- "NOT USED"

temp <- "Custom list of files provided as file1 and file2 will be used\n\n"

message(temp)
cat(file="log.txt", temp, append=T, sep="\n")
} else {


if (files[1]=="latest"){
source(paste(folder, "/robots.txt", sep=""))
file1 <- list.files(paste(last_data, "/_data", sep=""), full.names=T, pattern="_[rR]1")
file2 <- list.files(paste(last_data, "/_data", sep=""), full.names=T, pattern="_[rR]2")} else {

if (length(files)>1 ){
file1 <- files[grep("_[rR]1", files)]
file2 <- files[grep("_[rR]2", files)]
}
}

merge_identical <- sub(".*_data/(.*)_[rR]1.*", "\\1", file1)==sub(".*_data/(.*)_[rR]2.*", "\\1",file2)
# merging not identical reads
if(!sum(merge_identical)==length(merge_identical)){
warning("There is a problem with the files you want to merge. Not all fastq files have a matchign pair with identical name. Please check. Package stopped.")
setwd("../")
stop()
}

}






if(length(grep(".*N_debris_r1..*", file1))==1){message("N_debris are excluded and not merged.")}

file1 <- file1[!grepl(".*N_debris_r1..*", file1)] # remove debres from list
file2 <- file2[!grepl(".*N_debris_r2..*", file2)] # remove debres from list


message(paste("Starting to PE merge ", length(file1), " samples.", sep=""))
message(" ")

# new file names

new_names <- sub(".*(/.*)", "\\1", file1)
if(fastq){new_names <- sub("[rR]1.*", "PE.fastq", new_names)} else {new_names <- sub("[rR]1.*", "PE.fasta", new_names)}

new_names <- paste(folder, "/_data", new_names, sep="")


dir.create(paste(folder, "/_stats/merge_stats", sep=""))
log_names <- sub("_data", "_stats/merge_stats", new_names)
log_names <- sub("_PE.fast[aq]", "_PE_log.txt", log_names)

cmd <- paste(" -fastq_mergepairs \"", file1, "\" -reverse \"", file2,  "\" ", if(fastq){"-fastqout"} else {"-fastaout"}, " \"", new_names, "\"", " -report ", log_names, " -fastq_maxdiffs ", fastq_maxdiffs , " -fastq_pctid ", fastq_pctid, " -fastq_trunctail 0 -fastq_minovlen ", fastq_minovlen, sep="")

tab_exp <- NULL
for (i in 1:length(cmd)){
system2("usearch", cmd[i], stdout=F, stderr=F)
temp <- readLines(log_names[i])

# save cmd in log name!
cat(file=log_names[i], c(paste("usearch", cmd[i]), temp), sep="\n")


# table export
merged <- as.numeric(sub(".*merged \\((.*)..", "\\1",temp[6]))
rawdata_counts <- as.numeric(sub(".* / (.*) pairs merged.*", "\\1", temp[6]))
median_length <- as.numeric(sub("(.*)Median", "\\1",temp[11]))
temp_count <- Count_sequences(new_names[i], fastq)
short_name <- sub(".*_data/(.*)_PE.fast.", "\\1", new_names[i])
tab_exp <- rbind(tab_exp, c(short_name, rawdata_counts, temp_count, merged, median_length))

meep <- paste(short_name, ": ", merged, "% merged - median length: ", median_length, sep="")
message(meep)
cat(file="log.txt", meep, append=T, sep="\n")
}

cat(file="log.txt", "\n", append=T, sep="\n")


tab_exp <- data.frame(tab_exp)
names(tab_exp) <- c("Sample", "Sequ_count_in", "Sequ_count_out", "percent_merged", "median_length")

write.csv(tab_exp, paste(folder, "/_stats/sequ_length_abund.csv", sep=""))



# make some plots
temp <- read.csv(paste(folder, "/_stats/sequ_length_abund.csv", sep=""), stringsAsFactors=F)

Sequences_lost(temp$Sequ_count_in, temp$Sequ_count_out, temp$Sample, out=paste(folder, "/_stats/Sequences_merged.pdf", sep=""), main=paste(folder, ": Proportion of reads merged", sep=""))
Sequences_lost(temp$Sequ_count_in, temp$Sequ_count_out, temp$Sample, rel=T, out=paste(folder, "/_stats/Sequences_merged_rel.pdf", sep=""), main=paste(folder, ": Proportion of reads merged", sep=""))

merged_message <- paste("\nOn average ", round(mean(temp$percent_merged), 2), "% sequences merged (SD = ", round(sd(temp$percent_merged), 2), "%).\n", sep="")
message(merged_message)
cat(file="log.txt", merged_message, append=T, sep="\n")


#make length distribution plots
if(LDist){

dir.create(paste(folder, "_stats/length distribution", sep="/"))

message("Generating length distribution plots. If this takes to long you can turn this option off with setting \"LDist=F\".")

for (i in 1:length(new_names)){

pdfname <- sub("/_data/", "/_stats/length distribution/", new_names[i])
pdfname <- sub(".fast.", ".pdf", pdfname)

message(paste("Plotting ", sub(".*distribution/(.*)_PE_.pdf","\\1", pdfname), sep=""))
Length_distribution(new_names[i], pdfname)
}
message(" ")
}# Ldist end


message("Done with PE merging!")

cat(file="log.txt", paste(Sys.time(), "Done with PE merging", "", "*** Module completed!\n\n", sep="\n"), append=T, sep="\n")

}

