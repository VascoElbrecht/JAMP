# U_max_ee v0.1

U_max_ee <- function(files="latest", max_ee=1, fastq=F){

folder <- Core(module="U_max_ee")
cat(file="log.txt", c("\n","Version v0.2", "\n"), append=T, sep="\n")
message(" ")

if (files[1]=="latest"){
source(paste(folder, "/robots.txt", sep=""))
files <- list.files(paste(last_data, "/_data", sep=""), full.names=T)
}

temp <- paste("Starting to quality filter (max expected errors = " , max_ee, ") in ", length(files), " samples.", sep="")
message(temp)
message(" ")
cat(file="log.txt", temp, append=T, sep="\n")

# new file names
new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- paste(folder, "/", new_names, sep="")


if(fastq){
new_names <- sub(".fastq", "_ee.fastq", new_names) # keep fastq
} else {
new_names <- sub(".fastq", "_ee.fasta", new_names) # convert to fasta
}
new_names <- sub("_ee\\.", paste("_ee", max_ee, ".", sep=""), new_names)



dir.create(paste(folder, "/_stats/merge_stats", sep=""))
log_names <- sub("_data", "_stats/merge_stats", new_names)
log_names <- sub("_ee.fast.", "_ee.txt", log_names)

# cmd max EE
cmd <- paste("-fastq_filter \"", files, "\"", if(fastq){" -fastqout "} else {" -fastaout "}, "\"", new_names, "\" -fastq_maxee ", max_ee, " -fastq_qmax 60", sep="")

i <- 11
tab_exp <- NULL
for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
cat(A, file=log_names[i], sep="\n")

imput <- A[grep("Reads \\(", A)]
imput <- sub("(.*)Reads.*", "\\1", imput)
imput <- as.numeric(imput)

pass <- A[grep("Filtered reads", A)]
pass <- sub("(.*)Filtered reads .*", "\\1", pass)
pass <- as.numeric(pass)

pct <- round(pass/imput*100, digits=2)

# deal with empty files
if(length(imput)==0){
imput <- 0
pct <- 0
}
if(length(pass)==0){ # can still contain files but 0 passed
pass <- 0
pct <- 0
}



short_name <- sub(".*_data/(.*)_PE.*", "\\1", new_names[i])
tab_exp <- rbind(tab_exp, c(short_name, imput, pass, pct))

if(pass>0){
meep <- paste(short_name, ": ", pct, "% pass EE ", sep="")} else {
meep <- paste("WARNING:", short_name, ": ", if(imput>0){"0 sequences passed quality filtering!"} else {"Imput file did not contain sequences!"}, sep="")
}

message(meep)
cat(file="log.txt", meep, append=T, sep="\n")
}
cat(file="log.txt", "\n", append=T, sep="\n")

tab_exp <- data.frame(tab_exp)
names(tab_exp) <- c("Sample", "Sequ_count_in", "Sequ_count_out", "pct_pass")

write.csv(tab_exp, paste(folder, "/_stats/max_ee_stats.csv", sep=""))

# make some plots

temp <- read.csv(paste(folder, "/_stats/max_ee_stats.csv", sep=""), stringsAsFactors=F)

Sequences_lost(temp$Sequ_count_in, temp$Sequ_count_out, sub("_PE.*", "", temp$Sample), out=paste(folder, "/_stats/LowQuality_sequences.pdf", sep=""), main=paste(folder, ": max expected errors ", max_ee, sep=""))
Sequences_lost(Reads_in=temp$Sequ_count_in, Reads_out=temp$Sequ_count_out, Sample_names=sub("_PE.*", "", temp$Sample), rel=T, out=paste(folder, "/_stats/LowQuality_sequences_rel.pdf", sep=""), main=paste(folder, ": max expected errors ", max_ee, sep=""))

merged_message <- paste("\nSequences with sufficient read quality: ", round(mean(temp$pct_pass), 2), "% on average (SD = ", round(sd(temp$pct_pass), 2), "%).\n", sep="")
message(merged_message)
cat(file="log.txt", merged_message, append=T, sep="\n")




#make length distribution plots
if(F){

dir.create(paste(folder, "_stats/length distribution", sep="/"))

message("Generating length distribution plots. If this takes to long you can turn this option off with setting \"LDist=F\".")

for (i in 1:length(new_names)){

pdfname <- sub("/_data/", "/_stats/length distribution/", new_names[i])
pdfname <- sub(".fast.", ".pdf", pdfname)

message(paste("Plotting ", sub(".*distribution/(.*)_PE_.*pdf","\\1", pdfname), sep=""))
Length_distribution(new_names[i], pdfname)
}
message(" ")
}# Ldist end


message(" ")
message("Module completed!")

cat(file="log.txt", paste(Sys.time(), "*** Module completed!", "", sep="\n"), append=T, sep="\n")


}

