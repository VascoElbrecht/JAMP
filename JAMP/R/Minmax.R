# Minmax v0.2

Minmax <- function(files="latest", min=NA, max=NA, plusminus=c(NA, 10), fastq=T, LDist=F, delete_data=T, exe="cutadapt"){

folder <- Core(module="Minmax", delete_data=delete_data)
cat(file="../log.txt", c("Module Version: v0.2", "\n"), append=T, sep="\n")

files_to_delete <- NULL

# cutadapt version
temp <- paste("Version: ", "Cutadapt v", system2(exe, "--v", stdout=T, stderr=T), sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")
message(" ")


if (files[1]=="latest"){
source(paste(folder, "/robots.txt", sep=""))
files <- list.files(paste(last_data, "/_data", sep=""), full.names=T)
}

# new file names
new_names <- sub(".*/(.*)", "\\1", files)
if(fastq){new_names <- sub(".fastq$", "_minmax.fastq", new_names)}else{
new_names <- sub(".fasta$", "_minmax.fasta", new_names)}
new_names <- paste("_data/", new_names, sep="")

#log names
dir.create(paste(folder, "/_stats/_cutadapt_logs", sep=""))
log_names <- sub("_data", "_stats/_cutadapt_logs", new_names)
log_names <- sub(".fast[aq]$", ".txt", log_names)

# NEED TO REDO THIS, for now only plusminus works, no individual scores, no only onesided
# calculate min max
if(!is.na(plusminus[1])){
min <- plusminus[1] - plusminus[2]
max <- plusminus[1] + plusminus[2]
}


# make cmd
cmd <- paste("\"", files, "\" -o \"", folder, "/", new_names, "\"", if(!is.na(min[1])){" -m "}, if(!is.na(min[1])){min}, if(!is.na(min[1])){" -M "}, if(!is.na(min[1])){max}, sep="")

files_to_delete <- c(files_to_delete, paste(folder, "/", new_names, sep=""))

# report
temp <- paste("Starting to discard reads that don't fit the target length in ", length(cmd), " files:", sep="")
cat(file="log.txt", temp , append=T, sep="\n")
message(temp)


exp <- NULL
temp <- new_names
for (i in 1:length(cmd)){
A <- system2(exe, cmd[i], stdout=T, stderr=T)
cat(file=paste(folder, "/", log_names[i], sep=""), A, append=T, sep="\n")

# reporting
stats <- readLines(paste(folder, "/", log_names[i], sep=""))

reads_in <- stats[grep("Total reads processed:", stats)[1]]
reads_in <- sub(".* processed: +", "", reads_in)
reads_in <- as.numeric(gsub(",", "", reads_in))

reads_out <- stats[grep("Reads written \\(passing filters\\):", stats)[1]]
reads_out <- sub(".* filters.: +", "", reads_out)
reads_out <- sub(" .*", "", reads_out)
reads_out <- as.numeric(gsub(",", "", reads_out))

keep <- round(reads_out/reads_in*100, digits=2)
exp <- rbind(exp, c(sub(".*_data/(.*)", "\\1", temp[i]), reads_in, reads_out, keep))

meep <- paste(sub(".*_data/(.*)_PE.*", "\\1", temp[i]), " - ", keep, "% reads passed", sep="")
cat(file="log.txt", meep, append=T, sep="\n")
message(meep)
}

exp <- data.frame(exp)
names(exp) <- c("Sample", "Sequ_count_in", "Sequ_count_out", "pct_pass")
write.csv(exp, paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_minmax_pass.csv", sep=""))

# make plots


# make some plots
temp <- read.csv(paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_minmax_pass.csv", sep=""), stringsAsFactors=F)

Sequences_lost(temp$Sequ_count_in, temp$Sequ_count_out, sub("_PE.*", "", temp$Sample), out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_minmax_discarded.pdf", sep=""),  main=paste(folder, ": Proportion of reads with ", min, "-", max, "bp length", sep=""))
Sequences_lost(temp$Sequ_count_in, temp$Sequ_count_out, sub("_PE.*", "", temp$Sample), rel=T, out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_minmax_discarded_rel.pdf", sep=""), main=paste(folder, ": Proportion of reads with ", min, "-", max, "bp length", sep=""))

merged_message <- paste("\nProportion of sequences of expected length: ", round(mean(temp$pct_pass), 2), "% on average (SD = ", round(sd(temp$pct_pass), 2), "%,  Median = ", median(temp$pct_pass), "%).\n", sep="")
message(merged_message)
cat(file="log.txt", merged_message, append=T, sep="\n")


#make length distribution plots
if(LDist){

dir.create(paste(folder, "_stats/length distribution", sep="/"))

message("Generating length distribution plots. If this takes to long you can turn this option off with setting \"LDist=T\".")

for (i in 1:length(new_names)){

pdfname <- sub("_data/", "_stats/length distribution/", new_names[i])
pdfname <- sub(".fast.$", ".pdf", pdfname)

message(paste("Plotting ", sub(".*distribution/(.*)_PE_.*pdf","\\1", pdfname), sep=""))
Length_distribution(paste(folder, "/", new_names[i], sep=""), paste(folder, "/", pdfname, sep=""), fastq=fastq)
}
message(" ")
}# Ldist end

message(" ")
message(" Module completed!")

cat(file=paste(folder, "/robots.txt", sep=""), "\n# DELETE_START", files_to_delete, "# DELETE_END", append=T, sep="\n")

cat(file="log.txt", paste(Sys.time(), "\n", "*** Module completed!", "", sep="\n"), append=T, sep="\n")
}

