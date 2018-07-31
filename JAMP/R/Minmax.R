# Minmax v0.1

Minmax <- function(files="latest", min=NA, max=NA, plusminus=c(NA, 10), fastq=T, LDist=F){

folder <- Core(module="Minmax")
cat(file="../log.txt", c("Module Version: v0.2", "\n"), append=T, sep="\n")

# cutadapt version
temp <- paste("Version: ", "Cutadapt v", system2("cutadapt", "--v", stdout=T, stderr=T), sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")
message(" ")


if (files[1]=="latest"){
source(paste(folder, "/robots.txt", sep=""))
files <- list.files(paste(last_data, "/_data", sep=""), full.names=T)
}

# new file names
new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub(".fast", "_minmax.fast", new_names)

#log names
dir.create(paste(folder, "/_stats/_cutadapt_logs", sep=""))
log_names <- sub("_data", "_stats/_cutadapt_logs", new_names)
log_names <- sub(".fast[aq]", ".txt", log_names)

# NEED TO REDO THIS, for now only plusminus works, no individual scores, no only onesided
# calculate min max
if(!is.na(plusminus[1])){
min <- plusminus[1] - plusminus[2]
max <- plusminus[1] + plusminus[2]
}


# make cmd
cmd <- paste("\"", files, "\" -o \"", folder, "/", new_names, "\" -f ", if(fastq){"fastq"}else{"fasta"}, if(!is.na(min)){" -m "}, if(!is.na(min)){min}, if(!is.na(min)){" -M "}, if(!is.na(min)){max}, sep="")

# report
temp <- paste("Starting to discard reads that don't fit the target length in ", length(cmd), " files:", sep="")
cat(file="log.txt", temp , append=T, sep="\n")
message(temp)


exp <- NULL
temp <- new_names
for (i in 1:length(cmd)){
A <- system2("cutadapt", cmd[i], stdout=T, stderr=T)
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
write.csv(exp, paste(folder, "/_stats/minmax_pass.csv", sep=""))

# make plots


# make some plots
temp <- read.csv(paste(folder, "/_stats/minmax_pass.csv", sep=""), stringsAsFactors=F)

Sequences_lost(temp$Sequ_count_in, temp$Sequ_count_out, sub("_PE.*", "", temp$Sample), out=paste(folder, "/_stats/Sequences_discarded.pdf", sep=""),  main=paste(folder, ": Proportion of reads with ", min, "-", max, "bp length", sep=""))
Sequences_lost(temp$Sequ_count_in, temp$Sequ_count_out, sub("_PE.*", "", temp$Sample), rel=T, out=paste(folder, "/_stats/Sequences_discarded_rel.pdf", sep=""), main=paste(folder, ": Proportion of reads with ", min, "-", max, "bp length", sep=""))

merged_message <- paste("\nProportion of sequences of expected length: ", round(mean(temp$pct_pass), 2), "% on average (SD = ", round(sd(temp$pct_pass), 2), "%,  Median = ", median(temp$pct_pass), "%).\n", sep="")
message(merged_message)
cat(file="log.txt", merged_message, append=T, sep="\n")


#make length distribution plots
if(LDist){

dir.create(paste(folder, "_stats/length distribution", sep="/"))

message("Generating length distribution plots. If this takes to long you can turn this option off with setting \"LDist=T\".")

for (i in 1:length(new_names)){

pdfname <- sub("_data/", "_stats/length distribution/", new_names[i])
pdfname <- sub(".fast.", ".pdf", pdfname)

message(paste("Plotting ", sub(".*distribution/(.*)_PE_.*pdf","\\1", pdfname), sep=""))
Length_distribution(paste(folder, "/", new_names[i], sep=""), paste(folder, "/", pdfname, sep=""), fastq=fastq)
}
message(" ")
}# Ldist end

message(" ")
message(" Module completed!")

cat(file="log.txt", paste(Sys.time(), "\n", "*** Module completed!", "", sep="\n"), append=T, sep="\n")


}

