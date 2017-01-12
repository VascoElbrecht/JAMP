# Minmax v0.1

Minmax <- function(files="latest", min=NA, max=NA, plusminus=c(NA, 10), fastq=T){

Core(module="Minmax")
cat(file="../log.txt", c("Module Version: v0.1", "\n"), append=T, sep="\n")

# cutadapt version
temp <- paste("Util Version: ", "Cutadapt v", system2("cutadapt", "--v", stdout=T, stderr=T), sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")
message(" ")


if (files[1]=="latest"){
source("robots.txt")
files <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T)
}

# new file names
new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub(".fast", "_minmax.fast", new_names)

#log names
dir.create("_stats/_cutadapt_logs")
log_names <- sub("_data", "_stats/_cutadapt_logs", new_names)
log_names <- sub(".fast[aq]", ".txt", log_names)

# NEED TO REDO THIS, for now only plusminus works, no individual scores, no only onesided
# calculate min max
if(!is.na(plusminus[1])){
min <- plusminus[1] - 10
max <- plusminus[1] + 10
}


# make cmd
cmd <- paste("\"", files, "\" -o \"", new_names, "\" -f ", if(fastq){"fastq"}else{"fasta"}, if(!is.na(min)){" -m "}, if(!is.na(min)){min}, if(!is.na(min)){" -M "}, if(!is.na(min)){max}, sep="")

# report
temp <- paste("Starting to discard reads that don't fit the target length in ", length(cmd), " files:", sep="")
cat(file="../log.txt", temp , append=T, sep="\n")
message(temp)


exp <- NULL
temp <- new_names
for (i in 1:length(cmd)){
A <- system2("cutadapt", cmd[i], stdout=T, stderr=T)
cat(file=log_names[i], A, append=T, sep="\n")

# reporting
stats <- readLines(log_names[i])

reads_in <- stats[grep("Total reads processed:", stats)[1]]
reads_in <- sub(".* processed: +", "", reads_in)
reads_in <- as.numeric(gsub(",", "", reads_in))

reads_out <- stats[grep("Reads written \\(passing filters\\):", stats)[1]]
reads_out <- sub(".* filters.: +", "", reads_out)
reads_out <- sub(" .*", "", reads_out)
reads_out <- as.numeric(gsub(",", "", reads_out))

keep <- round(reads_out/reads_in*100, digits=2)
exp <- rbind(exp, c(sub(".*_data/(.*)", "\\1", temp[i]), reads_out, keep))

meep <- paste(sub(".*_data/(.*)_PE.*", "\\1", temp[i]), " - ", keep, "% reads passed", sep="")
cat(file="../log.txt", meep, append=T, sep="\n")
message(meep)
}

exp <- data.frame(exp)
names(exp) <- c("Sample", "Abundance", "pct_pass")
write.csv(exp, "_stats/minmax_pass.csv")

message(" ")
message(" Module completed!")

cat(file="../log.txt", paste(Sys.time(), "\n", "Module completed!", "", sep="\n"), append=T, sep="\n")

setwd("../")
}

