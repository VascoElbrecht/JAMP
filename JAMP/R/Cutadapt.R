# Cutadapt v0.1

Cutadapt <- function(files="latest", forward=NA, reverse=NA, fastq=T){

FW_only <- F
if (is.na(reverse)){FW_only <- T}


Core(module="Cutadapt")
cat(file="../log.txt", c("Module Version: v0.1", "\n"), append=T, sep="\n")

# cutadapt version
temp <- paste("Using Version: ", "Cutadapt v", system2("cutadapt", "--v", stdout=T, stderr=T), sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")
message(" ")


if (files[1]=="latest"){
source("robots.txt")
files <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T)
}

# new file names
new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub(".fast", "_cut.fast", new_names)

if(!fastq){ #rename to fasta if fastq=F
new_names <- sub(".fastq", ".fasta", new_names)
}

# get primer sequences / names
primers <- read.csv(paste(system.file(package="JAMP"), "/primers.csv", sep=""), stringsAsFactors=F)


# match primer number to number of files!

if(length(files)>1 & length(forward)==1){
forward <- rep(forward, length(files))
}
if(length(files)>1 & length(reverse)==1){
reverse <- rep(reverse, length(files))
}



# replace primer name with sequence
if(length(forward)==1){
fw <- primers$Primer_Sequence[primers$Primer_Name==forward]
if(length(fw)==0){fw <- forward}
}
if(length(reverse)==1){
rw <- primers$Primer_Sequence[primers$Primer_Name==reverse]
if(length(fw)==0){rw <- reverse}
}

# match multiple primer names
if(length(forward)>1){
fw <- primers$Primer_Sequence[match(forward, primers$Primer_Name)]
fw[is.na(fw)] <- forward[is.na(fw)]
}
if(length(reverse)>1){
rw <- primers$Primer_Sequence[match(reverse, primers$Primer_Name)]
rw[is.na(rw)] <- reverse[is.na(rw)]
}

#build rev comp of rw, using "seqinr"
for (i in which(!is.na(rw))){
rw[i] <- paste(rev(comp(unlist(strsplit(rw[i], "")), forceToLower=F, ambiguous=T)), collapse="")
}


# add: write down used primers in log!
if (FW_only){
cmd1 <- paste("-g ^", fw, " -o \"", new_names, "\" \"", files, "\"", " -f ", if(fastq){"fastq"}else{"fasta"}, " --discard-untrimmed", sep="") # forward adapter

temp <- paste("Starting to remove adapters (primers) on forward direction only ", length(cmd1), " files:", sep="")
cat(file="../log.txt", temp , append=T, sep="\n")
message(temp)
}else{
cmd1 <- paste("-g ^", fw, " -o _data/temp.txt \"", files, "\"", " -f ", if(fastq){"fastq"}else{"fasta"}, " --discard-untrimmed", sep="") # forward adapter
cmd2 <- paste("-a ", rw, "$ -o \"", new_names, "\" _data/temp.txt -f ", if(fastq){"fastq"}else{"fasta"}, " --discard-untrimmed", sep="") #rverse adapter

temp <- paste("Starting to remove adapters (primers) on both ends in ", length(cmd1), " files:", sep="")
cat(file="../log.txt", temp , append=T, sep="\n")
message(temp)
}


dir.create("_stats/_cutadapt_logs")
log_names <- sub("_data", "_stats/_cutadapt_logs", new_names)
log_names <- sub(".fast[aq]", ".txt", log_names)

exp <- NULL
temp <- new_names
for (i in 1:length(cmd1)){
A <- system2("cutadapt", cmd1[i], stdout=T, stderr=T)
cat(file=log_names[i], A, append=T, sep="\n")
if(!FW_only){
A <- system2("cutadapt", cmd2[i], stdout=T, stderr=T)
cat(file=log_names[i], A, append=T, sep="\n")
}# reporting

stats <- readLines(log_names[i])

reads_in <- stats[grep("Total reads processed:", stats)[1]]
reads_in <- sub(".* processed: +", "", reads_in)
reads_in <- as.numeric(gsub(",", "", reads_in))

reads_out <- stats[grep("Reads written \\(passing filters\\):", stats)[if(FW_only){1}else{2}]]
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
write.csv(exp, "_stats/cut_pass.csv")

if(!FW_only){file.remove("_data/temp.txt")} # remove temporary file
message(" ")
message(" Module completed!")

cat(file="../log.txt", paste(Sys.time(), "\n", "Module completed!", "", sep="\n"), append=T, sep="\n")

setwd("../")
}

