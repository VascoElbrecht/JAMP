# Cutadapt v0.1

Cutadapt <- function(files="latest", forward=NA, reverse=NA, fastq=T){

Core(module="Cutadapt")
cat(file="../log.txt", c("Module Version: v0.1", "\n"), append=T, sep="\n")

# cutadapt version
temp <- paste("Util Version: ", "Cutadapt v", system2("cutadapt", "--v", stdout=T, stderr=T), sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")
message(" ")


if (files=="latest"){
source("robots.txt")
files <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T)
}

# new file names
new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub(".fast", "_cut.fast", new_names)

# get primer sequences / names
primers <- read.csv(paste(system.file(package="JAMP"), "/primers.csv", sep=""), stringsAsFactors=F)


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
}
if(length(reverse)>1){
rw <- primers$Primer_Sequence[match(reverse, primers$Primer_Name)]
}

#build rev comp of rw, using "seqinr"
for (i in 1:length(rw)){
rw[i] <- paste(rev(comp(unlist(strsplit(rw[i], "")), forceToLower=F, ambiguous=T)), collapse="")
}


# add: write down used primers in log!

cmd1 <- paste("-g ^", fw, " -o _data/temp.txt \"", files, "\"", " -f ", if(fastq){"fastq"}else{"fasta"}, " --discard-untrimmed", sep="") # forward adapter
cmd2 <- paste("-a ", rw, "$ -o ", new_names, " _data/temp.txt -f ", if(fastq){"fastq"}else{"fasta"}, " --discard-untrimmed", sep="") #rverse adapter

temp <- paste("Starting to remove adapters (primers) on both ends in ", length(cmd1), " files:", sep="")
cat(file="../log.txt", temp , append=T, sep="\n")
message(temp)


temp <- new_names
for (i in 1:length(cmd1)){
system2("cutadapt", cmd1[i], stdout=T, stderr=T)
system2("cutadapt", cmd2[i], stdout=T, stderr=T)

# reporting
meep <- sub(".*_data/(.*)", "\\1", temp[i])
cat(file="../log.txt", meep, append=T, sep="\n")
message(meep)
}

# delete temp primer file!!!



message(" ")
message(" Module completed!")

cat(file="../log.txt", paste(Sys.time(), "\n", "Module completed!", "", sep="\n"), append=T, sep="\n")

setwd("../")
}

