# Cutadapt v0.1

Cutadapt <- function(files="latest", forward=NA, reverse=NA, bothsides=F, anchoring=T, LDist=F, cores=1, delete_data=T, ambsequ=F, error=0.1, exe="cutadapt", exe2="vsearch", fastq=T){

# do checks
FW_only <- F
if (is.na(reverse[1])){FW_only <- T}



folder <- Core(module="Cutadapt", delete_data=delete_data)
cat(file="log.txt", c("Module Version: v0.2", "\n"), append=T, sep="\n")

files_to_delete <- NULL

# cutadapt version
temp <- paste("Using Version: ", "Cutadapt v", system2(exe, "--v", stdout=T, stderr=T), sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")
message(" ")


if (files[1]=="latest"){
source(paste(folder, "/robots.txt", sep=""))
files <- list.files(paste(last_data, "/_data", sep=""), full.names=T)
}

# new file names
new_names <- sub(".*/(.*)", "\\1", files)
new_names <- sub(".fast", "_cut.fast", new_names)
new_names <- paste(folder, "/_data/", new_names, sep="")

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
if(length(rw)==0){rw <- reverse}
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

#replace I (inosine) with N
fw <- gsub("I", "N", fw)
rw <- gsub("I", "N", rw)


#build rev comp of rw, using "seqinr"
for (i in which(!is.na(rw))){
rw[i] <- paste(rev(comp(unlist(strsplit(rw[i], "")), forceToLower=F, ambiguous=T)), collapse="")
}


# add: write down used primers in log!
if (FW_only){
cmd1 <- paste("-g ", if(anchoring){"^"}, fw, " -o \"", new_names, "\" \"", files, "\"", " --discard-untrimmed --cores=", cores, if(ambsequ){" --match-read-wildcards"}, " --error-rate ", error, sep="") # forward adapter
files_to_delete <- c(files_to_delete, new_names)

temp <- paste("Starting to remove adapters (primers) on forward direction only ", length(cmd1), " files:", sep="")
cat(file="log.txt", temp , append=T, sep="\n")
message(temp)
}else{
cmd1 <- paste("-g ", if(anchoring){"^"}, fw, " -o ", folder, "/_data/temp.txt \"", files, "\"", " --discard-untrimmed --cores=", cores, if(ambsequ){" --match-read-wildcards"}, " --error-rate ", error, sep="") # forward adapter
cmd2 <- paste("-a ", rw, if(anchoring){"$"}," -o \"", new_names, "\" ", folder, "/_data/temp.txt ", " --discard-untrimmed --cores=", cores, if(ambsequ){" --match-read-wildcards"}, " --error-rate ", error, sep="") #rverse adapter
files_to_delete <- c(files_to_delete, new_names)

temp <- paste("Starting to remove adapters (primers) on both ends in ", length(cmd1), " files:", sep="")
cat(file="log.txt", temp , append=T, sep="\n")
message(temp)
}

if(bothsides){
temp <- paste("Notice: Option \"bothsides\" enabled, will also trim adapters on reverse complements and orientate all sequences in the same forward direction.", sep="")
cat(file="log.txt", temp , append=T, sep="\n")
message(temp)
}


dir.create(paste(folder, "/_stats/_cutadapt_logs", sep=""))
log_names <- sub("_data", "_stats/_cutadapt_logs", new_names)
log_names <- sub(".fast[aq]", ".txt", log_names)

exp <- NULL
temp <- new_names

for (i in 1:length(cmd1)){
A <- system2(exe, cmd1[i], stdout=T, stderr=T)
cat(file=log_names[i], A, append=T, sep="\n")
if(!FW_only){
if(!bothsides){ # not both sided
A <- system2(exe, cmd2[i], stdout=T, stderr=T)
cat(file=log_names[i], A, append=T, sep="\n")
} else {
# trimm FW
#A <- system2("cutadapt", cmd1[i], stdout=T, stderr=T)
#cat(file=log_names[i], A, append=T, sep="\n")

# trimm RW
rev_primer <- paste("-a ", rw[i], if(anchoring){"$"}," -o ", folder, "/_data/temp_A.txt ", folder, "/_data/temp.txt ", " --discard-untrimmed --cores=", cores, sep="", if(ambsequ){" --match-read-wildcards"}, " --error-rate ", error) #rverse adapter

A <- system2(exe, rev_primer, stdout=T, stderr=T)
cat(file=log_names[i], A, append=T, sep="\n")

# build revcomp

revcom_cmd <- paste("-fastx_revcomp ", files[i],  if(fastq){" -fastqout "} else {" -fastaout "}, folder, "/_data/temp2_RC.txt  -label_suffix _RC",  sep="")
A <- system2(exe2, revcom_cmd, stdout=T, stderr=T)

cat(file=log_names[i], paste(exe2, revcom_cmd, sep=""), append=T, sep="\n")
cat(file=log_names[i], A, append=T, sep="\n")

# cut primers again! simmillar commands

fw_primer_cmd <- paste("-g ", if(anchoring){"^"}, fw[i], " -o ", folder, "/_data/temp.txt ", folder, "/_data/temp2_RC.txt", " --discard-untrimmed --cores=", cores, if(ambsequ){" --match-read-wildcards"}, sep="", " --error-rate ", error) # forward adapter

A <- system2(exe, fw_primer_cmd, stdout=T, stderr=T)
cat(file=log_names[i], A, append=T, sep="\n")


rev_primer <- paste("-a ", rw[i], if(anchoring){"$"}," -o ", folder, "/_data/temp_B.txt ", folder, "/_data/temp.txt --discard-untrimmed --cores=", cores, if(ambsequ){" --match-read-wildcards"}, sep="", " --error-rate ", error) #rverse adapter

A <- system2(exe, rev_primer, stdout=T, stderr=T)
cat(file=log_names[i], A, append=T, sep="\n")

# combine files!
combine_cmd <- paste("", folder, "/_data/temp_A.txt ", folder, "/_data/temp_B.txt > ", new_names[i], sep="")

system2("cat", combine_cmd)
cat(file=log_names[i], paste("\ncat ", combine_cmd, sep=""), append=T, sep="\n")

# remove temp files
file.remove(paste(folder, "/_data/temp_A.txt", sep=""))
file.remove(paste(folder, "/_data/temp_B.txt", sep=""))
file.remove(paste(folder, "/_data/temp2_RC.txt", sep=""))
}
}
# reporting

stats <- readLines(log_names[i])
reads_in <- stats[grep("Total reads processed:", stats)[1]]
reads_in <- sub(".* processed: +", "", reads_in)
reads_in <- as.numeric(gsub(",", "", reads_in))

if(!bothsides){

reads_out <- stats[grep("Reads written \\(passing filters\\):", stats)[if(FW_only){1}else{2}]]
reads_out <- sub(".* filters.: +", "", reads_out)
reads_out <- sub(" .*", "", reads_out)
reads_out <- as.numeric(gsub(",", "", reads_out))

keep <- round(reads_out/reads_in*100, digits=2)
exp <- rbind(exp, c(sub(".*_data/(.*)", "\\1", temp[i]), reads_in, reads_out, keep))

meep <- paste(sub(".*_data/(.*)_PE.*", "\\1", temp[i]), " - ", keep, "% reads passed", sep="")
cat(file="log.txt", meep, append=T, sep="\n")
message(meep)
} else { # processing for both sided

reads_outA <- stats[grep("Reads written \\(passing filters\\):", stats)[2]]
reads_outA <- sub(".* filters.: +", "", reads_outA)
reads_outA <- sub(" .*", "", reads_outA)
reads_outA <- as.numeric(gsub(",", "", reads_outA))

reads_outB <- stats[grep("Reads written \\(passing filters\\):", stats)[4]]
reads_outB <- sub(".* filters.: +", "", reads_outB)
reads_outB <- sub(" .*", "", reads_outB)
reads_outB <- as.numeric(gsub(",", "", reads_outB))

reads_out <- reads_outA+ reads_outB

keep <- round(reads_out/reads_in*100, digits=2)
exp <- rbind(exp, c(sub(".*/(.*)", "\\1", temp[i]), reads_in, reads_out, keep))

meep <- paste(sub(".*_data/(.*)_PE.*", "\\1", temp[i]), " - ", keep, "% reads passed (", round(reads_outA/reads_in*100, digits=2), "% forward + ", round(reads_outB/reads_in*100, digits=2), "% reverse orientation)", sep="")
cat(file="log.txt", meep, append=T, sep="\n")
message(meep)

}
}

exp <- data.frame(exp)
names(exp) <- c("Sample", "Sequ_count_in", "Sequ_count_out", "pct_pass")
write.csv(exp, paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_cut_pass.csv", sep=""))

if(!FW_only){file.remove(paste(folder, "/_data/temp.txt", sep=""))} # remove temporary file




# make some plots
temp <- read.csv(paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_cut_pass.csv", sep=""), stringsAsFactors=F)

Sequences_lost(temp$Sequ_count_in, temp$Sequ_count_out, sub("_PE.*", "", temp$Sample), out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_Primers_trimmed.pdf", sep=""), main=paste(folder, ": Proportion of reads with no primer detected", sep=""))
Sequences_lost(temp$Sequ_count_in, temp$Sequ_count_out, sub("_PE.*", "", temp$Sample), rel=T, out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_Primers_trimmed_rel.pdf", sep=""),  main=paste(folder, ": Proportion of reads with no primer detected", sep=""))

merged_message <- paste("\nPrimers were successfully removed from ", round(mean(temp$pct_pass), 2), "% on average (SD = ", round(sd(temp$pct_pass), 2), "%).\n", sep="")
message(merged_message)
cat(file="log.txt", merged_message, append=T, sep="\n")




#make length distribution plots
if(LDist){

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
message(" Module completed!")

cat(file=paste(folder, "/robots.txt", sep=""), "\n# DELETE_START", files_to_delete, "# DELETE_END", append=T, sep="\n")

cat(file="log.txt", paste(Sys.time(), "\n", "*** Module completed!\n\n", sep="\n"), append=T, sep="\n")


}

