# U_truncate v0.1

U_truncate <- function(files="latest", left=0, right=0, trunclen=NA, fastq=T, rename=T, keepshort=F, exe="usearch", delete_data=T){

folder <- Core(module="U_truncate", delete_data=delete_data)
cat(file="log.txt", c("\n","Version v0.2", "\n"), append=T, sep="\n")
message(" ")

files_to_delete <- NULL

if (files[1]=="latest"){
source(paste(folder, "/robots.txt", sep=""))
files <- list.files(paste(last_data, "/_data", sep=""), full.names=T)
}

temp <- paste("Starting to truncate reads of ", length(files), " samples.", sep="")
message(temp)
message(" ")
cat(file="log.txt", temp, append=T, sep="\n")



# new file names

new_names <- sub(".*(_data/.*)", "\\1", files)
if(rename){ # rename files to indicate trimming
new_names <- sub(".fast", "_trunc.fast", new_names)
}
new_names <- paste(folder, "/", new_names, sep="")


# allow for vector truncation (in case of r1 or r2)
if(sub(".*/", "", exe)=="usearch"){
cmd <- paste("-fastx_truncate \"", files,"\"", " -stripleft ", left, " -stripright ", right, if(!is.na(trunclen)){paste(" -trunclen ", trunclen, sep="")}, if(keepshort){paste(" -padlen ", trunclen, sep="")}, if(fastq){" -fastqout "} else {" -fastaout "}, "\"", new_names, "\"", sep="")
} else { # vsearch version
cmd <- paste("-fastx_filter \"", files,"\"", " -fastq_stripleft ", left, " -fastq_stripright ", right, if(!is.na(trunclen)){paste(" -fastq_trunclen ", trunclen, sep="")}, if(keepshort){paste(" -padlen ", trunclen, sep="")}, if(fastq){" -fastqout "} else {" -fastaout "}, "\"", new_names, "\"", sep="")
}


files_to_delete <- c(files_to_delete, new_names)

tab_exp <- NULL
for (i in 1:length(cmd)){
system2(exe, cmd[i], stdout=T, stderr=T)

new_count <- Count_sequences(new_names[i], fastq= fastq)
old_count <- Count_sequences(files[i], fastq= fastq)
passed <- round(new_count/old_count*100, digits=2)

A <- system2(exe, paste("-fastx_info \"", new_names[i], "\" -secs 5", sep=""), stdout=T, stderr=T)
medianL <- as.numeric(sub(".*median (.*), hi.*", "\\1", A[grep("Lengths min ", A)]))

cat(file=paste(folder, "/_stats/log_length.txt", sep=""), new_names[i], "\n", A,"\n\n", append=T, sep="\n")


tab_exp <- rbind(tab_exp, c(sub(".*_data/(.*)", "\\1", new_names[i]), new_count, passed, medianL))

meep <- paste(sub(".*_data/(.*)", "\\1", new_names[i]), ": ", passed, "% passed - medianL: ", medianL, sep="")
message(meep)
cat(file="log.txt", meep, append=T, sep="\n")
}

cat(file="log.txt", "\n", append=T, sep="\n")

tab_exp <- data.frame(tab_exp)
names(tab_exp) <- c("Sample", "Abundance", "pct_pass", "medianL")
write.csv(tab_exp, paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_truncate_pass.csv", sep=""))


message(" ")
message("Module completed!")

cat(file=paste(folder, "/robots.txt", sep=""), "\n# DELETE_START", files_to_delete, "# DELETE_END", append=T, sep="\n")

cat(file="log.txt", paste(Sys.time(), "*** Module completed!", "", sep="\n"), append=T, sep="\n")
}

