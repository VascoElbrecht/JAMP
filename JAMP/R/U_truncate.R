# U_truncate v0.1

U_truncate <- function(files="latest", left=0, right=0, fastq=T){

Core(module="U_truncate")
cat(file="../log.txt", c("\n","Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files=="latest"){
source("robots.txt")
files <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T)
}

temp <- paste("Starting to truncate reads of ", length(files), " samples.", sep="")
message(temp)
message(" ")
cat(file="../log.txt", temp, append=T, sep="\n")



# new file names

new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub(".fast", "_trunc.fast", new_names)


cmd <- paste("-fastx_truncate \"", files,"\"", " -stripleft ", left, " -stripright ", right, if(fastq){" -fastqout "} else {" -fastaout "}, "\"", new_names, "\"", sep="")


tab_exp <- NULL
for (i in 1:length(cmd)){
system2("usearch", cmd[i], stdout=T, stderr=T)

new_count <- Count_sequences(new_names[i])
old_count <- Count_sequences(files[i])
passed <- round(new_count/old_count*100, digits=2)

A <- system2("usearch", paste("-fastx_info \"", new_names[i], "\" -secs 5", sep=""), stdout=T, stderr=T)
medianL <- as.numeric(sub(".*med (.*), hi.*", "\\1", A[grep("Lengths min ", A)]))

cat(file="_stats/log_length.txt", new_names[i], "\n", A,"\n\n", append=T, sep="\n")


tab_exp <- rbind(tab_exp, c(sub("_data/(.*)_PE_.*", "\\1", new_names[i]), new_count, passed, medianL))

meep <- paste(sub("_data/(.*)_PE_.*", "\\1", new_names[i]), ": ", passed, "% passed - medianL: ", medianL, sep="")
message(meep)
cat(file="../log.txt", meep, append=T, sep="\n")
}

cat(file="../log.txt", "\n", append=T, sep="\n")

tab_exp <- data.frame(tab_exp)
names(tab_exp) <- c("Sample", "Abundance", "pct_pass", "medianL")
write.csv(tab_exp, "_stats/truncate_pass.csv")


message(" ")
message("Module completed!")

cat(file="../log.txt", paste(Sys.time(), "Module completed!", "", sep="\n"), append=T, sep="\n")

setwd("../")
}

