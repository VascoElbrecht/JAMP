# U_cluster_otus v0.1

U_cluster_otus() <- function(files="latest", RC=c(T,F), fastq=T, copy_unchanged=T){

Core(module="U_cluster_otus()")
cat(file="../log.txt", c("Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files=="latest"){
source("robots.txt")
files <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T)
}

# CREATE CLUSTERING FOLDERS IN DATA





# new file names
new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub("_PE.", "_PE_RC.", new_names)


cmd <- paste("-fastx_revcomp \"", files[RC], if(fastq){"\" -fastqout \""} else {" -fastaout \""}, new_names[RC], "\" -label_suffix _RC",  sep="")

temp <- paste("RevComp is generated for the following ", length(files[RC]), " files:", sep="")
cat(file="../log.txt", temp , append=T, sep="\n")
message(temp)


temp <- new_names[RC]
for (i in 1:length(cmd)){
system2("usearch", cmd[i], stdout=T, stderr=T)
meep <- sub(".*_data/(.*)", "\\1", temp[i])
cat(file="../log.txt", meep, append=T, sep="\n")
message(meep)
}




message(" ")
message(" Module completed!")

cat(file="../log.txt", paste(Sys.time(), "\n", "Module completed!", "", sep="\n"), append=T, sep="\n")

setwd("../")
}

