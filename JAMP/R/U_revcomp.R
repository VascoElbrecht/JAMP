# U_merge_PE v0.1
# maybe add option to split and merge large files automatically?

U_revcomp <- function(files="latest", fastq=T){

Core(module="U_revcomp")
cat(file="../log.txt", c("Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files=="latest"){
source("robots.txt")
files <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T)
}

# copy over files which are not rev comp


# new file names
new_names <- sub(".*(_data/.*)", "\\1", file1)

cmd <- paste("",  sep="")


"-fastx_revcomp ", temp.txt, " -fastqout ", 


tab_exp <- NULL
for (i in 1:length(cmd)){
system2("usearch", cmd[i], stdout=F, stderr=F)
temp <- readLines(log_names[i])
cat(file="../log.txt", meep, append=T, sep="\n")
}

cat(file="../log.txt", "\n", append=T, sep="\n")



write.csv(tab_exp, "_stats/sequ_length_abund.csv")

# make some plots?

message(" ")
message(" Done with PE merging")

cat(file="../log.txt", paste(Sys.time(), "Done with PE merging", "", sep="\n"), append=T, sep="\n")

setwd("../")
}

