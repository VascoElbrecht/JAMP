# U_subset v0.1

U_subset <- function(files="latest", sample_size=100000, fastq_out=F, ranseed=NA){

folder <- Core(module="U_subset")
cat(file="log.txt", c("\n","Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files[1]=="latest"){
source(paste(folder, "/robots.txt", sep=""))
files <- list.files(paste(last_data, "/_data", sep=""), full.names=T)
}

temp <- paste("Starting to subsample to " , sample_size, " sequences) in ", length(files), " samples in ", if(fastq_out){"fastq"}else{"fasta"}, " format.\n\nFiles processed:", sep="")
message(temp)
#message(" ")
cat(file="log.txt", temp, append=T, sep="\n")


# write full number!
temp <- options()$scipen
options("scipen" = 100)
sample_size <- as.character(sample_size)
options("scipen" = temp)



# new file names
new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- paste(folder, "/", new_names, sep="")


if(fastq_out){
new_names <- sub(".fastq", "_sub.fastq", new_names) # keep fastq
} else {
new_names <- sub(".fast.", "_sub.fasta", new_names) # keep fasta
}
new_names <- sub("_sub", paste("_N", sample_size, sep=""), new_names)

dir.create(paste(folder, "/_stats/subset_stats", sep=""))
log_names <- sub("_data", "_stats/subset_stats", new_names)
log_names <- sub("_sub.fast.", "_sub.txt", log_names)


# cmd max EE
cmd <- paste("-fastx_subsample \"", files, "\"", if(fastq_out){" -fastqout "} else {" -fastaout "}, "\"", new_names, "\" -sample_size ", sample_size, " -sizein -sizeout", if(!is.na(ranseed)){paste(" ranseed ", ranseed)}, sep="")


tab_exp <- NULL
for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
cat(paste("usearch", cmd[i], "\n"), file=log_names[i], sep="\n", append=F)
cat(A, file=log_names[i], sep="\n", append=T)

message(sub("_data/","", new_names[i]))
cat(file="log.txt", sub("_data/","", new_names[i]), append=T, sep="\n")

}
cat(file="log.txt", "\n", append=T, sep="\n")

# make some plots?

message(" ")
message("Module completed!")

cat(file="log.txt", paste(Sys.time(), "*** Module completed!", "", sep="\n"), append=T, sep="\n")


}

