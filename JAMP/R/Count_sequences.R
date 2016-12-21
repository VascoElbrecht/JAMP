# Count_sequences v0.1

Count_sequences <- function(files=NA, fastq=T){


if(fastq){ #count fastq
cmd <- paste("", files, " | wc -l", sep="")

A <- NULL
for (i in 1:length(files)){
A[i] <- system2("cat", cmd[i], stdout=T)
}

abundance <- as.numeric(sub(" ", "", A))/4

} else { # count fasta
cmd <- paste(" -o \">\" ", files, " | wc -l", sep="")

A <- NULL
for (i in 1:length(files)){
A[i] <- system2("fgrep", cmd[i], stdout=T)
}

abundance <- as.numeric(sub(" ", "", A))


}

return(abundance)
}
