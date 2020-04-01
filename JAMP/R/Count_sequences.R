# Count_sequences v0.1

Count_sequences <- function(files=NA, fastq=T, count_size=F){

# add autodetect of fastq or not!



if(!count_size){
if(fastq){ #count fastq
cmd <- paste("", files, " | wc -l", sep="")

A <- NULL
for (i in 1:length(files)){
A[i] <- system2("cat", cmd[i], stdout=T)
}

abundance <- as.numeric(sub(" ", "", A))/4

}

if(!fastq){ # count fasta
cmd <- paste(" -c \"^>\" \"", files, "\"", sep="")
#cmd <- paste(" -o \"^>\" \"", files, "\" | wc -l", sep="") # old

A <- NULL
for (i in 1:length(files)){
A[i] <- system2("grep", cmd[i], stdout=T)
}

abundance <- as.numeric(sub(" ", "", A))

}
}


if(count_size){ # count fasta - size annotations

A <- NULL
for (i in 1:length(files)){
fasta_sequ <- readLines(files[i])
fasta_sequ <- fasta_sequ[grep(">", fasta_sequ)]
A[i] <- sum(as.numeric(sub(".*;size=(.*);?", "\\1", fasta_sequ)))
}

abundance <- as.numeric(sub(" ", "", A))

}

return(abundance)
}
