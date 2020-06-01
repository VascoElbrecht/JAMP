# Count_sequences v0.1

Count_sequences <- function(files=NA, fastq=T, count_size=F){

# add autodetect of fastq or not!

empty <- file.info(files)$size==0

if(!count_size){
if(fastq){ #count fastq
cmd <- paste("", files, " | wc -l", sep="")

A <- rep(0, length(files))
for (i in which(!empty)){
A[i] <- system2("cat", cmd[i], stdout=T)
}

abundance <- as.numeric(sub(" ", "", A))/4

}

if(!fastq){ # count fasta
cmd <- paste(" -c \"^>\" \"", files, "\"", sep="")
#cmd <- paste(" -o \"^>\" \"", files, "\" | wc -l", sep="") # old

A <- rep(0, length(files))
for (i in which(!empty)){
A[i] <- system2("grep", cmd[i], stdout=T)
}

abundance <- as.numeric(sub(" ", "", A))

}
}


if(count_size){ # count fasta - size annotations

A <- rep(0, length(files))
for (i in which(!empty)){
fasta_sequ <- readLines(files[i])
fasta_sequ <- fasta_sequ[grep(">", fasta_sequ)]

if(grepl(";$", fasta_sequ[1])){
A[i] <- sum(as.numeric(sub(".*;size=(.*);", "\\1", fasta_sequ)))
} else {
A[i] <- sum(as.numeric(sub(".*;size=(.*)", "\\1", fasta_sequ)))
}

}

abundance <- A
}
return(abundance)

}
# 200601
# Usearch
# fasta_sequ <- "NAME;size=327;"
#
# Vsearch
# fasta_sequ <- "NAME;size=327"