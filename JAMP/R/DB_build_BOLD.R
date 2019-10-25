# 191025 convert BOLD ref data to DB

# setwd("~/Desktop/bold")
tsv="bold_data.tsv"

function <- buildBOLDdb(tsv="ref_tsv_fro_bold.txt", fasta="Save_file_as.fasta", csv="BOLD_taxonomy.csv", minlength=500, maxlength=NA){

message("Building BOLD based refference database using: ", tsv)

data <- read.csv(tsv, sep="\t", stringsAsFactors=F)

message("Sequences detected: ", nrow(data), "\nRemoving gaps from sequences!")
data$nucleotides <- gsub("-", "", data$nucleotides)

message("Removing terminal Ns fromm sequences!")
data$nucleotides <- sub("N*(.*[ACGTUWSMKRYBDHV])N*", "\\1", data$nucleotides)

if(is.na(minlength)&is.na(maxlength)){
message("Not applying any length filtering, minlength and maxlength are both set to NA.")
} else {sequlength <- nchar(data$nucleotides)}




write.csv(data, "test.csv")



#table(A-B)
#table(A)
#length(A)
#table(B)
#head(data[A==0,])

# keep only sequences of 658 bp


data2 <- data[B>=658-10&B<=658+10,]

table(nchar(data2$nucleotides))

nrow(data)
nrow(data2)



length(unique(data$bin_uri))
length(unique(data2$bin_uri))


#remaining bins after filtering: 27970 of 32921

head(data2)

# write ref DB

cat(paste(paste(">", data2$processid, "__", data2$bin_uri, "\n", sep=""), data2$nucleotides, "\n", sep=""), file="Ontraio_edited_ref_v2.fasta", sep="")

nrow(data2)


}












