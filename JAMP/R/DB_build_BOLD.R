# 191025 convert BOLD ref data to DB


buildBOLDdb <- function(tsv="ref_tsv_fro_bold.txt", savefasta=NA, savetaxonomy=NA, minlength=500, maxlength=NA, blacklist=NA, GenBank=F){

message("Building BOLD based refference database using: ", tsv)
data <- read.csv(tsv, sep="\t", stringsAsFactors=F)

if(!GenBank){
GB <- data$institution_storing=="Mined from GenBank, NCBI"
message("Excluding sequences mined from Genbank: ", sum(GB), "/", nrow(data), " (", round(sum(GB)/nrow(data)*100, 2), "% removed).")
data <- data[!GB,]
}



message("Removing unneeded coulumns from table.")
data <- data[,c(1, 8:24, 27, 72)]


if(!is.na(blacklist)){
rm <- data$processid %in% blacklist
data <- data[!rm,]
message("Removing ", sum(rm), " sequences based on ", length(blacklist), " processids provided in blacklist (", round(100-nrow(data)/length(rm)*100,2), "% sequences removed)")
}


message("Sequences detected: ", nrow(data), "\nRemoving gaps from sequences!")
data$nucleotides <- gsub("-", "", data$nucleotides)

message("Removing terminal Ns fromm sequences!")
data$nucleotides <- sub("N*(.*[ACGTUWSMKRYBDHV])N*", "\\1", data$nucleotides)

sequlength <- nchar(data$nucleotides)
message("Removing ", sum(sequlength==0), "/", nrow(data), " specimens that don't come with a sequence (", round(sum(sequlength==0)/nrow(data)*100, 2), "% removed).")
data <- data[sequlength!=0,]


if(is.na(minlength)&is.na(maxlength)){
message("Not applying any length filtering, minlength and maxlength are both set to NA.")
}


if(!is.na(minlength)){
sequlength <- nchar(data$nucleotides)
message("Applying length filtering: ", sum(sequlength<minlength), "/", nrow(data), " reads are shorter than ", minlength, " bp are beeing removed (", round(sum(sequlength<minlength)/nrow(data)*100, 2), "% of sequences removed).")
data <- data[sequlength>=minlength,]
}

if(!is.na(maxlength)){
sequlength <- nchar(data$nucleotides)
message("Applying length filtering: ", sum(sequlength>maxlength), "/", nrow(data), " reads are shorter than ", minlength, " bp are beeing removed (", round(sum(sequlength>maxlength)/nrow(data)*100, 2), "% of sequences removed).")
data <- data[sequlength<maxlength,]
}


data[data==""] <- NA


# catch taxonomic ID

taxon <- rep(0, nrow(data))

taxon[which(!is.na(data$phylum_taxID))] <- 3
taxon[which(!is.na(data$class_taxID))] <- 5
taxon[which(!is.na(data$order_taxID))] <- 7
taxon[which(!is.na(data$family_taxID))] <- 9
taxon[which(!is.na(data$subfamily_taxID))] <- 11
taxon[which(!is.na(data$genus_taxID))] <- 13
taxon[which(!is.na(data$species_taxID))] <- 15
taxon[which(!is.na(data$subspecies_taxID))] <- 17


message("Preparing taxonomy reference files...")
if(is.na(savefatsa)){
savefatsa <- paste("DB_", sub("(.*)\\..*$", "\\1", tsv), ".fasta", sep="")
}

if(is.na(savetaxonomy)){
savetaxonomy <- paste("DB_", sub("(.*)\\..*$", "\\1", tsv), "_taxonomy.csv", sep="")
}


ID <- NULL
for (i in 1:nrow(data)){
ID[i] <- paste(data$processid[i], "_", data$bin_uri[i], "_taxID=BOLD", data[i, taxon[i]], sep="")
}


meep <- rep(NA, nrow(data))
taxonomyDB <- data.frame("taxID"=meep, "phylum"=meep, "class"=meep, "order"=meep, "family"=meep, "genus"=meep, "species"=meep, "subspecies"=meep, "taxRef"=meep)

for (i in 1:nrow(data)){
taxonomyDB[i,] <- c(paste("BOLD:", data[i, taxon[i]], sep=""), as.vector(unlist(data[i, c(4,6,8,10, 14, 16, 18, 19)])))
}

taxonomyDB <- taxonomyDB[!duplicated(taxonomyDB$taxID),]
message(nrow(taxonomyDB), " unique taxonomic IDs detected (from ", nrow(data), " specimens).")
taxonomyDB$taxRef <- gsub("[()]", "", taxonomyDB$taxRef)

write.csv(taxonomyDB, savetaxonomy, row.names=F)

cat(paste(paste(">", ID, "\n", sep=""), data$nucleotides, "\n", sep=""), file= savefatsa, sep="")


message("Module complete! use ", savefatsa, " as a refference database, in combination with the taxonomy table ", savetaxonomy)
}

