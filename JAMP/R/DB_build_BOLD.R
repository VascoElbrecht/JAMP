# 191025 convert BOLD ref data to DB


buildBOLDdb <- function(tsv="ref_tsv_fro_bold.txt", savefasta=NA, savetaxonomy=NA, minlength=500, maxlength=NA, blacklist=NA, GenBank=F, stripptsv=T){



if(file.size(tsv)/1000000>1000 & stripptsv){
message("tsv is larger tan 1 GB, only the relavant columns will be extracted before importing into R!")
tsv_stripped <- sub("\\.[tc].*$", "_stripped.tsv", tsv)
system2("awk", paste("-v FS='\t' -v OFS='\t' ' {print $1, $6, $8, $ 9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $27, $72}' ", tsv, " > ", tsv_stripped, sep=""))
tsv <- tsv_stripped
}


message("Building BOLD based refference database using: ", tsv)
data <- read.csv(tsv, sep="\t", stringsAsFactors=F)

if(!exists("tsv_stripped")){
message("Removing unneeded coulumns from table.")
data <- data[,c(1, 6, 8:24, 27, 72)]
}


if(!GenBank){
GB <- data$institution_storing=="Mined from GenBank, NCBI"
message("Excluding sequences mined from Genbank: ", sum(GB), "/", nrow(data), " (", round(sum(GB)/nrow(data)*100, 2), "% removed).")
data <- data[!GB,]
}

# rm GB info
data <- data[,c(-2)]



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
dataID <- data[,grep("ID", names(data))]
taxon <- rep(0, nrow(data))

taxon[which(!is.na(dataID$phylum_taxID))] <- 1
taxon[which(!is.na(dataID$class_taxID))] <- 2
taxon[which(!is.na(dataID$order_taxID))] <- 3
taxon[which(!is.na(dataID$family_taxID))] <- 4
taxon[which(!is.na(dataID$subfamily_taxID))] <- 5
taxon[which(!is.na(dataID$genus_taxID))] <- 6
taxon[which(!is.na(dataID$species_taxID))] <- 7
taxon[which(!is.na(dataID$subspecies_taxID))] <- 8


message("Preparing taxonomy reference files...")
if(is.na(savefasta)){
savefasta <- paste("DB_", sub("(.*)\\..*$", "\\1", tsv), ".fasta", sep="")
}

if(is.na(savetaxonomy)){
savetaxonomy <- paste("DB_", sub("(.*)\\..*$", "\\1", tsv), "_taxonomy.csv", sep="")
}

# make sequence ID
taxID <- as.vector(unlist(t(dataID)))[taxon+seq(0, nrow(dataID)*8-1, 8)]
ID <- paste(data$processid, "_", data$bin_uri, "_taxID=BOLD:", taxID, sep="")


# make taxonomy table
taxonomyDB <- data.frame("taxID"=paste("BOLD:", taxID, sep=""), "phylum"=data$phylum_name, "class"=data$class_name, "order"=data$order_name, "family"=data$family_name, "genus"=data$genus_name, "species"=data$species_name, "subspecies"=data$subspecies_name, "taxRef"=data$identification_reference, stringsAsFactors=F)

# remove taxID duplicates
taxonomyDB <- taxonomyDB[!duplicated(taxonomyDB$taxID),]
message(nrow(taxonomyDB), " unique taxonomic IDs detected (from ", nrow(data), " specimens).")
taxonomyDB$taxRef <- gsub("[()]", "", taxonomyDB$taxRef)

write.csv(taxonomyDB, savetaxonomy, row.names=F)

cat(paste(paste(">", ID, "\n", sep=""), data$nucleotides, "\n", sep=""), file= savefasta, sep="")


message("Module complete! use ", savefasta, " as a refference database, in combination with the taxonomy table ", savetaxonomy)
}

