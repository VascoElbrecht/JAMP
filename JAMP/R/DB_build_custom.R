# 191025 convert BOLD ref data to DB
# taxID not required (can be auto generated)
# ID = sequence ID, required, as well as sequences them self

buildCustomdb <- function(savefasta=NA, savetaxonomy=NA, minlength=500, maxlength=NA, taxIDprefix="JAMP", fasta=NA, ID=NA, phylum=NA, class=NA, order=NA, family=NA, genus=NA, species=NA, subspecies=NA, taxRef=NA){

message("Building refference database using custom data!")

if(is.na(fasta)){
stop(call="No fasta file or vector of sequences is provided. You need sequences if you would like to build a reference database. Function stopped!")
}
if(is.na(ID)){
stop(call="No sequence IDs are provided in \"ID\". Having a unique identifyer for each sequence is reequired. Function stopped!")
}
if(is.na(savefasta)){
stop(call="No file name for saving fasta files is given in \"savefasta\". Please provide a file name to save the fasta files.")
}
if(is.na(savetaxonomy)){
stop(call="No file name for saving taxonomy csv is given in \"savetaxonomy\". Please provide a file name to save the taxonomy csv.")
}
# apply length filtering


# NEED TO BE ADDED



# Build taxonomy

sumtax <- paste(phylum, class, order, family, genus, species, subspecies, taxRef, sep="")
sumtaxU <- unique(sumtax)
sumtaxU <- match(sumtax, sumtaxU)

taxonomy <- data.frame("taxID"=paste(taxIDprefix, sumtaxU, sep=":"), "phylum"=phylum, "class"=class, "order"=order, "family"=family, "genus"=genus, "species"=species, "subspecies"=subspecies, "taxRef"=taxRef, stringsAsFactors=F)

taxonomy <- taxonomy[!duplicated(taxonomy$taxID),]

write.csv(taxonomy, savetaxonomy, row.names=F)

#build fasta


fastaname <- paste(">", ID, "_", taxIDprefix, ":", sumtaxU, sep="")

cat(paste(fastaname, "\n", fasta, sep=""), sep="\n", file=savefasta)



message("Module complete! use ", savefasta, " as a refference database, in combination with the taxonomy table ", savetaxonomy)
}

