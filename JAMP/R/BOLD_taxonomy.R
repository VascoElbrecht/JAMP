# Bold taxonomy

Bold_taxonomy <- function(file=NA, folder="", MM=c(0.98, 0.95, 0.90, 0.85)){

#cat(file="../log.txt", c("\n","Version v0.1", "\n"), append=T, sep="\n")
message("Getting taxonomy via BOLD api.")

# read table
data <- read.csv(file, stringsAsFactors=F)
csv <- T

if(ncol(data)==1){ # fasta
data <- read.fasta(file, forceDNAtolower=F, as.string=T)
data <- data.frame("ID"=names(data), sequ=as.vector(unlist(data)), stringsAsFactors=F)
scv<-F
}

sequ <- which(!is.na(data$sequ)) # sequences exist

dir.create("_BOLD_data", showWarnings = FALSE)
dir.create("_BOLD_data/_stats", showWarnings = FALSE)

#based on ID_engine_bp_v2_anna_NovemberR.R
# query sequences to bold API

# check if data already downloaded

sequ_inital_DL <- sequ

# skipp if already downloaded!
temp_files <- paste(data$ID[sequ_inital_DL], ".csv", sep="")
sequ_inital_DL <- sequ_inital_DL[!temp_files%in%list.files("_BOLD_data/_stats", "csv")]

for (i in sequ_inital_DL){

hittable <- NULL

hittable <- bold_identify(data$sequ[i], db="COX1")

if(is.null(hittable[[1]])){ # no hits, check rev comp
revcomp <- paste(rev(comp(strsplit(data$sequ[i], "")[[1]])), collapse="")
hittable <- bold_identify(revcomp, db="COX1")[[1]]}

if(!is.null(hittable)){ # if matched, get taxonomy
hittable <- bold_identify_parents(hittable, wide=T)
hittable <- hittable[[1]]
hittable <- hittable[order(hittable$similarity, decreasing=T),]
}

write.table(hittable, file=paste("_BOLD_data/_stats/", data$ID[i], ".csv", collapse="", sep=""), sep="\t")
message(i)

}


new_taxonomy <- NULL
for (i in 1:length(temp_files)){

if(readLines(paste("_BOLD_data/_stats/", data$ID[i], ".csv", collapse="", sep=""))[1]=="\"\""){new_taxonomy <- rbind(new_taxonomy, c(temp_files[i] , NA, NA, NA, NA, NA, NA, NA))} else {#no hit



data_bold <- read.csv(paste("_BOLD_data/_stats/", temp_files[i], sep=""), sep="\t", stringsAsFactors=F)

if(is.null(data_bold$family)){data_bold <- data.frame(data_bold, "family"=rep(NA, nrow(data_bold)))}
if(is.null(data_bold$genus)){data_bold <- data.frame(data_bold, "genus"=rep(NA, nrow(data_bold)))}
if(is.null(data_bold$species)){data_bold <- data.frame(data_bold, "species"=rep(NA, nrow(data_bold)))}

data_bold <- data.frame("ID"=data_bold$ID, "taxitenti"=data_bold$taxonomicidentification, "similarity"=data_bold$similarity, "order"=data_bold$order, "family"=data_bold$family, "genus"=data_bold$genus, "species"=data_bold$species, stringsAsFactors=F)

data_bold$species[grep("sp.", data_bold$species)] <- data_bold$species[grep("sp.", data_bold$species)] <- NA

data_bold$species[!data_bold$similarity>= MM[1]] <- NA
data_bold$genus[!data_bold$similarity>= MM[2]] <- NA
data_bold$family[!data_bold$similarity>= MM[3]] <- NA
data_bold$order[!data_bold$similarity>= MM[4]] <- NA

data_bold[,4:7][is.na(data_bold[,4:7])] <- "NA"

# species
tab_species <- as.data.frame(table(data_bold$species[data_bold$similarity>= MM[1]]), stringsAsFactors=F)
tab_species <- tab_species[order(tab_species$Freq, decreasing=T),]
if(length(tab_species)==0){adj_species <- "NA"} else {
tab_species <- tab_species[!tab_species$Var1=="NA",]
adj_species <- tab_species[1,1]
}

# genus
tab_genus <- as.data.frame(table(data_bold$genus[data_bold$similarity>= MM[2]]), stringsAsFactors=F)
tab_genus <- tab_genus[order(tab_genus$Freq, decreasing=T),]
if(length(tab_genus)==0){adj_genus <- "NA"} else {
tab_genus <- tab_genus[! tab_genus$Var1=="NA",]
adj_genus <- tab_genus[1,1]
if(!is.null(nrow(tab_species))){
if(nrow(tab_species)>0){
tab_genus <- tab_species[1,]
tab_genus$Var1 <- sub("(.*) .*", "\\1",tab_genus$Var1)
adj_genus <- tab_genus[1,1]}
}
}


# family
tab_family <- as.data.frame(table(data_bold$family[data_bold$similarity>= MM[3]]), stringsAsFactors=F)
tab_family <- tab_family[order(tab_family$Freq, decreasing=T),]
if(length(tab_family)==0){adj_family <- "NA"} else {

tab_family <- tab_family[! tab_family $Var1=="NA",]
adj_family <- tab_family[1,1]
}


# order
tab_order <- as.data.frame(table(data_bold$order[data_bold$similarity>= MM[4]]), stringsAsFactors=F)
tab_order <- tab_order[order(tab_order$Freq, decreasing=T),]
if(length(tab_order)==0){adj_order <- "NA"} else {
tab_order <- tab_order[! tab_order $Var1=="NA",]
adj_order <- tab_order[1,1]
}



# Add taxonomic level
# add taxon name!

bold_taxonomy <- if(!grepl("NA", adj_species)){
bold_level <- "Species"
bold_taxonomy <- tab_species[1,1]
} else if(!grepl("NA", adj_genus)){
bold_level <- "Genus"
bold_taxonomy <- tab_genus[1,1]
} else if(!grepl("NA", adj_family)){
bold_level <- "Family"
bold_taxonomy <- tab_family[1,1]
} else if(!grepl("NA", adj_order)){
bold_level <- "Order"
bold_taxonomy <- tab_order[1,1]} else {
bold_level <- NA
bold_taxonomy <- NA
}


new_taxonomy <- rbind(new_taxonomy, c(temp_files[i] , adj_order, adj_family, adj_genus, adj_species, data_bold$similarity[1], bold_taxonomy, bold_level))

bold_level <- NA
bold_taxonomy <- NA
}

} # end loop

new_taxonomy <- data.frame(new_taxonomy, stringsAsFactors=F)
names(new_taxonomy) <- c("BOLD_OTU", "BOLD_Order", "BOLD_Family", "BOLD_Genus", "BOLD_Species", "BOLD_best_match", "bold_taxon", "bold_level")

if(csv){
new_taxonomy <- rbind(new_taxonomy, NA)
}

export <- cbind(data, new_taxonomy)

write.csv(export, file=paste("_BOLD_data/BOLD_API_taxonomy.csv", sep=""), row.names=F)

}
