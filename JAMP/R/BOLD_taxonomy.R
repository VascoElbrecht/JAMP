# Bold taxonomy

Bold_taxonomy <- function(file=NA, folder="", MM=c(0.98, 0.95, 0.90, 0.85)){

oldwd <- getwd()
setwd(folder)
getwd()
#cat(file="../log.txt", c("\n","Version v0.1", "\n"), append=T, sep="\n")
message("Getting taxonomy via BOLD api.")

data <- read.csv(file, stringsAsFactors=F)
sequ <- which(!is.na(data$sequ)) # sequences exist

dir.create("_BOLD_data", showWarnings = FALSE)
setwd("_BOLD_data")
dir.create("_stats", showWarnings = FALSE)

#based on ID_engine_bp_v2_anna_NovemberR.R
# query sequences to bold API

# check if data already downloaded

sequ_inital_DL <- sequ

# skipp if already downloaded!
temp_files <- paste(data$ID[sequ_inital_DL], ".csv", sep="")
sequ_inital_DL <- sequ_inital_DL[!temp_files%in%list.files("_stats", "csv")]

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

write.table(hittable, file=paste("_stats/", data$ID[i], ".csv", collapse="", sep=""), sep="\t")
message(i)

}


i <- 1
new_taxonomy <- NULL
for (i in 1:length(temp_files)){

data_bold <- read.csv(paste("_stats/", temp_files[i], sep=""), sep="\t", stringsAsFactors=F)

data_bold <- data.frame("ID"=data_bold$ID, "taxitenti"=data_bold$taxonomicidentification, "similarity"=data_bold$similarity, "order"=data_bold$order, "family"=data_bold$family, "genus"=data_bold$genus, "species"=data_bold$species, stringsAsFactors=F)
#data_bold <- data_bold[order(data_bold$similarity, decreasing=T),]

data_bold$species[grep("sp.", data_bold$species)] <- data_bold$species[grep("sp.", data_bold$species)] <- NA

data_bold$species[!data_bold$similarity>= MM[1]] <- NA
data_bold$genus[!data_bold$similarity>= MM[2]] <- NA
data_bold$family[!data_bold$similarity>= MM[3]] <- NA
data_bold$order[!data_bold$similarity>= MM[4]] <- NA

data_bold[,4:7][is.na(data_bold[,4:7])] <- "NA"

# species
tab_species <- as.data.frame(table(data_bold$species[data_bold$similarity>= MM[1]]), stringsAsFactors=F)
tab_species <- tab_species[order(tab_species$Freq, decreasing=T),]
adj_species <- paste(tab_species[1,1], " ", tab_species$Freq[1], "/", sum(tab_species$Freq), sep="")

# genus
tab_genus <- as.data.frame(table(data_bold$genus[data_bold$similarity>= MM[2]]), stringsAsFactors=F)
tab_genus <- tab_genus[order(tab_genus$Freq, decreasing=T),]
adj_genus <- paste(tab_genus[1,1], " ", tab_genus$Freq[1], "/", sum(tab_genus$Freq), sep="")

# family
tab_family <- as.data.frame(table(data_bold$family[data_bold$similarity>= MM[3]]), stringsAsFactors=F)
tab_family <- tab_family[order(tab_family$Freq, decreasing=T),]
adj_family <- paste(tab_family[1,1], " ", tab_family$Freq[1], "/", sum(tab_family$Freq), sep="")

# order
tab_order <- as.data.frame(table(data_bold$order[data_bold$similarity>= MM[4]]), stringsAsFactors=F)
tab_order <- tab_order[order(tab_order$Freq, decreasing=T),]
adj_order <- paste(tab_order[1,1], " ", tab_order$Freq[1], "/", sum(tab_order$Freq), sep="")

# to do:
# Add taxonomic level
# add taxon name!

new_taxonomy <- rbind(new_taxonomy, c(temp_files[i] , adj_order, adj_family, adj_genus, adj_species, data_bold$similarity[1]))

}

new_taxonomy <- data.frame(new_taxonomy, stringsAsFactors=F)
names(new_taxonomy) <- c("BOLD_OTU", "BOLD_Order", "BOLD_Family", "BOLD_Genus", "BOLD_Species", "BOLD_best_match")
new_taxonomy <- rbind(new_taxonomy, NA)


export <- cbind(data, new_taxonomy)

write.csv(export, file=paste("../", sub(".csv", "_BOLD.csv", file), sep=""), row.names=F)

setwd(oldwd)
}
