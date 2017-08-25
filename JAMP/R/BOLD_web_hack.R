# Bold web hack


Bold_web_hack <- function(file=NA, folder=NULL, MM=c(0.98, 0.95, 0.90, 0.85)){

oldwd <- getwd()
if(!is.null(folder)){
dir.create(folder, showWarnings = FALSE)
setwd(folder)
}
getwd()
#cat(file="../log.txt", c("\n","Version v0.1", "\n"), append=T, sep="\n")

dir.create("_BOLD_web_hack", showWarnings = FALSE)
dir.create("_BOLD_web_hack/OTUs", showWarnings = FALSE)
setwd("_BOLD_web_hack")

if(!is.null(folder)){
data <- readLines(paste("../../", file, sep=""), warn=F)
} else {
data <- readLines(paste("../", file, sep=""), warn=F)
}



OTU_start <- grep("Query: ", data)
#OTU_end <- which("Sampling Sites For Top Hits (>98% Match)"==data)
#OTU_end <- sort(c(which("Unable to match any records in the selected database. "==data), OTU_end))
OTU_end <- c(OTU_start[-1]-1, length(data))


for (i in 1:length(OTU_start)){


temp <- data[OTU_start[i]:OTU_end[i]]


OTU <- sub("Query: (.*) ", "\\1", data[OTU_start[i]])

if(temp[5]=="Unable to match any records in the selected database. "){
temp_tab <- data.frame(Phylum, Class, Order, Family, Genus, Species, Similarity, Status, stringsAsFactors=F)
write.csv(temp_tab, paste("OTUs/", OTU, ".csv", sep=""), row.names=F)
} else {

similarity <- suppressWarnings(as.numeric(temp))
whois <- which(!is.na(similarity))
whois_num <- whois

Similarity <- similarity[whois]
Status <- temp[whois+2]

# get phylum

whois <- c(which(temp =="Phylum\tClass\tOrder\tFamily\tGenus\tSpecies\tSubspecies\tSimilarity (%)\tStatus")+1, whois[-length(whois)]+4)


# Phylum
Phylum <- temp[whois]
Class <- temp[whois+2]
Order <- temp[whois+4]
Family <- temp[whois+6]
Genus <- temp[whois+8]
Species <- temp[whois+10]
Subspecies <- temp[whois+12] # ignore subspecies for now

temp_tab <- data.frame(Phylum, Class, Order, Family, Genus, Species, Similarity, Status, stringsAsFactors=F)

for (k in 1:nrow(temp_tab)){
kill <- which(suppressWarnings(as.numeric(as.vector(unlist(temp_tab[k, 1:6]))))>0)[1]
if(!is.na(kill)){
temp_tab[k, kill:6] <- ""

}
}

write.csv(temp_tab, paste("OTUs/", OTU, ".csv", sep=""), row.names=F)
}
}


# subset taxonomy stuff

files <- sub("Query: (.*) ", "\\1", data[OTU_start])

BOLD_tab <- NULL

i <- 5

for (i in 1:length(files)){
temp <- read.csv(paste("OTUs/", files[i], ".csv", sep=""), stringsAsFactors=F)

if(nrow(temp)==0){ # no hits!
BOLD_tab <- rbind(BOLD_tab, cbind(files[i], temp[1,]))

} else {
temp[temp==""] <- NA

#remove taxa lower than treshhold!
A <- temp$Similarity<MM[1]*100
temp$Species[A] <- NA
A <- temp$Similarity<MM[2]*100
temp$Genus[A] <- NA
A <- temp$Similarity<MM[3]*100
temp$Family[A] <- NA
A <- temp$Similarity<MM[4]*100
temp$Order[A] <- NA

# get taxon
END <- T # stop if written in table
who_sp <- which(temp$Similarity>= MM[1]*100)
who_sp <- which(!is.na(temp$Species[who_sp])) #rm NA

if(length(who_sp)>0){ # if species present
tab_species <- as.data.frame(table(temp$Species[who_sp]), stringsAsFactors=F)
tab_species <- tab_species[order(tab_species$Freq, decreasing=T),]
tab_species <- if(is.vector(tab_species)){tab_species[1]} else{tab_species[1,1]}

exp <- temp[who_sp[temp$Species[who_sp]==tab_species][1],]
BOLD_tab <- rbind(BOLD_tab, cbind(files[i], exp))
END <- F
}

who_gen <- which(temp$Similarity>= MM[2]*100)
who_gen <- which(!is.na(temp$Genus[who_gen])) #rm NA

if(length(who_gen)>0 & END){

tab_genus <- as.data.frame(table(temp$Genus[who_gen]), stringsAsFactors=F)
tab_genus <- tab_genus[order(tab_genus$Freq, decreasing=T),]
tab_genus <- if(is.vector(tab_genus)){tab_genus[1]} else{tab_genus[1,1]}

exp <- temp[who_gen[temp$Genus[who_gen]==tab_genus][1],]
BOLD_tab <- rbind(BOLD_tab, cbind(files[i], exp))
END <- F

}

who_fam <- which(temp$Similarity>= MM[3]*100)
who_fam <- which(!is.na(temp$Family[who_fam])) #rm NA

if(length(who_fam)>0 & END){

tab_family <- as.data.frame(table(temp$Family[who_fam]), stringsAsFactors=F)
tab_family <- tab_family[order(tab_family $Freq, decreasing=T),]
tab_family <- if(is.vector(tab_family)){tab_family[1]} else{tab_family[1,1]}


exp <- temp[who_fam[temp$Family[who_fam]==tab_family][1],]
exp[5:6] <- NA
BOLD_tab <- rbind(BOLD_tab, cbind(files[i], exp))
END <- F
}

who_ord <- which(temp$Similarity>= MM[4]*100)
who_ord <- which(!is.na(temp$Order[who_ord])) #rm NA

if(length(who_fam)>0 & END){
	
tab_order <- as.data.frame(table(temp$Order[who_ord]), stringsAsFactors=F)
tab_order <- tab_order[order(tab_order$Freq, decreasing=T),]
tab_order <- if(is.vector(tab_order)){tab_order[1]} else{tab_order[1,1]}

exp <- temp[who_ord[temp$Order[who_ord]==tab_order][1],]
exp[4:6] <- NA
BOLD_tab <- rbind(BOLD_tab, cbind(files[i], exp))
END <- F
}

if(END){
exp <- temp[1,]
exp[3:6] <- NA
BOLD_tab <- rbind(BOLD_tab, cbind(files[i], temp[1,]))
}
} # "skipp no hits"
} # taxa tab loop end


BOLD_tab2 <- data.frame(BOLD_tab, stringsAsFactors=F)
names(BOLD_tab2)[1] <- "OTU_ID"


write.csv(BOLD_tab2, file=sub(".txt", "_taxonomy.csv", file), row.names=F)

setwd(oldwd)
}


