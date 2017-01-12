# Bold taxonomy

Bold_taxonomy <- function(file=NA, folder=""){

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
temp <- paste(data$ID[sequ_inital_DL], ".csv", sep="")
sequ_inital_DL <- sequ_inital_DL[!temp%in%list.files("_stats", "csv")]

for (i in sequ_inital_DL){

hittable <- NULL

hittable <- bold_identify(data$sequ[i], db="COX1")
revcomp <- paste(rev(comp(strsplit(data$sequ[i], "")[[1]])), collapse="")
if(is.null(hittable[[1]])){ # no hits, check rev comp
hittable <- bold_identify(revcomp, db="COX1")[[1]]} #else { # hits FW, add rev com hits
	#temp <- bold_identify(revcomp, db="COX1")[[1]]
	#if(!is.null(temp[[1]])){hittable <- rbind(temp[[1]], hittable)}
#}
# write bold results in file!

write.table(hittable, file=paste("_stats/", data$ID[i], ".csv", collapse="", sep=""), sep="\t")
message(i)

}

# download additional taxonomic information from BIN information
if(!is.null(which(!sequ%in%sequ_inital_DL))){
message(paste("OTUs Taxonomy identified, thus skipped:\n", paste(which(!sequ%in%sequ_inital_DL), collapse=", "), sep=""))
}


# read in hit tables and steal taxonomic information from bold! 
otu_files <- list.files("_stats", full.names=T, pattern=".csv")


sequ_inital_DL <- sequ

# skipp if already downloaded!
temp <- paste("_stats/", data$ID[sequ_inital_DL], ".csv", sep="")
sequ_inital_DL <- sequ_inital_DL[!temp%in%list.files("_stats", "_hacked.txt")]
message(paste("hacked taxonomy availabl, thus skipped for OTUs:\n", which(!sequ%in%sequ_inital_DL), sep=""))


for(s in sequ_inital_DL){ # adjust numbers here after crash

if(readLines(otu_files[s])[1]!="\"\""){

tab <- read.table(otu_files[s], stringsAsFactors=F)

bins <- bold_specimens(ids=tab$ID)# get bin uri
tab_bins <- merge(tab, bins, by.x="ID", by.y="processid")
tab_bins <- tab_bins[order(tab_bins$similarity, decreasing=T),]


IDtab <- data.frame("ID"=tab_bins$ID, "bin_uri"=tab_bins$bin_uri, "similarity"=tab_bins$similarity, "country"=tab_bins$specimen_country, "orig_tax"=tab_bins$taxonomicidentification, "Order"=NA, "Family"=NA, "Subfamily"=NA, "Genus"=NA, "Species"=NA, stringsAsFactors=F)

unique_bins <- unique(IDtab$bin_uri) # get uinique bins

unique_bins <- unique_bins[unique_bins!=" "] # remove missing data?

if (is.na(unique_bins[1])){}else{

for (j in 1:length(unique_bins)){
html <- htmlParse(paste("http://www.boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=", unique_bins[j], sep=""), encoding="UCS-2LE")

html <- readHTMLTable(html, stringsAsFactors=F) 
length(html)

if (length(html)==2){
html <- readLines(paste("http://www.boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=", unique_bins[j], sep=""), encoding="UCS-2LE", warn=F)
Sys.sleep(1) # adding pause to not get kicked out by bold server

newID <- html[grep("has been synonymized", html)]
newID <- sub(".*clusteruri=(.*)\".*", "\\1", newID)
html <- htmlParse(paste("http://www.boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=", newID, sep=""), encoding="UCS-2LE")
html <- readHTMLTable(html, stringsAsFactors=F) 

}

if(length(html)!=2){
order <- html[[16]]$V3[which(html[[16]]$V2=="Order:")]
IDtab[IDtab$bin_uri==unique_bins[j], 6] <- sub(" \\[.*\\]", "", order)
family <- html[[16]]$V3[which(html[[16]]$V2=="Family:")]
IDtab[IDtab$bin_uri==unique_bins[j], 7] <- sub(" \\[.*\\]", "", family)
subfamily <- html[[16]]$V3[which(html[[16]]$V2=="Subfamily:")]
IDtab[IDtab$bin_uri==unique_bins[j], 8] <- sub(" \\[.*\\]", "", subfamily)
genus <- html[[16]]$V3[which(html[[16]]$V2=="Genus:")]
IDtab[IDtab$bin_uri==unique_bins[j], 9] <- sub(" \\[.*\\]", "", genus)
species <- html[[16]]$V3[which(html[[16]]$V2=="Species:")]
IDtab[IDtab$bin_uri==unique_bins[j], 10] <- sub(" \\[.*\\]", "", species)
}
}
write.table(IDtab, paste(sub(".csv", "", otu_files[s]), "_hacked.txt", sep="", collapse=""), sep="\t")}
}

}









setwd(oldwd)
}
