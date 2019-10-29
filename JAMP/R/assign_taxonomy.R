# 191027 - assign taxonomy

assign_taxonomy <- function(query=NA, reffasta=NA, reftaxonomy=NA, id=0.8, maxaccepts=32, maxrejects=256, strands="both", threads=12, query_cov=1, exactmatch=T, bracket=1, MM=c(0.85, 0.90, 0.95, 0.98), exe="usearch"){

# check that variables are not empty
if(is.na(query)){
stop(call="No query fasta file provided! Make sure to give path to the fasta file in \"query\"! Function stopped!")
}
if(is.na(reffasta)){
stop(call="No refference fasta file provided! Make sure to give path to the fasta file in \"reffasta\"! Function stopped!")
}
if(is.na(reftaxonomy)){
stop(call="No taxonomy table provided! Make sure to give path to the taxonomy csv in \"reftaxonomy\"! Function stopped!")
}

# check that imput files exist
if(!file.exists(query)){
stop(call=c("File \"", query, "\" does not exist. Please check file path! Function stopped!"))
}
if(!file.exists(reffasta)){
stop(call=c("File \"", reffasta, "\" does not exist. Please check file path! Function stopped!"))
}
if(!file.exists(reftaxonomy)){
stop(call=c("File \"", reftaxonomy, "\" does not exist. Please check file path! Function stopped!"))
}

fasta <- readLines(query) # read for later IDs

hitfile <- paste(sub("(.*)\\..*$", "\\1", query), "_hits.txt", sep="")

system2(exe,  paste("-usearch_global ", query, " -db ", reffasta," -id ", id, " -blast6out ", hitfile, " -strand ", strands, " -maxaccepts ", maxaccepts, " -maxrejects ", maxrejects, " -query_cov ", query_cov, " -maxhits 100 -threads 16", sep=""), stdout="", stderr="")


# reformatting 
hits <- read.csv(hitfile, stringsAsFactors=F, sep="\t", header=F)

names(hits) <- c("Query", "DB", "id", "alignlength", "missmatches", "gaps", "qstart", "qend", "DBstart", "DBend", "eval", "bitscore")
hits <- hits[,-c(11,12)]
hits <- data.frame(hits, "seqID"=sub("(.*)_taxID=.*", "\\1", hits$DB), "taxID"=sub(".*taxID=(.*)", "\\1", hits$DB), stringsAsFactors=F)

hits <- hits[c(1, 11, 3:10, 12)]


# processing individual queries
subset <- hits[0,]
query <- unique(hits$Query)

for(i in 1:length(query)){
temp <- hits[hits$Query==query[i],]
if(temp$id[1]==100&exactmatch){
subset <- rbind(subset, temp[temp$id==100,])
} else {
subset <- rbind(subset, temp[temp$id>=c(max(temp$id)-bracket),])
}
}

# save top hits
taxonomy <- read.csv(reftaxonomy, stringsAsFactors=F)
tophits <- merge(subset, taxonomy, by="taxID")

tophits <- tophits[order(tophits$id, decreasing=T),]
tophits <- tophits[order(as.numeric(sub(".*_", "", tophits$Query)), decreasing=F),]

tophits <- tophits[,c(2, 12:19, 4, 1, 3, 5:11)]
write.csv(tophits, sub("hits.txt", "tophits.csv", hitfile), row.names=F)


# extract main top hit and number of taxa matched
meep <- rep(NA, length(query))
taxassigned <- data.frame("Query"=meep, "id"=meep, "taxID"=meep, "taxaMatched"=meep)

for(i in 1:length(query)){
temp <- subset[subset$Query==query[i],]
temp2 <- as.data.frame(table(temp$taxID), stringsAsFactors=F)
taxassigned[i,] <- c(as.vector(unlist(temp[which(max(temp2$Freq)==temp2$Freq)[1], c(1, 3, 11)])), nrow(temp2))
}

# add taxonomic information + filter %tage wise


taxassigned <- merge(taxassigned, taxonomy, by="taxID")

taxassigned$id <- as.numeric(taxassigned$id)

taxassigned$taxRef[taxassigned$id<MM[4]*100] <- NA
taxassigned$subspecies[taxassigned$id<MM[4]*100] <- NA
taxassigned$species[taxassigned$id<MM[4]*100] <- NA
taxassigned$genus[taxassigned$id<MM[3]*100] <- NA
taxassigned$family[taxassigned$id<MM[2]*100] <- NA
taxassigned$order[taxassigned$id<MM[1]*100] <- NA


# add no hits in and match sorting order
fasta <- fasta[grep(">", fasta)]
fasta <- sub("^>", "", fasta)

fasta <- data.frame(fasta, stringsAsFactors=F)

exp <- merge(fasta, taxassigned, by.x="fasta", by.y="Query", all=T, sort=F)
exp <- exp[order(as.numeric(sub(".*_", "", exp$fasta))),]

exp <- exp[,c(1, 5:12, 3, 4)]
write.csv(exp, sub("hits.txt", "HitTable.csv", hitfile), row.names=F)


}



