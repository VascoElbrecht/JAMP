# U_cluster_otus v0.1

Map2ref <- function(files="latest", refDB=NULL, id=0.97, strand="plus", onlykeephits=F, filter=0.01){


folder <- Core(module="Map2ref")
cat(file="log.txt", c("Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files[1]=="latest"){
source(paste(folder, "/robots.txt", sep=""))
files <- list.files(paste(last_data, "/_data", sep=""), full.names=T)
}





# Dereplicate files using USEARCH
dir.create(paste(folder, "/_data/1_derep", sep=""))

new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub("_PE.*", "_PE_derep.fasta", new_names)
new_names <- sub("_data", "_data/1_derep", new_names)
new_names <- paste(folder, "/", new_names, sep="")

cmd <- paste("-fastx_uniques \"", files, "\" -fastaout \"", new_names, "\" -sizeout",  sep="")

temp <- paste(length(files), " files are dereplicated (incl. singletons):", sep="")
cat(file="log.txt", temp , append=T, sep="\n")
message(temp)


temp <- new_names
for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
meep <- sub(".*_data/(.*)", "\\1", temp[i])
cat(file="log.txt", meep, append=T, sep="\n")
cat(file=paste(folder, "/_stats/1_derep_logs.txt", sep=""), meep, A, "\n", append=T, sep="\n")
message(meep)
}

# add saving unused sequences as extra files!!!


# Mapp to refDB
dir.create(paste(folder, "/_data/2_mapping", sep=""))
dir.create(paste(folder, "/_stats/map_logs", sep=""))

blast_names <- sub("1_derep", "2_mapping", new_names)
blast_names <- sub("_derep.fasta", ".txt", blast_names)

log_names <- sub("_data/2_mapping/", "_stats/map_logs/", blast_names)


cmd <- paste("-usearch_global ", new_names, " -db \"", refDB, "\" -strand ", strand, " -id 0.97 -blast6out \"", blast_names, "\" -maxhits 1", sep="")


temp <- paste("Comparing ", length(cmd)," files with dereplicated reads (incl. singletons) against refDB: \"", sub(".*/(.*)", "\\1", refDB), "\" using \"usearch_global\" and Usearch.\n", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

exp <- NULL
temp <- new_names
for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
cat(file= log_names[i], paste("usearch ", cmd[i], sep=""), "\n", A, append=F, sep="\n")

meep <- sub("_data/.*/(.*)", "\\1", temp[i])
pass <- sub(".*, (.*)% matched\r", "\\1", A[grep("matched\r", A)])
exp <- rbind(exp, c(meep, pass))
glumanda <- paste(meep," - ", pass, "% reads matched", sep="")
cat(file="log.txt", glumanda, append=T, sep="\n")
message(glumanda)
}





# condensing hit tables!
files <- blast_names

tab <- c("NULL")
tab <- as.data.frame(tab, stringsAsFactors=F)
names(tab) <- "ID"

for (i in 1:length(files)){
data <- read.csv(files[i], sep="\t", header=F, stringsAsFactors=F)

names(data) <- c("query", "ref", "ident", "length", "mism", "gap", "qstart", "qend", "target_s", "target_e", "e.value", "bitscore")

data <- data[,c(-11,-12)]

data <- cbind(data, "abund"=as.numeric(sub(".*size=(.*);", "\\1", data$query)), stringsAsFactors=F)

#head(data)

temp <- aggregate(data$abund, by=list(data$ref), FUN="sum")
tab <- merge(tab , temp, by.x="ID", by.y="Group.1", all=T, sort=T)
names(tab)[i+1] <- sub(".*2_mapping/(.*).txt", "\\1", files[i])
}

tab <- tab[-1,] # remove NULL entry in the beginning
tab[is.na(tab)] <- 0



sequ <- read.fasta(refDB, forceDNAtolower=F, as.string=T)

# KEEP ONLY HITS OR KEEP ALL IN DB
if(onlykeephits){
temp2 <- match(attr(sequ, "name"), tab$ID)
tab2 <- cbind(tab, "sequ"=as.vector(unlist(sequ[temp2])))

} else {
temp2 <- data.frame("ID"=attr(sequ, "name"))

tab2 <- merge(tab, temp2, "ID", all=T)
tab2[is.na(tab2)] <- 0

#add sequences
temp2 <- match(attr(sequ, "name"), tab2$ID)
tab2 <- cbind(tab2, "sequ"=as.vector(unlist(sequ[temp2])))

}



# filter to relative abundance
rel_abund <- tab2

sampleabundance <- colSums(rel_abund[,2:(ncol(rel_abund)-1)])
for (i in 2:(ncol(rel_abund)-1)){
rel_abund[i] <- rel_abund[i]/sampleabundance[i-1]*100
rel_abund[i][rel_abund[i]<filter] <- 0
}
# write rel abundance tab
rel_abund <- rel_abund[order(rowSums(rel_abund[-c(1, ncol(rel_abund))]), decreasing=T),] # sort table by row sums
write.csv(file=paste(folder, "/3_rel_abundnace_ZEROs.csv", sep=""), rel_abund, row.names=F)




# write RAW table
tab2 <- tab2[order(rowSums(tab2[-c(1, ncol(tab2))]), decreasing=T),] # sort table by row sums

write.csv(file=paste(folder, "/3_Raw_hit_table.csv", sep=""), tab2, row.names=F)





# make plots!

pdf(paste(folder, "/rel_zero2.pdf", sep=""), height=(nrow(rel_abund)+20)/10, width=(ncol(rel_abund)-1)/2)

temp_heat <- rel_abund[,2:(ncol(rel_abund)-1)]
row.names(temp_heat) <- rel_abund[,1]

OTU_heatmap(temp_heat, abundance=F, col=rev(c("Red", "Orange", "Yellow", "Black")))
dev.off()








if(F){



temp <- "\n\nOTU table generated (including OTU sequences): 3_Raw_OTU_table.csv"
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

exp2 <- data.frame("ID"=exp[,1], "Abundance"=colSums(tab[-1]), "pct_pass"=exp[,2], row.names=1:length(exp[,1]))

write.csv(exp2, file=paste(folder, "/_stats/3_pct_matched.csv", sep=""))

#### end raw data table

### make abundance filtering
if(!is.na(filter)){
#tab2 <- read.csv(file="3_Raw_OTU_table.csv", stringsAsFactors=F)

start <- which(names(tab2)=="ID")+1
stop <- which(names(tab2)=="sequ")-1

temp <- tab2[, start:stop]

meep <- paste("Discarding OTUs with below ", filter, "% abundance across at least ", filterN, " out of ", ncol(temp), " samples.", sep="")
message(meep)
cat(file="log.txt", meep, append=T, sep="\n")

temp2 <- temp
sampleabundance <- colSums(temp)
for (i in 1:ncol(temp)){
temp2[i] <- temp[i]/sampleabundance[i]*100
}

# subset OTUs
subset <- rowSums(temp2>=filter)
subset2 <- subset >= filterN

# reporting
meep <- paste("Discarded OTUs: ", sum(!subset2)," out of ",  length(subset2), " discarded (", round(100-sum(subset2)/length(subset2)*100, 2), "%)", sep="")
message(meep)
cat(file="log.txt", meep, append=T, sep="\n")

#write subsetted OTU table

exp <- tab2[subset2,]
exp <- rbind(exp, NA)



exp[nrow(exp), start:stop] <- colSums(tab2[!subset2, start:stop])
exp$ID[nrow(exp)] <- paste("below_", filter, sep="")
exp$sort[nrow(exp)] <- exp$sort[nrow(exp)-1]

tail(exp)


# make folder 
dir.create(paste(folder, "/_data/5_subset/", sep=""))


write.csv(exp, file=paste(folder, "/_data/5_subset/5_OTU_sub_", filter, "_not_rematched.csv", sep=""), row.names=F)

OTU_sub_filename <- paste(folder, "/_data/5_subset/5_OTU_sub_", filter, ".fasta", sep="")
write.fasta(as.list(exp$sequ[-nrow(exp)]), exp$ID[-nrow(exp)], file.out=OTU_sub_filename)

# remapping of reads against subsetted OTUs!
# compare against refernce sequences (including singletons)

dir.create(paste(folder, "/_data/5_subset/usearch_global", sep=""))
dir.create(paste(folder, "/_stats/5_subset/", sep=""))

#new_names <- list.files("_data/1_derep_inc_singletons", full.names=T)
blast_names <- sub("_PE_derep.*.fasta", ".txt", new_names)
blast_names <- sub("1_derep_inc_singletons", "5_subset/usearch_global", blast_names)
blast_names <- sub("1_derep_unoise3", "5_subset/usearch_global", blast_names)
log_names <- sub("_data/", "_stats/", blast_names)
log_names <- sub("/usearch_global", "", log_names)


cmd <- paste("-usearch_global ", new_names, " -db ", OTU_sub_filename, " -strand plus -id 0.97 -blast6out ", blast_names, " -maxhits 1", sep="")


temp <- paste("\n\nRemapping ", length(cmd)," files (incl. singletons) against subsetted OTUs \"", OTU_sub_filename, "\" using \"usearch_global\".\n", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

exp <- NULL
temp <- new_names
for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
cat(file= log_names[i], paste("usearch ", cmd[i], sep=""), "\n", A, append=F, sep="\n")

meep <- sub(".*singletons/(.*)", "\\1", temp[i])
pass <- sub(".*, (.*)% matched\r", "\\1", A[grep("matched\r", A)])
exp <- rbind(exp, c(meep, pass))
glumanda <- paste(meep," - ", pass, "% reads matched", sep="")
cat(file="log.txt", glumanda, append=T, sep="\n")
message(glumanda)
}

#Remapping end

# Writing subsetted & remapped OTU table!
# Write raw data OTU table! incl OTU sequences
files <- blast_names

tab <- c("NULL")
tab <- as.data.frame(tab, stringsAsFactors=F)
names(tab) <- "ID"

for (i in 1:length(files)){
data <- read.csv(files[i], sep="\t", header=F, stringsAsFactors=F)

names(data) <- c("query", "otu", "ident", "length", "mism", "gap", "qstart", "qend", "target_s", "target_e", "e.value", "bitscore")

data <- data[,c(-11,-12)]

data <- cbind(data, "abund"=as.numeric(sub(".*size=(.*);", "\\1", data$query)), "otu_no"=sub("(.*);size.*", "\\1", data$otu), stringsAsFactors=F)

head(data)

temp <- aggregate(data$abund, by=list(data$otu_no), FUN="sum")
tab <- merge(tab , temp, by.x="ID", by.y="Group.1", all=T, sort=T)
names(tab)[i+1] <- sub(".*/(.*)_PE.*", "\\1", files[i])
}

head(tab)

tab <- tab[-1,] # remove NULL entry in the beginning
tab[is.na(tab)] <- 0

mrew <- tab$ID
mrew <- as.numeric(gsub("OTU_(.*)", "\\1", mrew))
tab <- tab[order(as.numeric(mrew)),]

# add sequences!

sequ <- read.fasta(OTU_sub_filename, forceDNAtolower=F, as.string=T)

tab2 <- cbind("sort"=as.numeric(sub("OTU_", "", tab[,1])), tab, "sequ"=unlist(sequ))

names(tab2) <- sub(".*_data/5_subset/usearch_global/(.*).txt", "\\1", names(tab2)) # SUBSET HERE

# add below OTUs

subSums <- read.csv(paste(folder, "/3_Raw_OTU_table.csv", sep=""))
subSums <- as.vector(c(subSums[nrow(subSums)-1, 1]+1, paste("below_", filter, sep=""), colSums(subSums[,-c(1,2, ncol(subSums))])-colSums(tab2[-c(1,2, ncol(tab2))]), NA))


tab3 <- rbind(tab2, subSums)


# set values 2 zero

expZERO <- tab3
#d <- 8
for (d in 3:(ncol(expZERO)-1)){

ZERO <- as.numeric(expZERO[,d])/sampleabundance[d-2]*100>=filter
discarded <- sum(as.numeric(expZERO[-nrow(expZERO),d][!ZERO[-nrow(expZERO)]]))
expZERO[-nrow(expZERO),d][!ZERO[-nrow(expZERO)]] <- 0 # set zero

expZERO[nrow(expZERO),d] <- as.numeric(expZERO[nrow(expZERO),d])+ discarded # add counts to discarded

}

#relative abundance table
expZEROrel <- expZERO
#d <- 8
for (d in 3:(ncol(expZERO)-1)){

expZEROrel[,d] <- as.numeric(expZERO[,d])/sampleabundance[d-2]*100


}



write.csv(file=paste(folder, "/5_OTU_table_", filter,".csv", sep=""), tab3, row.names=F)
write.csv(file=paste(folder, "/5_OTU_table_", filter,"_ZERO.csv", sep=""), expZERO, row.names=F)
write.csv(file=paste(folder, "/5_OTU_table_", filter,"_ZERO_rel.csv", sep=""), expZEROrel, row.names=F)

temp <- paste("\n\nSubsetted OTU table generated (", filter, "% abundance in at least ", filterN," sample): ", sub(paste(folder, "/_data/5_subset/", sep=""), "", OTU_sub_filename), sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")
#### end subsetted OTU table

} # end subsetting



# make plots

RAW <- read.csv(paste(folder, "/3_Raw_OTU_table.csv", sep=""), stringsAsFactors=F)
KEEP <- read.csv(paste(folder, "/5_OTU_table_", filter,".csv", sep=""), stringsAsFactors=F)

higlight <- rev(!RAW$ID %in% KEEP$ID)

pdf(paste(folder, "/_stats/OTU_plot_3_RAW.pdf", sep=""), height=(nrow(RAW)+20)/10, width=(ncol(RAW)-1)/2)

OTU_heatmap(paste(folder, "/3_Raw_OTU_table.csv", sep=""), abundance=T)


pos <- ncol(RAW) - 2.5
for (i in 1:length(higlight)){

if(higlight[i]){
rect(pos, i-0.5, pos+1, i+0.5, col="Lightgray", border=F)
text(pos+0.1, i, rev(RAW$ID)[i], adj=0, cex=0.5)
}

}
dev.off()

OTU_heatmap(paste(folder, "/5_OTU_table_", filter,".csv", sep=""), out=paste(folder, "/_stats/OTU_plot_5_", filter, ".pdf", sep=""), abundance=T)
OTU_heatmap(paste(folder, "/5_OTU_table_", filter,"_ZERO.csv", sep=""), out=paste(folder, "/_stats/OTU_plot_5_", filter, "_ZERO.pdf", sep=""), abundance=T)
OTU_heatmap(paste(folder, "/5_OTU_table_", filter,"_ZERO_rel.csv", sep=""), out=paste(folder, "/_stats/OTU_plot_5_", filter, "_ZERO_rel.pdf", sep=""), abundance=T, rel=T)

} # temp stop

temp <- "\nModule completed!"
message(temp)
cat(file="log.txt", paste(Sys.time(), "*** Module completed!", "", sep="\n"), append=T, sep="\n")


}

