# U_cluster_otus v0.1

Map2ref <- function(files="latest", refDB=NULL, id=0.97, strand="plus", onlykeephits=T, filter=0.01, maxaccepts=1, maxrejects=32, exe="usearch", heatmap=T, delete_data=T, JV=F){



A <- system2(exe, stdout=T)
version <- as.numeric(sub("usearch v(.+)\\.+.*\\..*_.*", "\\1", A[1]))


folder <- Core(module="Map2ref", delete_data=delete_data)
cat(file="log.txt", c("Version v0.1", "\n"), append=T, sep="\n")
message(" ")

files_to_delete <- NULL

if (files[1]=="latest"){
source(paste(folder, "/robots.txt", sep=""))
files <- list.files(paste(last_data, "/_data", sep=""), full.names=T)
}



# check for empty files!
empty <- !file.info(files)$size>0

if(sum(empty)>0){
temp <- paste("WARNING! ", sum(empty), " file", if(sum(empty>1)){"s do"} else {" does"}, " NOT contain sequences and are not included in clustering:", sep="", "\n", paste(files[empty], collapse="\n"), "\n")
message(temp)
cat(file="log.txt", temp , append=T, sep="\n")
}
excluded <- NULL
excluded <- files[empty]




# Dereplicate files using USEARCH
dir.create(paste(folder, "/_data/1_derep", sep=""))

new_names <- sub(".*(/.*)", "_data\\1", files[!empty])
new_names <- sub("_PE.*", "_PE_derep.fasta", new_names)
new_names <- sub("_data", "_data/1_derep", new_names)
new_names <- paste(folder, "/", new_names, sep="")

cmd <- paste(if(version<9){"-derep_fulllength"}else{"-fastx_uniques"}, " \"", files[!empty], "\" -fastaout \"", new_names, "\" -sizeout -sizein",  sep="")

files_to_delete <- c(files_to_delete, new_names)

temp <- paste(length(cmd), " files are dereplicated (incl. singletons):", sep="")
cat(file="log.txt", temp , append=T, sep="\n")
message(temp)


temp <- new_names
sequ_in <- NULL
for (i in 1:length(cmd)){
A <- system2(exe, cmd[i], stdout=T, stderr=T)
meep <- sub(".*_data/(.*)", "\\1", temp[i])
cat(file="log.txt", meep, append=T, sep="\n")
cat(file=paste(folder, "/_stats/1_derep_logs.txt", sep=""), meep, A, "\n", append=T, sep="\n")
message(meep)

sequ_in[i] <- as.numeric(sub(".* (.*) seqs.*", "\\1", A[grep("seqs", A)]))
}

temp <- sub(".*/(.*)", "\\1", files)
expTab <- data.frame("ID"=sub("(.*)_PE.*", "\\1", temp), stringsAsFactors=F, "sequ_in"=0)
expTab$sequ_in[!empty] <- sequ_in
temp <- suppressWarnings(as.numeric(substr(expTab$ID,1,1)))
temp <- which(!is.na(temp))

for(i in temp){
expTab$ID[i] <- paste("X", expTab$ID[i], sep="")
}

expTab <- expTab[order(expTab$ID),]

# add saving unused sequences as extra files!!!
# Mapp to refDB
dir.create(paste(folder, "/_data/2_mapping", sep=""))
dir.create(paste(folder, "/_stats/map_logs", sep=""))
dir.create(paste(folder, "/_data/3_nohit_fasta", sep=""))



blast_names <- sub("1_derep", "2_mapping", new_names)
blast_names <- sub("_derep.fasta", ".txt", blast_names)

nohit <- sub("1_derep", "3_nohit_fasta", new_names)


log_names <- sub("_data/2_mapping/", "_stats/map_logs/", blast_names)

if(!JV){
cmd <- paste("-usearch_global ", new_names, " -db \"", refDB, "\" -strand ", strand, " -id ", id, " -blast6out \"", blast_names, "\" -maxhits 1", " -notmatched \"", nohit, "\" -maxaccepts ", maxaccepts, " -maxrejects ", maxrejects, sep="")
} else { #legacy version, remove at some point
cmd <- paste("-usearch_global ", new_names, " -db \"", refDB, "\" -strand ", strand, " -id ", id, " -blast6out \"", blast_names, "\" -notmatched \"", nohit, "\" -maxaccepts ", maxaccepts, " -maxrejects ", maxrejects, sep="")
}
files_to_delete <- c(files_to_delete, blast_names)

temp <- paste("Comparing ", length(cmd)," files with dereplicated reads (incl. singletons) against refDB: \"", sub(".*/(.*)", "\\1", refDB), "\" using \"usearch_global\" and Usearch. Minimum identity (id) is ", id, ".\n", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

exp <- NULL
temp <- new_names
for (i in 1:length(cmd)){
A <- system2(exe, cmd[i], stdout=T, stderr=T)
cat(file= log_names[i], paste("usearch ", cmd[i], sep=""), "\n", A, append=F, sep="\n")

meep <- sub("_data/.*/(.*)", "\\1", temp[i])
pass <- sub(".*, (.*)% matched\r", "\\1", A[grep("matched\r", A)])
exp <- rbind(exp, c(meep, pass))
glumanda <- paste(meep," - ", pass, "% reads matched", sep="")
cat(file="log.txt", glumanda, append=T, sep="\n")
message(glumanda)
}



# checking for empty files

# check for empty files
empty <- !file.info(blast_names)$size>0

if(sum(empty)>0){
temp <- paste("\nWARNING! ", sum(empty), " file", if(sum(empty>1)){"s do"} else {" does"}, " not contain any reads matching the reference sequences!", sep="", "\n", paste(files[empty], collapse="\n"))
message(temp)
cat(file="log.txt", temp , append=T, sep="\n")
}

excluded2 <- c(excluded, blast_names[empty])



# condensing hit tables!
files <- blast_names[!empty]

tab <- c("NULL")
tab <- as.data.frame(tab, stringsAsFactors=F)
names(tab) <- "ID"

for (i in 1:length(files)){
data <- read.csv(files[i], sep="\t", header=F, stringsAsFactors=F)

names(data) <- c("query", "ref", "ident", "length", "mism", "gap", "qstart", "qend", "target_s", "target_e", "e.value", "bitscore")

data <- data[,c(-11,-12)]

if(JV){ # legacy version, OTU sorting by OTU number, not recommended any where else.

data <- data.frame(data, "OTUnum"=as.numeric(sub("OTU", "", data$ref)), stringsAsFactors=F)
temp <- data[order(data$OTUnum, decreasing=F),]
temp2 <- temp[order(temp$ident, decreasing=T),]
temp3 <- temp2[order(temp2$query),]

#temp4[temp4$query=="M00517:488:000000000-BV32D:1:1101:19803:2074;size=9961;",]

temp4 <- temp3[!duplicated(temp3$query),]

data <- temp4

}


data <- cbind(data, "abund"=as.numeric(sub(".*size=(.*);", "\\1", data$query)), stringsAsFactors=F)

#head(data)

temp <- aggregate(data$abund, by=list(data$ref), FUN="sum")
tab <- merge(tab , temp, by.x="ID", by.y="Group.1", all=T, sort=T)
names(tab)[i+1] <- sub(".*2_mapping/(.*).txt", "\\1", files[i])
}
names(tab) <- sub("_PE$", "", names(tab))

tab <- tab[-1,] # remove NULL entry in the beginning
tab[is.na(tab)] <- 0

# re add empty files

if(length(excluded2)>0){
# add excluded sequences
excluded2 <- sub(".*/(.*)", "\\1", excluded2)
excluded2 <- sub("(.*).txt", "\\1", excluded2)
excluded2 <- sub("(.*)_PE_.*", "\\1", excluded2)

if(length(excluded2)>0){
zero <- NULL
for(j in 1:length(excluded2)){
zero <- cbind(zero, rep(0, nrow(tab)))
}
}

# add zero to table
zero <- data.frame(zero, stringsAsFactors=F)
names(zero) <- excluded2

# order table by file names!
tab2 <- data.frame(tab, zero, stringsAsFactors=F)
tab <- tab2[,c(1, order(names(tab2)[-1])+1)]
}
# end readd




sequ <- read.fasta(refDB, forceDNAtolower=F, as.string=T)

# KEEP ONLY HITS OR KEEP ALL IN DB
if(onlykeephits){
temp2 <- match(tab$ID, attr(sequ, "name"))
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
rel_abund[is.na(rel_abund)] <- 0 # empty cols
# write rel abundance tab
rel_abund <- rel_abund[order(rowSums(rel_abund[-c(1, ncol(rel_abund))]), decreasing=T),] # sort table by row sums
write.csv(file=paste(folder, "/3_rel_abundnace_ZEROs.csv", sep=""), rel_abund, row.names=F)



# write RAW table
tab2 <- tab2[order(rowSums(tab2[-c(1, ncol(tab2))]), decreasing=T),] # sort table by row sums

write.csv(file=paste(folder, "/3_Raw_hit_table.csv", sep=""), tab2, row.names=F)


# write stats file
temp <- colSums(tab2[-c(1, ncol(tab2))])
expTab <- data.frame(expTab, "reads_mapped"=temp, stringsAsFactors=F)
temp <- round(expTab$reads_mapped/expTab$sequ_in*100, 2)
temp[is.na(temp)] <- 0
expTab <- data.frame(expTab, "pct_pass"=temp, stringsAsFactors=F)

row.names(expTab) <- 1:nrow(expTab)

write.csv(expTab, file=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_3_pct_matched.csv", sep=""))


# make plots!
# % matched

Sequences_lost(expTab$sequ_in, expTab$reads_mapped, expTab$ID, out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_Reads_mapped.pdf", sep=""), main=paste(folder, ": Reads mapped (with ", id, ")", sep=""))
Sequences_lost(expTab$sequ_in, expTab$reads_mapped, expTab$ID, out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_Reads_mapped_rel.pdf", sep=""), main=paste(folder, ": Reads mapped (with ", id, ")", sep=""), rel=T)




# heatmap
if(heatmap){
pdf(paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_rel_zero2.pdf", sep=""), height=(nrow(rel_abund)+20)/10, width=(ncol(rel_abund)-1)/2)

temp_heat <- rel_abund[,2:(ncol(rel_abund)-1)]
row.names(temp_heat) <- rel_abund[,1]

OTU_heatmap(temp_heat, abundance=F, col=rev(c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")))
dev.off()
}




cat(file=paste(folder, "/robots.txt", sep=""), "\n# DELETE_START", files_to_delete, "# DELETE_END", append=T, sep="\n")

temp <- "\nModule completed!"
message(temp)
cat(file="log.txt", paste(Sys.time(), "*** Module completed!", "", sep="\n"), append=T, sep="\n")


}

