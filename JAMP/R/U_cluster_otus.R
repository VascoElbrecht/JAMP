# U_cluster_otus v0.1

U_cluster_otus <- function(files="latest", minuniquesize=2, strand="plus", filter=0.01, filterN=1){
#, unoise_min=NA - unoise denoising removed, no longer supported!

folder <- Core(module="U_cluster_otus")
cat(file="log.txt", c("Version v0.2", "\n"), append=T, sep="\n")
message(" ")

if (files[1]=="latest"){
source(paste(folder, "/robots.txt", sep=""))
files <- list.files(paste(last_data, "/_data", sep=""), full.names=T)
}

# Dereplicate files using USEARCH
dir.create(paste(folder, "/_data/1_derep_inc_singletons", sep=""))

new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub("_PE.*", "_PE_derep.fasta", new_names)
new_names <- sub("_data", "_data/1_derep_inc_singletons", new_names)
new_names <- paste(folder, "/", new_names, sep="")

cmd <- paste("-fastx_uniques \"", files, "\" -fastaout \"", new_names, "\" -sizeout",  sep="")

temp <- paste(length(files), " files are dereplicated (incl. singletons!):", sep="")
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

# unoise filtering of individual samples
if (F){
dir.create(paste(folder, "/_data/1_derep_unoise3", sep=""))

# change names
denoised_names <- new_names
denoised_names <- sub("1_derep_inc_singletons", "1_derep_unoise3", denoised_names)
denoised_names <- sub(".fasta", "_unoise3.fasta", denoised_names)

cmd <- paste("-unoise3 ", new_names, " -zotus ", denoised_names, " -minsize ", unoise_min, sep="")

temp <- paste("\nDenoising ", length(cmd), " files using unoise3 wiht a minimum cluster size of ", unoise_min, ":", sep="")
message(temp)
cat(file="log.txt", "\n", temp, append=T, sep="")


for(i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)

A <- c(paste("usearch", cmd[i], sep=" "), A, "\n\n\n")
cat(A, file=paste(folder, "/_stats/1_unoise3_logs.txt", sep=""), sep="\n", append=T)

temp <- paste(sub("_data/1_derep_unoise3/(.*)_PE_derep_unoise3.fasta", "\\1", denoised_names[i]), ": ", sub(".*100.0% (.*) good.*", "\\1 Amplicons keept", A[grep("100.0% .* good", A)]), sep="")
message(temp)
cat(file="log.txt", "\n", temp, append=T, sep="")

}

new_names <- denoised_names # use denoised data for OTU clustering

} # end unoise


# 2 make OTUs!
# merge all files into one

dir.create(paste(folder, "/_data/2_OTU_clustering", sep=""))

cmd <- paste(paste(new_names, collapse=" "), " > ", folder, "/_data/2_OTU_clustering/A_all_files_united.fasta", collapse="", sep="")
A <- system2("cat", cmd, stdout=T, stderr=T)

#check <- readLines("_data/2_OTU_clustering/A_all_files_united.fasta")
#count <- as.numeric(sub(".*size=(.*);", "\\1", check))
#sum(count, na.rm=T)

# write logs

temp <- paste(length(files), " dereplicated files where merged (inc singleotns) into file:\n\"",folder, "/_data/2_OTU_clustering/A_all_files_united.fasta\"", sep="")
message("\n", temp)
cat(file="log.txt", "\n", temp, append=T, sep="\n")


cat(file=paste(folder, "/_stats/2_OTU_clustering_log.txt", sep=""), temp, "", paste("cat", cmd), append=T, sep="\n")

# dereplicate "A_all_files_united.fasta" using Vsearch!
cmd <- paste("-derep_fulllength ", folder, "/_data/2_OTU_clustering/A_all_files_united.fasta -output ", folder, "/_data/2_OTU_clustering/B_all_derep_min", minuniquesize, ".fasta -sizein -sizeout -minuniquesize ", minuniquesize, sep="")

filename_all_unique <- paste("B_all_derep_min", minuniquesize, ".fasta", sep="")

A <- system2("vsearch", cmd, stdout=T, stderr=T)

temp <- paste("Total number of sequences (not dereplicated): ", sub(".*nt in (.*) seqs.*", "\\1", A[grep("seqs, min", A)]), "\n", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

temp <- paste("United sequences are dereplicated with minuniquesize = ", minuniquesize , " into a total of ", sub("(.*) unique sequences.*", "\\1", A[grep(" unique sequences", A)]), " unique sequences.", "\n", "File prepared for OTU clustering: \"", filename_all_unique, "\"", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

# derep log
cat(file=paste(folder, "/_stats/2_OTU_clustering_log.txt", sep=""), "\n", A, "", paste("cat", cmd), append=T, sep="\n")

# Actual clustering of dereplicated file 
OTU_file <- sub(".fasta", "_OTUs.fasta", filename_all_unique)
OTU_file <- sub("B_", "C_", filename_all_unique)

cmd <- paste(" -cluster_otus ", folder, "/_data/2_OTU_clustering/", filename_all_unique, " -otus ", folder, "/_data/2_OTU_clustering/", OTU_file, " -uparseout ", folder, "/_data/2_OTU_clustering/", sub(".fasta", "_OTUtab.txt", OTU_file), " -relabel OTU_ -strand ", strand, sep="")

A <- system2("usearch", cmd, stdout=T, stderr=T) # cluster OTUs!

chimeras <- sub(".*OTUs, (.*) chimeras\r", "\\1", A[grep("chimeras\r", A)])
OTUs <- sub(".*100.0% (.*) OTUs, .* chimeras\r", "\\1", A[grep("chimeras\r", A)])

temp <- paste("\n", "Clustering reads from \"", filename_all_unique, "\nminuniquesize = ", minuniquesize, "\nstrand = ", strand, "\nChimeras discarded: ", chimeras, "\nOTUs written: ", OTUs, " -> file \"", OTU_file, "\"\n", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")


# compare against refernce sequences (including singletons)

dir.create(paste(folder, "/_data/3_Compare_OTU_derep/", sep=""))
dir.create(paste(folder, "/_stats/3_Compare_OTU_derep/", sep=""))

blast_names <- sub("1_derep_inc_singletons", "3_Compare_OTU_derep", new_names)
blast_names <- sub("1_derep_unoise3", "3_Compare_OTU_derep", blast_names)
blast_names <- sub("_PE_derep.*.fasta", ".txt", blast_names)
log_names <- sub("_data", "_stats", blast_names)


cmd <- paste("-usearch_global ", new_names, " -db ", "\"", folder, "/_data/2_OTU_clustering/", OTU_file, "\"", " -strand plus -id 0.97 -blast6out \"", blast_names, "\" -maxhits 1", sep="")


temp <- paste("Comparing ", length(cmd)," files with dereplicated reads (incl. singletons) against OTUs \"", folder, "/", OTU_file, "\" using \"usearch_global\".\n", sep="")
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
names(tab)[i+1] <- sub(".*derep/(.*).txt", "\\1", files[i])
}

head(tab)

tab <- tab[-1,] # remove NULL entry in the beginning
tab[is.na(tab)] <- 0

mrew <- tab$ID
mrew <- as.numeric(gsub("OTU_(.*)", "\\1", mrew))
tab <- tab[order(as.numeric(mrew)),]

# add sequences!

sequ <- read.fasta(paste(folder, "/_data/2_OTU_clustering/", OTU_file, sep=""), forceDNAtolower=F, as.string=T)

tab2 <- cbind("sort"=sub("OTU_", "", tab[,1]), tab, "sequ"=unlist(sequ))


write.csv(file=paste(folder, "/3_Raw_OTU_table.csv", sep=""), tab2, row.names=F)

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



temp <- "\nModule completed!"
message(temp)
cat(file="log.txt", paste(Sys.time(), "*** Module completed!", "", sep="\n"), append=T, sep="\n")


}

