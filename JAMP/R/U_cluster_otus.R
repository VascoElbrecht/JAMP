# U_cluster_otus v0.1

U_cluster_otus <- function(files="latest", minuniquesize=2, otu_radius_pct=3, strand="plus", filter=0.01, filterN=1){

Core(module="U_cluster_otus")
cat(file="../log.txt", c("Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files[1]=="latest"){
source("robots.txt")
files <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T)
}

# Dereplicate files using USEARCH
dir.create("_data/1_derep_inc_singletons")

new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub("_PE.*", "_PE_derep.fasta", new_names)
new_names <- sub("_data", "_data/1_derep_inc_singletons", new_names)

cmd <- paste("-fastx_uniques \"", files, "\" -fastaout \"", new_names, "\" -sizeout",  sep="")

temp <- paste(length(files), " files are dereplicated (incl. singletons!):", sep="")
cat(file="../log.txt", temp , append=T, sep="\n")
message(temp)


temp <- new_names
for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
meep <- sub(".*_data/(.*)", "\\1", temp[i])
cat(file="../log.txt", meep, append=T, sep="\n")
cat(file="_stats/1_derep_logs.txt", meep, A, "\n", append=T, sep="\n")
message(meep)
}

# 2 make OTUs!
# merge all files into one

dir.create("_data/2_OTU_clustering")

cmd <- paste(paste(new_names, collapse=" "), "> _data/2_OTU_clustering/A_all_files_united.fasta", collapse=" ")
A <- system2("cat", cmd, stdout=T, stderr=T)

#check <- readLines("_data/2_OTU_clustering/A_all_files_united.fasta")
#count <- as.numeric(sub(".*size=(.*);", "\\1", check))
#sum(count, na.rm=T)

# write logs

temp <- paste(length(files), " dereplicated files where merged (inc singleotns) into file:\n\"_data/2_OTU_clustering/A_all_files_united.fasta\"", sep="")
message("\n", temp)
cat(file="../log.txt", "\n", temp, append=T, sep="\n")


cat(file="_stats/2_OTU_clustering_log.txt", temp, "", paste("cat", cmd), append=T, sep="\n")

# dereplicate "A_all_files_united.fasta" using Vsearch!
cmd <- paste("-derep_fulllength _data/2_OTU_clustering/A_all_files_united.fasta -output _data/2_OTU_clustering/B_all_derep_min", minuniquesize, ".fasta -sizein -sizeout -minuniquesize ", minuniquesize, sep="")

filename_all_unique <- paste("B_all_derep_min", minuniquesize, ".fasta", sep="")

A <- system2("vsearch", cmd, stdout=T, stderr=T)

temp <- paste("Total number of sequences (not dereplicated): ", sub(".*nt in (.*) seqs.*", "\\1", A[grep("seqs, min", A)]), "\n", sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")

temp <- paste("United sequences are dereplicated with minuniquesize = ", minuniquesize , " into a total of ", sub("(.*)uniques written.*", "\\1", A[grep("uniques written", A)]), "unique sequences.", "\n", "File prepared for OTU clustering: \"", filename_all_unique, "\"", sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")

# derep log
cat(file="_stats/2_OTU_clustering_log.txt", "\n", A, "", paste("cat", cmd), append=T, sep="\n")

# Actual clustering of dereplicated file 
OTU_file <- sub(".fasta", "_OTUs.fasta", filename_all_unique)
OTU_file <- sub("B_", "C_", filename_all_unique)

cmd <- paste(" -cluster_otus _data/2_OTU_clustering/", filename_all_unique, " -otus _data/2_OTU_clustering/", OTU_file, " -uparseout _data/2_OTU_clustering/", sub(".fasta", "_OTUtab.txt", OTU_file), " -relabel OTU_ -otu_radius_pct ", otu_radius_pct, " -strand ", strand, sep="")

A <- system2("usearch", cmd, stdout=T, stderr=T) # cluster OTUs!

chimeras <- sub(".*OTUs, (.*) chimeras\r", "\\1", A[grep("chimeras\r", A)])
OTUs <- sub(".*100.0% (.*) OTUs, .* chimeras\r", "\\1", A[grep("chimeras\r", A)])

temp <- paste("\n", "Clustering reads from \"", filename_all_unique, "\nminuniquesize = ", minuniquesize, "\notu_radius_pct = ", otu_radius_pct, "\nstrand = ", strand, "\nChimeras discarded: ", chimeras, "\nOTUs written: ", OTUs, " -> file \"", OTU_file, "\"\n", sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")


# compare against refernce sequences (including singletons)

dir.create("_data/3_Compare_OTU_derep/")
dir.create("_stats/3_Compare_OTU_derep/")

blast_names <- sub("1_derep_inc_singletons", "3_Compare_OTU_derep", new_names)
blast_names <- sub("_PE_derep.fasta", ".txt", blast_names)
log_names <- sub("_data", "_stats", blast_names)


cmd <- paste("-usearch_global ", new_names, " -db ", "\"_data/2_OTU_clustering/", OTU_file, "\"", " -strand plus -id 0.97 -blast6out \"", blast_names, "\" -maxhits 1", sep="")


temp <- paste("Comparing ", length(cmd)," files with dereplicated reads (incl. singletons) against OTUs \"", OTU_file, "\" using \"usearch_global\".\n", sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")

exp <- NULL
temp <- new_names
for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
cat(file= log_names[i], paste("usearch ", cmd[i], sep=""), "\n", A, append=F, sep="\n")

meep <- sub(".*singletons/(.*)", "\\1", temp[i])
pass <- sub(".*, (.*)% matched\r", "\\1", A[grep("matched\r", A)])
exp <- rbind(exp, c(meep, pass))
glumanda <- paste(meep," - ", pass, "% reads matched", sep="")
cat(file="../log.txt", glumanda, append=T, sep="\n")
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

sequ <- read.fasta(paste("_data/2_OTU_clustering/", OTU_file, sep=""), forceDNAtolower=F, as.string=T)

tab2 <- cbind("sort"=sub("OTU_", "", tab[,1]), tab, "sequ"=unlist(sequ))


write.csv(file="3_Raw_OTU_table.csv", tab2, row.names=F)

temp <- "\n\nOTU table generated (including OTU sequences): 3_Raw_OTU_table.csv"
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")

exp2 <- data.frame("ID"=exp[,1], "Abundance"=colSums(tab[-1]), "pct_pass"=exp[,2], row.names=1:length(exp[,1]))

write.csv(exp2, file="_stats/3_pct_matched.csv")

#### end raw data table

### make abundance filtering
if(!is.na(filter)){
#tab2 <- read.csv(file="3_Raw_OTU_table.csv", stringsAsFactors=F)

start <- which(names(tab2)=="ID")+1
stop <- which(names(tab2)=="sequ")-1

temp <- tab2[, start:stop]

meep <- paste("Discarding OTUs with below ", filter, "% abundance across at least ", filterN, " out of ", ncol(temp), " samples.", sep="")
message(meep)
cat(file="../log.txt", meep, append=T, sep="\n")

temp2 <- temp
sampleabundance <- colSums(temp)
for (i in 1:ncol(temp)){
temp2[i] <- temp[i]/sampleabundance[i]*100
}

# subset OTUs
subset <- rowSums(temp2>=filter)
subset2 <- subset >= filterN

# reporting
meep <- paste("Discarded OTUs: ", sum(!subset2)," / ",  length(subset2), " (", round(100-sum(subset2)/length(subset2)*100, 2), "% discarded)", sep="")
message(meep)
cat(file="../log.txt", meep, append=T, sep="\n")

#write subsetted OTU table

exp <- tab2[subset2,]
exp <- rbind(exp, NA)

exp[nrow(exp), start:stop] <- colSums(tab2[!subset2, start:stop])
exp$ID[nrow(exp)] <- paste("below_", filter, sep="")
exp$sort[nrow(exp)] <- exp$sort[nrow(exp)-1]+1


# make folder 
dir.create("_data/4_subset/")


write.csv(exp, file=paste("_data/4_subset/4_OTU_sub_", filter, "_not_rematched.csv", sep=""), row.names=F)

write.fasta(as.list(exp$sequ[-nrow(exp)]), exp$ID[-nrow(exp)], file.out=paste("_data/4_subset/4_OTU_sub_", filter, ".fasta", sep=""))

} # end subsetting



temp <- "\nModule completed!"
message(temp)
cat(file="../log.txt", paste(Sys.time(), temp, "", sep="\n"), append=T, sep="\n")

setwd("../")
}

