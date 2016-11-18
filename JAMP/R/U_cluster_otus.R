# U_cluster_otus v0.1

U_cluster_otus <- function(files="latest", minuniquesize=2, otu_radius_pct=3, strand="plus"){

Core(module="U_cluster_otus")
cat(file="../log.txt", c("Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files=="latest"){
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

A <- system2("Vsearch", cmd, stdout=T, stderr=T)

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

A <- system2("Usearch", cmd, stdout=T, stderr=T)




message(" ")
message(" Module completed!")

cat(file="../log.txt", paste(Sys.time(), "\n", "Module completed!", "", sep="\n"), append=T, sep="\n")

setwd("../")
}

