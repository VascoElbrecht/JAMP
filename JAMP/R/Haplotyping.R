# Haplotyping v0.1

haplotyping <- function(folder=NA, ampliconLength=NA, minsize=5, minrelsize=0.001, minOTUabund=0.1, AbundInside=1, otu_radius_pct=3, strand="plus"){



Core(module="Haplotyping")
cat(file="../log.txt", c("Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files[1]=="latest"){
source("robots.txt")
files <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T)
}

# Dereplicate files using Usearch
dir.create("_data/1_derep")

# count sequences in each file
counts <- Count_sequences(files, fastq=F)
size <- round(counts* minrelsize/100) # get nim abundance
size[size<minsize] <- minrelsize # min size!


new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub("_PE.*", "_PE_derep", new_names)
new_names <- paste(new_names, "_", size, ".txt", sep="")
new_names <- sub("_data", "_data/1_derep", new_names)

cmd <- paste("-fastx_uniques \"", files, "\" -fastaout \"", new_names, "\" -sizeout", " -minuniquesize ", size,  sep="")

temp <- paste(length(files), " files are dereplicated and sequences bwelow ", minrelsize, "% (or minuniqesize of ", minsize, ") are beeing discarded:", sep="")
cat(file="../log.txt", temp , append=T, sep="\n")
message(temp)

temp <- new_names
for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
meep <- sub(".*_data/(.*)", "\\1", temp[i])
cat(file="_stats/1_derep_logs.txt", meep, A, "\n", append=T, sep="\n")

log_count <- Count_sequences(paste("_data/", meep, sep=""))
log <- paste(sub(".*_data/1_derep/(.*)", "\\1", temp[i]), ": ", log_count, " of ", counts[i], " keept (", round((log_count/counts[i])*100, digits=4), "%, min size: ", size[i],")", sep="")
cat(file="../log.txt", log , append=T, sep="\n")
message(log)
}

# 2 make OTUs!
# merge all files into one

dir.create("_data/2_OTU_clustering")

cmd <- paste(paste(new_names, collapse=" "), "> _data/2_OTU_clustering/A_all_files_united.fasta", collapse=" ")
A <- system2("cat", cmd, stdout=T, stderr=T)

temp <- paste(length(files), " dereplicated files where merged into file:\n\"_data/2_OTU_clustering/A_all_files_united.fasta\"", sep="")
message("\n", temp)
cat(file="../log.txt", "\n", temp, append=T, sep="\n")


cat(file="_stats/2_OTU_clustering_log.txt", temp, "", paste("cat", cmd), append=T, sep="\n")

# dereplicate "A_all_files_united.fasta" using Vsearch!
cmd <- paste("-derep_fulllength _data/2_OTU_clustering/A_all_files_united.fasta -output _data/2_OTU_clustering/B_all_derep.fasta -sizein -sizeout -relabel Uniq", sep="")

A <- system2("vsearch", cmd, stdout=T, stderr=T)

temp <- paste("Total number of sequences (not dereplicated): ", sub(".*nt in (.*) seqs.*", "\\1", A[grep("seqs, min", A)]), "\n", sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")

temp <- paste("United sequences are dereplicated + size filtered into a total of ", sub("(.*) unique sequences.*", "\\1", A[grep(" unique sequences", A)]), " unique sequences.", "\n", "File prepared for OTU clustering: B_all_derep.fasta", sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")

# derep log
cat(file="_stats/2_OTU_clustering_log.txt", "\n", A, "", paste("cat", cmd), append=T, sep="\n")

# min max sequence length (cutadapt)

cmd <- paste("_data/2_OTU_clustering/B_all_derep.fasta -o _data/2_OTU_clustering/C_all_derep_MM.fasta -m ", ampliconLength, " -M ", ampliconLength, sep="")
A <- system2("cutadapt", cmd, stdout=T, stderr=T)
getwd()

cat(file="_stats/2_OTU_clustering_log.txt", "\n", A, "", paste("cutadapt", cmd), append=T, sep="\n")

stats <- A
reads_in <- stats[grep("Total reads processed:", stats)[1]]
reads_in <- sub(".* processed: +", "", reads_in)
reads_in <- as.numeric(gsub(",", "", reads_in))

reads_out <- stats[grep("Reads written \\(passing filters\\):", stats)[1]]
reads_out <- sub(".* filters.: +", "", reads_out)
reads_out <- sub(" .*", "", reads_out)
reads_out <- as.numeric(gsub(",", "", reads_out))

keep <- round(reads_out/reads_in*100, digits=2)


meep <- paste("Filtering ", reads_in, " reads with MIN MAX ", ampliconLength, " bp: keep ", reads_out, " (", keep, "%)", "\nWritten in file: C_all_derep_MM.fasta", sep="")
cat(file="../log.txt", meep, append=T, sep="\n")
message(meep)



# Actual OTU clustering of dereplicated filtered file! 

cmd <- paste(" -cluster_otus _data/2_OTU_clustering/C_all_derep_MM.fasta -otus _data/2_OTU_clustering/D_OTUs.fasta -uparseout _data/2_OTU_clustering/D_OTU_table.txt -relabel OTU_ -otu_radius_pct ", otu_radius_pct, " -strand ", strand, sep="")

A <- system2("usearch", cmd, stdout=T, stderr=T) # cluster OTUs!

# cluster log
cat(file="_stats/2_OTU_clustering_log.txt", "\n", paste("usearch", cmd), "", A, "", append=T, sep="\n")

chimeras <- sub(".*OTUs, (.*) chimeras\r", "\\1", A[grep("chimeras\r", A)])
OTUs <- sub(".*100.0% (.*) OTUs, .* chimeras\r", "\\1", A[grep("chimeras\r", A)])

temp <- paste("\n", "Clustering reads from\n\"C_all_derep_MM.fasta\" \notu_radius_pct = ", otu_radius_pct, "\nstrand = ", strand, "\nChimeras discarded: ", chimeras, "\nOTUs written: ", OTUs, " -> file \"", OTU_file, "\"\n", sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")


# compare reads against dereplicated and rename sequences!

dir.create("_data/3_rename")

DNA_master <- read.fasta("_data/2_OTU_clustering/C_all_derep_MM.fasta", as.string=T, forceDNAtolower=F)



names(sample) <- names(DNA_master)[sample%in%DNA_master]


match(sample, DNA_master)

sample <- read.fasta(new_names[i], as.string=T, forceDNAtolower=F)









# Mapp reads (filtered abund & MM) against OTUs
blast_names <- sub("1_derep", "4_mapp", new_names)
log_names <- sub("_data", "_stats", blast_names)
dir.create("_data/4_mapp")

cmd <- paste("-usearch_global ", new_names, " -db ", "\"_data/2_OTU_clustering/D_OTUs.fasta\"", " -strand plus -id ", (100-otu_radius_pct)/100, " -blast6out \"", blast_names, "\" -maxhits 1", sep="")


for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
cat(file= "_stats/3_mapping.txt", paste("vsae rch", cmd[i], sep=""), A, "\n\n\n", append=T, sep="\n")
}
message("Reads remapped!")

i <- 1
# subset haplotypes
for (i in 1:length(blast_names)){

data <- read.csv(blast_names[i], header=F, sep="\t")

head(data)




}











temp <- "\nModule completed!"
message(temp)
cat(file="../log.txt", paste(Sys.time(), temp, "", sep="\n"), append=T, sep="\n")

setwd("../")


}
