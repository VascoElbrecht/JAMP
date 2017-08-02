# Haplotyping v0.1

Denoise <- function(files="latest",  strategy="unoise", unoise_alpha=5, minsize=10, minrelsize=0.001){



Core(module="Denoising")
cat(file="../log.txt", c("Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files[1]=="latest"){
source("robots.txt")
files <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T)
}



# count sequences in each file
counts <- Count_sequences(files, fastq=F)
size <- round(counts* minrelsize/100) # get nim abundance
size[size<minsize] <- minsize # min size!



# Dereplicate files using Usearch
dir.create("_data/1_derep")


new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub("_PE.*", "_PE_derep", new_names)
new_names <- paste(new_names, "_size_", size, ".txt", sep="")
new_names <- sub("_data", "_data/1_derep", new_names)

cmd <- paste("-fastx_uniques \"", files, "\" -fastaout \"", new_names, "\" -sizeout", " -minuniquesize ", size,  sep="")

temp <- paste(length(files), " files are dereplicated and sequences in each sample below ", minrelsize, "% (or minuniqesize of ", minsize,")  are beeing discarded:", sep="")
cat(file="../log.txt", temp , append=T, sep="\n")
message(temp)


temp <- new_names
for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
meep <- sub(".*_data/(.*)", "\\1", temp[i])
cat(file="_stats/1_derep_logs.txt", paste("usearch ", cmd[i], sep="") , append=T, sep="\n")
cat(file="_stats/1_derep_logs.txt", meep, A, "\n", append=T, sep="\n")


log_count <- Count_sequences(new_names[i], count_size=T)
log <- paste(sub(".*_data/1_derep/(.*)", "\\1", temp[i]), ": ", log_count, " of ", counts[i], " keept (", round((log_count/counts[i])*100, digits=4), "%, min size: ", size[i],")", sep="")
cat(file="../log.txt", log , append=T, sep="\n")
message(log)
}


# merge all files into one!


cat(file="_stats/1_derep_logs.txt", paste("\nCombining all files in a single file (samples_pooled.txt):\n", paste("cmd", cmd, collapse="", sep=""), collapse="", sep="") , append=T, sep="\n")
cat(file="../log.txt", "\nCombining all files in a single file (samples_pooled.txt)\n", append=T, sep="\n")

# dereplicating pooled file
message("Combining all files in a single file (samples_pooled.txt)")
cmd <- paste(paste(new_names, collapse=" "), "> _data/1_derep/samples_pooled.txt")
system2("cat", cmd)

# dereplicating files
info <- "Dereplicating pooled sequences! (no min size)"
message(info)
cat(file="../log.txt", info, append=T, sep="\n")

cmd <- "-fastx_uniques \"_data/1_derep/samples_pooled.txt\" -fastaout \"_data/1_derep/samples_pooled_derep.txt\" -sizein -sizeout"
A <- system2("usearch", cmd, stdout=T, stderr=T)

cat(file="_stats/1_derep_logs.txt", paste("usearch", cmd, sep=""), append=T, sep="\n")
cat(file="_stats/1_derep_logs.txt", A, append=T, sep="\n")

# renaming all sequences!
info <- "Renaming pooled sequences, and applying same names to the dereplicated files."
message(info)
cat(file="../log.txt", info, append=T, sep="\n")


haplo <- read.fasta("_data/1_derep/samples_pooled_derep.txt", forceDNAtolower=F, as.string=T)

temp <- sub(".*(;size=.*;)", "\\1", names(haplo))
temp2 <- paste("haplo_", 1:length(haplo), temp, sep="")

write.fasta(haplo, temp2, "_data/1_derep/samples_pooled_derep_renamed.txt")


# rename single files!
dir.create("_data/2_renamed")
renamed <- sub(".txt", "_renamed.txt", new_names)
renamed <- sub("/1_derep/", "/2_renamed/", renamed)


for (i in 1:length(new_names)){
sample <- read.fasta(new_names[i], as.string=T, forceDNAtolower=F)
matched <- match(sample, haplo)

new_sample <- haplo[matched] # DNA sequences
new_haplo_seque_names <- paste("haplo_", matched, sub(".*(;size=.*;)", "\\1", names(sample)), sep="") # sizes

write.fasta(new_sample, new_haplo_seque_names, renamed[i])
}

# UNOISE3
# Apply denoising on the POOLED dereplicated renamed file!

info <- paste("Denoising the file 1_derep/samples_pooled_derep_renamed.txt (containing", length(renamed), "samples).")
message(info)
cat(file="../log.txt", info, append=T, sep="\n")

cmd <- paste("-unoise3 \"_data/1_derep/samples_pooled_derep_renamed.txt\" -zotus \"_data/1_derep/samples_pooled_+_denoised.txt\" -unoise_alpha ", unoise_alpha,  sep="")

A <- system2("usearch", cmd, stdout=T, stderr=T)
cat(file="_stats/2_unoise.txt", c(info, "", paste("usearch", cmd), "", A), append=T, sep="\n")

info <- paste("Denoising compelte! ", Count_sequences("_data/1_derep/samples_pooled_derep_renamed.txt", fastq=F), " sequences were denoised using ", strategy, ".", "\nA total of ", sub(".*100.0% (.*) good, .* chimeras\r", "\\1", A[length(A)-1]), " haplotypes remained after denoising!", sep="")
message(info)
cat(file="../log.txt", info, append=T, sep="\n")



#Zotus, get old names back (original haplotypes)!

Zotus <- read.fasta("_data/1_derep/samples_pooled_+_denoised.txt", as.string=T, forceDNAtolower=F)
renamed_sequ <- read.fasta("_data/1_derep/samples_pooled_derep_renamed.txt", as.string=T, forceDNAtolower=F)

matched <- match(Zotus, renamed_sequ)
new_sample <- renamed_sequ[matched] # DNA sequences

write.fasta(new_sample, names(new_sample), "_data/1_derep/samples_pooled_+_denoised_renamed.txt")


names(new_sample) <- sub(";size(.*);", "", names(new_sample))
haplotypes <- new_sample

dir.create("_data/3_unoise")

denoised_sequences <- sub("2_renamed", "3_unoise", renamed)
denoised_sequences <- sub("PE_derep_size_", "", denoised_sequences)
denoised_sequences <- sub("_renamed.txt", "_denoised.txt", denoised_sequences)

i<-1
# chack dereplicated files agains the list of haplotypes (unoise3 all files)
for (i in 1:length(cmd)){

sample <- read.fasta(renamed[i], as.string=T, forceDNAtolower=F)
sample_keep <- sample[sample%in%haplotypes]

write.fasta(sample_keep, names(sample_keep), denoised_sequences[i])


info <- paste(sub("_data/3_unoise/(.*)_denoised.txt", "\\1", denoised_sequences[i]), ": ", length(sample_keep), " of ", length(sample), " sequences remained after denoising (", round(length(sample_keep)/length(sample)*100, 2), "%)", sep="")
message(info)
cat(file="../log.txt", info, append=T, sep="\n")

}


# bringing size back!
















# 2 make OTUs!
# merge all files into one

dir.create("_data/3_OTU_clustering")

cmd <- paste(paste(new_names2, collapse=" "), "> _data/3_OTU_clustering/A_all_files_united.fasta", collapse=" ")
A <- system2("cat", cmd, stdout=T, stderr=T)

temp <- paste(length(files), " dereplicated files where merged into file:\n\"_data/3_OTU_clustering/A_all_files_united.fasta\"", sep="")
message("\n", temp)
cat(file="../log.txt", "\n", temp, append=T, sep="\n")


cat(file="_stats/3_OTU_clustering_log.txt", temp, "", paste("cat", cmd), append=T, sep="\n")

# dereplicate "A_all_files_united.fasta" using Vsearch!
cmd <- paste("-derep_fulllength _data/3_OTU_clustering/A_all_files_united.fasta -output _data/3_OTU_clustering/B_all_derep.fasta -sizein -sizeout -relabel Uniq", sep="")

A <- system2("vsearch", cmd, stdout=T, stderr=T)

temp <- paste("Total number of sequences (not dereplicated): ", sub(".*nt in (.*) seqs.*", "\\1", A[grep("seqs, min", A)]), "\n", sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")

temp <- paste("United sequences are dereplicated + size filtered into a total of ", sub("(.*) unique sequences.*", "\\1", A[grep(" unique sequences", A)]), " unique sequences.", "\n", "File prepared for OTU clustering: B_all_derep.fasta", sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")

# derep log
cat(file="_stats/3_OTU_clustering_log.txt", "\n", A, "", paste("cat", cmd), append=T, sep="\n")


# Actual OTU clustering of dereplicated filtered file! 

cmd <- paste(" -cluster_otus _data/3_OTU_clustering/B_all_derep.fasta -otus _data/3_OTU_clustering/C_OTUs.fasta -uparseout _data/3_OTU_clustering/C_OTU_table.txt -relabel OTU_ -otu_radius_pct ", otu_radius_pct, " -strand ", strand, sep="")

A <- system2("usearch", cmd, stdout=T, stderr=T) # cluster OTUs!

# cluster log
cat(file="_stats/3_OTU_clustering_log.txt", "\n", paste("usearch", cmd), "", A, "", append=T, sep="\n")

chimeras <- sub(".*OTUs, (.*) chimeras\r", "\\1", A[grep("chimeras\r", A)])
OTUs <- sub(".*100.0% (.*) OTUs, .* chimeras\r", "\\1", A[grep("chimeras\r", A)])

temp <- paste("\n", "Clustering reads from\n\"B_all_derep.fasta\" \notu_radius_pct = ", otu_radius_pct, "\nstrand = ", strand, "\nChimeras discarded: ", chimeras, "\nOTUs written: ", OTUs, " -> file \"C_OTUs.fasta\"\n", sep="")
message(temp)
cat(file="../log.txt", temp, append=T, sep="\n")

# RNENAME
# compare reads against dereplicated and RENAME sequences!

dir.create("_data/4_rename")

DNA_master <- read.fasta("_data/3_OTU_clustering/B_all_derep.fasta", as.string=T, forceDNAtolower=F)
names(DNA_master) <- sub("(.*);size=.*", "\\1", names(DNA_master))

new_names3 <- sub("2_MinMax", "4_rename", new_names2)

for (i in 1:length(new_names2)){

sample <- read.fasta(new_names2[i], as.string=T, forceDNAtolower=F)

size <- sub(".*(;size=.*;)", "\\1", names(sample))
sequnames <- paste(names(DNA_master)[match(sample, DNA_master)], size, sep="")

write.fasta(sample, names= sequnames, new_names3[i])
}
message("read renamed! ame as in \"B_all_derep.fasta\"")
# END renaming


# Mapp reads (filtered abund & MM) against OTUs
blast_names <- sub("4_rename", "5_mapp", new_names3)
log_names <- sub("_data", "_stats", blast_names)
dir.create("_data/5_mapp")

cmd <- paste("-usearch_global ", new_names3, " -db ", "\"_data/3_OTU_clustering/C_OTUs.fasta\"", " -strand plus -id ", (100-otu_radius_pct)/100, " -blast6out \"", blast_names, "\" -maxhits 1", sep="")


for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
cat(file= "_stats/5_mapping.txt", paste("vsearch", cmd[i], sep=""), A, "\n\n\n", append=T, sep="\n")
}
message("Reads remapped!")


# HAPLO TYPING STUFF
# subset haplotypes + writing individual files

dir.create("_data/6_haplotypes") # make folders!
folder <- paste("_data/6_haplotypes/", sub("_data/5_mapp/(.*)_PE_derep.*", "\\1", blast_names), sep="")


i <- 1
for (i in 1:length(blast_names)){

data <- read.csv(blast_names[i], header=F, sep="\t", stringsAsFactors=F)

data$V11 <- sub("(.*);size=.*", "\\1", data$V1)
data$V12 <- as.numeric(sub(".*;size=(.*);", "\\1", data$V1))

# subset low abundant OTUs
temp <- aggregate(data$V12, list(data$V2), "sum")
OTUsum <- sum(temp$x)

temp <- cbind(temp, "relabund"=temp$x/OTUsum*100)

# subset OTUs to keep!
temp <- temp[temp$relabund>=minOTUabund,]

data <- data[data$V2%in%temp$Group.1,]

report <- paste("Subsetting OTUs with ", minOTUabund, " % anundance; Keeping ", nrow(temp), " OTUs", sep="")
message(report)

temp_foldername <- folder[i]
dir.create(temp_foldername)

#Processing of individual OTUs
# temp = list of OTUs with high enough abundnace (minOTUabund)
k <- 15
for(k in 1:nrow(temp)){

dir.create(paste(temp_foldername, "/", temp$Group.1[k], sep=""), showWarnings=F)

meep <- data[data$V2==temp$Group.1[k],]

# convert to rel abundance
meep$V12 <- meep$V12/sum(meep$V12)*100

meep <- cbind(meep, "keep"=meep$V12 >= AbundInside)

#recalculate rel abundance of OTUs left
meep <- cbind(meep, "keeprel"=meep$V12)
meep$keeprel[meep$keep] <- round(meep$keeprel[meep$keep]/sum(meep$keeprel[meep$keep])*100, 2)
meep$keeprel[!meep$keep] <- NA

# write table as csv!
write.csv(meep, paste(temp_foldername, "/", temp$Group.1[k], "/", temp$Group.1[k], "_tab.csv", sep=""), row.names=F)


# save fasta files
fasta_save <- DNA_master[names(DNA_master)%in%meep$V11[meep$keep]]
glumanda_names <- meep$V1[meep$keep][match(sub("(.*);size.*", "\\1", meep$V1[meep$keep]), names(fasta_save))] # same sorting

write.fasta(fasta_save, names=glumanda_names, paste(temp_foldername, "/", temp$Group.1[k], "/", temp$Group.1[k], "_sequ.txt", sep=""))



# make plot!
pdf(file= paste(temp_foldername, "/", temp$Group.1[k], "/", temp$Group.1[k], "_plot.pdf", sep=""), width=6, height=6, useDingbats=F)
plot(meep$V12, ylim=c(0.01, 100), log="y", yaxt="n", ylab="rel. proportions within OTU", main=temp$Group.1[k])
axis(side=2, at=c(100, 10, 1, 0.1, 0.01, 0.001), labels=c("100", "10", "1", "0.1","0.01", "0.001"), las=2)
axis(side=2, at=c(seq(20,90,10), seq(2,9,1), seq(0.2,0.9,0.1), seq(0.02,0.09,0.01)), labels=F, las=2, tck=-0.01)
lines(c(0,nrow(meep)), c(AbundInside, AbundInside), col="Red")
dev.off()

# make plot, write fasta of sequ to keep

}
}

i<-1
# unite haplotypes into single table

master_tab <- data.frame("V11"= "Uniq86")

for (i in 1:length(folder)){

OTUs <- list.files(folder[i])

exp <- NULL # extract individual haplotypes
for (k in 1:length(OTUs)){

otu <- read.csv(paste(folder[i], "/", OTUs[k], "/", OTUs[k], "_tab.csv", sep=""), stringsAsFactors=F)

# count sequences in OTU
otu_abund <- sum(as.numeric(sub(".*size=(.*);", "\\1", otu$V1)))

otu <- otu[otu$keep,] # keep high abund haplotypes
otu <- cbind(otu[c(2, 11, 14)], "indiv"=as.numeric(sub(".*size=(.*);", "\\1", otu$V1)), otu_abund) # extract useful info

exp <- rbind(exp, otu)
}

temp_name <- sub("_data/6_haplotypes/", "", folder[i])
names(exp) <- c(paste(temp_name, "_OTU", sep=""), "V11", paste(temp_name, "_rel", sep=""), paste(temp_name, "_abund", sep=""), paste(temp_name, "_OTU_abund", sep=""))

# merge exp haplotype tables from samples

master_tab <- merge(master_tab, exp, by="V11", all=T)

}

# reorganise haplo table 

master_tab <- master_tab[c(1, (1:length(folder))*4-2, (1:length(folder))*4-1, (1:length(folder))*4,(1:length(folder))*4+1)]

master_tab$V11 <- as.character(master_tab$V11)

OTU <- NULL

for (i in 1:length(folder)){
keep <- as.vector(!is.na(master_tab[(1+i)]))
OTU[keep] <- master_tab[keep, (1+i)]
}

master_tab <- cbind(OTU, master_tab[-c(2:(length(folder)+1))])


master_tab <- master_tab[order(as.numeric(sub("OTU_", "", master_tab$OTU)), as.numeric(sub("Uniq", "", master_tab$V11)), decreasing=F),]
names(master_tab)[2] <- "Haplotypes"

# add sequences to table
master_tab <- cbind(sort= c(1:nrow(master_tab)), master_tab, "sequ"=unlist(DNA_master[match(master_tab$Haplotypes, names(DNA_master))]))

write.csv(master_tab, file="haplo_tab.csv", row.names=F)


write.fasta(DNA_master[match(master_tab$Haplotypes, names(DNA_master))], names=paste(master_tab$OTU, master_tab$Haplotypes, sep="_"), "haplo_fasta.txt")




temp <- "\nModule completed!"
message(temp)
cat(file="../log.txt", paste(Sys.time(), temp, "", sep="\n"), append=T, sep="\n")

setwd("../")


}
