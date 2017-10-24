# Haplotyping v0.1

Denoise <- function(files="latest",  strategy="unoise", unoise_alpha=5, minsize=10, minrelsize=0.0001, OTUmin=0.01, minhaplosize=0.003, withinOTU=5, eachsampleOTUmin=NULL, minHaploPresence=1, minOTUPresence=1, renameSamples="(.*_.*)_cut.*"){



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
new_names <- sub(".fasta", "_derep.fasta", new_names)
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
message("\nCombining all files in a single file (samples_pooled.txt)")
cmd <- paste(paste(paste("\"", new_names, "\"", sep=""), collapse=" "), "> _data/1_derep/samples_pooled.txt")
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
info <- "Renaming pooled sequences, and applying same names to the dereplicated files.\n"
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

info <- paste("\nDenoising the file 1_derep/samples_pooled_derep_renamed.txt (containing", length(renamed), "samples).")
message(info)
cat(file="../log.txt", info, append=T, sep="\n")

cmd <- paste("-unoise3 \"_data/1_derep/samples_pooled_derep_renamed.txt\" -zotus \"_data/1_derep/samples_pooled_+_denoised.txt\" -unoise_alpha ", unoise_alpha,  sep="")

A <- system2("usearch", cmd, stdout=T, stderr=T)
cat(file="_stats/2_unoise.txt", c(info, "", paste("usearch", cmd), "", A), append=T, sep="\n")

info <- paste("Denoising compelte! ", Count_sequences("_data/1_derep/samples_pooled_derep_renamed.txt", fastq=F), " sequences were denoised using ", strategy, ".", "\nA total of ", sub(".*100.0% (.*) good, .* chimeras\r", "\\1", A[length(A)-1]), " haplotypes remained after denoising!\n", sep="")
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


# check dereplicated files agains the list of haplotypes (unoise3 all files)
for (i in 1:length(denoised_sequences)){

sample <- read.fasta(renamed[i], as.string=T, forceDNAtolower=F)
sample_keep <- sample[sample%in%haplotypes]

write.fasta(sample_keep, names(sample_keep), denoised_sequences[i])


info <- paste(sub("_data/3_unoise/(.*)_denoised.txt", "\\1", denoised_sequences[i]), ": ", length(sample_keep), " of ", length(sample), " sequences remained after denoising (", round(length(sample_keep)/length(sample)*100, 2), "%)", sep="")
message(info)
cat(file="../log.txt", info, append=T, sep="\n")

}


# Cluster into OTUs (for OTU table information)

cmd <- paste(" -cluster_otus _data/1_derep/samples_pooled_+_denoised_renamed.txt -otus _data/1_derep/samples_pooled_+_denoised_renamed_OTUsequ.txt -uparseout _data/1_derep/samples_pooled_+_denoised_renamed_OTUtable.txt -relabel OTU_ -strand plus", sep="")

A <- system2("usearch", cmd, stdout=T, stderr=T) # cluster OTUs!

cat(file="_stats/2_unoise.txt", c("Clustering haplotypes into OTUs for OTU table!", "", paste("usearch", cmd), "", A), append=T, sep="\n")

chimeras <- as.numeric(sub(".*100.0% .* OTUs, (.*) chimeras\r", "\\1", A[grep("chimeras\r", A)]))
OTUs <- as.numeric(sub(".*100.0% (.*) OTUs, .* chimeras\r", "\\1", A[grep("chimeras\r", A)]))
if(is.na(chimeras)){chimeras<-0}

info <- paste("Clustered ", length(haplotypes), " haplotype sequences (cluster_otus, 3% simmilarity) into ", OTUs, " OTUs (+", chimeras, " chimeras).\nOTUs and (potentially) chimeric sequences will be included in the Haplotype table!\n", sep="" )
message(info)
cat(file="../log.txt", info, append=T, sep="\n")


# generate one united haplotype table!

OTUs <- read.csv("_data/1_derep/samples_pooled_+_denoised_renamed_OTUtable.txt", stringsAsFactors=F, sep="\t", header=F)

k <- 1
OTU_list <- NULL
for (i in 1:nrow(OTUs)){

if(OTUs$V2[i]=="OTU"){OTU_list[i] <- paste("OTU_", k, sep="")
k <- k+1} else
if(OTUs$V2[i]=="match"){OTU_list[i] <- sub(".*;top=(OTU_.*)\\(.*", "\\1", OTUs$V3[i])} else {OTU_list[i] <- OTUs$V2[i]}

}


data <- data.frame("haplotype"=names(haplotypes), "OTU"=OTU_list, stringsAsFactors=F)

for (i in 1:length(denoised_sequences)){
sample <- names(read.fasta(denoised_sequences[i]))
matched <- match(sub(";size=.*;", "", sample), data$haplotype)
abundance <- rep(0, nrow(data))
abundance[matched] <- as.numeric(sub(".*;size=(.*);", "\\1", sample))

data <- cbind(data, abundance)
names(data)[i+2] <- sub("_data/3_unoise/(.*)_denoised.txt", "\\1", denoised_sequences[i])
}

data <- cbind(data, "sequences"=unlist(haplotypes), stringsAsFactors=F)
# sort by OTUs
data <- data[order(suppressWarnings(as.numeric(sub("OTU_", "", data$OTU)))),]
data <- cbind("sort"=1:nrow(data), data)

data <- rbind(data, c(nrow(data)+1, NA, "rm_bydenoising",  counts-colSums(data[4:(ncol(data)-1)]), NA))

dir.create("_data/4_denoised")

#right place?
names(data) <- sub(renameSamples, "\\1", names(data))

write.csv(file="_data/4_denoised/A_Raw_haplotable.csv", data, row.names=F)

write.fasta(as.list(data$sequences[-nrow(data)]), paste(data$OTU[-nrow(data)], data$haplotype[-nrow(data)], sep="__"), "_data/4_denoised/A_Raw_haplo_sequ_byOTU.txt")
write.fasta(as.list(data$sequences[-nrow(data)]), data$haplotype[-nrow(data)], "_data/4_denoised/A_Raw_haplo_sequ.txt")

info <- paste("\nHaplotype table generated!\nRaw data and fasta files are available in _data/4_denoised (no subsetting)\n\nNow appling subsetting to the dataset!\n\nHaplotypes below ", minhaplosize, "% abundance (GLOBAL ABUNDANCE, across all haplotypes in a sample) in at least one sample are beeing discarded. The relative abundance is based on the number of sequences available before denoising (imput files).\nWaringing: All abundances in the table below ", minhaplosize, "% are set to 0. See Raw_haplotable.csv tble for orignial data without subsetting!\n\n", sep="")
message(info)
cat(file="../log.txt", info, append=T, sep="\n")


# apply subsetting!

data <- read.csv("_data/4_denoised/A_Raw_haplotable.csv", stringsAsFactors=F)

# minhaplosize (on each haplotype)
for (i in 1:(ncol(data)-4)){ # set to 0!
temp <- data[i+3]/sum(data[i+3])*100
data[nrow(data), i+3] <- data[nrow(data), i+3] + sum(data[i+3][temp < minhaplosize])
data[i+3][temp < minhaplosize] <- 0
}


# remove all rows with 0!
data <- data[rowSums(data[5:ncol(data)-1])!=0,]

data[nrow(data), 3] <- paste("denoised+below", minhaplosize, sep="")

write.csv(file=paste("_data/4_denoised/B_haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, ".csv", sep=""), data, row.names=F)

# ADD OTU subsetting based on OTUmin

info <- paste("\nSubsetting OTUs based on minimum relative abundance. OTUs with not at least ", OTUmin, "% abundance in ONE of the ", length(files)," samples are removed from the dataset!\n\n", sep="")
message(info)
cat(file="../log.txt", info, append=T, sep="\n")


data <- read.csv(file=paste("_data/4_denoised/B_haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, ".csv", sep=""), stringsAsFactors=F)

temp <- data[5:ncol(data)-1]
sums <- colSums(temp)
exp <- aggregate(unlist(temp[1]), list(data$OTU), "sum")[1]

i<-1
for (i in 1:ncol(temp)){
exp <- cbind(exp, aggregate(unlist(temp[i]), list(data$OTU), "sum")[2])
exp[i+1] <- exp[i+1]/sums[i]*100
exp[i+1] <- exp[i+1] >= OTUmin
}

keep <- rowSums(exp[-1])
keep2 <- exp$Group.1[keep>0]

discarded <- data[!data$OTU %in% keep2,] # adddiscarded reads to count
data <- data[data$OTU %in% keep2,]
data[nrow(data), (5:ncol(data)-1)] <- data[nrow(data), (5:ncol(data)-1)] + colSums(discarded[5:ncol(data)-1])



info <- paste(sum(keep>0), " of ", length(exp$Group.1), " OTUs remain in the dataset. ", sum(keep==0), " OTUs were discarded (", round(sum(keep==0)/length(exp$Group.1)*100, 2), "%)\n", sep="", "\nDiscarded OTUs are:\n", paste(exp$Group.1[keep==0], sep="", collapse=", "), "\n\n")
message(info)
cat(file="../log.txt", info, append=T, sep="\n")

round(sum(keep==0)/length(exp$Group.1)*100, 2)


#### end OTU subsetting, write final fasta files

#write.csv(file=paste("haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, ".csv", sep=""), data, row.names=F)
write.csv(file=paste("_data/4_denoised/C_haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, ".csv", sep=""), data, row.names=F)


write.fasta(as.list(data$sequences[-nrow(data)]), paste(data$OTU[-nrow(data)], data$haplotype[-nrow(data)], sep="__"), paste("_data/4_denoised/C_haplo_sequ_byOTU_", minhaplosize, "_minOTU_", OTUmin, ".txt", sep=""))
write.fasta(as.list(data$sequences[-nrow(data)]), data$haplotype[-nrow(data)], paste("_data/4_denoised/C_haplo_sequ_", minhaplosize, "_minOTU_", OTUmin, ".txt", sep=""))

centroids <- which(!duplicated(data$OTU))
centroids <- centroids[-length(centroids)]
write.fasta(as.list(data$sequences[centroids]), data$OTU[centroids], paste("_data/4_denoised/C_haplo_OTU_Centroids_", minhaplosize, "_minOTU_", OTUmin, ".txt", sep=""))


# "D" additionall filtering! remove haplotypes within OTUs and OTUs in indiv. samples with below X number of reads.
# withinOTU
# eachsampleOTUmin


otus <- unique(data$OTU)
otus <- otus[1:(length(otus)-1)] # don't considder last row (stats only)

# poduce table indicating what was deleted
PA_tab <- data

colSumsAbund <- colSums(data[4:(ncol(data)-1)]) # tracke lost reads
count_haplo <- 0
count_OTU <- 0

for(i in 4:(ncol(data)-1)){

for(k in 1:length(otus)){
select <- which(data$OTU== otus[k])
temp <- data[select, i]/sum(data[select, i])*100
temp[is.na(temp)] <- 0  # replace if NA because 0

#PA_tab
PA_tab[select[which(temp!=0&temp<withinOTU)], i] <- "low_Haplo"
count_haplo <- sum(temp!=0&temp<withinOTU)+count_haplo

temp[temp<withinOTU] <- 0 # delete haplotypes below treshhold within each OTU
temp[temp>0] <- data[select, i][temp>0]
data[select, i] <- temp


if(!is.null(eachsampleOTUmin)){ #delete OTU with to low abundance 
if(sum(data[select, i])<eachsampleOTUmin&sum(data[select, i])!=0){
data[select, i] <- 0
PA_tab[select, i] <- "low_OTU"
count_OTU <- count_OTU+1
}
} # end delete low abundance OTUs

}
}

# add lost sequences
data[nrow(data), 4:(ncol(data)-1)] <- colSumsAbund - colSums(data[-nrow(data), 4:(ncol(data)-1)])



# remove empty haplotypes!
# report some stats!
previous_haplo <- nrow(data)

data <- data[rowSums(data[4:(ncol(data)-1)])>0,]
PA_tab <- PA_tab[rowSums(data[4:(ncol(data)-1)])>0,]

report <- paste("\nRemoved ", previous_haplo - nrow(data), " of ", previous_haplo, " haplotypes (", round((previous_haplo - nrow(data))/previous_haplo*100, 2), "%), based on min haplo abundance withing OTU = ", withinOTU, " (withinOTU) and min number of sequences within each OTU in each sample = ", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin}, " (eachsampleOTUmin).\nnumber of intances where a haplotype in a sample was set to zero: ", count_haplo, "\nNumber of OTUs below ", " in individual samples: ", count_OTU,   sep="")
message(report)
cat(file="../log.txt", report, append=T, sep="\n")

# end additional filtering

# generate relative abuncnance haplotype table!
data_rel <- data

for(i in 4:(ncol(data)-1)){

for(k in 1:length(otus)){
select <- which(data$OTU== otus[k])
temp <- data[select, i]/sum(data[select, i])*100
temp[is.na(temp)] <- 0  # replace if NA because 0
data_rel[select, i] <- temp

}
}

# save files!

write.csv(file=paste("_data/4_denoised/C_haplotable_HIGHLIGHT_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},  ".csv", sep=""), PA_tab, row.names=F)

#write.csv(file=paste("D_haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},  ".csv", sep=""), data, row.names=F)
write.csv(file=paste("_data/4_denoised/D_haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},  ".csv", sep=""), data, row.names=F)

#write.csv(file=paste("D_haplotable_RELabund_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},  ".csv", sep=""), data_rel, row.names=F)
write.csv(file=paste("_data/4_denoised/D_haplotable_RELabund_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},  ".csv", sep=""), data_rel, row.names=F)


write.fasta(as.list(data$sequences[-nrow(data)]), paste(data$OTU[-nrow(data)], data$haplotype[-nrow(data)], sep="__"), paste("_data/4_denoised/D_haplo_sequ_byOTU_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin}, ".txt", sep=""))
write.fasta(as.list(data$sequences[-nrow(data)]), data$haplotype[-nrow(data)], paste("_data/4_denoised/D_haplo_sequ_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},".txt", sep=""))

centroids <- which(!duplicated(data$OTU))
centroids <- centroids[-length(centroids)]
write.fasta(as.list(data$sequences[centroids]), data$OTU[centroids], paste("_data/4_denoised/D_haplo_OTU_Centroids_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin}, ".txt", sep=""))


###########
###########
# additional subsetting recommended for LARGE datasets, OTUs and haplotypes have to be present in at least XXX samples.
info <- paste("\nCarryng out additional presence based subsetting (on OTUs / haplotypes). This is useful for large data sets to reduce noise!\nminHaploPresence = ", minHaploPresence, " (All haplotypes which are not present on more than ", minHaploPresence, " are discarded).\n", sep="")
message(info)
cat(file="../log.txt", info, append=T, sep="\n")


# cound and subset - haplo presence
temp <- data[,-ncol(data)] # rm sequences
temp[, -c(1:3)] <- temp[, -c(1:3)]>0 # T/F table
temp <- cbind(temp, "Nhaplo"=rowSums(temp[, -c(1:3)]))

ond_nrow <- nrow(data)-1
# remove all haplotypes present on only minHaploPresence locations

addme <- colSums(data[temp$Nhaplo< minHaploPresence, -c(1:3, ncol(data))]) # add discarded counts

# discarded sequences!
data <- data[temp$Nhaplo>= minHaploPresence,]

data[nrow(data), -c(1:3, ncol(data))] <- data[nrow(data), -c(1:3, ncol(data))]+ addme # add discarded counts to table


# add discarded haplotypes to highlight list!
# match highlight list
PA_tab$haplotype[is.na(PA_tab$haplotype)] <- "NA"
temp$haplotype[is.na(temp$haplotype)] <- "NA"

haplo_list <- rep("NA", length(PA_tab$haplotype))
haplo_low <- match(temp$haplotype, PA_tab$haplotype)
haplo_list[haplo_low] <- c(temp$Nhaplo>= minHaploPresence)


PA_tab2 <- data.frame(PA_tab[,-ncol(PA_tab)], "minHaploPresence"= haplo_list, PA_tab[,ncol(PA_tab)])
#write.csv(file="PA_tab2.csv", PA_tab2, row.names=F)


# haplo frequency
X <- data.frame(table(temp$Nhaplo[-length(temp$Nhaplo)]))
X$Var1 <- as.numeric(as.character(X$Var1))
Y <- data.frame("y"=c(1:c(ncol(data)-4)))
XY <- merge(X, Y, by.x="Var1", by.y="y", sort=T, all=T)
XY$Freq[is.na(XY$Freq)] <- 0

write.csv(file="_stats/E_Haplotypes_frequency.csv", XY, row.names=F)

pdf(file="_stats/E_Haplotype_presence.pdf")
hist(temp$Nhaplo[-length(temp$Nhaplo)], breaks=ncol(data)-4, xlim=c(1,ncol(data)), col="Gray", border=NA, xlab="Number of samples in which\nthe respective haplotype is present", main=paste(100-round(c(nrow(data)-1)/ond_nrow*100, 2), "% of haplotypes are discared!", sep=""))
lines(c(minHaploPresence, minHaploPresence), c(0, length(temp$Nhaplo)), col="Red")
dev.off()

info <- paste(ond_nrow-nrow(data), " of ", ond_nrow-1, " haplotypes discarded (", 100-round(c(nrow(data)-1)/(ond_nrow-1)*100, 2), "%), because they are present in less than ", minHaploPresence, " samples.\nA histogram plot showing the haplotype distribution across samples was written in the stats folder (_stats/E_Haplotype_presence.pdf).\n\n", sep="")
message(info)
cat(file="../log.txt", info, append=T, sep="\n")



# next subset
# remove all OTUs present in only 5 locations

# remove noisy_chimera / perfect_chimera

info <- "Discarding chimeras! Sequences discarded:"
Nold <- nrow(data)
data <- data[data$OTU!="noisy_chimera",]
info <- c(info, paste("noisy_chimera: ", Nold-nrow(data), sep=""))
Nold <- nrow(data)
data <- data[data$OTU!="perfect_chimera",]
info <- c(info, paste("perfect_chimera: ", Nold-nrow(data), sep=""))

message(paste(info, collapse="\n"), "\n")
cat(file="../log.txt", info, "\n", append=T, sep="\n")




OTU_list <- unique(data$OTU)
keep_list <- rep(length(OTU_list), F) # save OTUs to keep


for (i in 1:length(OTU_list)){

keep <- which(data$OTU==OTU_list[i])
mycols <- colSums(data[keep,-c(1:3,ncol(data))]>0)
OTUsum <- sum(mycols>0) # was at 1 before

keep_list <- c(keep_list, rep(OTUsum, length(keep)))
}

#nrow(data)

#data$OTU[!duplicated(data$OTU)]
pdf(file="_stats/E_OTU_presence_across_samples.pdf")
plot(keep_list[!duplicated(data$OTU)][-length(keep_list[!duplicated(data$OTU)])], xlab="OTU", ylab="Presence in N samples", ylim=c(0,43))
lines(c(-100, length(keep_list[!duplicated(data$OTU)])+100), c(minOTUPresence-0.5, minOTUPresence-0.5), col="Red")
dev.off()


# write highlight file!

PA_tab2$haplotype[is.na(PA_tab$haplotype)] <- "NA"
data$haplotype[is.na(data$haplotype)] <- "NA"



# OTUs keept
keepme <- keep_list>=minOTUPresence

OTU_list <- rep("NA", length(PA_tab2$haplotype))
OTU_low <- match(data$haplotype, PA_tab2$haplotype)
OTU_list[OTU_low] <- keepme

PA_tab3 <- data.frame(PA_tab2[,-ncol(PA_tab2)], "minOTUPresence"= OTU_list, "sequences"=PA_tab2[,ncol(PA_tab2)])

# write highlight table
write.csv(PA_tab3, "_data/4_denoised/E_highlight.csv", row.names=F)

# do actual OTU subsetting
oldN <- nrow(data)


#count number of discarded sequences


addme <- colSums(data[keep_list<minOTUPresence, -c(1:3, ncol(data))]) # add discarded counts

data <- data[keep_list>=minOTUPresence,] # OTU subsetting!

data[nrow(data), -c(1:3, ncol(data))] <- data[nrow(data), -c(1:3, ncol(data))]+ addme # add discarded counts to table



info <- paste("Discarting OUTs present in not at least ", minOTUPresence, " out of ", ncol(data)-4, " samples.\n", oldN - nrow(data), " of ", oldN, " OTUs discarded (", round(c(oldN - nrow(data))/oldN*100, 2), "%).\n", sep="")

message(info)
cat(file="../log.txt", info, "\n", append=T, sep="\n")

# write files

data_rel <- data

for(i in 4:(ncol(data)-1)){

for(k in 1:length(otus)){
select <- which(data$OTU== otus[k])
temp <- data[select, i]/sum(data[select, i])*100
temp[is.na(temp)] <- 0  # replace if NA because 0
data_rel[select, i] <- temp

}
}


write.csv(data, "_data/4_denoised/E_haplo_table.csv", row.names=F)
write.csv(data, "E_haplo_table.csv", row.names=F)

write.csv(data_rel, "_data/4_denoised/E_haplo_table_rel.csv", row.names=F)
write.csv(data_rel, "E_haplo_table_rel.csv", row.names=F)



write.fasta(as.list(data$sequences[-nrow(data)]), paste(data$OTU[-nrow(data)], data$haplotype[-nrow(data)], sep="__"), paste("_data/4_denoised/E_haplo_sequ_by_OTU.txt", sep=""))
write.fasta(as.list(data$sequences[-nrow(data)]), data$haplotype[-nrow(data)], paste("_data/4_denoised/E_haplo_sequ.txt", sep=""))

centroids <- which(!duplicated(data$OTU))
centroids <- centroids[-length(centroids)]
write.fasta(as.list(data$sequences[centroids]), data$OTU[centroids], paste("_data/4_denoised/E_haplo_OTU_Centroids.txt", sep=""))




temp <- "\nModule completed!"
message(temp)
cat(file="../log.txt", paste(Sys.time(), temp, "", sep="\n"), append=T, sep="\n")

setwd("../")


}
