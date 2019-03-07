# Haplotyping v0.1

Denoise <- function(files="latest",  strategy="unoise", unoise_alpha=5, minsize=10, minrelsize=0.0001, poolsamples=F, OTUmin=0.01, minhaplosize=0.003, withinOTU=5, eachsampleOTUmin=NULL, minHaploPresence=1, minOTUPresence=1, renameSamples="(.*_.*)_cut.*", exe="usearch"){



folder <- Core(module="Denoising")
cat(file="log.txt", c("Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files[1]=="latest"){
source(paste(folder, "/robots.txt", sep=""))
files <- list.files(paste(last_data, "/_data", sep=""), full.names=T)
}



# count sequences in each file
counts <- Count_sequences(files, fastq=F)
size <- round(counts* minrelsize/100) # get nim abundance
size[size<minsize] <- minsize # min size!



# Dereplicate files using Usearch
dir.create(paste(folder, "/_data/1_derep", sep=""))


new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub(".fasta", "_derep.fasta", new_names)
new_names <- paste(new_names, "_size_", size, ".txt", sep="")
new_names <- sub("_data", "_data/1_derep", new_names)
new_names <- paste(folder, "/", new_names, sep="")

cmd <- paste("-fastx_uniques \"", files, "\" -fastaout \"", new_names, "\" -sizeout", " -minuniquesize ", size,  sep="")

temp <- paste(length(files), " files are dereplicated and sequences in each sample below ", minrelsize, "% (or minuniqesize of ", minsize,")  are beeing discarded:", sep="")
cat(file="log.txt", temp , append=T, sep="\n")
message(temp)


temp <- new_names
for (i in 1:length(cmd)){
A <- system2(exe, cmd[i], stdout=T, stderr=T)
meep <- sub(".*_data/(.*)", "\\1", temp[i])
cat(file=paste(folder, "/_stats/1_derep_logs.txt", sep=""), paste("usearch ", cmd[i], sep="") , append=T, sep="\n")
cat(file=paste(folder, "/_stats/1_derep_logs.txt", sep=""), meep, A, "\n", append=T, sep="\n")


log_count <- Count_sequences(new_names[i], count_size=T)
log <- paste(sub(".*_data/1_derep/(.*)", "\\1", temp[i]), ": ", log_count, " of ", counts[i], " keept (", round((log_count/counts[i])*100, digits=4), "%, min size: ", size[i],")", sep="")
cat(file="log.txt", log , append=T, sep="\n")
message(log)
}


# pooled processing
# merge all files into one!

if(poolsamples){
cat(file=paste(folder, "/_stats/1_derep_logs.txt", sep=""), paste("\nCombining all files in a single file (samples_pooled.txt):\n", paste("cmd", cmd, collapse="", sep=""), collapse="", sep="") , append=T, sep="\n")
cat(file="log.txt", "\nCombining all files in a single file (samples_pooled.txt)\n", append=T, sep="\n")

# dereplicating pooled file
message("\nCombining all files in a single file (samples_pooled.txt)")
cmd <- paste(paste(paste("\"", new_names, "\"", sep=""), collapse=" "), "> ", folder, "/_data/1_derep/samples_pooled.txt", sep="")
system2("cat", cmd)

# dereplicating files
info <- "Dereplicating pooled sequences! (no min size)"
message(info)
cat(file="log.txt", info, append=T, sep="\n")

cmd <- paste("-fastx_uniques \"", folder, "/_data/1_derep/samples_pooled.txt\" -fastaout \"", folder, "/_data/1_derep/samples_pooled_derep.txt\" -sizein -sizeout", sep="")
A <- system2(exe, cmd, stdout=T, stderr=T)

cat(file=paste(folder, "/_stats/1_derep_logs.txt", sep=""), paste("usearch", cmd, sep=""), append=T, sep="\n")
cat(file=paste(folder, "/_stats/1_derep_logs.txt", sep=""), A, append=T, sep="\n")

# renaming all sequences!
info <- "Renaming pooled sequences, and applying same names to the dereplicated files.\n"
message(info)
cat(file="log.txt", info, append=T, sep="\n")


haplo <- read.fasta(paste(folder, "/_data/1_derep/samples_pooled_derep.txt", sep=""), forceDNAtolower=F, as.string=T)

temp <- sub(".*(;size=.*;)", "\\1", names(haplo))
temp2 <- paste("haplo_", 1:length(haplo), temp, sep="")

write.fasta(haplo, temp2, paste(folder, "/_data/1_derep/samples_pooled_derep_renamed.txt", sep=""))


# rename single files!
dir.create(paste(folder, "/_data/2_renamed", sep=""))
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
cat(file="log.txt", info, append=T, sep="\n")

cmd <- paste("-unoise3 \"", folder, "/_data/1_derep/samples_pooled_derep_renamed.txt\" -zotus \"", folder, "/_data/1_derep/samples_pooled_+_denoised.txt\" -unoise_alpha ", unoise_alpha,  sep="")

A <- system2(exe, cmd, stdout=T, stderr=T)
cat(file=paste(folder, "/_stats/2_unoise.txt", sep=""), c(info, "", paste("usearch", cmd), "", A), append=T, sep="\n")

info <- paste("Denoising compelte! ", Count_sequences(paste(folder, "/_data/1_derep/samples_pooled_derep_renamed.txt", sep=""), fastq=F), " sequences were denoised using ", strategy, ".", "\nA total of ", sub(".*100.0% (.*) good, .* chimeras\r", "\\1", A[length(A)-1]), " haplotypes remained after denoising!\n", sep="")
message(info)
cat(file="log.txt", info, append=T, sep="\n")



#Zotus, get old names back (original haplotypes)!

Zotus <- read.fasta(paste(folder, "/_data/1_derep/samples_pooled_+_denoised.txt", sep=""), as.string=T, forceDNAtolower=F)
renamed_sequ <- read.fasta(paste(folder, "/_data/1_derep/samples_pooled_derep_renamed.txt", sep=""), as.string=T, forceDNAtolower=F)

matched <- match(Zotus, renamed_sequ)
new_sample <- renamed_sequ[matched] # DNA sequences

write.fasta(new_sample, names(new_sample), paste(folder, "/_data/1_derep/samples_pooled_+_denoised_renamed.txt", sep=""))


names(new_sample) <- sub(";size(.*);", "", names(new_sample))
haplotypes <- new_sample

dir.create(paste(folder, "/_data/3_unoise", sep=""))

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
cat(file="log.txt", info, append=T, sep="\n")

}
}
# end denoising POOLED samples


# If denoising on individual sampels!
if(!poolsamples){

dir.create(paste(folder, "/_data/2_denoised", sep=""))

denoised <- new_names
denoised <- sub("1_derep", "2_denoised", denoised)


cmd <- paste("-unoise3 \"", new_names, "\" -zotus \"", denoised,"\" -unoise_alpha ", unoise_alpha, " -sizein -sizeout", sep="")

# check for empty dereplicated files
empty_filesTF <- file.info(new_names)[,1]>0

tempM <- paste("A toatal of ", sum(!empty_filesTF)," of ", length(empty_filesTF), " files are empty and will not be denoised!\n", paste(sub(".*_data/1_derep/", "", new_names)[!empty_filesTF], collapse="\n"), sep="")
message(tempM)
cat(file="log.txt", tempM, append=T, sep="\n")

empty_files <- new_names[!empty_filesTF]


for (i in which(empty_filesTF)){

A <- system2(exe, cmd[i], stdout=T, stderr=T)

cat(file=paste(folder, "/_stats/2_denoise_logs.txt", sep=""), paste("usearch", cmd[i], sep=""), append=T, sep="\n")
cat(file=paste(folder, "/_stats/2_denoise_logs.txt", sep=""), A, append=T, sep="\n")

seqin <- Count_sequences(new_names[i], fastq=F)
seqout <- Count_sequences(denoised[i], fastq=F)

temp <- paste("Sample ", sub(".*/(.*)_PE_.*", "\\1", new_names[i]), " denoised ", seqin, " sequences to ", seqout, " ESVs (", round(seqout/seqin*100, 2), "% keeped)", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")
}

# Include abundance information again

for (i in which(empty_filesTF)){
de <- read.fasta(denoised[i], as.string=T, forceDNAtolower=F)
fast <- read.fasta(new_names[i], as.string=T, forceDNAtolower=F)
matched <- match(fast, de)
matched <- matched[!is.na(matched)]

names(de) <- names(fast[matched])
write.fasta(de, names(de), file.out= denoised[i])
}

} # Done with indiv procesing.


# Cluster into OTUs (for OTU table information)
if(poolsamples){
cmd <- paste(" -cluster_otus ", folder, "/_data/1_derep/samples_pooled_+_denoised_renamed.txt -otus ", folder, "/_data/1_derep/samples_pooled_+_denoised_renamed_OTUsequ.txt -uparseout ", folder, "/_data/1_derep/samples_pooled_+_denoised_renamed_OTUtable.txt -relabel OTU_ -strand plus", sep="")

A <- system2(exe, cmd, stdout=T, stderr=T) # cluster OTUs!

cat(file=paste(folder, "/_stats/2_unoise.txt", sep=""), c("Clustering haplotypes into OTUs for OTU table!", "", paste("usearch", cmd), "", A), append=T, sep="\n")

chimeras <- as.numeric(sub(".*100.0% .* OTUs, (.*) chimeras\r", "\\1", A[grep("chimeras\r", A)]))
OTUs <- as.numeric(sub(".*100.0% (.*) OTUs, .* chimeras\r", "\\1", A[grep("chimeras\r", A)]))
if(is.na(chimeras)){chimeras<-0}

info <- paste("Clustered ", length(haplotypes), " haplotype sequences (cluster_otus, 3% simmilarity) into ", OTUs, " OTUs (+", chimeras, " chimeras).\nOTUs and (potentially) chimeric sequences will be included in the Haplotype table!\n", sep="" )
message(info)
cat(file="log.txt", info, append=T, sep="\n")
} # end pooled processing


# OTU clustering for indiv samples
if(!poolsamples){

dir.create(paste(folder, "/_data/3_pooledESV", sep=""))




cat(file=paste(folder, "/_stats/1_derep_logs.txt", sep=""), paste("\nCombining all files in a single file (samples_pooled.txt):\n", paste("cmd", cmd, collapse="", sep=""), collapse="", sep="") , append=T, sep="\n")
cat(file="log.txt", "\nCombining all files in a single file (samples_pooled.txt)\n", append=T, sep="\n")

# dereplicating pooled file
message("\nCombining all individually denoised files in a single file (samples_pooled.txt)")
cmd <- paste(paste(paste("\"", denoised[empty_filesTF], "\"", sep=""), collapse=" "), "> ", folder, "/_data/3_pooledESV/1_samples_pooled.txt", sep="")
system2("cat", cmd)

# dereplicating files
info <- "Dereplicating pooled sequences!"
message(info)
cat(file="log.txt", info, append=T, sep="\n")

cmd <- paste("-fastx_uniques \"", folder, "/_data/3_pooledESV/1_samples_pooled.txt\" -fastaout \"", folder, "/_data/3_pooledESV/2_samples_pooled_derep.txt\" -sizein -sizeout", sep="")
A <- system2(exe, cmd, stdout=T, stderr=T)

cat(file=paste(folder, "/_stats/1_derep_logs.txt", sep=""), paste("usearch", cmd, sep=""), append=T, sep="\n")
cat(file=paste(folder, "/_stats/1_derep_logs.txt", sep=""), A, append=T, sep="\n")


# rename files for ESV
ESVs <- read.fasta(paste(folder, "/_data/3_pooledESV/2_samples_pooled_derep.txt", sep=""), forceDNAtolower=F, as.string=T)



names(ESVs) <- paste("haplo_", 1:length(ESVs), sub(".*(;size=.*;)", "\\1", names(ESVs)), sep="")
write.fasta(ESVs, names=names(ESVs), paste(folder, "/_data/3_pooledESV/3_ESV_list.txt", sep=""))


# rename ESVs in dereplicated files
ESV_list <- read.fasta(paste(folder, "/_data/3_pooledESV/3_ESV_list.txt", sep=""), forceDNAtolower=F, as.string=T)


FileListDenoised <- list.files(paste(folder, "/_data/2_denoised", sep=""), full.names=T)

dir.create(paste(folder, "/_data/3_unoise", sep=""))

for (i in 1:length(FileListDenoised)){

temp <- read.fasta(FileListDenoised[i], forceDNAtolower=F, as.string=T)

temp_size <- sub(".*(;size=.*;)", "\\1", names(temp))

names(temp) <- sub(";size=.*;", "", names(ESV_list)[match(temp, ESV_list)])

names(temp) <- paste(names(temp), temp_size, sep="")

temp_filename <- sub("2_denoised", "3_unoise", FileListDenoised[i])
temp_filename <- sub(".txt$", "_denoised.txt", temp_filename)

write.fasta(temp, names=names(temp), temp_filename)

} # end renaming



cmd <- paste(" -cluster_otus ", folder, "/_data/3_pooledESV/3_ESV_list.txt -otus ", folder, "/_data/1_derep/samples_pooled_+_denoised_renamed_OTUsequ.txt -uparseout ", folder, "/_data/1_derep/samples_pooled_+_denoised_renamed_OTUtable.txt -relabel OTU_ -strand plus", sep="")

A <- system2(exe, cmd, stdout=T, stderr=T) # cluster OTUs!

cat(file=paste(folder, "/_stats/2_unoise.txt", sep=""), c("Clustering haplotypes into OTUs for OTU table!", "", paste("usearch", cmd), "", A), append=T, sep="\n")

chimeras <- as.numeric(sub(".*100.0% .* OTUs, (.*) chimeras\r", "\\1", A[grep("chimeras\r", A)]))
OTUs <- as.numeric(sub(".*100.0% (.*) OTUs, .* chimeras\r", "\\1", A[grep("chimeras\r", A)]))
if(is.na(chimeras)){chimeras<-0}

info <- paste("Clustered ", length(ESV_list), " haplotype sequences (cluster_otus, 3% simmilarity) into ", OTUs, " OTUs (+", chimeras, " chimeras).\nOTUs and (potentially) chimeric sequences will be included in the Haplotype table!\n", sep="" )
message(info)
cat(file="log.txt", info, append=T, sep="\n")

haplotypes <- ESV_list
denoised_sequences <- list.files(paste(folder, "/_data/3_unoise", sep=""), full.names=T)


# need to keep working on this

} # end processing indiv denoised files




# generate one united haplotype table!

OTUs <- read.csv(paste(folder, "/_data/1_derep/samples_pooled_+_denoised_renamed_OTUtable.txt", sep=""), stringsAsFactors=F, sep="\t", header=F)

k <- 1
OTU_list <- NULL
for (i in 1:nrow(OTUs)){


if(grepl("otu", OTUs$V2[i])){OTU_list[i] <- paste("OTU_", k, sep="")
k <- k+1} else
if(OTUs$V2[i]=="match"){OTU_list[i] <- sub(".*;top=(OTU_.*)\\(.*", "\\1", OTUs$V3[i])} else {OTU_list[i] <- OTUs$V2[i]}

}


data <- data.frame("haplotype"=sub(";size=.*;", "",names(haplotypes)), "OTU"=OTU_list, stringsAsFactors=F)

for (i in 1:length(denoised_sequences)){
sample <- names(read.fasta(denoised_sequences[i]))
matched <- match(sub(";size=.*;", "", sample), data$haplotype)
abundance <- rep(0, nrow(data))
abundance[matched] <- as.numeric(sub(".*;size=(.*);", "\\1", sample))

data <- cbind(data, abundance)
names(data)[i+2] <- sub(".*_data/3_unoise/(.*)_denoised.txt", "\\1", denoised_sequences[i])
}
#head(data)
data <- cbind(data, "sequences"=unlist(haplotypes), stringsAsFactors=F)
# sort by OTUs
data <- data[order(suppressWarnings(as.numeric(sub("OTU_", "", data$OTU)))),]
data <- cbind("sort"=1:nrow(data), data)

data <- rbind(data, c(nrow(data)+1, NA, "rm_bydenoising",  counts[which(empty_filesTF)]-colSums(data[4:(ncol(data)-1)]), NA))

dir.create(paste(folder, "/_data/4_denoised", sep=""))

#right place?
names(data) <- sub(renameSamples, "\\1", names(data))

# Add back in empty cells
tempname <- sub(".*_data/1_derep/", "", new_names[which(!empty_filesTF)])
tempname <- sub(renameSamples, "\\1", tempname)

data2 <- data[,-c(ncol(data))]

for (i in 1:length(tempname)){
data2 <- data.frame(data2, rep(0, nrow(data2)), stringsAsFactors=F)
names(data2)[ncol(data2)] <- tempname[i]
}

data3 <- data2[,-c(1:3)]
data4 <-data3[,order(names(data3))]

data5 <- cbind(data[,c(1:3)], data4, data[,ncol(data)])

names(data5)[ncol(data5)] <- "sequences"
data <- data5

write.csv(file=paste(folder, "/_data/4_denoised/A_Raw_haplotable.csv", sep=""), data, row.names=F)

write.fasta(as.list(data$sequences[-nrow(data)]), paste(data$OTU[-nrow(data)], data$haplotype[-nrow(data)], sep="__"), paste(folder, "/_data/4_denoised/A_Raw_haplo_sequ_byOTU.txt", sep=""))
write.fasta(as.list(data$sequences[-nrow(data)]), data$haplotype[-nrow(data)], paste(folder, "/_data/4_denoised/A_Raw_haplo_sequ.txt", sep=""))

info <- paste("\nHaplotype table generated!\nRaw data and fasta files are available in _data/4_denoised (no subsetting)\n\nNow appling subsetting to the dataset!\n\nHaplotypes below ", minhaplosize, "% abundance (GLOBAL ABUNDANCE, across all haplotypes in a sample) in at least one sample are beeing discarded. The relative abundance is based on the number of sequences available before denoising (imput files).\nWaringing: All abundances in the table below ", minhaplosize, "% are set to 0. See Raw_haplotable.csv tble for orignial data without subsetting!\n\n", sep="")
message(info)
cat(file="log.txt", info, append=T, sep="\n")


# apply subsetting!

data <- read.csv(paste(folder, "/_data/4_denoised/A_Raw_haplotable.csv", sep=""), stringsAsFactors=F)

# minhaplosize (on each haplotype)
for (i in 1:(ncol(data)-4)){ # set to 0!
temp <- data[i+3]/sum(data[i+3])*100
data[nrow(data), i+3] <- data[nrow(data), i+3] + sum(data[i+3][temp < minhaplosize])
data[i+3][temp < minhaplosize] <- 0
}


# remove all rows with 0!
data <- data[rowSums(data[5:ncol(data)-1])!=0,]

data[nrow(data), 3] <- paste("denoised+below", minhaplosize, sep="")

write.csv(file=paste(folder,"/_data/4_denoised/B_haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, ".csv", sep=""), data, row.names=F)

# ADD OTU subsetting based on OTUmin

info <- paste("\nSubsetting OTUs based on minimum relative abundance. OTUs with not at least ", OTUmin, "% abundance in ONE of the ", length(files)," samples are removed from the dataset!\n\n", sep="")
message(info)
cat(file="log.txt", info, append=T, sep="\n")


data <- read.csv(file=paste(folder, "/_data/4_denoised/B_haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, ".csv", sep=""), stringsAsFactors=F)

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
cat(file="log.txt", info, append=T, sep="\n")

round(sum(keep==0)/length(exp$Group.1)*100, 2)


#### end OTU subsetting, write final fasta files

#write.csv(file=paste("haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, ".csv", sep=""), data, row.names=F)
write.csv(file=paste(folder, "/_data/4_denoised/C_haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, ".csv", sep=""), data, row.names=F)


write.fasta(as.list(data$sequences[-nrow(data)]), paste(data$OTU[-nrow(data)], data$haplotype[-nrow(data)], sep="__"), paste(folder, "/_data/4_denoised/C_haplo_sequ_byOTU_", minhaplosize, "_minOTU_", OTUmin, ".txt", sep=""))
write.fasta(as.list(data$sequences[-nrow(data)]), data$haplotype[-nrow(data)], paste(folder, "/_data/4_denoised/C_haplo_sequ_", minhaplosize, "_minOTU_", OTUmin, ".txt", sep=""))

centroids <- which(!duplicated(data$OTU))
centroids <- centroids[-length(centroids)]
write.fasta(as.list(data$sequences[centroids]), data$OTU[centroids], paste(folder, "/_data/4_denoised/C_haplo_OTU_Centroids_", minhaplosize, "_minOTU_", OTUmin, ".txt", sep=""))


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

report <- paste("\nRemoved ", c(previous_haplo-1) - c(nrow(data)-1), " of ", previous_haplo, " haplotypes (", round((c(previous_haplo-1) - (nrow(data)-1))/c(previous_haplo-1)*100, 2), "%), based on min haplo abundance within OTU = ", withinOTU, " (withinOTU) and min number of sequences within each OTU in each sample = ", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin}, " (eachsampleOTUmin).\nnumber of intances where a haplotype in a sample was set to zero: ", count_haplo, "\nNumber of OTUs below ", " in individual samples: ", count_OTU,   sep="")
message(report)
cat(file="log.txt", report, append=T, sep="\n")

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

write.csv(file=paste(folder, "/_data/4_denoised/C_haplotable_HIGHLIGHT_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},  ".csv", sep=""), PA_tab, row.names=F)

#write.csv(file=paste("D_haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},  ".csv", sep=""), data, row.names=F)
write.csv(file=paste(folder, "/_data/4_denoised/D_haplotable_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},  ".csv", sep=""), data, row.names=F)

#write.csv(file=paste("D_haplotable_RELabund_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},  ".csv", sep=""), data_rel, row.names=F)
write.csv(file=paste(folder, "/_data/4_denoised/D_haplotable_RELabund_alpha_", unoise_alpha, "_haplosize_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},  ".csv", sep=""), data_rel, row.names=F)


write.fasta(as.list(data$sequences[-nrow(data)]), paste(data$OTU[-nrow(data)], data$haplotype[-nrow(data)], sep="__"), paste(folder, "/_data/4_denoised/D_haplo_sequ_byOTU_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin}, ".txt", sep=""))
write.fasta(as.list(data$sequences[-nrow(data)]), data$haplotype[-nrow(data)], paste(folder, "/_data/4_denoised/D_haplo_sequ_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin},".txt", sep=""))

centroids <- which(!duplicated(data$OTU))
centroids <- centroids[-length(centroids)]
write.fasta(as.list(data$sequences[centroids]), data$OTU[centroids], paste(folder, "/_data/4_denoised/D_haplo_OTU_Centroids_", minhaplosize, "_minOTU_", OTUmin, "_withinOTU_", withinOTU, "_eachsampleOTUmin_", if(is.null(eachsampleOTUmin)){"NULL"}else{eachsampleOTUmin}, ".txt", sep=""))


###########
###########
# additional subsetting recommended for LARGE datasets, OTUs and haplotypes have to be present in at least XXX samples.
info <- paste("\nCarryng out additional presence based subsetting (on OTUs / haplotypes). This is useful for large data sets to reduce noise!\nminHaploPresence = ", minHaploPresence, " (All haplotypes which are not present on more than ", minHaploPresence, " are discarded).\n", sep="")
message(info)
cat(file="log.txt", info, append=T, sep="\n")


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

write.csv(file=paste(folder, "/_stats/E_Haplotypes_frequency.csv", sep=""), XY, row.names=F)


pdf(file=paste(folder, "/_stats/E_Haplotype_presence.pdf", sep=""))
plot(NULL, ylim=c(0, max(XY$Freq)*1.05), xlim=c(0.5,ncol(data)-3.5), xlab="Number of samples in which\nthe respective haplotype is present", main=paste(100-round(c(nrow(data)-1)/ond_nrow*100, 2), "% of haplotypes are discared!", sep=""), xaxt="n", ylab="Haplotype abundance")
axis(1, c(1:c(ncol(data)-4)), labels=c(1:c(ncol(data)-4)))
for (i in 1:length(XY$Freq)){
rect(i-0.5, 0, i+0.5, XY$Freq[i], col="gray", border=NA)
}

lines(c(minHaploPresence-0.5, minHaploPresence-0.5), c(0, length(temp$Nhaplo)), col="Red")
dev.off()

info <- paste((ond_nrow)-(nrow(data)-1), " of ", ond_nrow-1, " haplotypes discarded (", round(c((ond_nrow)-(nrow(data)-1))/(ond_nrow)*100, 2), "%), because they are present in less than ", minHaploPresence, " samples.\nA histogram plot showing the haplotype distribution across samples was written in the stats folder (_stats/E_Haplotype_presence.pdf).\n\n", sep="")
message(info)
cat(file="log.txt", info, append=T, sep="\n")



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
cat(file="log.txt", info, "\n", append=T, sep="\n")




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
pdf(file=paste(folder, "/_stats/E_OTU_presence_across_samples.pdf", sep=""))
plot(keep_list[!duplicated(data$OTU)][-length(keep_list[!duplicated(data$OTU)])], xlab="OTU", ylab="Presence in N samples", ylim=c(0,ncol(data)-4))
lines(c(-100, length(keep_list[!duplicated(data$OTU)])+100), c(minOTUPresence-0.5, minOTUPresence-0.5), col="Red")
dev.off()
tail(data)

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
write.csv(PA_tab3, paste(folder, "/_data/4_denoised/E_highlight.csv", sep=""), row.names=F)

# do actual OTU subsetting
oldN <- nrow(data)


#count number of discarded sequences


addme <- colSums(data[keep_list<minOTUPresence, -c(1:3, ncol(data))]) # add discarded counts

data <- data[keep_list>=minOTUPresence,] # OTU subsetting!

data[nrow(data), -c(1:3, ncol(data))] <- data[nrow(data), -c(1:3, ncol(data))]+ addme # add discarded counts to table



info <- paste("Discarting OUTs present in not at least ", minOTUPresence, " out of ", ncol(data)-4, " samples.\n", oldN - nrow(data), " of ", oldN, " OTUs discarded (", round(c(oldN - nrow(data))/oldN*100, 2), "%).\n", sep="")

message(info)
cat(file="log.txt", info, "\n", append=T, sep="\n")

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


write.csv(data, paste(folder, "/_data/4_denoised/E_haplo_table.csv", sep=""), row.names=F)
write.csv(data, paste(folder, "/E_haplo_table.csv", sep=""), row.names=F)

write.csv(data_rel, paste(folder, "/_data/4_denoised/E_haplo_table_rel.csv", sep=""), row.names=F)
write.csv(data_rel, paste(folder, "/E_haplo_table_rel.csv", sep=""), row.names=F)



write.fasta(as.list(data$sequences[-nrow(data)]), paste(data$OTU[-nrow(data)], data$haplotype[-nrow(data)], sep="__"), paste(folder, "/_data/4_denoised/E_haplo_sequ_by_OTU.txt", sep=""))
write.fasta(as.list(data$sequences[-nrow(data)]), data$haplotype[-nrow(data)], paste(folder, "/_data/4_denoised/E_haplo_sequ.txt", sep=""))

centroids <- which(!duplicated(data$OTU))
centroids <- centroids[-length(centroids)]
write.fasta(as.list(data$sequences[centroids]), data$OTU[centroids], paste(folder, "/_data/4_denoised/E_haplo_OTU_Centroids.txt", sep=""))


Denoise_barplot(paste(folder, "/E_haplo_table.csv", sep=""), out=paste(folder, "/_stats/Haplotype_barplot_E.pdf", sep=""), emptyOTUs=T)
message(paste(sep="Haplotype distribution plot generated:\n", folder, "/_stats/Haplotype_barplot_E.pdf\n\nYou may use the fuction \"Denoise_barplot\"(",folder, "/E_haplo_table.csv) to further customise that plot.\n\n"))


temp <- "\nModule completed!"
message(temp)
cat(file="log.txt", paste(Sys.time(), temp, "", "*** Module completed!", sep="\n"), append=T, sep="\n")


}
