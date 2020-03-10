# U_cluster_otus v0.1

Cluster_otus <- function(files="latest", minuniquesize=2, strand="plus", filter=0.01, filterN=1, exe="vsearch", otu_radius_pct=3, mapp_singletons=T, maxaccepts=1, maxrejects=32, delete_data=T, heatmap=T, threads=NA){
#, unoise_min=NA - unoise denoising removed, no longer supported!

folder <- Core(module="Cluster_otus", delete_data=delete_data)
cat(file="log.txt", c("Version v0.2", "\n"), append=T, sep="\n")
message(" ")

files_to_delete <- NULL

if(exe=="usearch"){
A <- system2(exe, stdout=T, stderr=T)
version <- as.numeric(sub("usearch v(.+)\\.+.*\\..*_.*", "\\1", A[1]))

if(otu_radius_pct!=3){
if(version > 8){
temp <-  paste("WARNING: You did set a custom clustering treshold of ", otu_radius_pct, " but are using Usearch version ", version, "! The custom treshold was removed with version 9 (nothing I can do about this), thus please provide a path in exe to Usearch8, to make this work. You can download the older version from the Usearch website or use vsearch. The script will stop now.", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")
stop("Script stopped!")
}
}
}


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

# Dereplicate files using Vsearch
if(mapp_singletons){
dir.create(paste(folder, "/_data/1_derep_inc_singletons", sep=""))
} else {
dir.create(paste(folder, "/_data/1_derep_minsize_", minuniquesize, sep=""))
}

new_names <- sub(".*(_data/.*)", "\\1", files[!empty])
new_names <- sub("_PE.*", "_PE_derep.fasta", new_names)
if(mapp_singletons){
new_names <- sub("_data", "_data/1_derep_inc_singletons", new_names)
}else{
new_names <- sub("_data", paste("_data/1_derep_minsize_", minuniquesize, sep=""), new_names)
}
new_names <- paste(folder, "/", new_names, sep="")

if(exe=="usearch"){
# might need to be changed!
cmd <- paste("-derep_fulllength \"", files[!empty], "\" -output \"", new_names, "\" -sizeout -sizein", if(! mapp_singletons){paste(" -minuniquesize ", minuniquesize, sep="")}, sep="")
}else{
cmd <- paste("-derep_fulllength \"", files[!empty], "\" -output \"", new_names, "\" -sizeout -sizein", if(! mapp_singletons){paste(" -minuniquesize ", minuniquesize, sep="")}, if(!is.na(threads)){paste(" -threads ", threads, sep="")}, sep="")
}

files_to_delete <- c(files_to_delete, new_names)

if(!mapp_singletons){
temp <- paste(length(cmd), " files are dereplicated using vsearch, but since mapp_singletons = F, sequences below minuniquesize = ", minuniquesize, " are discarded in each sample;", sep="")
} else {
temp <- paste(length(cmd), " files are dereplicated (incl. singletons!):", sep="")}
cat(file="log.txt", temp , append=T, sep="\n")
message(temp)


temp <- new_names
for (i in 1:length(cmd)){
A <- system2(exe, cmd[i], stdout=T, stderr=T)
meep <- sub(".*_data/(.*)", "\\1", temp[i])
cat(file="log.txt", meep, append=T, sep="\n")
cat(file=paste(folder, "/_stats/1_derep_logs.txt", sep=""), meep, A, "\n", append=T, sep="\n")
message(meep)
}


# 2 make OTUs!
# merge all files into one

dir.create(paste(folder, "/_data/2_OTU_clustering", sep=""))

# HOT FIX 190227 - split cat
if(length(new_names)<800){

cmd <- paste(paste(new_names, collapse=" "), " > ", folder, "/_data/2_OTU_clustering/A_all_files_united.fasta", collapse="", sep="")
A <- system2("cat", cmd, stdout=T, stderr=T)


files_to_delete <- c(files_to_delete, paste(folder, "/_data/2_OTU_clustering/A_all_files_united.fasta", sep="")) } else {
# SPLIT
# TEMP A
cmd <- paste(paste(new_names[1:700], collapse=" "), " > ", folder, "/_data/2_OTU_clustering/A_all_files_TEMP_A.fasta", collapse="", sep="")
A <- system2("cat", cmd, stdout=T, stderr=T)
# TEMP B
cmd <- paste(paste(new_names[701:length(new_names)], collapse=" "), " > ", folder, "/_data/2_OTU_clustering/A_all_files_TEMP_B.fasta", collapse="", sep="")
A <- system2("cat", cmd, stdout=T, stderr=T)

message(length(new_names[1:700]))
message(length(new_names[701:length(new_names)]))

cmd <- paste(folder, "/_data/2_OTU_clustering/A_all_files_TEMP_A.fasta ", folder, "/_data/2_OTU_clustering/A_all_files_TEMP_B.fasta", " > ", folder, "/_data/2_OTU_clustering/A_all_files_united.fasta", collapse="", sep="")
A <- system2("cat", cmd, stdout=T, stderr=T)



} # end hot fix



#check <- readLines("_data/2_OTU_clustering/A_all_files_united.fasta")
#count <- as.numeric(sub(".*size=(.*)", "\\1", check))
#sum(count, na.rm=T)

# write logs

temp <- paste(length(new_names), " dereplicated files where merged ", if(mapp_singletons==T){"(inc singleotns) "}, "into file:\n\"",folder, "/_data/2_OTU_clustering/A_all_files_united.fasta\"", sep="")
message("\n", temp)
cat(file="log.txt", "\n", temp, append=T, sep="\n")


cat(file=paste(folder, "/_stats/2_OTU_clustering_log.txt", sep=""), temp, "", paste("cat", cmd), append=T, sep="\n")

# dereplicate "A_all_files_united.fasta" using Vsearch!

if(exe=="usearch"){
# Needs to be updated
cmd <- paste("-derep_fulllength ", folder, "/_data/2_OTU_clustering/A_all_files_united.fasta -output ", folder, "/_data/2_OTU_clustering/B_all_derep_min", minuniquesize, ".fasta -sizein -sizeout -minuniquesize ", minuniquesize, sep="")
} else {
cmd <- paste("-derep_fulllength ", folder, "/_data/2_OTU_clustering/A_all_files_united.fasta -output ", folder, "/_data/2_OTU_clustering/B_all_derep_min", minuniquesize, ".fasta -sizein -sizeout -minuniquesize ", minuniquesize, sep="")
}

filename_all_unique <- paste("B_all_derep_min", minuniquesize, ".fasta", sep="")

files_to_delete <- c(files_to_delete, paste(folder, "/_data/2_OTU_clustering/B_all_derep_min2.fasta", sep=""))


A <- system2(exe, cmd, stdout=T, stderr=T)

temp <- paste("Total number of sequences (not dereplicated): ", sub(".*nt in (.*) seqs.*", "\\1", A[grep("seqs, min", A)]), "\n", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

temp <- paste("Pooled sequences are dereplicated with minuniquesize = ", minuniquesize , " into a total of ", sub("(.*) unique sequences.*", "\\1", A[grep(" unique sequences", A)]), " unique sequences.", "\n", "File prepared for OTU clustering: \"", filename_all_unique, "\"", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

# derep log
cat(file=paste(folder, "/_stats/2_OTU_clustering_log.txt", sep=""), "\n", paste("cat", cmd), "\n", A, "", append=T, sep="\n")

# Actual clustering of dereplicated file 
OTU_file <- sub(".fasta", "_OTUs.fasta", filename_all_unique)
OTU_file <- sub("B_", "C_", OTU_file)

# add option to cluster OTUs at other tahn 3% (only wrks with usearch8 or lower)

if(exe=="usearch"){
cmd <- paste(" -cluster_otus ", folder, "/_data/2_OTU_clustering/", filename_all_unique, " -otus ", folder, "/_data/2_OTU_clustering/", OTU_file, " -uparseout ", folder, "/_data/2_OTU_clustering/", sub(".fasta", "_OTUtab.txt", OTU_file), " -relabel OTU_ -strand ", strand, if(otu_radius_pct!=3){paste(" -otu_radius_pct ", otu_radius_pct, sep="")}, sep="")

A <- system2(exe, cmd, stdout=T, stderr=T) # cluster OTUs!

cat(file=paste(folder, "/_stats/2_OTU_clustering_log.txt", sep=""), "\n", paste("usearch", cmd), "\n", A, "", append=T, sep="\n")


#files_to_delete <- c(files_to_delete, paste(folder, "/_data/2_OTU_clustering/", OTU_file, sep=""))
files_to_delete <- c(files_to_delete, paste(folder, "/_data/2_OTU_clustering/", sub(".fasta", "_OTUtab.txt", OTU_file), sep=""))


if(version!=8){
chimeras <- sub(".*OTUs, (.*) chimeras\r", "\\1", A[grep("chimeras\r", A)])
OTUs <- sub(".*100.0% (.*) OTUs, .* chimeras\r", "\\1", A[grep("chimeras\r", A)])
} else {
chimeras <- sub(".*OTUs, (.*) chimeras.*", "\\1", A[grep("chimeras.*", A)])
OTUs <- sub(".*100.0% (.*) OTUs, .* chimeras.*", "\\1", A[grep("chimeras.*", A)])
}

temp <- paste("\n", "Clustering reads from \"", filename_all_unique, "\nminuniquesize = ", minuniquesize, "\nstrand = ", strand, "\nChimeras discarded: ", chimeras, "\nOTUs written: ", OTUs, " -> file \"", OTU_file, "\"\n", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")
} else {
# Clustering with vsearch! Needs extra denoising. 200309
# sequencesa re already sorted by derelication before
cmd <- paste(" -cluster_smallmem ", folder, "/_data/2_OTU_clustering/", filename_all_unique, " -centroids ", folder, "/_data/2_OTU_clustering/", sub("OTUs.fasta", "OTUs+chimeras.fasta", OTU_file), " -strand ", strand, " -usersort -id ", (100-otu_radius_pct)*0.01, sep="", " -uc ", folder, "/_data/2_OTU_clustering/", sub(".fasta", "_OTUtab.txt", OTU_file))

A <- system2(exe, cmd, stdout=T, stderr=T) # cluster OTUs!

cat(file=paste(folder, "/_stats/2_OTU_clustering_log.txt", sep=""), "\n", paste("vsearch", cmd), "\n", A, "", append=T, sep="\n")

input <- as.numeric(sub(".* nt in (.*) seqs, min.*", "\\1", A[grep(" nt in .* seqs, min", A)]))
OTUs <- as.numeric(sub("Clusters: (.*) Size .*", "\\1", A[grep("^Clusters: ", A)]))

temp <- paste("Clustered ", input, " sequences into ", OTUs, " OTU clusters using a treshold of ", otu_radius_pct, "%.\n\nNext step: Denovo chimera removal.\n")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

# vsearch chimera removal
cmd <- paste(" -uchime_denovo ", folder, "/_data/2_OTU_clustering/", sub("OTUs.fasta", "OTUs+chimeras.fasta", OTU_file), " -nonchimeras ",  folder, "/_data/2_OTU_clustering/", OTU_file, sep="", " -sizein -sizeout -fasta_width 0")
# -fasta_width 0 -> no wrapping
A <- system2(exe, cmd, stdout=T, stderr=T) # remove chimeras!

cat(file=paste(folder, "/_stats/2_OTU_clustering_log.txt", sep=""), "\n", paste("vsearch", cmd), "\n", A, "", append=T, sep="\n")


chimeras <- as.numeric(sub("^Found (.*) .* chimeras, .*", "\\1", A[grep("^Found 13902 .* chimeras, .*", A)]))
OTUs <- as.numeric(sub(".* chimeras, (.*) .* non-chimeras.*", "\\1", A[grep("^Found 13902 .* chimeras, .*", A)]))

temp <- paste("Chimeras removed de novo: ", chimeras, " (", round(chimeras/(chimeras+OTUs)*100, 2), "%)\n          OTUs remaining: ", OTUs , " (", 100-round(chimeras/(chimeras+OTUs)*100, 2), "%)", sep="", "\n")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

}

# compare against refernce sequences (including singletons)

dir.create(paste(folder, "/_data/3_Compare_OTU_derep/", sep=""))
dir.create(paste(folder, "/_stats/3_Compare_OTU_derep/", sep=""))

blast_names <- sub("1_derep_inc_singletons", "3_Compare_OTU_derep", new_names)

if(mapp_singletons){
blast_names <- sub("1_derep_inc_singletons", "3_Compare_OTU_derep", new_names)
}else{
blast_names <- sub(paste("1_derep_minsize_", minuniquesize, sep=""), "3_Compare_OTU_derep", new_names)
}

blast_names <- sub("1_derep_unoise3", "3_Compare_OTU_derep", blast_names)
blast_names <- sub("_PE_derep.*.fasta", ".txt", blast_names)
log_names <- sub("_data", "_stats", blast_names)

if(exe=="usearch"){
cmd <- paste("-usearch_global ", new_names, " -db ", "\"", folder, "/_data/2_OTU_clustering/", OTU_file, "\"", " -strand plus -id ", (100-otu_radius_pct)/100, " -blast6out \"", blast_names, "\" -maxhits 1 -maxaccepts ", maxaccepts, " -maxrejects ", maxrejects, sep="")
} else {
# vserch
cmd <- paste("-usearch_global ", new_names, " -db ", "\"", folder, "/_data/2_OTU_clustering/", OTU_file, "\"", " -strand plus -id ", (100-otu_radius_pct)/100, " -blast6out \"", blast_names, "\" -maxhits 1 -maxaccepts ", maxaccepts, " -maxrejects ", maxrejects, sep="", if(!is.na(threads)){paste(" -threads ", threads, sep="")})
}


files_to_delete <- c(files_to_delete, blast_names)

temp <- paste("Comparing ", length(cmd)," files with dereplicated reads (incl. singletons) against OTUs \"", folder, "/", OTU_file, "\" using \"usearch_global\" with id=", (100-otu_radius_pct)/100," as implemented in ", if(exe=="usearch"){"usearch"}else{"vsearch"}, ".\n", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

exp <- NULL
temp <- new_names
for (i in 1:length(cmd)){
A <- system2(exe, cmd[i], stdout=T, stderr=T)
cat(file= log_names[i], paste("usearch ", cmd[i], sep=""), "\n", A, append=F, sep="\n")

meep <- sub("_data/.*/(.*)", "\\1", temp[i])

if(exe=="usearch"){
pass <- sub(".*, (.*)% matched\r", "\\1", A[grep("matched\r", A)])
} else {
pass <- sub(".*, (.*)% matched\r", "\\1", A[grep("matched\r", A)])
}
exp <- rbind(exp, c(meep, pass))
glumanda <- paste(meep," - ", pass, "% reads matched", sep="")
cat(file="log.txt", glumanda, append=T, sep="\n")
message(glumanda)
}


# check for empty files
empty <- !file.info(blast_names)$size>0

if(sum(empty)>0){
temp <- paste("\nWARNING! ", sum(empty), " file", if(sum(empty>1)){"s do"} else {" does"}, " not contain any sequences matching the OTU sequences!", sep="", "\n", paste(files[empty], collapse="\n"))
message(temp)
cat(file="log.txt", temp , append=T, sep="\n")
}

excluded2 <- c(excluded, blast_names[empty])



# Write raw data OTU table! incl OTU sequences
files <- blast_names[!empty]

tab <- c("NULL")
tab <- as.data.frame(tab, stringsAsFactors=F)
names(tab) <- "ID"

for (i in 1:length(files)){
data <- read.csv(files[i], sep="\t", header=F, stringsAsFactors=F)

names(data) <- c("query", "otu", "ident", "length", "mism", "gap", "qstart", "qend", "target_s", "target_e", "e.value", "bitscore")

data <- data[,c(-11,-12)]

data <- cbind(data, "abund"=as.numeric(sub(".*size=(.*)", "\\1", data$query)), "otu_no"=sub("(.*)size.*", "\\1", data$otu), stringsAsFactors=F)

head(data)

temp <- aggregate(data$abund, by=list(data$otu_no), FUN="sum")
tab <- merge(tab , temp, by.x="ID", by.y="Group.1", all=T, sort=T)
names(tab)[i+1] <- sub(".*derep/(.*).txt", "\\1", files[i])
}



tab <- tab[-1,] # remove NULL entry in the beginning
tab[is.na(tab)] <- 0

mrew <- tab$ID
mrew <- as.numeric(gsub("OTU_(.*)", "\\1", mrew))
tab <- tab[order(as.numeric(mrew)),]



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



# add sequences!
sequ <- read.fasta(paste(folder, "/_data/2_OTU_clustering/", OTU_file, sep=""), forceDNAtolower=F, as.string=T)

# check for different numbers
temp <- sort(table(unlist(c(names(sequ), tab$ID))), decreasing=T)
temp <- temp[temp==1]
if(length(temp)>0){
report <- paste("\n\nWARNING: The following OTU", if(length(temp)>1){"s"}, " did not get any sequences assigned in read mapping. This might happen if they are very low abundant and closely related OTUs are near by. Still maybe take a look.\nAffected sequences: ", paste(names(temp), collapse=" "), sep="")
message(report)
cat(file="log.txt", report, append=T, sep="\n")
}


tab2 <- data.frame("sort"=as.numeric(sub("OTU_", "", tab[,1])), tab, "sequ"=unlist(sequ[match(tab$ID, names(sequ))]), stringsAsFactors=F)


write.csv(file=paste(folder, "/3_Raw_OTU_table.csv", sep=""), tab2, row.names=F)



temp <- "\n\nOTU table generated (including OTU sequences): 3_Raw_OTU_table.csv"
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")


# make raw fasta file with size info
Rsums <- rowSums(tab2[-c(1,2, ncol(tab2))])
names_fasta <- paste(">", tab2$ID, ";size=", Rsums, sep="")
cat(file=paste(folder, "/3_Raw_OTU.fasta", sep=""), sep="\n", paste(names_fasta, tab2$sequ, sep="\n"))

temp <- "OTU files are written also as fasta file in 3_Raw_OTU.fasta"
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")



if(length(excluded)>0){
exp <- rbind(exp, cbind(cbind(excluded), "empty fasta file"))
exp[,1] <- sub(".*/(.*)", "\\1", exp[,1])
exp[,1] <- sub("(.*).txt", "\\1", exp[,1])
exp[,1] <- sub("(.*)_PE_.*", "\\1", exp[,1])

exp <- exp[order(exp[,1]),]
}

exp2 <- data.frame("ID"=exp[,1], "Abundance"=colSums(tab[-1]), "pct_pass"=exp[,2], row.names=1:length(exp[,1]))


write.csv(exp2, file=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_3_pct_matched.csv", sep=""))

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

temp2[is.na(temp2)] <- 0

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
exp$sort[nrow(exp)] <- exp$sort[nrow(exp)-1]+1




# make folder 
dir.create(paste(folder, "/_data/5_subset/", sep=""))


write.csv(exp, file=paste(folder, "/_data/5_subset/5_OTU_sub_", filter, "_not_rematched.csv", sep=""), row.names=F)

OTU_sub_filename <- paste(folder, "/_data/5_subset/5_OTU_sub_", filter, ".fasta", sep="")
write.fasta(as.list(exp$sequ[-nrow(exp)]), exp$ID[-nrow(exp)], file.out=OTU_sub_filename)
file.copy(OTU_sub_filename, sub("_data/5_subset/", "", OTU_sub_filename))


# remapping of reads against subsetted OTUs!
# compare against refernce sequences (including singletons)

dir.create(paste(folder, "/_data/5_subset/usearch_global", sep=""))
dir.create(paste(folder, "/_stats/5_subset/", sep=""))

blast_names <- sub("_PE_derep.*.fasta", ".txt", new_names)

if(mapp_singletons){
blast_names <- sub("1_derep_inc_singletons", "5_subset/usearch_global", blast_names)
}else{
blast_names <- sub(paste("1_derep_minsize_", minuniquesize, sep=""), "5_subset/usearch_global", blast_names)
}

blast_names <- sub("1_derep_unoise3", "5_subset/usearch_global", blast_names)
log_names <- sub("_data/", "_stats/", blast_names)
log_names <- sub("/usearch_global", "", log_names)


cmd <- paste("-usearch_global ", new_names, " -db ", OTU_sub_filename, " -strand plus -id ", (100-otu_radius_pct)/100, " -blast6out ", blast_names, " -maxhits 1 -maxaccepts ", maxaccepts, " -maxrejects ", maxrejects, sep="")

files_to_delete <- c(files_to_delete, blast_names)

temp <- paste("\n\nRemapping ", length(cmd)," files (incl. singletons) against subsetted OTUs \"", OTU_sub_filename, "\" using \"usearch_global\" using id=",(100-otu_radius_pct)/100,".\n", sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")

exp <- NULL
temp <- new_names
for (i in 1:length(cmd)){
A <- system2(exe, cmd[i], stdout=T, stderr=T)
cat(file= log_names[i], paste(exe, " ", cmd[i], sep=""), "\n", A, append=F, sep="\n")

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



empty <- !file.info(blast_names)$size>0

if(sum(empty)>0){
temp <- paste("\nWARNING! ", sum(empty), " file", if(sum(empty>1)){"s do"} else {" does"}, " not contain any sequences matching the OTU sequences!", sep="", "\n", paste(blast_names[empty], collapse="\n"))
message(temp)
cat(file="log.txt", temp , append=T, sep="\n")
}

excluded2 <- c(excluded, blast_names[empty])




files <- blast_names[!empty]

tab <- c("NULL")
tab <- as.data.frame(tab, stringsAsFactors=F)
names(tab) <- "ID"

for (i in 1:length(files)){
data <- read.csv(files[i], sep="\t", header=F, stringsAsFactors=F)

names(data) <- c("query", "otu", "ident", "length", "mism", "gap", "qstart", "qend", "target_s", "target_e", "e.value", "bitscore")

data <- data[,c(-11,-12)]

data <- cbind(data, "abund"=as.numeric(sub(".*size=(.*)", "\\1", data$query)), "otu_no"=sub("(.*)size.*", "\\1", data$otu), stringsAsFactors=F)

head(data)

temp <- aggregate(data$abund, by=list(data$otu_no), FUN="sum")
tab <- merge(tab , temp, by.x="ID", by.y="Group.1", all=T, sort=T)
names(tab)[i+1] <- sub(".*usearch_global/(.*).txt", "\\1", files[i])
}

head(tab)

tab <- tab[-1,] # remove NULL entry in the beginning
tab[is.na(tab)] <- 0

mrew <- tab$ID
mrew <- as.numeric(gsub("OTU_(.*)", "\\1", mrew))
tab <- tab[order(as.numeric(mrew)),]



# add sequences!

sequ <- read.fasta(OTU_sub_filename, forceDNAtolower=F, as.string=T)

tab2 <- data.frame("sort"=as.numeric(sub("OTU_", "", tab[,1])), tab, "sequ"=unlist(sequ), stringsAsFactors=F)

names(tab2) <- sub(".*_data/5_subset/usearch_global/(.*).txt", "\\1", names(tab2)) # SUBSET HERE



# add below treshold OTU counts!

subSums <- read.csv(paste(folder, "/3_Raw_OTU_table.csv", sep=""))
sums <- colSums(subSums[,-c(1,2, ncol(subSums))])
subSums <- as.vector(c(subSums[nrow(subSums)-1, 1]+1, paste("below_", filter, sep=""), sums[sums>0]-colSums(tab2[-c(1,2, ncol(tab2))]), NA))


tab3 <- rbind(tab2, subSums)


# add samples with no hits!
# add excluded samples
if(length(excluded2)>0){
excluded2 <- sub(".*/(.*)", "\\1", excluded2)
excluded2 <- sub("(.*).txt", "\\1", excluded2)
excluded2 <- sub("(.*)_PE_.*", "\\1", excluded2)

if(length(excluded2)>0){
zero <- NULL
for(j in 1:length(excluded2)){
zero <- cbind(zero, rep(0, nrow(tab3)))
}
}

# add zero to table
zero <- data.frame(zero)
names(zero) <- excluded2



# order table by file names!
tab4 <- cbind(tab3[,-ncol(tab3)], zero)
tab4 <- tab4[,c(1, 2, order(names(tab4)[c(-1, -2)])+2)]
tab4 <- data.frame(tab4, "sequ"=tab3[,ncol(tab3)])
} else {tab4 <- tab3}

#head(tab4)

# recalculate sample abundance after sorting table
sampleabundance2 <- NULL
for (f in 3:(ncol(tab4)-1)){
sampleabundance2[f-2] <- sum(as.numeric(tab4[,f]))
}




# set values 2 zero

expZERO <- tab4
#d <- 8
for (d in 3:(ncol(expZERO)-1)){

ZERO <- as.numeric(expZERO[,d])/sampleabundance2[d-2]*100>=filter
discarded <- sum(as.numeric(expZERO[-nrow(expZERO),d][!ZERO[-nrow(expZERO)]]))
expZERO[-nrow(expZERO),d][!ZERO[-nrow(expZERO)]] <- 0 # set zero

expZERO[nrow(expZERO),d] <- as.numeric(expZERO[nrow(expZERO),d])+ discarded # add counts to discarded

}

expZERO[,-ncol(expZERO)][is.na(expZERO[,-ncol(expZERO)])] <- 0




#relative abundance table
expZEROrel <- expZERO
#d <- 8
for (d in 3:(ncol(expZERO)-1)){

expZEROrel[,d] <- as.numeric(expZERO[,d])/sampleabundance2[d-2]*100


}
expZEROrel[,-ncol(expZEROrel)][is.na(expZEROrel[,-ncol(expZEROrel)])] <- 0



write.csv(file=paste(folder, "/5_OTU_table_", filter,".csv", sep=""), tab4, row.names=F)
write.csv(file=paste(folder, "/5_OTU_table_", filter,"_ZERO.csv", sep=""), expZERO, row.names=F)
write.csv(file=paste(folder, "/5_OTU_table_", filter,"_ZERO_rel.csv", sep=""), expZEROrel, row.names=F)

temp <- paste("\n\nSubsetted OTU table generated (", filter, "% abundance in at least ", filterN," sample): ", sub(paste(folder, "/_data/5_subset/", sep=""), "", OTU_sub_filename), sep="")
message(temp)
cat(file="log.txt", temp, append=T, sep="\n")
#### end subsetted OTU table

} # end subsetting



# make plots

RAW <- read.csv(paste(folder, "/3_Raw_OTU_table.csv", sep=""), stringsAsFactors=F)

if(is.na(filter)){  # no filtering, plot only raw data
if(heatmap){
OTU_heatmap(paste(folder, "/3_Raw_OTU_table.csv", sep=""), out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_OTU_plot_raw.pdf", sep=""), abundance=T, col=rev(c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")))
}
}

if(!is.na(filter)){ # only highlight when subsetting
KEEP <- read.csv(paste(folder, "/5_OTU_table_", filter,".csv", sep=""), stringsAsFactors=F)

higlight <- rev(!RAW$ID %in% KEEP$ID)

if(heatmap){
pdf(paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_OTU_plot_3_RAW.pdf", sep=""), height=(nrow(RAW)+20)/10, width=(ncol(RAW)-1)/2)

OTU_heatmap(paste(folder, "/3_Raw_OTU_table.csv", sep=""), abundance=T, col=rev(c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")))


pos <- ncol(RAW) - 2.5
for (i in 1:length(higlight)){

if(higlight[i]){
rect(pos, i-0.5, pos+1, i+0.5, col="Lightgray", border=F)
text(pos+0.1, i, rev(RAW$ID)[i], adj=0, cex=0.5)
}

}
dev.off()
}

if(heatmap){
OTU_heatmap(paste(folder, "/5_OTU_table_", filter,".csv", sep=""), out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_OTU_plot_5_", filter, ".pdf", sep=""), abundance=T, col=rev(c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")))
OTU_heatmap(paste(folder, "/5_OTU_table_", filter,"_ZERO.csv", sep=""), out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_OTU_plot_5_", filter, "_ZERO.pdf", sep=""), abundance=T, col=rev(c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")))
OTU_heatmap(file=paste(folder, "/5_OTU_table_", filter,"_ZERO_rel.csv", sep=""), out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_OTU_plot_5_", filter, "_ZERO_rel.pdf", sep=""), abundance=T, rel=T, col=rev(c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")))
} else { # only plot unfiltered data if zero filtering is applied!
OTU_heatmap(paste(folder, "/3_Raw_OTU_table.csv", sep=""), out=paste(folder, "/_stats/", sub("(.)_.*", "\\1", folder), "_3_Raw_OTU_table.pdf", sep=""), abundance=T, col=rev(c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")))
}
}
if(!heatmap){message("Heatmap generation skipped!")}









cat(file=paste(folder, "/robots.txt", sep=""), "\n# DELETE_START", files_to_delete, "# DELETE_END", append=T, sep="\n")

temp <- "\nModule completed!"
message(temp)
cat(file="log.txt", paste(Sys.time(), "*** Module completed!", "", sep="\n"), append=T, sep="\n")
}
