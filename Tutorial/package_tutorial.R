# 170201 jamp pipeline, short tutorial
#
# download the raw data from and place it into the tutorial folder (Study we are looking at https://peerj.com/articles/3006/) 
# MiSeq Run1, R1 direction:
# https://dx.doi.org/10.6084/m9.figshare.4039821.v1
# MiSeq Run1, R2 direction:
# https://dx.doi.org/10.6084/m9.figshare.4039860.v1
# 
# Installing dependencies needed fro JAMP
install.packages(c("bold", "XML", "seqinr", "devtools"), dependencies=T)
# Load devtools and install package directly from GitHub
library("devtools")
install_github("VascoElbrecht/JAMP", subdir="JAMP")



setwd("~/Desktop/JAMP_pipeline")
list.files()

library("JAMP")

Demultiplexing_shifted(file1="../16_S10_L001_R1_001_run1.fastq", file2="../16_S10_L001_R2_001_run1.fastq", tags="../_converter/indexe_1.csv", combinations="../_converter/combos_1.csv")

# Paired end merging
U_merge_PE(fastq_pctid=75)

# check for PhiX
system2("usearch", "-usearch_global A_Demultiplexing_shifted/_data/N_debris_r1.txt -db PhiX.fasta -id 0.9 -strand both -blast6out PhiX_table.txt")



# trimm primers
Cutadapt(forward="GGWACWGGWTGAACWGTWTAYCCYCC", # mlCOIintF
reverse="TANACYTCNGGRTGNCCRAARAAYCA") # jgHCO, I (inosin) replaced with N

# trim reads in different orrientation
PEmerged <- list.files("~/Desktop/JAMP_pipeline/B_U_merge_PE/_data", full.names=T)

Cutadapt(files= PEmerged, forward="TANACYTCNGGRTGNCCRAARAAYCA", # jgHCO, I (inosin) replaced with N
reverse="GGWACWGGWTGAACWGTWTAYCCYCC") # mlCOIintF


U_revcomp(RC=T)


# merge forward and now reverse complement reverse reads!

dir.create("F_merge/_data", recursive=T)

FW <- list.files("C_Cutadapt/_data", full.names=T)
RC <- list.files("E_U_revcomp/_data", full.names=T)

i <- 1

for (i in 1:length(FW)){
system2("cat", paste(FW[i], " ", RC[i], " > F_merge/_data/", sub("E_U_revcomp/_data/(.*_).*", "\\1", RC[i]), "merged.txt", sep=""))
}

cat(file="log.txt", append=T, c("\nPROCESSING MODULE:", "F_merge"), sep="\n")


# discard with non target length
Minmax(min=(313-10), max=(313+10))

# discard reads above 1 expected error
U_max_ee(max_ee=1)

# subsample to lowest sample size, should be done if samples are widely different in sequencing depth (as one starts with)
U_subset(sample_size=60000)


#cluster OTUs
U_cluster_otus(filter=0.01)
file.rename("J_U_cluster_otus", "J_U_cluster_otus - 60k")

#cluster OTUs (without subsetting)
no_subset <- list.files("~/Desktop/JAMP_pipeline/H_U_max_ee/_data", full.names=T)

U_cluster_otus(files= no_subset, filter=0.01)


# assign taxonomy to OTUs without sub setting! K_U_cluster_otus
Bold_web_hack(file="K_bold_results.txt")


# setting read counts below 0.01% to 0
data <- read.csv("K_U_cluster_otus/5_OTU_table_0.01.csv", stringsAsFactors=F)
colSum <- colSums(data[,3:9])
i <- 3
for(i in 3:9){
keep <- data[-nrow(data),i]/colSum[i-2]*100>=0.01 # find values blow 0.01
data[nrow(data),i] <- data[nrow(data),i]+sum(data[,i][!keep]) # below 0.01
data[,i][!keep] <- 0 # set lavues to 0
}
write.csv(data, "K_U_cluster_otus/data_cleaned.csv", row.names=F)

data2 <- data
#relative abundance
for(i in 3:9){
data2[,i] <- round(data2[,i]/colSum[i-2]*100, 4) #converte data into relative abundance
}
write.csv(data2, "K_U_cluster_otus/data_cleaned_rel.csv", row.names=F)




# haplotyping
# from merged data:
no_subset <- list.files("~/Desktop/JAMP_pipeline/G_Minmax/_data", full.names=T)

#Keep only sequences of 313 bp length
Minmax(file=no_subset, min=313, max=313)

# stricter EE filtering
U_max_ee(max_ee=0.2)

#rename files to fasta files (Will fix this in the future)
file.rename(list.files("M_U_max_ee/_data", full.names=T), sub("txt", "fasta", list.files("M_U_max_ee/_data",  full.names=T)))

Denoise(minsize=5, minrelsize=0.001, OTUmin=0.1, minHaploPresence=2)


