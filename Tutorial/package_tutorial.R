# 12052021 jamp pipeline, short tutorial
#

# 
# Installing dependencies needed fro JAMP
install.packages(c("bold", "XML", "seqinr", "devtools", "fastqcr"), dependencies=T)
# Load devtools and install package directly from GitHub
library("devtools")
install_github("VascoElbrecht/JAMP", subdir="JAMP")


# set workin directory
setwd("~/Desktop/JAMP_pipeline/Tutorial")
# download the raw data from and place it into the tutorial folder (Study we are looking at https://peerj.com/articles/3006/) 
# MiSeq Run1, R1 direction:
# https://dx.doi.org/10.6084/m9.figshare.4039821.v1
# MiSeq Run1, R2 direction:
# https://dx.doi.org/10.6084/m9.figshare.4039860.v1

download.file("https://ndownloader.figshare.com/files/6503952", "16_S10_L001_R1_001_run1.fastq.gz")
download.file("https://ndownloader.figshare.com/files/6503991", "16_S10_L001_R2_001_run1.fastq.gz")
system2("gunzip", "16_S10_L001_R1_001_run1.fastq.gz")
system2("gunzip", "16_S10_L001_R2_001_run1.fastq.gz")



list.files()

library("JAMP")

# JAMP generates a new folder for each processing step! Should your files alreay be demultiplexed or you want to start somwhere else in the pipeline with preprocessed reads, you can generate an empty folder and place your files to be processed in "_data"
Empty_folder()

#To delete the last generated folder, run
Remove_last_folder()



# In this example we are dealing with sequence raw data that is not jet demultiplexed. To demultiplex run:
Demultiplexing_shifted(file1="16_S10_L001_R1_001_run1.fastq", file2="16_S10_L001_R2_001_run1.fastq", tags="_converter/indexe_1.csv", combinations="_converter/combos_1.csv")



# check for PhiX
# only subsample 10000 reads
if(T){
system2("vsearch", "-fastx_subsample A_Demultiplexing_shifted/_data/N_debris_R1.fastq -fastaout A_Demultiplexing_shifted/_data/N_debris_R1.fasta -sample_size 10000")
} else { # check all reads in read 1 for PhiX
system2("paste", " - - - - < A_Demultiplexing_shifted/_data/N_debris_R1.fastq | cut -f 1,2 | sed 's/^@/>/' | tr \"\t\" \"\n\" > A_Demultiplexing_shifted/_data/N_debris_R1.fasta")
}
system2("vsearch", "-usearch_global A_Demultiplexing_shifted/_data/N_debris_R1.fasta -db PhiX.fasta -id 0.9 -strand both -blast6out PhiX_table.txt -maxrejects 1 -maxaccepts 1", stdout=T, stderr=T)


# Paired end merging
Merge_PE()


# trimm primers (mlCOIintF and jgHCO)
Cutadapt(forward="GGWACWGGWTGAACWGTWTAYCCYCC", reverse="TAIACYTCIGGRTGICCRAARAAYCA", bothsides=T)
#by using "bothsides=T", forward or reverse primers are detected on both ends. This is not nessesary for fusion primers.

# discard with non target length
Minmax(min=(313-10), max=(313+10))

# discard reads above 1 expected error
Max_ee(max_ee=1)

# subsample to lowest sample size, should be done if samples are widely different in sequencing depth (as one starts with)
# needs to be umdated, right now still on userch! Vserch not supported right now
U_subset(sample_size=60000)




#cluster OTUs
Cluster_otus(filter=0.01)
file.rename("G_U_cluster_otus", "G_U_cluster_otus - 60k")

#cluster OTUs (without subsetting)
no_subset <- list.files("E_U_max_ee/_data", full.names=T)

Cluster_otus(files= no_subset, filter=0.01)

# assign taxonomy to OTUs without sub setting! K_U_cluster_otus
Bold_web_hack(file="K_bold_results.txt")





# haplotyping
# from merged data:
no_subset <- list.files("D_Minmax/_data", full.names=T)

# Keep only sequences of 313 bp length
Minmax(file=no_subset, min=313, max=313)

# Stricter EE filtering
U_max_ee(max_ee=0.2)

# Extracting haplotypes from the metabarcoding data
Denoise(minsize=5, minrelsize=0.001, OTUmin=0.01, minHaploPresence=1)





