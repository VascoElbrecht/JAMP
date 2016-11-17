# 161111 - JAMP tutorial!
setwd("~/Documents/UNI_und_VORLESUNGEN/GitHub/JAMP/") # set the path to the PrimerMinder folder you just downloaded

#install.packages("seqinr", dependencies=T)

# install the PrimerMiner package icl dependencies
install.packages("JAMP", repos = NULL, type="source")



# load the package into R
library("JAMP")

# base directory
setwd("~/Desktop/package_test2/")

Demultiplexing_shifted("../Ak15_1.fastq", "../Ak15_2.fastq", tags="BF_BR", combinations="../AK_demulti_used.csv")

U_merge_PE() # merge PE

# select sequ to make revcomp of!
revcomp <- list.files("B_U_merge_PE/_data")
revcomp_tf <- grepl(".*_.*_.R.*_.*_.*", revcomp)
cbind(revcomp, revcomp_tf)

U_revcomp(RC= revcomp_tf) # make RevComp of selected reads

Cutadapt(forward="BF2", reverse="BR1")


setwd("D_Cutadapt")


