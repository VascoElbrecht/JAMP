# 161111 - JAMP tutorial!
setwd("~/Documents/UNI_und_VORLESUNGEN/GitHub/JAMP/") # set the path to the PrimerMinder folder you just downloaded

# install the PrimerMiner package icl dependencies
install.packages("JAMP", repos = NULL, type="source", dependencies=T)



# load the package into R
library("JAMP")

# base directory
setwd("~/Desktop/package_test2/")

Demultiplexing_shifted("../Ak15_1.fastq", "../Ak15_2.fastq", tags="BF_BR", combinations="../AK_demulti_used.csv")

U_merge_PE() # merge PE



meep <- Count_sequences(list.files("A_Demultiplexing_shifted/_data", full.names=T))

getwd()
setwd("B_U_merge_PE")


