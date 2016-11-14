# 161111 - JAMP tutorial!
setwd("~/Documents/UNI_und_VORLESUNGEN/GitHub/JAMP/") # set the path to the PrimerMinder folder you just downloaded

# install the PrimerMiner package icl dependencies
install.packages("JAMP", repos = NULL, type="source", dependencies=T)

# load the package into R
library("JAMP")

setwd("~/Desktop/package_test/")
# creating configuration file and batch downloading reads


Demultiplexing_shifted("~/Documents/UNI_und_VORLESUNGEN/11 phd projects/1 Meta ANNA chiro/2 Ak15 DATA/Ak15_1.fastq.gz", "~/Documents/UNI_und_VORLESUNGEN/11 phd projects/1 Meta ANNA chiro/2 Ak15 DATA/Ak15_2.fastq.gz", tags="BF_BR", combinations="../AK_demulti_used.csv")

getwd()





Count_sequences(list.files("A_Demultiplexing_shifted/_data", full.names=T)[1:10])


