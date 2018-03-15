# 180125 convert table into files for demultiplexing

# set directory
setwd("~/Documents/UNI_und_VORLESUNGEN/14 Guelph/1 TEACHING/2018 metabarcoding course/1 JAMP/_converter")

#table
data <- read.csv("Leray2017.txt", stringsAsFactor=F, sep="\t")

data <- data[1:7,]

ind <- substr(data$Primer.sequence, 1, 6)


indexe <- data.frame("barcode"=ind, "rm"=nchar(ind), "ID"=paste("i", 1:7, sep=""), "comment"=data$Primer.sequence)
indexe

write.csv(indexe, row.names=F, file="indexe_1.csv")

combos <- data.frame("ID"=paste(indexe$ID, "_", indexe$ID, sep=""), "File1"=paste("L", 1:7, "_r1.txt", sep=""), "File2"=paste("L", 1:7, "_r2.txt", sep=""))

combos

write.csv(combos, row.names=F, file="combos_1.csv")

