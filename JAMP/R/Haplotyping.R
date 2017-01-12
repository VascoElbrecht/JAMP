# Haplotyping v0.1

haplotyping <- function(folder=NA, replicate_one=NA, replicate_two=NA, new_name=NA, min_OTU=0.1, min_haplo=1){


#cat(file="../log.txt", c("\n","Version v0.1", "\n"), append=T, sep="\n")
message("making haplo / OTU mixes")


setwd(folder)
dir.create("_haplotypes", showWarnings = FALSE)
setwd("_haplotypes")
dir.create("_fasta_OTU_files", showWarnings = FALSE)

denoised <- list.files("../_data/1_derep_unoise2", full.names=T)
mapping <- list.files("../_data/5_subset/usearch_global", full.names=T)


for(i in 1:length(mapping)){

DNA <- read.fasta(denoised[i], as.string=T, forceDNAtolower=F)
DNA <- data.frame("name"=names(DNA), "sequ"=as.vector(unlist(DNA)), stringsAsFactors=F)
OTUtab <- read.csv(mapping[i], sep="\t", stringsAsFactors=F, header=F)

folder_sample <- sub("../_data/1_derep_unoise2/(.*)_PE.*", "\\1", denoised[i])
dir.create(paste("_fasta_OTU_files/", folder_sample, sep=""))
matchedOTUs <- unique(OTUtab$V2)


for(k in 1:length(matchedOTUs)){

temp <- OTUtab[OTUtab$V2== matchedOTUs[k],]

cat(paste(DNA$name[DNA$name%in%temp$V1], "\n", DNA$sequ[DNA$name%in%temp$V1], collapse="\n", sep=""), file=paste("_fasta_OTU_files/", folder_sample, "/", matchedOTUs[k], ".txt", sep=""))
}

}


}
