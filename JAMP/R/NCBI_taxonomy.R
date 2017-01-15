# Bold taxonomy

Ncbi_taxonomy <- function(file=NA, folder=""){

oldwd <- getwd()
setwd(folder)
getwd()
#cat(file="../log.txt", c("\n","Version v0.1", "\n"), append=T, sep="\n")
message("Getting taxonomy via NCBI api.")

data <- read.csv(file, stringsAsFactors=F)
sequ <- which(!is.na(data$sequ)) # sequences exist

dir.create("_NCBI_data", showWarnings = FALSE)
setwd("_NCBI_data")
dir.create("_stats", showWarnings = FALSE)


# query sequences to bold API

# check if data already downloaded

sequ_inital_DL <- sequ

# skipp if already downloaded!
temp <- paste(data$ID[sequ_inital_DL], ".csv", sep="")
sequ_inital_DL <- sequ_inital_DL[!temp%in%list.files("_stats", "csv")]

for (i in sequ_inital_DL){



}






setwd(oldwd)
}
