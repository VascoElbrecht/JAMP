# FastQC

FastQC <- function(files="latest", ){

# get folder path!
if(file.exists("log.txt")){
log <- readLines("log.txt") # load log
step <- which(log=="PROCESSING MODULE:")+1 # detect 

folder <- log[step[length(step)]]
message("Runing module: FastQC")

} else { # if nothing has been run jet
message("Runing module: FastQC")
dir.create("_FastQC", showWarnings=F)
folder <- "_FastQC"
}



if(files=="latest"){
files <- list.files(paste(folder, "/_data", sep=""), full.name=T)
}


dir.create(paste(folder, "/FastQC/_raw_reports/", sep=""), showWarnings=F)

# run fastQC module
cmd <- paste("-o ", folder, "/FastQC/_raw_reports/ ", sep="", paste("\"", files, "\" ", collapse="", sep="") )

A <- system2("fastQC", cmd)


system2("vsearch", paste("--fastq_eestats ", files[1], sep="", " --output test.txt"))


# run fastqcr
QCpath <- paste(folder, "/FastQC/_raw_reports/", sep="")
QCfiles <- list.files(QCpath, full.names=T, pattern="fastqc.zip")

QClist <- list()
for(i in 1:length(QCfiles)){
QClist[[i]] <- qc_read(QCfiles[i])
}

A <- qc_read(QCfiles[i])
qc_plot(QClist[[5]], "Per sequence quality scores")


# aggregate




message(" ")
message(" Module completed!")

cat(file="log.txt", paste(Sys.time(), "\n", "*** Module completed!", "", sep="\n"), append=T, sep="\n")


}

