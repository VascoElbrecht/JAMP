# Core v0.1

Core <- function(files="last", module=NA, letter="A"){

if(files=="last"){ # look in log what last action was to load files from data folder!

}

ABC <- unlist(strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ", ""))

if(file.exists("log.txt")){
log <- readLines("log.txt") # load log
step <- which(log=="PROCESSING MODULE:") # detect processing modules
} else {step <- NULL}

# ADD routine to remove deleted folders from log!


new_folder <- paste(ABC[length(step)+1], ") ", module, sep="")
dir.create(new_folder)
dir.create(paste(new_folder, "_data", sep="/"))
dir.create(paste(new_folder, "_stats", sep="/"))
setwd(new_folder)

# write mudule name in log
temp <- paste("##########", Sys.time(), "PROCESSING MODULE:", new_folder, sep="\n")

cat(file="../log.txt", temp, append=T, sep="\n")





warning("Hello World")


setwd("../") # remove and implement into modules!!!

}







