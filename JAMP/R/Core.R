# Core v0.1

Core <- function(files="last", module=NA, letter="A", delete_data=F){

if(files=="last"){ # look in log what last action was to load files from data folder!
}

#numbering for folders
ABC <- unlist(strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ", ""))
ABC <- c(ABC, paste("Z0", 1:9, sep=""), paste("Z", 10:74, sep=""))

if(file.exists("log.txt")){
log <- readLines("log.txt") # load log
step <- which(log=="PROCESSING MODULE:") # detect processing modules
} else {step <- NULL}

# ADD routine to remove deleted folders from log!


new_folder <- paste(ABC[length(step)+1], "_", module, sep="")
dir.create(new_folder)
dir.create(paste(new_folder, "_data", sep="/"))
dir.create(paste(new_folder, "_stats", sep="/"))
setwd(new_folder)

# write mudule name in log
temp <- paste("##########", Sys.time(), "PROCESSING MODULE:", new_folder, "", sep="\n")

cat(file="../log.txt", temp, append=T, sep="")

if(is.null(step)){temp <- NA} else {temp <- log[step[length(step)]+1]}
cat(file="robots.txt", "delete_data=", delete_data, "\n", "last_data=\"", temp, "\"\n", append=T, sep="")






message(paste("Runing module:", module))

 
}







