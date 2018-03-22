# Core v0.2

Core <- function(files="last", module=NA, letter="A", delete_data=F){


#numbering for folders
ABC <- unlist(strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ", ""))
ABC <- c(ABC, paste("Z0", 1:9, sep=""), paste("Z", 10:74, sep=""))

if(file.exists("log.txt")){
log <- readLines("log.txt") # load log
step <- which(log=="PROCESSING MODULE:") # detect processing modules
} else {step <- NULL
}

# Ask to delete incomplete folders!

if(!is.null(step)){
completed <- "*** Module completed!" %in% log[step[length(step)]:length(log)]

if(!completed){
message(paste("WARNING: It looks like the previously run module did crash! Do you want to delete the prvious folder:\n", getwd(), "/", log[step[length(step)]+1], "\n\n(recommended, but confirm the pasth is correct!)", sep=""))
Sys.sleep(0.1)
}
if(!completed){
temp <- menu(c("y", "n"), title="delete the prvious data? y = Yes, n = No")


if(temp==1){

message(paste("Deleteing the folder ", log[step[length(step)]+1], "!", sep=""))
system2("rm", paste("-fr", log[step[length(step)]+1]))
if(length(step)!=1){cat(file="log.txt", log[1:(step[length(step)]-3)], append=F, sep="\n")} else {cat(file="log.txt", sep="")}

} else {"Did not delete anything! but please check your imput and log files as JAMP will likely not run correctly!"}
}
}


if(!file.exists("log.txt")){
cat(file="log.txt", sep="")
}

log <- readLines("log.txt") # reload log
step <- which(log=="PROCESSING MODULE:") 

if(!is.null(step)){new_folder <- paste(ABC[length(step)+1], "_", module, sep="")} else {new_folder <- paste(ABC[1], "_", module, sep="")} 
dir.create(new_folder)
dir.create(paste(new_folder, "_data", sep="/"))
dir.create(paste(new_folder, "_stats", sep="/"))


# write mudule name in log
temp <- paste("##########", Sys.time(), "PROCESSING MODULE:", new_folder, "", sep="\n")

cat(file="log.txt", temp, append=T, sep="\n")

if(is.null(step)){temp <- NA} else {temp <- log[step[length(step)]+1]}
cat(file=paste(new_folder, "/robots.txt", sep=""), "delete_data=", delete_data, "\n", "last_data=\"", temp, "\"\n", append=T, sep="")


message(paste("Runing module:", module))
return(new_folder)
}







