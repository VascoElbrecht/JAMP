# delete_data

delete_data <- function(ask=T){

cat(file="log_deleted_files.txt", "Running module \"delete_data\"\n", append=T, sep="\n")


# detect folders from log file
if(file.exists("log.txt")){
log <- readLines("log.txt") # load log
step <- which(log=="PROCESSING MODULE:")+1 # detect 
folder <- log[step]
} else { # 
stop("No log file found, module stopped!")
}

# detecting files to delete!

files_to_delete <- NULL
number_of_files_per_folder <- NULL

for (i in 1:length(folder)){
data <- readLines(paste(folder[i], "/robots.txt", sep=""))

TF <- sum(data=="delete_data=TRUE")
if(TF==0){ # don't delete folder
number_of_files_per_folder[i] <- "FALSE"
} else {

A <- which(data=="# DELETE_START")+1
B <- which(data=="# DELETE_END")-1

if(!A>B){
files_to_delete <- c(files_to_delete, data[A:B])
number_of_files_per_folder[i] <- length(data[A:B])
} else {
number_of_files_per_folder[i] <- "EMPTY!"
}
}
}

# reporting
temp_num <- suppressWarnings(as.numeric(number_of_files_per_folder))


cat(file="log_deleted_files.txt", files_to_delete, append=T, sep="\n")

if(ask){
temp <- menu(c("y", "n"), title=paste("Will delete a total of ", sum(temp_num, na.rm=T), " files in ", length(!is.na(temp_num)), " folders.\nSee table below for details and verify that his is as supposed. You can also take a look in the \"log_deleted_files.txt\" to verify wich files shoudl be deleted!", sep="", "\n\n", paste(number_of_files_per_folder ," - ", folder, collapse="\n"), "\n\nWould you like to permanently DELETE these files? y = Yes, n = No"))
} else {temp=2}

if(temp==2){
message("ABORTED! - No files deleted")
}

# actually delete files!
if(temp==1){
for (i in 1:length(files_to_delete)){

if(!file.exists(files_to_delete[i])){
message(paste(files_to_delete[i], " does NOT exist!", sep=""))

updatefile <- readLines("log_deleted_files.txt")
updatefile[updatefile ==files_to_delete[i]] <- paste(updatefile[updatefile ==files_to_delete[i]], " - FILE DOES NOT EXIST!")
cat(file="log_deleted_files.txt", updatefile, append=F, sep="\n")

} else { # file exists, let's delete it
file.remove(files_to_delete[i])

updatefile <- readLines("log_deleted_files.txt")
updatefile[updatefile ==files_to_delete[i]] <- paste(updatefile[updatefile ==files_to_delete[i]], " - DELETED!")
cat(file="log_deleted_files.txt", updatefile, append=F, sep="\n")
}

}

}

message("Files deleted!")


}

