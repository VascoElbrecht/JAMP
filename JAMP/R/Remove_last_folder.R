# Core v0.2

Remove_last_folder <- function(DoAsk=T){

if(file.exists("log.txt")){
log <- readLines("log.txt") # load log
step <- which(log=="PROCESSING MODULE:") # detect processing modules
} else {stop()
}

if(DoAsk){
message(paste("WARNING: Are you SURE you want to remove the last module \"", log[step[length(step)]+1], "\" and delete everything in folder and it's entries in the log file?\nThis action can not be undone!\n", getwd(), "/", log[step[length(step)]+1], "\n", sep=""))
Sys.sleep(0.1)

temp <- menu(c("y", "n"), title="delete the prvious data? y = Yes, n = No")
} else {temp <- 1}

if(temp==1){

message(paste("Deleteing the folder ", log[step[length(step)]+1], "!", sep=""))
system2("rm", paste("-fr", log[step[length(step)]+1]))
if(length(step)!=1){
cat(file="log.txt", log[1:(step[length(step)]-3)], append=F, sep="\n")
} else {
jup <- file.remove("log.txt")
}


} else {
message("Stopped! Did not delete anything!")
}

}







