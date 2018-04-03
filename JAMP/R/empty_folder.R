# U_max_ee v0.1

Empty_folder <- function(){

folder <- Core(module="Empty_Folder")
cat(file="log.txt", c("\n","Version v0.1", "\n"), append=T, sep="\n")

temp <- paste("Folder \"", folder, "\" generated!\nPlease place files for further processing in \"", folder, "/_data\"", sep="", "\n\nPlease have a look at the JAMP wiki on github for assistance \nhttps://github.com/VascoElbrecht/JAMP/wiki\n")

message(temp)

cat(file="log.txt",temp, append=T, sep="")


message(" ")
message("Module completed!")

cat(file="log.txt", paste(Sys.time(), "*** Module completed!", "", sep="\n"), append=T, sep="\n")


}

