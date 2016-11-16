# U_merge_PE v0.1
# maybe add option to split and merge large files automatically?

U_revcomp <- function(files="latest", RC=c(T,F), fastq=T, copy_unchanged=T){

Core(module="U_revcomp")
cat(file="../log.txt", c("Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files=="latest"){
source("robots.txt")
files <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T)
}

# new file names
new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub("_PE.", "_PE_RC.", new_names)

# copy over files which are not rev comp
if(copy_unchanged){
#copy_name <- sub(".*/(_data/.*)", "\\1", files)
temp <- file.copy(files[!RC], new_names[!RC])

meep <- c(paste(length(temp), " files where copied to ", sub(".*/(.*)", "\\1", getwd()), " without generating the RevComp:", sep=""), sub(".*_data/", "", files[!RC]))
for(i in 1:length(meep)){
cat(file="../log.txt", meep[i] , append=T, sep="\n")
message(meep[i])
}
cat(file="../log.txt", "" , append=T, sep="\n")
}
message(" ")


cmd <- paste("-fastx_revcomp \"", files[RC], if(fastq){"\" -fastqout \""} else {" -fastaout \""}, new_names[RC], "\" -label_suffix _RC",  sep="")

temp <- paste("RevComp is generated for the following ", length(files[RC]), " files:", sep="")
cat(file="../log.txt", temp , append=T, sep="\n")
message(temp)


temp <- new_names[RC]
for (i in 1:length(cmd)){
system2("usearch", cmd[i], stdout=T, stderr=T)
meep <- sub(".*_data/(.*)", "\\1", temp[i])
cat(file="../log.txt", meep, append=T, sep="\n")
message(meep)
}




message(" ")
message(" Module completed!")

cat(file="../log.txt", paste(Sys.time(), "\n", "Module completed!", "", sep="\n"), append=T, sep="\n")

setwd("../")
}

