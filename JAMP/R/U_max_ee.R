# U_max_ee v0.1

U_max_ee <- function(files="latest", max_ee=0.5){

Core(module="U_max_ee")
cat(file="../log.txt", c("\n","Version v0.1", "\n"), append=T, sep="\n")
message(" ")

if (files=="latest"){
source("robots.txt")
files <- list.files(paste("../", last_data, "/_data", sep=""), full.names=T)
}

temp <- paste("Starting to quality filter (max expected errors = " , max_ee, ") in ", length(files), " samples.", sep="")
message(temp)
message(" ")
cat(file="../log.txt", temp, append=T, sep="\n")

# new file names
new_names <- sub(".*(_data/.*)", "\\1", files)
new_names <- sub(".fastq", "_ee.fasta", new_names)

dir.create("_stats/merge_stats")
log_names <- sub("_data", "_stats/merge_stats", new_names)
log_names <- sub("_ee.fasta", "_ee.txt", log_names)

# cmd max EE
cmd <- paste("-fastq_filter \"", files, "\" -fastaout \"", new_names, "\" -fastq_maxee ", max_ee, sep="")


tab_exp <- NULL
for (i in 1:length(cmd)){
A <- system2("usearch", cmd[i], stdout=T, stderr=T)
cat(A, file=log_names[i], sep="\n")

imput <- A[grep("FASTQ recs", A)]
imput <- sub(" +(.*) +FASTQ.*", "\\1", imput)
imput <- as.numeric(imput)

pass <- A[grep("Converted ", A)]
pass <- sub(" +(.*) +Converted.*", "\\1", pass)
pass <- as.numeric(pass)

pct <- round(pass/imput*100, digits=2)

short_name <- sub("_data/(.*)_PE.*", "\\1", new_names[i])
tab_exp <- rbind(tab_exp, c(short_name, pass, pct))

meep <- paste(short_name, ": ", pct, "% pass EE ", sep="")
message(meep)
cat(file="../log.txt", meep, append=T, sep="\n")
}
cat(file="../log.txt", "\n", append=T, sep="\n")

tab_exp <- data.frame(tab_exp)
names(tab_exp) <- c("Sample", "Sequ_count", "percent_merged")

write.csv(tab_exp, "_stats/max_ee_stats.csv")

# make some plots?

message(" ")
message("Module completed!")

cat(file="../log.txt", paste(Sys.time(), "Module completed!", "", sep="\n"), append=T, sep="\n")

setwd("../")
}

