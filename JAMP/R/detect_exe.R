# detect v0.1

detect_exe <- function(exe){
A <- system2(exe, "-version", stdout=T, stderr=T)


if(grepl("usearch", A)[1]){
return("usearch")
}

if(grepl("vsearch", A)[1]){
return("vsearch")
}
}

