# Plot illumina tiles from fastq files
# A, G, C, T, N - color


PlotTiles <- function(file=NA, subset=100000, cycles="ALL", folder="TilePlots", color=c("Red", "Gold", "Blue", "Green2", "Gray"), format="png", width=2, height=19, Nalert=0.10, cex=0.75, exe="vsearch"){



# make folder

if(!file.exists(folder)){
message("Making folder ", folder)
dir.create(folder)
dir.create(paste(folder, "/plots", sep=""))
dir.create(paste(folder, "/stats", sep=""))

} else {
temp <- menu(c("y", "n"), title=paste("Folder \"", folder ,"\" already exists. Things inside might be overwritten. Is that ok?\n\ny = Yes, over write exisiting files!\nn = No, stop function!", sep=""))
if(temp==2){
stop("Fuction stopped, nothing overwritten. Feel free to specify a different folder with \"folder\".")
}
}





# Import file
compression <- gsub(".*\\.(.*)$", "\\1", file)

#if(compression=="gz"){
#} else if(compression=="bz2"){
#} else if(compression=="fastq"){
#} else {
#warning("Raw data has to be compressed in \".gz\" or \".bz2\" or uncompressed in \".fastq\" format! Processing was stopped, please check your raw data format.")
#stop()
#}


#data <- readLines(con, subset*4) # read chunck
#close(con)

# apply subsampling using vsearch

suppressWarnings(A <- system2(exe, paste("-fastx_subsample \"", file, "\" -fastqout \"", folder, "/subset_", format(subset, scientific=F), ".txt\" -sample_size ", format(subset, scientific=F), sep=""), stdout=T, stderr=T))


if(length(grep("Fatal error: Cannot subsample more reads than in the original sample", A))>0){
count <- sub(".* nt in (.*) seqs, min .*", "\\1", A[grep("seqs, min", A)])
message("File ", file, "contains only ", count, " sequences, thus no subsampling (subset=", format(subset, scientific=F), ") will be applyed!\nReading in file...")

data <- readLines(file)

} else {
message(paste(A, collapse="\n"))
message("subsetting was applyed!")
data <- readLines(paste(folder, "/subset_", format(subset, scientific=F), ".txt", sep=""))
}


# get cycles to plot

if(cycles=="ALL"){
message("Automatically detecting the number of cycles in file!")
cycles <- 1:max(nchar(data[seq(2, length(data), 4)]))

message("Will generate plots for possition 1 to ", max(cycles), ". Specify the numbers of the cycles in \"cycle\" if you would like to plot a specific cycle.")
} else {
message("Plotting user defined possitions: ", paste(cycles, collapse=" "))
}


#data <- c("@FS10000474:11:BPG80716-2107:1:1101:1930:1000", "@M00307:138:000000000-CTLD2:1:1101:10100:1464 2:N:0:1")

data[seq(1, length(data), 4)] <- sub(" .:.:.:.", "", data[seq(1, length(data), 4)])

sub(".+:.+:.+:.+:(.+):.+:.+ ?.?:?.?:?.?:?.?", "\\1", data)
sub(".+:.+:.+:.+:(.+):.+:.+ .:.:.:.", "\\1", data)

# extract relevant information from fastq

data2 <- data.frame("tile"=as.numeric(sub(".*:(.*):.*:.*", "\\1", data[seq(1, length(data), 4)])), "x"=as.numeric(sub(".*:.*:(.*):.*", "\\1", data[seq(1, length(data), 4)])), "y"=as.numeric(sub(".*:.*:.*:(.*)", "\\1", data[seq(1, length(data), 4)])),  stringsAsFactors=F, "sequ"=data[seq(2, length(data), 4)])

unique <- unique(data2$tile)
unique <- sort(unique)

message("A total of ", length(unique), " tiles where detected!" )
message("Tiles: ", paste(unique, collapse=", "))

# make plots
savename <- sub("[.*/]?(.*)\\..*", "\\1", file)


cat("Cycle\tTile\tSpots\tA\tT\tG\tC\tN\tIssues", file=paste(folder, "/warnings.tsv", sep=""), sep="\n", append=F)



for (k in cycles){

message("processing possition: ", k)
# Extract base composition at specified cycle!
extacted <- substr(data2$sequ, k, k)

stats <- data.frame("Cycle"=NA, "Tile"=NA, "Spots"=NA, "A"=NA, "T"=NA, "G"=NA, "C"=NA, "N"=NA, "Issues"=NA)



if(format=="png"){
png(paste(folder, "/plots/", savename, "_", k, ".png", sep=""), width=width*440, height=height*400)
}
if(format=="pdf"){
pdf(paste(folder, "/plots/", savename, "_", k, ".pdf", sep=""), width=width*4.4, height=height*4, useDingbats=F)
}

par(mfcol=c(height, width))

for(i in 1:length(unique)){

# select tile
temp <- extacted[data2$tile==unique[i]]

if(i==1){
stats[1,] <- c(k, i, length(temp[temp!=""]), sum(temp=="A"), sum(temp=="T"), sum(temp=="G"), sum(temp=="C"), sum(temp=="N"), NA)
} else {
stats <- rbind(stats, c(k, i, length(temp[temp!=""]), sum(temp=="A"), sum(temp=="T"), sum(temp=="G"), sum(temp=="C"), sum(temp=="N"), NA))
}


if(stats[i,3]==stats[i,8]){
stats[i,9] <- "No Call"
message("WARNING: No Call in cycle ", k, " tile ", i, "!")
cat(paste(stats[i,], collapse="\t"), file=paste(folder, "/warnings.tsv", sep=""), sep="\n", append=T)
} else if(stats[i,8]/stats[i,3]>=Nalert){
stats[i,9] <- paste("Above ", Nalert*100, "% Ns!")
message("WARNING: Tile ", i, " in cycle ", k, " contains ", round(stats[i,8]/stats[i,3]*100), "% Ns!")
cat(paste(stats[i,], collapse="\t"), file=paste(folder, "/warnings.tsv", sep=""), sep="\n", append=T)
}


# add colors
temp[temp=="A"] <- color[1]
temp[temp=="G"] <- color[2]
temp[temp=="C"] <- color[3]
temp[temp=="T"] <- color[4]
temp[temp=="N"] <- color[5]

par(mar=c(0,0,0,0))
# convert short sequences
x <- data2$x[data2$tile==unique[i]][temp!=""]
y <- data2$y[data2$tile==unique[i]][temp!=""]
temp <- temp[temp!=""]

plot(x, y, main="", col= temp, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", pch=20, cex=cex)
text(min(x), max(y), paste("Tile: ", unique[i], sep="", " - clusters: ", length(x), pos=4))



sum(temp=="Gray")/length(temp)

}

dev.off()

# write stats file

write.table(stats, file=paste(folder, "/stats/", savename, "_", k, ".tsv", sep=""), sep="\t", row.names=F)

# Add way to detect mixed tiles




}






}

