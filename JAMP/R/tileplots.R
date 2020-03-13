# Plot illumina tiles from fastq files
# A, G, C, T, N - color


PlotTiles <- function(color=c("Red", "Yellow", "Blue", "Green", "Gray")){













data <- readLines("~/Desktop/primer_test/JAMP/A_Empty_Folder/_data/_101_4_r2.fastq")


for (k in 1:301){


base <- k

data2 <- data.frame("tile"=as.numeric(sub("@M05074:192:000000000-J2J3L:1:(.*):.*:.* 2:N:0:0", "\\1", data[seq(1, length(data), 4)])), "x"=as.numeric(sub("@M05074:192:000000000-J2J3L:1:.*:(.*):.* 2:N:0:0", "\\1", data[seq(1, length(data), 4)])), "y"=as.numeric(sub("@M05074:192:000000000-J2J3L:1:.*:.*:(.*) 2:N:0:0", "\\1", data[seq(1, length(data), 4)])),  stringsAsFactors=F, "base"=substr(data[seq(2, length(data), 4)], base, base))



unique <- unique(data2$tile)


pdf(paste("R2/tile_R2_", k, ".pdf", sep=""))

for(i in 1:length(unique)){

temp <- data2$base[data2$tile==unique[i]]
temp[temp=="A"] <- color[1]
temp[temp=="G"] <- color[2]
temp[temp=="C"] <- color[3]
temp[temp=="T"] <- color[4]
temp[temp=="N"] <- color[5]

plot(data2$x[data2$tile==unique[i]], data2$y[data2$tile==unique[i]], main=paste("tile: ", unique[i], sep="", " - clusters: ", length(data2$x[data2$tile==unique[i]])), col= temp)
}

dev.off()

}






}

