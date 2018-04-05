Denoise_barplot <- function(table="E_haplo_table.csv", out=NA, height=6, width=7, mfrow=c(5, 40), emptyOTUs=T){


if(is.character(table)){ # standard JAMP table
data <- read.csv(table, stringsAsFactors=F)
data <- data[-nrow(data), -c(1,2, ncol(data))]
} else {data <- table}



# "white",
colors <- c("#ffffff", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", rep("gray", 100))



head(data)
#  make large plot showing haplotype compossition per sample

OTUs <- unique(data[1])

i <- 31
k <- 4

if(!is.na(out)){
pdf(file= out, height, width)
}
par(mfrow=c(mfrow[1], mfrow[2]), oma=c(0,0,0,0), mar=c(0,0,1,0), lwd=0.5)

for (i in 1:length(OTUs) ){
temp <- data[data[1]==OTUs[i], c(2:ncol(data))]

for (k in 1:ncol(temp)){ # alculate relative abundance
temp[k] <- suppressWarnings(as.numeric(unlist(temp[k])))/sum(suppressWarnings(as.numeric(unlist(temp[k]))))*100
}
temp[is.na(temp)] <- 0

cols <- length(temp)*100

# relative proportions
#temp <- rowSums(temp)/sum(rowSums(temp))*100
temp <- rowSums(temp) # count coloumns
temp <- temp[temp!=0]

if(emptyOTUs){
barplot(cbind(c(cols-sum(temp), sort(temp, decreasing=T))), col= colors, main = sub("OTU_", "", OTUs[i]), xlab="", ylab="", yaxt="n", axisnames=F, cex.main=0.7)
} else {
barplot(cbind(sort(temp, decreasing=T)), col= colors[-1], main = sub("OTU_", "", OTUs[i]), xlab="", ylab="", yaxt="n", axisnames=F, cex.main=0.7)
}

}
if(!is.na(out)){
dev.off()
}

}



