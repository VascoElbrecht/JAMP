
Denoise_plots <- function(file="E_highlight.csv", out="")





data <- read.csv(file, stringsAsFactors=F)
head(data)
data[data=="low_Haplo"] <- 0



pdf(paste(out, "/", "Haplo_per_OTU.pdf",sep=""), height=6, width=7)
barplot(table(table(data$OTU)), yaxt="n", ylab="Number of OTUs", xlab="Number of haplotypes within respective OTU")
axis(2, las=1)
dev.off()

sum(table(table(data$OTU)))


# "white",
colors <- c( "#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B", "#9E67AB", "#CE7058", "#D77FB4", "#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B", "#9E67AB", "#CE7058", "#D77FB4", "#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B", "#9E67AB", "#CE7058", "#D77FB4", "#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B", "#9E67AB", "#CE7058", "#D77FB4", "#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B", "#9E67AB", "#CE7058", "#D77FB4", "#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B", "#9E67AB", "#CE7058", "#D77FB4")


head(data[1:10])
#  make large plot showing haplotype compossition per sample

OTUs <- unique(data$OTU)

OTUs[i]


i <- 31
k <- 4
pdf(paste(out, "/", "Haplo_distribution_BETA_needs improvement.pdf",sep=""), height=7, width=6)
par(mfrow=c(5, 40), oma=c(0,0,0,0), mar=c(0,0,1,0), lwd=0.5)

for (i in 1:(length(OTUs)-1) ){
temp <- data[data$OTU==OTUs[i], c(4:(ncol(data)-3))]

for (k in 1:ncol(temp)){ # alculate relative abundance
temp[k] <- as.numeric(unlist(temp[k]))/sum(as.numeric(unlist(temp[k])))*100
}
temp[is.na(temp)] <- 0

# relative proportions
temp <- rowSums(temp)/sum(rowSums(temp))*100
temp <- temp[temp!=0]


barplot(cbind(sort(temp, decreasing=T)), col= colors, main = sub("OTU_", "", OTUs[i]), xlab="", ylab="", yaxt="n", axisnames=F, cex.main=0.7)

}
dev.off()

















