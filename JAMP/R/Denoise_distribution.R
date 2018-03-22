# 171031

setwd("~/Documents/UNI_und_VORLESUNGEN/14 Guelph/1 Parks /1 R scripts/1 haplo per OTU/finland")



data <- read.csv("~/Documents/UNI_und_VORLESUNGEN/11 phd projects/1 Meta HAPLOTYPING/2 JAMP finland/r maps/_haplomaps/highligh_orange/170825_BF2+BR2_haplotable_alpha_5_haplosize_0.01_minOTU_0.1_orange.csv", stringsAsFactors=F)

data[data=="orange"] <- 0

pdf("hist_haplo_finland_BF2+BR2_180128.pdf", height=6, width=7)
barplot(table(table(data$OTU)), yaxt="n", ylab="Number of OTUs", xlab="Number of haplotypes within respective OTU")
axis(2, las=1)
dev.off()

sum(table(table(data$OTU)))


# add.pies
colors <- c("#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B", "#9E67AB", "#CE7058", "#D77FB4")



#  make large plot showing haplotype compossition per sample

OTUs <- unique(data$OTU)

i <- 9
k <- 1
pdf("haplo_composition_finland_BF2+BR2_180128+2.pdf", height=7, width=6)
par(mfrow=c(5, 40), oma=c(0,0,0,0), mar=c(0,0,1,0), lwd=0.5)

for (i in 1:length(OTUs) ){
temp <- data[data$OTU==OTUs[i], -c(1:4, ncol(data))]

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

















