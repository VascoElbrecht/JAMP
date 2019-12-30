# beta
OTU_heatmap <- function(file=table, out="", abundance=F, plot0=F, rel=F, col=c("#2b83ba", "#abdda4", "#ffffbf", "#fdae61", "#d7191c")){

ORIGscipen <- getOption("scipen")
options(scipen=10)

if(is.character(file)){
data <- read.csv(file, stringsAsFactors=F)
row.names(data) <- data[,2]
data <- data[,-c(1,2, ncol(data))]
} else {
data <- file
row.names(data) <- data[,2]
data <- data[-c(1,2)]}

backup <- data


# make relatuve abundance
sums <- colSums(data)

for (k in 1:ncol(data)){
data[k] <- data[k]/sums[k]*100
if(rel){ # save as relative abundance
backup[k] <- round(data[k], 4)
backup[k][is.na(backup[k])] <- 0
}
data[k] <- suppressWarnings(log10(data[k])) 			# log10
data[,k][!is.finite(data[,k])] <- NA	# remove inf
data[,k][data[,k]< -3] <- -3
}


#invert
data <- data[nrow(data):1,]
backup <- backup[nrow(data):1,]




mycol <- colorRampPalette(col)
ramp <- data.frame("ID"=rev(seq(-3, 2, 0.1)), "col"=mycol(52)[1:51], stringsAsFactors=F)


##### MAKE PLOT

if(out!=""){
pdf(out, height=(nrow(data)+20)/10, width=(ncol(data)+2)/2)
}

par(mar=c(0,0,0,0))
plot(NULL, ylim=c(c(nrow(data)+20)*0.035, c(nrow(data)+20)*0.965), xlim=c(0.75, ncol(data)+2), xlab="", ylab="", xaxt="n", yaxt="n", bty="n")


for (m in 1:ncol(data)){
text(m, nrow(data)+2, names(data)[m], srt=45, adj=0)
for(i in 1:nrow(data)){
wichColor <- which(data[i,m]>ramp$ID)[1]
rect(m-0.5, i-0.5, m+0.5, i+0.5, col=ramp$col[wichColor], border=F)
if(abundance){
	if(plot0){text(m, i, backup[i,m], cex=0.5)} else if(backup[i,m]>0){text(m, i, backup[i,m], cex=0.5)}
	}
}
}

text(ncol(data)+0.6, 1:nrow(data), row.names(data), adj=0, cex=0.5)

# add legend

for(i in 1:51){
rect(0.45+i*0.05, nrow(data)+17, 0.5+i*0.05, nrow(data)+20, col=ramp$col[i], border=F)
}

MYlabels <- c("100%", "10%", "1%", "0.1%", "0.01%", "0.001%")

k <- 1

for (i in c(0.525, 1.025, 1.525, 2.025, 2.525, 3.025)){
lines(c(i, i), c(nrow(data)+17, nrow(data)+16.5))
text(i, nrow(data)+16.3, labels=MYlabels[k], srt=90, adj=1, cex=0.5)
k <- k+1
}




i <- 1



0.45+0.025+0.05*11

if(out!=""){
dev.off()
}
options(scipen= ORIGscipen)
}# heatmap end






if(F){


subsample <- function(sequences, N){
exp <- NULL

#sampling of OTU table by sequences depth defined below
for (i in 1:N){
data <- sample(data$BIN.URI, sequences, prob=data$Sequences, replace=T)
exp[i] <- length(sort(table(data)))
}
return(exp)
}

#number of subreads sampled cannot exceed the number of availble reads
steps <- 10^c(seq(0,7, 0.1))
steps[steps<sum(data$Sequences)]

#table will contain number of sequences samples, mean bin count with st dev for all replicates
tab <- data.frame("seqNumber"=1, "meanBIN"=1, "SDbin"=1)
tab <- tab[-1,]
#k is variable for sequences called by the function above
for (k in steps[steps<sum(data$Sequences)]){
#N=number of desierd replicates
subset <- subsample(k, N=50)
tab <- rbind(tab, cbind(k, mean(subset), sd(subset)))
}

write.csv(tab, file=paste(g,"_",sub("path","",files[g]),sep=""))

}



