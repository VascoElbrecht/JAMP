# work in progres - not working jet
OTU_plots <- function(file=table, out=""){



subsample <- function(sequences, N){
exp <- NULL

#sampling of OTU table by sequences depth defined below
for (i in 1:N){
temp <- sample(data$BIN.URI, sequences, prob=data$Sequences, replace=T)
exp[i] <- length(sort(table(temp)))
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

write.csv(tab, file=paste(g,"_",sub("/Users/tbraukma/Desktop/mBRAVE_test/","",files[g]),sep=""))

}




}

















