# 200402 - lulu processing

LULU <- function(table=NA, folder="Lulu", ratio_type = "min", min_ratio = 1, min_match = 84, cooccurence = 0.95, exe="vsearch", heatmap=T){

# check if folder exists
if(!file.exists(folder)){ # folder does not exist
dir.create(folder)
message(paste("Cerated folder \"", folder, "\"\n", sep=""))
} else { #Folder exists, ask to delete it.


tempSelect <- menu(c("y", "n"), title=paste("Folder \"", folder, "\" elready exists. Do you wnat to delete and over write it? y = Yes, n = No", sep=""))

if(tempSelect==1){
system2("rm", paste(" -fr ", folder, sep=""))
message("folder deleted!\n")
message(paste("Creating folder \"", folder, "\"\n", sep=""))
dir.create(folder)
}
if(tempSelect==2){
message("Nothing deleted! Function stopped!")
stop("User aborted!")
}
}



# read OTU table
data <- read.csv(table, stringsAsFactors=F)

data$sequ <- toupper(data$sequ)

# need to implement for ESVs as well!
row.names(data) <- data$ID

keep <- !names(data) %in% c("ID", "sequ", "sort", "ESV")

# write table OTU_IDs + abundance in each sample
otuonly <- data[, keep]

#write.table(otuonly, paste(folder, "/OTUonly.tsv", sep=""), quote=F, sep="\t")


# write OTUs as CSV

cat(file=paste(folder, "/raw_OTU.fasta", sep=""),paste(paste(">", data$ID, sep=""), toupper(data$sequ), sep="\n"), sep="\n")


# Match OTUs against them self
system2("vsearch", paste("--usearch_global ", folder, "/raw_OTU.fasta --db  ", folder, "/raw_OTU.fasta --self --id .84 --iddef 1 --userout  ", folder, "/match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10", sep=""))



# Run lulu
matchlist <- read.csv(paste(folder, "/match_list.txt", sep=""), sep="\t", stringsAsFactors=F)


message("applying lulu filtering")

curated_result <- lulu(otuonly, matchlist, minimum_ratio_type = ratio_type, minimum_ratio = min_ratio, minimum_match = min_match, minimum_relative_cooccurence = cooccurence)


# merge tables for statistics
message("Done!")

# move lulu log
system2("mv", paste(list.files(pattern="lulu.log_"), " ", sep="", folder))




message("OTUs total:     ", nrow(data))
message("OTUs curated:   ", curated_result$curated_count, " (", round(curated_result$curated_count/nrow(data)*100, 2), "%).")
message("OTUs discarded: ", curated_result$discarded_count, " (", round(curated_result$discarded_count/nrow(data)*100, 2), "%)." )

message(paste("Writing raw table with lulu scoring! Can be used for more accurate processing once taxonomy is included.\nFile: ", folder, "/1 Raw_lulu_flags.csv", sep=""))

originaltable <- data.frame("ID"=row.names(curated_result$original_table), curated_result$original_table, "sequ"=data$sequ[match(row.names(curated_result$original_table), data$ID)], curated_result$otu_map, "keept" = row.names(curated_result$original_table) %in% curated_result$curated_otus, stringsAsFactors=F)


# sort by OTU number
originaltable <- originaltable[order(as.numeric(sub("OTU_", "", originaltable$ID)), decreasing=F),]

originaltable <- data.frame("sort"=1:nrow(originaltable), originaltable, stringsAsFactors=F)

# write OTU flaggs of raw data!
write.csv(originaltable, paste(folder, "/1 Raw_lulu_flags.csv", sep=""), row.names=F)


# make fasta file for each currated OTU cluster
message("Making individual fasta files", if(heatmap){" and heatmaps"}, " for collapsed OTUs.")

dir.create(paste(folder, "/1 OTU_cluster", sep=""))
if(heatmap){dir.create(paste(folder, "/1 OTU_cluster_plots", sep=""))}

uniOTU <- unique(originaltable$parent_id)
for(i in 1:length(uniOTU)){
temp <- originaltable[uniOTU[i]==originaltable$parent_id,]

cat(file=paste(folder, "/1 OTU_cluster/", uniOTU[i], ".fasta", sep=""),paste(paste(">", temp$ID, ";size=", temp$total, sep=""), toupper(temp$sequ), sep="\n"), sep="\n")

if(heatmap){
pdf(paste(folder, "/1 OTU_cluster_plots/", uniOTU[i], ".pdf", sep=""), height=(nrow(temp)+40)/10, width=(ncol(temp)-1)/3)
OTU_heatmap(temp[, -c(1,2,(ncol(temp)-6):ncol(temp))], abundance=T, rel=F, col=rev(c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")))
dev.off()
}
} # loop end!

if(heatmap){
message("Making heatmap for raw data file, sorted by OTUs merged together")

temp <- originaltable[order(originaltable$total, decreasing=T),]
temp <- temp[order(as.numeric(sub("OTU_", "", temp$parent_id))),]
row.names(temp) <- paste(temp$parent_id, temp$ID, temp$total, sep="   ")
temp <- temp[, -c(1,2,(ncol(temp)-6):ncol(temp))]

pdf(paste(folder, "/1 Raw_lulu_flags_plot.pdf", sep=""), height=(nrow(temp)+20)/10, width=(ncol(temp)-1)/3)
OTU_heatmap(temp, abundance=T, rel=F, col=rev(c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")))
dev.off()
}

# writing condensed OTU table!
message("Writing condensed OTU table!")


data2 <- curated_result$curated_table
data2 <- data2[order(as.numeric(sub("OTU_", "", row.names(data2))), decreasing=F),]



if(heatmap){
message("Making heatmap for condensed OTUs.")

pdf(paste(folder, "/2 condensed_OTUs_plot.pdf", sep=""), height=(nrow(temp)+20)/10, width=(ncol(temp)-1)/3)
OTU_heatmap(data2, abundance=T, rel=F, col=rev(c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba")))
dev.off()
}


data3 <- data.frame("sort"=1:nrow(data2), "ID"= row.names(data2), data2, "sequ"=data$sequ[match(row.names(data2), data$ID)], stringsAsFactors=F)


write.csv(data3, paste(folder, "/2 condensed_OTUs.csv", sep=""), row.names=F)

cat(file=paste(folder, "/2 condensed_OTUs.fasta", sep=""),paste(paste(">", data3$ID, sep=""), toupper(data3$sequ), sep="\n"), sep="\n")

message("Condensed OTU table and fasta file written. Processing completed!")

}














