library(reshape)
library(Biostrings)
library(ggplot2)
library(gridExtra)
library(ggseqlogo)
library(msa)
library(stringr)
setwd("~/Desktop/Data/2024/Bioinformatics/T3PKS_msa/")
#loading plant.faa as reference and all.fa for fungi t3pks
plant <- readAAStringSet("plant-pks-reviewed.fasta", format = "fasta")
t3pks <- readAAStringSet("all_unique.fa", format = "fasta")
#read alignment
msa <- readAAMultipleAlignment("all_unique_rev_plannt_msa", format = "clustal")
msa_seq <- readAAStringSet(msa)
#important residues in the first plant seq C164,H303,N336, F215, F265, T132, S133, T194 ,S338 ,P375 ,G256
#MsCHS AAA02824.1
MsCHS_res <- c(164, 303,336,215,265,132,133,194,338,375,256)
MsCHS_res <- sort(MsCHS_res)
#indexing msa finding  AAA02824.1
MsCHS_msa <- unlist(str_split(toString(msa@unmasked[which(msa@unmasked@ranges@NAMES == "AAA02824.1")]), ""))
index <- 0
index_msa <- 0
res_msa <- c()
for( i in 1:length(MsCHS_msa)){
  index_msa <- index_msa + 1
  if(MsCHS_msa[i] != "-"){
    index <- index +1
    if(index %in% MsCHS_res){
      res_msa<- c(res_msa, index_msa)
    }
  }
}
#msa <- chartr("-", "X", msa)
#making simple seqlogo
#crashes if use entire sequence - slice according positions first, then pass to seqlogo
sliced_seq <- c()
for( i in 1:length(msa@unmasked)){
  seq <- unlist(str_split(toString(msa@unmasked[i]), ""))
  sliced_seq <- c(sliced_seq, str_c(seq[res_msa], collapse = ""))
}
#seqlogo
#making labels
labels <- c()
MsCHS <- unlist(str_split(str_replace_all( toString(msa@unmasked[which(msa@unmasked@ranges@NAMES == "AAA02824.1")]), "-", ""), ""))
for( i in 1:length(MsCHS_res)){
  labels <- c(labels, paste(MsCHS[MsCHS_res[i]], MsCHS_res[i], sep = ""))
}

plant_logo <- ggseqlogo(sliced_seq[1:length(plant)+1], method = "p") + 
  scale_x_discrete(breaks = 1:11,labels = labels, limits = factor(1:11))+labs(title = "Plant")
fungi_logo <-  ggseqlogo(sliced_seq[185:length(sliced_seq)], method = "p")+ 
  scale_x_discrete(breaks = 1:11,labels = labels, limits = factor(1:11))+ labs(title = "Fungi")
gplant <- ggplotGrob(plant_logo)
gfungi <- ggplotGrob(fungi_logo)
grid.arrange(gplant, gfungi, nrow = 2)


cairo_ps(filename = "seqlogo.eps",
         width = 50, height = 20, pointsize = 12,
         fallback_resolution = 600)
print(grid.arrange(gplant, gfungi, nrow = 2)
)
dev.off()

#writing sliced fasta files for pLogo
for(i in 1:27){
  header <- paste('>', msa@unmasked@ranges@NAMES[i])  
  write(header, file = "plant_sliced.fa", append=TRUE)
  write(sliced_seq[i], file =  "plant_sliced.fa", append=TRUE)}

for(i in 28:length(sliced_seq)){
  header <- paste('>', msa@unmasked@ranges@NAMES[i])  
  write(header, file = "fungi_sliced.fa", append=TRUE)
  write(sliced_seq[i], file =  "fungi_sliced.fa", append=TRUE)}


writeXStringSet(sliced_seq[1:27], "plant_align.fa", format = "fasta")
writeXStringSet(msa@unmasked[28:length(sliced_seq)], "fungi_align.fa", format = "fasta")


