library(reshape)
library(Biostrings)
library(ggplot2)
library(gridExtra)
library(stringr)
library(dplyr)
library(googlesheets4)
setwd("~/Desktop/Data/2023/ML/T3PKS/final/")


#function to load input_df and melt it to long format,replace y/n. values to 0/1
melt_input_df <- function(input_df, range_1){
  df_long <- melt(input_df, id.vars = colnames(input_df[, c(1:range_1)]), measure.vars = colnames(input_df[, (range_1+1):length(input_df)]))
  df_long <- subset(df_long, is.na(value) == F)
  for( i in 1:length(df_long$value)){
    if(df_long$value[i] == 'n'){
      df_long$value[i] <- 0
    }
    else{df_long$value[i] <- 1}
    
  }
  return(df_long)
}
#writing fasta from df
write_fasta <- function(input_df, output_fasta){
  for(i in 1:length(input_df$Gene)){
    header <- paste('>', input_df$Gene[i])  
    write(header, file = output_fasta, append=TRUE)
    write(input_df$seq[i], file = output_fasta, append=TRUE)
  } 
}
#plotting heatmap
plot_heatmap <- function(input_df, x, y, fill){
  g <- ggplot(input_df)+ geom_tile(aes(x = x, y = y, fill = fill), linewidth = 0.4, color = "black")+
    coord_fixed()+ theme( legend.position = "none",panel.background  = element_blank(),
                          axis.ticks = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
                          text = element_text( face = "plain", size=14), axis.text.x  = element_text(angle = 90, vjust = 0.5))
  return(g)
}

make_df_for_prediction <- function(input, substrates){
  seqs <- readAAStringSet(input, format = "fasta")
  duplicated_names <- seqs@ranges@NAMES[which(duplicated(seqs@ranges@NAMES) == T)] 
  mines_names <- unique(seqs@ranges@NAMES)
  Gene <- c()
  Substrate <- c()
  for(i in 1:length(mines_names)){
    for(j in 1:length(substrates$Substrate)){
      Gene <- c(Gene, mines_names[i])
      Substrate <- c(Substrate, substrates$Substrate[j])
    }
  }
  df_mined <- as.data.frame(cbind(Gene, Substrate))
  df_mined <- subset(df_mined, !(Gene %in% duplicated_names))
  return(df_mined)
  
}

make_df_substrate_prediction <- function(genes, smiles){
  df <- as.data.frame(c())
  for(i in genes){
    for(j in smiles){
      df <- rbind(df, c(i, j))
    }
  }
  colnames(df) <- c("Gene", "Substrate")
  return(df)
}

#data_set with 
input_df <- as.data.frame(read_sheet("https://docs.google.com/spreadsheets/d/17XtMcsRywdFOYnWKTbYyUnTUB1Qz6EFiCep7AoY1Bag/edit?usp=sharing"))
#separating training set
training_df <- input_df[,1:16]
training_long <- melt_input_df(training_df, 5)
colnames(training_long) <- c(colnames(training_df[1:5]), "Substrate", "value")
#separatinge xtra cell_free reaction
extra_substrate <- as.data.frame(read_sheet("https://docs.google.com/spreadsheets/d/1l966bZB2yfbYZ85Q2toriq7F18Y1eDrjlP0Z1BYVhMQ/edit?usp=sharing"))
extra_long <- melt_input_df(extra_substrate, 4)
colnames(extra_long) <- c(colnames(extra_long[1:4]), "Substrate", "value")
#reading published data
published_df <- as.data.frame(read_sheet("https://docs.google.com/spreadsheets/d/1-U2FxTf1-KSYlZcB8p4dUCpUgbVMoHdNphNiSqhwgtM/edit?usp=sharing"))
long_published_df <- melt_input_df(published_df[,1:32], 3)
colnames(long_published_df) <- c(colnames(long_published_df[1:3]), "Substrate", "value")

#make dataframe for substrate prediction 
substrates_prediction_df <- make_df_substrate_prediction(unique(training_long$Gene), unique(extra_long$Substrate))



#wirting dataset for ML
write.csv(training_long, "input_data/T3PKS_cellfree_df.csv", sep = ",")
write.csv(extra_long, "input_data/xtra_substrates.csv", sep = ",")
write.csv(long_published_df, "input_data/published.csv", sep = ",")
write.csv(substrates_prediction_df, "input_data/substrate_prediction_df.csv", sep = ",")
#write fasta for published
write_fasta(published_df, "published.fasta")
