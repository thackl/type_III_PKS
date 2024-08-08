library(ggplot2)
library(gridExtra)
library(stringr)
library(dplyr)
library(scales)
setwd("~/Desktop/Data/2023/ML/T3PKS/final/")

gene_levels <- c("PhCHS", "AtamPKS2", "AserPKS1", "AserPKS2", "HargPKS1", "XacuPKS1",
                 "AwakPKS", "AtriPKS", "AtamPKS1", "AsesPKS", "AastPKS", "AcosPKS",
                 "PtriPKS", "XylPKS", "AthePKS", "TtonPKS", "PflaPKS", "AlupPKS2", "AlupPKS1",
                 "AiizPKS", "AneoPKS", "MpolPKS", "FerePKS", "CadPKS", "PverPKS", "SinsPKS",
                 "CgloPKS", "MoryPKS", "VmalPKS", "DhelPKS", "DliqPKS", "FmanPKS", "PficPKS",
                 "HargPKS2", "XacuPKS2", "BiscPKS", "DdecPKS", "HypPKS")

substrates_levels <- c('4-pentynoyl-coa',"cyclohexanoyl-coa"  ,"4-cyclohexyl-4-oxobutyryl-coa",
                       "4-cyclopentyl-4-oxobutyryl-coa" ,"thiazole-4-carboxylyl-coa","2-furoyl-coa" ,
                       "trans-2-Phenylcyclopropane-1-carboxylyl-coa","cinnamoyl-coa","2-chlorophenylacetyl-coa",
                       "3-fluorobenzoyl-coa","4-fluorophenylacetyl-coa","2-amino-5-chlorobenzoyl-coa" )


#box plot for model scores
scores_training <- as.data.frame(read.csv2("output_data/model_scores_training.csv", sep = ","))
scores_training$set <- "training"

#new box plot for SI with training models

models_si <- ggplot(subset(scores_training, Model != "Extra Trees"), aes(x = Model, y = as.numeric(value), fill = Embedding))+
  stat_summary(fun = "mean",fun.max = function(x) mean(x) + sd(x),  
               fun.min = function(x) mean(x) - sd(x),
               geom = "errorbar", position="dodge", size = 0.5)+
  stat_summary(fun = "mean",  size = 3, geom = "bar", position="dodge")+
  scale_y_continuous(limits = c(0.0, 1.05),oob = rescale_none, expand = c(0,0), 
                     breaks = seq(0,1.05, by = 0.2))+
  theme(panel.background  = element_blank(), axis.line = element_line(),
        axis.title = element_blank(), axis.text = element_text(size = 20), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "right")+
  theme(panel.background  = element_blank(), axis.line = element_line(),
                       axis.title = element_blank(), axis.text = element_text(size = 14), 
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
                       legend.position = "right")+
  scale_fill_manual(values = c("#B05670", "#005b96", "#005f43"), name = "Embeding")+
  facet_grid( cols = vars(Metric))
models_si

published_score <- as.data.frame(read.csv2("output_data/published_scores.csv", sep = ","))
published_score$set <- "published"
published_score$value <- as.numeric(published_score$value)
substrate_score <- as.data.frame(read.csv2("output_data/extra_substrates_scores.csv", sep = ","))
substrate_score$set <- "substrates"
substrate_score$value <- as.numeric(substrate_score$value)
plant_score <- as.data.frame(read.csv2("output_data/published_plant_scores.csv", sep =","))
plant_score$value <- as.numeric(plant_score$value)


score_df <- rbind(scores_training,published_score, substrate_score)


score_df$value <- as.numeric(score_df$value)
score_df$Metric <- factor(score_df$Metric, levels = c("Accuracy",  "ROC AUC", "Precision","Recall","F1"))
score_df$set <- factor(score_df$set, levels = c("training", "published", "substrates"))

#plot for the main figure
model_score <- ggplot(subset(score_df, Model == "Multi-layer Perceptron" & Embedding == "MACCS/ProtTrans"), aes(x = Metric, y = value, fill = set)) + 
  geom_boxplot()+
  scale_y_continuous(limits = c(0.5, 1.05),oob = rescale_none, expand = c(0,0))+
  theme(panel.background  = element_blank(), axis.line = element_line(),
        axis.title = element_blank(), axis.text = element_text(size = 20), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "right")+
  scale_fill_manual(values = c("#B05670", "#005b96", "#005f43"), name = "Set", labels = c("Test", "Extra substrates", "Published fungal"))
model_score

model_score <- ggplot(subset(score_df, Model == "Multi-layer Perceptron" & Embedding == "MACCS/ProtTrans"),
                      aes(x = Metric, y = value, fill = set, group = set)) +
  stat_summary(fun = "mean",fun.max = function(x) mean(x) + sd(x),  
               fun.min = function(x) mean(x) - sd(x),
               geom = "errorbar", position="dodge", size = 0.5)+
  stat_summary(fun = "mean",  size = 3, geom = "bar", position="dodge")+
  scale_y_continuous(limits = c(0.0, 1.05),oob = rescale_none, expand = c(0,0), 
                     breaks = seq(0,1.05, by = 0.2))+
  theme(panel.background  = element_blank(), axis.line = element_line(),
        axis.title = element_blank(), axis.text = element_text(size = 20), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "right")+
  scale_fill_manual(values = c("#B05670", "#005f43", "#005b96"), name = "Set", 
                    labels = c("Test", "Published fungal", "Extra substrates"))
model_score

cairo_ps(filename = "figures/ML_scores.eps",
         width = 5.5, height = 4.9, pointsize = 12,
         fallback_resolution = 600)
print(model_score)
dev.off()
gModel <- ggplotGrob(model_score)
#stat test for training dataset
library(tidyverse)
library(rstatix)
library(ggpubr)

stats_df <- subset(scores_training, Model != "Extra Trees")
stats_df$Model_embedding <- paste(stats_df$Embedding, stats_df$Model, sep ="_")
stats_df$value <- as.numeric(stats_df$value)
stats_df$Model_embedding <- factor(stats_df$Model_embedding, 
                                   levels = c("Morgan/ProtTrans_Decision Tree" , "RDKit/ProtTrans_Decision Tree" ,
                                              "MACCS/ProtTrans_Decision Tree","Morgan/ProtTrans_Random Forest"  ,
                                              "RDKit/ProtTrans_Random Forest","MACCS/ProtTrans_Random Forest" ,
                                              "Morgan/ProtTrans_Multi-layer Perceptron" ,"RDKit/ProtTrans_Multi-layer Perceptron" ,
                                              "MACCS/ProtTrans_Multi-layer Perceptron" ) )


#colnames(data) <- c("V1", "V2", "V3")
stat.test <- stats_df %>%
  group_by(Metric) %>%
  t_test(value ~ Model_embedding) %>%
  #adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

myplot <- ggboxplot(
  stats_df, x = "Metric", y = "value",
  fill = "Model_embedding", palette = "npg",legend = "right",
  ggtheme = theme_pubr(border = TRUE),
  repel = T) + rotate_x_text(angle = 90, hjust = NULL, vjust = NULL)
myplot
stat.test <- stat.test %>% add_xy_position(x = "Metric")
myplot + stat_pvalue_manual(hide.ns = TRUE,stat.test, label = "p.adj.signif", step.increase = 0)
#descriptive statistic
library(psych)
statistics_training <- describeBy(stats_df, list(stats_df$Model_embedding,stats_df$Metric))
statistics_published  <- describeBy(published_score, list(published_score$Embedding, published_score$Model,
                                                          published_score$Metric))
statistics_plant <- describeBy(plant_score, list(plant_score$Embedding, plant_score$Model,
                                                 plant_score$Metric))
statistics_substrate <- describeBy(substrate_score, list(substrate_score$Embedding, substrate_score$Model,
                                                         substrate_score$Metric))



statistics_published

ggsave("figures/ML_models_training.png", 
       myplot + stat_pvalue_manual(hide.ns = TRUE,stat.test, label = "p.adj.signif"),
       width = 50, height = 20,units = "cm", dpi = 600)


#heatmap plot for main figures
plot_heatmap <- function(input_df, x, y, fill){
  g <- ggplot(input_df)+ geom_tile(aes(x = x, y = y, fill = fill), linewidth = 0.4, color = "black", size = 5)+
    coord_fixed()+ theme( legend.position = "right",panel.background  = element_blank(),
                          axis.ticks = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
                          text = element_text( face = "plain", size=14), axis.text.x  = element_text(angle = 90, vjust = 0.2, hjust = 0.95))
  return(g)
}
#reading results for published enzymes and adding a columnm for substrates in training data set
published_heatmap_df <- as.data.frame(read.csv2("output_data/published_prediction.csv"))
in_training <- c()
for( i in published_heatmap_df$Substrate){
  if( i %in% unique(training_long$Substrate)){
    in_training <- c(in_training, "X")
  }
  else{
    in_training <- c(in_training, "")
  }
}
published_heatmap_df$intraining <- in_training

published_exp <- plot_heatmap(published_heatmap_df, published_heatmap_df$Substrate, published_heatmap_df$Gene,
                              as.numeric(published_heatmap_df$value))+ scale_fill_gradient(low = "white", high = "#005f43")+
  theme(legend.position = "right")+labs(fill = "Probability")
published_exp
published_pred <- plot_heatmap(published_heatmap_df, published_heatmap_df$Substrate, published_heatmap_df$Gene,
                               as.numeric(published_heatmap_df$MACCS.ProtTrans_Multi.layer.Perceptron))+ scale_fill_gradient(low = "white", high = "#B05670")+
  theme(legend.position = "right")+labs(fill = "Probability")+
  geom_text(aes(x = published_heatmap_df$Substrate, y = published_heatmap_df$Gene,label = published_heatmap_df$intraining))
published_pred
gpublished_exp <- ggplotGrob(published_exp)
gpublished_pred <- ggplotGrob(published_pred)
#substrates
substrate_heatmap_df <- as.data.frame(read.csv2("output_data/extra_substrates_prediction.csv"))
substrate_heatmap_df$Substrate <- factor(substrate_heatmap_df$Substrate, levels = substrates_levels)
substrate_exp <- plot_heatmap(substrate_heatmap_df, substrate_heatmap_df$Substrate, substrate_heatmap_df$Gene,
                              substrate_heatmap_df$value)+ scale_fill_gradient(low = "white", high = "#005b96")+
                  scale_x_discrete(labels = seq(1,length(unique(substrate_heatmap_df$Substrate)),1))+theme(axis.text.x = element_text(angle = 0))+
  labs(fill = "Probability")
substrate_pred <- plot_heatmap(substrate_heatmap_df, substrate_heatmap_df$Substrate, substrate_heatmap_df$Gene,
                               as.numeric(substrate_heatmap_df$MACCS.ProtTrans_Multi.layer.Perceptron))+ scale_fill_gradient(low = "white", high = "#B05670")+
  scale_x_discrete(labels = seq(1,length(unique(substrate_heatmap_df$Substrate)),1))+theme(axis.text.x = element_text(angle = 0))+
  labs(fill = "Probability")

gsubstrate_exp <- ggplotGrob(substrate_exp)
gsubstrate_pred <- ggplotGrob(substrate_pred)

lay <- rbind(c(1,2,3),
             c(4,5,9))
grid.arrange(gModel,gsubstrate_exp, gsubstrate_pred,gpublished_exp, gpublished_pred,  layout_matrix = lay )
cairo_ps(filename = "figures/ML_fig1.eps",
         width = 20, height = 20, pointsize = 12,
         fallback_resolution = 600)
print(grid.arrange(gModel,gsubstrate_exp, gsubstrate_pred,gpublished_exp, gpublished_pred,  layout_matrix = lay )
)
dev.off()

#plants/bacterial prediction for SI figure
plant_score <- as.data.frame(read.csv2("output_data/published_plant_scores.csv", sep =","))

plant_score_plot <- ggplot(subset(plant_score, Model == "Multi-layer Perceptron" & Embedding == "MACCS/ProtTrans"), 
                           aes(x = Metric, y = as.numeric(value))) + 
  stat_summary(fun = "mean",fun.max = function(x) mean(x) + sd(x),  
               fun.min = function(x) mean(x) - sd(x),
               geom = "errorbar", position="dodge", size = 0.5)+
  stat_summary(fun = "mean",  size = 3, geom = "bar", position="dodge", fill = "#B05670")+
  scale_y_continuous(limits = c(0, 1.05),oob = rescale_none, expand = c(0,0))+
  theme(panel.background  = element_blank(), axis.line = element_line(),
        axis.title = element_blank(), axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "right")
plant_score_plot
gPlantScore <- ggplotGrob(plant_score_plot)


published_plant_df <- as.data.frame(read.csv2("output_data/published_plant_prediction.csv"))
in_training2 <- c()
for( i in published_plant_df$Substrate){
  if( i %in% unique(training_long$Substrate)){
    in_training2 <- c(in_training2, "X")
  }
  else{
    in_training2 <- c(in_training2, "")
  }
}
published_plant_df$intraining <- in_training2

plants_level <- c("ARAS1","ARAS2", "QNS-Marmelos",  "QNS-Microcarpa", "ANS-Microcarpa", "BNS-Palmatum","ArsB",
                  "ArsC","BpsA","PhlD","gcs","Sg-RppA",
                  "Sc-RppA","Ncs", "Se-RppA" )

published_plant_df$Gene <- factor(published_plant_df$Gene, levels = rev(plants_level))

published_plant <- plot_heatmap(published_plant_df, published_plant_df$Substrate, published_plant_df$Gene,
                                published_plant_df$value)+ scale_fill_gradient(low = "white", high = "#005f43")+
  theme(legend.position = "right",  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+labs(fill = "Probability")
published_plant
predicted_plant <- plot_heatmap(published_plant_df, published_plant_df$Substrate, published_plant_df$Gene,
                                as.numeric(published_plant_df$MACCS.ProtTrans_Multi.layer.Perceptron))+ scale_fill_gradient(low = "white", high = "#B05670")+
  geom_text(aes(x = published_plant_df$Substrate, y = published_plant_df$Gene,label = published_plant_df$intraining))+
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) )+labs(fill = "Probability")
predicted_plant
gpublished_plant <- ggplotGrob(published_plant)
gpredicted_plant <- ggplotGrob(predicted_plant)
lay3 <- rbind(c(1,2,2,2,3,3,3), 
              c(4,2,2,2,3,3,3))
grid.arrange(gPlantScore,gpublished_plant,gpredicted_plant, layout_matrix = lay3)
ggsave("figures/figSx_plant_bacteria.png", 
       grid.arrange(gPlantScore,gpredicted_plant, gpublished_plant,layout_matrix = lay3),
       width = 50, height = 20,units = "cm", dpi = 600)

cairo_ps(filename = "figures/figSx_plant_bacteria.eps",
         width = 50, height = 20, pointsize = 12,
         fallback_resolution = 600)
print(grid.arrange(gPlantScore,gpublished_plant,gpredicted_plant, layout_matrix = lay3))
dev.off()


#plot prediction for all enzymes and substrates
all_prediction_df <- as.data.frame(read.csv2("output_data/all_substrate_prediction.csv"))
all_prediction_df$Gene <- factor(all_prediction_df$Gene, levels = gene_levels)

all_prediction_df$Substrate <- factor(all_prediction_df$Substrate, levels = rev(substrates_levels))
all_prediction_plot <- plot_heatmap(all_prediction_df,all_prediction_df$Gene, all_prediction_df$Substrate, 
                                     as.numeric(all_prediction_df$pred))+
                        scale_fill_gradient(low = "white", high = "#B05670", name = "Probability")+
  scale_y_discrete(labels = seq(24, 13, -1))+theme(axis.text.y = element_text(face = "bold"))

all_prediction_plot

cairo_ps(filename = "figures/all_predictions.eps",
         width = 50, height = 20, pointsize = 12,
         fallback_resolution = 600)
print(all_prediction_plot)
dev.off()

#plots for descriptive model
descriptive_score <- as.data.frame(read.csv2("output_data/decision_tree_score.csv"))
descriptive_score <- melt(descriptive_score, id.vars = "X")


descriptivel_score_plot<- ggplot(descriptive_score, aes(x = variable, y = as.numeric(value)))+
  geom_boxplot(fill = "#B05670")+theme(panel.background  = element_blank(), axis.line = element_line(),
                       axis.title = element_blank(), axis.text = element_text(size = 16), 
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
                       legend.position = "right")


gDescriptivel_score <- ggplotGrob(descriptivel_score_plot)
#violin plot for Decision tree model importances
importances_df <- as.data.frame(read.csv2("output_data/importances.csv", sep = ","))
importances_df <- melt(importances_df, id.vars = "X")
importances_df$value <- as.numeric(importances_df$value)

descriptive_boxplot <-  ggplot(importances_df, aes(x = variable, y = as.numeric(value)))+
  geom_boxplot(fill = "#B05670")+theme(panel.background  = element_blank(), axis.line = element_line(),
                                       axis.title = element_blank(), axis.text = element_text(size = 40),
                                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                                       legend.position = "right")
descriptive_boxplot
gdescriptive_boxplot <- ggplotGrob(descriptive_boxplot)

cairo_ps(filename = "figures/descriptive_model.eps",
         width = 28.37838, height = 10, pointsize = 20,
         fallback_resolution = 600)
print(descriptive_boxplot )
dev.off()


#making alignment plot
msa_df_input <- as.data.frame(read.csv2("output_data/msa_ligand_cys.csv", sep = ","))
msa_df <- melt(msa_df_input, id.vars = c("X", "Importance"))
msa_df <- subset(msa_df, !(variable %in% c("ORAS", "MsCHS")))
#clustal color scheme
AA <- c("A", "I", "L", "M","F", "W", "V", "R", "K",
        "N", "Q", "S", "T", "D", "E", "C", "G", "H", "Y", "-",
        "P")
AA_pallete <- c("#7fa0f0","#7fa0f0","#7fa0f0","#7fa0f0","#7fa0f0","#7fa0f0", "#7fa0f0","#f01405","#f01405", 
                "#02ff00","#02ff00","#02ff00","#02ff00", "#c048c0", "#c048c0", "#f08080", "#f09047", "#14a4a4", "#14a4a4", "white", 
                "#ffff00")

pallete_vec <- c()
for(i in msa_df$value){
  pallete_vec <- c(pallete_vec,AA_pallete[which(AA == i)] )
}
#gene levels from phylogenetic tree
gene_levels <- c("PhCHS", "AtamPKS2", "AserPKS1", "AserPKS2", "HargPKS1", "XacuPKS1",
                 "AwakPKS", "AtriPKS", "AtamPKS1", "AsesPKS", "AastPKS", "AcosPKS",
                 "PtriPKS", "XylPKS", "AthePKS", "TtonPKS", "PflaPKS", "AlupPKS2", "AlupPKS1",
                 "AiizPKS", "AneoPKS", "MpolPKS", "FerePKS", "CadPKS", "PverPKS", "SinsPKS",
                 "CgloPKS", "MoryPKS", "VmalPKS", "DhelPKS", "DliqPKS", "FmanPKS", "PficPKS",
                 "HargPKS2", "XacuPKS2", "BiscPKS", "DdecPKS", "HypPKS")
msa_df$variable <- factor(msa_df$variable, levels = gene_levels)
msa_df$X <- factor(msa_df$X, levels = msa_df_input$X)
msa_df$color <- pallete_vec
#filter only important residue

msa_plot <- ggplot(msa_df) + 
  geom_tile(aes(x = X, y = variable), fill = msa_df$color, alpha = 20*as.numeric(msa_df$Importance)) +
  #scale_fill_manual(values = pallete_vec)+
 geom_text(aes(x = X, y = variable, label = value))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.05),
        panel.background = element_blank(), axis.title.y = element_blank())+
  labs(x = "ORAS residues")+
  coord_fixed()
msa_plot
gmsa <- ggplotGrob(msa_plot)

#making seqlogo plot for residues from descriptive model 
#for each substrate separating active and non active substrtates
seqlogo_df <- as.data.frame(t(msa_df_input))
colnames(seqlogo_df) <- seqlogo_df[1,]
seqlogo_df <- seqlogo_df[3:length(seqlogo_df$L26),]
#leaving only important residues with importance >0 
important_residues <- which(colnames(seqlogo_df) %in% unique(subset(msa_df, as.numeric(msa_df$Importance) != 0 )$X)) 
important_residues_labels <- colnames(seqlogo_df)[important_residues]
seqlogo_df <- seqlogo_df[, important_residues]
seq <- c()
for(i in 1:length(seqlogo_df$L26)){
  seq <- c(seq, str_replace_all(toString(seqlogo_df[i,]), ", ", ""))
}
seqlogo_df <- as.data.frame(cbind(rownames(seqlogo_df), seq)) 
seqlogo_df <- subset(seqlogo_df, !(V1 %in% c("ORAS", "MsCHS")))
#loading activity dataframe
input_df <- as.data.frame(read_sheet("https://docs.google.com/spreadsheets/d/17XtMcsRywdFOYnWKTbYyUnTUB1Qz6EFiCep7AoY1Bag/edit?usp=sharing"))
#separating training set and use function from make_df.R
input_df <- input_df[,1:16]
activity_df <- as.data.frame(c())
for( i in seqlogo_df$V1){
  activity_df <- rbind(activity_df, input_df[which(input_df$Gene == i), 6:length(input_df)])
}
seqlogo_df <- cbind(seqlogo_df, activity_df)
seqlogo_df <- melt_input_df(seqlogo_df, 2)
colnames(seqlogo_df) <- c("Gene", "seq", "Substrate", "value")
plot_names <- c()
for( i in unique(seqlogo_df$Substrate)){
  df <-subset(seqlogo_df, Substrate == i)
  name <- str_replace_all(i, "-", "")
  assign(paste(name, "_active",sep = ""), ggplotGrob(ggseqlogo(subset(df, value == 1)$seq, method = "p")+
      scale_x_discrete(breaks = 1:length(important_residues),labels = important_residues_labels, limits = factor(1:length(important_residues)))+
        labs(title = paste(i, "active (on",length(subset(df, value == 1)$seq), "sequences)",  sep = " "))+
        theme(legend.position = "none")))
  assign(paste(name, "_inactive",sep = ""), ggplotGrob(ggseqlogo(subset(df, value == 0)$seq, method = "p")+
  scale_x_discrete(breaks = 1:length(important_residues),labels = important_residues_labels, limits = factor(1:length(important_residues)))+
  labs(title = paste(i, "inactive  (on",length(subset(df, value == 0)$seq), "sequences)", sep = " "))+
    theme(legend.position = "none")))
  plot_names <- c(plot_names, paste(name, "_active",sep = ""), paste(name, "_inactive",sep = ""))
  
}

grid.arrange( coumaroylcoa_active, coumaroylcoa_inactive,benzoylcoa_active,
              benzoylcoa_inactive,hexanoylcoa_active,hexanoylcoa_inactive,
              decanoylcoa_active,decanoylcoa_inactive,myristoylcoa_active,
              myristoylcoa_inactive,oleoylcoa_active,oleoylcoa_inactive,
              phytanoylcoa_active,phytanoylcoa_inactive,NCH3anthraniloylcoa_active,
              NCH3anthraniloylcoa_inactive,phenylacetylcoa_active,phenylacetylcoa_inactive,
              acetylcoa_active,acetylcoa_inactive,methylcrotonylcoa_active,
              methylcrotonylcoa_inactive,
              ncol=2)
grid.arrange(phytanoylcoa_active,phytanoylcoa_inactive,NCH3anthraniloylcoa_active,
            NCH3anthraniloylcoa_inactive,phenylacetylcoa_active,phenylacetylcoa_inactive,
            acetylcoa_active,acetylcoa_inactive,methylcrotonylcoa_active,
            methylcrotonylcoa_inactive, ncol=2)
#making seqlogo dividing msa in half up to CadPKS
seq_half_df <-  as.data.frame(t(msa_df_input))
colnames(seq_half_df) <- seq_half_df[1,]
res_labels2 <- colnames(seq_half_df)
seq_half_df <- seq_half_df[3:length(seq_half_df$L26),]

seq <- c()
for(i in 1:length(seq_half_df$L26)){
  seq <- c(seq, str_replace_all(toString(seq_half_df[i,]), ", ", ""))
}
seq_half_df <- as.data.frame(cbind(rownames(seq_half_df), seq)) 
seq_half_df <- subset(seq_half_df, !(V1 %in% c("ORAS", "MsCHS")))

#plotting half msa
seq1 <- ggseqlogo(subset(seq_half_df, V1 %in% gene_levels[which(gene_levels == "CadPKS"):length(gene_levels)])$seq, method = "p")+
  scale_x_discrete(breaks = 1:length(res_labels2),labels = res_labels2, limits = factor(1:length(res_labels2)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.1)) + labs(title = "HypPKS - CadPKS")

seq2 <- ggseqlogo(subset(seq_half_df, V1 %in% gene_levels[1:which(gene_levels == "CadPKS")-1])$seq, method = "p")+
  scale_x_discrete(breaks = 1:length(res_labels2),labels = res_labels2, limits = factor(1:length(res_labels2)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.1)) + labs(title = "FerePKS - PhCHS")

gseq1 <- ggplotGrob(seq1)
gseq2 <- ggplotGrob(seq2)

pyrone <- c('SinsPKS','AthePKS','AastPKS',
            'AsesPKS','XacuPKS', 'HargPKS')

pyron_seq <-  ggseqlogo(subset(seq_half_df, V1 %in% pyrone)$seq, method = "p")+
  scale_x_discrete(breaks = 1:length(res_labels2),labels = res_labels2, limits = factor(1:length(res_labels2)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.1)) + labs(title = "pyron")
resorcynol_seq <-  ggseqlogo(subset(seq_half_df, !(V1 %in% pyrone))$seq, method = "p")+
  scale_x_discrete(breaks = 1:length(res_labels2),labels = res_labels2, limits = factor(1:length(res_labels2)))+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.1)) + labs(title = "resorcynol")

gpyrone <- ggplotGrob(pyron_seq)
gresorcynol <- ggplotGrob(resorcynol_seq)
grid.arrange(gpyrone, gresorcynol,nrow =2)

#arranging figure for descriptive model
lay4 <- rbind(c(1,2, 2,2,2,2),
              c(3,3,3,4,4,4),
              c(3,3,3,5,5,5))

grid.arrange(gDescriptivel_score, gViolin_plot, gmsa,gseq1, gseq2,layout_matrix = lay4)

cairo_ps(filename = "figures/figXXX_descriptive_model.eps",
         width = 50, height = 20, pointsize = 12,
         fallback_resolution = 600)
print(grid.arrange(gDescriptivel_score, gViolin_plot, gmsa,gseq1, gseq2,layout_matrix = lay4)
)
dev.off()



gExp_heatmap <- ggplotGrob(experimental)
grid.arrange(gmsa, gExp_heatmap, nrow = 1)


lay2 <- rbind(c(1,2,2,2),
              c(1,2,2,2),
              c(3,3,3,4),
              c(3,3,3,4),
              c(3,3,3,4))
grid.arrange(gDescriptivel_score, gViolin_plot, gmsa,layout_matrix = lay2)

#plotting simmilarity substarate generated from RDkit
#plotting simmilarity substarate generated from RDkit
substrates_similarity <- as.data.frame(read.csv2("output_data//substrate_simmilarity_df.csv", sep = ","))
training_subs <- colnames(input_df)[6:16]
validation_subs <- unique(substrate_heatmap_df$Substrate)
published_subs <- unique(c(unique(published_heatmap_df$Substrate), unique(published_plant_df$Substrate)))

subs_df <- as.data.frame(read.csv2('input_data/sub_space_with_valid.csv', sep = ","))
substrates_similarity <- subset(substrates_similarity, Substrate %in% subs_df$Substrate)
training_set <- c()
validation_set <- c()
published_set <- c()

for( i in substrates_similarity$Substrate){
  if(i %in% validation_subs){
    validation_set <- c(validation_set,1)
  }
  if(!(i %in% validation_subs)){
    validation_set <- c(validation_set,0)
  }
  
  
  if(i %in% published_subs){
    published_set <- c(published_set,1)
  }
  if(!(i %in% published_subs)){
    published_set <- c(published_set,0)
  }
  
  if(i %in% training_subs){
    training_set <- c(training_set,1)
  }
  if(!(i %in% training_subs)){
    training_set <- c(training_set,0)
  }
  
}
substrates_similarity$train <- training_set
substrates_similarity$valid <- validation_set
substrates_similarity$pub <- published_set

pal <- c("#005f43", "#B05670","#7fafa1", "lightblue", "#e6e6e6")

library(ggrepel)
substrates_plot <- ggplot() + geom_point(data = substrates_similarity, aes( x =as.numeric(tsne1), y = as.numeric(tsne2)), 
                      size = 12, color = "#e6e6e6")+
  geom_text_repel(data = substrates_similarity, aes(x =as.numeric(tsne1), y = as.numeric(tsne2),label = Substrate), size = 9)+
  geom_point(data = subset(substrates_similarity, train == 1), aes( x =as.numeric(tsne1), y = as.numeric(tsne2)), 
             size = 12, color = "#B05670")+
  geom_point(data = subset(substrates_similarity, valid == 1), aes( x =as.numeric(tsne1), y = as.numeric(tsne2)), 
             size = 12, color = "lightblue")+
  geom_point(data = subset(substrates_similarity, pub == 1 & train == 1), aes( x =as.numeric(tsne1), y = as.numeric(tsne2)), 
             size = 6, colour = "#005f43")+
  geom_point(data = subset(substrates_similarity, pub == 1 & valid == 1), aes( x =as.numeric(tsne1), y = as.numeric(tsne2)), 
             size = 6, colour = "#005f43")+
  geom_point(data = subset(substrates_similarity, pub == 1 & valid == 0 & train == 0 ),aes( x =as.numeric(tsne1), y = as.numeric(tsne2)), 
             size = 12, colour = "#005f43")+
   theme(legend.position  = "bottom",panel.background = element_blank(), axis.line = element_line(colour = 'black'),
        text = element_text( face = "plain", size=22))+ labs(x = "tsne1", y = "tsne2")
  

cairo_ps(filename = "figures/substrates_space.eps",
         width = 20, height = 20, pointsize = 12,
         fallback_resolution = 600)
print(substrates_plot)
dev.off()

