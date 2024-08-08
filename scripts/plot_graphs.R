library(ggplot2)
library(gridExtra)
library(stringr)
library(dplyr)
library(scales)

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
scores_training <- as.data.frame(read.csv2("../output_data/model_scores_training.csv", sep = ","))
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

published_score <- as.data.frame(read.csv2("../output_data/published_scores.csv", sep = ","))
published_score$set <- "published"
published_score$value <- as.numeric(published_score$value)
substrate_score <- as.data.frame(read.csv2("../output_data/extra_substrates_scores.csv", sep = ","))
substrate_score$set <- "substrates"
substrate_score$value <- as.numeric(substrate_score$value)
plant_score <- as.data.frame(read.csv2("../output_data/published_plant_scores.csv", sep =","))
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
ggsave("../figures/ML_scores.png", model_score,
       width = 50, height = 20,units = "cm", dpi = 600)

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


ggsave("../figures/ML_models_training.png", 
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

ggsave("../figures/ML_fig1.png", grid.arrange(gModel,gsubstrate_exp, 
       gsubstrate_pred,gpublished_exp, gpublished_pred,  layout_matrix = lay ),
       width = 50, height = 20,units = "cm", dpi = 600)

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


published_plant_df <- as.data.frame(read.csv2("../output_data/published_plant_prediction.csv"))
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
ggsave("../figures/figSx_plant_bacteria.png", 
       grid.arrange(gPlantScore,gpredicted_plant, gpublished_plant,layout_matrix = lay3),
       width = 50, height = 20,units = "cm", dpi = 600)

#plot prediction for all enzymes and substrates
all_prediction_df <- as.data.frame(read.csv2("../output_data/all_substrate_prediction.csv"))
all_prediction_df$Gene <- factor(all_prediction_df$Gene, levels = gene_levels)

all_prediction_df$Substrate <- factor(all_prediction_df$Substrate, levels = rev(substrates_levels))
all_prediction_plot <- plot_heatmap(all_prediction_df,all_prediction_df$Gene, all_prediction_df$Substrate, 
                                     as.numeric(all_prediction_df$pred))+
                        scale_fill_gradient(low = "white", high = "#B05670", name = "Probability")+
  scale_y_discrete(labels = seq(24, 13, -1))+theme(axis.text.y = element_text(face = "bold"))

all_prediction_plot
ggsave("../figures/all_predictions.png", all_prediction_plot,
       width = 50, height = 20,units = "cm", dpi = 600)

#plots for descriptive model
descriptive_score <- as.data.frame(read.csv2("../output_data/decision_tree_score.csv"))
descriptive_score <- melt(descriptive_score, id.vars = "X")


descriptivel_score_plot<- ggplot(descriptive_score, aes(x = variable, y = as.numeric(value)))+
  geom_boxplot(fill = "#B05670")+theme(panel.background  = element_blank(), axis.line = element_line(),
                       axis.title = element_blank(), axis.text = element_text(size = 16), 
                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
                       legend.position = "right")


gDescriptivel_score <- ggplotGrob(descriptivel_score_plot)
#box plot for Decision tree model importances
importances_df <- as.data.frame(read.csv2("../output_data/importances.csv", sep = ","))
importances_df <- melt(importances_df, id.vars = "X")
importances_df$value <- as.numeric(importances_df$value)

descriptive_boxplot <-  ggplot(importances_df, aes(x = variable, y = as.numeric(value)))+
  geom_boxplot(fill = "#B05670")+theme(panel.background  = element_blank(), axis.line = element_line(),
                                       axis.title = element_blank(), axis.text = element_text(size = 40),
                                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                                       legend.position = "right")
descriptive_boxplot
gdescriptive_boxplot <- ggplotGrob(descriptive_boxplot)
ggsave("../figures/descriptive_model.eps", descriptive_boxplot,
       width = 28, height = 10,units = "cm", dpi = 600)
