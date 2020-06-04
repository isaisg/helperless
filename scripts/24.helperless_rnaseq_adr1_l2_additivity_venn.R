library(ohchibi)
library(DESeq2)
library(paletteer)
library(scales)
library(ggtree)
library(egg)
library(clusterProfiler)
library(org.At.tair.db)
library(ComplexHeatmap)
library(extrafont)
loadfonts(device = "pdf")
library(Rmisc)
library(emmeans)
library(ggvenn)


#Graphic parameters


setwd('/home/isai/Documents/results/helperless//scripts/')
set.seed(130816)


### Upregulated ###
Res <- readRDS(file = "../cleandata/plotting_structures.RDS")


mup <- read.table(file = "../rawdata/ADR1-L2_DV_up.csv",header = T,sep = "\t") %$% ADR1.L2_DV_up %>%
  as.character

Dat <- readRDS(file = "../cleandata/res_contrasts.RDS")
df <- Dat$df_contrasts
padj_thres <- 0.05
lfc_thres <- 1

s_8h_col_ev <- df$Contrast %>% grep(pattern = "Col.8h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


load(file = "../cleandata/helperless_contrasts_results.RData")

col_eti_8h_s4_up  
 
#Select mean object to plot 
df <- Res$mean

# melted <- df %>% 
#   subset(Genotype == "Col-0" & (Treatment == "avrRps4" | Treatment == "EV" )) %>% droplevels

melted <- df %>% 
  subset(Genotype == "Col-0") %>% droplevels



#Create a list
s4_up <- union(col_eti_4h_s4_up,col_eti_8h_s4_up)
t2_up  <- union(col_eti_4h_t2_up,col_eti_8h_t2_up)
m1_up <- union(col_eti_4h_m1_up,col_eti_8h_m1_up)

mlista <- list(
  EV = s_8h_col_ev,
  avrRps4 = s4_up,
  avrRpt2 = t2_up,
  avrRpm1 = m1_up,
  avr1L2DV = mup
)

sets <- names(mlista) 

Res <- NULL
for(set in sets){
  mgenes <- mlista[[which(names(mlista) == set)]]
  melted_sub <- melted[melted$Gene %in% mgenes,] %>%
    droplevels
  df_sum <- summarySE(data = melted_sub,measurevar = "zscore",groupvars = c("Genotype","Time","Treatment"))
  df_sum$Set <- set
  #Now compute comparisons here
  mtime <- set %>% gsub(pattern = ".*_",replacement = "")
  melted_sub <- melted_sub %>% subset(Time == mtime) %>% droplevels
  df_sum$UId <- paste(df_sum$Genotype,df_sum$Treatment,df_sum$Time,sep= "_")

  Res <- rbind(Res,df_sum)
  
}


#Compute average
melted_sum <- dcast(data = Res,formula = Genotype~Time+Treatment,
                    fun.aggregate = mean,value.var = "zscore") %>% melt

melted_sum$Time <- melted_sum$variable %>% gsub(pattern = "_.*",replacement = "") %>%
  factor(levels = c("0h","0.5h","4h","8h"))
melted_sum$Treatment <- melted_sum$variable %>% gsub(pattern = ".*_",replacement = "") %>%
  factor(levels = c("EV","avrRps4","avrRpt2","avrRpm1"))

colnames(melted_sum)[3] <- "zscore"

Res$Set <- Res$Set %>% 
  factor(levels = sets)

paleta_set <- c("#FFDA73","#B200CC","#0EFFBC","#39B387","#FF866B")
names(paleta_set) <- c("EV","avrRps4","avr1L2DV","avrRpt2","avrRpm1")

p <- Res %>% 
  #subset(Treatment == "avrRps4") %>% droplevels %>%
  ggplot(data = .,aes(Time,zscore)) +
    geom_smooth(method = "loess",aes(group = Set),color = "black",size = 0.2,se = T) +
    geom_point(aes(fill = Set),shape = 21,size = 8,alpha = 0.8) +
    facet_grid(.~Treatment,space = "free",scales = "free") +
    theme_ohchibi(size_panel_border = 2) +
    ylab(label = "z-score") +
    theme(
      strip.background.x = element_blank(),
      strip.background.y = element_blank(),
      strip.text.x = element_text(size = 20,family = "Arial",face = "bold"),
      strip.text.y = element_text(size = 20,family = "Arial",face = "bold"),
      plot.title = element_text(hjust = 0.5,size = 20,family = "Arial",face = "bold")
      
    ) +
    xlab(label = "Time (hpi)") +
  scale_fill_manual(values = paleta_set)
oh.save.pdf(p = p,outname = "avr1ld2_additivity_up.pdf",outdir = "../figures/",width = 16,height = 8)  



# paleta <- c("#FFDA73","#B200CC","#0EFFBC")
# p <- ggvenn(data = mlista,fill_alpha = 0.6,
#                fill_color = paleta,text_size = 6,set_name_size = 10) 
# oh.save.pdf(p = p,outname = "avr1ld2_venn_up.pdf",outdir = "../figures/",width = 8,height = 8)


### Downregulated ###
rm(list=ls())
Res <- readRDS(file = "../cleandata/plotting_structures.RDS")


mdown <- read.table(file = "../rawdata/ADR1-L2_DV_down.csv",header = T,sep = "\t") %$% ADR1.L2_DV_down %>%
  as.character

Dat <- readRDS(file = "../cleandata/res_contrasts.RDS")
df <- Dat$df_contrasts
padj_thres <- 0.05
lfc_thres <- -1

s_8h_col_ev <- df$Contrast %>% grep(pattern = "Col.8h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange <(lfc_thres))%$% Gene %>% as.character


load(file = "../cleandata/helperless_contrasts_results.RData")


#Select mean object to plot 
df <- Res$mean

# melted <- df %>% 
#   subset(Genotype == "Col-0" & (Treatment == "avrRps4" | Treatment == "EV" )) %>% droplevels

melted <- df %>% 
  subset(Genotype == "Col-0" ) %>% droplevels




#Create a list
s4_down <- union(col_eti_4h_s4_down,col_eti_8h_s4_down)
t2_down <- union(col_eti_4h_t2_down,col_eti_8h_t2_down)
m1_down <- union(col_eti_4h_m1_down,col_eti_8h_m1_down)

mlista <- list(
  EV = s_8h_col_ev,
  avrRps4 = s4_down,
  avrRpt2 = t2_down,
  avrRpm1 = m1_down,
  avr1L2DV = mdown
)

sets <- names(mlista) 

Res <- NULL
for(set in sets){
  mgenes <- mlista[[which(names(mlista) == set)]]
  melted_sub <- melted[melted$Gene %in% mgenes,] %>%
    droplevels
  df_sum <- summarySE(data = melted_sub,measurevar = "zscore",groupvars = c("Genotype","Time","Treatment"))
  df_sum$Set <- set
  #Now compute comparisons here
  mtime <- set %>% gsub(pattern = ".*_",replacement = "")
  melted_sub <- melted_sub %>% subset(Time == mtime) %>% droplevels
  df_sum$UId <- paste(df_sum$Genotype,df_sum$Treatment,df_sum$Time,sep= "_")
  
  Res <- rbind(Res,df_sum)
  
}


#Compute average
melted_sum <- dcast(data = Res,formula = Genotype~Time+Treatment,
                    fun.aggregate = mean,value.var = "zscore") %>% melt

melted_sum$Time <- melted_sum$variable %>% gsub(pattern = "_.*",replacement = "") %>%
  factor(levels = c("0h","0.5h","4h","8h"))
melted_sum$Treatment <- melted_sum$variable %>% gsub(pattern = ".*_",replacement = "") %>%
  factor(levels = c("EV","avrRps4","avrRpt2","avrRpm1"))

colnames(melted_sum)[3] <- "zscore"

Res$Set <- Res$Set %>% 
  factor(levels = sets)

paleta_set <- c("#FFDA73","#B200CC","#0EFFBC","#39B387","#FF866B")
names(paleta_set) <- c("EV","avrRps4","avr1L2DV","avrRpt2","avrRpm1")


p <- Res %>% 
  #subset(Treatment == "avrRps4") %>% droplevels %>%
  ggplot(data = .,aes(Time,zscore)) +
  geom_smooth(method = "loess",aes(group = Set),color = "black",size = 0.2,se = T) +
  geom_point(aes(fill = Set),shape = 21,size = 8,alpha = 0.8) +
  facet_grid(.~Treatment,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 2) +
  ylab(label = "z-score") +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(size = 20,family = "Arial",face = "bold"),
    strip.text.y = element_text(size = 20,family = "Arial",face = "bold"),
    plot.title = element_text(hjust = 0.5,size = 20,family = "Arial",face = "bold")
    
  ) +
  xlab(label = "Time (hpi)") +
  scale_fill_manual(values =paleta_set)
oh.save.pdf(p = p,outname = "avr1ld2_additivity_down.pdf",outdir = "../figures/",width = 16,height = 8)  



# paleta <- c("#FFDA73","#B200CC","#0EFFBC")
# p <- ggvenn(data = mlista,fill_alpha = 0.6,
#             fill_color = paleta,text_size = 6,set_name_size = 10) 
# oh.save.pdf(p = p,outname = "avr1ld2_venn_down.pdf",outdir = "../figures/",width = 8,height = 8)

