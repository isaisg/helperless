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




setwd('/home/isai/Documents/results/helperless//scripts/')
set.seed(130816)

Dat <- readRDS(file = "../cleandata/res_contrasts.RDS")
df <- Dat$df_contrasts
padj_thres <- 0.05
lfc_thres <- 1
s_8h_col_ev <- df$Contrast %>% grep(pattern = "Col.8h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character



### Compute the desired contrast
s_8h_col_s4 <- df$Contrast %>% grep(pattern = "Col.8h.S4_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_col_ev %>% length
s_8h_col_s4 %>% length

pti_col <- intersect(s_8h_col_s4,s_8h_col_ev) 




mlista <- list(
  pticol_8h = pti_col
  
)



Res <- readRDS(file = "../cleandata/plotting_structures.RDS")
#Select mean object to plot 
melted <- Res$mean

### Up ###
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
  df_stat <- NULL
  for(treat in melted_sub$Treatment %>% unique){
    df_temp <- melted_sub %>% subset(Treatment == treat) %>% droplevels
    m1 <- aov(formula = zscore ~Genotype,data = df_temp)
    m1_em <- emmeans(m1,specs = "Genotype") %>% CLD %>%
      as.data.frame
    m1_em$.group <- m1_em$.group %>% gsub(pattern = " ",replacement = "")
    colnames(m1_em)[which(colnames(m1_em) == ".group")] <- "group"
    Res_em <- m1_em
    
    #Adjust the letters of tukey
    Res_em$LLetters  <- Res_em$group %>% 
      strsplit(split = "") %>% 
      lapply(X = .,FUN = function(x)letters[as.numeric(x)])
    
    #Loop 
    Res_em$Letters <- rep(0,nrow(Res_em))
    for(i in 1:nrow(Res_em)){
      Res_em$Letters[i] <- Res_em$LLetters[[i]] %>% 
        unlist %>% toupper  %>% paste0(collapse = "")
    }
    Res_em <- Res_em[,c(1,9)]
    Res_em$Treatment <- treat
    Res_em$Time <- mtime
    df_stat <- rbind(df_stat,Res_em)
  }
  df_sum$UId <- paste(df_sum$Genotype,df_sum$Treatment,df_sum$Time,sep= "_")
  df_stat$UId <- paste(df_stat$Genotype,df_stat$Treatment,df_stat$Time,sep= "_")
  df_sum$group <- match(df_sum$UId,df_stat$UId) %>%
    df_stat$Letters[.]
  Res <- rbind(Res,df_sum)
  
}
Res$group <- Res$group %>% factor


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



Res$SetEffector <- Res$Set %>% gsub(pattern = "_.*",replacement = "") %>%
  factor()

Res$SetTime <- Res$Set %>% gsub(pattern = ".*_",replacement = "") %>%
  factor(levels = c("4h","8h"))



paleta_sig <- paletteer_d(package = "ggthemes",palette = "calc",n = 6)
names(paleta_sig) <- c("A","B","AB","C","BC","D")

#Plot the palettes
p <- ggplot(data = Res,aes(Time,zscore,fill = group)) + 
  geom_point(shape = 21,size = 6) +
  scale_fill_manual(values = paleta_sig,na.value = "#D9D9D9") +
  theme_ohchibi()

Res_end <- Res

p_up <- Res_end %>% 
  #subset(SetTime =="4h") %>%
  droplevels  %>% 
  ggplot(data = .,aes(Time,zscore)) +
  geom_smooth(method = "loess",aes(group = Genotype),color = "black",size = 0.2,se = T) +
  geom_point(aes(shape = Genotype,fill = group,color = group),color = "black",size = 8,alpha = 0.7) +
  facet_grid(SetEffector~Treatment,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 2) +
  ylab(label = "z-score") +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(size = 20,family = "Arial",face = "bold"),
    strip.text.y = element_text(size = 20,family = "Arial",face = "bold"),
    plot.title = element_text(hjust = 0.5,size = 20,family = "Arial",face = "bold")
    
  ) +
  scale_shape_manual(values = c(21,22,23,24))  +
  scale_fill_manual(values = paleta_sig,na.value = "#D9D9D9")+
  scale_color_manual(values = paleta_sig,na.value = "#D9D9D9")+
  xlab(label = "Time (hpi)")


oh.save.pdf(p = p_up,outname = "pti_8h_s4_hpk.pdf",
            outdir = "../figures/",width = 16,height = 8)


### Write table ###
colnames(Res_end)
mfinal <- Res_end[,c(1,2,3,5,9,11,12,13)]

write.table(x = mfinal,file = "../cleandata/suptable_figs3.csv",append = F,quote = F,sep = ",",col.names = T,row.names = F)
