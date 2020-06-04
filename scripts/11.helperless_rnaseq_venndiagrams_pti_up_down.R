library(ohchibi)
library(ggvenn)
library(UpSetR)

library(extrafont)
loadfonts(device = "pdf")



setwd('/home/isai/Documents/results/pierre/scripts/')
set.seed(130816)

#Read structures
Dat <- readRDS(file = "../cleandata/res_contrasts.RDS")
df <- Dat$df_contrasts
melted <- Dat$melted



####################### Up #########################################
####################################################################
####################################################################

#Define thresholds to call enrichments
padj_thres <- 0.05
lfc_thres <- 1

### Perform all contrasts ####
### 30m hours ###
s_30m_col_ev <- df$Contrast %>% grep(pattern = "Col.30m.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_30m_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.30m.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_30m_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.30m.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_30m_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.30m.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

### 4 hours ###
s_4h_col_ev <- df$Contrast %>% grep(pattern = "Col.4h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_4h_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.4h.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_4h_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.4h.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_4h_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.4h.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

### 8 hours ###
s_8h_col_ev <- df$Contrast %>% grep(pattern = "Col.8h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.8h.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.8h.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.8h.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


### Venn diagrams ####

## Inside 30 mins ##
mlist <- list(
  col = s_30m_col_ev,
  adr1 = s_30m_adr1_ev,
  nrg1 = s_30m_nrg1_ev,
  hpk = s_30m_hpk_ev
)


p_30 <- ggvenn(data = mlist,fill_alpha = 0.6)

p_30_upset <- upset(data = fromList(mlist),set_size.numbers_size = T,text.scale = 3)


## Inside 4h  ##
mlist <- list(
  col = s_4h_col_ev,
  adr1 = s_4h_adr1_ev,
  nrg1 = s_4h_nrg1_ev,
  hpk = s_4h_hpk_ev
)


p_4h <- ggvenn(data = mlist,fill_alpha = 0.6)

p_4h_upset <- upset(data = fromList(mlist),set_size.numbers_size = T,text.scale = 3)



## Inside 8h  ##
mlist <- list(
  col = s_8h_col_ev,
  adr1 = s_8h_adr1_ev,
  nrg1 = s_8h_nrg1_ev,
  hpk = s_8h_hpk_ev
)


p_8h <- ggvenn(data = mlist,fill_alpha = 0.6)


p_8h_upset <- upset(data = fromList(mlist),set_size.numbers_size = T,text.scale = 3)


composition <- egg::ggarrange(p_30,p_4h,p_8h,nrow = 1)
oh.save.pdf(p = composition,outname = "venns_fast_pti_up.pdf",outdir = "../figures/",width = 20,height = 8)


pdf(file = "../figures/new_upset_up_30m_pti.pdf",width = 12,height = 12)
p_30_upset
dev.off()

pdf(file = "../figures/new_upset_up_4h_pti.pdf",width = 12,height = 12)
p_4h_upset
dev.off()


pdf(file = "../figures/new_upset_up_8h_pti.pdf",width = 12,height = 12)
p_8h_upset
dev.off()





####################### Down #######################################
####################################################################
####################################################################

#Define thresholds to call enrichments
padj_thres <- 0.05
lfc_thres <- -1

### Perform all contrasts ####
### 30m hours ###
s_30m_col_ev <- df$Contrast %>% grep(pattern = "Col.30m.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_30m_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.30m.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_30m_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.30m.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_30m_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.30m.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

### 4 hours ###
s_4h_col_ev <- df$Contrast %>% grep(pattern = "Col.4h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_4h_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.4h.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_4h_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.4h.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_4h_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.4h.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

### 8 hours ###
s_8h_col_ev <- df$Contrast %>% grep(pattern = "Col.8h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.8h.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.8h.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.8h.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


### Venn diagrams ####

## Inside 30 mins ##
mlist <- list(
  col = s_30m_col_ev,
  adr1 = s_30m_adr1_ev,
  nrg1 = s_30m_nrg1_ev,
  hpk = s_30m_hpk_ev
)


p_30 <- ggvenn(data = mlist,fill_alpha = 0.6)
p_30_upset <- upset(data = fromList(mlist),set_size.numbers_size = T,text.scale = 3)


## Inside 4h  ##
mlist <- list(
  col = s_4h_col_ev,
  adr1 = s_4h_adr1_ev,
  nrg1 = s_4h_nrg1_ev,
  hpk = s_4h_hpk_ev
)


p_4h <- ggvenn(data = mlist,fill_alpha = 0.6)
p_4h_upset <- upset(data = fromList(mlist),set_size.numbers_size = T,text.scale = 3)


## Inside 8h  ##
mlist <- list(
  col = s_8h_col_ev,
  adr1 = s_8h_adr1_ev,
  nrg1 = s_8h_nrg1_ev,
  hpk = s_8h_hpk_ev
)


p_8h <- ggvenn(data = mlist,fill_alpha = 0.6)

p_8h_upset <- upset(data = fromList(mlist),set_size.numbers_size = T,text.scale = 3)


composition <- egg::ggarrange(p_30,p_4h,p_8h,nrow = 1)
oh.save.pdf(p = composition,outname = "venns_fast_pti_down.pdf",outdir = "../figures/",width = 20,height = 8)


pdf(file = "../figures/new_upset_down_30m_pti.pdf",width = 12,height = 12)
p_30_upset
dev.off()

pdf(file = "../figures/new_upset_down_4h_pti.pdf",width = 12,height = 12)
p_4h_upset
dev.off()


pdf(file = "../figures/new_upset_down_8h_pti.pdf",width = 12,height = 12)
p_8h_upset
dev.off()

