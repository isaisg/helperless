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



setwd('/home/isai/Documents/results/helperless/scripts/')
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

s_30m_col_s4 <- df$Contrast %>% grep(pattern = "Col.30m.S4_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_30m_col_t2 <- df$Contrast %>% grep(pattern = "Col.30m.T2_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_30m_col_m1 <- df$Contrast %>% grep(pattern = "Col.30m.M1_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_30m_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.30m.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_30m_adr1_s4 <- df$Contrast %>% grep(pattern = "ADR1.30m.S4_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_30m_adr1_t2 <- df$Contrast %>% grep(pattern = "ADR1.30m.T2_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_30m_adr1_m1 <- df$Contrast %>% grep(pattern = "ADR1.30m.M1_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_30m_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.30m.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_30m_hpk_s4 <- df$Contrast %>% grep(pattern = "HPK.30m.S4_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_30m_hpk_t2 <-  df$Contrast %>% grep(pattern = "HPK.30m.T2_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_30m_hpk_m1 <- df$Contrast %>% grep(pattern = "HPK.30m.M1_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character



s_30m_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.30m.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_30m_nrg1_s4 <- df$Contrast %>% grep(pattern = "nrg1.30m.S4_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_30m_nrg1_t2 <- df$Contrast %>% grep(pattern = "nrg1.30m.T2_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_30m_nrg1_m1 <- df$Contrast %>% grep(pattern = "nrg1.30m.M1_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character



### 4 hours ###
s_4h_col_ev <- df$Contrast %>% grep(pattern = "Col.4h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_4h_col_s4 <- df$Contrast %>% grep(pattern = "Col.4h.S4_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_4h_col_t2 <- df$Contrast %>% grep(pattern = "Col.4h.T2_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_4h_col_m1 <- df$Contrast %>% grep(pattern = "Col.4h.M1_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_4h_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.4h.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_4h_adr1_s4 <- df$Contrast %>% grep(pattern = "ADR1.4h.S4_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_4h_adr1_t2 <- df$Contrast %>% grep(pattern = "ADR1.4h.T2_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_4h_adr1_m1 <- df$Contrast %>% grep(pattern = "ADR1.4h.M1_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_4h_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.4h.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_4h_hpk_s4 <- df$Contrast %>% grep(pattern = "HPK.4h.S4_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_4h_hpk_t2 <-  df$Contrast %>% grep(pattern = "HPK.4h.T2_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_4h_hpk_m1 <- df$Contrast %>% grep(pattern = "HPK.4h.M1_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character



s_4h_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.4h.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_4h_nrg1_s4 <- df$Contrast %>% grep(pattern = "nrg1.4h.S4_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_4h_nrg1_t2 <- df$Contrast %>% grep(pattern = "nrg1.4h.T2_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_4h_nrg1_m1 <- df$Contrast %>% grep(pattern = "nrg1.4h.M1_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


### 8 hours ###
s_8h_col_ev <- df$Contrast %>% grep(pattern = "Col.8h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_col_s4 <- df$Contrast %>% grep(pattern = "Col.8h.S4_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_col_t2 <- df$Contrast %>% grep(pattern = "Col.8h.T2_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_col_m1 <- df$Contrast %>% grep(pattern = "Col.8h.M1_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_8h_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.8h.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_adr1_s4 <- df$Contrast %>% grep(pattern = "ADR1.8h.S4_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_adr1_t2 <- df$Contrast %>% grep(pattern = "ADR1.8h.T2_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_adr1_m1 <- df$Contrast %>% grep(pattern = "ADR1.8h.M1_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_8h_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.8h.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_hpk_s4 <- df$Contrast %>% grep(pattern = "HPK.8h.S4_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_8h_hpk_t2 <-  df$Contrast %>% grep(pattern = "HPK.8h.T2_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_8h_hpk_m1 <- df$Contrast %>% grep(pattern = "HPK.8h.M1_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character



s_8h_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.8h.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character


s_8h_nrg1_s4 <- df$Contrast %>% grep(pattern = "nrg1.8h.S4_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_nrg1_t2 <- df$Contrast %>% grep(pattern = "nrg1.8h.T2_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_nrg1_m1 <- df$Contrast %>% grep(pattern = "nrg1.8h.M1_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character



#### ETI Contrasts ####
### 30 minutes ###

col_s4_30m <- s_30m_col_s4[!(s_30m_col_s4 %in% s_30m_col_ev)]
col_t2_30m <- s_30m_col_t2[!(s_30m_col_t2 %in% s_30m_col_ev)]
col_m1_30m <- s_30m_col_m1[!(s_30m_col_m1 %in% s_30m_col_ev)]


adr1_s4_30m <- s_30m_adr1_s4[!(s_30m_adr1_s4 %in% s_30m_adr1_ev)]
adr1_t2_30m <- s_30m_adr1_t2[!(s_30m_adr1_t2 %in% s_30m_adr1_ev)]
adr1_m1_30m <- s_30m_adr1_m1[!(s_30m_adr1_m1 %in% s_30m_adr1_ev)]

hpk_s4_30m <- s_30m_hpk_s4[!(s_30m_hpk_s4 %in% s_30m_hpk_ev)]
hpk_t2_30m <- s_30m_hpk_t2[!(s_30m_hpk_t2 %in% s_30m_hpk_ev)]
hpk_m1_30m <- s_30m_hpk_m1[!(s_30m_hpk_m1 %in% s_30m_hpk_ev)]

nrg1_s4_30m <- s_30m_nrg1_s4[!(s_30m_nrg1_s4 %in% s_30m_nrg1_ev)]
nrg1_t2_30m <- s_30m_nrg1_t2[!(s_30m_nrg1_t2 %in% s_30m_nrg1_ev)]
nrg1_m1_30m <- s_30m_nrg1_m1[!(s_30m_nrg1_m1 %in% s_30m_nrg1_ev)]

### 4 hours ###

col_s4_4h <- s_4h_col_s4[!(s_4h_col_s4 %in% s_4h_col_ev)]
col_t2_4h <- s_4h_col_t2[!(s_4h_col_t2 %in% s_4h_col_ev)]
col_m1_4h <- s_4h_col_m1[!(s_4h_col_m1 %in% s_4h_col_ev)]


adr1_s4_4h <- s_4h_adr1_s4[!(s_4h_adr1_s4 %in% s_4h_adr1_ev)]
adr1_t2_4h <- s_4h_adr1_t2[!(s_4h_adr1_t2 %in% s_4h_adr1_ev)]
adr1_m1_4h <- s_4h_adr1_m1[!(s_4h_adr1_m1 %in% s_4h_adr1_ev)]

hpk_s4_4h <- s_4h_hpk_s4[!(s_4h_hpk_s4 %in% s_4h_hpk_ev)]
hpk_t2_4h <- s_4h_hpk_t2[!(s_4h_hpk_t2 %in% s_4h_hpk_ev)]
hpk_m1_4h <- s_4h_hpk_m1[!(s_4h_hpk_m1 %in% s_4h_hpk_ev)]

nrg1_s4_4h <- s_4h_nrg1_s4[!(s_4h_nrg1_s4 %in% s_4h_nrg1_ev)]
nrg1_t2_4h <- s_4h_nrg1_t2[!(s_4h_nrg1_t2 %in% s_4h_nrg1_ev)]
nrg1_m1_4h <- s_4h_nrg1_m1[!(s_4h_nrg1_m1 %in% s_4h_nrg1_ev)]


### 8 hours ###

col_s4_8h <- s_8h_col_s4[!(s_8h_col_s4 %in% s_8h_col_ev)]
col_t2_8h <- s_8h_col_t2[!(s_8h_col_t2 %in% s_8h_col_ev)]
col_m1_8h <- s_8h_col_m1[!(s_8h_col_m1 %in% s_8h_col_ev)]


adr1_s4_8h <- s_8h_adr1_s4[!(s_8h_adr1_s4 %in% s_8h_adr1_ev)]
adr1_t2_8h <- s_8h_adr1_t2[!(s_8h_adr1_t2 %in% s_8h_adr1_ev)]
adr1_m1_8h <- s_8h_adr1_m1[!(s_8h_adr1_m1 %in% s_8h_adr1_ev)]

hpk_s4_8h <- s_8h_hpk_s4[!(s_8h_hpk_s4 %in% s_8h_hpk_ev)]
hpk_t2_8h <- s_8h_hpk_t2[!(s_8h_hpk_t2 %in% s_8h_hpk_ev)]
hpk_m1_8h <- s_8h_hpk_m1[!(s_8h_hpk_m1 %in% s_8h_hpk_ev)]

nrg1_s4_8h <- s_8h_nrg1_s4[!(s_8h_nrg1_s4 %in% s_8h_nrg1_ev)]
nrg1_t2_8h <- s_8h_nrg1_t2[!(s_8h_nrg1_t2 %in% s_8h_nrg1_ev)]
nrg1_m1_8h <- s_8h_nrg1_m1[!(s_8h_nrg1_m1 %in% s_8h_nrg1_ev)]


########
mvariables <- ls(pattern = "^col|^adr1|^nrg1|^hpk")

lista_up <- list()
for(var in mvariables){
  nom <- paste0(var,"_up")
  lista_up[[nom]] <- get(var)
  
}

rm(list=mvariables)
rm(list=ls(pattern = "^s"))
rm(global_ev)
rm(var)
rm(mvariables)



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
  subset(padj < padj_thres & log2FoldChange <(lfc_thres))%$% Gene %>% as.character

s_30m_col_s4 <- df$Contrast %>% grep(pattern = "Col.30m.S4_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_30m_col_t2 <- df$Contrast %>% grep(pattern = "Col.30m.T2_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_30m_col_m1 <- df$Contrast %>% grep(pattern = "Col.30m.M1_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_30m_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.30m.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_30m_adr1_s4 <- df$Contrast %>% grep(pattern = "ADR1.30m.S4_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_30m_adr1_t2 <- df$Contrast %>% grep(pattern = "ADR1.30m.T2_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_30m_adr1_m1 <- df$Contrast %>% grep(pattern = "ADR1.30m.M1_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_30m_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.30m.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_30m_hpk_s4 <- df$Contrast %>% grep(pattern = "HPK.30m.S4_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_30m_hpk_t2 <-  df$Contrast %>% grep(pattern = "HPK.30m.T2_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_30m_hpk_m1 <- df$Contrast %>% grep(pattern = "HPK.30m.M1_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character



s_30m_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.30m.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_30m_nrg1_s4 <- df$Contrast %>% grep(pattern = "nrg1.30m.S4_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_30m_nrg1_t2 <- df$Contrast %>% grep(pattern = "nrg1.30m.T2_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_30m_nrg1_m1 <- df$Contrast %>% grep(pattern = "nrg1.30m.M1_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character



### 4 hours ###
s_4h_col_ev <- df$Contrast %>% grep(pattern = "Col.4h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_4h_col_s4 <- df$Contrast %>% grep(pattern = "Col.4h.S4_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_4h_col_t2 <- df$Contrast %>% grep(pattern = "Col.4h.T2_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_4h_col_m1 <- df$Contrast %>% grep(pattern = "Col.4h.M1_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_4h_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.4h.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_4h_adr1_s4 <- df$Contrast %>% grep(pattern = "ADR1.4h.S4_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_4h_adr1_t2 <- df$Contrast %>% grep(pattern = "ADR1.4h.T2_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_4h_adr1_m1 <- df$Contrast %>% grep(pattern = "ADR1.4h.M1_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_4h_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.4h.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_4h_hpk_s4 <- df$Contrast %>% grep(pattern = "HPK.4h.S4_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_4h_hpk_t2 <-  df$Contrast %>% grep(pattern = "HPK.4h.T2_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_4h_hpk_m1 <- df$Contrast %>% grep(pattern = "HPK.4h.M1_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character



s_4h_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.4h.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_4h_nrg1_s4 <- df$Contrast %>% grep(pattern = "nrg1.4h.S4_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_4h_nrg1_t2 <- df$Contrast %>% grep(pattern = "nrg1.4h.T2_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_4h_nrg1_m1 <- df$Contrast %>% grep(pattern = "nrg1.4h.M1_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


### 8 hours ###
s_8h_col_ev <- df$Contrast %>% grep(pattern = "Col.8h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_col_s4 <- df$Contrast %>% grep(pattern = "Col.8h.S4_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_col_t2 <- df$Contrast %>% grep(pattern = "Col.8h.T2_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_col_m1 <- df$Contrast %>% grep(pattern = "Col.8h.M1_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_8h_adr1_ev <- df$Contrast %>% grep(pattern = "ADR1.8h.EV_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_adr1_s4 <- df$Contrast %>% grep(pattern = "ADR1.8h.S4_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_adr1_t2 <- df$Contrast %>% grep(pattern = "ADR1.8h.T2_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_adr1_m1 <- df$Contrast %>% grep(pattern = "ADR1.8h.M1_vs_ADR1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_8h_hpk_ev <- df$Contrast %>% grep(pattern = "HPK.8h.EV_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_hpk_s4 <- df$Contrast %>% grep(pattern = "HPK.8h.S4_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_8h_hpk_t2 <-  df$Contrast %>% grep(pattern = "HPK.8h.T2_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_8h_hpk_m1 <- df$Contrast %>% grep(pattern = "HPK.8h.M1_vs_HPK.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character



s_8h_nrg1_ev <- df$Contrast %>% grep(pattern = "nrg1.8h.EV_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character


s_8h_nrg1_s4 <- df$Contrast %>% grep(pattern = "nrg1.8h.S4_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_nrg1_t2 <- df$Contrast %>% grep(pattern = "nrg1.8h.T2_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character

s_8h_nrg1_m1 <- df$Contrast %>% grep(pattern = "nrg1.8h.M1_vs_nrg1.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange < (lfc_thres))%$% Gene %>% as.character



#### ETI Contrasts ####
### 30 minutes ###


col_s4_30m <- s_30m_col_s4[!(s_30m_col_s4 %in% s_30m_col_ev)]
col_t2_30m <- s_30m_col_t2[!(s_30m_col_t2 %in% s_30m_col_ev)]
col_m1_30m <- s_30m_col_m1[!(s_30m_col_m1 %in% s_30m_col_ev)]


adr1_s4_30m <- s_30m_adr1_s4[!(s_30m_adr1_s4 %in% s_30m_adr1_ev)]
adr1_t2_30m <- s_30m_adr1_t2[!(s_30m_adr1_t2 %in% s_30m_adr1_ev)]
adr1_m1_30m <- s_30m_adr1_m1[!(s_30m_adr1_m1 %in% s_30m_adr1_ev)]

hpk_s4_30m <- s_30m_hpk_s4[!(s_30m_hpk_s4 %in% s_30m_hpk_ev)]
hpk_t2_30m <- s_30m_hpk_t2[!(s_30m_hpk_t2 %in% s_30m_hpk_ev)]
hpk_m1_30m <- s_30m_hpk_m1[!(s_30m_hpk_m1 %in% s_30m_hpk_ev)]

nrg1_s4_30m <- s_30m_nrg1_s4[!(s_30m_nrg1_s4 %in% s_30m_nrg1_ev)]
nrg1_t2_30m <- s_30m_nrg1_t2[!(s_30m_nrg1_t2 %in% s_30m_nrg1_ev)]
nrg1_m1_30m <- s_30m_nrg1_m1[!(s_30m_nrg1_m1 %in% s_30m_nrg1_ev)]

### 4 hours ###


col_s4_4h <- s_4h_col_s4[!(s_4h_col_s4 %in% s_4h_col_ev)]
col_t2_4h <- s_4h_col_t2[!(s_4h_col_t2 %in% s_4h_col_ev)]
col_m1_4h <- s_4h_col_m1[!(s_4h_col_m1 %in% s_4h_col_ev)]


adr1_s4_4h <- s_4h_adr1_s4[!(s_4h_adr1_s4 %in% s_4h_adr1_ev)]
adr1_t2_4h <- s_4h_adr1_t2[!(s_4h_adr1_t2 %in% s_4h_adr1_ev)]
adr1_m1_4h <- s_4h_adr1_m1[!(s_4h_adr1_m1 %in% s_4h_adr1_ev)]

hpk_s4_4h <- s_4h_hpk_s4[!(s_4h_hpk_s4 %in% s_4h_hpk_ev)]
hpk_t2_4h <- s_4h_hpk_t2[!(s_4h_hpk_t2 %in% s_4h_hpk_ev)]
hpk_m1_4h <- s_4h_hpk_m1[!(s_4h_hpk_m1 %in% s_4h_hpk_ev)]

nrg1_s4_4h <- s_4h_nrg1_s4[!(s_4h_nrg1_s4 %in% s_4h_nrg1_ev)]
nrg1_t2_4h <- s_4h_nrg1_t2[!(s_4h_nrg1_t2 %in% s_4h_nrg1_ev)]
nrg1_m1_4h <- s_4h_nrg1_m1[!(s_4h_nrg1_m1 %in% s_4h_nrg1_ev)]


### 8 hours ###


col_s4_8h <- s_8h_col_s4[!(s_8h_col_s4 %in% s_8h_col_ev)]
col_t2_8h <- s_8h_col_t2[!(s_8h_col_t2 %in% s_8h_col_ev)]
col_m1_8h <- s_8h_col_m1[!(s_8h_col_m1 %in% s_8h_col_ev)]


adr1_s4_8h <- s_8h_adr1_s4[!(s_8h_adr1_s4 %in% s_8h_adr1_ev)]
adr1_t2_8h <- s_8h_adr1_t2[!(s_8h_adr1_t2 %in% s_8h_adr1_ev)]
adr1_m1_8h <- s_8h_adr1_m1[!(s_8h_adr1_m1 %in% s_8h_adr1_ev)]

hpk_s4_8h <- s_8h_hpk_s4[!(s_8h_hpk_s4 %in% s_8h_hpk_ev)]
hpk_t2_8h <- s_8h_hpk_t2[!(s_8h_hpk_t2 %in% s_8h_hpk_ev)]
hpk_m1_8h <- s_8h_hpk_m1[!(s_8h_hpk_m1 %in% s_8h_hpk_ev)]

nrg1_s4_8h <- s_8h_nrg1_s4[!(s_8h_nrg1_s4 %in% s_8h_nrg1_ev)]
nrg1_t2_8h <- s_8h_nrg1_t2[!(s_8h_nrg1_t2 %in% s_8h_nrg1_ev)]
nrg1_m1_8h <- s_8h_nrg1_m1[!(s_8h_nrg1_m1 %in% s_8h_nrg1_ev)]


########
mvariables <- ls(pattern = "^col|^adr1|^nrg1|^hpk")

lista_down <- list()
for(var in mvariables){
  nom <- paste0(var,"_down")
  lista_down[[nom]] <- get(var)
  
}

rm(list=mvariables)
rm(list=ls(pattern = "^s"))
rm(global_ev)
rm(var)
rm(mvariables)


mlists <- list(Up = lista_up,Down = lista_down)
saveRDS(object = mlists,file = "../cleandata/list_contrasts.RDS")

rm(list=ls())
