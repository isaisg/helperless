library(ohchibi)
library(ggvenn)
library(UpSetR)

library(extrafont)
loadfonts(device = "pdf")



setwd('/home/isai/Documents/results/helperless//scripts/')
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


### 30 minutes ###
r <- intersect(s_30m_col_ev,s_30m_nrg1_ev)

adr1_sufficient_30m <- r[ ! ( r %in% s_30m_hpk_ev ) ]

adr1_requiring_30m <- adr1_sufficient_30m[!(adr1_sufficient_30m %in% s_30m_adr1_ev)] 

nrg1_adr1_redundant_30m <- adr1_sufficient_30m[(adr1_sufficient_30m %in% s_30m_adr1_ev)] 


s <- intersect(s_30m_col_ev,s_30m_adr1_ev) 

#The contrary

nrg1_sufficient_30m <- s[!(s %in% s_30m_hpk_ev)] 


nrg1_requiring_30m <- nrg1_sufficient_30m[!(nrg1_sufficient_30m %in%s_30m_nrg1_ev)] 


helper_dependent_30m <-s_30m_col_ev[!(s_30m_col_ev %in% s_30m_hpk_ev)]

col_eti_30m =  s_30m_col_ev


### 4h  ###
r <- intersect(s_4h_col_ev,s_4h_nrg1_ev)

adr1_sufficient_4h <- r[ ! ( r %in% s_4h_hpk_ev ) ]

adr1_requiring_4h <- adr1_sufficient_4h[!(adr1_sufficient_4h %in% s_4h_adr1_ev)] 

nrg1_adr1_redundant_4h <- adr1_sufficient_4h[(adr1_sufficient_4h %in% s_4h_adr1_ev)] 


s <- intersect(s_4h_col_ev,s_4h_adr1_ev) 

#The contrary

nrg1_sufficient_4h <- s[!(s %in% s_4h_hpk_ev)] 


nrg1_requiring_4h <- nrg1_sufficient_4h[!(nrg1_sufficient_4h %in%s_4h_nrg1_ev)] 


helper_dependent_4h <-s_4h_col_ev[!(s_4h_col_ev %in% s_4h_hpk_ev)]

col_eti_4h =  s_4h_col_ev


### 8h  ###
r <- intersect(s_8h_col_ev,s_8h_nrg1_ev)

adr1_sufficient_8h <- r[ ! ( r %in% s_8h_hpk_ev ) ]

adr1_requiring_8h <- adr1_sufficient_8h[!(adr1_sufficient_8h %in% s_8h_adr1_ev)] 

nrg1_adr1_redundant_8h <- adr1_sufficient_8h[(adr1_sufficient_8h %in% s_8h_adr1_ev)] 


s <- intersect(s_8h_col_ev,s_8h_adr1_ev) 

#The contrary

nrg1_sufficient_8h <- s[!(s %in% s_8h_hpk_ev)] 


nrg1_requiring_8h <- nrg1_sufficient_8h[!(nrg1_sufficient_8h %in%s_8h_nrg1_ev)] 


helper_dependent_8h <-s_8h_col_ev[!(s_8h_col_ev %in% s_8h_hpk_ev)]

col_eti_8h =  s_8h_col_ev



adr1_requiring_30m_up <- adr1_requiring_30m
adr1_requiring_4h_up <- adr1_requiring_4h
adr1_requiring_8h_up <- adr1_requiring_8h
adr1_sufficient_30m_up <- adr1_sufficient_30m
adr1_sufficient_4h_up <- adr1_sufficient_4h
adr1_sufficient_8h_up <- adr1_sufficient_8h
col_eti_30m_up <- col_eti_30m
col_eti_4h_up <- col_eti_4h
col_eti_8h_up <- col_eti_8h
helper_dependent_30m_up <- helper_dependent_30m
helper_dependent_4h_up <- helper_dependent_4h
helper_dependent_8h_up <- helper_dependent_8h
nrg1_adr1_redundant_30m_up <- nrg1_adr1_redundant_30m
nrg1_adr1_redundant_4h_up <- nrg1_adr1_redundant_4h
nrg1_adr1_redundant_8h_up <- nrg1_adr1_redundant_8h
nrg1_requiring_30m_up <- nrg1_requiring_30m
nrg1_requiring_4h_up <- nrg1_requiring_4h
nrg1_requiring_8h_up <- nrg1_requiring_8h
nrg1_sufficient_30m_up <- nrg1_sufficient_30m
nrg1_sufficient_4h_up <- nrg1_sufficient_4h
nrg1_sufficient_8h_up <- nrg1_sufficient_8h




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

### 30 minutes ###
r <- intersect(s_30m_col_ev,s_30m_nrg1_ev)

adr1_sufficient_30m <- r[ ! ( r %in% s_30m_hpk_ev ) ]

adr1_requiring_30m <- adr1_sufficient_30m[!(adr1_sufficient_30m %in% s_30m_adr1_ev)] 

nrg1_adr1_redundant_30m <- adr1_sufficient_30m[(adr1_sufficient_30m %in% s_30m_adr1_ev)] 


s <- intersect(s_30m_col_ev,s_30m_adr1_ev) 

#The contrary

nrg1_sufficient_30m <- s[!(s %in% s_30m_hpk_ev)] 


nrg1_requiring_30m <- nrg1_sufficient_30m[!(nrg1_sufficient_30m %in%s_30m_nrg1_ev)] 


helper_dependent_30m <-s_30m_col_ev[!(s_30m_col_ev %in% s_30m_hpk_ev)]

col_eti_30m =  s_30m_col_ev


### 4h  ###
r <- intersect(s_4h_col_ev,s_4h_nrg1_ev)

adr1_sufficient_4h <- r[ ! ( r %in% s_4h_hpk_ev ) ]

adr1_requiring_4h <- adr1_sufficient_4h[!(adr1_sufficient_4h %in% s_4h_adr1_ev)] 

nrg1_adr1_redundant_4h <- adr1_sufficient_4h[(adr1_sufficient_4h %in% s_4h_adr1_ev)] 


s <- intersect(s_4h_col_ev,s_4h_adr1_ev) 

#The contrary

nrg1_sufficient_4h <- s[!(s %in% s_4h_hpk_ev)] 


nrg1_requiring_4h <- nrg1_sufficient_4h[!(nrg1_sufficient_4h %in%s_4h_nrg1_ev)] 


helper_dependent_4h <-s_4h_col_ev[!(s_4h_col_ev %in% s_4h_hpk_ev)]

col_eti_4h =  s_4h_col_ev


### 8h  ###
r <- intersect(s_8h_col_ev,s_8h_nrg1_ev)

adr1_sufficient_8h <- r[ ! ( r %in% s_8h_hpk_ev ) ]

adr1_requiring_8h <- adr1_sufficient_8h[!(adr1_sufficient_8h %in% s_8h_adr1_ev)] 

nrg1_adr1_redundant_8h <- adr1_sufficient_8h[(adr1_sufficient_8h %in% s_8h_adr1_ev)] 


s <- intersect(s_8h_col_ev,s_8h_adr1_ev) 

#The contrary

nrg1_sufficient_8h <- s[!(s %in% s_8h_hpk_ev)] 


nrg1_requiring_8h <- nrg1_sufficient_8h[!(nrg1_sufficient_8h %in%s_8h_nrg1_ev)] 


helper_dependent_8h <-s_8h_col_ev[!(s_8h_col_ev %in% s_8h_hpk_ev)]

col_eti_8h =  s_8h_col_ev


adr1_requiring_30m_down <- adr1_requiring_30m
adr1_requiring_4h_down <- adr1_requiring_4h
adr1_requiring_8h_down <- adr1_requiring_8h
adr1_sufficient_30m_down <- adr1_sufficient_30m
adr1_sufficient_4h_down <- adr1_sufficient_4h
adr1_sufficient_8h_down <- adr1_sufficient_8h
col_eti_30m_down <- col_eti_30m
col_eti_4h_down <- col_eti_4h
col_eti_8h_down <- col_eti_8h
helper_dependent_30m_down <- helper_dependent_30m
helper_dependent_4h_down <- helper_dependent_4h
helper_dependent_8h_down <- helper_dependent_8h
nrg1_adr1_redundant_30m_down <- nrg1_adr1_redundant_30m
nrg1_adr1_redundant_4h_down <- nrg1_adr1_redundant_4h
nrg1_adr1_redundant_8h_down <- nrg1_adr1_redundant_8h
nrg1_requiring_30m_down <- nrg1_requiring_30m
nrg1_requiring_4h_down <- nrg1_requiring_4h
nrg1_requiring_8h_down <- nrg1_requiring_8h
nrg1_sufficient_30m_down <- nrg1_sufficient_30m
nrg1_sufficient_4h_down <- nrg1_sufficient_4h
nrg1_sufficient_8h_down <- nrg1_sufficient_8h




################ We have all the comparison lists ################
rm(r)
rm(s)
rm(Dat)
rm(df)
rm(melted)
rm(lfc_thres)
rm(padj_thres)
rm(list = ls(pattern = "^s_"))

toremove <- ls(pattern = "^nrg1|^helper|^adr1|^col") %>% 
  grep(pattern = "up$|down$",value = T,invert = T) 
rm(list = toremove)


save.image(file = "../cleandata/helperless_contrasts_results_pti.RData")
