library(ohchibi)
library(DESeq2)
library(paletteer)
library(scales)
library(ggtree)
library(egg)
library(clusterProfiler)
library(org.At.tair.db)

setwd('/home/isai/Documents/results/helperless/scripts/')

Dat <- readRDS(file ="../cleandata/res_contrasts.RDS")
df <- Dat$df_contrasts 


#Adjust the levels so it keeps the same structure
df$Genotype <- df$Genotype %>% factor(levels = c("Col","ADR1","nrg1","HPK"))

df <- which(!(is.na(df$Genotype))) %>%
  df[.,] %>% droplevels

#Create structure to share
df_log <- dcast(data = df,formula = Gene~Genotype+Treatment+Time,value.var = "log2FoldChange")
df_qval <- dcast(data = df,formula = Gene~Genotype+Treatment+Time,value.var = "padj")



colnames(df_log)[2:ncol(df_log)] <- paste0("LFC_",colnames(df_log)[2:ncol(df_log)] )
colnames(df_qval)[2:ncol(df_qval)] <- paste0("qval_",colnames(df_qval)[2:ncol(df_qval)] )


all <- merge(df_log,df_qval, by = "Gene")



write.table(x = all,file = "../cleandata/res_contrasts_all.wide.csv",
            append = F,quote = F,sep = ",",row.names = F,col.names = T)
