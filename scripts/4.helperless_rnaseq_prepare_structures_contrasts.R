library(ohchibi)
library(DESeq2)

setwd('/home/isai/Documents/results/helperless/scripts/')
set.seed(130816)

#### Create matrix for plotting #####
### Tab_sum ###
Dat <- readRDS(file = "../cleandata/dat_rnaseq.RDS")
dds <- readRDS(file = "../cleandata/dds_rnaseq.RDS")

Tab_z <- dds %>% vst %>% assay %>% t %>% scale %>% t
Dat_z <- create_dataset(Tab = Tab_z,Map = Dat$Map)
Dat_z$Map$Treatment <- Dat_z$Map$Treatment %>% as.character
Dat_z$Map$Treatment[is.na(Dat_z$Map$Treatment)] <- "T0"
Dat_z$Map$Treatment <- Dat_z$Map$Treatment %>% factor

Dat_z$Map$Treatment <- Dat_z$Map$Treatment %>% 
  factor(levels = c("T0","EV","M1","S4","T2"))
Dat_z$Map$Genotype <- Dat_z$Map$Genotype %>% 
  factor(levels = c("Col","ADR1","nrg1","HPK"))
Dat_z$Map$Time <- Dat_z$Map$Time %>% factor(levels = c("T0","30m","4h","8h"))

melted <- Dat_z$Tab %>% melt
colnames(melted) <- c("Gene","Sample_Id","zscore")

melted <- merge(melted,Dat_z$Map,by = "Sample_Id") 

Tab_sum <- acast(data = melted,formula = Gene~group,
                 fun.aggregate = mean,value.var = "zscore")

melted <- Tab_sum %>% melt
colnames(melted) <- c("Gene","group","zscore")
melted$Genotype <- melted$group %>% gsub(pattern = "\\..*",replacement = "")
melted$Treatment <- melted$group %>% gsub(pattern = ".*\\.",replacement = "")
melted$Time <- melted$group %>% gsub(pattern = "^\\S+?\\.",replacement = "") %>%
  gsub(pattern = "\\..*",replacement = "") %>% 
  factor(levels = c("T0","30m","4h","8h"))


###Read intra and inter structures derived from dds objects
df_intra <- read.table(file = "../cleandata/deseq2_res_contrasts.intra.tsv",header = T,sep = ",")
df_inter <- read.table(file = "../cleandata/deseq2_res_contrasts.inter.tsv",header = T,sep = "\t")

## Combine both structures in a dataframe where all tests are contained
df <- rbind(df_intra,df_inter)


mlist <- list(melted = melted,df_contrasts = df)
saveRDS(object = mlist,file = "../cleandata/res_contrasts.RDS")
