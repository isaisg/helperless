library(ohchibi)
library(DESeq2)


#Create structure with new contrasts
setwd('/home/isai/Documents/results/helperless/scripts/')

#Load data
load(file = "../cleandata/deseq_space.first.RData")
load(file = "../cleandata/deseq_space.second.RData")
load(file = "../cleandata/deseq_space.third.RData")

#List all variables with the given profile

mdfs <- ls(pattern = "df.*")
Res <- NULL


for(d in mdfs){
  nom <- d %>% gsub(pattern = "df_",replacement = "")
  nom_vec <- nom %>% strsplit(split = "_") %>% unlist
  nom <- paste0(nom_vec[2],".",nom_vec[4],".",nom_vec[1],"_vs_",
                nom_vec[3],".",nom_vec[4],".",nom_vec[1])
  df <- d %>% get
  df$Gene <- rownames(df)
  df$Contrast <- nom
  df$Genotype <- "Inter"
  df$Treatment <- toupper(nom_vec[1])
  df$Time <- nom_vec[4]
  rownames(df) <- NULL
  Res <- rbind(Res,df)
}

Res$Contrast <-Res$Contrast %>%
  gsub(pattern = "adr1",replacement = "ADR1") %>%
  gsub(pattern = "col",replacement = "Col") %>%
  gsub(pattern = "hpk",replacement = "HPK") %>%
  gsub(pattern = "m1",replacement = "M1") %>%
  gsub(pattern = "s4",replacement = "S4") %>%
  gsub(pattern = "t2",replacement = "T2")



write.table(x = Res,file = "../cleandata/deseq2_res_contrasts.inter.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
