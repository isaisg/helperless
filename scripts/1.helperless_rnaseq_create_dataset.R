library(ohchibi)
library(DESeq2)

setwd('/home/isai/Documents/results/helperless/scripts/')

df <- read.table(file = "../rawdata/rnaseq_pipeline_mRNA_counts.txt",
                  header = T,row.names = 1,sep = "\t",skip = 1,
                  comment.char = "",quote = "")

#Create Tab
Tab <- df[,6:ncol(df)] %>% as.matrix
colnames(Tab) <- colnames(Tab) %>%
  gsub(pattern = "...mapped.",replacement = "") %>%
  gsub(pattern ="_.*",replacement = "" )
Tab <- Tab %>% t

Map <-data.frame(Sample_Id = rownames(Tab))
Map$Id <- Map$Sample_Id %>% 
  gsub(pattern = "^[A-H][0-9]{1,2}",replacement = "")

#Go manually and correct some of the Ids
Map$Id[110] <- "26.6.1"
Map$Id[114] <- "26.6.3"
Map$Id[151] <- "30.4.1"
Map$Id[155] <- "30.4.3"
Map$Id[272] <- "26.6.2"
Map$Id[312] <- "30.4.2"


Res <- NULL
for(id in Map$Id){
  ns <- id %>% strsplit(split = "\\.") %>%
    unlist
  if(length(ns) == 4){
    res <- data.frame(Genotype = ns[1],Time = ns[2],
               Treatment = ns[3],TechRep = ns[4])
  } else{
    res <- data.frame(Genotype = ns[1],Time = ns[2],
                      Treatment = NA,TechRep = ns[3])
  }
  Res <- rbind(Res,res)
}

Map <- cbind(Map,Res)
rownames(Map) <- Map$Sample_Id
Map$group <- paste(Map$Genotype,Map$Time,Map$Treatment,sep = ".") %>%
  factor

Dat <- create_dataset(Tab = Tab %>% t,Map = Map)

Dat <- Dat %>% subset.Dataset(subset = Genotype == "Col" |
                         Genotype == "HPK" | 
                         Genotype == "nrg1" |
                         Genotype == "ADR1",drop = T,clean = T)

Dat <- clean(Dat)


#Create dds object 
dds <- DESeqDataSetFromMatrix(countData = Dat$Tab,
                       colData = Dat$Map,design = ~ TechRep + group)

#EXplore for possible outliers
Tab_original <- vst(object = dds,blind = F) %>% assay %>% t 
Tab_v <- Tab_original %>% scale

#Perform PCA
mpca <- prcomp(x = Tab_v,center = F,scale. = F)
scores <- mpca$x %>% as.data.frame
scores$Sample_Id <- rownames(scores)
scores <- merge(scores,Dat$Map,by = "Sample_Id")

#Plot general PCA
ggplot(data = scores,aes(PC1,PC2)) +
  geom_point(size = 5,aes(color = Time))
toremove <- scores %>% subset(PC2 > 200) %$% Sample_Id


ggplot(data = scores,aes(PC1,PC2)) +
  geom_point(size = 5,aes(color = Time)) +
  scale_y_continuous(limits = c(-150,150))

toremove <- scores %>% subset(Time == "4h" & PC1 < 0 & PC2 < -80) %$%
  Sample_Id %>% c(toremove,.)

toremove <- scores %>% subset(Time == "30m" & PC1 >100 & PC2 < -40) %$%
  Sample_Id %>% c(toremove,.)

toremove <- scores %>% subset(Time == "4h" & PC1 > 90 & PC2 < -50) %$%
  Sample_Id %>% c(toremove,.)

#Redo pca with the removed samples
Tab_v <- Tab_original[!(rownames(Tab_original) %in% toremove) ,]
mgenes <- which(apply(X = Tab_v,MARGIN = 2,FUN = var)  == 0) %>%
  names
Tab_v <- Tab_v[,which(!(colnames(Tab_v) %in% mgenes) )]
Tab_v <- Tab_v %>% scale

mpca <- prcomp(x = Tab_v,center = F,scale. = F)
scores <- mpca$x %>% as.data.frame
scores$Sample_Id <- rownames(scores)
scores <- merge(scores,Dat$Map,by = "Sample_Id")

#Plot general PCA
ggplot(data = scores,aes(PC1,PC2)) +
  geom_point(size = 5,aes(color = Time))

#Do intra timepoint exploration
scores %>% 
  subset(Time == "T0") %>% droplevels %>%
  ggplot(aes(PC1,PC2)) +
  geom_point(size = 5,aes(color = Genotype))

scores %>% 
  subset(Time == "T0") %>% droplevels %>%
  ggplot(aes(PC1,PC2)) +
  geom_point(size = 5,aes(color = Genotype)) +
  scale_y_continuous(limits = c(-150,-90))

#Identify that final outlier
toremove <- scores %>% 
  subset(Time == "T0" & Genotype == "ADR1" & PC2 > -90) %$%
  Sample_Id %>% c(toremove,.)

scores %>% 
  subset(Time == "30m") %>% droplevels %>%
  ggplot(aes(PC1,PC2)) +
  geom_point(size = 5,aes(color = Genotype,shape = Treatment))

scores %>% 
  subset(Time == "4h") %>% droplevels %>%
  ggplot(aes(PC1,PC2)) +
  geom_point(size = 5,aes(color = Genotype,shape = Treatment))

scores %>% 
  subset(Time == "8h") %>% droplevels %>%
  ggplot(aes(PC1,PC2)) +
  geom_point(size = 5,aes(color = Genotype,shape = Treatment))


#### Create final structure with removed outliers
Dat_clean <- remove_samples(Dat = Dat,samples = toremove)
Dat_clean <- clean(Dat_clean)
Dat_clean$Map$Rep <- paste0("Rep",Dat_clean$Map$TechRep) %>%
  factor

saveRDS(object = Dat_clean,file = "../cleandata/dat_rnaseq.RDS")


#Create dds object 
dds <- DESeqDataSetFromMatrix(countData = Dat_clean$Tab,
                              colData = Dat_clean$Map,
                              design = ~ Rep + group)


dds <- DESeq(dds)

saveRDS(object = dds,file = "../cleandata/dds_rnaseq.RDS")

