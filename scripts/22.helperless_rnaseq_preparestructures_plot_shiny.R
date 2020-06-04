library(ohchibi)
library(DESeq2)
setwd('/home/isai/Documents/results/helperless/scripts/')


### Prepare matrix for heatmap###
Dat <- readRDS(file = "../cleandata/dat_rnaseq.RDS")
dds <- readRDS(file = "../cleandata/dds_rnaseq.RDS")

Tab_z <- dds %>% vst %>% assay %>% t %>% scale %>% t

### Create Dat_z 
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


#Add timepoint  to each one
melted_t0 <- melted %>% subset(Treatment == "T0") %>% droplevels

new_melted <- NULL
treatments <- melted$Treatment %>% levels %>% grep(pattern = "T0",invert = T,value = T)
for(geno in unique(melted_t0$Genotype)){
  melted_temp <- melted_t0 %>% subset(Genotype == geno) %>% droplevels
  for(treat in treatments){
    mtemp <- melted_temp
    mtemp$Treatment <- treat
    new_melted <- rbind(new_melted,mtemp)
  }
}


melted <- melted %>% subset(Treatment != "T0") %>% droplevels
melted <- rbind(melted,new_melted)





melted$Genotype <- melted$Genotype %>% gsub(pattern = "Col",replacement = "Col-0") %>%
  gsub(pattern = "HPK",replacement = "helperless") %>%
  gsub(pattern = "ADR1",replacement = "adr1 triple") %>%
  gsub(pattern = "nrg1",replacement = "nrg1.1 nrg1.2") %>%
  factor(levels = c("Col-0","adr1 triple","nrg1.1 nrg1.2","helperless"))

melted$Time <- melted$Time %>% gsub(pattern = "T0",replacement = "0h") %>%
  gsub(pattern = "30m",replacement = "0.5h")  %>%
  factor(levels = c("0h","0.5h","4h","8h"))

melted$Treatment <- melted$Treatment %>%
  gsub(pattern = "S4",replacement = "avrRps4") %>%
  gsub(pattern = "T2",replacement = "avrRpt2") %>%
  gsub(pattern = "M1",replacement = "avrRpm1") %>%
  factor(levels = c("EV","avrRps4","avrRpt2","avrRpm1") )



head(melted)
melted <- melted[,c(2,3,5,6,7,8,10)]

melted$group <- paste(melted$Treatment,melted$Genotype,melted$Time,sep = " ")

order_groups <- with(melted,order(Treatment,Genotype,Time)) %>%
  melted$group[.] %>% unique

melted$group <- melted$group %>% factor(levels = order_groups)

Tab_sum <- acast(data = melted,formula = Gene~Treatment+Genotype+Time,
                 fun.aggregate = mean,value.var = "zscore")


melted_sum <- Tab_sum %>% melt
melted_sum$Treatment <- melted_sum$Var2 %>% gsub(pattern = "_.*",replacement = "") %>%
  factor(levels = c("EV","avrRps4","avrRpt2","avrRpm1") )
melted_sum$Time <- melted_sum$Var2 %>% gsub(pattern = ".*_",replacement = "") %>%
  factor(levels = c("0h","0.5h","4h","8h"))
mtemp <- melted_sum$Var2 %>%as.character %>%strsplit(split = "_") %>% unlist %>%
  matrix(data = .,ncol = 3,byrow = T)
melted_sum$Genotype <- mtemp[,2] %>% 
  factor(levels = c("Col-0","adr1 triple","nrg1.1 nrg1.2","helperless"))
melted_sum$group <- paste(melted_sum$Treatment,melted_sum$Genotype,melted_sum$Time,sep = " ")
melted_sum$group <- melted_sum$group %>% factor(levels = order_groups)
melted_sum <- melted_sum[,-2]
colnames(melted_sum)[c(1,2)] <- c("Gene","zscore")

Res <- list(
  raw = melted,
  mean = melted_sum
)
saveRDS(object =Res,file = "../cleandata/plotting_structures.RDS" )
