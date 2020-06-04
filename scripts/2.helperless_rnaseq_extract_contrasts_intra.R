library(ohchibi)
library(DESeq2)

setwd('/home/isai/Documents/results/pierre/scripts/')

Dat <- readRDS(file = "../cleandata/dat_rnaseq.RDS")
dds <- readRDS(file = "../cleandata/dds_rnaseq.RDS")

### Compute all possible contrasts ###
Dat$Map$Time <- Dat$Map$Time %>% factor %>% relevel(ref = "T0")
Dat$Map$Treatment <- Dat$Map$Treatment %>% factor %>% relevel(ref = "EV")

genotypes <- Dat$Map$Genotype %>% as.character %>% unique
treatments <- Dat$Map$Treatment %>% 
  as.character %>% unique %>% na.omit %>%
  as.character

timepoints <- c("30m","4h","8h")
reftime <- "T0"

### Loop over all genotypes to compute the results from the model ###
Res <- NULL
for(geno in genotypes){
  for(treat in treatments){
    for(time in timepoints){
      cat("Working on ",geno," ",treat," ",time,"\n")
      #Extract results from the model
      mgroup <- Dat$Map %>% 
        subset(Genotype == geno & Treatment == treat & Time == time) %$%
        group %>% as.character %>% unique
      mref <- Dat$Map %>% 
        subset(Genotype == geno & Time == reftime) %$%
        group %>% as.character %>% unique
      mcontrast <- c("group",mgroup,mref)
      cat("Contrast ",mcontrast,"\n")
      res <- results(object = dds,contrast = mcontrast,parallel = TRUE) %>% as.data.frame
      res$Gene <- rownames(res)
      res$Contrast <- paste0(mcontrast[2],"_vs_",mcontrast[3])
      res$Genotype <- geno
      res$Treatment <- treat
      res$Time <- time
      rownames(res) <- NULL
      Res <- rbind(Res,res)
    }
  }
}
