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

#Graphic parameters
#Define colors
paleta_time <- paletteer_d(package = "rcartocolor",palette = "Safe")[1:4]

color_4h <- paleta_time[3]
color_8h <- paleta_time[4]
color_4h <- "black"
color_8h <- "black"

#Parameters for upset plots
size_titles_upset <- 30
size_text_upset <- 25
size_sets_upset <- 15



setwd('/home/isai/Documents/results/helperless//scripts/')
set.seed(130816)

#Load dataset
load(file = "../cleandata/helperless_contrasts_results.RData")

#Compute the intersections

######### UP #################

######## 4 hours #############

####### ETI genes across different effectors ####
lista <- list(
  col_eti_4h_s4 = col_eti_4h_s4_up,
  col_eti_4h_t2 = col_eti_4h_t2_up,
  col_eti_4h_m1 = col_eti_4h_m1_up
)


m <- lista %>%list_to_matrix %>%
  make_comb_mat



m_eti_4h_up <- m


##### Helper dependent genes #####
#S4#
lista <- list(
  col_eti_4h_s4 = col_eti_4h_s4_up,
  adr1_sufficient_4h_s4 = adr1_sufficient_4h_s4_up,
  adr1_requiring_4h_s4 = adr1_requiring_4h_s4_up,
  nrg1_sufficient_4h_s4 = nrg1_sufficient_4h_s4_up,
  nrg1_requiring_4h_s4 = nrg1_requiring_4h_s4_up,
  helper_dependent_4h_s4 = helper_dependent_4h_s4_up
)

m_s4 <- list_to_matrix(lista) %>%
  make_comb_mat()


#T2#
lista <- list(
  col_eti_4h_t2 = col_eti_4h_t2_up,
  adr1_sufficient_4h_t2 = adr1_sufficient_4h_t2_up,
  adr1_requiring_4h_t2 = adr1_requiring_4h_t2_up,
  nrg1_sufficient_4h_t2 = nrg1_sufficient_4h_t2_up,
  nrg1_requiring_4h_t2 = nrg1_requiring_4h_t2_up,
  helper_dependent_4h_t2 = helper_dependent_4h_t2_up
)

m_t2 <- list_to_matrix(lista) %>%
  make_comb_mat()



#M1#
lista <- list(
  col_eti_4h_m1 = col_eti_4h_m1_up,
  adr1_sufficient_4h_m1 = adr1_sufficient_4h_m1_up,
  adr1_requiring_4h_m1 = adr1_requiring_4h_m1_up,
  nrg1_sufficient_4h_m1 = nrg1_sufficient_4h_m1_up,
  nrg1_requiring_4h_m1 = nrg1_requiring_4h_m1_up,
  helper_dependent_4h_m1 = helper_dependent_4h_m1_up
)

m_m1 <- list_to_matrix(lista) %>%
  make_comb_mat()


#Here we need to compute the gene 
#categories in list format to input it to gene ontology
noms <- comb_size(m_s4) %>% names
matrices <- c("m_s4","m_t2","m_m1")
list_helper <- list()
for(mat in matrices){
  temp <- get(mat)
  mgenes <- NULL
  for(cm in noms[-1]){
    cat("Working on",mat,"\t",cm,"\n")
    mg <- extract_comb(m = temp,comb_name = cm)
    mnom <- paste0(mat,"_",cm)
    list_helper[[mnom]] <- mg
    mgenes <- c(mgenes,extract_comb(m = temp,comb_name = cm))
  }
  mpr <- paste0(mat,"_","helperdep")
  list_helper[[mpr]] <-mgenes 
  
}

names(list_helper) <- paste0(names(list_helper),"_4h")
list_4h_up <- list_helper



######## 8 hours #############

####### ETI genes across different effectors ####
lista <- list(
  col_eti_8h_s4 = col_eti_8h_s4_up,
  col_eti_8h_t2 = col_eti_8h_t2_up,
  col_eti_8h_m1 = col_eti_8h_m1_up
)


m <- lista %>%list_to_matrix %>%
  make_comb_mat


m_eti_8h_up <- m

##### Helper dependent genes #####

#S4#
lista <- list(
  col_eti_8h_s4 = col_eti_8h_s4_up,
  adr1_sufficient_8h_s4 = adr1_sufficient_8h_s4_up,
  adr1_requiring_8h_s4 = adr1_requiring_8h_s4_up,
  nrg1_sufficient_8h_s4 = nrg1_sufficient_8h_s4_up,
  nrg1_requiring_8h_s4 = nrg1_requiring_8h_s4_up,
  helper_dependent_8h_s4 = helper_dependent_8h_s4_up
)

m_s4 <- list_to_matrix(lista) %>%
  make_comb_mat()


#T2#
lista <- list(
  col_eti_8h_t2 = col_eti_8h_t2_up,
  adr1_sufficient_8h_t2 = adr1_sufficient_8h_t2_up,
  adr1_requiring_8h_t2 = adr1_requiring_8h_t2_up,
  nrg1_sufficient_8h_t2 = nrg1_sufficient_8h_t2_up,
  nrg1_requiring_8h_t2 = nrg1_requiring_8h_t2_up,
  helper_dependent_8h_t2 = helper_dependent_8h_t2_up
)

m_t2 <- list_to_matrix(lista) %>%
  make_comb_mat()



#M1#
lista <- list(
  col_eti_8h_m1 = col_eti_8h_m1_up,
  adr1_sufficient_8h_m1 = adr1_sufficient_8h_m1_up,
  adr1_requiring_8h_m1 = adr1_requiring_8h_m1_up,
  nrg1_sufficient_8h_m1 = nrg1_sufficient_8h_m1_up,
  nrg1_requiring_8h_m1 = nrg1_requiring_8h_m1_up,
  helper_dependent_8h_m1 = helper_dependent_8h_m1_up
)

m_m1 <- list_to_matrix(lista) %>%
  make_comb_mat()


#Here we need to compute the gene 
#categories in list format to input it to gene ontology
noms <- comb_size(m_s4) %>% names
matrices <- c("m_s4","m_t2","m_m1")
list_helper <- list()
for(mat in matrices){
  temp <- get(mat)
  mgenes <- NULL
  for(cm in noms[-1]){
    cat("Working on",mat,"\t",cm,"\n")
    mg <- extract_comb(m = temp,comb_name = cm)
    mnom <- paste0(mat,"_",cm)
    list_helper[[mnom]] <- mg
    mgenes <- c(mgenes,extract_comb(m = temp,comb_name = cm))
  }
  mpr <- paste0(mat,"_","helperdep")
  list_helper[[mpr]] <-mgenes 
  
}

names(list_helper) <- paste0(names(list_helper),"_8h")
list_8h_up <- list_helper




######### Down #################

######## 4 hours #############

####### ETI genes across different effectors ####
lista <- list(
  col_eti_4h_s4 = col_eti_4h_s4_down,
  col_eti_4h_t2 = col_eti_4h_t2_down,
  col_eti_4h_m1 = col_eti_4h_m1_down
)


m <- lista %>%list_to_matrix %>%
  make_comb_mat


m_eti_4h_down <- m


##### Helper dependent genes #####
#S4#
lista <- list(
  col_eti_4h_s4 = col_eti_4h_s4_down,
  adr1_sufficient_4h_s4 = adr1_sufficient_4h_s4_down,
  adr1_requiring_4h_s4 = adr1_requiring_4h_s4_down,
  nrg1_sufficient_4h_s4 = nrg1_sufficient_4h_s4_down,
  nrg1_requiring_4h_s4 = nrg1_requiring_4h_s4_down,
  helper_dependent_4h_s4 = helper_dependent_4h_s4_down
)

m_s4 <- list_to_matrix(lista) %>%
  make_comb_mat()


#T2#
lista <- list(
  col_eti_4h_t2 = col_eti_4h_t2_down,
  adr1_sufficient_4h_t2 = adr1_sufficient_4h_t2_down,
  adr1_requiring_4h_t2 = adr1_requiring_4h_t2_down,
  nrg1_sufficient_4h_t2 = nrg1_sufficient_4h_t2_down,
  nrg1_requiring_4h_t2 = nrg1_requiring_4h_t2_down,
  helper_dependent_4h_t2 = helper_dependent_4h_t2_down
)

m_t2 <- list_to_matrix(lista) %>%
  make_comb_mat()



#M1#
lista <- list(
  col_eti_4h_m1 = col_eti_4h_m1_down,
  adr1_sufficient_4h_m1 = adr1_sufficient_4h_m1_down,
  adr1_requiring_4h_m1 = adr1_requiring_4h_m1_down,
  nrg1_sufficient_4h_m1 = nrg1_sufficient_4h_m1_down,
  nrg1_requiring_4h_m1 = nrg1_requiring_4h_m1_down,
  helper_dependent_4h_m1 = helper_dependent_4h_m1_down
)

m_m1 <- list_to_matrix(lista) %>%
  make_comb_mat()


#Here we need to compute the gene 
#categories in list format to input it to gene ontology
noms <- comb_size(m_s4) %>% names
matrices <- c("m_s4","m_t2","m_m1")
list_helper <- list()
for(mat in matrices){
  temp <- get(mat)
  mgenes <- NULL
  for(cm in noms[-1]){
    cat("Working on",mat,"\t",cm,"\n")
    mg <- extract_comb(m = temp,comb_name = cm)
    mnom <- paste0(mat,"_",cm)
    list_helper[[mnom]] <- mg
    mgenes <- c(mgenes,extract_comb(m = temp,comb_name = cm))
  }
  mpr <- paste0(mat,"_","helperdep")
  list_helper[[mpr]] <-mgenes 
  
}

names(list_helper) <- paste0(names(list_helper),"_4h")
list_4h_down <- list_helper



######## 8 hours #############

####### ETI genes across different effectors ####
lista <- list(
  col_eti_8h_s4 = col_eti_8h_s4_down,
  col_eti_8h_t2 = col_eti_8h_t2_down,
  col_eti_8h_m1 = col_eti_8h_m1_down
)


m <- lista %>%list_to_matrix %>%
  make_comb_mat


m_eti_8h_down <- m

##### Helper dependent genes #####

#S4#
lista <- list(
  col_eti_8h_s4 = col_eti_8h_s4_down,
  adr1_sufficient_8h_s4 = adr1_sufficient_8h_s4_down,
  adr1_requiring_8h_s4 = adr1_requiring_8h_s4_down,
  nrg1_sufficient_8h_s4 = nrg1_sufficient_8h_s4_down,
  nrg1_requiring_8h_s4 = nrg1_requiring_8h_s4_down,
  helper_dependent_8h_s4 = helper_dependent_8h_s4_down
)

m_s4 <- list_to_matrix(lista) %>%
  make_comb_mat()


#T2#
lista <- list(
  col_eti_8h_t2 = col_eti_8h_t2_down,
  adr1_sufficient_8h_t2 = adr1_sufficient_8h_t2_down,
  adr1_requiring_8h_t2 = adr1_requiring_8h_t2_down,
  nrg1_sufficient_8h_t2 = nrg1_sufficient_8h_t2_down,
  nrg1_requiring_8h_t2 = nrg1_requiring_8h_t2_down,
  helper_dependent_8h_t2 = helper_dependent_8h_t2_down
)

m_t2 <- list_to_matrix(lista) %>%
  make_comb_mat()



#M1#
lista <- list(
  col_eti_8h_m1 = col_eti_8h_m1_down,
  adr1_sufficient_8h_m1 = adr1_sufficient_8h_m1_down,
  adr1_requiring_8h_m1 = adr1_requiring_8h_m1_down,
  nrg1_sufficient_8h_m1 = nrg1_sufficient_8h_m1_down,
  nrg1_requiring_8h_m1 = nrg1_requiring_8h_m1_down,
  helper_dependent_8h_m1 = helper_dependent_8h_m1_down
)

m_m1 <- list_to_matrix(lista) %>%
  make_comb_mat()


#Here we need to compute the gene 
#categories in list format to input it to gene ontology
noms <- comb_size(m_s4) %>% names
matrices <- c("m_s4","m_t2","m_m1")
list_helper <- list()
for(mat in matrices){
  temp <- get(mat)
  mgenes <- NULL
  for(cm in noms[-1]){
    cat("Working on",mat,"\t",cm,"\n")
    mg <- extract_comb(m = temp,comb_name = cm)
    mnom <- paste0(mat,"_",cm)
    list_helper[[mnom]] <- mg
    mgenes <- c(mgenes,extract_comb(m = temp,comb_name = cm))
  }
  mpr <- paste0(mat,"_","helperdep")
  list_helper[[mpr]] <-mgenes 
  
}

names(list_helper) <- paste0(names(list_helper),"_8h")
list_8h_down <- list_helper



### merge lists
list_up <- c(list_4h_up,list_8h_up)
list_down <- c(list_4h_down,list_8h_down)

#Append the col eti information
lista_col_up <- list(
  col_eti_4h_s4 = col_eti_4h_s4_up,
  col_eti_4h_t2 = col_eti_4h_t2_up,
  col_eti_4h_m1 = col_eti_4h_m1_up,
  col_eti_8h_s4 = col_eti_8h_s4_up,
  col_eti_8h_t2 = col_eti_8h_t2_up,
  col_eti_8h_m1 = col_eti_8h_m1_up
)



lista_col_down <- list(
  col_eti_4h_s4 = col_eti_4h_s4_down,
  col_eti_4h_t2 = col_eti_4h_t2_down,
  col_eti_4h_m1 = col_eti_4h_m1_down,
  col_eti_8h_s4 = col_eti_8h_s4_down,
  col_eti_8h_t2 = col_eti_8h_t2_down,
  col_eti_8h_m1 = col_eti_8h_m1_down
)


list_up <- c(list_up,lista_col_up)
list_down <- c(list_down,lista_col_down)

#Put both lists together
names(list_up) <- paste0(names(list_up),"_up")
names(list_down) <- paste0(names(list_down),"_down")

lista <- c(list_up,list_down)


### Load the melted structure
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

Tab_sum <- acast(data = melted,formula = Gene~group,
                 fun.aggregate = mean,value.var = "zscore")



melted <- Tab_sum %>% melt
colnames(melted) <- c("Gene","group","zscore")
melted <- merge(melted,Dat_z$Map[,c(3,4,5,7)] %>% unique,by = "group") 

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

melted$Time <- melted$Time %>% factor(levels = c("T0","30m","4h","8h"))
melted$Treatment <- melted$Treatment %>% factor(levels = c("EV","S4","T2","M1"))

paleta <- c("#3d98d3","#F99D1E","#EC0B88")
names(paleta) <- c("s4","t2","m1")
### Create unions ###
s4 <- union(lista$col_eti_4h_s4_up,lista$col_eti_8h_s4_up) 
t2 <- union(lista$col_eti_4h_t2_up,lista$col_eti_8h_t2_up) 
m1 <- union(lista$col_eti_4h_m1_up,lista$col_eti_8h_m1_up) 

mlista <- list(
  s4 = s4,
  t2 = t2,
  m1 = m1
)

### Up ###
sets <- names(mlista) 

Res <- NULL
for(set in sets){
  mgenes <- mlista[[which(names(mlista) == set)]]
  melted_sub <- melted[melted$Gene %in% mgenes,] %>%
    droplevels
  df_sum <- summarySE(data = melted_sub,measurevar = "zscore",groupvars = c("Genotype","Time","Treatment"))
  df_sum$Set <- set
  Res <- rbind(Res,df_sum)
}


#Compute average
melted_sum <- dcast(data = Res,formula = Genotype~Time+Treatment,
                    fun.aggregate = mean,value.var = "zscore") %>% melt

melted_sum$Time <- melted_sum$variable %>% gsub(pattern = "_.*",replacement = "") %>%
  factor(levels = c("T0","30m","4h","8h"))
melted_sum$Treatment <- melted_sum$variable %>% gsub(pattern = ".*_",replacement = "") %>%
  factor(levels = c("EV","S4","T2","M1"))

colnames(melted_sum)[3] <- "zscore"

Res$Set <- Res$Set %>% 
  factor(levels = c("s4","t2","m1"))


p_up <- Res %>% 
  subset(Genotype == "Col" | Genotype == "HPK") %>%
  droplevels  %>% 
  ggplot(data = .,aes(Time,zscore)) +
  geom_smooth(method = "loess",aes(group = Treatment),color = "black",size = 2,se = T) +
  geom_point(shape = 21,aes(fill = Set),size = 6,alpha = 0.9) +
  facet_grid(.~Genotype+Treatment,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 2) +
  ylab(label = "z-score") +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = 20,family = "Arial",face = "bold"),
    plot.title = element_text(hjust = 0.5,size = 20,family = "Arial",face = "bold")
    
  ) +
  scale_fill_manual(values = paleta) +
  ggtitle(label = "Upregulated")

Res_up <- Res
melted_sum_up <- melted_sum

### Down ###
### Create unions ###
s4 <- union(lista$col_eti_4h_s4_down,lista$col_eti_8h_s4_down) 
t2 <- union(lista$col_eti_4h_t2_down,lista$col_eti_8h_t2_down) 
m1 <- union(lista$col_eti_4h_m1_down,lista$col_eti_8h_m1_down) 

mlista <- list(
  s4 = s4,
  t2 = t2,
  m1 = m1
)

sets <- names(mlista) 

Res <- NULL
for(set in sets){
  mgenes <- mlista[[which(names(mlista) == set)]]
  melted_sub <- melted[melted$Gene %in% mgenes,] %>%
    droplevels
  df_sum <- summarySE(data = melted_sub,measurevar = "zscore",groupvars = c("Genotype","Time","Treatment"))
  df_sum$Set <- set
  Res <- rbind(Res,df_sum)
}



#Compute average
melted_sum <- dcast(data = Res,formula = Genotype~Time+Treatment,
                    fun.aggregate = mean,value.var = "zscore") %>% melt

melted_sum$Time <- melted_sum$variable %>% gsub(pattern = "_.*",replacement = "") %>%
  factor(levels = c("T0","30m","4h","8h"))
melted_sum$Treatment <- melted_sum$variable %>% gsub(pattern = ".*_",replacement = "") %>%
  factor(levels = c("EV","S4","T2","M1"))

colnames(melted_sum)[3] <- "zscore"

Res$Set <- Res$Set %>% 
  factor(levels = c("s4","t2","m1"))



p_down <- Res %>% 
  subset(Genotype == "Col" | Genotype == "HPK") %>%
  droplevels  %>% 
  ggplot(data = .,aes(Time,zscore)) +
  geom_smooth(method = "loess",aes(group = Treatment),color = "black",size = 2,se = T) +
  geom_point(shape = 21,aes(fill = Set),size = 6,alpha = 0.9) +
  facet_grid(.~Genotype+Treatment,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 2) +
  ylab(label = "z-score") +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_text(size = 20,family = "Arial",face = "bold"),
    plot.title = element_text(hjust = 0.5,size = 20,family = "Arial",face = "bold")
  ) +
  scale_fill_manual(values = paleta) +
  ggtitle(label = "Downregulated")

Res_down <- Res
melted_sum_down <- melted_sum


composition <- egg::ggarrange(p_up + theme(legend.position = "none"),
                              p_down + theme(legend.position = "none"),nrow = 1)
oh.save.pdf(p = composition,outname = "eti_directionality.end.pdf",outdir = "../figures/",width = 32,height = 8)

oh.save.pdf(p = p_up,outname = "eti_up.pdf",outdir = "../figures/",width = 14,height = 8)
oh.save.pdf(p = p_down,outname = "eti_down.pdf",outdir = "../figures/",width = 14,height = 8)

### Prepare table
a <- Res_up[,c(1,2,3,5,9)]
b <- Res_down[,c(1,2,3,5,9)]

a$Direction <- "Upregulated"
b$Direction <- "Downregulated"

mfinal <- rbind(a,b)
write.table(x = mfinal,file = "../cleandata/suptable_fig5bd.csv",append = F,quote = F,sep = ",",col.names = T,row.names = F)
paleta_treatment <- c(paleta_time[1],
                      paletteer_d(package = "awtools",
                                  palette = "ppalette")[c(7,3,5,2)])
