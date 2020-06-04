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


mup <- read.table(file = "../rawdata/ADR1-L2_DV_up.csv",header = T,sep = "\t") %$% ADR1.L2_DV_up %>%
  as.character



mnoms <- names(list_up) %>% grep(pattern = "s4",value = T) %>%
  grep(pattern = "8h",value = T)


list_up <- list_up[which(names(list_up) %in% mnoms)]
list_up[["ADR1_L2_DV_8h"]] <- mup



intersect(list_up$col_eti_8h_s4_up,list_up$ADR1_L2_DV_8h) %>% length

Dat <- readRDS(file = "../cleandata/res_contrasts.RDS")
df <- Dat$df_contrasts
padj_thres <- 0.05
lfc_thres <- 1
s_8h_col_ev <- df$Contrast %>% grep(pattern = "Col.8h.EV_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

s_8h_col_s4 <- df$Contrast %>% grep(pattern = "Col.8h.S4_vs_Col.T0.NA") %>%
  df[.,] %>% droplevels %>% 
  subset(padj < padj_thres & log2FoldChange > (lfc_thres))%$% Gene %>% as.character

col_s4_8h <- s_8h_col_s4[!(s_8h_col_s4 %in% s_8h_col_ev)]


intersect(list_up$ADR1_L2_DV_8h,s_8h_col_ev) %>% length()
length(intersect(list_up$ADR1_L2_DV_8h,list_up$col_eti_8h_s4_up))

### Up regulated cluster profiler
## Gene ontology analysis
cg <- compareCluster(geneCluster=list_up,
                     fun="enrichGO",
                     keyType       = "TAIR",
                     OrgDb         = org.At.tair.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 1,
                     qvalueCutoff  = 1)

df_plot_raw <- cg %>% as.data.frame


#Loop over each category
#Display only first 25 categories per cluster
chosen_description <- NULL
df_num <- NULL
chosen_id <- NULL
for(clust in df_plot_raw$Cluster %>% levels){
  df_sub <- df_plot_raw %>% subset(Cluster  == clust) %>% droplevels
  df_sub <- df_sub %>% subset(p.adjust < 0.05) %>% droplevels
  
  #Order to determine top 25
  df_sub <- with(df_sub,order(Cluster,p.adjust)) %>%
    df_sub[.,]
  
  chosen_description <- c(chosen_description,
                          df_sub$Description[1:25] %>% na.omit %>% as.character)
  
  chosen_id <- c(chosen_id,df_sub$ID[1:25] %>% na.omit %>% as.character)
  
  #I still need to determine the total number of genes mapped per cluster and it
  numgenes <- df_sub %$% geneID %>%
    as.character %>% strsplit(split = "\\/") %>%
    unlist %>%
    #unique%>%
    length
  df_num <- rbind(df_num,data.frame(Cluster = clust,NumGenes = numgenes))
}

chosen_id <- chosen_id %>% unique
chosen_description <- chosen_description %>% unique



#Subset the structure based on p.adjust cutoff

df_plot <- df_plot_raw[df_plot_raw$ID %in% chosen_id,] %>% droplevels

#Merge with total of genes mapped per cluster

df_plot <- merge(df_plot,df_num,by = "Cluster")

### Modify the figure accordingly
#Order the descriptions based on cluster  and then p.adjustment to recapitulate
#the figure constructed by dotplot
order_description <- with(df_plot,order(Cluster,p.adjust)) %>%
  df_plot[.,] %$% Description %>% as.character %>% unique

df_plot$Description <- df_plot$Description %>% factor(levels = order_description %>% rev)


df_plot$Num <- as.numeric(df_plot$GeneRatio %>% gsub(pattern = "\\/.*",replacement = ""))/
  as.numeric(df_plot$GeneRatio %>% gsub(pattern = ".*\\/",replacement = ""))

#Order the clusters based on time and effector
df_plot$Time <- "4h"
df_plot$Time[df_plot$Cluster %>% grep(pattern = "8h")] <- "8h"


df_plot$Treatment <- "S4"
df_plot$Treatment[df_plot$Cluster %>% grep(pattern = "t2")] <- "T2"
df_plot$Treatment[df_plot$Cluster %>% grep(pattern = "m1")] <- "M1"
df_plot$Treatment <- df_plot$Treatment %>% factor(levels = c("S4","T2","M1"))


#Type is gonna be the x axis
df_plot$Type <- "Col_eti"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "helperdep")] <- "helperdep"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "100001")] <- "adr1nrg1syn"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "111001")] <- "adr1spe"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "110101")] <- "adr1nrg1red"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "100111")] <- "nrg1spe"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "ADR1_L2_DV")] <- "ADR1_L2_DV"



df_plot$Type <- df_plot$Type %>% 
  factor(levels = c("Col_eti","helperdep","adr1nrg1syn","adr1nrg1red","adr1spe","nrg1spe","ADR1_L2_DV"))


df_plot$Direction <- "Down"
df_plot$Direction[df_plot$Cluster %>% grep(pattern = "up")] <- "Up"
df_plot$Direction <- df_plot$Direction %>% factor(levels = c("Up","Down"))

#Color NA by other thing
df_plot$p.adjust[df_plot$p.adjust > 0.05] <- NA
df_plot$Significance <- FALSE
df_plot$Significance[which(!(is.na(df_plot$p.adjust)))] <- TRUE

df_plot$Treatment <- df_plot$Treatment %>% 
  as.character %>%
  gsub(pattern = "S4",replacement = "avrRps4") %>%
  gsub(pattern = "T2",replacement = "avrRpt2") %>%
  gsub(pattern = "M1",replacement = "acrRpm1") %>% 
  factor

df_plot$Type <- df_plot$Type %>% as.character %>%
  gsub(pattern =  "Col_eti",replacement = "Col-0 ETI") %>%
  gsub(pattern = "helperdep",replacement = "RNL-dependent") %>%
  gsub(pattern = "adr1nrg1syn",replacement = "adr1 nrg1 codependent") %>%
  gsub(pattern = "adr1nrg1red",replacement = "adr1 nrg1 redundant") %>%
  gsub(pattern = "adr1spe",replacement = "adr1 specific") %>%
  gsub(pattern = "nrg1spe",replacement = "nrg1 specific") %>%
  factor(levels = c("Col-0 ETI","RNL-dependent","adr1 nrg1 codependent","adr1 nrg1 redundant","adr1 specific","nrg1 specific","ADR1_L2_DV"))


p <- ggplot(data = df_plot,aes(Type,Description)) +
  geom_point(data = df_plot %>% subset(Significance == FALSE),aes(color = p.adjust,size = Count)) +
  geom_point(data = df_plot %>% subset(Significance == TRUE),aes(color = p.adjust,size = Count))+ 
  facet_grid(.~Time+Treatment,space = "free",scale = "free") +
  scale_color_paletteer_c(package = "viridis",palette = "plasma",na.value = "#BFBFBF") +
  theme(
    axis.text.x = element_text(family = "Arial",face ="bold",color = "black",angle = 45,vjust = 1,hjust = 1,size = 10),
    axis.text.y = element_text(family = "Arial",face = "bold",size = 10, color = "black"),
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_line(colour = "#D9D9D9",size = 0.5)
  )

oh.save.pdf(p = p,outname = "geneontology_eti_up_ADR1_L2_DV.pdf",
            outdir = "../figures/",width = 10,height = 14,family = "Arial")




### Down regulated cluster profiler
## Gene ontology analysis
mdown<- read.table(file = "../rawdata/ADR1-L2_DV_down.csv",header = T,sep = "\t") %$% ADR1.L2_DV_down %>%
  as.character



mnoms <- names(list_down) %>% grep(pattern = "s4",value = T) %>%
  grep(pattern = "8h",value = T)


list_down <- list_down[which(names(list_down) %in% mnoms)]
list_down[["ADR1_L2_DV_8h"]] <- mdown


cg <- compareCluster(geneCluster=list_down,
                     fun="enrichGO",
                     keyType       = "TAIR",
                     OrgDb         = org.At.tair.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 1,
                     qvalueCutoff  = 1)

df_plot_raw <- cg %>% as.data.frame


#Loop over each category
#Display only first 25 categories per cluster
chosen_description <- NULL
df_num <- NULL
chosen_id <- NULL
for(clust in df_plot_raw$Cluster %>% levels){
  df_sub <- df_plot_raw %>% subset(Cluster  == clust) %>% droplevels
  df_sub <- df_sub %>% subset(p.adjust < 0.05) %>% droplevels
  
  #Order to determine top 25
  df_sub <- with(df_sub,order(Cluster,p.adjust)) %>%
    df_sub[.,]
  
  chosen_description <- c(chosen_description,
                          df_sub$Description[1:25] %>% na.omit %>% as.character)
  
  chosen_id <- c(chosen_id,df_sub$ID[1:25] %>% na.omit %>% as.character)
  
  #I still need to determine the total number of genes mapped per cluster and it
  numgenes <- df_sub %$% geneID %>%
    as.character %>% strsplit(split = "\\/") %>%
    unlist %>%
    #unique%>%
    length
  df_num <- rbind(df_num,data.frame(Cluster = clust,NumGenes = numgenes))
}

chosen_id <- chosen_id %>% unique
chosen_description <- chosen_description %>% unique



#Subset the structure based on p.adjust cutoff

df_plot <- df_plot_raw[df_plot_raw$ID %in% chosen_id,] %>% droplevels

#Merge with total of genes mapped per cluster

df_plot <- merge(df_plot,df_num,by = "Cluster")

### Modify the figure accordingly
#Order the descriptions based on cluster  and then p.adjustment to recapitulate
#the figure constructed by dotplot
order_description <- with(df_plot,order(Cluster,p.adjust)) %>%
  df_plot[.,] %$% Description %>% as.character %>% unique

df_plot$Description <- df_plot$Description %>% factor(levels = order_description %>% rev)


df_plot$Num <- as.numeric(df_plot$GeneRatio %>% gsub(pattern = "\\/.*",replacement = ""))/
  as.numeric(df_plot$GeneRatio %>% gsub(pattern = ".*\\/",replacement = ""))

#Order the clusters based on time and effector
df_plot$Time <- "4h"
df_plot$Time[df_plot$Cluster %>% grep(pattern = "8h")] <- "8h"


df_plot$Treatment <- "S4"
df_plot$Treatment[df_plot$Cluster %>% grep(pattern = "t2")] <- "T2"
df_plot$Treatment[df_plot$Cluster %>% grep(pattern = "m1")] <- "M1"
df_plot$Treatment <- df_plot$Treatment %>% factor(levels = c("S4","T2","M1"))


#Type is gonna be the x axis
df_plot$Type <- "Col_eti"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "helperdep")] <- "helperdep"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "100001")] <- "adr1nrg1syn"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "111001")] <- "adr1spe"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "110101")] <- "adr1nrg1red"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "100111")] <- "nrg1spe"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "ADR1_L2_DV")] <- "ADR1_L2_DV"



df_plot$Type <- df_plot$Type %>% 
  factor(levels = c("Col_eti","helperdep","adr1nrg1syn","adr1nrg1red","adr1spe","nrg1spe","ADR1_L2_DV"))


df_plot$Direction <- "Down"
df_plot$Direction[df_plot$Cluster %>% grep(pattern = "up")] <- "Up"
df_plot$Direction <- df_plot$Direction %>% factor(levels = c("Up","Down"))

#Color NA by other thing
df_plot$p.adjust[df_plot$p.adjust > 0.05] <- NA
df_plot$Significance <- FALSE
df_plot$Significance[which(!(is.na(df_plot$p.adjust)))] <- TRUE


df_plot$Treatment <- df_plot$Treatment %>% 
  as.character %>%
  gsub(pattern = "S4",replacement = "avrRps4") %>%
  gsub(pattern = "T2",replacement = "avrRpt2") %>%
  gsub(pattern = "M1",replacement = "acrRpm1") %>% 
  factor

df_plot$Type <- df_plot$Type %>% as.character %>%
  gsub(pattern =  "Col_eti",replacement = "Col-0 ETI") %>%
  gsub(pattern = "helperdep",replacement = "RNL-dependent") %>%
  gsub(pattern = "adr1nrg1syn",replacement = "adr1 nrg1 codependent") %>%
  gsub(pattern = "adr1nrg1red",replacement = "adr1 nrg1 redundant") %>%
  gsub(pattern = "adr1spe",replacement = "adr1 specific") %>%
  gsub(pattern = "nrg1spe",replacement = "nrg1 specific") %>%
  factor(levels = c("Col-0 ETI","RNL-dependent","adr1 nrg1 codependent","adr1 nrg1 redundant","adr1 specific","nrg1 specific","ADR1_L2_DV"))


p <- ggplot(data = df_plot,aes(Type,Description)) +
  geom_point(data = df_plot %>% subset(Significance == FALSE),aes(color = p.adjust,size = Count)) +
  geom_point(data = df_plot %>% subset(Significance == TRUE),aes(color = p.adjust,size = Count))+ 
  facet_grid(.~Time+Treatment,space = "free",scale = "free") +
  scale_color_paletteer_c(package = "viridis",palette = "plasma",na.value = "#BFBFBF") +
  theme(
    axis.text.x = element_text(family = "Arial",face ="bold",color = "black",angle = 45,vjust = 1,hjust = 1,size = 10),
    axis.text.y = element_text(family = "Arial",face = "bold",size = 10, color = "black"),
    strip.background.x = element_blank(),
    strip.text.x = element_text(family = "Arial",face = "bold",size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(),
    panel.grid.major = element_line(colour = "#D9D9D9",size = 0.5)
  )

oh.save.pdf(p = p,outname = "geneontology_eti_down_ADR1_L2_DV.pdf",
            outdir = "../figures/",width = 10,height = 14,family = "Arial")
