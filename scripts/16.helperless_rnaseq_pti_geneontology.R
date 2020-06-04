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
load(file = "../cleandata/helperless_contrasts_results_pti.RData")

#Compute the intersections

######### UP #################

##### Helper dependent genes #####


#30minutes#
lista <- list(
  col_eti_30m = col_eti_30m_up,
  adr1_sufficient_30m = adr1_sufficient_30m_up,
  adr1_requiring_30m = adr1_requiring_30m_up,
  nrg1_sufficient_30m = nrg1_sufficient_30m_up,
  nrg1_requiring_30m = nrg1_requiring_30m_up,
  helper_dependent_30m = helper_dependent_30m_up
)


m_30 <- list_to_matrix(lista) %>%
  make_comb_mat()


#4 hours#
lista <- list(
  col_eti_4h = col_eti_4h_up,
  adr1_sufficient_4h = adr1_sufficient_4h_up,
  adr1_requiring_4h = adr1_requiring_4h_up,
  nrg1_sufficient_4h = nrg1_sufficient_4h_up,
  nrg1_requiring_4h = nrg1_requiring_4h_up,
  helper_dependent_4h = helper_dependent_4h_up
)


#Try complext heatmap implementation
#Solution here  to add number of interesections
#https://support.bioconductor.org/p/118557/
m_4h <- list_to_matrix(lista) %>%
  make_comb_mat()


#8 hours#
lista <- list(
  col_eti_8h = col_eti_8h_up,
  adr1_sufficient_8h = adr1_sufficient_8h_up,
  adr1_requiring_8h = adr1_requiring_8h_up,
  nrg1_sufficient_8h = nrg1_sufficient_8h_up,
  nrg1_requiring_8h = nrg1_requiring_8h_up,
  helper_dependent_8h = helper_dependent_8h_up
)


#Try complext heatmap implementation
#Solution here  to add number of interesections
#https://support.bioconductor.org/p/118557/
m_8h <- list_to_matrix(lista) %>%
  make_comb_mat()


#Here we need to compute the gene 
#categories in list format to input it to gene ontology
noms <- comb_size(m_30) %>% names
matrices <- c("m_30","m_4h","m_8h")
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

names(list_helper) <- paste0(names(list_helper),"_up")
list_up <- list_helper


######### Down #################


#30minutes#
lista <- list(
  col_eti_30m = col_eti_30m_down,
  adr1_sufficient_30m = adr1_sufficient_30m_down,
  adr1_requiring_30m = adr1_requiring_30m_down,
  nrg1_sufficient_30m = nrg1_sufficient_30m_down,
  nrg1_requiring_30m = nrg1_requiring_30m_down,
  helper_dependent_30m = helper_dependent_30m_down
)


#Try complext heatmap implementation
#Solution here  to add number of interesections
#https://sdownport.bioconductor.org/p/118557/
m_30 <- list_to_matrix(lista) %>%
  make_comb_mat()


#4 hours#
lista <- list(
  col_eti_4h = col_eti_4h_down,
  adr1_sufficient_4h = adr1_sufficient_4h_down,
  adr1_requiring_4h = adr1_requiring_4h_down,
  nrg1_sufficient_4h = nrg1_sufficient_4h_down,
  nrg1_requiring_4h = nrg1_requiring_4h_down,
  helper_dependent_4h = helper_dependent_4h_down
)


#Try complext heatmap implementation
#Solution here  to add number of interesections
#https://sdownport.bioconductor.org/p/118557/
m_4h <- list_to_matrix(lista) %>%
  make_comb_mat()


#8 hours#
lista <- list(
  col_eti_8h = col_eti_8h_down,
  adr1_sufficient_8h = adr1_sufficient_8h_down,
  adr1_requiring_8h = adr1_requiring_8h_down,
  nrg1_sufficient_8h = nrg1_sufficient_8h_down,
  nrg1_requiring_8h = nrg1_requiring_8h_down,
  helper_dependent_8h = helper_dependent_8h_down
)


#Try complext heatmap implementation
#Solution here  to add number of interesections
#https://sdownport.bioconductor.org/p/118557/
m_8h <- list_to_matrix(lista) %>%
  make_comb_mat()

#Here we need to compute the gene 
#categories in list format to input it to gene ontology
noms <- comb_size(m_30) %>% names
matrices <- c("m_30","m_4h","m_8h")
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

names(list_helper) <- paste0(names(list_helper),"_up")
list_down <- list_helper

#Add columbia genes
lista_col_up <- list(
  col_eti_30_up = col_eti_30m_up,
  col_eti_4h_up = col_eti_4h_up,
  col_eti_8h_up = col_eti_8h_up
)



lista_col_down <- list(
  col_eti_30_down = col_eti_30m_down,
  col_eti_4h_down = col_eti_4h_down,
  col_eti_8h_down = col_eti_8h_down
)


list_up <- c(list_up,lista_col_up)
list_down <- c(list_down,lista_col_down)

lista <- c(list_up,list_down)




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
df_plot$Time <- "30m"
df_plot$Time[df_plot$Cluster %>% grep(pattern = "4h")] <- "4h"
df_plot$Time[df_plot$Cluster %>% grep(pattern = "8h")] <- "8h"


#Type is gonna be the x axis
df_plot$Type <- "Col_eti"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "helperdep")] <- "helperdep"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "100001")] <- "adr1nrg1syn"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "111001")] <- "adr1spe"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "110101")] <- "adr1nrg1red"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "100111")] <- "nrg1spe"



df_plot$Type <- df_plot$Type %>% 
  factor(levels = c("Col_eti","helperdep","adr1nrg1syn","adr1nrg1red","adr1spe","nrg1spe"))


df_plot$Direction <- "Down"
df_plot$Direction[df_plot$Cluster %>% grep(pattern = "up")] <- "Up"
df_plot$Direction <- df_plot$Direction %>% factor(levels = c("Up","Down"))

#Color NA by other thing
df_plot$p.adjust[df_plot$p.adjust > 0.05] <- NA
df_plot$Significance <- FALSE
df_plot$Significance[which(!(is.na(df_plot$p.adjust)))] <- TRUE

p <- ggplot(data = df_plot,aes(Type,Description)) +
  geom_point(data = df_plot %>% subset(Significance == FALSE),aes(color = p.adjust,size = Count)) +
  geom_point(data = df_plot %>% subset(Significance == TRUE),aes(color = p.adjust,size = Count))+ 
  facet_grid(.~Time,space = "free",scale = "free") +
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

oh.save.pdf(p = p,outname = "geneontology_pti_up.pdf",
            outdir = "../figures/",width = 16,height = 21,family = "Arial")



### Down regulated cluster profiler
## Gene ontology analysis
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
df_plot$Time <- "30m"
df_plot$Time[df_plot$Cluster %>% grep(pattern = "4h")] <- "4h"
df_plot$Time[df_plot$Cluster %>% grep(pattern = "8h")] <- "8h"


#Type is gonna be the x axis
df_plot$Type <- "Col_eti"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "helperdep")] <- "helperdep"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "100001")] <- "adr1nrg1syn"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "111001")] <- "adr1spe"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "110101")] <- "adr1nrg1red"
df_plot$Type[df_plot$Cluster  %>% grep(pattern = "100111")] <- "nrg1spe"



df_plot$Type <- df_plot$Type %>% 
  factor(levels = c("Col_eti","helperdep","adr1nrg1syn","adr1nrg1red","adr1spe","nrg1spe"))


df_plot$Direction <- "Down"
df_plot$Direction[df_plot$Cluster %>% grep(pattern = "up")] <- "Up"
df_plot$Direction <- df_plot$Direction %>% factor(levels = c("Up","Down"))

#Color NA by other thing
df_plot$p.adjust[df_plot$p.adjust > 0.05] <- NA
df_plot$Significance <- FALSE
df_plot$Significance[which(!(is.na(df_plot$p.adjust)))] <- TRUE

p <- ggplot(data = df_plot,aes(Type,Description)) +
  geom_point(data = df_plot %>% subset(Significance == FALSE),aes(color = p.adjust,size = Count)) +
  geom_point(data = df_plot %>% subset(Significance == TRUE),aes(color = p.adjust,size = Count))+ 
  facet_grid(.~Time,space = "free",scale = "free") +
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

oh.save.pdf(p = p,outname = "geneontology_pti_down.pdf",
            outdir = "../figures/",width = 16,height = 21,family = "Arial")



