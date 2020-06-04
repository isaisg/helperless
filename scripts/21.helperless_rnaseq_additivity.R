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
library(emmeans)

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


s4_4h<- lista$col_eti_4h_s4_up
t2_4h <- lista$col_eti_4h_t2_up
m1_4h <- lista$col_eti_4h_m1_up

s4_8h<- lista$col_eti_8h_s4_up
t2_8h <- lista$col_eti_8h_t2_up
m1_8h <- lista$col_eti_8h_m1_up

mlista <- list(
  s4_4h = s4_4h,
  t2_4h = t2_4h,
  m1_4h = m1_4h,
  s4_8h = s4_8h,
  t2_8h = t2_8h,
  m1_8h = m1_8h
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
  #Now compute comparisons here
  mtime <- set %>% gsub(pattern = ".*_",replacement = "")
  melted_sub <- melted_sub %>% subset(Time == mtime) %>% droplevels
  df_stat <- NULL
  for(treat in melted_sub$Treatment %>% unique){
    df_temp <- melted_sub %>% subset(Treatment == treat) %>% droplevels
    m1 <- aov(formula = zscore ~Genotype,data = df_temp)
    m1_em <- emmeans(m1,specs = "Genotype") %>% CLD %>%
      as.data.frame
    m1_em$.group <- m1_em$.group %>% gsub(pattern = " ",replacement = "")
    colnames(m1_em)[which(colnames(m1_em) == ".group")] <- "group"
    Res_em <- m1_em
    
    #Adjust the letters of tukey
    Res_em$LLetters  <- Res_em$group %>% 
      strsplit(split = "") %>% 
      lapply(X = .,FUN = function(x)letters[as.numeric(x)])
    
    #Loop 
    Res_em$Letters <- rep(0,nrow(Res_em))
    for(i in 1:nrow(Res_em)){
      Res_em$Letters[i] <- Res_em$LLetters[[i]] %>% 
        unlist %>% toupper  %>% paste0(collapse = "")
    }
    Res_em <- Res_em[,c(1,9)]
    Res_em$Treatment <- treat
    Res_em$Time <- mtime
    df_stat <- rbind(df_stat,Res_em)
  }
  df_sum$UId <- paste(df_sum$Genotype,df_sum$Treatment,df_sum$Time,sep= "_")
  df_stat$UId <- paste(df_stat$Genotype,df_stat$Treatment,df_stat$Time,sep= "_")
  df_sum$group <- match(df_sum$UId,df_stat$UId) %>%
    df_stat$Letters[.]
  Res <- rbind(Res,df_sum)
  
}
Res$group <- Res$group %>% factor

#Compute average
melted_sum <- dcast(data = Res,formula = Genotype~Time+Treatment,
                    fun.aggregate = mean,value.var = "zscore") %>% melt

melted_sum$Time <- melted_sum$variable %>% gsub(pattern = "_.*",replacement = "") %>%
  factor(levels = c("T0","30m","4h","8h"))
melted_sum$Treatment <- melted_sum$variable %>% gsub(pattern = ".*_",replacement = "") %>%
  factor(levels = c("EV","S4","T2","M1"))

colnames(melted_sum)[3] <- "zscore"

Res$Set <- Res$Set %>% 
  factor(levels = sets)

Res$SetEffector <- Res$Set %>% gsub(pattern = "_.*",replacement = "") %>%
  factor(levels = c("s4","t2","m1"))

Res$SetTime <- Res$Set %>% gsub(pattern = ".*_",replacement = "") %>%
  factor(levels = c("4h","8h"))

Res$Genotype <- Res$Genotype %>%   factor(levels = c("Col","ADR1","nrg1","HPK"))

Res$Treatment <- Res$Treatment %>% as.character %>% gsub(pattern = "S4",replacement = "avrRps4") %>%
  gsub(pattern = "T2",replacement = "avrRpt2") %>%
  gsub(pattern = "M1",replacement = "avrRpm1") %>%
  factor(levels = c("EV","avrRps4","avrRpt2","avrRpm1"))

Res$Time <- Res$Time %>% gsub(pattern = "T0",replacement = "0") %>%
  gsub(pattern = "30m",replacement = "0.5") %>% 
  gsub(pattern = "h",replacement = "") %>%
  factor

Res$Genotype <- Res$Genotype %>%
  gsub(pattern = "Col",replacement = "Col-0") %>%
  gsub(pattern = "ADR1",replacement = "adr1 triple") %>%
  gsub(pattern = "nrg1",replacement = "nrg1.1 nrg1.2") %>%
  gsub(pattern = "HPK",replacement = "helperless") %>%
  factor(levels= c("Col-0","adr1 triple","nrg1.1 nrg1.2","helperless"))


paleta_sig <- paletteer_d(package = "ggthemes",palette = "calc",n = 6)
names(paleta_sig) <- c("A","B","AB","C","BC","D")

#Plot the palettes
p <- ggplot(data = Res,aes(Time,zscore,fill = group)) + 
  geom_point(shape = 21,size = 6) +
  scale_fill_manual(values = paleta_sig,na.value = "#D9D9D9") +
  theme_ohchibi()
oh.save.pdf(p = p,outname = "legend_additivity.pdf",
            outdir = "../figures/",width = 16,height = 8)




p_up <- Res %>% 
  subset(SetTime =="4h") %>%
droplevels  %>% 
  ggplot(data = .,aes(Time,zscore)) +
  geom_smooth(method = "loess",aes(group = Genotype),color = "black",size = 0.2,se = T) +
  geom_point(aes(shape = Genotype,fill = group,color = group),color = "black",size = 6,alpha = 0.7) +
  facet_grid(SetTime~SetEffector+Treatment,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 2) +
  ylab(label = "z-score") +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(size = 20,family = "Arial",face = "bold"),
    strip.text.y = element_text(size = 20,family = "Arial",face = "bold"),
    plot.title = element_text(hjust = 0.5,size = 20,family = "Arial",face = "bold")
    
  ) +
  scale_shape_manual(values = c(21,22,23,24))  +
  scale_fill_manual(values = paleta_sig,na.value = "#D9D9D9")+
  scale_color_manual(values = paleta_sig,na.value = "#D9D9D9")+
  ggtitle(label = "Upregulated") + 
  xlab(label = "Time (hpi)")

oh.save.pdf(p = p_up,outname = "eti_additivity_upregulated_4h.pdf",
            outdir = "../figures/",width = 24,height = 8)



p_up <- Res %>% 
  subset(SetTime =="8h") %>%
  droplevels  %>% 
  ggplot(data = .,aes(Time,zscore)) +
  geom_smooth(method = "loess",aes(group = Genotype),color = "black",size = 0.2,se = T) +
  geom_point(aes(shape = Genotype,fill = group,color = group),color = "black",size = 6,alpha = 0.7) +
  facet_grid(SetTime~SetEffector+Treatment,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 2) +
  ylab(label = "z-score") +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(size = 20,family = "Arial",face = "bold"),
    strip.text.y = element_text(size = 20,family = "Arial",face = "bold"),
    plot.title = element_text(hjust = 0.5,size = 20,family = "Arial",face = "bold")
    
  ) +
  scale_shape_manual(values = c(21,22,23,24))  +
  scale_fill_manual(values = paleta_sig,na.value = "#D9D9D9")+
  scale_color_manual(values = paleta_sig,na.value = "#D9D9D9")+
  ggtitle(label = "Upregulated") + 
  xlab(label = "Time (hpi)")

oh.save.pdf(p = p_up,outname = "eti_additivity_upregulated_8h.pdf",
            outdir = "../figures/",width = 24,height = 8)

Res_up <- Res

######## Down ##########
s4_4h<- lista$col_eti_4h_s4_down
t2_4h <- lista$col_eti_4h_t2_down
m1_4h <- lista$col_eti_4h_m1_down

s4_8h<- lista$col_eti_8h_s4_down
t2_8h <- lista$col_eti_8h_t2_down
m1_8h <- lista$col_eti_8h_m1_down

mlista <- list(
  s4_4h = s4_4h,
  t2_4h = t2_4h,
  m1_4h = m1_4h,
  s4_8h = s4_8h,
  t2_8h = t2_8h,
  m1_8h = m1_8h
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
  #Now compute comparisons here
  mtime <- set %>% gsub(pattern = ".*_",replacement = "")
  melted_sub <- melted_sub %>% subset(Time == mtime) %>% droplevels
  df_stat <- NULL
  for(treat in melted_sub$Treatment %>% unique){
    df_temp <- melted_sub %>% subset(Treatment == treat) %>% droplevels
    m1 <- aov(formula = zscore ~Genotype,data = df_temp)
    m1_em <- emmeans(m1,specs = "Genotype") %>% CLD %>%
      as.data.frame
    m1_em$.group <- m1_em$.group %>% gsub(pattern = " ",replacement = "")
    colnames(m1_em)[which(colnames(m1_em) == ".group")] <- "group"
    Res_em <- m1_em
    
    #Adjust the letters of tukey
    Res_em$LLetters  <- Res_em$group %>% 
      strsplit(split = "") %>% 
      lapply(X = .,FUN = function(x)letters[as.numeric(x)])
    
    #Loop 
    Res_em$Letters <- rep(0,nrow(Res_em))
    for(i in 1:nrow(Res_em)){
      Res_em$Letters[i] <- Res_em$LLetters[[i]] %>% 
        unlist %>% toupper  %>% paste0(collapse = "")
    }
    Res_em <- Res_em[,c(1,9)]
    Res_em$Treatment <- treat
    Res_em$Time <- mtime
    df_stat <- rbind(df_stat,Res_em)
  }
  df_sum$UId <- paste(df_sum$Genotype,df_sum$Treatment,df_sum$Time,sep= "_")
  df_stat$UId <- paste(df_stat$Genotype,df_stat$Treatment,df_stat$Time,sep= "_")
  df_sum$group <- match(df_sum$UId,df_stat$UId) %>%
    df_stat$Letters[.]
  Res <- rbind(Res,df_sum)
  
}
Res$group <- Res$group %>% factor

#Compute average
melted_sum <- dcast(data = Res,formula = Genotype~Time+Treatment,
                    fun.aggregate = mean,value.var = "zscore") %>% melt

melted_sum$Time <- melted_sum$variable %>% gsub(pattern = "_.*",replacement = "") %>%
  factor(levels = c("T0","30m","4h","8h"))
melted_sum$Treatment <- melted_sum$variable %>% gsub(pattern = ".*_",replacement = "") %>%
  factor(levels = c("EV","S4","T2","M1"))

colnames(melted_sum)[3] <- "zscore"

Res$Set <- Res$Set %>% 
  factor(levels = sets)

Res$SetEffector <- Res$Set %>% gsub(pattern = "_.*",replacement = "") %>%
  factor(levels = c("s4","t2","m1"))

Res$SetTime <- Res$Set %>% gsub(pattern = ".*_",replacement = "") %>%
  factor(levels = c("4h","8h"))

Res$Genotype <- Res$Genotype %>%   factor(levels = c("Col","ADR1","nrg1","HPK"))

Res$Treatment <- Res$Treatment %>% as.character %>% gsub(pattern = "S4",replacement = "avrRps4") %>%
  gsub(pattern = "T2",replacement = "avrRpt2") %>%
  gsub(pattern = "M1",replacement = "avrRpm1") %>%
  factor(levels = c("EV","avrRps4","avrRpt2","avrRpm1"))

Res$Time <- Res$Time %>% gsub(pattern = "T0",replacement = "0") %>%
  gsub(pattern = "30m",replacement = "0.5") %>% 
  gsub(pattern = "h",replacement = "") %>%
  factor

Res$Genotype <- Res$Genotype %>%
  gsub(pattern = "Col",replacement = "Col-0") %>%
  gsub(pattern = "ADR1",replacement = "adr1 triple") %>%
  gsub(pattern = "nrg1",replacement = "nrg1.1 nrg1.2") %>%
  gsub(pattern = "HPK",replacement = "helperless") %>%
  factor(levels= c("Col-0","adr1 triple","nrg1.1 nrg1.2","helperless"))



p_down <- Res %>% 
  subset(SetTime =="4h") %>%
  droplevels  %>% 
  ggplot(data = .,aes(Time,zscore)) +
  geom_smooth(method = "loess",aes(group = Genotype),color = "black",size = 0.2,se = T) +
  geom_point(aes(shape = Genotype,fill = group),color = "black",size = 6,alpha = 0.7) +
  facet_grid(SetTime~SetEffector+Treatment,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 2) +
  ylab(label = "z-score") +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(size = 20,family = "Arial",face = "bold"),
    strip.text.y = element_text(size = 20,family = "Arial",face = "bold"),
    plot.title = element_text(hjust = 0.5,size = 20,family = "Arial",face = "bold")
    
  ) +
  scale_shape_manual(values = c(21,22,23,24))  +
  scale_fill_manual(values = paleta_sig,na.value = "#D9D9D9")+
  scale_color_manual(values = paleta_sig,na.value = "#D9D9D9")+
  ggtitle(label = "Downregulated") + 
  xlab(label = "Time (hpi)")

oh.save.pdf(p = p_down,outname = "eti_additivity_downregulated_4h.pdf",
            outdir = "../figures/",width = 24,height = 8)



p_down <- Res %>% 
  subset(SetTime =="8h") %>%
  droplevels  %>% 
  ggplot(data = .,aes(Time,zscore)) +
  geom_smooth(method = "loess",aes(group = Genotype),color = "black",size = 0.2,se = T) +
  geom_point(aes(shape = Genotype,fill = group,color = group),color = "black",size = 6,alpha = 0.7) +
  facet_grid(SetTime~SetEffector+Treatment,space = "free",scales = "free") +
  theme_ohchibi(size_panel_border = 2) +
  ylab(label = "z-score") +
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(size = 20,family = "Arial",face = "bold"),
    strip.text.y = element_text(size = 20,family = "Arial",face = "bold"),
    plot.title = element_text(hjust = 0.5,size = 20,family = "Arial",face = "bold")
    
  ) +
  scale_shape_manual(values = c(21,22,23,24))  +
  scale_fill_manual(values = paleta_sig,na.value = "#D9D9D9")+
  scale_color_manual(values = paleta_sig,na.value = "#D9D9D9")+
  ggtitle(label = "Downregulated") + 
  xlab(label = "Time (hpi)")

oh.save.pdf(p = p_down,outname = "eti_additivity_downregulated_8h.pdf",
            outdir = "../figures/",width = 24,height = 8)


Res_down <- Res

Res_up$Direction <- "Upregulated"
Res_down$Direction <- "Downregulated"

Res <- rbind(Res_up,Res_down)

mfinal <- Res[,c(1,2,3,5,9,11,12,13,14)]

write.table(x = mfinal,file = "../cleandata/suptable_figs6ab.csv",append = F,quote = F,sep = ",",col.names = T,row.names = F)


















