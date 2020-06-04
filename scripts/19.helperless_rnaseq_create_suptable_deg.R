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
  for(cm in noms){
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
  for(cm in noms){
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
  for(cm in noms){
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
  for(cm in noms){
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


#Create a table with all the genes 
genes <- NULL
noms <- names(lista) %>% grep(pattern = "helperdep",invert = T,value = T)

df <- lista[[noms[1]]] %>% as.character %>% as.data.frame
df$Sig <- "Significant"
colnames(df)[1] <- "Gene"

res <- df

for(i in 2:length(noms)){
  df <- lista[[noms[i]]] %>% as.character %>% as.data.frame
  df$Sig <- "Significant"
  colnames(df)[1] <- "Gene"
  res <- merge(res,df, by = "Gene",all = T)
}

colnames(res)[2:ncol(res)] <- noms
res[(is.na(res))] <- "NotSignificant"

#Read gene description
df_desc <- read.table(file = "../cleandata/mgene_description_20131231.txt",
                      header = F,sep = "\t",fill = NA,quote = "",comment.char = "")
df_desc$Gene <- df_desc$V1 %>% gsub(pattern = "\\..*",replacement = "")

Description <- match(res$Gene,df_desc$Gene) %>%
  df_desc$V2[.] %>% as.character

res$Description <- Description


#Fix the names
colnames(res) <- colnames(res) %>% gsub(pattern = "m_",replacement = "") %>%
  gsub(pattern = "down",replacement = "Downregulated") %>%
  gsub(pattern = "up",replacement = "Upregulated") %>%
  gsub(pattern = "_s4_",replacement = "_avrRps4_") %>%
  gsub(pattern = "_t2_",replacement = "_avrRpt2_") %>%
  gsub(pattern = "_m1_",replacement = "_avrRpm1_") %>%
gsub(pattern = "^s4_",replacement = "avrRps4_") %>%
  gsub(pattern = "^t2_",replacement = "avrRpt2_") %>%
  gsub(pattern = "^m1_",replacement = "avrRpm1_") %>%
  gsub(pattern = "col_eti_",replacement = "Col-0_ETI_") %>%
  gsub(pattern = "100000",replacement = "hNLR_independent") %>%
  gsub(pattern = "100001",replacement = "ADR1+NRG1_dependent") %>%
  gsub(pattern = "111001",replacement = "ADR1_specific") %>%
  gsub(pattern = "110101",replacement = "ADR1/NRG1_redundant") %>%
  gsub(pattern = "100111",replacement = "NRG1_specific")  %>%
  gsub(pattern = "Col-0_ETI_4h_avrRpt2",replacement = "avrRpt2_Col-0_ETI_4h") %>%
  gsub(pattern = "Col-0_ETI_8h_avrRpt2",replacement = "avrRpt2_Col-0_ETI_8h") %>%
  gsub(pattern = "Col-0_ETI_4h_avrRps4",replacement = "avrRps4_Col-0_ETI_4h") %>%
  gsub(pattern = "Col-0_ETI_8h_avrRps4",replacement = "avrRps4_Col-0_ETI_8h") %>%
  gsub(pattern = "Col-0_ETI_4h_avrRpm1",replacement = "avrRpm1_Col-0_ETI_4h") %>%
  gsub(pattern = "Col-0_ETI_8h_avrRpm1",replacement = "avrRpm1_Col-0_ETI_8h") 

res <- res[,c(1,32,2:6,33,7:11,34,12:16,35,17:21,36,22:26,37,27:31,68,38:42,69,43:47,70,48:52,71,53:57,72,58:62,73,63:67,74)]

write.table(x = res,file = "../cleandata/suptable_degs.csv",append = F,quote = F,sep = ",",row.names = F,col.names = T)
