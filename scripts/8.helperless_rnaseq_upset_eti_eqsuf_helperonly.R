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


cs = comb_size(m)
pdf(file = "../figures/rnaseq_upset_up_4h_eti.pdf",width = 16,height = 12)

ht <- UpSet(m, 
            set_order = c("col_eti_4h_s4","col_eti_4h_t2","col_eti_4h_m1"),
            comb_order = 1:7, 
            pt_size = unit(10, "mm"), lwd = 4,
            top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs)),
                                                  height = unit(20, "cm"),
                                                  annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                                  annotation_name_rot = 90,
                                                  annotation_name_offset = unit(2,"cm"),
                                                  axis_param = list(
                                                    gp = gpar(fontsize = size_text_upset))
            ),
            left_annotation = rowAnnotation(
              "Set size" = anno_barplot(set_size(m),
                                        annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                        axis_param = list(direction = "reverse",
                                                          gp = gpar(fontsize = size_text_upset),
                                                          labels_rot = 0
                                        ),
                                        border = FALSE, 
                                        gp = gpar(fill = "black"),
                                        width = unit(10, "cm")
              ),
              annotation_name_gp  = gpar(fontsize = size_titles_upset)
            ),
            right_annotation =NULL,
            row_names_side = "right",
            row_names_gp = gpar(fontsize = size_sets_upset),
)
ht = draw(ht)
co = column_order(ht)

nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = size_text_upset), 
            just = "bottom",
            default.units = "native")
})
dev.off()

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



### Helper intersections at 4 hours ###
noms <- comb_size(m_s4) %>% names
matrices <- c("m_s4","m_t2","m_m1")
list_helper <- list()
i <- 0
for(mat in matrices){
  i <- i + 1
  temp <- get(mat)
  mgenes <- NULL
  for(cm in noms[-1]){
    cat("Working on",mat,"\t",cm,"\n")
    mgenes <- c(mgenes,extract_comb(m = temp,comb_name = cm))
  }
  list_helper[[i]] <-mgenes 
  
}

names(list_helper) <- c("s4","t2","m1")



m <- list_helper %>%list_to_matrix %>%
  make_comb_mat


cs = comb_size(m)
pdf(file = "../figures/rnaseq_upset_up_4h_helperdependent.pdf",width = 16,height = 12)

ht <- UpSet(m, 
            set_order = c("s4","t2","m1"),
            comb_order = 1:7, 
            pt_size = unit(10, "mm"), lwd = 4,
            top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs)),
                                                  height = unit(20, "cm"),
                                                  annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                                  annotation_name_rot = 90,
                                                  annotation_name_offset = unit(2,"cm"),
                                                  axis_param = list(
                                                    gp = gpar(fontsize = size_text_upset))
            ),
            left_annotation = rowAnnotation(
              "Set size" = anno_barplot(set_size(m),
                                        annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                        axis_param = list(direction = "reverse",
                                                          gp = gpar(fontsize = size_text_upset),
                                                          labels_rot = 0
                                        ),
                                        border = FALSE, 
                                        gp = gpar(fill = "black"),
                                        width = unit(10, "cm")
              ),
              annotation_name_gp  = gpar(fontsize = size_titles_upset)
            ),
            right_annotation =NULL,
            row_names_side = "right",
            row_names_gp = gpar(fontsize = size_sets_upset),
)
ht = draw(ht)
co = column_order(ht)

nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = size_text_upset), 
            just = "bottom",
            default.units = "native")
})
dev.off()

m_helper_dependent_4h_up <- m

extract_comb(m = m_eti_4h_up,comb_name = "111") %>%
  intersect(
    extract_comb(m = m_helper_dependent_4h_up,comb_name = "111")
  )



######## 8 hours #############

####### ETI genes across different effectors ####
lista <- list(
  col_eti_8h_s4 = col_eti_8h_s4_up,
  col_eti_8h_t2 = col_eti_8h_t2_up,
  col_eti_8h_m1 = col_eti_8h_m1_up
)


m <- lista %>%list_to_matrix %>%
  make_comb_mat


cs = comb_size(m)
pdf(file = "../figures/rnaseq_upset_up_8h_eti.pdf",width = 16,height = 12)

ht <- UpSet(m, 
            set_order = c("col_eti_8h_s4","col_eti_8h_t2","col_eti_8h_m1"),
            comb_order = 1:7, 
            pt_size = unit(10, "mm"), lwd = 4,
            top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs)),
                                                  height = unit(20, "cm"),
                                                  annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                                  annotation_name_rot = 90,
                                                  annotation_name_offset = unit(2,"cm"),
                                                  axis_param = list(
                                                    gp = gpar(fontsize = size_text_upset))
            ),
            left_annotation = rowAnnotation(
              "Set size" = anno_barplot(set_size(m),
                                        annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                        axis_param = list(direction = "reverse",
                                                          gp = gpar(fontsize = size_text_upset),
                                                          labels_rot = 0
                                        ),
                                        border = FALSE, 
                                        gp = gpar(fill = "black"),
                                        width = unit(10, "cm")
              ),
              annotation_name_gp  = gpar(fontsize = size_titles_upset)
            ),
            right_annotation =NULL,
            row_names_side = "right",
            row_names_gp = gpar(fontsize = size_sets_upset),
)
ht = draw(ht)
co = column_order(ht)

nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = size_text_upset), 
            just = "bottom",
            default.units = "native")
})
dev.off()
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



#Helper intersections
noms <- comb_size(m_s4) %>% names
matrices <- c("m_s4","m_t2","m_m1")
list_helper <- list()
i <- 0
for(mat in matrices){
  i <- i + 1
  temp <- get(mat)
  mgenes <- NULL
  for(cm in noms[-1]){
    cat("Working on",mat,"\t",cm,"\n")
    mgenes <- c(mgenes,extract_comb(m = temp,comb_name = cm))
  }
  list_helper[[i]] <-mgenes 
  
}

names(list_helper) <- c("s4","t2","m1")


m <- list_helper %>%list_to_matrix %>%
  make_comb_mat


cs = comb_size(m)
pdf(file = "../figures/rnaseq_upset_up_8h_helperdependent.pdf",width = 16,height = 12)

ht <- UpSet(m, 
            set_order = c("s4","t2","m1"),
            comb_order = 1:7, 
            pt_size = unit(10, "mm"), lwd = 4,
            top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs)),
                                                  height = unit(20, "cm"),
                                                  annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                                  annotation_name_rot = 90,
                                                  annotation_name_offset = unit(2,"cm"),
                                                  axis_param = list(
                                                    gp = gpar(fontsize = size_text_upset))
            ),
            left_annotation = rowAnnotation(
              "Set size" = anno_barplot(set_size(m),
                                        annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                        axis_param = list(direction = "reverse",
                                                          gp = gpar(fontsize = size_text_upset),
                                                          labels_rot = 0
                                        ),
                                        border = FALSE, 
                                        gp = gpar(fill = "black"),
                                        width = unit(10, "cm")
              ),
              annotation_name_gp  = gpar(fontsize = size_titles_upset)
            ),
            right_annotation =NULL,
            row_names_side = "right",
            row_names_gp = gpar(fontsize = size_sets_upset),
)
ht = draw(ht)
co = column_order(ht)

nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = size_text_upset), 
            just = "bottom",
            default.units = "native")
})
dev.off()

m_helper_dependent_8h_up <- m


extract_comb(m = m_eti_8h_up,comb_name = "111") %>%
  intersect(
    extract_comb(m = m_helper_dependent_8h_up,comb_name = "111")
  )



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


cs = comb_size(m)
pdf(file = "../figures/rnaseq_upset_down_4h_eti.pdf",width = 16,height = 12)

ht <- UpSet(m, 
            set_order = c("col_eti_4h_s4","col_eti_4h_t2","col_eti_4h_m1"),
            comb_order = 1:7, 
            pt_size = unit(10, "mm"), lwd = 4,
            top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs)),
                                                  height = unit(20, "cm"),
                                                  annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                                  annotation_name_rot = 90,
                                                  annotation_name_offset = unit(2,"cm"),
                                                  axis_param = list(
                                                    gp = gpar(fontsize = size_text_upset))
            ),
            left_annotation = rowAnnotation(
              "Set size" = anno_barplot(set_size(m),
                                        annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                        axis_param = list(direction = "reverse",
                                                          gp = gpar(fontsize = size_text_upset),
                                                          labels_rot = 0
                                        ),
                                        border = FALSE, 
                                        gp = gpar(fill = "black"),
                                        width = unit(10, "cm")
              ),
              annotation_name_gp  = gpar(fontsize = size_titles_upset)
            ),
            right_annotation =NULL,
            row_names_side = "right",
            row_names_gp = gpar(fontsize = size_sets_upset),
)
ht = draw(ht)
co = column_order(ht)

nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = size_text_upset), 
            just = "bottom",
            default.units = "native")
})
dev.off()
dev.off()
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



### Helper intersections at 4 hours ###
noms <- comb_size(m_s4) %>% names
matrices <- c("m_s4","m_t2","m_m1")
list_helper <- list()
i <- 0
for(mat in matrices){
  i <- i + 1
  temp <- get(mat)
  mgenes <- NULL
  for(cm in noms[-1]){
    cat("Working on",mat,"\t",cm,"\n")
    mgenes <- c(mgenes,extract_comb(m = temp,comb_name = cm))
  }
  list_helper[[i]] <-mgenes 
  
}

names(list_helper) <- c("s4","t2","m1")



m <- list_helper %>%list_to_matrix %>%
  make_comb_mat


cs = comb_size(m)
pdf(file = "../figures/rnaseq_upset_down_4h_helperdependent.pdf",width = 16,height = 12)

ht <- UpSet(m, 
            set_order = c("s4","t2","m1"),
            comb_order = 1:7, 
            pt_size = unit(10, "mm"), lwd = 4,
            top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs)),
                                                  height = unit(20, "cm"),
                                                  annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                                  annotation_name_rot = 90,
                                                  annotation_name_offset = unit(2,"cm"),
                                                  axis_param = list(
                                                    gp = gpar(fontsize = size_text_upset))
            ),
            left_annotation = rowAnnotation(
              "Set size" = anno_barplot(set_size(m),
                                        annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                        axis_param = list(direction = "reverse",
                                                          gp = gpar(fontsize = size_text_upset),
                                                          labels_rot = 0
                                        ),
                                        border = FALSE, 
                                        gp = gpar(fill = "black"),
                                        width = unit(10, "cm")
              ),
              annotation_name_gp  = gpar(fontsize = size_titles_upset)
            ),
            right_annotation =NULL,
            row_names_side = "right",
            row_names_gp = gpar(fontsize = size_sets_upset),
)
ht = draw(ht)
co = column_order(ht)

nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = size_text_upset), 
            just = "bottom",
            default.units = "native")
})
dev.off()

m_helper_dependent_4h_down <- m

extract_comb(m = m_eti_4h_down,comb_name = "111") %>%
  intersect(
    extract_comb(m = m_helper_dependent_4h_down,comb_name = "111")
  )


######## 8 hours #############

####### ETI genes across different effectors ####
lista <- list(
  col_eti_8h_s4 = col_eti_8h_s4_down,
  col_eti_8h_t2 = col_eti_8h_t2_down,
  col_eti_8h_m1 = col_eti_8h_m1_down
)


m <- lista %>%list_to_matrix %>%
  make_comb_mat


cs = comb_size(m)
pdf(file = "../figures/rnaseq_upset_down_8h_eti.pdf",width = 16,height = 12)

ht <- UpSet(m, 
            set_order = c("col_eti_8h_s4","col_eti_8h_t2","col_eti_8h_m1"),
            comb_order = 1:7, 
            pt_size = unit(10, "mm"), lwd = 4,
            top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs)),
                                                  height = unit(20, "cm"),
                                                  annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                                  annotation_name_rot = 90,
                                                  annotation_name_offset = unit(2,"cm"),
                                                  axis_param = list(
                                                    gp = gpar(fontsize = size_text_upset))
            ),
            left_annotation = rowAnnotation(
              "Set size" = anno_barplot(set_size(m),
                                        annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                        axis_param = list(direction = "reverse",
                                                          gp = gpar(fontsize = size_text_upset),
                                                          labels_rot = 0
                                        ),
                                        border = FALSE, 
                                        gp = gpar(fill = "black"),
                                        width = unit(10, "cm")
              ),
              annotation_name_gp  = gpar(fontsize = size_titles_upset)
            ),
            right_annotation =NULL,
            row_names_side = "right",
            row_names_gp = gpar(fontsize = size_sets_upset),
)
ht = draw(ht)
co = column_order(ht)

nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = size_text_upset), 
            just = "bottom",
            default.units = "native")
})
dev.off()

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



#Helper intersections
noms <- comb_size(m_s4) %>% names
matrices <- c("m_s4","m_t2","m_m1")
list_helper <- list()
i <- 0
for(mat in matrices){
  i <- i + 1
  temp <- get(mat)
  mgenes <- NULL
  for(cm in noms[-1]){
    cat("Working on",mat,"\t",cm,"\n")
    mgenes <- c(mgenes,extract_comb(m = temp,comb_name = cm))
  }
  list_helper[[i]] <-mgenes 
  
}

names(list_helper) <- c("s4","t2","m1")


m <- list_helper %>%list_to_matrix %>%
  make_comb_mat


cs = comb_size(m)
pdf(file = "../figures/rnaseq_upset_down_8h_helperdependent.pdf",width = 16,height = 12)

ht <- UpSet(m, 
            set_order = c("s4","t2","m1"),
            comb_order = 1:7, 
            pt_size = unit(10, "mm"), lwd = 4,
            top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs)),
                                                  height = unit(20, "cm"),
                                                  annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                                  annotation_name_rot = 90,
                                                  annotation_name_offset = unit(2,"cm"),
                                                  axis_param = list(
                                                    gp = gpar(fontsize = size_text_upset))
            ),
            left_annotation = rowAnnotation(
              "Set size" = anno_barplot(set_size(m),
                                        annotation_name_gp  = gpar(fontsize = size_titles_upset),
                                        axis_param = list(direction = "reverse",
                                                          gp = gpar(fontsize = size_text_upset),
                                                          labels_rot = 0
                                        ),
                                        border = FALSE, 
                                        gp = gpar(fill = "black"),
                                        width = unit(10, "cm")
              ),
              annotation_name_gp  = gpar(fontsize = size_titles_upset)
            ),
            right_annotation =NULL,
            row_names_side = "right",
            row_names_gp = gpar(fontsize = size_sets_upset),
)
ht = draw(ht)
co = column_order(ht)

nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = size_text_upset), 
            just = "bottom",
            default.units = "native")
})
dev.off()


m_helper_dependent_8h_down <- m


extract_comb(m = m_eti_8h_down,comb_name = "111") %>%
  intersect(
    extract_comb(m = m_helper_dependent_8h_down,comb_name = "111")
  )


### Save matrices ###
mlists <- list(
  m_eti_4h_up = m_eti_4h_up,
  m_eti_8h_up = m_eti_8h_up,
  m_helper_dependent_4h_up = m_helper_dependent_4h_up,
  m_helper_dependent_8h_up = m_helper_dependent_8h_up,
  m_eti_4h_down = m_eti_4h_down,
  m_eti_8h_down = m_eti_8h_down,
  m_helper_dependent_4h_down = m_helper_dependent_4h_down,
  m_helper_dependent_8h_down = m_helper_dependent_8h_down
)


saveRDS(object = mlists,file = "../cleandata/helperless_contrasts_results_matrices.RDS")


