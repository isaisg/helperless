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

#Red ojbect
mlist <- readRDS(file = "../cleandata/list_contrasts.RDS")

#Load dataset
load(file = "../cleandata/helperless_contrasts_results.RData")

###### Up ##############


### S4 ###
lista <- list(
  col_eti_4h_s4 = mlist$Up$col_s4_4h_up,
  adr1_sufficient_4h_s4 = adr1_sufficient_4h_s4_up,
  adr1_requiring_4h_s4 = adr1_requiring_4h_s4_up,
  nrg1_sufficient_4h_s4 = nrg1_sufficient_4h_s4_up,
  nrg1_requiring_4h_s4 = nrg1_requiring_4h_s4_up,
  helper_dependent_4h_s4 = helper_dependent_4h_s4_up
)



#Try complext heatmap implementation
#Solution here  to add number of interesections
#https://support.bioconductor.org/p/118557/
m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_4h_s4 <- m

pdf(file = "../figures/upset_4h_s4_eti.up.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_4h_s4",
                          "nrg1_requiring_4h_s4","nrg1_sufficient_4h_s4",
                          "adr1_requiring_4h_s4","adr1_sufficient_4h_s4",
                          "helper_dependent_4h_s4"),
            comb_order = 1:5, 
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

### Analysis for T2


lista <- list(
  col_eti_4h_t2 = mlist$Up$col_t2_4h_up,
  adr1_sufficient_4h_t2 = adr1_sufficient_4h_t2_up,
  adr1_requiring_4h_t2 = adr1_requiring_4h_t2_up,
  nrg1_sufficient_4h_t2 = nrg1_sufficient_4h_t2_up,
  nrg1_requiring_4h_t2 = nrg1_requiring_4h_t2_up,
  helper_dependent_4h_t2 = helper_dependent_4h_t2_up
)


m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_4h_t2 <- m

pdf(file = "../figures/upset_4h_t2_eti.up.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_4h_t2",
                          "nrg1_requiring_4h_t2","nrg1_sufficient_4h_t2",
                          "adr1_requiring_4h_t2","adr1_sufficient_4h_t2",
                          "helper_dependent_4h_t2"),
            comb_order = 1:5, 
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


### Analysis for M1

lista <- list(
  col_eti_4h_m1 = mlist$Up$col_m1_4h_up,
  adr1_sufficient_4h_m1 = adr1_sufficient_4h_m1_up,
  adr1_requiring_4h_m1 = adr1_requiring_4h_m1_up,
  nrg1_sufficient_4h_m1 = nrg1_sufficient_4h_m1_up,
  nrg1_requiring_4h_m1 = nrg1_requiring_4h_m1_up,
  helper_dependent_4h_m1 = helper_dependent_4h_m1_up
)



m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_4h_m1 <- m

pdf(file = "../figures/upset_4h_m1_eti.up.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_4h_m1","nrg1_requiring_4h_m1","nrg1_sufficient_4h_m1",
                          "adr1_requiring_4h_m1","adr1_sufficient_4h_m1",
                          "helper_dependent_4h_m1"),
            comb_order = 1:5, 
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


### 8 hours ###
lista <- list(
  col_eti_8h_s4 = mlist$Up$col_s4_8h_up,
  adr1_sufficient_8h_s4 = adr1_sufficient_8h_s4_up,
  adr1_requiring_8h_s4 = adr1_requiring_8h_s4_up,
  nrg1_sufficient_8h_s4 = nrg1_sufficient_8h_s4_up,
  nrg1_requiring_8h_s4 = nrg1_requiring_8h_s4_up,
  helper_dependent_8h_s4 = helper_dependent_8h_s4_up
)



#Try complext heatmap implementation
#Solution here  to add number of interesections
#https://support.bioconductor.org/p/118557/
m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_8h_s4 <- m

pdf(file = "../figures/upset_8h_s4_eti.up.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_8h_s4",
                          "nrg1_requiring_8h_s4","nrg1_sufficient_8h_s4",
                          "adr1_requiring_8h_s4","adr1_sufficient_8h_s4",
                          "helper_dependent_8h_s4"),
            comb_order = 1:5, 
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

### Analysis for T2
lista <- list(
  col_eti_8h_t2 = mlist$Up$col_t2_8h_up,
  adr1_sufficient_8h_t2 = adr1_sufficient_8h_t2_up,
  adr1_requiring_8h_t2 = adr1_requiring_8h_t2_up,
  nrg1_sufficient_8h_t2 = nrg1_sufficient_8h_t2_up,
  nrg1_requiring_8h_t2 = nrg1_requiring_8h_t2_up,
  helper_dependent_8h_t2 = helper_dependent_8h_t2_up
)


m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_8h_t2 <- m

pdf(file = "../figures/upset_8h_t2_eti.up.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_8h_t2",
                          "nrg1_requiring_8h_t2","nrg1_sufficient_8h_t2",
                          "adr1_requiring_8h_t2","adr1_sufficient_8h_t2",
                          "helper_dependent_8h_t2"),
            comb_order = 1:5, 
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


### Analysis for M1

lista <- list(
  col_eti_8h_m1 = mlist$Up$col_m1_8h_up,
  adr1_sufficient_8h_m1 = adr1_sufficient_8h_m1_up,
  adr1_requiring_8h_m1 = adr1_requiring_8h_m1_up,
  nrg1_sufficient_8h_m1 = nrg1_sufficient_8h_m1_up,
  nrg1_requiring_8h_m1 = nrg1_requiring_8h_m1_up,
  helper_dependent_8h_m1 = helper_dependent_8h_m1_up
)



m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_8h_m1 <- m

pdf(file = "../figures/upset_8h_m1_eti.up.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_8h_m1","nrg1_requiring_8h_m1","nrg1_sufficient_8h_m1",
                          "adr1_requiring_8h_m1","adr1_sufficient_8h_m1",
                          "helper_dependent_8h_m1"),
            comb_order = 1:5, 
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

#### Down ###

### S4 ###
lista <- list(
  col_eti_4h_s4 = mlist$Down$col_s4_4h_down,
  adr1_sufficient_4h_s4 = adr1_sufficient_4h_s4_down,
  adr1_requiring_4h_s4 = adr1_requiring_4h_s4_down,
  nrg1_sufficient_4h_s4 = nrg1_sufficient_4h_s4_down,
  nrg1_requiring_4h_s4 = nrg1_requiring_4h_s4_down,
  helper_dependent_4h_s4 = helper_dependent_4h_s4_down
)



#Try complext heatmap implementation
#Solution here  to add number of interesections
#https://sdownport.bioconductor.org/p/118557/
m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_4h_s4 <- m

pdf(file = "../figures/upset_4h_s4_eti.down.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_4h_s4",
                          "nrg1_requiring_4h_s4","nrg1_sufficient_4h_s4",
                          "adr1_requiring_4h_s4","adr1_sufficient_4h_s4",
                          "helper_dependent_4h_s4"),
            comb_order = 1:5, 
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

### Analysis for T2


lista <- list(
  col_eti_4h_t2 = mlist$Down$col_t2_4h_down,
  adr1_sufficient_4h_t2 = adr1_sufficient_4h_t2_down,
  adr1_requiring_4h_t2 = adr1_requiring_4h_t2_down,
  nrg1_sufficient_4h_t2 = nrg1_sufficient_4h_t2_down,
  nrg1_requiring_4h_t2 = nrg1_requiring_4h_t2_down,
  helper_dependent_4h_t2 = helper_dependent_4h_t2_down
)


m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_4h_t2 <- m

pdf(file = "../figures/upset_4h_t2_eti.down.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_4h_t2",
                          "nrg1_requiring_4h_t2","nrg1_sufficient_4h_t2",
                          "adr1_requiring_4h_t2","adr1_sufficient_4h_t2",
                          "helper_dependent_4h_t2"),
            comb_order = 1:5, 
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


### Analysis for M1

lista <- list(
  col_eti_4h_m1 = mlist$Down$col_m1_4h_down,
  adr1_sufficient_4h_m1 = adr1_sufficient_4h_m1_down,
  adr1_requiring_4h_m1 = adr1_requiring_4h_m1_down,
  nrg1_sufficient_4h_m1 = nrg1_sufficient_4h_m1_down,
  nrg1_requiring_4h_m1 = nrg1_requiring_4h_m1_down,
  helper_dependent_4h_m1 = helper_dependent_4h_m1_down
)



m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_4h_m1 <- m

pdf(file = "../figures/upset_4h_m1_eti.down.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_4h_m1","nrg1_requiring_4h_m1","nrg1_sufficient_4h_m1",
                          "adr1_requiring_4h_m1","adr1_sufficient_4h_m1",
                          "helper_dependent_4h_m1"),
            comb_order = 1:5, 
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


### 8 hours ###
lista <- list(
  col_eti_8h_s4 = mlist$Down$col_s4_8h_down,
  adr1_sufficient_8h_s4 = adr1_sufficient_8h_s4_down,
  adr1_requiring_8h_s4 = adr1_requiring_8h_s4_down,
  nrg1_sufficient_8h_s4 = nrg1_sufficient_8h_s4_down,
  nrg1_requiring_8h_s4 = nrg1_requiring_8h_s4_down,
  helper_dependent_8h_s4 = helper_dependent_8h_s4_down
)



#Try complext heatmap implementation
#Solution here  to add number of interesections
#https://sdownport.bioconductor.org/p/118557/
m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_8h_s4 <- m

pdf(file = "../figures/upset_8h_s4_eti.down.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_8h_s4",
                          "nrg1_requiring_8h_s4","nrg1_sufficient_8h_s4",
                          "adr1_requiring_8h_s4","adr1_sufficient_8h_s4",
                          "helper_dependent_8h_s4"),
            comb_order = 1:5, 
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

### Analysis for T2
lista <- list(
  col_eti_8h_t2 = mlist$Down$col_t2_8h_down,
  adr1_sufficient_8h_t2 = adr1_sufficient_8h_t2_down,
  adr1_requiring_8h_t2 = adr1_requiring_8h_t2_down,
  nrg1_sufficient_8h_t2 = nrg1_sufficient_8h_t2_down,
  nrg1_requiring_8h_t2 = nrg1_requiring_8h_t2_down,
  helper_dependent_8h_t2 = helper_dependent_8h_t2_down
)


m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_8h_t2 <- m

pdf(file = "../figures/upset_8h_t2_eti.down.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_8h_t2",
                          "nrg1_requiring_8h_t2","nrg1_sufficient_8h_t2",
                          "adr1_requiring_8h_t2","adr1_sufficient_8h_t2",
                          "helper_dependent_8h_t2"),
            comb_order = 1:5, 
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


### Analysis for M1

lista <- list(
  col_eti_8h_m1 = mlist$Down$col_m1_8h_down,
  adr1_sufficient_8h_m1 = adr1_sufficient_8h_m1_down,
  adr1_requiring_8h_m1 = adr1_requiring_8h_m1_down,
  nrg1_sufficient_8h_m1 = nrg1_sufficient_8h_m1_down,
  nrg1_requiring_8h_m1 = nrg1_requiring_8h_m1_down,
  helper_dependent_8h_m1 = helper_dependent_8h_m1_down
)



m <- list_to_matrix(lista) %>%
  make_comb_mat()
mat_8h_m1 <- m

pdf(file = "../figures/upset_8h_m1_eti.down.pdf",width = 16,height = 12)

cs = comb_size(m)
ht <- UpSet(m, 
            set_order = c("col_eti_8h_m1","nrg1_requiring_8h_m1","nrg1_sufficient_8h_m1",
                          "adr1_requiring_8h_m1","adr1_sufficient_8h_m1",
                          "helper_dependent_8h_m1"),
            comb_order = 1:5, 
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

