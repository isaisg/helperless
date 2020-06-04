library(ohchibi)
library(UpSetR)
library(ggvenn)


setwd('/home/isai/Documents/results/helperless/scripts/')
set.seed(130816)


mlist <- readRDS(file = "../cleandata/list_contrasts.RDS")

#S4
genes_hpk_8h_s4 <- c(mlist$Up$hpk_s4_8h_up) %>% unique
genes_col_8h_s4 <- c(mlist$Up$col_s4_8h_up) %>% unique


genes_hpk_4h_s4 <- c(mlist$Up$hpk_s4_4h_up) %>% unique
genes_col_4h_s4 <- c(mlist$Up$col_s4_4h_up) %>% unique

#T2
genes_hpk_8h_t2 <- c(mlist$Up$hpk_t2_8h_up) %>% unique
genes_col_8h_t2 <- c(mlist$Up$col_t2_8h_up) %>% unique



genes_hpk_4h_t2 <- c(mlist$Up$hpk_t2_4h_up) %>% unique
genes_col_4h_t2 <- c(mlist$Up$col_t2_4h_up) %>% unique


#M1
genes_hpk_8h_m1 <- c(mlist$Up$hpk_m1_8h_up) %>% unique
genes_col_8h_m1 <- c(mlist$Up$col_m1_8h_up) %>% unique



genes_hpk_4h_m1 <- c(mlist$Up$hpk_m1_4h_up) %>% unique
genes_col_4h_m1 <- c(mlist$Up$col_m1_4h_up) %>% unique

paleta <- c("#B35D5D","#9EFFBA","#FFCB82","#549ECC")

#S4
mlist <- list(
  h4_hpk = genes_hpk_4h_s4,
  h4_col = genes_col_4h_s4,
  h8_col = genes_col_8h_s4,
  h8_hpk = genes_hpk_8h_s4
)



p_s4 <- ggvenn(data = mlist,fill_alpha = 0.6,
               fill_color = paleta) 

p_s4_upset <- upset(data = fromList(mlist),set_size.numbers_size = T,text.scale = 3)



#T2
mlist <- list(
  h4_hpk = genes_hpk_4h_t2,
  h4_col = genes_col_4h_t2,
  h8_col = genes_col_8h_t2,
  h8_hpk = genes_hpk_8h_t2
)


p_t2 <- ggvenn(data = mlist,fill_alpha = 0.6,
               fill_color = paleta)  

p_t2_upset <- upset(data = fromList(mlist),set_size.numbers_size = T,text.scale = 3)


#M1
mlist <- list(
  h4_hpk = genes_hpk_4h_m1,
  h4_col = genes_col_4h_m1,
  h8_col = genes_col_8h_m1,
  h8_hpk = genes_hpk_8h_m1
)


p_m1 <- ggvenn(data = mlist,fill_alpha = 0.6,
               fill_color =paleta) 


p_m1_upset <- upset(data = fromList(mlist),set_size.numbers_size = T,text.scale = 3)


composition <- egg::ggarrange(p_s4,p_t2,p_m1,nrow = 1)
oh.save.pdf(p = composition,outname = "venns_fast.pdf",outdir = "../figures/",width = 20,height = 8)


pdf(file = "../figures/new_upset_up_s4.pdf",width = 12,height = 12)
p_s4_upset
dev.off()

pdf(file = "../figures/new_upset_up_t2.pdf",width = 12,height = 12)
p_t2_upset
dev.off()


pdf(file = "../figures/new_upset_up_m1.pdf",width = 12,height = 12)
p_m1_upset
dev.off()
