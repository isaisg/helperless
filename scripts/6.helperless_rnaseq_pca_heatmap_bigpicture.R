library(ohchibi)
library(DESeq2)
library(paletteer)
library(scales)
library(ggtree)
library(egg)
library(clusterProfiler)
library(org.At.tair.db)
library(UpSetR)
library(extrafont)
loadfonts(device = "pdf")

setwd('/home/isai/Documents/results/helperless/scripts/')

df <- read.table(file = "../cleandata/deseq2_res_contrasts.intra.tsv",header = T,sep = ",")

#Adjust the levels so it keeps the same structure
df$Genotype <- df$Genotype %>% factor(levels = c("Col","ADR1","nrg1","HPK"))

### Prepare matrix for heatmap###
Dat <- readRDS(file = "../cleandata/dat_rnaseq.RDS")
dds <- readRDS(file = "../cleandata/dds_rnaseq.RDS")

Tab_z <- dds %>% vst %>% assay %>% t %>% scale %>% t

### Create Dat_z object to do permanova analysis
Dat_z <- create_dataset(Tab = Tab_z,Map = Dat$Map)
Dat_z$Map$Treatment <- Dat_z$Map$Treatment %>% as.character
Dat_z$Map$Treatment[is.na(Dat_z$Map$Treatment)] <- "T0"
Dat_z$Map$Treatment <- Dat_z$Map$Treatment %>% factor

Dat_z$Map$Treatment <- Dat_z$Map$Treatment %>% 
  factor(levels = c("T0","EV","M1","S4","T2"))
Dat_z$Map$Genotype <- Dat_z$Map$Genotype %>% 
  factor(levels = c("Col","ADR1","nrg1","HPK"))
Dat_z$Map$Time <- Dat_z$Map$Time %>% factor(levels = c("T0","30m","4h","8h"))

#Permanova using euclidean matrix
distfun <- function(x,method) vegan::vegdist(x = x, method = "euclidean")
Tab_euc <- distfun(t(Dat_z$Tab))


mypermanova <- adonis(Tab_euc ~  Time + Genotype +Treatment + Rep,
                      data = Dat_z$Map,
                      permutations = 9999)

mypermanova

res <- chibi.permanova(mypermanova = mypermanova,terms_exclude_plot = "Rep") 
paleta_perm <- c(paletteer_d(package = "ochRe",palette = "namatjira_qual")[c(1,3,5)],"white")

names(paleta_perm) <- c("Genotype","Treatment","Time","Residual")
p_perm <- res$p + scale_fill_manual(values = paleta_perm)  +
  theme(
    legend.position = "none"
  )

#Do global PCA 
mpca <- prcomp(x = Dat_z$Tab %>% t,center = F,scale. = F)
scores <- mpca$x %>% as.data.frame
scores$Sample_Id <- rownames(scores)
rownames(scores) <- NULL
scores <- merge(scores,Dat_z$Map, by = "Sample_Id")

paleta_time <- paletteer_d(package = "awtools",palette = "mpalette")[1:4]
#paleta_time <- paletteer_d(package = "rcartocolor",palette = "Safe")[1:4]



summary(mpca)
p_time <- ggplot(data = scores,aes(PC1,PC2)) +
  geom_vline(xintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
  geom_point(shape =21, aes(fill = Time),size = 10,stroke = 1) +
  xlab(label = "PC1 32.08%") + ylab(label = "PC2 10.40%") +
  theme_ohchibi(size_panel_border = 1) +
  scale_fill_manual(values = paleta_time) +
  theme(
    legend.position = "none"
  )


paleta_treatment <- c(paleta_time[1],
                      paletteer_d(package = "awtools",
                                  palette = "ppalette")[c(7,3,5,2)])

p_treatment <- ggplot(data = scores,aes(PC1,PC2)) +
  geom_vline(xintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
  geom_point(shape =21, aes(fill = Treatment),size = 10,stroke = 1) +
  xlab(label = "PC1 32.08%") + ylab(label = "PC2 10.40%") +
  theme_ohchibi(size_panel_border = 1) +
  scale_fill_manual(values = paleta_treatment) +
  theme(legend.position = "none")


p_geno <- ggplot(data = scores,aes(PC1,PC2)) +
  geom_vline(xintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
  geom_point(shape =21, aes(fill = Genotype),size = 10,stroke = 1) +
  xlab(label = "PC1 32.08%") + ylab(label = "PC2 10.40%") +
  theme_ohchibi(size_panel_border = 1) +
  theme(legend.position = "none") +
  scale_fill_paletteer_d(package = "ochRe",palette = "tasmania")



composition <- egg::ggarrange(p_perm,p_time,p_treatment,p_geno,
                              nrow = 1,widths = c(0.05,1,1,1))
oh.save.pdf(p = composition,outname = "pca_rnaseq_permanova.pdf",
            outdir = "../figures/",
            width = 30,height = 10,family = "Arial")

### Figure with shapes for genotypes ###
p <- ggplot(data = scores,aes(PC1,PC2)) +
  geom_vline(xintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "dashed",size = 4 , color = "#D9D9D9") +
  geom_point( aes(fill = Treatment,shape = Genotype),size = 8,stroke = 1) +
  xlab(label = "PC1 32.08%") + ylab(label = "PC2 10.40%") +
  theme_ohchibi(size_panel_border = 2) +
  scale_fill_manual(values = paleta_treatment) +
  scale_shape_manual(values = c(21,22,23,24)) 

oh.save.pdf(p = p,outname = "pca_rnaseq_coleffector_shapegenotype.pdf",
            outdir = "../figures/",
            width = 12,height = 10,family = "Arial")




## Create heatmap
#Do a summary per treatment and then create a matrix to cluste rgenes
melted <- Dat_z$Tab %>% melt
colnames(melted) <- c("Gene","Sample_Id","zscore")

melted <- merge(melted,Dat_z$Map,by = "Sample_Id") 

Tab_sum <- acast(data = melted,formula = Gene~group,
                 fun.aggregate = mean,value.var = "zscore")


#Subset only significant genes
chosen_genes <- df %>% subset(padj < 0.05 & abs(log2FoldChange) > 1) %$% Gene %>%
  as.character %>% unique 

Tab_sum <- Tab_sum[rownames(Tab_sum) %in% chosen_genes,]

write.table(x = Tab_sum,file = "../cleandata/suptable_fig4ab.csv",append = F,quote = F,sep = ",",row.names = T,col.names = T)


#Cluster genes
dist_genes <- as.dist(1-cor(Tab_sum %>% t))
mclust_genes <- hclust(d = dist_genes,method = "complete")


mclust_genes %>% plot
#8 clusters using ward.D2
df_cg <- mclust_genes %>% cutree(k = 19) %>% 
  data.frame(Gene = names(.),
             ClusterGene = paste0("CG",.),row.names = NULL)

df_cg <- df_cg[,-1]


order_genes <- mclust_genes$order %>% mclust_genes$labels[.]

melted <- Tab_sum %>% melt
colnames(melted) <- c("Gene","group","zscore")
melted <- merge(melted,Dat_z$Map[,c(3,4,5,7)] %>% unique,by = "group") 

#Merge with cluster of genes
melted <- merge(melted,df_cg, by = "Gene")

#Order the data.frame genes
melted$Gene <- melted$Gene %>% factor(levels = order_genes)

#Determine the order of the clusters
order_cg <- melted[with(melted,order(Gene)),] %$% ClusterGene %>% 
  as.character %>% unique %>% rev
melted$ClusterGene <- melted$ClusterGene %>% factor(levels = order_cg)


#Define limits of palette for scales argument
melted$zscore %>% sort %>% plot
#-3 and 3 are the limits

#Order the groups
order_group <- with(melted,order(Time,Treatment,Genotype)) %>%
  melted[.,] %$% group %>% unique

melted$group <- melted$group %>% factor(levels = order_group)
melted$Treatment <- melted$Treatment %>%
  factor(levels = c("T0","EV","S4","T2","M1"))

p_heatmap <- ggplot(data = melted,aes(group,Gene)) +
  geom_raster(aes(fill = zscore)) +  
  scale_fill_paletteer_c(package = "pals",palette = "ocean.balance",
                         limits = c(-3,3),oob = squish) +
  facet_grid(ClusterGene~Time+Treatment,space = "free",scales = "free") +
  theme_ohchibi(size_axis_text.x = 20,
                angle_text.x = 90,
                size_axis_text.y = 20,
                size_axis_title.x = 22,
                size_axis_title.y = 0,
                legend_proportion_size = 2,
                size_title_text = 30,
                size_legend_text = 25,
                size_panel_border = 1.5,
                size_lines_panel = 0) +
  theme(axis.ticks = element_blank(),
        panel.spacing = unit(0.075, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_text(family = "Arial",face = "bold",size = 30),
        strip.text.y = element_text(family = "Arial",face = "bold",size = 10,angle = 0),
        axis.text.y = element_blank(),
        axis.text.x = element_text(hjust = 1,vjust = 0.5,family = "Arial",face = "bold",size = 25),
        axis.title.x = element_blank()
  ) +
  scale_x_discrete(expand = c(0,0))


oh.save.pdf(p = p_heatmap,outname = "heatmap.bigpicture.npars.pdf",
            outdir = "../figures/",height = 16,width = 20)


