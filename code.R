library(DESeq2)
library(ggrepel)
library(pheatmap)
library(tidyverse)    
library(magrittr)    
library(WGCNA)
library(gprofiler2)

dds <- DESeqDataSetFromMatrix(countData=combined_counts, 
                              +                               colData=metadata, 
                              +                               design=~dex, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)

sig_genes <- res[which(abs(res$log2FoldChange) > 1 & res$padj < 0.01), ]
significant_gene_names <- rownames(sig_genes)
write.csv(significant_gene_names, file = "Deseq significant_genes.txt", quote = FALSE, row.names = FALSE)

#save up/downregulated genes lists
upregulated_genes <- subset(res, log2FoldChange > 1 & padj < 0.01)
downregulated_genes <- subset(res, log2FoldChange < -1 & padj < 0.01)
upregulated_gene_names <- rownames(upregulated_genes)
downregulated_gene_names <- rownames(downregulated_genes)
write.csv(upregulated_gene_names, file = "Deseq upregulated_genes.txt", quote = FALSE, row.names = FALSE)
write.csv(downregulated_gene_names, file = "Deseq downregulated_genes.txt", quote = FALSE, row.names = FALSE)


###########
plotMA(res, main="MA Plot", ylim=c(-10,10), xlab="Mean of normalised counts", ylab="Log2 fold change")
head(results(dds, tidy=TRUE))
res <- res[order(res$padj),]
head(res)

#use plotCounts fxn to compare the normalised counts
#between inoculated and control groups for top 6 genes
par(mfrow=c(2,3))
plotCounts(dds, gene="MTR_2g019250", intgroup="dex")
plotCounts(dds, gene="MTR_6g005360", intgroup="dex")
plotCounts(dds, gene="MTR_1g090667", intgroup="dex")
plotCounts(dds, gene="MTR_6g005340", intgroup="dex")
plotCounts(dds, gene="MTR_7g012070", intgroup="dex")
plotCounts(dds, gene="MTR_6g005330", intgroup="dex")
############

# Calculate the variance stabilising transformation
vsd <- vst(dds, blind=FALSE)

# Perform PCA
pca_data <- plotPCA(vsd, intgroup=c("dex"), returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Plot PCA
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = dex)) +
  geom_point(aes(shape = dex), size = 3) +  
  theme_bw(base_size = 14) +  
  theme(legend.position = "right",  
        legend.background = element_blank(), 
        legend.key = element_blank(), 
        plot.background = element_blank(),  
        panel.background = element_rect(fill = "grey90", color = NA), 
        panel.grid.major = element_line(color = "white"), 
        panel.grid.minor = element_blank(), 
        axis.title = element_text(face = "bold"),  
        plot.title = element_text(hjust = 0.5)) + 
  scale_shape_manual(values = c(16, 17)) +  
  scale_color_manual(values = c("Inoculated" = "red", "Control" = "blue")) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"), 
       y = paste0("PC2: ", percentVar[2], "% variance"), 
       color = "Group", 
       shape = "Group") +  
  ggtitle("PCA Plot") 
pca_plot <- pca_plot + geom_text_repel(aes(label=rownames(pca_data)),
                                                                               size=3,
                                                                               box.padding = unit(0.35, "lines"),
                                                                               point.padding = unit(0.5, "lines"))
print(pca_plot)

#plot sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
pheatmap(sampleDistMatrix ,
         +          clustering_distance_rows = sampleDists,
         +          clustering_distance_cols = sampleDists,
         +          color = colorRampPalette(c("darkblue", "white"))(100), main = "Sample-to-sample distances")
par(mfrow=c(1,1))
plotDispEsts(dds)


#########
#wgcna
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # 95 quantile to reduce dataset
expr_normalised <- wpn_vsd[ rv_wpn > q95_wpn, ]
inputmatrix = t(expr_normalised) #transpose matric for wgcna


allowWGCNAThreads()
powers <- c(c(1:10), seq(from = 10, to = 30, by = 5), seq(from = 30, to = 44, by = 2))
sft <- pickSoftThreshold(inputmatrix, powerVector = powers)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h = 0.90, col = "blue")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h = 0.0, col = "blue")

picked_power = 25
temp_cor <- cor       
cor <- WGCNA::cor
netwk <- blockwiseModules(inputmatrix,    
                          # == Adjacency Function ==
                          power = picked_power,               
                          networkType = "signed",
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor
mergedColors = labels2colors(netwk$colors)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

MEs0 <- moduleEigengenes(inputmatrix, mergedColors)$eigengenes
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
MEs0$treatment = row.names(MEs0)
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

new_labels <- c("Control 1", "Control 2", "Control 3", "Control 4", "Control 5", "Control 6", 
                "Inoculated 1", "Inoculated 2", "Inoculated 3", "Inoculated 4", "Inoculated 5", "Inoculated 6")

mME$Treatment <- factor(mME$treatment, levels = unique(mME$treatment), labels = new_labels)

# Module-trait heatmap
mME %>%
  ggplot(aes(x = Treatment, y = name, fill = value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


modules_of_interest = c("blue", "yellow", "turquoise")
submod = module_df %>%
  subset(colors %in% modules_of_interest)
row.names(module_df) = module_df$gene_id
subexpr = expr_normalised[submod$gene_id,]
submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )
submod_df$Treatment <- factor(submod_df$name, levels = submod_df$name, labels = new_labels)
submod_df %>%
  ggplot(aes(x = Treatment, y = value, group = gene_id)) +
  geom_line(aes(color = module), alpha = 0.2) +
  scale_color_manual(values = c("yellow" = "yellow", "turquoise" = "turquoise", "blue" = "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(rows = vars(module)) +
  labs(x = "Treatment", y = "Normalised expression")


#functional annotation with gprofiler2
results <- gconvert(query = all_gs$V1, organism = "mtruncatula", target = "ENSG")
matched_genes <- results[, -c(1:3)]

genes_turq <- module_df$gene_id[module_df$colors == "turquoise"]
genes_yellow <- module_df$gene_id[module_df$colors == "yellow"]
genes_blue <- module_df$gene_id[module_df$colors == "blue"]

ENSG_turq <- matched_genes[matched_genes$target %in% genes_turq, ]
ENSG_yellow <- matched_genes[matched_genes$target %in% genes_yellow, ]
ENSG_blue <- matched_genes[matched_genes$target %in% genes_blue, ]

write.table(ENSG_turq, file = "ENSG_turq.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ENSG_yellow, file = "ENSG_yellow.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ENSG_blue, file = "ENSG_blue.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

gost_turq <- gost(
  query = ENSG_turq$target,
  organism = "mtruncatula", 
  sources = c("GO:MF", "GO:BP", "GO:CC", "KEGG"), 
  user_threshold = 0.05, 
  ordered_query = FALSE)
gost_yellow <- gost(
  query = ENSG_yellow$target,
  organism = "mtruncatula", 
  sources = c("GO:MF", "GO:BP", "GO:CC", "KEGG"),  
  user_threshold = 0.05,  
  ordered_query = FALSE)
#no results for yellow at 0.05, so not a significant module
gost_blue <- gost(
  query = ENSG_blue$target,
  organism = "mtruncatula", 
  sources = c("GO:MF", "GO:BP", "GO:CC", "KEGG"), 
  user_threshold = 0.05,  
  ordered_query = FALSE)

pf1 <- gostplot(gost_turq, capped = FALSE, interactive = FALSE)
pf2 <- gostplot(gost_blue, capped = FALSE, interactive = FALSE)
publish_gosttable(gost_turq,
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)
publish_gosttable(gost_blue,
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)

publish_gostplot(pf1, highlight_terms = c(#terms of interest), 
                 width = 20, height = 5, filename = "gostplot+table_blue.pdf")
publish_gostplot(pf2, highlight_terms = c(#terms of interest), 
                 width = 20, height = 5, filename = "gostplot+table_blue.pdf")


#create gem file for cytoscape
gostgem_turq <- gost(
  query = ENSG_turq$target,
  organism = "mtruncatula", 
  evcodes = TRUE, 
  multi_query = FALSE, 
  sources = c("GO", "REAC", "MIRNA", "CORUM", "HP", "HPA", "WP"),
  user_threshold = 0.05,
  ordered_query = FALSE)
gem_turq <- gostgem_turq$result[,c("term_id", "term_name", "p_value", "intersection")]
colnames(gem_turq) <- c("GO.ID", "Description", "p.Val", "Genes")
gem_turq$FDR <- gem_turq$p.Val
gem_turq$Phenotype = "+1"
gem_turq <- gem_turq[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
write.table(gem_turq, file = "gProfiler_gem_turq.txt", sep = "\t", quote = F, row.names = F)

gostgem_blue <- gost(
  query = ENSG_blue$target,
  organism = "mtruncatula", 
  evcodes = TRUE, 
  multi_query = FALSE, 
  sources = c("GO", "REAC", "MIRNA", "CORUM", "HP", "HPA", "WP"),
  user_threshold = 0.05,
  ordered_query = FALSE)
gem_blue <- gostgem_blue$result[,c("term_id", "term_name", "p_value", "intersection")]
colnames(gem_blue) <- c("GO.ID", "Description", "p.Val", "Genes")
gem_blue$FDR <- gem_blue$p.Val
gem_blue$Phenotype = "+1"
gem_blue <- gem_blue[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
write.table(gem_blue, file = "gProfiler_gem_blue.txt", sep = "\t", quote = F, row.names = F)
