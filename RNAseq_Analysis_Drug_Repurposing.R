if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "hgu133a.db", "pheatmap"))

library(GEOquery)
library(limma)
library(hgu133a.db)
library(pheatmap)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "hgu133a.db", "pheatmap"))

library(GEOquery)
library(limma)
library(hgu133a.db)
library(pheatmap)
gse <- getGEO("GSE6740", GSEMatrix = TRUE)
data <- exprs(gse[[1]])     # Expression matrix
pdata <- pData(gse[[1]])    # Sample metadata
head(pdata$title)           # Check sample titles
cd4_samples <- grep("CD4", pdata$title)  # Get CD4+ sample indices
data_cd4 <- data[, cd4_samples]         # Subset expression matrix
pdata_cd4 <- pdata[cd4_samples, ]       # Subset metadata
head(pdata_cd4$title)                   # Verify CD4+ samples only
group <- ifelse(grepl("N", pdata_cd4$title), "Control", "HIV")
group <- factor(group, levels = c("Control", "HIV"))
table(group)  # Check sample distribution
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "HIV_vs_Control")

fit <- lmFit(data_cd4, design)
fit <- eBayes(fit)

deg_results <- topTable(fit, coef = "HIV_vs_Control", adjust = "fdr", number = Inf)
head(deg_results)
deg_filtered <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 1)
write.csv(deg_filtered, "GSE6740_CD4_DEGs.csv")
library(ggplot2)
deg_results$Significant <- ifelse(deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1,
                                  ifelse(deg_results$logFC > 1, "Up", "Down"), "NS")

# Volcano plot
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey"))

# Heatmap (Top 50 DEGs)
top50 <- rownames(deg_filtered)[1:50]
pheatmap(data_cd4[top50, ], annotation_col = data.frame(Group = group))
# Adjust top genes dynamically based on available DEGs
top_n <- min(50, nrow(deg_filtered))
top_genes <- rownames(deg_filtered)[1:top_n]

pheatmap(data_cd4[top_genes, ], annotation_col = data.frame(Group = group))
# Check how many DEGs passed the filter
nrow(deg_filtered)

# Check first few gene names in deg_filtered
head(rownames(deg_filtered))

# Check if they exist in data_cd4
head(rownames(data_cd4))
all(rownames(deg_filtered) %in% rownames(data_cd4))
deg_filtered <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 0.5)
top_n <- min(50, nrow(deg_filtered)) 
top_genes <- rownames(deg_filtered)[1:top_n]
mat <- data_cd4[top_genes, ]    # Subset matrix for selected genes
pheatmap(mat,
         annotation_col = data.frame(Group = group),
         scale = "row",         # Normalize expression per gene
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE)
length(top_genes)
head(top_genes)
pheatmap(mat,
         annotation_col = data.frame(Group = group),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE)
# Check matrix dimensions
dim(mat)

# Check if it contains any NA or empty values
any(is.na(mat))

# View first few rows and columns
mat[1:5, 1:5]
top_genes <- intersect(rownames(deg_filtered), rownames(data_cd4))
mat <- data_cd4[top_genes, , drop = FALSE]

pheatmap(mat,
         annotation_col = data.frame(Group = group),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE)
# Intersect genes between deg_filtered and data_cd4
top_genes <- intersect(rownames(deg_filtered), rownames(data_cd4))

# Check if we have any valid genes
if (length(top_genes) > 0) {
  mat <- data_cd4[top_genes, , drop = FALSE]
  
  # Ensure matrix is numeric
  mat <- as.matrix(mat)
  
  # Plot heatmap
  pheatmap(mat,
           annotation_col = data.frame(Group = group),
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = FALSE)
} else {
  message("No overlapping genes found between deg_filtered and data_cd4.")
}
deg_filtered <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 0.5)
# Intersect genes between deg_filtered and data_cd4
top_genes <- intersect(rownames(deg_filtered), rownames(data_cd4))

# Check if we have any valid genes
if (length(top_genes) > 0) {
  mat <- data_cd4[top_genes, , drop = FALSE]
  
  # Ensure matrix is numeric
  mat <- as.matrix(mat)
  
  # Plot heatmap
  pheatmap(mat,
           annotation_col = data.frame(Group = group),
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = FALSE)
} else {
  message("No overlapping genes found between deg_filtered and data_cd4.")
}
length(top_genes)
deg_filtered <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 0.2)
mat <- data_cd4[top_genes, , drop = FALSE]
dim(mat)
# Ensure matrix is numeric
mat <- as.matrix(mat)
mode(mat) <- "numeric"

# Ensure annotation column matches samples
annotation <- data.frame(Group = group)
rownames(annotation) <- colnames(mat)

# Plot heatmap
pheatmap(mat,
         annotation_col = annotation,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE)
# Re-generate and save Volcano plot
library(ggplot2)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point() +
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted P-Value") +
  ggsave("Volcano_Plot.png", width = 8, height = 6, dpi = 300)
library(ggplot2)

# Create volcano plot and save it to an object
volcano_plot <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point() +
  scale_color_manual(values = c("red", "blue", "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted P-Value")

# Display plot
print(volcano_plot)

# Save plot as PNG
ggsave("Volcano_Plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)
# Extract significant DEGs
deg_filtered <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 0.2)

# Extract gene symbols (assuming "Gene" column is present in deg_results)
gene_list <- deg_filtered$Gene
install.packages("clusterProfiler")
install.packages("org.Hs.eg.db")
install.packages("enrichplot")
install.packages("ggnewscale")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggnewscale)
gene_ids <- bitr(gene_list, fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)
library(clusterProfiler)
library(org.Hs.eg.db)
gene_ids <- clusterProfiler::bitr(gene_list,
                                  fromType = "SYMBOL",
                                  toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db)
# Gene vector for enrichment
entrez_ids <- gene_ids$ENTREZID

# GO Enrichment
go_results <- enrichGO(gene          = entrez_ids,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENTREZID",
                       ont           = "ALL",      # BP, CC, MF or ALL
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)

# KEGG Enrichment
kegg_results <- enrichKEGG(gene         = entrez_ids,
                           organism     = 'hsa',
                           pvalueCutoff = 0.05)
# GO barplot
barplot(go_results, showCategory = 10, title = "Top 10 GO Terms")

# KEGG barplot
barplot(kegg_results, showCategory = 10, title = "Top 10 KEGG Pathways")

# Alternatively, dotplots
dotplot(go_results, showCategory = 10, title = "GO Enrichment")
dotplot(kegg_results, showCategory = 10, title = "KEGG Pathway Enrichment")
library(enrichplot)
barplot(go_results, showCategory = 10, title = "Top 10 GO Terms")
dotplot(go_results, showCategory = 10, title = "GO Enrichment")
head(go_results)
go_results <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,    # relaxed
  qvalueCutoff  = 0.2
)
library(enrichplot)
if (nrow(as.data.frame(go_results)) > 0) {
  barplot(go_results, showCategory = 10, title = "Top 10 GO Terms")
  dotplot(go_results, showCategory = 10, title = "GO Enrichment")
} else {
  message("No GO enrichment results found.")
}
length(entrez_ids)        # Check how many mapped genes you have
head(entrez_ids)          # Preview some IDs
mapped_genes <- bitr(gene_list, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
head(mapped_genes)
# STEP 1: Relax DEG filters to get more genes
deg_filtered <- subset(deg_results, adj.P.Val < 0.1 & abs(logFC) > 0.1)

# STEP 2: Extract gene symbols and clean them
gene_list <- deg_filtered$Gene
gene_list <- toupper(trimws(gene_list))  # Uppercase + remove whitespace

# STEP 3: Map gene symbols to Entrez IDs
library(clusterProfiler)
library(org.Hs.eg.db)

mapped_genes <- bitr(gene_list,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",


                     OrgDb = org.Hs.eg.db)

# Check mapped genes
print(head(mapped_genes))
entrez_ids <- mapped_genes$ENTREZID
length(entrez_ids)
str(dataset)
head(dataset)

str(test_gene)
head(test_gene)
# Load biomaRt
library(biomaRt)

# Use Ensembl as the database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene symbol for Gene_ID column
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = dataset$Gene_ID,
  mart = mart
)

# Merge with your dataset
dataset_annotated <- merge(dataset, gene_map, by.x = "Gene_ID", by.y = "ensembl_gene_id")
library(biomaRt)

# Connect to Ensembl BioMart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map Ensembl Gene IDs to HGNC symbols and Entrez IDs
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = dataset$Gene_ID,
  mart = mart
)

# Merge with original dataset
dataset_annotated <- merge(dataset, gene_map, by.x = "Gene_ID", by.y = "ensembl_gene_id")

# Preview
head(dataset_annotated)
# Remove entries with missing Entrez IDs
gene_df <- dataset_annotated[!is.na(dataset_annotated$entrezgene_id), ]

# Optional: Filter based on expression level if needed
# For example, genes with expression > 1 (adjust as appropriate)
gene_df <- gene_df[gene_df$Expression > 1, ]

# Extract vector of Entrez IDs for enrichment
gene_entrez <- unique(gene_df$entrezgene_id)

# Check length
length(gene_entrez)
library(clusterProfiler)
library(org.Hs.eg.db)

ego <- enrichGO(
  gene         = gene_entrez,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",              # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

# View results
head(ego)
barplot(ego, showCategory = 20)
dotplot(ego, showCategory = 20)
head(gene_entrez)
length(gene_entrez)
class(gene_entrez)
gene_entrez <- as.character(gene_entrez)
library(org.Hs.eg.db)

# Check if IDs exist
valid_ids <- gene_entrez[gene_entrez %in% keys(org.Hs.eg.db, keytype = "ENTREZID")]
length(valid_ids)
ego <- enrichGO(
  gene          = valid_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)
head(ego)
barplot(ego, showCategory = 20)
head(gene_entrez)
str(gene_entrez)
length(gene_entrez)
head(mapped_genes)
head(gene_list)
str(gene_list)
library(clusterProfiler)
library(org.Hs.eg.db)

mapped_genes <- tryCatch({
  bitr(gene_list,
       fromType = "SYMBOL",
       toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)
}, error = function(e) {
  message("bitr() failed: ", e$message)
  NULL
})
gene_list[!gene_list %in% keys(org.Hs.eg.db, keytype = "SYMBOL")]
head(gene_list)
str(gene_list)
head(gene_list)
# Preview result
head(dataset_annotated)
sig_genes <- deg_results[deg_results$padj < 0.05 & abs(deg_results$log2FoldChange) > 1, ]
mapped_genes <- bitr(sig_genes$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_ids <- mapped_genes$ENTREZID
head(deg_results)
colnames(deg_results)[1] <- "gene"  # Adjust as per actual column name
deg_results <- read.csv("deg_results.csv", header = TRUE)
head(deg_results)
deg_results <- readRDS("deg_results.rds")
head(deg_results)
# Assuming you already have 'dds' created in Step 2
dds <- DESeq(dds)
res <- results(dds)
deg_results <- as.data.frame(res)
deg_results$gene <- rownames(deg_results)  # Add gene column for mapping
write.csv(deg_results, "deg_results.csv")
getwd()
setwd("C:/Users/pavan/OneDrive/Documents")
deg_results <- read.csv("deg_results.csv")
setwd('C:\Users\pavan\OneDrive\Desktop\Bioinformatics Research Projects\Projects\HIV-Drug Repurposing\Integrated Biological Networks for the Drug repurposing using ML\Step-Wise plan for project\GSE6740_CD4_DEGs.csv')
setwd("C:/Users/pavan/OneDrive/Desktop/Bioinformatics Research Projects/Project/$HIV-Drug Repurposing/...")
deg_results <- read.csv("C:/Users/pavan/OneDrive/Desktop/Bioinformatics Research Projects/Project/$HIV-Drug Repurposing/Integrated Biological Networks for the Drug repurposing using ML/Step-Wise plan for project/GSE6740_CD4_DEGS.csv")
deg_results <- read.csv("C:/Users/pavan/OneDrive/Desktop/Bioinformatics Research Projects/Project/$HIV-Drug Repurposing/Integrated Biological Networks for the Drug repurposing using ML/Step-Wise plan for project/GSE6740_CD4_DEGS.csv")
read.csv("C:/Users/pavan/OneDrive/Desktop/Bioinformatics Research Projects/Project$/HIV-Drug Repurposing/Integrated Biological Networks for the Drug repurposing using ML/Step-Wise plan for project/GSE6740_CD4_DEGS.csv")
setwd("C:/Users/pavan/OneDrive/Desktop")
a
a
library(hgu133a.db)
library(AnnotationDbi)
head(deg_results)
# Add gene symbols column using probe IDs (assumed in rownames)
deg_results$GeneSymbol <- mapIds(hgu133a.db,
                                 keys = rownames(deg_results),
                                 column = "SYMBOL",
                                 keytype = "PROBEID",
                                 multiVals = "first")
keytypes(hgu133a.db)
head(keys(hgu133a.db, keytype = "PROBEID"))
head(keys(hgu133a.db, keytype = "PROBEID"))
# Load required annotation libraries
library(hgu133a.db)
library(AnnotationDbi)

# Map probe IDs (row names of your deg_results) to gene symbols
deg_results$GeneSymbol <- mapIds(hgu133a.db,
                                 keys = rownames(deg_results),
                                 column = "SYMBOL",
                                 keytype = "PROBEID",
                                 multiVals = "first")

# Check if mapping worked
head(deg_results)
head(keys(hgu133a.db, keytype = "PROBEID"))
head(rownames(deg_results))
# Clean rownames (if they are numbers like "1", "2", "3", etc.)
head(rownames(deg_results))  # If you see "1", "2", etc. do this:
rownames(deg_results) <- deg_results$X  # Replace with actual column if it holds probe IDs
deg_results$X <- NULL  # remove extra column

# OR if you just need to trim spaces
rownames(deg_results) <- trimws(rownames(deg_results))
deg_results$GeneSymbol <- mapIds(hgu133a.db,
                                 keys = rownames(deg_results),
                                 column = "SYMBOL",
                                 keytype = "PROBEID",
                                 multiVals = "first")
head(deg_results)
# Step 5
upregulated <- deg_results[deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, ]
downregulated <- deg_results[deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, ]
up_genes <- unique(upregulated$GeneSymbol)
down_genes <- unique(downregulated$GeneSymbol)
library(clusterProfiler)
library(org.Hs.eg.db)
enrich_result <- enrichGO(gene = up_genes, 
                          OrgDb = org.Hs.eg.db, 
                          keyType = "SYMBOL", 
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.05)
barplot(enrich_result, showCategory=20)
library(stringr)

# Modify Description to wrap text at 40 characters
enrich_result@result$Description <- str_wrap(enrich_result@result$Description, width = 20)

# Plot again
barplot(enrich_result, showCategory = 20)
# Load necessary libraries
library(clusterProfiler)
library(ggplot2)
library(stringr)
library(dplyr)

# Extract the enrichment results
df <- as.data.frame(enrich_result)

# Wrap the Description text at 40 characters for better readability
df$Description <- str_wrap(df$Description, width = 20)

# Plot using ggplot2
ggplot(df[1:20, ], aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
  labs(
    title = "Top 20 GO Enrichment Terms",
    x = "GO Term",
    y = "Gene Count"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 2)
  )
library(stringr)
df$Description <- str_wrap(df$Description, width = 40)
ggplot(df[1:20, ], aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
  labs(
    title = "Top 20 GO Enrichment Terms",
    x = "GO Term",
    y = "Gene Count"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 5, lineheight = 0.8),  # key line added
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )
# Step 6
library(clusterProfiler)
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs
up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# KEGG Enrichment
kegg_result <- enrichKEGG(gene = up_entrez$ENTREZID,
                          organism = "hsa",
                          pvalueCutoff = 0.05)

# Visualize
barplot(kegg_result, showCategory = 20)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ReactomePA")
a
library(ReactomePA)
reactome_result <- enrichPathway(gene = up_entrez$ENTREZID,
                                 organism = "human",
                                 pvalueCutoff = 0.05)
barplot(reactome_result, showCategory = 20)
library(ggplot2)
library(enrichplot)  # Optional if needed

# Convert to ggplot and apply improvements
barplot(reactome_result, 
        showCategory = 20, 
        font.size = 10,        # Decrease font size
        title = "Reactome Pathway Enrichment",
        label_format = 40      # Wrap long labels after 40 characters
) + 
  theme(axis.text.y = element_text(size = 5))  # Optional fine-tuning
# Save GO enrichment
write.csv(df, "GO_Enrichment_Top20.csv", row.names = FALSE)

# Save KEGG result
write.csv(as.data.frame(kegg_result), "KEGG_Enrichment.csv", row.names = FALSE)

# Save Reactome result
write.csv(as.data.frame(reactome_result), "Reactome_Enrichment.csv", row.names = FALSE)
save.image("HIV_Drug_Repurposing_Enrichment.RData")


