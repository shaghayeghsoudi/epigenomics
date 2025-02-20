



# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "ReactomePA", "tximport", "edgeR"))

# Load libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)  # Annotation package for human genes
library(ReactomePA)  # Reactome Pathway Enrichment
library(tximport)  # For Salmon count files
library(edgeR)

### RNA-seq data

# Define the directory where quant.sf files are located
salmon_dir <- "path_to_salmon_quant/"

# List all quant.sf files
files <- list.files(path = salmon_dir, pattern = "quant.sf", full.names = TRUE, recursive = TRUE)

# Assign sample names based on directory structure
sample_names <- basename(dirname(files))
names(files) <- sample_names

# Import data using tximport
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Prepare metadata for DESeq2
sample_metadata <- data.frame(
  sample = sample_names,
  condition = c("Tumor", "Normal", "Tumor", "Normal"),  # Change as per experiment
  row.names = sample_names
)

# Run DESeq2 for Differential Expression Analysis
dds <- DESeqDataSetFromTximport(txi, colData = sample_metadata, design = ~ condition)
dds <- DESeq(dds)

# Extract results (log2FoldChange, p-values)
res <- results(dds, alpha = 0.05)  # FDR threshold 0.05
res <- as.data.frame(res)
res$gene <- rownames(res)

# Filter for significantly differentially expressed genes
sig_genes <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]  # Adjust cutoff as needed


##### 
## methylation significant

# Load significant hyper- and hypo-methylated genes
hyper_meth_genes <- c("NFIB", "SMARCA4", "SNX29", "STK11")  # Replace with your gene list
hypo_meth_genes <- c("MMP9", "PIK3R2", "GNA14", "ALKBH5")  # Replace with your gene list

# Convert gene names to uppercase for consistency
sig_genes$gene <- toupper(sig_genes$gene)
hyper_meth_genes <- toupper(hyper_meth_genes)
hypo_meth_genes <- toupper(hypo_meth_genes)

# Find overlapping genes
overlap_hyper <- sig_genes[sig_genes$gene %in% hyper_meth_genes, ]
overlap_hypo <- sig_genes[sig_genes$gene %in% hypo_meth_genes, ]

# Combine for pathway enrichment analysis
combined_genes <- unique(c(overlap_hyper$gene, overlap_hypo$gene))

print(length(combined_genes))  # Number of genes found in RNA-seq & methylation dat



#####

### pathway enrichment
# Convert Gene Symbols to Entrez IDs
gene_entrez <- bitr(combined_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Check if conversion was successful
head(gene_entrez)


# KEGG Pathway Enrichment
kegg_res <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",  # Human
  pvalueCutoff = 0.05
)

# View results
head(kegg_res@result)

# Plot KEGG Pathway Results
dotplot(kegg_res, showCategory = 10)  # Show top 10 pathways



# Reactome Pathway Enrichment
reactome_res <- enrichPathway(
  gene = gene_entrez$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.05
)

# View results
head(reactome_res@result)

# Plot Reactome Pathway Results
dotplot(reactome_res, showCategory = 10)


### Go analysis

# GO Enrichment Analysis
go_res <- enrichGO(
  gene = gene_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # "BP" = Biological Process, "MF" = Molecular Function, "CC" = Cellular Component
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# View results
head(go_res@result)

# Plot GO Enrichment Results
dotplot(go_res, showCategory = 10)