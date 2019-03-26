# We will be using DESeq2 for performing the differential expression analysis and additional R packages for data wrangling and plotting
# RNA-Seq analysis pipeline: 1. read raw counts. 2. Create a DESeq2 object 3. Create a DESeq2 object and Run the DESeq2 analysis 3. PCA and MA plot. 4. heatmap of differently regulated genes

library(DESeq2)

wt_rawcounts <- read.csv("fibrosis_wt_rawcounts.csv")
View(wt_rawcounts)

# Create genotype vector
genotype <- c("smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe")

# Create condition vector
condition <- c("fibrosis", "fibrosis", "fibrosis", "fibrosis", "normal", "normal", "normal")

# Create data frame
smoc2_metadata <- data.frame(genotype, condition)

# Assign the row names of the data frame
rownames(smoc2_metadata) <- c("smoc2_fibrosis1", "smoc2_fibrosis2", "smoc2_fibrosis3", "smoc2_fibrosis4", "smoc2_normal1", "smoc2_normal3", "smoc2_normal4")

# Use the match() function to reorder the columns of the raw counts
match(rownames(smoc2_metadata), colnames(smoc2_rawcounts))

# Reorder the columns of the count data
reordered_smoc2_rawcounts <- smoc2_rawcounts[ , match(rownames(smoc2_metadata), colnames(smoc2_rawcounts))]

# Create a DESeq2 object
dds_smoc2 <- DESeqDataSetFromMatrix(countData = reordered_smoc2_rawcounts,
                              colData = smoc2_metadata,
                              design = ~ condition)

# Determine the size factors to use for normalization
dds_smoc2 <- estimateSizeFactors(dds_smoc2)

# Extract the normalized counts
smoc2_normalized_counts <- counts(dds_smoc2, normalized=TRUE)


# Transform the normalized counts
vsd_smoc2 <- vst(dds_smoc2, blind = TRUE)

# Extract the matrix of transformed counts
vsd_mat_smoc2 <- assay(vsd_smoc2)

# Compute the correlation values between samples
vsd_cor_smoc2 <- cor(vsd_mat_smoc2)

# Plot the heatmap
pheatmap(vsd_cor_smoc2, annotation = select(smoc2_metadata, condition))

# Transform the normalized counts
vsd_smoc2 <- vst(dds_smoc2, blind = TRUE)

# Plot the PCA of PC1 and PC2
plotPCA(vsd_smoc2, intgroup="condition")

# Create DESeq2 object
dds_smoc2 <- DESeqDataSetFromMatrix(countData = reordered_smoc2_rawcounts,
                              colData = smoc2_metadata,
                              design = ~ condition)

# Run the DESeq2 analysis
dds_smoc2 <- DESeq(dds_smoc2)

# Plot dispersions
plotDispEsts(dds_smoc2)

# Extract the results of the differential expression analysis
smoc2_res <- results(dds_smoc2,
                contrast = c("condition", "fibrosis", "normal"),
                alpha = 0.05)


# Shrink the log2 fold change estimates to be more accurate
smoc2_res <- lfcShrink(dds_smoc2,
                    contrast = c("condition", "fibrosis", "normal"),
                    res = smoc2_res)


# Explore the results() function
?results

# Extract results
smoc2_res <- results(dds_smoc2,
            contrast = c("condition", "fibrosis", "normal"),
            alpha = 0.05,
            lfcThreshold = 0.32)

# Shrink the log2 fold changes
smoc2_res <- lfcShrink(dds_smoc2,
            contrast = c("condition", "fibrosis", "normal"),
            res = smoc2_res)

# Get an overview of the results
summary(smoc2_res)

# Save results as a data frame
smoc2_res_all <- data.frame(smoc2_res)

# Subset the results to only return the significant genes with p-adjusted values less than 0.05
smoc2_res_sig <- subset(smoc2_res_all, padj < 0.05)


# Create MA plot
plotMA(smoc2_res)

# Generate logical column
smoc2_res_all <- data.frame(smoc2_res) %>% mutate(threshold = padj < 0.05)

# Create the volcano plot
ggplot(smoc2_res_all) +
        geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
        xlab("log2 fold change") +
        ylab("-log10 adjusted p-value") +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)))


# Subset normalized counts to significant genes
sig_norm_counts_smoc2 <- normalized_counts_smoc2[rownames(smoc2_res_sig), ]

# Choose heatmap color palette
heat_colors <- brewer.pal(n = 6, name = "YlOrRd")

# Plot heatmap
pheatmap(sig_norm_counts_smoc2,
        color = heat_colors,
        cluster_rows = T,
        show_rownames = F,
        annotation = select(smoc2_metadata, condition),
        scale = "row")



####### Summary 1 ########
# Check that all of the samples are in the same order in the metadata and count data
all(rownames(all_metadata) %in% colnames(all_rawcounts))

# DESeq object to test for the effect of fibrosis regardless of genotype
dds_all <- DESeqDataSetFromMatrix(countData = all_rawcounts,
                        colData = all_metadata,
                        design = ~ genotype + condition)

# DESeq object to test for the effect of genotype on the effect of fibrosis
dds_complex <- DESeqDataSetFromMatrix(countData = all_rawcounts,
                                colData = all_metadata,
                                design = ~ genotype + condition + genotype:condition)



####### Summary 2 ########
# Log transform counts for QC
vsd_all <- vst(dds_all, blind = TRUE)

# Create heatmap of sample correlation values
vsd_all %>%
        assay() %>%
        cor() %>%
        pheatmap(annotation =  select(all_metadata, c("condition", "genotype")))

# Create the PCA plot for PC1 and PC2 and color by condition
plotPCA(vsd_all, intgroup = "condition")

# Create the PCA plot for PC1 and PC2 and color by genotype
plotPCA(vsd_all, intgroup = "genotype")


# Select significant genese with padj < 0.05
smoc2_sig <- subset(res_all, padj < 0.05) %>%
  				data.frame() %>%
  				rownames_to_column(var = "geneID")

# Extract the top 6 genes with padj values
smoc2_sig %>%
	arrange(padj) %>%
	select(geneID, padj) %>%
	head()
