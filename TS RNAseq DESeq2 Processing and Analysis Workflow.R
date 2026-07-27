##### RNAseq Processing, Normalization, and Statistical Analysis #####
# Method: DESeq
# Samples: WT and MAP3K4 Kinase Inactive (KI) Trophoblast Stem (TS) Cells
# Sample Number: 3 TS WT and 3 TS MAP3K4 KI
# Remove any other sample prior to analysis.

#Step 1: Load all required packages.
library(readxl) #Used for reading in excel files.
library(writexl) #Used for writing excel files.
library(tidyverse) #Multiple functions used for data processing.
library(DESeq2) #Used to analyze sequencing data.

#Step 2: Load all required data.
setwd("~/Library/CloudStorage/OneDrive-Personal/Trophoblast and Placenta Sequencing Datasets/Abell Mouse TSc RNAseq Data [GSE92425 (1), GSE148496 (2), Unpublished (3)]/Data") #Set working directory.

#A: Read-in datasets.
Proteins = read_xlsx('ProteinCodingGenes.xlsx') #Protein coding genes
Y_linked = read_xlsx('MouseYlinked.xlsx') #Y-linked genes

#B: Extract and Format Proteins That Are Not Y-linked.
Proteins = (Proteins %>% filter(!Symbol %in% Y_linked$Symbol))$Symbol

#C: Loading RNAseq Data Without the Non-Protein Coding Genes
M4_seq = (read_xlsx('TS_RNAseq.xlsx') %>% filter(GeneName %in% Proteins))[1:8]

#D: Create Dataframe for Count Data
M4_counts = data.frame(M4_seq[,3:8], row.names = M4_seq$GeneID)

#E: Create Dataframe for Metadata
M4_meta = data.frame(Cell.Type = colnames(M4_seq[,3:8]),
                     M4.status = as.factor(c('WT','WT','WT','KI','KI','KI')),
                     row.names = colnames(M4_seq[,3:8]))

#Step 3: Check to make sure colnames of count data match rownames of meta data.
all(colnames(M4_counts) %in% rownames(M4_meta)) #Should be TRUE
all(colnames(M4_counts) == rownames(M4_meta)) #Should be TRUE

#Step 4: Perform DESeq2 analysis pipeline.

#A: Create DESeq2 Object.
M4_deseq = DESeqDataSetFromMatrix(countData = M4_counts, colData = M4_meta, design = ~ M4.status)

#B: Normalize DESeq2 Data.
M4_norm = estimateSizeFactors(M4_deseq)

#C: Perform DESeq2 Analysis by Using the Wald Test.
M4_deseq = DESeq(M4_deseq, test = 'Wald')

#D: Extract Results From DESeq2 Analysis.
M4_results = results(M4_deseq, contrast = c('M4.status', 'KI', 'WT'), alpha = 0.05)

#E: Create Dataframe From Results.
DE_genes = as.data.frame(M4_results[order(M4_results$padj),]) 
DE_genes$GeneID = rownames(DE_genes)
DE_genes = DE_genes %>% 
  inner_join(M4_seq[,1:2], by = 'GeneID') %>%
  filter(padj < 0.05)
DE_genes = DE_genes[,c(7:8,1:6)]

TSRNAseq_DEseq_Results = as.data.frame(M4_results[order(M4_results$padj),])
TSRNAseq_DEseq_Results$GeneID = rownames(TSRNAseq_DEseq_Results)
TSRNAseq_DEseq_Results = TSRNAseq_DEseq_Results %>%
  inner_join(M4_seq[,1:2], by = 'GeneID')

# Step 5: Data Extraction

#A: Extract Results for All Genes.
setwd("~/Library/CloudStorage/OneDrive-Personal/Trophoblast and Placenta Sequencing Datasets/Abell Mouse TSc RNAseq Data [GSE92425 (1), GSE148496 (2), Unpublished (3)]/Results") #Set working directory.
write_xlsx(TSRNAseq_DEseq_Results, 'TSRNAseq_DEseq_Results.xlsx') #Save Results for All Genes as Excel File

#B: Extract Increased DEGs in MAP3K4 KI TS Cells.
up.genes = (DE_genes %>% filter(log2FoldChange > log2(1.5)))
write_xlsx(up.genes, 'MAP3K4 KI Up DEGs.xlsx') # Export for paper,

#C: Extract Decreased DEGs in MAP3K4 KI TS Cells.
down.genes = (DE_genes %>% filter(log2FoldChange < log2(0.75)))
write_xlsx(down.genes, 'MAP3K4 KI Down DEGs.xlsx') # Export for paper.

# Step 6: Data Visualization

#A: Format dataset using dplyr.
data = TSRNAseq_DEseq_Results %>%
  drop_na() %>% #Remove NAs from dataset.
  mutate(gene = GeneName) %>% # Change GeneName column to gene
  select(-c(GeneID, GeneName)) %>% #Remove unecessary columns.
  mutate(label = ifelse(gene == 'Hspb8', 'Hspb8', "")) # Label Hspb8 for paper (not necessary).

#B: Determine direction of Log2FC for each gene based on thresholds.
direction = c() # Assign empty direction vector and fill based on condition.

# Algorithm needed for filling empty direction vector based on conditions below:
# Increased gene threshold = gene log2FC > log2(1.5)
# Decreased gene threshold = gene log2FC < log2(0.75)
# Unchanged (neither) genes = log2(0.75) < gene log2FC < log2(1.5)

# Use for loop with if/else conditions:
for(i in 1:nrow(data)) {
  if(data$log2FoldChange[i] > log2(1.5)) { # Increased gene threshold = log2(1.5)
    direction = c(direction, 'up') # Increased genes (positive log2FC > threshold)
  } else if(data$log2FoldChange[i] < log2(0.75)) { # Decreased gene threshold = log2(0.75)
    direction = c(direction, 'down') # Decreased genes (negative log2FC < threshold)
  } else if (data$log2FoldChange[i] < log2(1.5) & 
             data$log2FoldChange[i] > log2(0.75)) { # Genes within thresholds (unchanged)
    direction = c(direction, 'neither') # Unchanged genes
  }
}

data$direction = direction # Assign filled direction vector to column in data.

#C: Determine significance of each gene based on thresholds.
significance = c() # Assign empty significance vector and fill based on condition.

# Algorithm needed for filling empty significance vector based on conditions below:
# Significant genes = gene padj (adjusted p-value) < 0.05
# Nonsignificant genes = gene padj (adjusted p-value) > 0.05

# Use for loop with if/else conditions:
for(i in 1:nrow(data)) {
  if(data$padj[i] < 0.05) { # Significance threshold
    significance = c(significance, 'Yes') # Significant gene (Yes)
  } else {
    significance = c(significance, 'No') # Nonsignificant gene (No)
  }
}

data$significance = significance # Assign significance vector to column in data.

#D: Export processed data.
write_xlsx(data, 'Processed_TSRNAseq_DEseq_Results.xlsx') #Save Processed data as excel file.

#E: Generate volcano plot using processed data.

# Graph format:
# Log2FC will be on the x-axis.
# -log(padj) will be on the y-axis.
# Significance threshold will be shown with lines.
# Significance and log2FC of genes will be demonstrated by a color scale:

colors = c() # Assign empty vector and fill with colors based on condition.

# Algorithm needed for filling empty significance vector based on conditions below:
# Increased DEGs (direction = up and significance = yes) will be orange red.
# Decreased DEGs (direction = down and significance = yes) will be royal blue.
# Non DEGs (genes that don't fit the criteria above) will be gray.
for(i in 1:nrow(data)) {
  if(data$direction[i] == 'up' & data$significance[i] == 'Yes') { # Increased DEG criteria.
    colors = c(colors, 'orangered') # Increased DEGs
  } else if(data$direction[i] == 'down' & data$significance[i] == 'Yes') { # Decreased DEG criteria.
    colors = c(colors, 'royalblue') # Decreased DEGs
  } else {
    colors = c(colors, 'gray35') # Non DEGs
  }
}

data$colors = colors # Assign column in data with filled colors vector.

# Generate plot using ggplot2 package.
plot = ggplot(data, 
              aes(x = log2FoldChange, 
                  y = -log2(padj))) +
  geom_point(shape = 21, stroke = 0, fill = colors) +
  geom_hline(
    yintercept = -log2(0.05), 
    linetype = "solid", 
    color = "black", 
    size = 1,
    alpha = 0.5) +
  geom_text(aes(label = label), vjust = -1.1, hjust = 0.2, fontface = 'italic') +
  scale_y_continuous(expand = c(0,5)) +
  scale_x_continuous(limits = c(-11,11)) +
  xlab('') +
  ylab('') +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1, lineend = 'square'),
        legend.title = element_blank(),
        aspect.ratio = 1)

# Save plot for figure 5 in paper (Mullins et al. 2026).
ggsave('TS_RNAseq_Volcano_Plot.png', plot, device = 'png', height = 60, width = 60, 
       dpi = 1200, units = 'mm', bg = 'transparent')
