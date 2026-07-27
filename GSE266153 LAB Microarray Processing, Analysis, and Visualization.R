##### Processing and Analysis of GSE266153 Microarray data for Figure 6 #####
# Data: Labyrinth differentiated TS cells.
# Method: DESeq2

#Step 1: Load required packages.
library(readxl) # Used for reading in excel files.
library(writexl) # Used for writing excel files.
library(tidyverse) # Multiple functions used for data processing.
library(DESeq2) # Used to analyze sequencing data.
library(stringr) # Used for character formatting.
library(reshape2) # Used for melting dfs.
library(ggsignif) # Used for adding stats to plots.
library(patchwork) # Used for figure generation.

#Step 2: Set working directory to folder containing data.
project_wd = getwd()
data_wd = paste(project_wd, '/Data/', sep = '')
setwd(data_wd)

#Step 3: Read-in all required data.

#A: Raw microarray data.
data = read_xlsx('GSE266153_raw_counts.xlsx')

#B: Protein-coding and y-linked genes.
Proteins = read_xlsx('ProteinCodingGenes.xlsx') #Protein coding genes
Y_linked = read_xlsx('MouseYlinked.xlsx') #Y-linked genes

#Step 4: Extract and Format Proteins That Are Not Y-linked.
Proteins = (Proteins %>% filter(!Symbol %in% Y_linked$Symbol))$Symbol

#Step 5: Process raw microarray data.

#A: Format data df
data = data %>%
  select(Probe, contains('Vector')) %>%
  group_by(Probe) %>%
  summarise(across(contains('Vector'), mean))

#B: Assign column names.
colnames(data) = c('Gene', 
                   paste(rep('T0', 5), c(1:5), sep = '_'),
                   paste(rep('T3', 5), c(1:5), sep = '_'),
                   paste(rep('T6', 5), c(1:5), sep = '_'))

#C: Filter data for protein coding genes only.
filtered_data = data %>%
  filter(Gene %in% Proteins)

#D: Create df for filtered microarray data.
counts = data.frame(base::sapply(filtered_data[,-c(1)], as.integer), 
                    row.names = filtered_data$Gene)

#Step 6: Formatting metadata for DESeq2 analysis.

# Create Dataframe for Metadata
meta = data.frame(diff_state = as.factor(
  str_split_fixed(colnames(filtered_data)[-c(1)],
                  '_', 2)[,1]),
  row.names = colnames(filtered_data)[-c(1)])

#Step 7: Check to make sure colnames of filtered micrarray data match rownames of meta data.
all(colnames(counts) %in% rownames(meta)) #Should be TRUE
all(colnames(counts) == rownames(meta)) #Should be TRUE

#Step 8: Create DESeq2 object.
deseq = DESeqDataSetFromMatrix(countData = counts, 
                               colData = meta, 
                               design = ~ diff_state)

#Step 9: Perform DESeq2 analysis by using the Wald test.
deseq = DESeq(deseq, test = 'Wald')

#Step 10: Normalize DESeq2 data for downstream visualization.

#A: Perform log2 normalization
log2_norm = rlog(deseq)
log2_norm_data = data.frame(Gene = rownames(assay(log2_norm)), 
                            round(assay(log2_norm), digits = 2))

#B: Set results wd
results_wd = paste(project_wd, '/Results/', sep = '')
setwd(results_wd)

#C: Save normalized data for later use.
write_xlsx(log2_norm_data, 'GSE266153_log2_normalized_counts.xlsx') 

#Step 11: Process and extract results from DESeq2 analysis.

#A: 0-day vs. 3-day results.

# Extraction of results.
results_0v3 = results(deseq, contrast = c('diff_state', 'T3', 'T0'), 
                      alpha = 0.05)

# Formatting results as a df.
de_0v3 = as.data.frame(results_0v3[order(results_0v3$padj),])
de_0v3 = de_0v3 %>%
  mutate(Gene = rownames(de_0v3), Sample = 'T3') %>%
  filter(padj <= 0.05) # Filtering for significant genes 
de_0v3 = de_0v3[,c(7:8,1:6)]

#B: 0-day vs. 6-day results.

# Extraction of results.
results_0v6 = results(deseq, contrast = c('diff_state', 'T6', 'T0'), 
                      alpha = 0.05)

# Formatting results as a df.
de_0v6 = as.data.frame(results_0v6[order(results_0v6$padj),])
de_0v6 = de_0v6 %>%
  mutate(Gene = rownames(de_0v6), Sample = 'T6') %>%
  filter(padj <= 0.05)
de_0v6 = de_0v6[,c(7:8,1:6)]

#C: Combining both results dfs together.
de_both = rbind(de_0v3, de_0v6) %>%
  arrange(padj) # Sort by significance

#D: Save combined df.
write_xlsx(de_both, 'GSE266153_T0vT3vT6_de.xlsx')

# Step 12: Create function for performing stats and displaying data as boxplots.

#A: display_data function
# data = df with normalized data
# gene = target gene for analysis
# stats = df with statistical results (Wald test)
# manual_limits = sets automatic limits for box plots unless numeric vector
# of length 2 (first number is min, second number is max) is assigned.
display_data = function(data, gene, stats, manual_limits = NULL) {
  
  # Generate plot data df.
  plot_data = data %>%
    filter(Gene == gene) %>%
    melt(., variable.name = 'Sample', value.name = 'Expression') %>%
    mutate(Treatment = str_split_fixed(Sample, '_', 2)[,1]) %>%
    mutate(Treatment = factor(Treatment, levels = unique(Treatment))) %>%
    select(Gene, Sample, Treatment, Expression)
  
  # Set min and max limits for box plots.
  if(is.null(manual_limits) == T) {
    min_limit = min(plot_data$Expression) * 0.925
    max_limit = max(plot_data$Expression) * 1.075
  } else {
    min_limit = manual_limits[1]
    max_limit = manual_limits[2]
  }
  
  
  # Generate statistical analysis df.
  stat_df = stats %>%
    filter(Gene == gene)
  
  # Set stat for each comparison to n.s. (not significant)
  t3 = 'n.s.'
  t6 = 'n.s.'
  
  # Use if/else statements to determine statistical significance.
  if('T3' %in% stat_df$Sample) {
    t3 = (stat_df %>% filter(Sample == 'T3'))$padj
    if(t3 < 0.05 & t3 >= 0.01) {
      t3 = '*'
    } else if(t3 < 0.01 & t3 >= 0.001) {
      t3 = '**'
    } else if(t3 < 0.001 & t3 >= 0.0001) {
      t3 = '***'
    } else if(t3 < 0.0001) {
      t3 = '****'
    }
  }
  # If t0 vs. t3 comparison is n.s., then t3 will remain as n.s.
  
  if('T6' %in% stat_df$Sample) {
    t6 = (stat_df %>% filter(Sample == 'T6'))$padj
    if(t6 < 0.05 & t6 >= 0.01) {
      t6 = '*'
    } else if(t6 < 0.01 & t6 >= 0.001) {
      t6 = '**'
    } else if(t6 < 0.001 & t6 >= 0.0001) {
      t6 = '***'
    } else if(t6 < 0.0001) {
      t6 = '****'
    }
  }
  # If t0 vs. t6 comparison is n.s., then t6 will remain as n.s.
  
  # Create annotation vector for plot
  anno = c(t3,t6)
  
  # Generate box plot of data using ggplot2 package.
  plot = ggplot(plot_data, aes(x = Treatment, y = Expression)) +
    stat_boxplot(coef = Inf, 
                 staplewidth = 0.5, 
                 linewidth = 1, 
                 color = 'black',
                 fill = 'white') +
    geom_point(shape = 21,
               fill = 'white', 
               size = 1, 
               stroke = 1,
               position = position_dodge2(width = 0.8),
               show.legend = F) +
    geom_signif(comparisons = list(c('T0','T3'), c('T0','T6')), 
                annotation = anno, size = 1, tip_length = 0, step_increase = 0.2,
                vjust = 0.5) +
    scale_y_continuous(limits = c(min_limit, max_limit)) +
    scale_color_manual(values = c('black','red')) +
    ggtitle(gene) +
    ylab(bquote(Log[2]~Expression)) +
    xlab('') +
    scale_x_discrete(labels = c('0\n-','3\n+','6\n+')) +
    theme_classic() +
    theme(plot.title = element_text(size = 12, color = 'black', hjust = 0.5, 
                                    vjust = 1, face = 'italic'),
          axis.text = element_text(size = 12, color = 'black'),
          axis.title = element_text(size = 12, color = 'black'),
          axis.line = element_line(linewidth = 1, color = 'black', 
                                   lineend = 'square'),
          axis.ticks = element_line(linewidth = 1),
          axis.ticks.length = unit(0.2, 'cm'),
          axis.ticks.x = element_blank(),
          aspect.ratio = 1.5)
  return(plot) # Function returns the plot
}

# Step 13: Display plots for genes used in Figure 6 of paper.
p1 = display_data(log2_norm_data, 'Hspb1', de_both) # Hspb1 boxplot
p2 = display_data(log2_norm_data, 'Hspb8', de_both) # Hspb8 boxplot

# Step 14: Paste both box plots together for figure.
fig = (p1 + theme(axis.title.x = element_blank())) | (p2 + theme(axis.title.y = element_blank()))

# Step 15: Save figure.
ggsave('Hsp_plot.png', fig, device = 'png', scale = 1, width = 4, height = 2.5, 
       dpi = 1200, units = 'in')

# Step 16: Figure details.
#A: Customization of figure occurred in Adobe Photoshop.

#B: Later versions of box plots for figure were set to same min/max limits
# and merged using Adobe photoshop.
