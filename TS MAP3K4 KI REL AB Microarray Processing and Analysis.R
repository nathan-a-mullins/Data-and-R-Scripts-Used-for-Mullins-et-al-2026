##### Visualization of Microarray Data for Mullins et al. 2026 #####
# Average Log2 fold change for treat (TS MAP3K4 KI with REL add back) and 
# control (TS MAP3K4 KI) along with the associated p-values from ANOVA test.

#Step 1: Load required packages.
library(readxl) # for reading in data.
library(tidyverse) # for formatting and displaying data.

#Step 2: Read-in and format data.
data = read_xlsx('GenelistwholeANOVA.xlsx') # Microarray data.

# Formatting df for downstream visualization.
data = data %>%
  select(`Gene Symbol`, 
         `p-value(Attribute)`, 
         `Fold-Change(treat vs. control)...9`) # Selecting relavent columns.

colnames(data) = c('Gene','p.value','fold.change') # Renaming columns.

#Step 3: Determine direction of fold change for each gene based on thresholds.
direction = c() # Assign empty direction vector and fill based on condition.

# Algorithm needed for filling empty direction vector based on conditions below:
# Increased gene threshold = gene log2FC > log2(1.5)
# Decreased gene threshold = gene log2FC < log2(0.75)
# Unchanged (neither) genes = log2(0.75) < gene log2FC < log2(1.5)

# Use for loop with if/else conditions:
for(i in 1:nrow(data)) {
  if(data$fold.change[i] > log2(1.5)) {
    direction = c(direction, 'up')
  } else if(data$fold.change[i] < log2(0.75)) {
    direction = c(direction, 'down')
  } else if (data$fold.change[i] < log2(1.5) & 
             data$fold.change[i] > log2(0.75)) {
    direction = c(direction, 'neither')
  }
}

data$direction = direction # Assign filled direction vector to column in data.

#Step 4: Determine significance of each gene based on thresholds.
significance = c() # Assign empty significance vector and fill based on condition.

# Algorithm needed for filling empty significance vector based on conditions below:
# Significant genes = gene padj (adjusted p-value) < 0.05
# Nonsignificant genes = gene padj (adjusted p-value) > 0.05

# Use for loop with if/else conditions:
for(i in 1:nrow(data)) {
  if(data$p.value[i] < 0.05) {
    significance = c(significance, 'Yes')
  } else {
    significance = c(significance, 'No')
  }
}

data$significance = significance # Assign filled direction vector to column in data.

#Step 5: Generate volcano plot using processed data (Figure 8A).

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

# For loop with if/else conditions used to generate color scale.
for(i in 1:nrow(data)) {
  if(data$direction[i] == 'up' & data$significance[i] == 'Yes') {
    colors = c(colors, 'orangered')
  } else if(data$direction[i] == 'down' & data$significance[i] == 'Yes') {
    colors = c(colors, 'royalblue')
  } else {
    colors = c(colors, 'gray35')
  }
}

data$colors = colors # Assign column in data with filled colors vector.

# Generate volcano plot using ggplot2 package.
plot = ggplot(data, 
              aes(x = fold.change, 
                  y = -log2(p.value))) +
  geom_point(shape = 21, stroke = 0, fill = data$colors) +
  geom_hline(
    yintercept = -log2(0.05), 
    linetype = "solid", 
    color = "black", 
    size = 1,
    alpha = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,30)) +
  scale_x_continuous(limits = c(-60,60)) +
  xlab('') +
  ylab('') +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1, lineend = 'square'),
        legend.title = element_blank(),
        aspect.ratio = 1)

# Save volcano plot.
ggsave('REL_Microarray_Volcano_Plot.png', plot, device = 'png', height = 60, width = 60, 
       dpi = 1200, units = 'mm', bg = 'transparent')

#Step 6: Generate bar plot showing number of DEGs using processed data (Figure 8B).

# Create separate df for increased DEGs.
up.degs = data %>%
  filter(p.value < 0.05 & 2^fold.change > 1.5)

# Create separate df for decreased DEGs.
down.degs = data %>%
  filter(p.value < 0.05 & 2^fold.change < 0.75)

# Format df for bar plot.
direction = c('Up', 'Down') # direction
n.genes = c(nrow(up.degs), nrow(down.degs)) # number of DEGs
n.gene.df = data.frame(direction, n.genes) # df used for plot

# Generate bar plot with number of DEGs using ggplot2 package.
plot = ggplot(n.gene.df, aes(x = factor(direction, levels = c('Up', 'Down')), 
                              y = n.genes)) +
  geom_col(fill = c('orangered','royalblue'), color = 'black', linewidth = 1) +
  scale_y_continuous(limits = c(0, 5000), expand = F) +
  scale_x_discrete(expand = c(1.2,0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1, lineend = 'square'),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        aspect.ratio = 2)

# Save DEG barplot.
ggsave('REL_Microarray_DEG_barplot.png', plot, height = 39, width = 67, 
       dpi = 1200, units = 'mm', bg = 'transparent')

#Step 6: Generate bar plot showing log2FC of Hspb1 and Hspb8 using processed data (Figure 8C-D).

# Create new df by filtering for Hspb1 and selecting relevant columns.
plot.data = data %>%
  filter(Gene %in% c('Hspb1')) %>%
  select(Gene, 
         `p.value`, 
         `fold.change`)

# Formatting plot data for ggplot2.
plot.data$Sample = c('RELAB')

# Creating dummy df for ggplot2.
control.data = data.frame(Gene = c('Hspb1'),
                          p.value = c(1),
                          fold.change = c(1),
                          Sample = c('KI'))

# Combining data df and dummy df into one df (plot data).
plot.data = rbind(plot.data, control.data) 

# Assigning desired aesthetics for geom_bar to bar variable.
bar = geom_bar(stat = 'summary',
               color = 'black',
               lwd = 1,
               position = position_dodge(width = 0.9),
               fill = 'red')

# Assigning desired themes to j22_theme variable.
j22_theme = theme_classic() +
  theme(axis.line = element_line(linewidth = 1, lineend = 'square'),
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(0.15, 'cm'),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        aspect.ratio = 1.75,
        plot.background = element_blank(),
        panel.background = element_blank()) 

# Generating bar plot using ggplot2 package to display Hspb1 data.
plot = ggplot(plot.data, 
              aes(x = factor(Sample, levels = c('KI','RELAB')), 
                  y = fold.change,
                  group = 1)) +
  bar +
  scale_x_discrete(expand = c(1.2,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = 'on', ylim = c(0,4)) +
  scale_fill_manual(values = c('red')) +
  guides(fill = 'none') +
  j22_theme

# Saving Hspb1 bar plot.
ggsave('Rel Hspb1 panel.png', plot, device = 'png', 
       height = 45, width = 45, dpi = 1200, units = 'mm', bg = 'transparent')

# Create new df by filtering for Hspb8 and selecting relevant columns.
plot.data = data %>%
  filter(Gene %in% c('Hspb8')) %>%
  select(Gene, 
         `p.value`, 
         `fold.change`)

# Formatting plot data for ggplot2.
plot.data$Sample = c('RELAB')

# Creating dummy df for ggplot2.
control.data = data.frame(Gene = c('Hspb8'),
                          p.value = c(1),
                          fold.change = c(1),
                          Sample = c('KI'))

# Combining data df and dummy df into one df (plot data).
plot.data = rbind(plot.data, control.data)

# Generating bar plot using ggplot2 package to display Hspb8 data.
plot = ggplot(plot.data, 
              aes(x = factor(Sample, levels = c('KI','RELAB')), 
                  y = fold.change,
                  group = 1)) +
  bar +
  scale_x_discrete(expand = c(1.2,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = 'on', ylim = c(0,10)) +
  scale_fill_manual(values = c('red')) +
  guides(fill = 'none') +
  j22_theme

# Saving Hspb8 bar plot.
ggsave('Rel Hspb8 panel.png', plot, device = 'png', 
       height = 45, width = 45, dpi = 1200, units = 'mm', bg = 'transparent')
