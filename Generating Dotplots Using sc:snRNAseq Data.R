##### Functions for Displaying sc/snRNAseq Data as Dotplots or Heatmaps #####

#Step 1: Load required packages.
library(tidyverse) # Used for data processing in functions.
library(ggplot2) # Used to plot data.

#Step 2: Read-in pre-formatted data from Jiang et al. study (GSE156125).
data = read.csv('jiang_average_melted.csv') # scRNAseq data used in Figure 6.

# Data format:
# Data: melted df of sc/snRNAseq data with following columns:
# Column 1: Cell Type
# Column 2: Gestational age (GA)
# Column 3: Gene
# Column 4: Mean of CLR expression across all cells within cell type and GA
# Column 5: Sum of CLR expression across all cells within cell type and GA
# Column 6: Median of CLR expression across all cells within cell type and GA
# Column 7: Percentage of cells within cell type and GA that express gene
# Column 8: Total amount of cells within cell type and GA
# Most sc/snRNAseq datasets can be formatted in this way.

#Step 3: Create function for processing and displaying data.
single_gene_dotplot = function(data, gene, stat = 'mean', color_set = 'inferno') {
  
  # Generate color scale for dot fills.
  if(color_set == 'inferno') { # Custom inferno color scale
    colors = c("#000004", # Highest expression.
               "#2b0c41",
               "#6a0c0b",
               "#781c3d",
               "#a73c60",
               "#cf4446",
               "#ed6925",
               "#fb9b06",
               "#f7d13d",
               "#fcffa4") # Lowest expression.
  } else if(color_set == 'viridis') { # Custom viridis color set
    colors = c("#440154", # Highest expression.
               "#482878",
               "#3e4989",
               "#31688e",
               "#26828e",
               "#1f9e89",
               "#35b779",
               "#6ece58",
               "#b5de2b",
               "#fde725") # Lowest expression.
  }
  
  # Reverse colors for plot.
  colors = rev(colors)
  
  # Create a nested function to modify alpha of color scale.
  modify_alpha = function(colors) {
    adj_colors = c() # Empty vector for color alpha adjustment.
    alpha = 1/length(colors) # Base alpha adjustment function.
    for(i in 1:length(colors)) { # Changes alpha based on low to high scale.
      adj_colors = c(adj_colors, adjustcolor(colors[i],alpha.f=alpha)) # Adjusting color based on alpha function.
      alpha = alpha + 1/length(colors) # Modify alpha function for next color.
    }
    return(adj_colors) # Returns adjusted color set.
  }
  
  # Creating adjusted color set for dot plot.
  new_colors = modify_alpha(colors)
  
  # Creating vector with original names of cell types in Marsh and Jiang datasets.
  if('TSC/ExE' %in% data$Cell_Type) { # Jiang data has TS cells.
    cells = c('TSC/ExE','LaTP','SynTII P','LaTP 2','SynTI P','S-TGC P','S-TGC',
              'EPC','EPC M','Gly-T','SpA-TGC','SpT','2* P-TGC P',
              '2* P-TGC','1* P-TGC')
    
    # Filtering target gene, cell types, and cell types with > than 5 total cells for Jiang data.
    df = data %>%
      filter(Gene == gene & Cell_Type %in% cells & Total_Cells > 5)
    
  } else { # Marsh data has these cell types.
    cells = c('LaTP','SynTII.Precursor','SynTII',
              'LaTP.2','SynTI.Precursor','SynTI')
    
    # Filtering target gene and cell types for Marsh data.
    df = data %>%
      filter(Gene == gene & Cell_Type %in% cells) 
    
    # Modifying names of cell types in df for simplicity and compatibility between the datasets.
    df$Cell_Type[df$Cell_Type == 'SynTII.Precursor'] = 'SynTII P'
    df$Cell_Type[df$Cell_Type == 'LaTP.2'] = 'LaTP 2'
    df$Cell_Type[df$Cell_Type == 'SynTI.Precursor'] = 'SynTI P'
    
    # Modifying names of cell types in cells vector for simplicity and compatibility between the datasets.
    cells[cells == 'SynTII.Precursor'] = 'SynTII P'
    cells[cells == 'LaTP.2'] = 'LaTP 2'
    cells[cells == 'SynTI.Precursor'] = 'SynTI P'
  }
  
  if(stat == 'mean') { # Argument for plotting mean expression.
    max_value = max(df$Mean_Expression) # max of mean expression across dataset.
    
    # Plot mean expression in df using ggplot2 package.
    plot = ggplot(df, aes(x = GA, # GA on x-axis.
                          y = factor(Cell_Type, levels = rev(cells)), # Cell types on y-axis.
                          fill = Mean_Expression, # Fill is scaled to mean expression.
                          size = Positive)) # Size is scaled to percent positive.
    
  } else if(stat == 'sum') { # Argument for plotting sum expression.
    max_value = max(df$Sum_Expression) # max of sum expression across dataset.
    
    # Plot sum expression in df using ggplot2 package.
    plot = ggplot(df, aes(x = GA, # GA on x-axis.
                          y = factor(Cell_Type, levels = rev(cells)), # Cell types on y-axis.
                          fill = Sum_Expression, # Fill is scaled to sum expression.
                          size = Positive)) # Size is scaled to percent positive.
    
  } else if(stat == 'median') { # Argument for plotting median expression
    max_value = max(df$Median_Expression) # max of median expression across dataset.
    
    # Plot median expression in df using ggplot2 package.
    plot = ggplot(df, aes(x = GA, # GA on x-axis.
                          y = factor(Cell_Type, levels = rev(cells)), # Cell types on y-axis.
                          fill = Median_Expression, # Fill is scaled to median expression.
                          size = Positive)) # Size is scaled to percent positive.
  }
  
  # Finish formatting the plot (same for mean, sum, and median expression).
  plot = plot +
    geom_point(shape = 21, stroke = 0) + # Formatting points ("dots").
    scale_fill_gradientn(colors = new_colors,
                         limits = c(0, max_value)) + # Formatting fill (expression) gradient with (0, max value) limits.
    guides(size = guide_legend(override.aes = list(fill = 'black'))) + # Formatting legend.
    theme_bw() + # Black and white theme
    theme(aspect.ratio = length(unique(df$Cell_Type))/length(unique(df$GA)), # Custom aspect ratio based on length of cell types and GA.
          axis.text = element_text(size = 18, color = 'black', face = 'bold'), # Formatting axis text.
          axis.text.x = element_text(size = 18, angle = 90, hjust = 1, 
                                     vjust = 0.5, color = 'black', face = 'bold'), # Formatting x-axis text.
          axis.text.y = element_text(size = 18, color = 'black', face = 'bold'), # Formatting y-axis text.
          axis.ticks = element_blank(), # Removing axis ticks.
          axis.title = element_blank(), # Removing axis titles.
          legend.key = element_rect(fill = "transparent", colour = NA), # Remove legend borders.
          legend.position = 'right', # Assigning legend position to right side.
          legend.text = element_text(size = 12, color = 'black', face = 'bold'), # Formatting legend text.
          legend.title = element_text(size = 18, color = 'black', face = 'bold'), # Formatting legend titles.
          panel.border = element_rect(color = "black", fill = NA, size = 3), # Formatting panel border.
          panel.grid = element_blank(), # Removing panel grid.
          plot.background = element_blank(), # Removing plot background.
          plot.title = element_text(hjust = 0.5, size = 24, color = 'black', face = 'bold.italic')) + # Formatting plot title.
    ggtitle(gene) # Setting plot title as target gene.
  return(plot) # Function returns the plot.
}

#Step 4: Create Hspb1 dotplot for Figure 6 Panel M in Mullins et al. 2026.
p1 = single_gene_dotplot(data, 'Hspb1') # Hspb1 dotplot using GSE156125 dataset.
ggsave('Hspb1.png', p1, device = 'png', dpi = 1200, height = 5, width = 7) # Save dotplot.

#Step 5: Create Hspb8 dotplot for Figure 6 Panel N in Mullins et al. 2026.
p2 = single_gene_dotplot(data, 'Hspb8') # Hspb8 dotplot using GSE156125 dataset.
ggsave('Hspb8.png', p2, device = 'png', dpi = 1200, height = 5, width = 7) # Save dotplot.

