# Packages ----------------------------------------------------------------

library(corrplot)
library(stringr)
library(microbiome)
library(hrbrthemes)
library(ggrepel)
library(rstatix)
library(doParallel)
library(DECIPHER); packageVersion("DECIPHER") 
library(phyloseq); packageVersion("phyloseq") 
library(Biostrings); packageVersion("Biostrings") 
library(ggplot2); packageVersion("ggplot2")
library(tidyverse); packageVersion("tidyverse") 
library(readxl)
library(pals)
library(ggpubr)
library(cowplot)
library(multcompView)
library(microbiomeMarker)
library(pheatmap)
library(pairwiseAdonis)
library(gplots)
library("SIAMCAT")
library(rlang)
library(DESeq2); packageVersion("DESeq2") 
library(vegan) 
library(magrittr)
library(dplyr)
library(ape)
library(biomformat); packageVersion("biomformat") 
library(ggforce)
library(Maaslin2); packageVersion("Maaslin2") 
library(reshape2)
library(picante)
library(rstatix)
library(ggpubr)
library(survival)
library(palmerpenguins)
library(multcomp)
library(agricolae)
library(knitr)
library(tibble)
library(dplyr)
library(RColorBrewer)
library(ggside)
library("ape")
library(viridis)
library(scales)
library(ggsignif)
library(metagenomeSeq); packageVersion("metagenomeSeq")  
library(dendextend); packageVersion("dendextend")   
library(rms); packageVersion("rms")
library(breakaway); packageVersion("breakaway")   
library(ggprism)
library(ggExtra)
library(FEAST)
library(phangorn); packageVersion("phangorn")
library(ggplot2)
library(dplyr)
library(ggpicrust2)
library(readr)
library(cowplot)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

# Import data and process -------------------------------------------------
# Read phyloseq object

outdir <- ""

data <- readRDS("dada2_phyloseq_Q30_no402.rds")

# Store metadata
metadata <- read.table("~/bioinformatics/oral_meta_datasets/oral_metadata_no402.tsv", sep = "\t", header = T, row.names = 1)

# create metadata phyloseq object
metadata <- sample_data(metadata)

# Update the metadata in the phyloseq object
sample_data(data.f) <- metadata

# Remove sample with <1000 reads
low_read_1000 <- c("CO40", "MGM23427")
data.f <- prune_samples(!(sample_names(data.f) %in% low_read_1000), data.f)

# Collapse at the genus level
data.f <- tax_glom(data.f, taxrank = "Genus")


# Demo barplot for all studies ---------------------------------------------

metadata %>%
  group_by(continent, study_group_name) %>%
  summarise(count = n()) %>%
  arrange(continent, desc(count))

facet_order <- c("study_group_name", "sample_type", "V_region", "continent", "Population")

# Reshape data to long format for faceting
metadata_long <- metadata %>%
  pivot_longer(cols = all_of(facet_order), names_to = "Facet_Variable", values_to = "Facet_Value") %>%
  mutate(
    Facet_Variable = factor(Facet_Variable, levels = facet_order),  # Order facets
    Facet_Value = fct_infreq(Facet_Value)  # Order bars within each facet by frequency (largest to smallest)
  )

# Create the bar plot faceted in one row
demo <- ggplot(metadata_long, aes(x = Facet_Value)) +
  geom_bar(width = 0.7, fill = "blue", color = "black", size = 1.2, position = position_dodge2(preserve = "single")) +  
  facet_grid(~ Facet_Variable, scales="free_x", space = "free") +  # All facets in one row
  theme_classic() +
  labs(x = NULL, y = "Number of samples") +  # Adjust labels
  theme_cowplot(20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +  # Minimal padding on y-axis
  coord_cartesian(ylim = c(0, max(table(metadata_long$Facet_Value))))  # Adjust y-axis limits

ggsave(filename = file.path(outdir, "fig1_demographs.png"), plot = demo, width = 14, height = 8, dpi = 300)

# Specify phyloseq objects for comparisons --------------------------------
# Store absolute abundances object
abs = data.f

# Convert into total sum scaling (relative abundance)
rel = transform_sample_counts(abs, function(x) x / sum(x))
# Remove taxa occuring <0.1%
rel <- filter_taxa(rel, function(x) sum(x) > .001, TRUE)
# Convert (AGAIN) into total sum scaling (relative abundance)
rel <- transform_sample_counts(rel, function(x) x / sum(x))

# sgn
sgn_rel <- subset_samples(rel, study_group_name=="CD" | study_group_name=="UC" | study_group_name=="healthy_control")

# Adults
rel_ad <- subset_samples(rel, study_alias != "SRP385133")
# Adults SGN
sgn_rel_ad <- subset_samples(rel_ad, study_group_name=="CD" | study_group_name=="UC" | study_group_name=="healthy_control")


# Bar plots ---------------------------------------------------------------
# ALL 
## Agglomerate at the genus level 
rel_ph <- tax_glom(sgn_rel, taxrank = 'Phylum')
rel_ph_df <- psmelt(rel_ph)

# Reorder Phylum based on total Abundance
rel_ph_df <- rel_ph_df %>%
  mutate(Phylum = fct_reorder(Phylum, Abundance, .fun = sum, .desc = TRUE))

# Select colors from Set1 and combine with Set3
set1_colors <- brewer.pal(8, "Set1")  # Maximum 9 colors
set3_colors <- brewer.pal(12, "Set3")[1:6]  # Select the first 3 colors from Set3

# Combine colors to create a custom palette
custom_colors <- c(set1_colors, set3_colors)

## Plot by AD status (abundances summarized by group)
rel_ph_SGN_plot.ap <- ggplot(rel_ph_df, aes(x=study_group_name, y=Abundance, fill=Phylum))+
  geom_bar(stat='identity', position='fill')+
  theme_classic()+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3))+
  #scale_x_discrete(labels=c('adult', 'pediatric')) +
  scale_fill_manual(values = custom_colors)
  #+ facet_wrap(~study_group_name)
 
# Construct the filename dynamically
output_filename <- sprintf("Figures/one/barplot_phy_sgn.png")

# Save the plot
ggsave(
  filename = output_filename,
  plot = rel_ph_SGN_plot.ap,
  width = 6.5,
  height = 7.5,
  dpi = 300
)


# FUNCTIONS ---------------------------------------------------------------
# Alpha diversity
generate_boxplot <- function(data, variable, population_col = "Population") {
  # Perform t-test
  t_test_result <- t.test(as.formula(paste(variable, "~", population_col)), data = data)
  
  # Extract p-value
  p_value <- t_test_result$p.value
  
  # Determine significance label
  significance_label <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ "ns"
  )
  
  # Create the boxplot
  plot <- ggplot(data, aes_string(x = population_col, y = variable, fill = population_col)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, size = 0.9) +
    geom_jitter(aes_string(color = population_col), size = 3, width = 0.1, height = 0, alpha = 0.7) +
    theme_cowplot(20) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(size = 0.5, color = 'grey'),
      axis.line = element_line(color = "black", size = 0.8),
      axis.title.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +  # Minimal padding
    coord_cartesian(ylim = c(0, max(data[[variable]]))) +  # Ensure no extra space below 0
    # Add significance label
    geom_text(
      aes(
        x = 1.5,
        y = max(data[[variable]]) - 0.05 * max(data[[variable]]),
        label = significance_label
      ),
      vjust = -0.5,
      size = 10
    )
  
  return(list(plot = plot, p_value = p_value, significance_label = significance_label))
}
generate_faceted_boxplot <- function(data, variable, grouping_col, facet_col, p_adjust_method = "BH") {
  # Step 1: Perform t-tests and adjust p-values
  t_test_results <- data %>%
    group_by(!!sym(facet_col)) %>%
    summarise(
      t_test = list(t.test(as.formula(paste(variable, "~", grouping_col)), data = cur_data())),  # Perform t-test
      .groups = "drop"
    ) %>%
    mutate(
      p_value = map_dbl(t_test, ~ .x$p.value),  # Extract p-values
      BH_corrected_p = p.adjust(p_value, method = p_adjust_method)  # Adjust p-values
    )
  
  # Step 2: Merge results back to the original data for annotation
  data_with_p <- data %>%
    left_join(t_test_results, by = facet_col)
  
  # Step 3: Create the faceted boxplot
  plot <- ggplot(data_with_p, aes_string(x = grouping_col, y = variable, fill = grouping_col)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, size = 0.9) +
    geom_jitter(aes_string(color = grouping_col), size = 3, width = 0.1, height = 0, alpha = 0.7) +
    theme_cowplot(20) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(size = 0.5, color = 'grey'),
      axis.line = element_line(color = "black", size = 0.8),
      axis.title.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +  # Minimal padding
    coord_cartesian(ylim = c(0, max(data[[variable]]))) +  # Ensure no extra space below 0
    # Add significance labels
    geom_text(
      data = t_test_results,
      aes(
        x = 1.5,
        y = max(data[[variable]]) - 0.05 * max(data[[variable]]),  # Calculate position
        label = case_when(
          BH_corrected_p < 0.001 ~ "***",
          BH_corrected_p < 0.01 ~ "**",
          BH_corrected_p < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      ),
      vjust = -0.5,
      size = 10,
      inherit.aes = FALSE  # Prevent ggplot from using default plot aesthetics
    ) +
    facet_wrap(as.formula(paste("~", facet_col)))  # Facet by the specified column
  
  return(list(plot = plot, t_test_results = t_test_results))
}

generate_boxplot_anova <- function(data, variable, grouping_col = "Population") {
  # Validate inputs
  if (!all(c(grouping_col, variable) %in% names(data))) {
    stop("Specified columns are not present in the data.")
  }
  
  # Perform ANOVA
  aov_result <- aov(as.formula(paste(variable, "~", grouping_col)), data = data)
  summary_aov <- summary(aov_result)
  
  # Perform Tukey's HSD test
  tukey_result <- TukeyHSD(aov_result)
  
  # Extract Tukey results into a dataframe and ensure ordering matches the comparisons
  tukey_df <- as.data.frame(tukey_result[[grouping_col]]) %>%
    mutate(
      comparison = rownames(.), 
      group1 = sub("-.*", "", comparison),  # Extract first group
      group2 = sub(".*-", "", comparison)   # Extract second group
    ) %>%
    arrange(group1, group2) %>%  # Ensure ordering matches comparison_list
    mutate(
      significance_label = case_when(
        `p adj` < 0.001 ~ "***",
        `p adj` < 0.01  ~ "**",
        `p adj` < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  # Generate y_positions for each comparison to avoid overlap
  max_value <- max(data[[variable]], na.rm = TRUE)
  
  # Ensure the comparison list matches Tukey's ordering
  comparison_list <- split(tukey_df[, c("group1", "group2")], seq(nrow(tukey_df)))
  comparison_list <- lapply(comparison_list, as.character)  # Ensure each pair is a character vector
  
  annotations <- tukey_df$significance_label  # Now aligned with comparisons
  
  y_positions <- seq(max_value + 0.1 * max_value, by = 0.1 * max_value, length.out = length(comparison_list))
  
  # Create the boxplot
  plot <- ggplot(data, aes_string(x = grouping_col, y = variable, fill = grouping_col)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, size = 0.9) +
    geom_jitter(aes_string(color = grouping_col), size = 3, width = 0.1, height = 0, alpha = 0.7) +
    theme_cowplot(20) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(size = 0.5, color = 'grey'),
      axis.line = element_line(color = "black", size = 0.8),
      axis.title.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    # Add pairwise significance annotations with corrected ordering
    ggsignif::geom_signif(
      comparisons = comparison_list,
      map_signif_level = FALSE,
      annotations = annotations,  # Use aligned significance labels
      y_position = y_positions,
      textsize = 10
    )
  
  return(list(
    plot = plot,
    anova_summary = summary_aov,
    tukey_result = tukey_result,
    tukey_df = tukey_df
  ))
}

# Beta diversity
# ALL
plot_beta_diversity <- function(physeq_object, comparator) {  
  # Step 1: Prune taxa
  physeq_object <- prune_taxa(taxa_sums(physeq_object) > 0, physeq_object)
  
  # Step 2: Ordination with Bray-Curtis distance
  bray <- ordinate(physeq_object, method = "PCoA", distance = "bray")
  
  # Step 3: Extract variance explained for the first two axes
  variance_axis1 <- bray$values$Relative_eig[1] * 100  # Percentage for Axis 1
  variance_axis2 <- bray$values$Relative_eig[2] * 100  # Percentage for Axis 2
  
  # Format axis labels with variance explained
  axis1_label <- sprintf("PCo1 (%.1f%%)", variance_axis1)
  axis2_label <- sprintf("PCo2 (%.1f%%)", variance_axis2)
  
  # Step 4: ADONIS test
  arc_dist_matrix <- phyloseq::distance(physeq_object, method = "bray")
  adonis_result <- vegan::adonis2(as.formula(paste("arc_dist_matrix ~", comparator)), 
                                  data = as(sample_data(physeq_object), "data.frame"))
  # Extract p-value from ADONIS test
  p_value <- adonis_result$`Pr(>F)`[1]
  title_text <- sprintf("PERMANOVA p = %.3f", p_value)
  
  print(adonis_result)
  
  # Optional: Pairwise ADONIS test
  pairwise_result <- pairwise.adonis(arc_dist_matrix, sample_data(physeq_object)[[comparator]])
  print(pairwise_result)
  
  # Step 5: Plot ordination
  p.s.nmds <- plot_ordination(physeq_object, bray, color = comparator, 
                              title="\nPERMANOVA p = 0.034")  # Customize the title if needed
  
  # Step 6: Customize the plot with dynamic axis labels
  beta <- p.s.nmds +
    geom_point(size = 2) +
    stat_ellipse(aes(group = !!sym(comparator)), linetype = 1) +
    theme_cowplot(20) +
    coord_cartesian(clip = "off") +
    labs(x = axis1_label, y = axis2_label, title = title_text) +  # Use dynamic labels for axes
    theme(
      legend.position = c(0.99, 0.99),  # Legend inside the plot box (top-right corner)
      legend.justification = c(1, 1),   # Align the legend by its top-right corner # Increase the legend text size
      legend.title = element_blank(),
      panel.border = element_rect(color = "black", size = 1, fill = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(
        fill = alpha("white", 0.7),  # Transparent white background
        color = "black",             # Border color
        size = 0.5                   # Border thickness
      ))
  
  # Step 7: Add marginal plots
  beta2 <- ggMarginal(beta, type = 'boxplot', groupFill = TRUE, size = 10)
  
  # Return the plot
  return(beta2)
}

# SGN
plot_beta_diversity_faceted <- function(physeq_object, comparator, facet_variable) {
  # Step 1: Prune taxa
  physeq_object <- prune_taxa(taxa_sums(physeq_object) > 0, physeq_object)
  
  # Step 2: Ordination with Bray-Curtis distance
  bray <- ordinate(physeq_object, method = "PCoA", distance = "bray")
  
  # Step 3: Extract variance explained for the first two axes
  variance_axis1 <- bray$values$Relative_eig[1] * 100  # Percentage for Axis 1
  variance_axis2 <- bray$values$Relative_eig[2] * 100  # Percentage for Axis 2
  
  # Format axis labels with variance explained
  axis1_label <- sprintf("PCo1 (%.1f%%)", variance_axis1)
  axis2_label <- sprintf("PCo2 (%.1f%%)", variance_axis2)
  
  # Step 4: ADONIS test
  arc_dist_matrix <- phyloseq::distance(physeq_object, method = "bray")
  adonis_result <- vegan::adonis2(as.formula(paste("arc_dist_matrix ~", comparator)), 
                                  data = as(sample_data(physeq_object), "data.frame"))
  print(adonis_result)
  
  # Optional: Pairwise ADONIS test
  pairwise_result <- pairwise.adonis(arc_dist_matrix, sample_data(physeq_object)[[comparator]])
  print(pairwise_result)
  
  # Extract p-value from ADONIS test
  p_value <- adonis_result$`Pr(>F)`[1]
  title_text <- sprintf("PERMANOVA p = %.3f", p_value)
  
  # Step 5: Plot ordination
  p.s.nmds <- plot_ordination(physeq_object, bray, color = comparator)
  
  # Step 6: Customize the plot with faceting and dynamic axis labels
  facetted_plot <- p.s.nmds +
    geom_point(size = 2) +
    stat_ellipse(aes(group = !!sym(comparator)), linetype = 1) +
    theme_cowplot(20) +
    coord_cartesian(clip = "off") +
    labs(
      x = axis1_label,
      y = axis2_label
    ) +
    facet_wrap(vars(!!sym(facet_variable)), scales = "fixed", ncol = 3) +  # Use fixed scales
    theme(
      legend.position = c(0.99, 0.99),  # Legend inside the plot box (top-right corner)
      legend.justification = c(1, 1),   # Align the legend by its top-right corner
      legend.text = element_text(size = 16),  # Increase the legend text size
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA), 
      # Remove the legend title
      legend.background = element_rect(
        fill = alpha("white", 0.7),  # Transparent white background
        color = "black",             # Border color
        size = 0.5                   # Border thickness
      ),
      plot.title = element_text(
        size = 16,                   # Title size
        margin = margin(t = 10, b = 5)  # Add space above (t) and below (b)
      ),
      strip.text = element_text(size = 20),  # Style for facet labels
      strip.background = element_rect(fill = "grey80", color = "grey80"),  # Standard grey background
      panel.spacing = unit(1, "lines")  # Add space between panels
    )
  
  # Return the facetted plot
  return(facetted_plot)
}
adonis_with_subsets <- function(physeq_object, comparator, p_adjust_method = "fdr") {
  # Define the subgroups for subsetting
  groups <- list(
    "Healthy Control" = subset_samples(physeq_object, IBD_group_name == "healthy_control"),
    "CD" = subset_samples(physeq_object, study_group_name == "CD"),
    "UC" = subset_samples(physeq_object, study_group_name == "UC")
  )
  
  # Initialize a list to store results
  adonis_results_list <- list()
  
  # Loop through each group and perform ADONIS
  for (group_name in names(groups)) {
    subset_physeq <- groups[[group_name]]  # Get the subset for the current group
    
    # Calculate Bray-Curtis distance
    dist_matrix <- phyloseq::distance(subset_physeq, method = "bray")
    
    # Perform ADONIS test
    adonis_result <- vegan::adonis2(
      as.formula(paste("dist_matrix ~", comparator)),
      data = as(sample_data(subset_physeq), "data.frame")
    )
    
    # Store the result
    adonis_results_list[[group_name]] <- data.frame(
      Group = group_name,
      Pr_F = adonis_result$`Pr(>F)`[1],
      R2 = adonis_result$R2[1],
      Note = "ADONIS performed successfully"
    )
  }
  
  # Combine all results into a single data frame
  adonis_results_df <- do.call(rbind, adonis_results_list)
  
  # Adjust p-values
  valid_p_values <- !is.na(adonis_results_df$Pr_F)  # Identify valid p-values
  adonis_results_df$Adjusted_P <- NA  # Initialize column for adjusted p-values
  adonis_results_df$Adjusted_P[valid_p_values] <- p.adjust(
    adonis_results_df$Pr_F[valid_p_values], 
    method = p_adjust_method
  )
  
  return(adonis_results_df)
}



# ALL samples -------------------------------------------------------------
# Alpha
alpha <- microbiome::alpha(rel, index = "all")
merge <- merge(metadata, alpha, by = "row.names")
alpha.2 <- microbiome::alpha(sgn_rel, index = "all")
mergey <- merge(metadata, alpha.2, by = "row.names")

# Generate boxplot
all_IBD_boxplot_shannon <- generate_boxplot(data = merge, variable = "diversity_shannon", population_col = "IBD_group_name")
all_IBD_boxplot_obs <- generate_boxplot(data = merge, variable = "observed", population_col = "IBD_group_name")

# Generate the boxplot with ANOVA and Tukey's HSD test
all_SGN_boxplot_shannon <- generate_boxplot_anova(data = mergey, variable = "diversity_shannon", grouping_col = "study_group_name")
all_SGN_boxplot_obs <- generate_boxplot_anova(data = mergey, variable = "observed", grouping_col = "study_group_name")

# Generate faceted boxplot
all_boxplot_faceted <- generate_faceted_boxplot(data = mergey, variable = "diversity_shannon", grouping_col = "Population", facet_col = "study_group_name")

# Beta
all_beta_IBD_plot <- plot_beta_diversity(rel, "IBD_group_name")
all_beta_sgn_plot <- plot_beta_diversity(sgn_rel, "study_group_name")
# SGN
AP_beta_sgn_plot <- plot_beta_diversity_faceted(
  physeq_object = sgn_rel,
  comparator = "Population",
  facet_variable = "study_group_name"
)
# Run the function with default p-value adjustment (FDR)
AP_sgn_adonis_results <- adonis_with_subsets(
  physeq_object = sgn_rel,
  comparator = "study_group_name"
)

# Fig 2
rel_ph_SGN_plot.ap
all_IBD_boxplot_shannon <- all_IBD_boxplot_shannon$plot
all_IBD_boxplot_obs <- all_IBD_boxplot_obs$plot
rel_ph_SGN_plot 
all_SGN_boxplot_shannon <- all_SGN_boxplot_shannon$plot
all_SGN_boxplot_obs <- all_SGN_boxplot_obs$plot
all_beta_IBD_plot
all_beta_sgn_plot

# Combine the top row
fig_2_IBD.u <- plot_grid(
  rel_ph_SGN_plot.ap,
  all_beta_IBD_plot,
  labels = c('A', 'B'), 
  ncol = 2,
  label_size = 18,
  axis = "tblr",
  align = "h",
  rel_widths = c(1,1.5)
)

# Combine the top row
fig_2_IBD.l <- plot_grid(
  all_IBD_boxplot_shannon,
  all_IBD_boxplot_obs,
  labels = c('C', 'D'), 
  ncol = 2,
  label_size = 18,
  axis = "tblr",
  align = "h"
)

fig2_ibd <- plot_grid(fig_2_IBD.u, fig_2_IBD.l, ncol = 1)

fig_2_sgn.u <- plot_grid(
  rel_ph_SGN_plot.ap,
  all_beta_sgn_plot,
  labels = c('A', 'B'), 
  ncol = 2,
  label_size = 18,
  axis = "tblr",
  align = "h",
  rel_widths = c(1,1.5)
)

fig_2_sgn.l <- plot_grid(
  all_SGN_boxplot_shannon,
  all_SGN_boxplot_obs,
  labels = c('C', 'D'), 
  ncol = 2,
  label_size = 18,
  axis = "tblr",
  align = "h"
)

fig2_sgn <- plot_grid(fig_2_sgn.u, fig_2_sgn.l, ncol = 1, rel_heights = c(1,1.2))

ggsave(filename = file.path(outdir, "fig2_IBD_all.png"), plot = fig2_ibd, width = 12.5, height = 15.5, dpi = 300)

ggsave(filename = file.path(outdir, "fig2_sgn_all.png"), plot = fig2_sgn, width = 12.5, height = 17, dpi = 300)


# Adults vs peads ---------------------------------------------------------
# Adults vs pediatric
# alpha
alpha <- microbiome::alpha(rel, index = "all")
merge <- merge(metadata, alpha, by = "row.names")
alpha.2 <- microbiome::alpha(sgn_rel, index = "all")
mergey <- merge(metadata, alpha.2, by = "row.names")

# Generate boxplot
AP_boxplot_shannon <- generate_boxplot(data = merge, variable = "diversity_shannon", population_col = "Population")
AP_boxplot_obs <- generate_boxplot(data = merge, variable = "observed", population_col = "Population")
# Generate faceted boxplot
AP_boxplot_faceted <- generate_faceted_boxplot(data = mergey, variable = "diversity_shannon", grouping_col = "Population", facet_col = "study_group_name")

# Beta
AP_beta_all_plot <- plot_beta_diversity(rel, "Population")
## SGN
AP_beta_sgn_plot <- plot_beta_diversity_faceted(
  physeq_object = sgn_rel,
  comparator = "Population",
  facet_variable = "study_group_name"
)
# Run the function with default p-value adjustment (FDR)
AP_sgn_adonis_results <- adonis_with_subsets(
  physeq_object = sgn_rel,
  comparator = "Population"
)


AP_boxplot_shannon <- AP_boxplot_shannon$plot
AP_boxplot_obs <- AP_boxplot_obs$plot
AP_boxplot_faceted <- AP_boxplot_faceted$plot
AP_beta_all_plot
AP_beta_sgn_plot
AP_sgn_adonis_results

## COWPLOT

fig_2.1 <- plot_grid(
  rel_ph_SGN_plot.ap,
  AP_boxplot_shannon,
  AP_boxplot_obs,
  labels = c('B', 'C', 'D'), 
  label_size = 18,
  ncol = 3, 
  axis = "tblr",
  align = "h"
)

fig_2.2 <- plot_grid(
  AP_beta_all_plot,
  AP_beta_sgn_plot, 
  labels = c('A', 'B', 'C', 'D'), 
  label_size = 18,
  ncol = 2, 
  axis = "tblr",
  rel_widths = c(1, 1),
  align = "v"
)

fig_2.3 <- plot_grid(
  AP_beta_sgn_plot,
  rel_ph_SGN_plot.ap,
  AP_boxplot_faceted, 
  labels = c('A', 'B', 'C'), 
  label_size = 18,
  ncol = 1, 
  axis = "tblr",
  align = "v"
)

AP_beta_all_plot <- plot_grid(
  AP_beta_all_plot, 
  label_size = 18,
  labels = c('A'))
  
ggsave(filename = file.path(outdir, "adults_peads_beta_all.png"), plot = AP_beta_all_plot, width = 8, height = 8, dpi = 300 )
ggsave(filename = file.path(outdir, "adults_peads_alpha_all.png"), plot = fig_2.1, width = 14, height = 7.5, dpi = 300 )

ggsave(filename = file.path(outdir, "adults_peads_SGN_facets.png"), plot = fig_2.3, width = 15, height = 25, dpi = 500 )

# Combine the top row
facet_all <- plot_grid(
  facetted_plot.a, 
  AP_beta_sgn_plot,
  labels = c('A', 'B'), 
  ncol = 1, 
  label_size = 18,
  align = "v",               # Align vertically
  axis = "l"                # Align the left and right axes
)

ggsave(filename = file.path(outdir, "adults_peads_facet_beta_alpha.pdf"), plot = facet_all, width = 14, height = 15, dpi = 300)


# Adults continent -------------------------------------------------------------
# Continent
# alpha
alpha <- microbiome::alpha(rel_ad, index = "all")
merge <- merge(metadata, alpha, by = "row.names")
alpha.2 <- microbiome::alpha(sgn_rel_ad, index = "all")
mergey <- merge(metadata, alpha.2, by = "row.names")

# Generate boxplot
con_boxplot_shannon <- generate_boxplot(data = merge, variable = "diversity_shannon", population_col = "continent")
con_boxplot_obs <- generate_boxplot(data = merge, variable = "observed", population_col = "continent")
# Generate faceted boxplot
con_boxplot_faceted <- generate_faceted_boxplot(data = mergey, variable = "diversity_shannon", grouping_col = "continent", facet_col = "study_group_name")

# Beta
con_beta_all_plot <- plot_beta_diversity(rel_ad, "continent")
# SGN
con_beta_sgn_plot <- plot_beta_diversity_faceted(
  physeq_object = sgn_rel_ad,
  comparator = "continent",
  facet_variable = "study_group_name"
)
# Run the function with default p-value adjustment (FDR)
con_sgn_adonis_results <- adonis_with_subsets(
  physeq_object = sgn_rel_ad,
  comparator = "continent"
)

con_boxplot_shannon
con_boxplot_obs
con_boxplot_faceted
con_beta_all_plot
con_beta_sgn_plot
con_sgn_adonis_results

# Combine the top row
A <- plot_grid(
  con_beta_all_plot, 
  labels = c('A'), 
  label_size = 18
)

# Combine the top row
B <- plot_grid(
  con_beta_sgn_plot, 
  labels = c('B'), 
  label_size = 18
)

ggsave(filename = file.path(outdir, "adults_beta_all_con.png"), plot = A, width = 8.5, height = 8, dpi = 300 )
ggsave(filename = file.path(outdir, "adults_beta_SGN_con.png"), plot = B, width = 14, height = 7.5, dpi = 300 )


# Adults ALL + SGN --------------------------------------------------------
# Agglomerate at the genus level 
rel_ph <- tax_glom(rel_ad, taxrank = 'Phylum')
rel_ph_df <- psmelt(rel_ph)

# Reorder Phylum based on total Abundance
rel_ph_df <- rel_ph_df %>%
  mutate(Phylum = fct_reorder(Phylum, Abundance, .fun = sum, .desc = TRUE))

# Select colors from Set1 and combine with Set3
set1_colors <- brewer.pal(8, "Set1")  # Maximum 9 colors
set3_colors <- brewer.pal(12, "Set3")[1:6]  # Select the first 3 colors from Set3

# Combine colors to create a custom palette
custom_colors <- c(set1_colors, set3_colors)

## Plot by AD status (abundances summarized by group)
rel_ph_IBD_plot_ad <- ggplot(rel_ph_df, aes(x=IBD_group_name, y=Abundance, fill=Phylum))+
  geom_bar(stat='identity', position='fill')+
  theme_classic()+
  theme(axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=10),
        plot.title = element_text(size=15),
        panel.border = element_rect(fill=NA, colour="black", size=3)) +
  scale_fill_manual(values = custom_colors)

# Construct the filename dynamically
output_filename <- sprintf("Figures/one/barplot_phy_SGN_adult.png")

# Save the plot
ggsave(
  filename = output_filename,
  plot = rel_ph_IBD_plot_ad,
  width = 6.5,
  height = 7.5,
  dpi = 300
)

# alpha
alpha <- microbiome::alpha(rel_ad, index = "all")
merge <- merge(metadata, alpha, by = "row.names")
alpha.2 <- microbiome::alpha(sgn_rel_ad, index = "all")
mergey <- merge(metadata, alpha.2, by = "row.names")

# alpha
# Generate boxplot
AD_all_IBD_boxplot_shannon <- generate_boxplot(data = merge, variable = "diversity_shannon", population_col = "IBD_group_name")
AD_all_IBD_boxplot_obs <- generate_boxplot(data = merge, variable = "observed", population_col = "IBD_group_name")

# Generate the boxplot with ANOVA and Tukey's HSD test
AD_all_SGN_boxplot_shannon <- generate_boxplot_anova(data = mergey, variable = "diversity_shannon", grouping_col = "study_group_name")
AD_all_SGN_boxplot_obs <- generate_boxplot_anova(data = mergey, variable = "observed", grouping_col = "study_group_name")

# Beta
AD_all_beta_IBD_plot <- plot_beta_diversity(rel_ad, "IBD_group_name")
AD_all_beta_sgn_plot <- plot_beta_diversity(sgn_rel_ad, "study_group_name")


# Fig 3
rel_ph_IBD_plot_ad
AD_all_IBD_boxplot_shannon <- AD_all_IBD_boxplot_shannon$plot
AD_all_IBD_boxplot_obs <- AD_all_IBD_boxplot_obs$plot
rel_ph_SGN_plot_ad
AD_all_SGN_boxplot_shannon <- AD_all_SGN_boxplot_shannon$plot
AD_all_SGN_boxplot_obs <- AD_all_SGN_boxplot_obs$plot
AD_all_beta_IBD_plot
AD_all_beta_sgn_plot

# Combine the top row
fig_2_IBD.u <- plot_grid(
  rel_ph_IBD_plot_ad,
  AD_all_beta_IBD_plot,
  labels = c('A', 'B'), 
  ncol = 2,
  label_size = 18,
  axis = "tblr",
  align = "h",
  rel_widths = c(1,1.5)
)

# Combine the top row
fig_2_IBD.l <- plot_grid(
  AD_all_IBD_boxplot_shannon,
  AD_all_IBD_boxplot_obs,
  labels = c('C', 'D'), 
  ncol = 2,
  label_size = 18,
  axis = "tblr",
  align = "h"
)

fig2_ibd <- plot_grid(fig_2_IBD.u, fig_2_IBD.l, ncol = 1)

fig_2_sgn.u <- plot_grid(
  rel_ph_SGN_plot_ad,
  AD_all_beta_sgn_plot,
  labels = c('A', 'B'), 
  ncol = 2,
  label_size = 18,
  axis = "tblr",
  align = "h",
  rel_widths = c(1,1.5)
)

fig_2_sgn.l <- plot_grid(
  AD_all_SGN_boxplot_shannon,
  AD_all_SGN_boxplot_obs,
  labels = c('C', 'D'), 
  ncol = 2,
  label_size = 18,
  axis = "tblr",
  align = "h"
)

fig2_sgn <- plot_grid(fig_2_sgn.u, fig_2_sgn.l, ncol = 1)

ggsave(filename = file.path(outdir, "fig3_IBD_all_ad.png"), plot = fig2_ibd, width = 12.5, height = 15.5, dpi = 300)

ggsave(filename = file.path(outdir, "fig3_sgn_all_ad.png"), plot = fig2_sgn, width = 15, height = 20, dpi = 300)

# Maaslin2 ----------------------------------------------------------------
# Adults only IBD

# Adults only SGN
# Subtypes
# Extract metadata
cd <- subset_samples(sgn_rel_ad, study_group_name == "CD" | study_group_name == "healthy_control")
uc <- subset_samples(sgn_rel_ad, study_group_name == "UC" | study_group_name == "healthy_control")

metadata_df <- as.data.frame(sample_data(cd))
metadata_df <- as.matrix(metadata_df)
metadata_df <- as.data.frame(metadata_df)

# Extract OTU table
otu <- otu_table(cd)
feat <- as.data.frame(otu)
feat <- t(feat)

## Tax table
tax <- tax_table(cd)
tax <- as.data.frame(tax)

# Ensure row names of df1 match column names of df2
matching_indices <- match(colnames(feat), rownames(tax))

# Replace column names in df2 with values from the "Genus" column in df1
colnames(feat)[!is.na(matching_indices)] <- tax$Genus[matching_indices[!is.na(matching_indices)]]

fit_data_random = Maaslin2(input_data     = feat, 
                           input_metadata = metadata_df, 
                           min_prevalence = 0.1,
                           max_significance = 0.05,
                           normalization  = "NONE",
                           output         = "IBD_rel_ad",
                           fixed_effects  = c("IBD_group_name","continent", "V_region"),
                           reference = "IBD_group_name,healthy_control")

fit_data_random = Maaslin2(input_data     = feat, 
                           input_metadata = metadata_df, 
                           min_prevalence = 0.1,
                           max_significance = 0.05,
                           normalization  = "NONE",
                           output         = "new_cd_rel_ad_v", 
                           fixed_effects  = c("study_group_name","continent", "V_region"),
                           reference = "study_group_name,healthy_control")

#maaslin output visualisation
df <- read.table("Chapter_3/Analysis/IBD_rel_ad/significant_results.tsv", 
                 sep = "\t", header = TRUE)
df <- read.table("Chapter_3/Analysis/cd_rel_ad_v/significant_results.tsv", 
                 sep = "\t", header = TRUE)
df <- read.table("Chapter_3/Analysis/uc_rel_ad_v/significant_results.tsv", 
                 sep = "\t", header = TRUE)

# Filter significant feature associated with IBD
df <- df %>%
  dplyr::filter(value == "CD")

# Keep features q < 0.05
df<- df[df$qval <= 0.05, ]
df$log_qval <- -log10(df$qval)  # Compute -log10(q_val) for size scaling

# Add a significance threshold line (e.g., q < 0.05)
significance_threshold <- -log10(0.05)

# Create a new column to classify positive and negative coefficients
df$effect_direction <- ifelse(df$coef > 0, "Positive", "Negative")

# Define colors
positive_color <- "#00AFBB"
negative_color <- "#FC4E07"

# Volcano Plot with Thicker Horizontal Arrows
volcano_plot <- ggplot(df, aes(x = coef, y = log_qval, label = feature)) +
  geom_point(aes(fill = effect_direction), 
             color = "black",       # Outline color
             shape = 21,            # Filled circles with border
             size = 5,              # Point size
             stroke = 1.2) +        # Outline thickness
  geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "red") +
  geom_text_repel(data = subset(df, qval < 0.05), 
                  size = 5, 
                  max.overlaps = 10, 
                  box.padding = 0.4) +  # Labels for significant points
  scale_fill_manual(values = c("Positive" = positive_color, "Negative" = negative_color)) +  # Fill colors
  labs(x = "Model Coefficient", 
       y = "-log10(q-value)", 
       fill = "") +
  theme_minimal(base_size = 20) +
  theme(legend.position = "top",
        panel.border = element_rect(color = "black", fill = NA, size = 1))  # Plot border

# Add Annotations (IBD and HC) with Outward-Pointing Arrows and Text at the Arrow Tail
volcano_plot <- volcano_plot +
  
  # IBD (top right, positive color) - Arrow pointing OUTWARD (right)
  annotate("text", x = max(df$coef) - 0.1, y = max(df$log_qval) + 0.2, 
           label = "CD", size = 10, fontface = "bold", color = positive_color, hjust = 1) +
  geom_segment(aes(x = max(df$coef), y = max(df$log_qval), 
                   xend = max(df$coef) + 0.4, yend = max(df$log_qval)), 
               arrow = arrow(length = unit(0.3, "inches"), type = "closed"), 
               color = positive_color, size = 3, lineend = "round") +  # Thicker arrow tail
  
  # HC (top left, negative color) - Arrow pointing OUTWARD (left)
  annotate("text", x = min(df$coef) + 0.1, y = max(df$log_qval) + 0.2, 
           label = "HC", size = 10, fontface = "bold", color = negative_color, hjust = 0) +
  geom_segment(aes(x = min(df$coef), y = max(df$log_qval), 
                   xend = min(df$coef) - 0.4, yend = max(df$log_qval)), 
               arrow = arrow(length = unit(0.3, "inches"), type = "closed"), 
               color = negative_color, size = 3, lineend = "round")  # Thicker arrow tail

# OFFSETTING #######
# Adjust the vertical offset
offset_y <- 0.5  # Adjust this value to control how much you move it upwards

# Add Annotations (IBD and HC) with Adjustments
volcano_plot <- volcano_plot +
  # CD (top right, positive color) - Arrow pointing OUTWARD (right)
  annotate("text", x = max(df$coef) - 0.1, y = max(df$log_qval) + 0.2 + offset_y, 
           label = "CD", size = 10, fontface = "bold", color = positive_color, hjust = 1) +
  geom_segment(aes(x = max(df$coef), y = max(df$log_qval) + offset_y, 
                   xend = max(df$coef) + 0.4, yend = max(df$log_qval) + offset_y), 
               arrow = arrow(length = unit(0.3, "inches"), type = "closed"), 
               color = positive_color, size = 3, lineend = "round") +  # Thicker arrow tail
  
  # HC (top left, negative color) - Arrow pointing OUTWARD (left)
  annotate("text", x = min(df$coef) + 0.1, y = max(df$log_qval) + 0.2 + offset_y, 
           label = "HC", size = 10, fontface = "bold", color = negative_color, hjust = 0) +
  geom_segment(aes(x = min(df$coef), y = max(df$log_qval) + offset_y, 
                   xend = min(df$coef) - 0.4, yend = max(df$log_qval) + offset_y), 
               arrow = arrow(length = unit(0.3, "inches"), type = "closed"), 
               color = negative_color, size = 3, lineend = "round")  # Thicker arrow tail


ggsave(filename = file.path(outdir, "fig4_maaslin2_IBD_ad_vol.png"), plot = volcano_plot, width = 11, height = 10, dpi = 300)
ggsave(filename = file.path(outdir, "fig4_maaslin2_cd_ad_vol.png"), plot = volcano_plot, width = 11, height = 10.5, dpi = 300)
