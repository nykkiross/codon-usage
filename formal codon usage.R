# Total GC, GC1, GC2, and GC3 calculations
# followed by RSCU calculations and plotting functions

# calculating GC content -------------------
# Load necessary libraries
install.packages("seqinr")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")

library(seqinr)
library(dplyr)
library(ggplot2)
library(tidyr)

# load your files

organism1 <- read.fasta("your/file/path.fasta")
organism2 <- read.fasta("your/file/path2.fasta")
# etc; can add as many organisms as you want to compare; can also compare groups of genes instead
gene_group1 <- read.fasta("your/file/path.fasta")
gene_group2 <- read.fasta("your/file/path2.fasta")

# Function to calculate GC, GC1, GC2, and GC3 per gene
calculate_gc_per_gene <- function(fasta_data) {
  gc_values <- sapply(fasta_data, function(seq) {
    seq_vec <- unlist(seq)  # Convert list element to character vector
    c(GC(seq_vec), GC1(seq_vec), GC2(seq_vec), GC3(seq_vec))
  })
  return(t(gc_values))  # Transpose to ensure columns represent GC types
}

gc_data <- list(
  "Organism 1" = calculate_gc_per_gene(organism1),
  "Organism 2" = calculate_gc_per_gene(organism2),
  "Gene group 1" = calculate_gc_per_gene(gene_group1),
  "Gene group 2" = calculate_gc_per_gene(gene_group2)
)

gc_results <- do.call(rbind, lapply(names(gc_data), function(group) {
  df <- as.data.frame(gc_data[[group]])  # Convert matrix to data frame
  colnames(df) <- c("GC", "GC1", "GC2", "GC3")  # Name columns
  df$Organism <- group  # Add group label
  return(df)
}))

# Compute and display average GC, GC1, GC2, and GC3 for each group
gc_averages <- gc_results %>%
  group_by(Organism) %>%
  summarise(across(GC:GC3, mean, na.rm = TRUE))
print(gc_averages)
# write csv
write.csv(gc_averages, "your/file/path/GC_content_averages.csv", row.names = FALSE)

# can set a custom order for plotting
custom_order <- c(
  "Organism 1", "Gene Group 1", "Organism 2", "Gene Group 2"
)

# Convert 'Organism' column to a factor with the specified order
gc_long$Organism <- factor(gc_long$Organism, levels = custom_order)

# Plot the average GC content for each group with the specified order
ggplot(gc_long, aes(x = Organism, y = Mean_GC, fill = GC_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Average GC Content per Organism/Gene Group", x = "Organism/Gene Group", y = "Mean GC Content (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# calculation and plotting of RSCU values -----------------------
# Install and load necessary packages
install.packages(c("seqinr", "ggplot2", "reshape2", "pheatmap", "factoextra"))

library(seqinr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(factoextra)

# Function to calculate codon usage
calculate_codon_usage <- function(sequence, seq_name) {
  sequence <- unlist(sequence)
  codon_usage <- uco(sequence, index = "rscu")
  codon_usage_df <- as.data.frame(t(codon_usage))
  codon_usage_df$Sequence <- seq_name
  return(codon_usage_df)
}

# Function to process a FASTA file and calculate codon usage for each gene
process_fasta_file <- function(file, group_name) {
  sequences <- read.fasta(file = file, seqtype = "DNA")
  codon_usage_list <- lapply(sequences, function(seq) {
    seq_name <- attr(seq, "name")
    calculate_codon_usage(seq, seq_name)
  })
  combined_df <- do.call(rbind, codon_usage_list)
  combined_df$Group <- group_name
  return(combined_df)
}

# OPTIONAL - set codon order
codon_order <- c("aaa", "aac", "aag", "aca", "aga", "agc", "atc", "atg", "caa",
                 "cac", "cag", "cca", "cgt", "cta", "gaa", "gac", "gca", "gga",
                 "gta", "tac", "tca", "tgc", "tta", "ttc", "ttg", "aat", "acc",
                 "acg", "act", "agg", "agt", "ata", "att", "cat", "ccc", "ccg",
                 "cct", "cga", "cgc", "cgg", "ctc", "ctg", "ctt", "gag", "gat",
                 "gcc", "gcg", "gct", "ggc", "ggg", "ggt", "gtc", "gtg", "gtt",
                 "taa", "tag", "tat", "tcc", "tcg", "tct", "tga", "tgg", "tgt",
                 "ttt")

# define file paths and group names
organism1_files <- c("your/file/path1.fasta", "your/file/path2.fasta", "your/file/path3.fasta", ..., "your/file/pathN.fasta")
organism2_files <- c("your/file/path1.fasta", "your/file/path2.fasta", "your/file/path3.fasta", ..., "your/file/pathN.fasta")
organism1_groups <- c("Fasta 1", "Fasta 2", "Fasta 3", ..., "Fasta N")
organism2_groups <- c("Fasta 1", "Fasta 2", "Fasta 3", ..., "Fasta N")

# Process each FASTA file and save as .csv files
organism1_codon_usage <- do.call(rbind, mapply(process_fasta_file, organism1_files, organism1_groups, SIMPLIFY = FALSE))
write.csv(organism1_codon_usage, file = "/file_path/organism1_codon_usage_data.csv")
organism2_codon_usage <- do.call(rbind, mapply(process_fasta_file, organism2_files, organism2_groups, SIMPLIFY = FALSE))
write.csv(organism2_codon_usage, file = "/file_path/organism2_codon_usage_data.csv")

# Combine both groups into a single data frame
combined_codon_usage <- rbind(organism1_codon_usage, organism2_codon_usage)

# Replace NA values with zeros
combined_codon_usage[is.na(combined_codon_usage)] <- 0

# additional checks to make sure your data is ready for plotting and analysis
# Ensure all values are numeric
codon_usage_matrix <- as.matrix(combined_codon_usage[, -c(ncol(combined_codon_usage)-1, ncol(combined_codon_usage))])
codon_usage_matrix <- apply(codon_usage_matrix, 2, function(x) suppressWarnings(as.numeric(x)))

# Replace any remaining NA values with zeros
codon_usage_matrix[is.na(codon_usage_matrix)] <- 0

# Ensure no infinite values are present
codon_usage_matrix[is.infinite(codon_usage_matrix)] <- 0

# Set row names for the matrix
row_names <- paste(combined_codon_usage$Group, combined_codon_usage$Sequence, sep = "-")
row.names(codon_usage_matrix) <- row_names

# Check for any remaining NA/NaN/Inf values
if (any(is.na(codon_usage_matrix)) || any(is.nan(codon_usage_matrix)) || any(is.infinite(codon_usage_matrix))) {
  print("Warning: NA, NaN, or Inf values detected in the codon usage matrix.")
}

# plot using PCA -----------
# Remove columns with zero variance
codon_usage_matrix <- codon_usage_matrix[, apply(codon_usage_matrix, 2, var) != 0]

# Perform PCA
codon_usage_pca <- prcomp(codon_usage_matrix, scale. = TRUE)

# Prepare data for PCA plot
pca_data <- as.data.frame(codon_usage_pca$x)
pca_data$Group <- combined_codon_usage$Group
pca_data$Sequence <- combined_codon_usage$Sequence

# Set the order of the factor levels for Group
pca_data$Group <- factor(pca_data$Group, levels = c("Organism 1", "Gene Group 1", "Organism 2", "Gene Group 2"))

# Define custom colors for the groups
group_colors <- c("Organism 1" = "red",
                  "Gene Group 1" = "purple",
                  "Organism 2" = "green",
                  "Gene Group 2" = "blue")


# Plot PCA without labels
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "Relative Synonymous Codon Usage in Organism 1 and Organism 2", x = "Principal Component 1", y = "Principal Component 2") +
  scale_color_manual(values = group_colors) +  # Apply custom colors
  theme_minimal()

# alternatively, plot PCA WITH labels
install.packages("ggrepel")
library(ggrepel)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sequence)) + 
  geom_point(size = 3) + 
  geom_text_repel(size = 3, max.overlaps = 15) +  # Adjust max.overlaps as needed
  labs(title = "Relative Synonymous Codon Usage in S. flexneri 2457T and phage Sf14",
       x = "Principal Component 1",
       y = "Principal Component 2") + 
  scale_color_manual(values = group_colors) +  
  theme_minimal()

# plot using heatmap -----------------
# this heatmap assumes multiple levels of groups, i.e. 2 organisms with multiple gene groups each
# if multiple levels are not needed, skip to line 223 and replace "annotation_df_ordered" with "annotation_df"
# also replace "codon_usage_matrix_ordered" with "codon_usage_matrix"

# Custom clustering function
custom_clustering <- function(data, annotation) {
  # Split the data by groups
  data_split <- split(as.data.frame(data), annotation$Group)
  # Cluster each group separately
  clustered_data <- lapply(data_split, function(group_data) {
    # Ensure group_data is a matrix
    group_data_matrix <- as.matrix(group_data)
    # Perform hierarchical clustering
    dist_matrix <- dist(group_data_matrix)
    hc <- hclust(dist_matrix)
    return(group_data_matrix[hc$order, , drop = FALSE])
  })
  # Combine the clustered groups
  clustered_data_combined <- do.call(rbind, clustered_data)
  # Return both the clustered data and the new row order
  new_order <- unlist(lapply(clustered_data, rownames))
  return(list(data = clustered_data_combined, order = new_order))
}

# Use the custom clustering function to order the rows
clustered_result <- custom_clustering(codon_usage_matrix, annotation_df)
codon_usage_matrix_ordered <- clustered_result$data

# Update the annotation data frame to match the new order
annotation_df_ordered <- annotation_df[clustered_result$order, , drop = FALSE]

# Ensure that the row names in the ordered matrix match the annotation data frame
row.names(annotation_df_ordered) <- row.names(codon_usage_matrix_ordered)

# Ensure all codons in the dataset are in the defined order
existing_codons <- colnames(codon_usage_matrix)
codon_order_filtered <- codon_order[codon_order %in% existing_codons]

# Reorder the matrix columns to match codon_order
codon_usage_matrix <- codon_usage_matrix[, codon_order_filtered]

# Generate the heatmap with custom clustering
pheatmap(codon_usage_matrix,
         cluster_rows = FALSE,   # Disable clustering of rows
         cluster_cols = FALSE,   # Disable clustering of columns
         show_rownames = FALSE,  # Option to show/hide row names
         annotation_row = annotation_df,
         annotation_colors = list(Group = c("Organism 1" = "red",
                                            "Gene Group 1" = "purple",
                                            "Organism 2" = "green",
                                            "Gene Group 2" = "blue")))

# STATISTICAL ANALYSES -----------------
# analyses of GC content ------------------

# Kruskal-Wallis tests for each GC type
gc_types <- c("GC", "GC1", "GC2", "GC3")
kruskal_results <- list()
dunn_results_list <- list()

for (gc_type in gc_types) {
  cat("\nPerforming Kruskal-Wallis test for", gc_type, "...\n")
  kruskal_test <- kruskal.test(GC_Content ~ Organism, data = gc_long %>% filter(GC_Type == gc_type))
  print(kruskal_test)
  kruskal_results[[gc_type]] <- kruskal_test
  # If Kruskal-Wallis test is significant, perform Dunn’s post-hoc test
  if (kruskal_test$p.value < 0.05) {
    dunn_test <- dunnTest(GC_Content ~ Organism, 
                          data = gc_long %>% filter(GC_Type == gc_type), 
                          method = "bonferroni")
    # Convert results to a data frame
    dunn_results_df <- as.data.frame(dunn_test$res)
    dunn_results_list[[gc_type]] <- dunn_results_df
    # Define file path for saving
    file_path <- paste0("/your/file/path/GC_", gc_type, "_DunnTest_Results.tsv")
    # Save as tab-delimited file
    write.table(dunn_results_df, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Dunn’s test results for", gc_type, "saved to:", file_path, "\n")
  }
}

# Boxplot with significance annotation for each GC Type
ggplot(gc_long, aes(x = Organism, y = GC_Content, fill = Organism)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test", label = "p.format") +  # Kruskal-Wallis p-value
  facet_wrap(~GC_Type, scales = "free") +  # Separate plots for GC, GC1, GC2, GC3
  labs(title = "GC Content Comparison by Gene Group", y = "GC Content", x = "Gene Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# RSCU stats ------------------------------------------------------

library(tidyr)
library(dplyr)

# Replace with your actual RSCU objects
rscu_org1 <- organism1_codon_usage
rscu_org2 <- organism2_codon_usage

# Add a column identifying the organism
rscu_org1$Organism <- "Organism 1"
rscu_org2$Organism <- "Organism 2"

# Combine both data frames
rscu_combined <- bind_rows(rscu_org1, rscu_org2)

numeric_cols <- sapply(rscu_combined, is.numeric)  # Identify numeric columns
codon_cols <- names(rscu_combined)[numeric_cols]  # Get codon column names

# Reshape to long format, keeping only numeric codon columns
rscu_long <- pivot_longer(
  rscu_combined,
  cols = all_of(codon_cols),  # Use only numeric codon columns
  names_to = "Codon",        # Column for codon names
  values_to = "RSCU"         # Column for RSCU values
)

# perform Wilcoxon rank-sum test for each codon
results <- rscu_long %>%
  group_by(Codon) %>%
  summarize(
    p_value = wilcox.test(RSCU[Organism == "Organism 1"], RSCU[Organism == "Organism 2"])$p.value,
    median_org1 = median(RSCU[Organism == "Organism 1"], na.rm = TRUE),
    median_org2 = median(RSCU[Organism == "Organism 2"], na.rm = TRUE),
    .groups = "drop"
  )
results <- results %>%
  mutate(Significant = ifelse(p_value < 0.05, "Yes", "No"))

# View the results
print(results)

# Save the results to a CSV file
write.csv(results, "/your/file/path/rscu_comparison_results.csv", row.names = FALSE)

# analysis of PCA data -----------------------
# PERMANOVA (Permutational Multivariate Analysis of Variance)
install.packages("vegan")
library(vegan)

distance_matrix <- dist(codon_usage_matrix, method = "euclidean")

permanova_result <- adonis2(distance_matrix ~ Group, data = combined_codon_usage, permutations = 999)
print(permanova_result)

# Pairwise PERMANOVA
codon_data <- combined_codon_usage[, !(colnames(combined_codon_usage) %in% c("Sequence", "Group"))]

anyNA(codon_data)
anyNA(combined_codon_usage$Group)

table(combined_codon_usage$Group)
nrow(codon_data) == sum(table(combined_codon_usage$Group))

distance_matrix <- vegdist(codon_data, method = "euclidean")

sub_data <- codon_data[combined_codon_usage$Group %in% c("Group 1", "Group 2"), ]
sub_factors <- combined_codon_usage$Group[combined_codon_usage$Group %in% c("Group 1", "Group 2")]

cat("Subset dimensions:", dim(sub_data), "\n")
cat("Subset factor levels:", unique(sub_factors), "\n")
dist_matrix <- vegdist(sub_data, method = "euclidean")
cat("Distance matrix dimensions:", dim(as.matrix(dist_matrix)), "\n")

adonis_result <- try(adonis2(dist_matrix ~ sub_factors))
adonis_result






