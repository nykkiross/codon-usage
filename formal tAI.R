# tAI Calculation using RSCU - Genome-Specific tRNA Adaptation Index Analysis ------------

# PART 1: Obtaining and uploading data from tRNAscan-SE for your organism(s)
# Load your organism into tRNAscan-SE and download the results into excel
# If you have your tRNAscan data saved in another format, extract your data accordingly
install.packages("readxl")
install.packages("dplyr")

library(readxl)
library(dplyr)

# function to process multiple excel files - assuming multiple organisms
process_multiple_excel_files <- function(tRNAscan_paths) {
  # Loop through each tRNAscan file
  for (file_path in tRNAscan_paths) {
    # Get the base file name without extension (for unique naming)
    file_name <- tools::file_path_sans_ext(basename(file_path))
    # Get all sheet names in the Excel file
    sheet_names <- excel_sheets(file_path)
    # Loop through each sheet in the file for tRNAscan-SE data
    for (sheet in sheet_names) {
      # Read the sheet data into a data frame
      sheet_data <- read_excel(file_path, sheet = sheet)
      # Clean and structure the tRNAscan-SE data
      tRNA_data <- sheet_data %>%
        select(Name, `tRNA_Type`, `Anti_Codon`, Score)
      # Create a unique name for the tRNAscan-SE data frame
      df_name <- paste0(file_name, "_tRNAscan_", sheet)
      # Assign the tRNAscan-SE data frame to the global environment
      assign(df_name, tRNA_data, envir = .GlobalEnv)
      # Print a message for each data frame created
      print(paste("Created tRNAscan data frame:", df_name))
    }
  }
}

tRNAscan_paths <- c("/file_path_1.xlxs", "/file_path_2.xlxs", ..., "/file_path_n.xlxs")
process_multiple_excel_files(tRNAscan_paths)

# rename your objects to more concise names
organism_tRNA <- file_tRNAscan_organism

# PART 2: Calculate Relative Synonymous Codon Usage (RSCU) for each organism ---------------
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


# Define file paths and group names
organism_genes <- ("/file_path.fasta")

# Calculate and save RSCU values for each organism
organism_rscu <- process_fasta_file(organism_genes, "Organism")
write.csv(organism_rscu, file = "/file_path/organism_rscu_data.csv")

# Reformat RSCU for proper alignment of codons and anticodons during tAI calculations ------
library(tidyverse)

# Codon Case Change
case_change <- c(
  "aaa" = "AAA", "aac" = "AAC", "aag" = "AAG", "aat" = "AAT",
  "aca" = "ACA", "acc" = "ACC", "acg" = "ACG", "act" = "ACT",
  "aga" = "AGA", "agc" = "AGC", "agg" = "AGG", "agt" = "AGT",
  "ata" = "ATA", "atc" = "ATC", "atg" = "ATG", "att" = "ATT",
  "caa" = "CAA", "cac" = "CAC", "cag" = "CAG", "cat" = "CAT",
  "cca" = "CCA", "ccc" = "CCC", "ccg" = "CCG", "cct" = "CCT",
  "cga" = "CGA", "cgc" = "CGC", "cgg" = "CGG", "cgt" = "CGT",
  "cta" = "CTA", "ctc" = "CTC", "ctg" = "CTG", "ctt" = "CTT",
  "gaa" = "GAA", "gac" = "GAC", "gag" = "GAG", "gat" = "GAT",
  "gca" = "GCA", "gcc" = "GCC", "gcg" = "GCG", "gct" = "GCT",
  "gga" = "GGA", "ggc" = "GGC", "ggg" = "GGG", "ggt" = "GGT",
  "gta" = "GTA", "gtc" = "GTC", "gtg" = "GTG", "gtt" = "GTT",
  "taa" = "TAA", "tac" = "TAC", "tag" = "TAG", "tat" = "TAT",
  "tca" = "TCA", "tcc" = "TCC", "tcg" = "TCG", "tct" = "TCT",
  "tga" = "TGA", "tgc" = "TGC", "tgg" = "TGG", "tgt" = "TGT",
  "tta" = "TTA", "ttc" = "TTC", "ttg" = "TTG", "ttt" = "TTT"
)

organism_rscu_long <- organism_rscu %>%
  pivot_longer(
    cols = -c(Sequence, Group),
    names_to = "Codon",
    values_to = "RSCU"
  ) %>%
  # Add AmAcid column by matching each codon with its amino acid in lowercase
  mutate(AmAcid = case_when(
    Codon %in% c("gca", "gcc", "gct", "gcg") ~ "Ala",
    Codon %in% c("aga", "agg", "cga", "cgc", "cgg", "cgt") ~ "Arg",
    Codon %in% c("aac", "aat") ~ "Asn",
    Codon %in% c("gac", "gat") ~ "Asp",
    Codon %in% c("tgc", "tgt") ~ "Cys",
    Codon %in% c("gaa", "gag") ~ "Glu",
    Codon %in% c("ttt", "ttc") ~ "Phe",
    Codon %in% c("gga", "ggc", "ggg", "ggt") ~ "Gly",
    Codon %in% c("cac", "cat") ~ "His",
    Codon %in% c("ata", "atc", "att") ~ "Ile",
    Codon %in% c("aaa", "aag") ~ "Lys",
    Codon %in% c("cta", "ctc", "ctg", "ctt", "tta", "ttg") ~ "Leu",
    Codon == "atg" ~ "Met",
    Codon %in% c("cca", "ccc", "ccg", "cct") ~ "Pro",
    Codon %in% c("caa", "cag") ~ "Gln",
    Codon %in% c("agc", "agt", "tca", "tcc", "tcg", "tct") ~ "Ser",
    Codon %in% c("aca", "acc", "acg", "act") ~ "Thr",
    Codon %in% c("gta", "gtc", "gtg", "gtt") ~ "Val",
    Codon == "tgg" ~ "Trp",
    Codon %in% c("tac", "tat") ~ "Tyr",
    TRUE ~ NA_character_
  )) %>%
  # Map Codon to uppercase Anti-Codon values
  mutate(Codon = case_change[Codon])
# Remove redundant columns and write to CSV
organism_rscu_long <- organism_rscu_long %>%
  select(Sequence, Group, Codon, RSCU, AmAcid)
# View the modified data
head(organism_rscu_long)
write.csv(organism_rscu_long, file = "/file_path/organism_rscu.csv")


# PART 3: Calculate tRNA Adaptation Index using RSCU and tRNAscan-SE data ------------
library(seqinr)
library(dplyr)

# Step 1: set data and dictionary
codon_dictionary <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT",
                      "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", 
                      "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", 
                      "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", 
                      "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT",
                      "TTA", "TTC", "TTG", "TTT")

# Ensure anticodon_dictionary is a vector of anticodons
anticodon_dictionary <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT",
                          "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", 
                          "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", 
                          "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", 
                          "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT",
                          "TTA", "TTC", "TTG", "TTT")

# organism data
codon_data <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT",
                "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", 
                "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", 
                "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", 
                "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT",
                "TTA", "TTC", "TTG", "TTT")

anticodon_data <- c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT",
                    "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", 
                    "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", 
                    "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", 
                    "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT",
                    "TTA", "TTC", "TTG", "TTT")

# Step 2: Check for codons and anticodons not in dictionary
missing_codons <- codon_data[!codon_data %in% codon_dictionary]
missing_anticodons <- anticodon_data[!anticodon_data %in% anticodon_dictionary]

# Display results
if (length(missing_codons) > 0) {
  cat("Codons not found in dictionary:\n")
  print(missing_codons)
} else {
  cat("All codons found in dictionary.\n")
}

if (length(missing_anticodons) > 0) {
  cat("Anti-Codons not found in dictionary:\n")
  print(missing_anticodons)
} else {
  cat("All anti-codons found in dictionary.\n")
}

for (anticodon in anticodon_data) {
  if (!(anticodon %in% anticodon_dictionary)) {
    print(paste("Missing anticodon:", anticodon))
  }
}

# Step 3: Generate organism-specific anticodon dictionary from general dictionary
# This step is necessary to turn the anticodon_dictionary object from an atomic
# vector into a dataframe that includes anticodon scores from tRNAscan-SE data

# Step 3a: function to accept raw tRNAscan-SE data and return anticodon data frame
generate_anticodon_dictionary <- function(trna_data) {
  # Ensure the input is a data frame with the expected columns
  if (!all(c("Anti_Codon", "Score") %in% colnames(trna_data))) {
    stop("Input data must contain 'Anti_Codon' and 'Score' columns.")
  }
  
  # Remove any rows with NA values in 'Anti_Codon' or 'Score'
  trna_data <- trna_data[!is.na(trna_data$Anti_Codon) & !is.na(trna_data$Score), ]
  
  # Create a named vector of scores using the Anti_Codon column as names
  anticodon_dictionary <- setNames(trna_data$Score, trna_data$Anti_Codon)
  
  return(anticodon_dictionary)
}

organism_tRNA_dict <- generate_anticodon_dictionary(organism_tRNA)
# if performing analysis with multiple organisms, generate a list of all anticodon libraries
all_tRNA_dictionaries <- list(org1, org2, org3, ..., orgN)


# Step 3b: reverse complement for matching to codons
reverse_complement <- function(codon) {
  # Check for NA input and return NA if codon is missing
  if (is.na(codon)) {
    return(NA)
  }
  
  # Define complementary base pairs
  complement_map <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G")
  
  # Convert codon to uppercase
  codon <- toupper(codon)
  
  # Split codon into individual bases, find complements, and reverse order
  rev_comp <- paste0(rev(complement_map[strsplit(codon, NULL)[[1]]]), collapse = "")
  
  return(rev_comp)
}

# OPTIONAL DEBUGGING ----------------------------------------------------------
# Step 1: Check unique codons in RSCU data
unique_codons <- unique(organism_rscu_long$Codon)
print("Unique codons in RSCU data:")
print(unique_codons)

# Step 2: Check unique anti-codons in tRNA dictionary
unique_anticodons <- names(organism_tRNA_dict)
print("Unique anticodons in tRNA dictionary:")
print(unique_anticodons)

# Step 3: Identify missing anticodon matches
missing_anticodons <- unique_codons[!sapply(unique_codons, function(codon) {
  anti_codon <- reverse_complement(codon)
  anti_codon %in% unique_anticodons
})]

if (length(missing_anticodons) > 0) {
  print("Codons without matching anticodons in the tRNA dictionary:")
  print(missing_anticodons)
} else {
  print("All codons have matching anticodons in the tRNA dictionary.")
}

unmatched <- setdiff(sapply(unique_codons, reverse_complement), unique_anticodons)
print("Unmatched codons after reverse complement:")
print(unmatched)

# END DEBUGGING
# RESUME tAI CALCULATION -----------------------------------------------------

# Step 4: calculation of tAI
calculate_tai_with_rscu <- function(rscu_data, anticodon_dictionary, full_genome_tai = NULL, min_weight = 0.01) {
  required_columns <- c("Codon", "RSCU", "Sequence")
  if (!all(required_columns %in% colnames(rscu_data))) {
    stop(paste("RSCU data must contain the following columns:", paste(required_columns, collapse = ", ")))
  }
  
  # Compute tAI weights for each codon
  rscu_data$tAI_weight <- mapply(function(codon, rscu_value) {
    anticodon <- reverse_complement(codon)
    
    # Find all matching anticodon scores
    matching_scores <- anticodon_dictionary[names(anticodon_dictionary) == anticodon]
    
    if (length(matching_scores) > 0) {
      # Sum all matching scores (allow for multiple tRNA genes for the same anticodon)
      norm_tRNA_score <- sum(matching_scores)
      tAI_weight <- (rscu_value + 1e-6) * pmax(norm_tRNA_score, min_weight)
    } else {
      # Assign fallback weight if anticodon is missing
      tAI_weight <- min_weight
    }
    
    return(tAI_weight)
  }, rscu_data$Codon, rscu_data$RSCU, SIMPLIFY = TRUE)
  
  # Ensure no zero values (avoid log(0))
  rscu_data$tAI_weight[rscu_data$tAI_weight == 0] <- min_weight
  
  # Calculate tAI for each gene
  gene_tai <- rscu_data %>%
    group_by(Sequence) %>%
    summarize(
      num_matched_codons = sum(tAI_weight > min_weight, na.rm = TRUE),  # Count matched codons
      total_codons = n(),  # Count total codons
      raw_tAI = exp(sum(log(pmax(tAI_weight, min_weight)), na.rm = TRUE) / num_matched_codons),  # Adjust log scaling
      .groups = "drop"
    )
  
  # **Adjusting Normalization Strategy**
  if (!is.null(full_genome_tai)) {
    # If a full genome tAI dataset is provided, normalize using its max value
    max_tAI_value <- max(full_genome_tai$raw_tAI, na.rm = TRUE)
  } else {
    # Otherwise, use the max within the current dataset (same as before)
    max_tAI_value <- max(gene_tai$raw_tAI, na.rm = TRUE)
  }
  
  # Normalize tAI values
  gene_tai <- gene_tai %>%
    mutate(tAI = raw_tAI / max_tAI_value)
  
  return(gene_tai)
}

# Run the function
organism_tAI <- calculate_tai_with_rscu(organism_rscu_long, organism_tRNA_dict, all_tRNA_dictionaries)
write.csv(organism_tAI, file = "/file_path/organism_tai.csv")

# ALTERNATIVE FUNCTION if comparing multiple organisms, assures scores are normalized within genomes -----------
calculate_tai_with_rscu_per_genome <- function(rscu_data, anticodon_dictionary, full_genome_tai = NULL, min_weight = 0.01) {
  required_columns <- c("Codon", "RSCU", "Sequence")
  if (!all(required_columns %in% colnames(rscu_data))) {
    stop(paste("RSCU data must contain the following columns:", paste(required_columns, collapse = ", ")))
  }
  
  # Compute tAI weights for each codon
  rscu_data$tAI_weight <- mapply(function(codon, rscu_value) {
    anticodon <- reverse_complement(codon)
    
    # Find all matching anticodon scores
    matching_scores <- anticodon_dictionary[names(anticodon_dictionary) == anticodon]
    
    if (length(matching_scores) > 0) {
      # Sum all matching scores (allow for multiple tRNA genes for the same anticodon)
      norm_tRNA_score <- sum(matching_scores)
      tAI_weight <- (rscu_value + 1e-6) * pmax(norm_tRNA_score, min_weight)
    } else {
      # Assign fallback weight if anticodon is missing
      tAI_weight <- min_weight
    }
    
    return(tAI_weight)
  }, rscu_data$Codon, rscu_data$RSCU, SIMPLIFY = TRUE)
  
  # Ensure no zero values (avoid log(0))
  rscu_data$tAI_weight[rscu_data$tAI_weight == 0] <- min_weight
  
  # Calculate tAI for each gene
  gene_tai <- rscu_data %>%
    group_by(Sequence) %>%
    summarize(
      num_matched_codons = sum(tAI_weight > min_weight, na.rm = TRUE),  # Count matched codons
      total_codons = n(),  # Count total codons
      raw_tAI = exp(sum(log(pmax(tAI_weight, min_weight)), na.rm = TRUE) / num_matched_codons),  # Adjust log scaling
      .groups = "drop"
    )
  
  # **Adjusting Normalization Strategy**
  if (!is.null(full_genome_tai)) {
    # If a full genome tAI dataset is provided, normalize using its max value
    max_tAI_value <- max(full_genome_tai$raw_tAI, na.rm = TRUE)
  } else {
    # Otherwise, use the max within the current dataset (same as before)
    max_tAI_value <- max(gene_tai$raw_tAI, na.rm = TRUE)
  }
  
  # Normalize tAI values
  gene_tai <- gene_tai %>%
    mutate(tAI = raw_tAI / max_tAI_value)
  
  return(gene_tai)
}

# OPTIONAL: Average tAI values for an organism to compare tAI across different organisms ---------------
average_tai <- function(tai_values) {
  mean_value <- mean(tai_values[[4]], na.rm = TRUE)
  return(mean_value)
}

organism_avg <- average_tai(organism_tAI)

# plot with simple bar graph
library(ggplot2)

ggplot(organism_tAI, aes(x = Sequence, y = tAI)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Organism tAI Values", x = "Sequence", y = "tAI Score") +
  scale_y_continuous(limits = c(0, 2000)) + #adjust upper limit according to your data
  theme_minimal()

# STATISTICAL ANALYSES -----------------
# step 1 - load tAI files

organism1_tAI <- read.csv("/your/file/path/organism1_tAI.csv")
organism2_tAI <- read.csv("/your/file/path/organism2_tAI.csv")

# step 2 - check normality of data with the Shapiro-Wilk test
shapiro.test(organism1_tAI$tAI)
shapiro.test(organism2_tAI$tAI)


# step 3 - define comparisons:
# compare two groups at the genome level - organism 1 and 2
# multiple-group comparisons - compare all gene groups
# cross-genome comparisons - compare at the group level between organism 1 and organism 2


# step 4 - compare two groups using Mann-Whitney U test (Wilcoxon Rank-Sum test) --------
# format data
organism1$Genome <- "Organism 1"
organism2$Genome <- "Organism 2"
combined_genome_tAI <- rbind(organism1_tAI, organism2_tAI)
combined_genome_tAI$Genome <- factor(combined_genome_tAI$Genome)
str(combined_genome_tAI) #check that data is in a data frame
summary(combined_genome_tAI$tAI) #quick summary statistics

# perform Wilcoxon Rank-Sum Test
wilcox.test(tAI ~ Genome, data = combined_genome_tAI)

boxplot(tAI ~ Genome, data = combined_genome_tAI, col = c("slateblue", "pink2"),
        main = "tAI Comparison: Organism 1 vs. Organism 2", xlab = "Genome", ylab = "tAI")

# step 5 - identify the effect size
# goes beyond p-values to help quantify the difference between the organisms
install.packages("effsize")
library(effsize)
cliff.delta(organism1_tAI$tAI, organism2_tAI$tAI)

# step 6 - group comparisons
# formatting
group1_tAI$GeneGroup <- "Group 1"
group2_tAI$GeneGroup <- "Group 2"
group3_tAI$GeneGroup <- "Group 3"
groupN_tAI$GeneGroup <- "Group N"

group_tAI <- rbind(group1_tAI, group2_tAI, group3_tAI, groupN_tAI)
group_tAI$GeneGroup <- factor(group_tAI$GeneGroup)

# check distribution of tAI across all gene groups
library(dplyr)
group_tAI %>%
  group_by(GeneGroup) %>%
  summarize(n = n(), mean_tAI = mean(tAI), median_tAI = median(tAI), sd_tAI = sd(tAI))

group_tAI$GeneGroup <- factor(
  group_tAI$GeneGroup,
  levels = c("Group 1", "Group 2", "Group 3", "Group N")
)
boxplot(tAI ~ GeneGroup, data = group_tAI, col = rainbow(length(unique(group_tAI$GeneGroup))),
        main = "tAI Distribution by Gene Group", xlab = "Gene Group", ylab = "tAI", las = 1)

# assuming samples are not normally distributed, Kruskal-Wallis test will be used
# (if your samples have a normal distribution, adjust accordingly)
kruskal.test(tAI ~ GeneGroup, data = group_tAI)

# if test shows significance between groups, proceed with pairwise analysis

# non-parametric post-hoc test: Dunn's test
install.packages("FSA")
library(FSA)
dunnTest(tAI ~ GeneGroup, data = group_tAI, method = "bonferroni")

# check effect size
library(effsize)

cliff.delta(group_tAI$tAI[group_tAI$GeneGroup == "Group 1"], 
            group_tAI$tAI[group_tAI$GeneGroup == "Group 2"])

cliff.delta(group_tAI$tAI[group_tAI$GeneGroup == "Group 1"], 
            group_tAI$tAI[group_tAI$GeneGroup == "Group 3"])

cliff.delta(group_tAI$tAI[group_tAI$GeneGroup == "Group 1"], 
            group_tAI$tAI[group_tAI$GeneGroup == "Group N"])

cliff.delta(group_tAI$tAI[group_tAI$GeneGroup == "Group 2"], 
            group_tAI$tAI[group_tAI$GeneGroup == "Group 3"])
# etc; repeat as many times as necessary








