# analysis of tRNA mutations - comparing host and phage matched tRNAs ----------
# make sure .csv file is formatted with the following columns
# "tRNA_ID"	"Pair_ID"	"Organism"	"Isotype"	"Anticodon"	"Anticodon_Loop"	"Full_Sequence"	"Consistent"
# "Consistent" is a Yes or No column that is based on output data from tRNAscan-SE 2.0

# please be sure to install RNAfold (ViennaRNA) FIRST
# if using a Windows system, ViennaRNA can be downloaded from "https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/install.html"
# if using a Mac or Linux system, open the terminal and load Homebrew using the commands: brew install viennarna

# Checking for mutations and saving analysis results
# make sure each tRNA you are analyzing is part of a pair, this script compares tRNAs of the same isotype/anticodon to identify mutations
# Load necessary libraries
library(readr)
library(dplyr)
library(stringr)

# Load your table
tRNA_data <- read_csv("/your/file/path/tRNAmutations.csv")

# Filter only rows with *organism 1* and *organism 2* tRNA pairs
paired_tRNAs <- tRNA_data %>%
  group_by(Pair_ID) %>%
  filter(n() == 2) %>%  # Keeps only properly paired rows
  ungroup()

# Split the dataset for comparison
org1_tRNA <- paired_tRNAs %>% filter(Organism == "Organism 1")
org2_tRNA <- paired_tRNAs %>% filter(Organism == "Organism 2")

# Ensure both datasets are aligned properly
comparison_data <- org1_tRNA %>%
  select(Pair_ID, Isotype, Anticodon, Anticodon_Loop) %>%
  rename(org1_Anticodon_Loop = Anticodon_Loop) %>%
  left_join(org2_tRNA %>% 
              select(Pair_ID, Anticodon_Loop) %>%
              rename(org2_Anticodon_Loop = Anticodon_Loop),
            by = "Pair_ID")

# Function to compare anticodon loops and find mutations --------------
compare_sequences <- function(seq1, seq2) {
  # Check if either sequence is missing
  if (is.na(seq1) | is.na(seq2)) {
    return("Missing sequence")
  }
  # Handle length mismatches
  if (nchar(seq1) != nchar(seq2)) {
    return("Length mismatch")
  }
  # Identify mutation positions
  mismatches <- which(strsplit(seq1, "")[[1]] != strsplit(seq2, "")[[1]])
  if (length(mismatches) == 0) {
    return("No mutations")
  } else {
    # Show positions and changes
    changes <- paste0("Pos ", mismatches, ": ", substr(seq1, mismatches, mismatches), 
                      "→", substr(seq2, mismatches, mismatches), collapse = "; ")
    return(changes)
  }
}

# check for NAs and remove them
comparison_data <- comparison_data %>%
  filter(!is.na(org1_Anticodon_Loop) & !is.na(org2_Anticodon_Loop))
# check for NULL values and remove them
comparison_data <- comparison_data %>%
  mutate(
    org1_Anticodon_Loop = as.character(org1_Anticodon_Loop),
    org2_Anticodon_Loop = as.character(org2_Anticodon_Loop)
  )
#check data
head(comparison_data$org1_Anticodon_Loop)
head(comparison_data$org2_Anticodon_Loop)

# Apply the function to find mutations
comparison_data <- comparison_data %>%
  rowwise() %>%
  mutate(Mutation_Details = compare_sequences(org1_Anticodon_Loop, org2_Anticodon_Loop)) %>%
  ungroup()

# Save the results as a CSV
write_csv(comparison_data, "/your/file/path/anticodon_loop_mutations.csv")

# Show the final mutation table
print(comparison_data)

# determining biological significance - identifying patterns in mutations ------------
comparison_data %>%
  group_by(Isotype) %>%
  summarise(Mutation_Count = sum(Mutation_Details != "No mutations")) %>%
  arrange(desc(Mutation_Count)) %>%
  print()

# find the most frequently mutated positions
mut_pos <- str_extract_all(comparison_data$Mutation_Details, "Pos \\d+") %>%
  unlist() %>%
  table() %>%
  sort(decreasing = TRUE)
print(mut_pos)

# visualize mutation frequency with a bar plot
library(ggplot2)
mutation_summary <- comparison_data %>%
  filter(Mutation_Details != "No mutations") %>%
  mutate(Mutation_Count = str_count(Mutation_Details, "Pos")) %>%
  group_by(Isotype) %>%
  summarise(Total_Mutations = sum(Mutation_Count))

ggplot(mutation_summary, aes(x = Isotype, y = Total_Mutations)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Mutation Frequency by Isotype", x = "tRNA Isotype", y = "Total Mutations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#plotting mutations by nucleotide position ---------------
library(tidyverse)

# Manually define all positions (1-7) and set known values
mut_pos_df <- tibble(
  Mutation_Position = c(1, 2, 3, 4, 5, 6, 7),  # Explicitly define positions 1-7 of the anticodon loop
  Frequency = c(0, 0, 0, 0, 0, 0, 0)  # replace 0s with your actual values
)

# Print to verify
print(mut_pos_df)

# Bar plot of mutation positions
ggplot(mut_pos_df, aes(x = as.factor(Mutation_Position), y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Mutation Frequency Across Anticodon Loop Positions",
       x = "Mutation Position",
       y = "Frequency of Mutations") +
  theme_minimal() +
  scale_x_discrete(limits = as.character(1:7)) +  # Ensure positions 1-7 are displayed
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability


# time to run RNAfold ---------------
# this will calculate the free energy (∆G) of your tRNAs and calculate the difference in ∆G (∆∆G) between pairs of tRNAs
# RNAfold should already be loaded to your computer
# Load necessary libraries
library(dplyr)
library(readr)

# Define function to run RNAfold and extract free energy (ΔG)
run_rnafold <- function(sequence) {
  temp_file <- tempfile(fileext = ".fa")
  writeLines(paste0(">seq\n", sequence), temp_file)
  # Run RNAfold with full path and capture output
  output <- system(paste("/usr/local/bin/RNAfold < ", temp_file), intern = TRUE)
  # Print the raw output for debugging
  print("Full RNAfold Output:")
  print(output)
  # Ensure output has at least 3 lines
  if (length(output) < 3) {
    return(NA)
  }
  # Extract last line containing structure and ΔG
  last_line <- output[length(output)]
  # Split by " (" to isolate ΔG
  split_line <- strsplit(last_line, " \\(")[[1]]
  # Get last element and remove closing parenthesis
  free_energy <- as.numeric(gsub("\\)", "", split_line[length(split_line)]))
  return(free_energy)
}

# optional - Test with a sample sequence
test_sequence <- "GGGAAAUCC"
run_rnafold(test_sequence)

# Load your tRNA data (update with correct file path)
tRNA_data <- read_csv("/your/file/path/tRNAmutations.csv")

# Run RNAfold for each tRNA and store ΔG values
tRNA_data <- tRNA_data %>%
  rowwise() %>%
  mutate(RNAfold_DeltaG = run_rnafold(Full_Sequence)) %>%
  ungroup()

# Save results
write_csv(tRNA_data, "/your/file/path/tRNA_RNAfold_results.csv")

# View results
print(tRNA_data)

# now time to visualize differences in ∆G between organisms/tRNA pairs - this is ∆∆G -----------
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)

# Load tRNA ΔG results
tRNA_data <- read_csv("/your/file/path/tRNA_RNAfold_Results.csv")

paired_tRNAs <- tRNA_data %>%
  select(Pair_ID, Organism, Isotype, Anticodon, RNAfold_DeltaG) %>%
  pivot_wider(names_from = Organism, values_from = RNAfold_DeltaG) %>%
  rename(org1_DeltaG = Organism1, org2_DeltaG = Organism2)

# Calculate ΔΔG
paired_tRNAs <- paired_tRNAs %>%
  mutate(Delta_DeltaG = org1_DeltaG - org2_DeltaG)

# Save the results
write_csv(paired_tRNAs, "/your/file/path/tRNA_DeltaG_comparison.csv")

# Print results
print(paired_tRNAs)

# visualize and run statistics -----------------
# Boxplot of ΔG by Organism**
ggplot(tRNA_data, aes(x = Organism, y = RNAfold_DeltaG, fill = Organism)) +
  geom_boxplot() +
  labs(title = "Comparison of tRNA ΔG Between Org1 and Org2",
       y = "ΔG (kcal/mol)", x = "Organism") +
  theme_minimal()

t_test_pair <- t.test(RNAfold_DeltaG ~ Organism, data = tRNA_data)
print(t_test_pair)

# Scatterplot of Paired tRNA ΔG**
ggplot(paired_tRNAs, aes(x = Sflex_DeltaG, y = Sf14_DeltaG, label = Isotype)) +
  geom_point(aes(color = Isotype), size = 4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +  # Reference line (y = x)
  geom_text(vjust = 1.5, size = 3) +
  labs(title = "Comparison of tRNA ΔG Values: Organism 1 vs. Organism 2",
       x = "Organism 1 ΔG (kcal/mol)", y = "Organism 2 ΔG (kcal/mol)") +
  theme_minimal()

# visualizing tRNAs from the RNAfold mutation data ------------
# run RNAfold and extract the dot-bracket structure to retrieve ∆G value
run_rnafold_structure <- function(sequence) {
  temp_file <- tempfile(fileext = ".fa")
  writeLines(paste0(">seq\n", sequence), temp_file)
  # Run RNAfold to get secondary structure
  output <- system(paste("/usr/local/bin/RNAfold < ", temp_file), intern = TRUE)
  # Print raw output for debugging
  print("Full RNAfold Output:")
  print(output)
  # Extract the secondary structure line (always last)
  structure_line <- output[length(output)]
  # Extract the structure and ΔG
  structure <- gsub("\\s*\\(.*\\)", "", structure_line)  # Remove ΔG value
  deltaG <- as.numeric(gsub(".*\\((-?[0-9.]+)\\).*", "\\1", structure_line))  # Extract ΔG
  return(list(structure = structure, deltaG = deltaG))
}

# **Run for an example tRNA sequence**
test_sequence <- "GCGGGTATAGAGAAAGGGTGTCTCACATGTCTCATTAGCATGTTATCGGTAGGTTCGACTCCTACACCCGCCTCCA"
test_structure <- run_rnafold_structure(test_sequence)

# repeat for your sequences
seq1 <- "insert tRNA sequence here"
seq1_structure <- run_rnafold_structure(seq1)

seq2 <- "insert tRNA sequence here"
seq2_structure <- run_rnafold_structure(seq2)

# continue for however many sequences you wish to investigate
seqN <- "insert tRNA sequence here"
seqN_structure <- run_rnafold_structure(seqN)


# analyzing impact of mutations on tRNA stability ------
# mutation table containing information on: mutation type, mutation position, isotype, anticodon, and ∆G/∆∆G values
library(tidyverse)

# Load CSV file (update with your actual file path)
mutation_data <- read_csv("/your/file/path/mutation_details.csv")

# Convert Mutation_Position to numeric (ignore empty cells)
mutation_data <- mutation_data %>%
  mutate(Mutation_Position = as.numeric(Mutation_Position))

print(mutation_data)

# plot mutation frequency by anticodon loop position -------
# Count mutations at each position
mutation_counts <- mutation_data %>%
  drop_na(Mutation_Position) %>%
  count(Mutation_Position)

# Bar plot of mutation frequency
ggplot(mutation_counts, aes(x = Mutation_Position, y = n)) +
  geom_bar(stat = "identity", fill = "blue", alpha = 0.7) +
  geom_text(aes(label = n), vjust = -0.5, size = 5) +
  scale_x_continuous(breaks = seq(min(mutation_counts$Mutation_Position, na.rm = TRUE), 
                                  max(mutation_counts$Mutation_Position, na.rm = TRUE), 1)) +
  labs(title = "Frequency of Mutations by Position",
       x = "Mutation Position",
       y = "Mutation Count") +
  theme_minimal()

# run comparison of mutation position to changes in ∆G
mutation_data <- mutation_data %>%
  mutate(Mutation_Position = as.numeric(Mutation_Position))

mutation_data <- mutation_data %>%
  mutate(Mutation_Position = ifelse(is.na(Mutation_Position), 0, Mutation_Position))

mutation_data_clean <- mutation_data %>%
  drop_na(Mutation_Position, delta_deltaG)

library(ggrepel)
ggplot(mutation_data_clean, aes(x = Mutation_Position, y = delta_deltaG, label = Isotype)) +
  geom_point(size = 4, color = "slateblue") +
  geom_text_repel(size = 5, max.overlaps = 10) +
  labs(title = "Impact of Mutation Position on ΔΔG",
       x = "Mutation Position",
       y = "ΔΔG (Organism 1 - Organism 2)") +
  theme_minimal()

# statistical analysis -----------
anova_result <- aov(delta_deltaG ~ factor(Mutation_Position), data = mutation_data)
summary(anova_result)

# further statistical analyses comparing mutated to non-mutated tRNAs ---------
# also compares tRNAscan-SE 2.0 consistent vs inconsistent tRNAs
# this comparison focuses on comparing different tRNAs for one organism, does NOT focus on analyzing differences between organisms

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)

# Load tRNA ΔG results
# .csv file should have the following columns: "Pair_ID", "Isotype", "Anticodon", "RNAfold_DeltaG", "Mutation", "Consistent"
tRNA_data <- read_csv("/your/file/path/org1_RNAfold.csv")

# Ensure we are only analyzing one organism's tRNAs
org1_tRNAs <- tRNA_data %>%
  filter(Organism == "Organism 1") %>%
  select(Pair_ID, Isotype, Anticodon, RNAfold_DeltaG, Mutation, Consistent)

# Group 1: Compare ΔG of mutated vs. non-mutated tRNAs
org1_tRNAs <- org1_tRNAs %>%
  mutate(Mutation = ifelse(Mutation == "Yes", "Mutation", "No Mutation"))

# Boxplot of ΔG for mutated vs. non-mutated tRNAs
ggplot(org1_tRNAs, aes(x = Mutation, y = RNAfold_DeltaG, fill = Mutation)) +
  geom_boxplot() +
  labs(title = "Comparison of ΔG for Mutated vs. Non-Mutated tRNAs",
       x = "Mutation Status",
       y = "ΔG (kcal/mol)") +
  theme_minimal()

# Perform a statistical test (t-test or Wilcoxon test based on distribution)
t_test_mutation <- t.test(RNAfold_DeltaG ~ Mutation, data = org1_tRNAs)
print(t_test_mutation)


# Group 2: Compare ΔG of tRNAs based on tRNAscan-SE consistency
org1_tRNAs <- org1_tRNAs %>%
  mutate(Consistency = ifelse(Consistent == "Yes", "Consistent", "Inconsistent"))

# Boxplot of ΔG for consistent vs. inconsistent tRNAs
ggplot(org1_tRNAs, aes(x = Consistency, y = RNAfold_DeltaG, fill = Consistency)) +
  geom_boxplot() +
  labs(title = "Comparison of ΔG for tRNAscan-SE Consistent vs. Inconsistent tRNAs",
       x = "tRNAscan-SE Consistency",
       y = "ΔG (kcal/mol)") +
  theme_minimal()

# Perform a statistical test (t-test or Wilcoxon test based on distribution)
t_test_consistency <- t.test(RNAfold_DeltaG ~ Consistency, data = org1_tRNAs)
print(t_test_consistency)

# Save results
write_csv(org1_tRNAs, "/your/file/path/org1_tRNA_DeltaG_analysis.csv")


