# codon-usage
Relative synonymous codon usage (RSCU) calculation, tRNA adaptation index (tAI) calculation, and analysis of tRNA mutations

This repository contains 3 scripts for computational analyses of codon usage and tRNA adaptation between organisms and gene groups. The scripts are general and should be modified to fit your own data.

The "formal codon usage" script contains functions for extracting FASTA file data, calculating total GC, GC1, GC2, and GC3 content at the gene level, plotting these values, and running statistical analyses on the resulting data. It also contains functions for calculating relative synonymous codon usage (RSCU) at the gene level, plotting the data in heatmaps and PCA plots, and for statistical analyses of the resulting RSCU data.

The "formal tAI" script contains functions to extract tRNAscan-SE 2.0 data from excel files, calculate codon usage as RSCU, and reformat the RSCU data to be appropriate for analysis of tRNA adaptation index (tAI). tAI is calculated from the tRNAscan data and the RSCU data, and is organism-specific so multiple types of organisms can be analyzed together. This script also includes functions to plot the data, an option to average the values, and statistical analyses for the resulting data.

tAI calculations are as follows:
Calculation of codon weights
tAIweight = (RSCU + ϵ) x norm tRNA score![image](https://github.com/user-attachments/assets/d191876f-c342-4a5e-841b-106569e6682a)
where ϵ is a small value (1 x 10-6) to prevent division by zero during log calculations (Schluter 1988), RSCU is the relative synonymous codon usage value for a given codon in a given genome, and norm tRNA score is the normalized tRNA score for the complementary anticodon. Any missing anticodons are assigned a value of 0.01

Calculating geometric means of all tAI weights
tAI = exp(log(tAIweight)/N)![image](https://github.com/user-attachments/assets/f000099b-3da0-45fd-8ed7-bcf5f98a62d9)
where N is the number of matched codons for a given anticodon.

Normalizing tAI values across a genome
tAIfinal = tAI/max(tAI)![image](https://github.com/user-attachments/assets/a074920a-f6e9-407d-813c-1f7b6b2291cc)
where tAI is the calculated tRNA adaptation index value for a given gene and max(tAI) is the maximum tAI across all genes in a given genome

Finally, the "formal mutation" script contains functions to identify anticodon loop mutations, run RNAfold to calculate ∆G values for tRNAs, and plot mutation data by tRNA isotype, mutation position, or ∆G values. Statistical analyses functions are also included.

If you have any questions or run into any issues, please feel free to email me at nicole.ross@ufl.edu
