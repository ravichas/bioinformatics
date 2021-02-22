getwd() # Get the working directory, and use setwd() to
# change it to any preferred location

# Note the installation instructions for BioConductor has changed
# Refer to the documentation pages for Bioconductor here,
# https://bioconductor.org/install/ and here
# https://bioconductor.org/packages/release/bioc/html/Biostrings.html

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

library(Biostrings)# Install the Biostrings library

data(BLOSUM50) # load the data for the BLOSUM50 matrix
BLOSUM50[1:4,1:4] # view the first four rows and columns

# of this matrix
nw <- pairwiseAlignment(AAString("PAWHEAE"),
                          AAString("HEAGAWGHEE"), substitutionMatrix = "BLOSUM50",
                          gapOpening = 0, gapExtension = -8) 

nw_p30 <- pairwiseAlignment(AAString("PAWHEAE"),
                        AAString("HEAGAWGHEE"), substitutionMatrix = "PAM30",
                        gapOpening = 0, gapExtension = -8) 

# create object nw aligning
# two amino acid strings with the specified matrix and gap
# penalties

nw # view the result. Try repeating this alignment with

# different gap penalties and scoring matrices. Biostrings
# includes 10 matrices (PAM30 PAM40, PAM70, PAM120, PAM250,
# BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, and BLOSUM100).

# compareStrings(nw) # view the consensus sequence
# consensusMatrix(nw)
