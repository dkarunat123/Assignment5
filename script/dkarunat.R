#-----------------------------------------------------------------------------#
## 1. Loading Dependencies: ======
#-----------------------------------------------------------------------------#

# List of CRAN packages
cran_packages <- c("progress", "beepr", "tidyverse", "stringr", "ape", "phytools")

# List of Bioconductor packages
bioc_packages <- c("Biostrings", "DECIPHER", "msa", "phangorn", "ips")

# Load CRAN packages
for (pkg in cran_packages) {
  library(pkg, character.only = TRUE)
}

# Install Bioconductor packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}


#-----------------------------------------------------------------------------#
## 2. NCBI Data Acquisition of Mustelidae COI and RAG1 Genes: ======
#-----------------------------------------------------------------------------#

# Create a function to fetch COI and RAG1 gene sequences for Mustelidae

fetch_sequences <- function(gene_name, length_range, batch_size = 50, output_file = NULL) {
  length_filter <- paste0(" AND ", length_range[1], ":", length_range[2], "[SLEN]")
  search_query <- entrez_search(
    db = "nuccore",
    term = paste0("Mustelidae[ORGN] AND ", gene_name, "[GENE]", length_filter),
    use_history = TRUE
  )
  
  total_hits <- search_query$count
  if (is.null(total_hits) || total_hits == 0) stop("No sequences found for the query.")
  
  message(paste("Total sequences found for", gene_name, ":", total_hits))
  batch_count <- ceiling(total_hits / batch_size)
  pb <- progress_bar$new(
    format = paste0("Fetching ", gene_name, " sequences [:bar] :percent in :elapsed"),
    total = batch_count,
    clear = FALSE,
    width = 60
  )
  
  if (!is.null(output_file)) {
    file_conn <- file(output_file, "w")
  }
  
  for (start in seq(0, total_hits - 1, by = batch_size)) {
    pb$tick()
    batch_data <- tryCatch({
      entrez_fetch(
        db = "nuccore",
        web_history = search_query$web_history,
        rettype = "fasta",
        retmax = batch_size,
        retstart = start
      )
    }, error = function(e) {
      message("Error fetching batch starting at ", start, ": ", e$message)
      return("")
    })
    
    if (!is.null(output_file) && nzchar(batch_data)) {
      write(batch_data, file_conn)
    }
  }
  
  if (!is.null(output_file)) {
    close(file_conn)
    dna_sequences <- readDNAStringSet(output_file, format = "fasta")
  } else {
    dna_sequences <- NULL
  }
  
  beepr::beep("ping")
  return(dna_sequences)
}


#-----------------------------------------------------------------------------#

coi_file <- "../data/Combined_FASTA_COI.fasta"
rag1_file <- "../data/Combined_FASTA_RAG1.fasta"


#Chosen Mustelidae sequence length ranges for COI & RAG1:

    # COI: This is the primary marker for DNA barcoding, known to have a typical length of 657 bps (Nugent et al., 2020). Therefore, a length range of 600-700 will likely be appropriate for the sequence search.
    # RAG1: While the gene itself is longer,  shorter regions are typically amplified and sequenced for phylogenetic studies due to their conserved regions. These shorter regions typically seem to be around the 1000-1100 bp length, evidenced by several studies using this gene: (Liu et al., 2023), (Sato et al., 2004).


# # Fetch COI sequences (600-700 bp)
# Mustelidae_COI <- fetch_sequences("COI", length_range = c(600, 700), output_file = coi_file)
# 187 sequences for the Mustelidae COI gene were retrieved

# # Fetch RAG1 sequences (900-1200 bp)
# Mustelidae_RAG1 <- fetch_sequences("RAG1", length_range = c(1000, 1100), output_file = rag1_file)
# 100 sequences for the Mustelidae RAG1 gene were retrieved


#-----------------------------------------------------------------------------#
## 3. Analyze Initial Sequence Lengths: ======
#-----------------------------------------------------------------------------#

# Create function to analyze initial sequence lengths by outputting a summary and a histogram for each gene

analyze_sequence_lengths <- function(fasta_files, titles, colors) {
  if (length(fasta_files) != length(titles) || length(fasta_files) != length(colors)) {
    stop("Lengths of fasta_files, titles, and colors must match.")
  }
  
  # Setup for side-by-side plots
  par(mfrow = c(1, length(fasta_files)), mar = c(4, 4, 2, 1))  # Adjust margins
  
  # Initialize a list to store sequence lengths
  sequence_data <- list()
  
  for (i in seq_along(fasta_files)) {
    # Read and analyze each file
    dna_sequences <- readDNAStringSet(fasta_files[i], format = "fasta")
    sequence_lengths <- nchar(as.character(dna_sequences))
    sequence_data[[i]] <- sequence_lengths
    
    # Print summary
    cat("\nSummary of", titles[i], "Lengths:\n")
    print(summary(sequence_lengths))
    
    # Plot histogram
    hist(
      sequence_lengths,
      main = paste(titles[i], "Lengths"),
      xlab = "Length (bp)",
      breaks = 20,
      col = colors[i],
      border = "white",
      ylim = c(0, 100)
    )
  }
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
  
  return(sequence_data)
}


# File paths for COI and RAG1
fasta_files <- c(coi_file, rag1_file)
titles <- c("Figure 1a: COI", "Figure 1b: RAG1")
colors <- c("#689ACD", "#EC5F67")

# Analyze sequences and display histograms side by side
sequence_lengths <- analyze_sequence_lengths(fasta_files, titles, colors)

#Figure 1a:
    # The COI sequences show more variability in their lengths. Most sequences are clustered between 660-700 bp, with a clear peak around the upper end of this range. A smaller number of sequences fall below this cluster.
    # The presence of shorter sequences (e.g., around 620–640 bp) suggests that some of the COI sequences might be incomplete or represent partial coding sequences.
    # This variability is expected because COI is a mitochondrial gene, which often evolves rapidly and may contain insertions, deletions, or incomplete sequences in public datasets.

# Figure 1b:
    # The RAG1 sequences are much longer (approximately 1030–1090 bp) and display a narrower length range, reflecting their more conserved nature as a nuclear gene.
    # The relatively uniform distribution of sequence lengths highlights the conserved nature of RAG1 across species, with fewer partial or incomplete sequences compared to COI.
    # The lower overall counts in some bins suggest fewer retrieved sequences or higher variability in data availability for this gene.

#-----------------------------------------------------------------------------#
## 4. Cleaning and Filtering: Mustelidae Gene Data ======
#-----------------------------------------------------------------------------#

# Function to create a tidy data frame from a FASTA file
create_sequence_dataframe <- function(fasta_file) {
  # Read the FASTA file
  dna_sequences <- readDNAStringSet(fasta_file, format = "fasta")
  
  # Convert to a tibble using tidyverse
  df <- tibble(
    Header = names(dna_sequences),      # Extract headers
    Sequence = as.character(dna_sequences)  # Extract sequences
  )
  
  return(df)
}


# Create data frames for COI and RAG1 genes
df_COI <- create_sequence_dataframe(coi_file)
df_RAG1 <- create_sequence_dataframe(rag1_file)

# View the data frames
print(df_COI)
print(df_RAG1)

#-----------------------------------------------------------------------------#

# Remove number/letter sequence and a space from the start of the Header columns for each gene's dataframe
df_COI <- df_COI %>%
  mutate(Header = str_remove(Header, "^[A-Za-z0-9\\.]+ "))

df_RAG1 <- df_RAG1 %>%
  mutate(Header = str_remove(Header, "^[A-Za-z0-9\\.]+ "))

# View the cleaned data frames
print(df_COI)
print(df_RAG1)

#-----------------------------------------------------------------------------#

# Split the Header column into Species and Details for COI
df_COI <- df_COI %>%
  mutate(
    Species = word(Header, 1, 2),   # Extract the first three words
    Details = word(Header, 3, -1)  # Extract everything after the first three words
  ) %>%
  select(Species, Details, Sequence)

# Split the Header column into Species and Details for RAG1
df_RAG1 <- df_RAG1 %>%
  mutate(
    Species = word(Header, 1, 2),   # Extract the first three words
    Details = word(Header, 3, -1)  # Extract everything after the first three words
  ) %>%
  select(Species, Details, Sequence)

# View the transformed data frames
head(df_COI)
head(df_RAG1)

#-----------------------------------------------------------------------------#

# Examine the Details column of both dataframes to check for important information that could impact cleaning, otherwise, this column can be disregarded from now on, as it will not be necessary for downstream analysis.

# View the the Details column from df_COI
#head(df_COI$Details, 187)

# View the Details column from df_RAG1
#head(df_RAG1$Details, 100)


# Project Relevance Criteria for Sequences:
    # The sequences must belong to the Mustelidae family.
    # For COI, only mitochondrial COI sequences are relevant. 
    # For RAG1, only nuclear RAG1 sequences are relevant. 
    # Non-target gene sequences, nuclear copies of mitochondrial genes (NUMTs), or ambiguous designations should be excluded.
    # Vouchers, isolates, or strains should explicitly correspond to Mustelidae species and not be generic or unidentified sequences.
    # Partial sequences or sequences marked as ambiguous or problematic may need to be excluded as well.

# Nuclear copies of mitochondrial genes (NUMTs) are fragments of mitochondrial DNA that have been transferred to the nuclear genome over evolutionary time (Xue et al., 2023). While they resemble mitochondrial sequences, they are no longer part of the mitochondrial genome and may have accumulated mutations or deletions, leading to incorrect phylogenetic inferences (Xue et al., 2023). Using NUMTs in mitochondrial phylogenies can introduce errors by misrepresenting the true evolutionary relationships derived from the mitochondrial genome.

# Using partial coding sequences (partial CDS) is generally acceptable for phylogenetic analysis, provided the sequences are sufficiently conserved and alignable across species. However, partial CDS sequences may limit the resolution of phylogenetic trees compared to full-length sequences. This is because they contain less genetic information, which might reduce the statistical power to detect subtle evolutionary patterns or relationships (Czech et al., 2022). Therefore, while partial CDS can be used, we need to ensure that the sequences are aligned properly and are representative of the regions being studied (COI/RAG1) to avoid biases in the analysis.

# So, I will remove rows that contain NUMTs in the COI dataframe as they would obstruct downstream analysis.

# There are 187 rows in the COI dataframe currently

# Filter out rows with "nuclear copy of mitochondrial gene" in the Details column for df_COI
df_COI <- df_COI %>%
  filter(!str_detect(Details, "nuclear copy of mitochondrial gene"))

# View the filtered data frames
head(df_COI)

# There are now 185 rows in the COI dataframe. Therefore, two sequences were removed.


#-----------------------------------------------------------------------------#

# Add Sequence_Length column to df_COI
df_COI <- df_COI %>%
  mutate(Sequence_Length = nchar(Sequence))

# Add Sequence_Length column to df_RAG1
df_RAG1 <- df_RAG1 %>%
  mutate(Sequence_Length = nchar(Sequence))

# View the updated data frames
print(df_COI)
print(df_RAG1)


#-----------------------------------------------------------------------------#

# Some species (in both data frames) have multiple associated sequences (i.e., there are multiple columns with the same Species name). I only want unique species in the Species column in order to do the downstream analysis (alignment, phylogenetic tree construction, etc). 

# I will be using the longest available sequence for each unique species because longer sequences contain more phylogenetic signals, improving the reliability of evolutionary analyses (Czech et al., 2022). This ensures more accurate placement of query sequences within the phylogenetic tree (Czech et al., 2022).

#This will remove duplicate sequences to retain one representative per species for each marker.

# For df_COI
df_COI <- df_COI %>%
  group_by(Species) %>%
  slice_max(order_by = Sequence_Length, n = 1, with_ties = FALSE) %>%
  ungroup()

# For df_RAG1
df_RAG1 <- df_RAG1 %>%
  group_by(Species) %>%
  slice_max(order_by = Sequence_Length, n = 1, with_ties = FALSE) %>%
  ungroup()

# Check the number of unique species retained
nrow(df_COI) # 21 unique species
nrow(df_RAG1) # 51 unique species


# Verify uniqueness of the Species column in df_COI
all_unique_COI <- nrow(df_COI) == n_distinct(df_COI$Species)
all_unique_COI # TRUE

# Verify uniqueness of the Species column in df_RAG1
all_unique_RAG1 <- nrow(df_RAG1) == n_distinct(df_RAG1$Species)
all_unique_RAG1 # TRUE

# Remove uneeded objects:
rm(all_unique_COI)
rm(all_unique_RAG1)


#-----------------------------------------------------------------------------#

# Visualizing the cleaned COI and RAG1 sequences using the analyze_sequence_lengths function (defined in 3.)

# # Save the cleaned sequences to FASTA files
# writeXStringSet(DNAStringSet(df_COI$Sequence), "../data/Cleaned_COI.fasta")
# writeXStringSet(DNAStringSet(df_RAG1$Sequence), "../data/Cleaned_RAG1.fasta")

# File paths for cleaned COI and RAG1
cleaned_fasta_files <- c("../data/Cleaned_COI.fasta", "../data/Cleaned_RAG1.fasta")
cleaned_titles <- c("Figure 2a: Cleaned COI", "Figure 2b: Cleaned RAG1")
cleaned_colors <- c("#689ACD", "#EC5F67")

# Visualize the cleaned COI and RAG1 sequences side by side
cleaned_sequence_lengths <- analyze_sequence_lengths(cleaned_fasta_files, cleaned_titles, cleaned_colors)


# Figure 2a:
    # The cleaned histogram shows significantly fewer sequences compared to the original. This indicates that sequences were removed during cleaning, likely due to exclusion criteria (e.g., partial sequences, NUMTs, or non-target sequences).
    # The range of sequence lengths is slightly reduced, clustering more tightly around 660–700 bp. This suggests that shorter, potentially incomplete or ambiguous sequences have been filtered out.
    # The cleaned dataset likely represents high-quality sequences that are more suitable for downstream analysis.

# Figure 2a:
    # The cleaned RAG1 dataset also contains fewer sequences compared to the original. Sequences may have been excluded due to being incomplete, ambiguous, or failing other project-specific criteria.
    # The cleaned sequences are tightly clustered around 1095 bp, reflecting a more conserved and complete set of sequences. This indicates that sequences with high variability or partial lengths (e.g., shorter than 1065 bp) were removed.

# The cleaning process has resulted in datasets with more uniform sequence lengths, indicative of higher quality and completeness. 


#-----------------------------------------------------------------------------#

# Implications of Results:

    # The COI dataset has a slightly more broad range of sequence lengths (628–697 bp), reflecting the variability in partial mitochondrial sequences. However, the lengths are consistent with typical COI sequences and suitable for alignment and 

    # The RAG1 dataset has a slightly more uniform distribution (1063–1095 bp), reflecting the conserved nature of this nuclear gene and its longer sequences, which is ideal for robust phylogenetic inference.


#-----------------------------------------------------------------------------#
## 5. Finding Overlapping Species in COI & RAG1 DataFrames: ======
#-----------------------------------------------------------------------------#

# Find overlapping species
overlapping_species <- intersect(df_COI$Species, df_RAG1$Species)

# View the overlapping species
cat("Overlapping species:\n")
print(overlapping_species)

# Count the number of overlapping species
cat("Number of overlapping species:", length(overlapping_species), "\n") 
# There are 19 overlapping species between the two dataframes.


#-----------------------------------------------------------------------------#

# To simplify the data we will be looking at for the rest of the analysis, a new dataframe will be created to hold only required information.

# Filter both data frames to include only overlapping species
df_COI_filtered <- df_COI %>%
  filter(Species %in% overlapping_species)

df_RAG1_filtered <- df_RAG1 %>%
  filter(Species %in% overlapping_species)

# Create the new dataframe with the specified columns
df_COI_RAG1_Overlap <- df_COI_filtered %>%
  inner_join(df_RAG1_filtered, by = "Species", suffix = c("_COI", "_RAG1")) %>%
  select(Species, Sequence_COI = Sequence_COI, Sequence_RAG1 = Sequence_RAG1)

# View the new dataframe & verify the number of rows is still 19 (indicating that no rows/species were lost)
print(df_COI_RAG1_Overlap) # 19 rows displayed


# Check that this new dataframe holds the correct information by comparing its contents to the rows that hold data for the 19 overlapping species in the df_COI and df_RAG1 dataframes

# Filter original data frames for only overlapping species
df_COI_filtered <- df_COI %>%
  filter(Species %in% overlapping_species)

df_RAG1_filtered <- df_RAG1 %>%
  filter(Species %in% overlapping_species)

# Check that Species columns match
species_match <- all(df_COI_RAG1_Overlap$Species == df_COI_filtered$Species) &&
  all(df_COI_RAG1_Overlap$Species == df_RAG1_filtered$Species)

# Check that Sequence_COI matches df_COI_filtered$Sequence
sequence_coi_match <- all(df_COI_RAG1_Overlap$Sequence_COI == df_COI_filtered$Sequence)

# Check that Sequence_RAG1 matches df_RAG1_filtered$Sequence
sequence_rag1_match <- all(df_COI_RAG1_Overlap$Sequence_RAG1 == df_RAG1_filtered$Sequence)

# Display results
cat("Do Species columns match across all data frames? ", species_match, "\n")
# TRUE
cat("Do Sequence_COI columns match? ", sequence_coi_match, "\n")
# TRUE
cat("Do Sequence_RAG1 columns match? ", sequence_rag1_match, "\n")
# TRUE


#-----------------------------------------------------------------------------#

# Examine the Sequence columns to make sure that only sequences with standard nucleotide bases remain:

    #In sequence alignments, bases such as Y (indicating a pyrimidine) and R (indicating a purine) are examples of IUPAC nucleotide codes, which represent ambiguous or uncertain nucleotide identities (Johnson, 2010). They are commonly used in sequences where multiple nucleotides are possible in that specific position (Johnson, 2010). 

# Define the set of standard nucleotide bases
standard_bases <- c("A", "T", "G", "C")

# Function to check for non-standard bases
check_non_standard_bases <- function(sequence_column) {
  # Identify species with non-standard bases
  non_standard_indices <- sapply(sequence_column, function(seq) {
    any(!strsplit(seq, "")[[1]] %in% standard_bases)
  })
  # Return the species names with non-standard bases
  df_COI_RAG1_Overlap$Species[non_standard_indices]
}


# Check for non-standard bases in Sequence_COI column
non_standard_species_COI <- check_non_standard_bases(df_COI_RAG1_Overlap$Sequence_COI)
non_standard_species_COI # All sequences contain standard bases

# Check for non-standard bases in Sequence_RAG1 column
non_standard_species_RAG1 <- check_non_standard_bases(df_COI_RAG1_Overlap$Sequence_RAG1)
non_standard_species_RAG1 # 6 Sequences contain non-standard bases


# Remove rows from the dataframe with species listed in non_standard_species_RAG1
df_COI_RAG1_Overlap <- df_COI_RAG1_Overlap %>%
  filter(!Species %in% non_standard_species_RAG1)

# View the updated dataframe
print(df_COI_RAG1_Overlap) 
cat("There are now only", nrow(df_COI_RAG1_Overlap), "overlapping species.", "\n") # 13 overlapping species


#-----------------------------------------------------------------------------#
## 6. Sequence Alignment of COI & RAG1 Genes: ======
#-----------------------------------------------------------------------------#

# Alignment is critical for comparative analyses in phylogenetics as it ensures that homologous positions are aligned across sequences (Wright, 2024). 

# For protein-coding genes like COI, codon-based alignment is preferred because it aligns sequences in a way that preserves reading frames (by introducing gaps in multiples of three nucleotides), which is vital for maintaining biological meaning and avoiding frame-shift errors (Wright, 2024). For this reason, the DECIPHER package's AlignTranslation function will be used for the COI sequences, ensuring that codons are properly aligned (Wright, 2024).

# Nuclear markers like RAG1, which are often non-coding or only partially coding, do not have the same strict requirement for reading frame preservation (Wright, 2024). Instead, these sequences benefit from nucleotide-based alignment methods, which maximize similarity across bases without codon constraints (Bonatesta, Kainrath, & Bodenhofer, 2024). For RAG1, the msa package and the MUSCLE algorithm were chosen for their efficiency and accuracy in aligning nucleotide sequences (Bonatesta, Kainrath, & Bodenhofer, 2024).

# The use of these tailored methods ensures that the alignment is both biologically accurate and optimized for the type of data. Misaligned sequences could lead to erroneous phylogenies, so using the appropriate method for each marker is essential for reliable downstream analyses.


#-----------------------------------------------------------------------------#

# Sequence Alignment for COI (Codon-based Alignment):

# Read the COI sequences for overlapping species
coi_sequences <- DNAStringSet(df_COI_RAG1_Overlap$Sequence_COI)

# Assign species names to the sequences
names(coi_sequences) <- df_COI_RAG1_Overlap$Species

# Perform codon-based alignment
alignment_COI <- AlignTranslation(coi_sequences)

# Convert alignment to DNAStringSet and save as a fasta file
aligned_COI_DNA <- as(alignment_COI, "DNAStringSet")
# writeXStringSet(aligned_COI_DNA, "../data/Aligned_COI.fasta")


#-----------------------------------------------------------------------------#

# Sequence Alignment for RAG1 (Nucleotide-based Alignment):

#>_< For large datasets with many sequences or very long sequences, aligning all sequences in a single step can exceed computational resources, leading to errors or program crashes. In this case, attempting to align the entire dataset for RAG1 using MUSCLE resulted in an error (MUSCLE finished by an unknown reason), likely due to memory constraints or limits on handling the large sequence set. 

# Batch Processing: Splitting the sequences into smaller, manageable batches allows each subset to be aligned individually, reducing the computational load at any given time. After aligning the batches, the results are merged into a single dataset, ensuring all sequences are included in the final alignment without overloading resources.


# Read the RAG1 sequences for overlapping species
rag1_sequences <- DNAStringSet(df_COI_RAG1_Overlap$Sequence_RAG1)

# Assign species names to the sequences before alignment
names(rag1_sequences) <- df_COI_RAG1_Overlap$Species

# Split sequences into batches for alignment
batch_size <- 10
sequence_batches <- split(rag1_sequences, ceiling(seq_along(rag1_sequences) / batch_size))

# Perform alignment for each batch
aligned_batches <- lapply(sequence_batches, function(batch) {
  msa(batch, method = "Muscle")
})

# Convert each batch to DNAStringSet
aligned_batches_DNA <- lapply(aligned_batches, function(alignment) {
  as(alignment, "DNAStringSet")
})

# Merge aligned batches into one DNAStringSet
aligned_RAG1 <- DNAStringSet(unlist(lapply(aligned_batches_DNA, as.character)))

# Reassign species names to the merged alignment & save as a fasta file
names(aligned_RAG1) <- df_COI_RAG1_Overlap$Species
# writeXStringSet(aligned_RAG1, "../data/Aligned_RAG1.fasta")


#-----------------------------------------------------------------------------#

# Viewing COI & RAG1 alignments:

# Load the aligned COI sequences
aligned_COI <- readDNAStringSet("../data/Aligned_COI.fasta", format = "fasta")

# View the alignment using BrowseSeqs
BrowseSeqs(aligned_COI)

# Load the aligned RAG1 sequences
aligned_RAG1 <- readDNAStringSet("../data/Aligned_RAG1.fasta", format = "fasta")

# View the alignment using BrowseSeqs
BrowseSeqs(aligned_RAG1)


# In BrowseSeqs alignments (Molecular Phylogenetic Techniques, 2023):
    # Regions with no gaps across sequences indicate conserved regions.
    # Regions with many gaps or variations suggest areas of divergence or lower conservation. Columns with many gaps suggest indels, common in evolutionary events.

# Using BrowseSeqs(), I was able to observe a significantly larger number of gaps in the COI alignment compared to the RAG1 alignment. This result was expected due to the inherent differences in the nature of these genes. The COI gene, being mitochondrial, evolves faster than nuclear genes like RAG1, leading to greater sequence variability across species (Li et al., 2024). 

# It can also be seen that within the COI alignment, certain regions exhibited few to no gaps, indicating conserved regions (Molecular Phylogenetic Techniques, 2023). These conserved areas are likely under strong selective pressure, preserving their sequences due to critical functional or structural roles in the gene's protein product (Liu et al., 2008). The presence of conserved regions amidst variable ones reflects the balance between evolutionary changes and the need to maintain essential gene functions. In contrast, RAG1, being a nuclear gene, evolves more slowly and uniformly across its sequence, resulting in fewer overall gaps and more conservation throughout the alignment. 


# Metrics like Robinson-Foulds distance (which will be used later on) rely on accurate and reliable tree structures. Poorly aligned or regions with a lot of gaps in the COI alignment can lead to artifacts in the tree (Du et al., 2019), which might obscure the true evolutionary signals and inflate incongruences between mitochondrial and nuclear markers.

# Testing evolutionary models requires robust alignments to accurately reflect homologous positions. Retaining noisy regions in COI could skew these analyses.

# Therefore, the COI alignment will be trimmed.


#-----------------------------------------------------------------------------#

#Trimming the COI Alignment:

# Convert COI alignment to a matrix for manual trimming
alignment_matrix <- as.matrix(aligned_COI)

# Calculate gap fractions per column
gap_fractions <- colMeans(alignment_matrix == "-")

# Retain columns with gap fractions below the threshold
threshold <- 0.3 # The threshold for trimming (relatively strict)
trimmed_matrix <- alignment_matrix[, gap_fractions <= threshold]

# Convert back to a DNAStringSet
trimmed_aligned_COI <- DNAStringSet(apply(trimmed_matrix, 1, paste0, collapse = ""))

# # Save the trimmed alignment
# writeXStringSet(trimmed_aligned_COI, "../data/Trimmed_Aligned_COI.fasta")


# The proportion of gaps in each column were calculated, identifying poorly aligned regions or overhangs. Then the columns with gap fractions below 0.5 (50%) were retained, ensuring that regions high in gaps were excluded. Finally, using the filtered columns, the alignment was recreated as a DNAStringSet. By keeping only relatively well-aligned regions, this was meant to increase the reliability of the phylogenetic analyses and reduce noise.


# Viewing New COI Alignment:

# Load the trimmed COI alignment
trimmed_aligned_COI <- readDNAStringSet("../data/Trimmed_Aligned_COI.fasta", format = "fasta")

# View the trimmed alignment using BrowseSeqs()
BrowseSeqs(trimmed_aligned_COI)

#I can see that regions with large numbers of gaps have been removed for the COI alignment compared to the previous COI alignment. 


#-----------------------------------------------------------------------------#
## 7. Model Testing & Constructing ML Phylogenetic Trees for COI & RAG1: ======
#-----------------------------------------------------------------------------#

# Maximum Likelihood (ML) trees are a type of phylogenetic tree constructed by identifying the tree topology that has the highest probability of producing the observed genetic data under a specified evolutionary model (Cho, 2012).

# Load the aligned COI and RAG1 sequences
aligned_COI <- readDNAStringSet("../data/Trimmed_Aligned_COI.fasta", format = "fasta")
aligned_RAG1 <- readDNAStringSet("../data/Aligned_RAG1.fasta", format = "fasta")

# Convert to DNAbin format
aligned_COI_DNAbin <- as.DNAbin(aligned_COI)
aligned_RAG1_DNAbin <- as.DNAbin(aligned_RAG1)

# Convert to phyDat format for ML tree construction
phyDat_COI <- phyDat(aligned_COI_DNAbin, type = "DNA")
phyDat_RAG1 <- phyDat(aligned_RAG1_DNAbin, type = "DNA")

# Generate random starting trees for COI and RAG1
random_tree_COI <- rtree(n = length(phyDat_COI), tip.label = names(phyDat_COI))
random_tree_RAG1 <- rtree(n = length(phyDat_RAG1), tip.label = names(phyDat_RAG1))


#-----------------------------------------------------------------------------#

# Model Testing:

# sink(tempfile()) sink() suppresses output on console

# Temporarily suppress console output
sink(tempfile())
model_test_COI <- modelTest(phyDat_COI, model = "all")
model_test_RAG1 <- modelTest(phyDat_RAG1, model = "all")
sink()  # Restore console output


# View the best models based on AIC and BIC
# print(model_test_COI)
# print(model_test_RAG1)

# Determine the model with the lowest AIC/BIC values:

# For COI
best_model_COI_AIC <- model_test_COI[which.min(model_test_COI$AIC), ]
best_model_COI_BIC <- model_test_COI[which.min(model_test_COI$BIC), ]

# For RAG1
best_model_RAG1_AIC <- model_test_RAG1[which.min(model_test_RAG1$AIC), ]
best_model_RAG1_BIC <- model_test_RAG1[which.min(model_test_RAG1$BIC), ]


# Summary of model testing results:
cat("Best model for COI based on AIC:", best_model_COI_AIC$Model, "\n")
cat("Best model for COI based on BIC:", best_model_COI_BIC$Model, "\n")

cat("Best model for RAG1 based on AIC:", best_model_RAG1_AIC$Model, "\n")
cat("Best model for RAG1 based on BIC:", best_model_RAG1_BIC$Model, "\n")


# Model testing is essential in phylogenetic analyses to identify the evolutionary model that best explains the data while balancing complexity and fit (Susko & Roger, 2020). Both the Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC) were evaluated, with AIC prioritizing model fit and BIC favoring simpler models (Susko & Roger, 2020). For the COI gene, the GTR (TIM2+G(4)) model was chosen as it had the lowest AIC and accounts for rate heterogeneity, crucial for mitochondrial genes that evolve rapidly and unevenly (Schliep & Bardel-Kahr, 2024). For the RAG1 gene, the SYM+I model was selected based on its lowest AIC and inclusion of invariant sites (Schliep et al., 2024) capturing the more uniform evolutionary dynamics of nuclear genes.

# While the BIC suggested alternative models, AIC's ability to account for data-specific complexity made it a better guide for this analysis (Susko & Roger, 2020). Using TIM2+G(4) for COI and SYM+I for RAG1 allows the phylogenetic analyses to reflect the distinct evolutionary patterns of each gene (Susko & Roger, 2020). This tailored approach improves the biological accuracy and reliability of the resulting phylogenetic trees, ensuring meaningful insights into evolutionary relationships.


#-----------------------------------------------------------------------------#

# Optimize ML tree for COI using TIM2+G(4)
ml_tree_COI <- pml(random_tree_COI, data = phyDat_COI)
optimized_ml_tree_COI <- optim.pml(
  ml_tree_COI, 
  model = "GTR", # TIM2 is a subset of GTR
  optGamma = TRUE, 
  optInv = FALSE
)

# Optimize ML tree for RAG1 using SYM+I
ml_tree_RAG1 <- pml(random_tree_RAG1, data = phyDat_RAG1)
optimized_ml_tree_RAG1 <- optim.pml(
  ml_tree_RAG1, 
  model = "SYM", 
  optGamma = FALSE, 
  optInv = TRUE
)


#-----------------------------------------------------------------------------#

# Plot the optimized ML trees:

# Calculate distances from the root for COI
distances_COI <- node.depth.edgelength(optimized_ml_tree_COI$tree)

# Calculate distances for only tip nodes (species)
distances_COI <- distances_COI[1:Ntip(optimized_ml_tree_COI$tree)]

# Normalize distances for COI
range_COI <- range(distances_COI)
normalized_COI <- (distances_COI - range_COI[1]) / (range_COI[2] - range_COI[1])
adjusted_COI <- 0.3 + 0.7 * normalized_COI 
tip_colors_COI <- colorRampPalette(c("#103F7C", "#689ACD"))(100)[as.numeric(cut(adjusted_COI, breaks = 100))]

# Calculate distances from the root for RAG1
distances_RAG1 <- node.depth.edgelength(optimized_ml_tree_RAG1$tree)

# Calculate distances for only tip nodes (species)
distances_RAG1 <- distances_RAG1[1:Ntip(optimized_ml_tree_RAG1$tree)]

# Normalize distances for RAG1
range_RAG1 <- range(distances_RAG1)
normalized_RAG1 <- (distances_RAG1 - range_RAG1[1]) / (range_RAG1[2] - range_RAG1[1])
adjusted_RAG1 <- 0.3 + 0.7 * normalized_RAG1 
tip_colors_RAG1 <- colorRampPalette(c("#5B2323", "#EC5F67"))(100)[as.numeric(cut(adjusted_RAG1, breaks = 100))]

# The different shade ranges of colors in these trees visually represent the evolutionary distance of each species from the root, with darker shades indicating closer proximity to the root and lighter shades reflecting greater evolutionary divergence.

# Plot COI tree
plot(
  optimized_ml_tree_COI$tree,
  main = "Figure 3a: Maximum Likelihood Tree (COI)",
  show.tip.label = FALSE
)
tiplabels(
  text = optimized_ml_tree_COI$tree$tip.label,
  col = tip_colors_COI,
  frame = "none",
  adj = 0.5,
  cex = 0.8
)

# Plot RAG1 tree
plot(
  optimized_ml_tree_RAG1$tree,
  main = "Figure 3b: Maximum Likelihood Tree (RAG1)",
  show.tip.label = FALSE
)
tiplabels(
  text = optimized_ml_tree_RAG1$tree$tip.label,
  col = tip_colors_RAG1,
  frame = "none",
  adj = 0.5,
  cex = 0.8
)


# The Maximum Likelihood (ML) tree based on the mitochondrial COI gene (Figure 3a) reveals the evolutionary relationships among Mustelidae species with a focus on more recent divergence events. Species like Neovison vison and Galictis vittata cluster closely, reflecting shared evolutionary histories, while others such as Pteronura brasiliensis and Gulo gulo display longer branches, suggesting earlier divergence or faster mutation rates. The variability in branch lengths and clustering patterns shows the rapid evolutionary rate of the COI gene, which was the expected outcome. These genes often evolve faster due to higher mutation rates and maternal inheritance, making them well-suited for studying recent evolutionary events but potentially less reliable for deeper phylogenies (Árnadóttir et al., 2024)

# The RAG1-based ML tree (Figure 3b), in contrast, depicts a more conserved view of Mustelidae evolution, with clusters like Neovison vison and Martes foina still evident but embedded within a more uniform structure. The shorter and more consistent branch lengths across the tree reflect the slower evolutionary rate of nuclear genes like RAG1, which evolve more uniformly due to their biparental inheritance and role in maintaining genetic stability (Árnadóttir et al., 2024). Species such as Gulo gulo and Galictis vittata appear as distinct outliers here as well.

# When compared, the COI and RAG1 trees highlight the expected differences between mitochondrial and nuclear markers in evolutionary studies. The faster mutation rate of the COI gene results in more variable branch lengths and clustering, capturing recent divergence but potentially introducing noise in deeper phylogenetic relationships (Árnadóttir et al., 2024). In contrast, the RAG1 tree provides a clearer and more stable view of deeper evolutionary splits due to its slower mutation rate and conserved nature (Árnadóttir et al., 2024).


#-----------------------------------------------------------------------------#
## 8. Examine Topological & Phylogenetic Congruence of ML Trees: ======
#-----------------------------------------------------------------------------#

# Examine Topological Congruence through Visualization - Cophylogenetic Plot:

# Set a random seed for reproducibility
set.seed(123)  

# Create a cophylogenetic plot (side-by-side comparison)
cophylo_obj <- cophylo(
  tr1 = optimized_ml_tree_COI$tree, # Tree 1: COI
  tr2 = optimized_ml_tree_RAG1$tree # Tree 2: RAG1
)

# Adjust margins to prevent cutoff
par(mar = c(1, 1, 1, 1))  # Adjust bottom, left, top, right margins

# Create and plot the cophylogenetic tree
plot(
  cophylo_obj,
  link.type = "curved",
  link.lwd = 2,
  link.col = "#927bd3",
  fsize = 0.7,
  main = "Figure 4: Cophylogenetic Comparison: COI vs. RAG1"
)

# This plot is a cophylogenetic plot, which visually compares the topologies of two phylogenetic trees—here, the COI-based tree (on the left) and the RAG1-based tree (on the right)—to highlight their congruence/incongruence (Charleston & Perkins, 2006). The dashed lines connecting species across the two trees represent how the same taxa are placed in each tree (Charleston & Perkins, 2006).

# This cophylogenetic plot highlights the incongruence between the COI (left) and RAG1 (right) phylogenetic trees, as evidenced by the many crossed dashed lines, reflecting discrepancies in species placements due to differences in evolutionary signals captured by the mitochondrial and nuclear markers. While some species, such as Gulo gulo, exhibit congruence with consistent placements across both trees, others, like Neovison vison, Martes foina, and Pteronura brasiliensis, show notable shifts in their positions, revealing distinct evolutionary histories. These differences align with the expected characteristics of the two genes: the mitochondrial COI gene evolves faster and captures recent divergences, whereas the nuclear RAG1 gene is more conserved and better represents deeper evolutionary splits. Overall, the plot underscores the complementary nature of mitochondrial and nuclear data, emphasizing the value of integrating multiple markers to construct a more comprehensive and robust understanding of phylogenetic relationships.


#-----------------------------------------------------------------------------#

# Calculate Metrics of Phylogenetic Congruence:

# The Robinson-Foulds (RF) distance is a metric to quantify topological differences between two phylogenetic trees (Llabrés et al., 2021). It compares the splits (bipartitions) in the trees and calculates the number of splits unique to each tree (Llabrés et al., 2021). 
      # A lower RF distance indicates greater congruence (Smith, n.d.):
            # 0: Perfect congruence (identical topologies).
            # 1: Complete incongruence (no shared splits).


# Compute Robinson-Foulds distance
rf_distance <- RF.dist(
  tree1 = optimized_ml_tree_COI$tree,
  tree2 = optimized_ml_tree_RAG1$tree,
  normalize = TRUE # Optional: Normalize the distance by the maximum possible RF distance
)
cat("Normalized Robinson-Foulds Distance:", rf_distance, "\n")

# The Normalized Robinson-Foulds (RF) Distance of 1 indicates that the two trees are completely incongruent in their topologies. Specifically, it means there is no overlap in the splits of the COI and RAG1 trees (Llabrés et al., 2021). This suggests that the evolutionary relationships inferred from the mitochondrial COI gene and the nuclear RAG1 gene differ entirely.


