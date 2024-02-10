library(Biostrings)

install.packages("seqinr")
library(seqinr)
library(msa)
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")
a
library(msa)
system.file("tex", "texshade.sty", package="msa")
?readDNAStringSet
mySequenceFile <- system.file("bioinformatics", "sequence-1fb", package="msa")
mySequenceFile <- system.file("examples", "exampleAA.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile)
mySequences
library(Biostrings)
?readDNAStringSet
dna_sequences <- readDNAStringSet("sequence-1fb")
dna_sequences <- readDNAStringSet("sequence-2fb")
dna_sequences <- readDNAStringSet("sequence-3fb")
dna_sequences <- readDNAStringSet("sequence-4fb")
dna_sequences <- readDNAStringSet("sequence-5fb")


sequence_1 <- readDNAStringSet("sequence-1fb")
sequence_2 <- readDNAStringSet("sequence-2fb")
sequence_3 <- readDNAStringSet("sequence-3fb")
sequence_4 <- readDNAStringSet("sequence-4fb")
sequence_5 <- readDNAStringSet("sequence-5fb")
combined_sequences <- DNAStringSet(c(sequence_1, sequence_2, sequence_3, sequence_4, sequence_5))


library(msa)


clustalw_alignment <- msa(combined_sequences, method = "ClustalW", type = "dna")

print(clustalw_alignment)

print(clustalw_alignment, show = "complete")
num_gaps <- nGaps(clustalw_alignment)
library(Biostrings)
alignment_matrix <- as.matrix(clustalw_alignment)
num_gaps <- colSums(alignment_matrix == "-")
total_gaps <- sum(num_gaps)
print(num_gaps)
print(total_gaps)

library(Biostrings)


alignment_length <- ncol(as.matrix(clustalw_alignment))

print(alignment_length)

library(Biostrings)

letter_freq <- alphabetFrequency(clustalw_alignment)


GC_content <- (letter_freq["G"] + letter_freq["C"]) / sum(letter_freq) * 100


print(GC_content)
print(letter_freq)


library(seqinr)


seqinr_alignment <- msaConvert(clustalw_alignment, type = "seqinr")


distance_matrix <- dist.alignment(seqinr_alignment)


print(distance_matrix)


print(ls())


myMaskedAlignment <- clustalw_alignment
colM <- IRanges(start = 1, end = 100)  
colmask(myMaskedAlignment) <- colM  
print(myMaskedAlignment)


# Convert the MsaDNAMultipleAlignment object to an alignment object
aligned_sequences <- msaConvert(myMaskedAlignment, type = "seqinr::alignment")

# Compute the distance matrix from the converted alignment object
distance_matrix <- dist.alignment(aligned_sequences)

# Find the indices of the pair with the smallest distance
min_distance_index <- which(distance_matrix == min(distance_matrix), arr.ind = TRUE)

# Get the names of the sequences corresponding to the indices
closest_samples <- rownames(myMaskedAlignment)[min_distance_index[, 1]]

# Print the names of the two closest samples
print(closest_samples)














