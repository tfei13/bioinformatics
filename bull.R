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


GC_content <- (letter_freq[,"G"] + letter_freq[,"C"]) / sum(letter_freq[,1:4]) * 100


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



aligned_sequences <- msaConvert(myMaskedAlignment, type = "seqinr::alignment")


distance_matrix <- dist.alignment(aligned_sequences)


min_distance_index <- which(distance_matrix == min(distance_matrix), arr.ind = TRUE)


closest_samples <- rownames(myMaskedAlignment)[min_distance_index[, 1]]

# Print the names of the two closest samples
print(closest_samples)




#answeers

library(seqinr)



install.packages("phangorn")
install.packages("tidyr")
install.packages("dplyr")

install.packages("phangorn") - the latest development version remotes::install_github("KlausVigo/phangorn")

# load in all of the libraries that you might need
# this should always be at the start of your script
library(msa)
library(Biostrings)
library(seqinr)
library(phangorn)
library(tidyr)
library(dplyr)

# set the working directory to the folder containing all of your scripts and data
# filepaths and files should always be in quotes. Variables in R should not.
setwd("/Users/ojohnson/Documents/GitHub/Bioinformatics_Spring2024")

# read in albatross Cytochrome B sequences 
# assign each one to a variable
# Note that these fasta files are contained in a folder called 'Diomedea_exulans'
seq_1 <- readDNAStringSet("sequence-1fb")
seq_2 <- readDNAStringSet("sequence-2fb")
seq_3 <- readDNAStringSet("sequence-3fb")
seq_4 <- readDNAStringSet("sequence-4fb")
seq_5 <- readDNAStringSet("sequence-5fb")


# combine samples into a single variable using the combine ('c') function
seqs <- c(seq_1, seq_2, seq_3, seq_4, seq_5)

# The default GenBank names are very very long, so lets
# rename the samples to something shorter and more legible
# we do this by assigning a list of characters (using the same 'c' function)
# to the 'names' of the combined seqs variable
# check what these names are by first running just the names() function
names(seqs)
names(seqs) <- c("sanfordi_U48946.1", "chionoptera_AF076048.1", "epomophora_AF076049.1",
                 "gibsoni_AF076050.1", "antipodensis_MH330008.1", "antipodensis_MH330010.1",
                 "exulans_U48947.1", "amsterdamensis_U48948.1")

# run the MSA! Assign it to a new variable
bulldog <- msa(seqs)

# check the alignment length, two different ways
nchar(bulldog)
print(bulldog, show="complete") 
# here, you can also calculate the number of gaps by hand, or use the next step

# Calculate the GC content. First, calculate the frequency of each nucleotide
alFreq <- alphabetFrequency(bulldog)
alFreq # here, it also gives you the number of dashes (-) in the alignment, which is the number of gaps

# now pull out the total number of G's and C's
# here, the 'sum' function takes the sum, as you might expect
# the square brackets are for accessing rows and columns of a matrix
# values before the comma access rows, those after the comma access columns
# we want the columns
GC <- sum(alFreq[,"C"]) + sum(alFreq[,"G"]) 
AT <- sum(alFreq[,"A"]) + sum(alFreq[,"T"]) 
# and calculate the percentage that are G or C (out of the total nucleotides)
GC / (GC + AT )


# calculate the GC content a different way, using the 'GC' function in the seqinr package
# we can only run this on one sample at a time, so let's read in the first sample using the 
# read.fasta() function in the seqinr package
# note that because there are multiple 'read.fasta()' functions, we need to specify that
# want to use the one in the seqinr package using the double colon
seq_1.seqinr <- seqinr::read.fasta("Diomedea_exulans/sequence_1.fasta")
# then select the sequence data from the variable. Use the '$' to access it 
seqinr::GC(seq_1.seqinr$U48946.1)

# calculate the identity matrix
# first, convert the alignment to the seqinr format using msaConvert
# because the dist.alignment() function is part of the seqinr package
seq_1 <- msaConvert(bulldog, type="seqinr::alignment")
d <- dist.alignment(bulldog, "identity")
d
# this is a fancy way to compare only my 'epomophora' sample to the other samples in the matrix
# and convert the numbers to a percentage 
100 - (round(as.matrix(d)[, "epomophora_AF076049.1", drop=FALSE], digits = 2) * 100)

# translate one sample to an amino acid sequence
# we again need to specify which package to use because the 'translate()' function 
# exists in both the Biostrings and seqinr packages
seq_1_AA <- Biostrings::translate(seq_1)
print(seq_1_AA)


library(Biostrings)

# Replace 'sequence.fasta' with the path to your file
dna_seq <- readDNAStringSet("sequence-1fb")

# Assuming dna_string is your DNA sequence
substring <- substr(dna_string, start = 1, stop = 15515)
translated_seq <- translate(Biostrings::DNAString(substring))


# Print the translated amino acid sequence
print(translated_seq)


# write the alignment to a fasta file (harder than I expected to figure this out) 
# there is a write function in the phangorn package, but not one that I could find in seqinr or Biostrings
# Biostrings has a write function, but not for fasta-formatted files
sequencefb_phyDat <- msaConvert(sequence-1fb, type="phangorn::phyDat")
write.phyDat(sequencefb_phyDat, "frenchbulldog/frenchie1alignment.fasta", format = "fasta")




library(seqinr)

# Assuming 'translated_seq' is your AAString or AAStringSet object containing amino acid sequences

# Convert AAStringSet to character vector
aa_seqs <- as.character(translated_seq)
aa_seqs
# Define the file name for the FASTA file
fasta_file <- "translated_sequences.fasta"

# Write amino acid sequences to a FASTA file
write.fasta(sequences = aa_seqs, names = paste0("Sequence_", seq_along(aa_seqs)), file.out = fasta_file)

# Read the FASTA file back into R
fasta_data <- read.fasta(fasta_file)

# Access sequences and names
sequences <- fasta_data$seq
names <- fasta_data$name

# Print the first sequence
print(sequences[[1]])

# Print the name of the first sequence
print(names[1])

#Creating FASTA File from aa seq ####
output_file <- "aa_sequence.fasta"
writeXStringsSet(aa_seqs, file = output_file,
              format = "fasta", width = 60)

write.fasta(aa_seqs, "")



