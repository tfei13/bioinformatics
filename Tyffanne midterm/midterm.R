#1.Import and align your DNA sequences ####
#Set the working directory 
setwd("/Users/kylecameron/Documents/GitHub/bioinformatics/Tyffanne midterm")

# Load in packages that I might need
library(seqinr)
library(msa)
library(Biostrings)
library(UniprotR)
library(protti)


# download the fasta files under bioinformatics
mySequences <- Biostrings::readDNAStringSet("midtermfasta")

# turn file into an alignment
midtermAlignment <- msa::msa(mySequences)

#show complete alignment
print(midtermAlignment, show="complete")

#2.Check to see how different your samples are from one another. Are any of them different from the rest? If so, what kinds of mutations do you observe in this individual (or individuals)? ####
# next i'm going to add color to my alignment to better visualize any differences between them
# Create a DNA String Set from the alignment
seqhs<-DNAMultipleAlignment(midtermAlignment)
print(seqhs)

#at this point you can see the differences in color. clearly Homo sapiens 6 and 10 are different from the rest especially 6.
#it has both delietion and substitutions.

#3. You suspect that an individual (or individuals) in this population might have some mutations in this gene, but you don’t know what this gene might be. Compare your sequences to a database to figure out what the gene is.####
#I performed a BLAST search on the GenBank website using the initial FASTA sequences file. 
#The top match for my gene is identified as "Homo sapiens hbb gene for beta globin, partial cds," with the accession number GenBank LC121775.

#4. Find the individual that is the most different from the rest of the individuals in your dataset. Translate that sequence to protein. Write it to a fasta file.####
#I'm confident from looking at the colorful alignment the most different one is homo_sapiens_6
#I'm also going to compute a distance matrix to confirm.
#next step is to translate it into a protein.
# Define the DNA sequence
# homosapiensall <- readDNAStringSet("midtermfasta") # already read in as mySequences
# allhomosapiensaligned <- msa(homosapiensall)
# allhomosapiensaligned
homosapiensAl <- msaConvert(midtermAlignment, type = "seqinr::alignment")
homosapiensAl
DistanceMatrix <- dist.alignment(homosapiensAl, "identity")
DistanceMatrix

# I assume this string below is for individual #6? If so, extract it from the original squence file: 
dna_string <- mySequences$Homo_sapiens_6
# dna_string <- DNAString("AATCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGGGCTGGGCATGAAAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCACCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGGTGAGTCTATGGGACCCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGG")

#translate in AA sequence
homo6AA <- Biostrings::translate(dna_string)
print(homo6AA)

# Define the sequence header
header <- ">Homo_sapiens_6"

# Define the amino acid sequence
sequence <- as.character(homo6AA)

# Write the sequence header and sequence to a FASTA file
writeLines(c(header, sequence), "homo_sapiens_6.fasta")

#5. Use a database to figure out what your protein matches to. Click on the record for the best match. What is the accession number of this entry? ####
#I used UniProt's BLAST tool to identify my protein's match. It's likely "Hemoglobin subunit beta, HBB, Homo sapiens (human)."
# Its UniProt accession number is A0A0J9YWK4.

#6. Either using R or by searching in the database, what disease(s) is this gene associated with? Does this person have the disease?####
#through research I found out HBB (hemoglobin subunit beta) is associated with several diseases in humans, including:
# Beta thalassemia (HBB/LCRB)
# Beta-zero thalassemia
# Beta-plus thalassemia
# Dominant-beta thalassemia
# Sickle cell disease
# Methemoglobinemia, beta-globin type
# Hemoglobin C disease
# I can't say with confidance what disease this person has, but my educated guess would be Beta thalassemia

#7. What is the 3-dimensional structure of this protein? You can include a screenshot or download of a photo of this structure in your GitHub repository.
 

