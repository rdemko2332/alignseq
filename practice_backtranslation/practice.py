from Bio import Seq

# We know the unaligned nucleotides, unaligned aa, and ALIGNED nucleotides
nuc_unaligned = list(SeqIO.read("nuc_unaligned.fasta", "fasta"))
aa_unaligned = list(SeqIO.read("aa_unaligned.fasta"))
nuc_aligned = list(SeqIO.read("nuc_aligned.fasta", "fasta"))

# We need to figure out ***aa_aligned**, and you can refer to the file to see what you should get.

# Use this script to get the hang of the for-loop you need to write to get to aa_aligned. Then, we can integrate into alignseq.py