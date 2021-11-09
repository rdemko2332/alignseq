from Bio import SeqIO
import shutil
import os

# We know the unaligned nucleotides, unaligned aa, and ALIGNED nucleotides
filename = "/home/demkor62/Desktop/alignseq/out.fasta"
format = "fasta"
def read_sequences_from_file(infile, format = "fasta"):
        """
            Function to read in a file of sequences as a list.
        """
        with open(infile, "r") as file_handle:
            records = list(SeqIO.parse(file_handle, format))
        return(records)
list_of_aligned_aa = read_sequences_from_file(filename, format)
def convert_sequence_list_to_dictionary(sequence_list):
        """
            Convert a list of sequences into a dictionary of sequences, where key is the sequence ID.
        """
        dictionary = {}
        for record in sequence_list:
            dictionary[record.id] = record.seq
        return(dictionary)
def read_in_aligned_data(filename, format = "fasta"):
        list_of_aligned_aa = read_sequences_from_file(filename, format)
        aligned_translated = convert_sequence_list_to_dictionary(list_of_aligned_aa)
        return(aligned_translated)
data = read_in_aligned_data("/home/demkor62/Desktop/alignseq/practice_backtranslation/aligned_nuc.fasta", format="fasta")
print(data)
program_path = shutil.which("mafft")
#os.system(program_path + " " + "--quiet --preservecase" + " /home/demkor62/Desktop/alignseq/practice_backtranslation/unaligned_nuc.fasta > /home/demkor62/Desktop/temp.fasta")
# PROBLEM! When the unaligned nucleotide fasta is placed in the aligner, it is not aligning it correctly. Check temp.fasta, aligned_nuc.fasta and unaligned_nuc.fasta
# If the aligner was working as intended and I recieved aligned_nuc.fasta, this is what I would do
table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        '---':'-'
    }
aa_aligned = {}
for id in data:
    protein_sequence = []
    for nuc in range(0,len(data[id]),3):
        codon = str(data[id][nuc:nuc+3])
        protein = table[codon]
        protein_sequence.append(protein)
    aa_aligned[str(id)] = protein_sequence
print(aa_aligned)

aa_aligned = {}
for id in data:
    protein_sequence = []
    for nuc in range(0,len(data[id]),3):
        codon = str(data[id][nuc:nuc+3])
        protein = table[codon]
        protein_sequence.append(protein)
    aa_aligned[str(id)] = protein_sequence
print(aa_aligned)

        
#def preform_backtranslation():







# We need to figure out ***aa_aligned**, and you can refer to the file to see what you should get.

# Use this script to get the hang of the for-loop you need to write to get to aa_aligned. Then, we can integrate into alignseq.py