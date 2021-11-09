from Bio import SeqIO
import shutil
import os
from Bio import AlignIO

def convert_sequence_list_to_dictionary(sequence_list):
    """
        Convert a list of sequences into a dictionary of sequences, where key is the sequence ID.
    """
    dictionary = {}
    for record in sequence_list:
        dictionary[record.id] = record.seq
    return(dictionary)
    
def write_translated_seqs_to_file(sequences, temporary_file):
    """
        Creates the file for aligned output and appends the translated data to it
    """
    with open(temporary_file, "w") as f:
        for record in sequences:
            f.write(">" + str(record) + "\n" + str(sequences[record]) + "\n")
    

nuc_unaligned = convert_sequence_list_to_dictionary(list(SeqIO.parse("unaligned_nuc.fasta", "fasta"))) # a user provides
#aa_unaligned = convert_sequence_list_to_dictionary(list(SeqIO.parse("unaligned_aa.fasta", "fasta"))    # your code translates
aa_aligned = convert_sequence_list_to_dictionary(list(SeqIO.parse("aligned_aa.fasta", "fasta")))  # mafft makes this!!
#nuc_aligned = convert_sequence_list_to_dictionary(list(SeqIO.parse("aligned_nuc.fasta", "fasta")))   # YOU NEED TO MAKE THIS


#{'stephanie': Seq('SRIL-RK'), 'rich': Seq('SM-LFR-')}
#{'stephanie': Seq('AGT CGA ATT CTG AGA AAA'), 'rich': Seq('AGTATGCTGTTCAGG')}

final = {}
for id in aa_aligned:
    aligned_nuc_sequence = ""
    aligned_aa_sequence = aa_aligned[id]
    unaligned_nuc_sequence = nuc_unaligned[id]
    nuc_index = 0
    for aa in aligned_aa_sequence:
        #print(aa)
        
        if aa == "-":
            new_codon = "---"
        else:
            new_codon = unaligned_nuc_sequence[nuc_index:nuc_index+3]
            nuc_index += 3
        
        aligned_nuc_sequence += new_codon
    print(aligned_nuc_sequence)
    final[id] = aligned_nuc_sequence
print(final) # THIS IS IT YAY!!!
    
    
    
# Write to file in Spielman hack!
write_translated_seqs_to_file(final, "hack.fasta")
AlignIO.convert("hack.fasta", "fasta", "hack.phy", "phylip-relaxed")
# AlignIO.convert(input, input format, output, output format)
    
    
    
    
    
    
    
    
    
    
    

