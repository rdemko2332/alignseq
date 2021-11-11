from Bio import SeqIO
from Bio import AlignIO
import argparse
import shutil
import os
from tempfile import mkstemp


def collect_input_arguments():
    """
        Collecting input arguments at the command line for use later.
    """

    parser = argparse.ArgumentParser(prog= 'Alignseq', description='Align Codon Sequences', usage='%(prog)s [options]',epilog="And that's how you make an Alignment!")
    parser.add_argument('--inf', metavar='Infile', action='store', help='A input file of codons')
    parser.add_argument('--outf', metavar='Outfile', action='store', help='An Output file (desired path) of codon Alignment')
    parser.add_argument('--prog', metavar='Program', action='store', help='Desired program to Align Sequences', default='mafft')
    parser.add_argument('--args', metavar='Arguments', action='store', help='Arguments for the program you are  running', default= "--quiet --preservecase")
    parser.add_argument('--outtranslated', metavar='Outfile for Translated Data', action='store', help='An Output file (desired path) for translated data')
    parser.add_argument('--outtransaligned', metavar='Outfile for Translated and Aligned Data', action='store', help='An Output file (desired path) for translated and aligned data')
    parser.add_argument('--outformat', metavar='Output Format', action='store', help='An Output Format', default = "fasta")

    return parser.parse_args()


class Settings:
    """
        Creating Settings class. This will take its arguments from the output of collect_input_arguments()
    """
    def __init__(self, infile, outfile, program, arguments, outtranslated, outtransaligned, outformat):
        self.infile = infile
        self.outfile = outfile
        self.program = program
        self.arguments = arguments
        self.outtranslated = outtranslated  # None
        self.outtransaligned = outtransaligned
        self.outformat = outformat
        
        self.check_infile_exists()
        self.check_program_exists()
    
    def check_infile_exists(self):
        """
            Assert that the input file containing nucleotide sequences exists.
        """
        assert os.path.exists(self.infile), "Path to infile does is incorrect."

    def check_program_exists(self):
        """
            Assert that the program being used to align the sequences exists
        """
        prog_path = shutil.which(self.program)
        assert prog_path is not None, "Cannot find alignment software."


class Sequences:
    """
        Creating Sequences class. 
    """
    def __init__(self):
        return None 

    def prepare_for_alignment(self, infile, temporary_file):
        self.translate_sequences(infile)
        self.write_translated_seqs_to_file(temporary_file)                

    def translate_sequences(self, infile): 
        """
            Translating sequences from nucleotides to amino acids. Creates and fills two dictionaries, self.unaligned and self.unaligned_translated
        """
        records = self.read_sequences_from_file(infile)
        self.unaligned = self.convert_sequence_list_to_dictionary(records)
        #print("self.unaligned -> ", self.unaligned)
        self.unaligned_translated = {}
        for record_id in self.unaligned:
            self.unaligned_translated[record_id] = self.unaligned[record_id].translate()
        #print("self.unaligned_translated ->", self.unaligned_translated)

    def write_translated_seqs_to_file(self, temporary_file):
        """
            Creates the file for aligned output and appends the translated data to it
        """
        with open(temporary_file, "w") as f:
            for record in self.unaligned_translated:
                f.write(">" + str(record) + "\n" + str(self.unaligned_translated[record]) + "\n")
    
    def read_sequences_from_file(self, infile, format = "fasta"):
        """
            Function to read in a file of sequences as a list.
        """
        with open(infile, "r") as file_handle:
            records = list(SeqIO.parse(file_handle, format))
        return(records)
    
    def convert_sequence_list_to_dictionary(self, sequence_list):
        """
            Convert a list of sequences into a dictionary of sequences, where key is the sequence ID.
        """
        dictionary = {}
        for record in sequence_list:
            dictionary[record.id] = record.seq
        return(dictionary)
    
    def read_in_aligned_data(self, filename, format = "fasta"):
        list_of_aligned_aa = self.read_sequences_from_file(filename, format)
        self.aligned_translated = self.convert_sequence_list_to_dictionary(list_of_aligned_aa)
        
    def backtranslate_sequences(self):
        final = {}
        for id in self.aligned_translated:
            aligned_nuc_sequence = ""
            aligned_aa_sequence = self.aligned_translated[id]
            unaligned_nuc_sequence = self.unaligned[id]
            nuc_index = 0
            for aa in aligned_aa_sequence:
                if aa == "-":
                    new_codon = "---"
                else:
                    new_codon = unaligned_nuc_sequence[nuc_index:nuc_index+3]
                    nuc_index += 3
                aligned_nuc_sequence += new_codon
            final[id] = aligned_nuc_sequence
        self.aligned_nuc = final
    
    def save_sequences_to_file(self, outfile, outputformat): #aa_aligned, nucfile
        temp_file_handle, temp_file_path = mkstemp()
        with open(temp_file_path, "w") as f:
            for record in self.aligned_nuc:
                f.write(">" + str(record) + "\n" + str(self.aligned_nuc[record]) + "\n")
        AlignIO.convert(temp_file_path, "fasta", outfile, outputformat)
        os.close(temp_file_handle)    


class Aligner:
    """
        Creating the Aligner class. Will also take arguments from the arguments specified from collect_input_arguments(). Generates what will be called in os.system to run the alignment.
    """
    def __init__(self, settings):
        '''
            Sets up Aligner instance.
        '''
        self.program_path = shutil.which(settings.program)
        #print(settings.outtransaligned)
        #print(settings.outfile)
        self.command = (self.program_path + " " + settings.arguments + " " +settings.outtranslated + " > " + settings.outtransaligned)
        # " ".join([self.program_path, settings.arguments, settings.outtranslated, ">", settings.outtransaligned])
    
    
    
    
    
    def __call__(self):
        os.system(self.command)


def main():
    # Collect input
    args = collect_input_arguments() 
    my_settings = Settings(args.inf, args.outf, args.prog, args.args, args.outtranslated, args.outtransaligned, args.outformat)
    
    # Initialize sequences and prepare for alignment
    my_sequences = Sequences() 
    my_sequences.prepare_for_alignment(my_settings.infile, my_settings.outtranslated) # writes to "my_settings.outtranslated"
    
    #Perform alignment
    my_aligner = Aligner(my_settings)
    my_aligner()
    
    # Backtranslate
    my_sequences.read_in_aligned_data(my_settings.outtransaligned)
    #print("aligned_translated ->", my_sequences.aligned_translated)
    my_sequences.backtranslate_sequences()
    #print(my_sequences.aligned_nuc)
    
    # Export
    my_sequences.save_sequences_to_file(my_settings.outfile, my_settings.outformat)

main()

# Example Usage For Me
#python3 alignseq.py --inf example.fasta --outtranslated outtranslated.fasta --outtransaligned outtransaligned.fasta --outf new_alignment.fasta 


"""

In Settings() make a function: check_outfiles()

If the files are None (self.outtranslated is None) then define as follows...
if this is my infile:
infile.fasta

# the outtranslated will be:
infile.fasta.translated # outtranslated default

infile.fasta.translated_aligned # outtransaligned default

infile.fasta.aligned # outfile default


ALSO DO THESE THINGS:
+ when you run without args like python3 alignseq.py, **the help should come up** not an error
+ make a MafftAligner child class (and ClustalOmegaAligner????) to inherit from Aligner parent. All the parent should define is self.program_path. Rest is in child.
+ optimize write_sequences_to_file (needs argument not self)
+ docstrings everywhere
"""

































