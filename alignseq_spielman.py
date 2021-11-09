from Bio import SeqIO
import argparse
import shutil
import os

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
    
    return parser.parse_args()


class Settings:
    """
        Creating Settings class. This will take its arguments from the output of collect_input_arguments()
    """
    def __init__(self, infile, outfile, program, arguments, outtranslated):
        self.infile = infile
        self.outfile = outfile
        self.program = program
        self.arguments = arguments
        self.outtranslated = outtranslated
        
        # SPIELMAN: Need to perform assertions by calling those methods
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
        Creating Sequences class. SPIELMAN: NOPE, DELETE THIS:  Will also take arguments from the output of collect_input_arguments()
    """
    def __init__(self):
        return None 

    def prepare_for_alignment(self, infile, temporary_file):
        self.translate_sequences(infile)
        self.write_translated_seqs_to_file(temporary_file)                

    def translate_sequences(self, infile): # SPIELMAN: See, we just pass infile as an argument. It's not a self.infile
        """
            Translating sequences from nucleotides to amino acids. Creates and fills two dictionaries, self.unaligned and self.unaligned_translated
        """
        records = self.read_sequences_from_file(infile)
        self.unaligned = self.convert_sequence_list_to_dictionary(records)
        self.unaligned_translated = {}
        for record_id in self.unaligned:
            self.unaligned_translated[record_id] = self.unaligned[record_id].translate()

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
        
    # def backtranslate_sequences(self, temporary_file2):
    #     read_in_aligned_data  <- will read in out.fasta
    
    # def save_sequences_to_file(self, aafile, nucfile)


class Aligner:
    """
        Creating the Aligner class. Will also take arguments from the arguments specified from collect_input_arguments(). Generates what will be called in os.system to run the alignment.
    """
    def __init__(self, settings):
        '''
            Sets up Aligner instance.
        '''
        
        self.program_path = shutil.which(settings.program)
        self.command = (self.program_path + " " + settings.arguments + " in.fasta > out.fasta")
    
    def __call__(self):
        os.system(self.command)


def main():
    # Collect input
    args = collect_input_arguments() 
    my_settings = Settings(args.inf, args.outf, args.prog, args.args, args.outtranslated)
    
    # Initialize sequences and prepare for alignment
    my_sequences = Sequences() 
    my_sequences.prepare_for_alignment(my_settings.infile, "in.fasta") # writes to "out.fasta"
    
    # Perform alignment
    my_aligner = Aligner(my_settings)
    my_aligner()
    
    # Backtranslate
    my_sequences.read_in_aligned_data("out.fasta")
    # my_sequences.perform_backtranslation()
    
    # Export
    # my_sequences.save_to_file()

main()



