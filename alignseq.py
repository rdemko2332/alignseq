from Bio import SeqIO
from Bio import AlignIO
import argparse
import shutil
import os
import sys
from tempfile import mkstemp

stop_codons = ["TAG", "TAA", "TGA"]

def collect_input_arguments():
    """
        Collecting input arguments at the command line for use later.
    """

    parser = argparse.ArgumentParser(prog= 'Alignseq', description='Align Codon Sequences', prefix_chars='+', usage='%(prog)s [options]', epilog="And that's how you make an Alignment!")
    parser.add_argument('+inf', metavar='Infile', action='store', help='A input file of codons')
    parser.add_argument('+outf', metavar='Outfile', action='store', help='An Output file (desired path) of codon Alignment')
    parser.add_argument('+prog', metavar='Program', action='store', help='Desired program to Align Sequences', default='mafft')
    parser.add_argument('+args', metavar='Arguments', action='store', nargs='*', help='Arguments for the program you are running. MUST BE PROVIDED IN QUOTES!!', default= "--quiet --preservecase")
    parser.add_argument('+outtranslated', metavar='Outfile for Translated Data', action='store', help='An Output file (desired path) for translated data')
    parser.add_argument('+outtransaligned', metavar='Outfile for Translated and Aligned Data', action='store', help='An Output file (desired path) for translated and aligned data')
    parser.add_argument('+outformat', metavar='Output Format', action='store', help='An Output Format', default = "fasta")
    
    # Will print the help menu if no arguments are passed to alignseq.py.
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    # Returns arguments for use in classes.
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
        self.outtranslated = outtranslated  
        self.outtransaligned = outtransaligned
        self.outformat = outformat
        
        self.check_infile_exists()
        self.check_program_exists()
        self.check_outfiles()
        self.check_arguments()
    
    def check_infile_exists(self):
        """
            Assert that the input file containing nucleotide sequences exists.
        """
        assert os.path.exists(self.infile), "Path to infile does not exist."

    def check_program_exists(self):
        """
            Assert that the program being used to align the sequences exists
        """
        prog_path = shutil.which(self.program)
        assert prog_path is not None, "Cannot find alignment software."
    
    def check_outfiles(self):
        """
            Specifies outfile, outtranslated and outtransligned files if they are not specified in collect_input_arguments().
        """
        if self.outfile is None:
            self.outfile = str(self.infile + '.aligned')
        if self.outtranslated is None:
            self.outtranslated = str(self.infile + '.translated')
        if self.outtransaligned is None:
            self.outtransaligned = str(self.infile + '.translated_aligned')
    
    def check_arguments(self):
        """
            Redefines default arguments in clustalo is specified as the program in collect_input_arguments().
        """
        if self.program == "clustalo" and self.arguments == "--quiet --preservecase":
            self.arguments = "-v"


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
        self.unaligned_translated = {}
        for record_id in self.unaligned:
        
            # Check sequence: Is its length a multiple of three?  len(sequence) % 3 == 0 is the goal
            assert len(self.unaligned[record_id]) % 3 == 0, "Error: Provides sequences are not codons. Lengths are not multiples of 3."
                
            translated = self.unaligned[record_id].translate()
            
            # If there's a stop, purge it from everything.
            if translated.endswith("*"):
                translated = translated.rstrip("*")
                self.unaligned[record_id] = self.unaligned[record_id][:-3]

            self.unaligned_translated[record_id] = translated





    def write_translated_seqs_to_file(self, temporary_file):
        """
            Creates the file for aligned output and writes the translated data to it
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
        """
            Generates self.aligned_translated from reading in the output from the specified aligner (outtransaligned)
        """
        list_of_aligned_aa = self.read_sequences_from_file(filename, format)
        self.aligned_translated = self.convert_sequence_list_to_dictionary(list_of_aligned_aa)
        
    def backtranslate_sequences(self):
        """
            Generates an aligned nucleotide sequence by comparing self.aligned_translated and self.unaligned
        """
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
        """
            Writes self.aligned_nuc to tempfile, converts the tempfile to the prefered outfile path and format, and closes the tmepfile.
        """
        temp_file_handle, temp_file_path = mkstemp()
        with open(temp_file_path, "w") as f:
            for record in self.aligned_nuc:
                f.write(">" + str(record) + "\n" + str(self.aligned_nuc[record]) + "\n")
        AlignIO.convert(temp_file_path, "fasta", outfile, outputformat)
        os.close(temp_file_handle)    


class Aligner:
    """
        Creating the Aligner class. 
    """
    def __init__(self, settings):
        '''
            Sets up Aligner instance. Defines the program path for the program specified in collect_input_arguments().
        '''
        self.program_path = shutil.which(settings.program)
        

class Mafft_aligner(Aligner):
    """
        Creates the Mafft_aligner subclass. 
    """
    def __init__(self, settings):
        """
            Sets up Mafft_aligner instance. super() allows program_path to be inherited from Aligner class. Defines self.command.
        """
        super().__init__(settings)
        self.command =  " ".join([self.program_path, settings.arguments, settings.outtranslated, ">", settings.outtransaligned])
        # " ".join([self.program_path, settings.arguments, settings.outtranslated, ">", settings.outtransaligned])
    def __call__(self):
        """
            Makes the Mafft_aligner instance callable, runs self.command when it is run.
        """
        os.system(self.command)

class ClustalOmega_aligner(Aligner):
    """
        Creates the ClustalOmega_aligner subclass.
    """
    def __init__(self, settings):
        """
            Sets up ClustalOmega_aligner instance. super() allows program_path to be inherited from Aligner class. Defines self.command.
        """
        super().__init__(settings)
        self.command = (self.program_path + " -i " + settings.outtranslated + " -o " + settings.outtransaligned + " " + " ".join(settings.arguments))
        # " ".join([self.program_path, settings.arguments, settings.outtranslated, ">", settings.outtransaligned])
    def __call__(self):
        """
            Makes the ClustalOmega_aligner instance callable, runs self.command when it is run.
        """
        os.system(self.command)



def main():
    # Collect input
    args = collect_input_arguments() 
    my_settings = Settings(args.inf, args.outf, args.prog, args.args, args.outtranslated, args.outtransaligned, args.outformat)
    
    # Initialize sequences and prepare for alignment
    my_sequences = Sequences() 
    my_sequences.prepare_for_alignment(my_settings.infile, my_settings.outtranslated) # writes to "my_settings.outtranslated"
    
    #Perform alignment
    if my_settings.program == "mafft":
        my_aligner = Mafft_aligner(my_settings)
    else: 
        my_aligner = ClustalOmega_aligner(my_settings)
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
#Mafft
#python3 alignseq.py +inf example.fasta 
#ClustalOmega
# python3 alignseq.py +inf example.fasta +outf out.fasta +args -v +prog clustalo


