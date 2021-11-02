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
        Creating Settings class. This will take it's arguments from the output of collect_input_arguments()
    """
    def __init__(self, infile, outfile, program, arguments, outtranslated):
        self.infile = infile
        self.outfile = outfile
        self.program = program
        self.arguments = arguments
        self.outtranslated = outtranslated
    
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
        Creating Sequences class. Will also take arguments from the output of collect_input_arguments()
    """
    def __init__(self, infile, outfile, outtranslated):
        self.infile = infile
        self.outfile = outfile
        self.outtranslated = outtranslated
        
    def translate_data(self, infile):
        """
            Translating data from nucleotides to amino acids. Creates and fills two dictionaries, self.unaligned and self.aligned
        """
        with open(infile, "r") as file_handle:
            records = list(SeqIO.parse(file_handle, "fasta"))
            self.unaligned = {}
            self.unaligned_translated = {}
            for record in records:
                self.unaligned[record.id] = record.seq
                self.unaligned_translated[record.id] = record.seq.translate()
    
    def write_translated_seqs_to_file(self):
        """
            Creates the file for aligned output as specified in collect_input_arguments() and appends the translated data to it
        """
        os.system("touch " + self.outtranslated)
        file = open(self.outtranslated, "a")  # append mode
        file.write(str(self.unaligned_translated))
        file.close()

    def read_in_aligned_data(self, outfile):
        '''
            Reads in the alignment data from the outfile specified in collect_input_arguments() that was created after the aligner instance was called
        '''
        self.aligned = {} 
        records = list(SeqIO.parse(outfile, "fasta"))
        for record in records:
            self.aligned[record.id] = record.seq

    def __call__(self):
        '''
            Runs both functions when an instance is called as a function.
        '''
        self.translate_data(self.infile)
        self.write_translated_seqs_to_file()


class Aligner:
    """
        Creating the Aligned class. Will also take arguments from the arguments specified from collect_input_arguments(). Generates what will be called in os.system to run the alignment.
    """
    def __init__(self, infile, outfile, program, arguments):
        self.infile = infile
        self.outfile = outfile
        self.program = program
        self.arguments = arguments
        '''
            Defines the path for the specified alignment program
        '''
        program_path = shutil.which(self.program)
        self.command = (program_path + " " + self.arguments + " " + self.infile + ">" + self.outfile)
    
    def __call__(self):
        os.system(self.command)


def main():
    args = collect_input_arguments() 
    my_settings = Settings(args.inf, args.outf, args.prog, args.args, args.outtranslated)
    my_sequences = Sequences(my_settings.infile, my_settings.outfile, my_settings.outtranslated)
    my_sequences()
    my_aligner = Aligner(my_settings.infile, my_settings.outfile, my_settings.program, my_settings.arguments)
    my_aligner()
    my_sequences.read_in_aligned_data(my_settings.outfile)   

main()



