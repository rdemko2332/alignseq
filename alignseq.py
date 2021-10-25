from Bio import SeqIO
import pprint
import argparse
import shutil
import os

def collect_input_arguments():
    """
    This stuff is documenting the function.
    """

    parser = argparse.ArgumentParser(prog= 'Alignseq', description='Align Codon Sequences', usage='%(prog)s [options]',epilog="And that's how you make an Alignment!")
    parser.add_argument('--inf', metavar='Infile', action='store', help='A input file of codons')
    parser.add_argument('--outf', metavar='Outfile', action='store', help='A Output file (desired path) of codon Alignment')
    parser.add_argument('--prog', metavar='Program', action='store', help='Desired program to Align Sequences', default='mafft')
    #parser.add_argument('--act', metavar='Action', action='store', help='Action you want to do')
    
    return parser.parse_args()

class Settings:
    def __init__(self, infile, outfile, program):
        self.infile = infile
        self.outfile = outfile
        self.program = program
        
        check_infile_exists()
        check_program_exists()
    
    
    def check_infile_exists(self):
        assert os.path.exists(self.infile), "Path to infile does is incorrect."

    def check_program_exists(self):
        prog_path = shutil.which(self.program)
        assert prog_path is not None, "Cannot find alignment software."


        
class Sequences:
    def __init__(self):
        # coming up eventually, not in init though!
        # self.aligned = {}
        # self.aligned_translated = {}
        
    def translate_data(self, infile):
        with open(infile, "r") as file_handle:
            records = list(SeqIO.parse(file_handle, "fasta"))
            self.unaligned_translated = {}
            for record in records:
                self.unaligned[record.id] = record.seq
                self.unaligned_translated[record.id] = record.seq.translate()
            
            
            
        
            
        # functions to do translation and backtranslation
# function to write final alignment to file  - very flexible with arguments!

#class Aligner:
# functions to perform the alignment
# functions to assert executables are there
#    def check_mafft(): SPIELMAN REMAINS UNSURE WHERE THIS SHOULD GO, you'll find the way by coding and something will make sense eventually
#        # This function should define the mafft path (shutil.which) and do the assert
    
#/home/demkor62/Downloads/nuc.fasta
      
def main():
    args = collect_input_arguments() 
    my_settings = Settings(args.inf, args.outf, args.prog)
    
    my_sequences = Sequences()
    my_sequences.translate_data(my_settings.infile)
    
    
    #os.system(mafft_path + " --quiet --preservecase " + my_settings.infile + " > " + my_settings.outfile)

main()



