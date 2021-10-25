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
        # before proceeding to define more stuff, make sure this definition is ok!
        # or, define a bunch, and then check a bunch.
        self.outfile = outfile
        self.program = program


def make_settings_class(infile, outfile, program):
    my_settings = Settings(infile,outfile, program)
    return(my_settings)

def check_infile_exists(instance_infile):
    result = os.path.exists(instance_infile)
    return(result)
        
class Sequences:
    def __init__(self, instance_infile, instance_outfile):
        self.infile = instance_infile
        self.outfile = instance_outfile
        self.translated = self.translate_data()
    def translate_data(self):
        with open(self.infile, "r") as file_handle:
            records = list(SeqIO.parse(file_handle, "fasta"))
            dictionary = {}
            for record in records:
                dictionary[record.id] = record.seq.translate()
            return(dictionary)
        # functions to do translation and backtranslation
# function to write final alignment to file  - very flexible with arguments!

#class Aligner:
# functions to perform the alignment
# functions to assert executables are there
#    def check_mafft(): SPIELMAN REMAINS UNSURE WHERE THIS SHOULD GO, you'll find the way by coding and something will make sense eventually
#        # This function should define the mafft path (shutil.which) and do the assert
    
#/home/demkor62/Downloads/nuc.fasta
      
def main():
    # Call a function to set up arguments and return an INSTANCE of the Settings class. Aka, a Settings object.
    args = collect_input_arguments() 
    my_settings = make_settings_class(args.inf, args.outf, args.prog)
    if check_infile_exists(my_settings.infile) == True:
        my_sequences = Sequences(my_settings.infile, my_settings.outfile)
        print(my_sequences.translated)
        if my_settings.program == 'mafft':
            mafft_path = shutil.which("mafft")
            os.system(mafft_path + " --quiet --preservecase " + my_settings.infile + " > " + my_settings.outfile)
    else:
        print('Input File Does Not Exist')

    # Call a function to set up an INSTANCE (or just define an instance here) of a Sequences class. Aka, a Sequences object.
    
main()



