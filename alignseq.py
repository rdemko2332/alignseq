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
    
    return parser
    
    
class Settings:
    def __init__(self, infile, outfile, program):
        self.infile = infile
        # before proceeding to define more stuff, make sure this definition is ok!
        # or, define a bunch, and then check a bunch.
        
        
        self.outfile = outfile
        self.program = program

    #def check_infile():
        # Make sure the input file exists: os.path.exists(self.infile)
    
    
    
#class Sequences:
# the sequence data, pre and post alignment (dictionaries encouraged)
# functions to do translation and backtranslation
# function to write final alignment to file  - very flexible with arguments!


#class Aligner:
# functions to perform the alignment
# functions to assert executables are there
#    def check_mafft(): SPIELMAN REMAINS UNSURE WHERE THIS SHOULD GO, you'll find the way by coding and something will make sense eventually
#        # This function should define the mafft path (shutil.which) and do the assert
    




def main():
    
    # Call a function to set up arguments and return an INSTANCE of the Settings class. Aka, a Settings object.
    args = collect_input_arguments() 
    #my_settings = Settings(args)

    # Call a function to set up an INSTANCE (or just define an instance here) of a Sequences class. Aka, a Sequences object.

      
main()



