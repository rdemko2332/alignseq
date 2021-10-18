import argparse
import shutil
import os

parser = argparse.ArgumentParser(prog= 'Alignseq', description='Align Nucleotide Sequences', usage='%(prog)s [options]',epilog="And that's how you make an Alignment!")
parser.add_argument('--inf', metavar='Infile', action='store', help='A input file of Nucleotides')
parser.add_argument('--ouf', metavar='Outfile', action='store', help='A Output file (desired path) of Nucleotide Alignment')
parser.add_argument('--prog', metavar='Program', action='store', help='Desired program to Align Sequences', default='mafft')
#parser.add_argument('--act', metavar='Action', action='store', help='Action you want to do')

args = parser.parse_args()
infi = str(args.inf)
outfi = str(args.ouf)
pro = str(args.prog)

class Vars:
  def __init__(self, infile, outfile, program):
      self.infile = infile
      self.outfile = outfile
      self.program = program

test = Vars(infi, outfi, pro)

def main():
    if test.program == "mafft":
        program_path = shutil.which("mafft")
        assert(program_path is not None) # make sure mafft is there!!
        variables = " --quiet --preservecase "
    if test.program == "placeholder":
        program_path = shutil.which("placeholder")
        assert(program_path is not None) # make sure placeholder is there!!
    os.system(program_path + variables + test.infile + " > " + test.outfile)
    os.system("cat "+ test.outfile)
      
main()



