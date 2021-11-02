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
        check_infile_exists()
        check_program_exists()
    
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
       # Spielman: NOPE! These are attributes of the Settings class. NOT SEQUENCES!
       #self.infile = infile
       # self.outfile = outfile
       # self.outtranslated = outtranslated


    def prepare_for_alignment(self, infile, temporary_file):
        self.translate_data(infile)
        self.write_translated_seqs_to_file(temporary_file)                
                
                
        
    def read_sequences_from_file(self, input_file, format = "fasta"):
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
                    
        
    def translate_sequences(self, infile): # SPIELMAN: See, we just pass infile as an argument. It's not a self.infile
        """
            Translating sequences from nucleotides to amino acids. Creates and fills two dictionaries, self.unaligned and self.unaligned_translated
        """
        records = self.read_sequences_from_file(infile)
        self.unaligned = self.convert_sequence_list_to_dictionary(records)
        self.unaligned_translated = {}
        for record_id in self.unaligned:
            self.unaligned[record_id] = self.unaligned.seq[record_id].translate()
    
    
    def write_translated_seqs_to_file(self, temporary_file):
        """
            Creates the file for aligned output and appends the translated data to it
        """
        os.system("touch " + temporary_file) # Spielman: ANY LUCK WITH TMP?
        # Spielman: BELOW DOES NOT WORK, and "append" mode is DANGER! DANGER!
        #file = open(self.unaligned_translated, "a")  # append mode
        #file.write(str(self.unaligned_translated))
        #file.close()
        with open(temporary_file, "w") as f:
            for record in self.unaligned_translated:
                f.write(">" + str(record) + "\n" + str(self.unaligned_translated[record]) + "\n")
        
    # def backtranslate_sequences(self, temporary_file)
    # def save_sequences_to_file(self, aafile, nucfile)
                





    # __call__ method here doesn't make "sense." What does it mean to <verb> a Sequence instance?
    # def __call__(self, infile):
    #     '''
    #         Runs both functions when an instance is called as a function.
    #     '''
    #     self.translate_data(infile) 
    #     self.write_translated_seqs_to_file()
   
   













class Aligner:
    """
        Creating the Aligner class. Will also take arguments from the arguments specified from collect_input_arguments(). Generates what will be called in os.system to run the alignment.
    """
    def __init__(self, settings):
        '''
            Sets up Aligner instance.
        '''

        self.program_path = shutil.which(settings.program)
        self.command = (self.program_path + " " + settings.arguments + " " + settings.infile + ">" + settings.outfile)
    
    def __call__(self):
        os.system(self.command)













def main():
    # Collect input
    args = collect_input_arguments() 
    my_settings = Settings(args.inf, args.outf, args.prog, args.args, args.outtranslated)
    
    # Initialize sequences and prepare for alignment
    my_sequences = Sequences() 
    my_sequences.prepare_for_alignment(my_settings.infile, "SOME TEMPORARY FILE GOES HERE????")
    
    # Perform alignment
    my_aligner = Aligner(my_settings)
    my_aligner()
    
    # Backtranslate
    my_sequences.read_in_aligned_data(my_settings.outfile)   
    # my_sequences.perform_backtranslation()
    
    # Export
    # my_sequences.save_to_file()

main()



