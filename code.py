#cd "C:\Users\Sarah Hesham\OneDrive\Documents\Python\Biopython_CommandLineInterpreter"

#Some important remarks regarding the getopt library:

#An option followed by a colon only means that it needs an argument.
# It doesn't mean that the option is enforced. You should write your own code 
# to enforce the existence of options/arguments.

#getopt.getopt(args, shortopts, longopts=[])
#Parses command line options and parameter list. 
#args is the argument list to be parsed, without the leading 
#reference to the running program. Typically, this means sys.argv[1:]. 
#shortopts is the string of option letters that the script wants to recognize, 
#with options that require an argument followed by a colon (':'; i.e., the same 
#format that Unix getopt() uses).

#exception getopt.GetoptError
#This is raised when an unrecognized option is found in the argument list 
#or when an option requiring an argument is given none. The argument to the exception 
#is a string indicating the cause of the error. For long options, an argument given 
#to an option which does not require one will also cause this exception to be raised. 
#The attributes msg and opt give the error message and related option; if there is no 
#specific option to which the exception relates, opt is an empty string.

#One of the most important remarks is:
#that the getopt function reads options and forms a dictionary of each option and its argument if and only if the first 
#input was an option ,rather than that it will consider the whole input as an argument even if there are option
#passed and the opts array will be empty in that case
 
#Example Usage

#Short Form

#python3 code.py -s 2000-1-2 -e 2002-1-1

#Long Form

#python3 code.py --start_date 2000-1-2 --end_date 2002-1-1

#Input Formats:

#In this code many types of inputs are allowed and also many types of errors are well handeled
#Permitted Input Formats:
#command name param1 param2 param3 ... -o option1 -x option2 -y option3 ...
#command name -o option1 -x option2 -y option3 ...
#command name param1 param2 param3 ...
#So, your input can be only options associated with their arguments parameters or
#Only arguments
#Or a mix of options and arguments parameters
#In all these cases, many other logical errors are handeled according to the function type and functionality
#So as to give user freedom in their input format

#All required imports needed for built-in function to be called
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Seq import transcribe 
from Bio import pairwise2 
from Bio.Blast import NCBIXML 
from Bio.Blast import NCBIWWW 
from Bio import SeqIO
import sys
import getopt

#Functions Implementation 12 functions

#Number One
#Calculate G+C content, return percentage (as float between 0 and 100).
#Copes mixed case sequences, and with the ambiguous nucleotide S (G or C) 
#when counting the G and C content. The percentage is calculated against the full length
#gc 
#Usage :gc seq 
#Parameters seq: a string represents the sequence 
#Description :This command takes a seq as a string and returns the gc percentage of it. 
#Sample run: gc AGCAT
def gc(seq):

     return('GC= {}'.format(GC(seq)))

#Number Two
#Transcription is the process of changing DNA sequence into RNA sequence. 
#The actual biological transcription process is performing a reverse complement 
#(TCAG → CUGA) to get the mRNA considering the DNA as template strand. 
#However, in bioinformatics and so in Biopython, we typically work directly with 
#the coding strand and we can get the mRNA sequence by changing the letter T to U.
#transcribe 
#Usage :transcribe seq 
#Options and arguments seq: a string represents the sequence 
#Description: This command takes a seq as a string and returns its transcription. 
#Sample run: transcribe AGCAT 
def TRANSCRIBE(seq):

     return('Transcribed= {}'.format(transcribe(seq)))

#Number Three
#The reverse_complement() method complements and reverses the resultant sequence from left to right.
#I could have used the library Bio.Alpha but, it needs a downgraded version for biopython library 
#reverse_complement 
#Usage: reverse_complement seq 
#Options and arguments seq: a string represents the sequence 
#Description: This command takes a seq as a string and returns its reverse complement. 
#Sample run : reverse_complement AGCAT 
def REVERSE_COMPLEMENT(seq):

    complement_table = seq.maketrans('acgtnACGTN', 'tgcanTGCAN')
    complemented_sequence = seq.translate(complement_table)
    reverse_complement_seq = complemented_sequence[::-1]
    seq = reverse_complement_seq
    return ('Reverse_Complement= {}'.format(seq))

#Number Four
#This function takes a sequence and calculates its nbases which are letters N or n
#calc_nbases 
#Usage: calc_nbases seq 
#Options and arguments seq: a string represents the sequence 
#Description: This command takes a seq and calculates its nbases 
#Sample run: calc_nbases NGCTN 
def calc_nbases(seq):

    return ('N_bases Count = {}'.format(seq.count('n') + seq.count('N')))

#Number Five
#is_valid 
#Usage: is_valid seq type 
#Options and arguments seq: a string represents the sequence type: a string that represents the 
#type of the sequence. It can be one of these keywords [protein, dna, rna] 
#Description: This command takes a seq and a type (protein, dna, rna) and returns a 
#Boolean value of whether it’s a valid type or not 
#Sample run: is_valid NGCTN protein
def is_valid(seq, type):

    for base in seq:
        if base not in "ACTGactg" and type in ["DNA", "dna"]:
            return False
        if base not in "ACUGacug" and type in ["RNA", "rna"]:
            return False
        if base not in "ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz" and type in ["Protein", "protein"]:
            return False
    return True

#Number Six
#filter_nbases 
#Usage: filter_nbases seq 
#Options and arguments seq: a string represents the sequence 
#Description: This command takes a seq and returns the Seq after removing n bases 
#Sample run: filter_nbases NGCTN
def filter_nbases(seq):

    valid_seq = ""
    for base in seq:
        if base not in "nN":
            valid_seq += base
    return ('N Bases Filtered = {}'.format(valid_seq))

#Number Seven
#seq_alignment 
#Usage :seq_alignment seq1 seq2 [-o file] 
#Options and arguments :seq1, seq2 : strings represents the sequence -o file : specifies the path of the output file 
#Description :This command takes 2 sequences as input and get all their alignments 
#along with the score. The -o is an optional parameter if we need the output 
#to be written on a file instead of the screen. 
#Sample run :seq_alignment AGCCT AGC -o output.txt 
def seq_alignment(seq1 , seq2 , outputFile=None):

        alignments = pairwise2.align.globalxx(Seq(seq1.upper()), Seq(seq2.upper())) #global alignment

        if(outputFile != None):
             sys.stdout = open(outputFile, "w")

        for alignment in alignments:  
                print(alignment)

#Number Eight
#seq_alignment_files 
#Usage: seq_alignment_files file1 file2 [-o file3] 
#Options and arguments : file1, file2: specifies the paths of the files to be aligned that contain the 
#sequences -o file3 : specifies the path of the output file 
#Description: This command takes 2 fasta files as input, each file contains a single 
#sequence. It reads the 2 sequences from files and get all their alignments 
#along with the score. The -o is an optional parameter if we need the output 
#to be written on a file instead of the screen. 
#Sample run: seq_alignment s1.fasta s1.fasta -o output.txt 
def seq_alignment_files(file1,file2,outputFile=None):

        try:
           if(file1[-6:] != ".fasta"):
               print("Error, file 1 is not fasta extension")
               sys.exit()
           file1 = open(file1)
        except FileNotFoundError or ValueError or IndexError:
           print("Error , file 1 not found")
           sys.exit()

        try:
           if(file2[-6:] != ".fasta"):
               print("Error, file 2 is not fasta extension")
               sys.exit()
           file2 = open(file2)
        except FileNotFoundError or ValueError or IndexError:
           print("Error , file 2 not found")
           sys.exit()

        seq1 = file1.read()
        if(seq1.count(">") >1):
            print("Error, file one has more than one sequence")
            sys.exit()

        seq2 = file2.read()
        if(seq2.count(">") >1):
            print("Error, file two has more than one sequence")
            sys.exit()

        file1.seek(0)
        seq1 = file1.readline()

        if(seq1[0] != ">"):
            print("Error, wrong format of data in fasta file 1 input , no > sign found")
            sys.exit()

        seq1 =file1.read().split("\n")
        seq1 = ''.join(seq1) 

        file2.seek(0)
        seq2 =file2.readline()

        if(seq2[0] != ">"):
            print("Error, wrong format of data in fasta file 2 input , no > sign found")
            sys.exit()

        seq2 = file2.read().split("\n")
        seq2 = ''.join(seq2)

        alignments = pairwise2.align.globalxx(Seq(seq1.upper()), Seq(seq2.upper())) #global alignment

        if(outputFile != None):
            sys.stdout = open(outputFile, "w")

        for alignment in alignments:  
                print(alignment)


#Number Nine
#online_alignment 
#Usage: online_align seq [-o file] 
#Options and arguments : seq : a string represents the sequence -o file : specifies the path of the output file 
#Description :This command takes a sequence and uses BLAST to search the internet for 
#its alignments. The output should be all the information in the resultant 
#BLAST record. The -o is an optional parameter if we need the output to be 
#written on a file instead of the screen. 
#Sample run: seq_alignment ACTGCCGTCAAGTCAG -o output.txt
def online_alignment(seq, outputfile=None):

     result_handle = NCBIWWW.qblast("blastn","nt",seq.upper()) 
     blast_record = NCBIXML.read(result_handle)

     if(outputfile != None):
         sys.stdout = open(outputfile, "w")

     for alignment in blast_record.alignments: 
          for hsp in alignment.hsps: 
            print('****Alignment****') 
            print('sequence:', alignment.title) 
            print('length:', alignment.length) 
            print('e value:', hsp.expect) 
            print(hsp.query) 
            print(hsp.match) 
            print(hsp.sbjct)


#Number Ten
#merge_fasta 
#Usage :merge_fasta file1 file2 [file3 …] [-o output_file] 
#Options and arguments : file1, file2, etc. : specifies the paths of the files to be merged -o output_file 
#: specifies the path of the output file 
#Description: This command takes any number of fasta files (at least two) and merge 
#their contents into one fasta output file. There is an option to write the 
#merge result in a file using -o option, otherwise the merge result will be 
#displayed on the console. 
#Hint: use variadic parameters in functions to handle the unknown number 
#of parameters. Don’t use lists. 
#You may be thinking that the solution with *args is very similar to the list solution. 
#That's because *args is internally a Tuple, which is an iterable sequence similar to lists.
#Sample run: merge_fasta f1.fasta f2.fasta f3.fasta f4.fasta
def merge_fasta(file1,file2,outputfile=None,*args):

        text = ""

        try:
           if(file1[-6:] != ".fasta"):
               print("Error, file 1 is not fasta extension")
               sys.exit()
           else:
              file1 = open(file1)
              text = text + file1.read()
        except FileNotFoundError or IndexError or ValueError:
           print("Error , file 1 not found")
           sys.exit()

        try:
           if(file2[-6:] != ".fasta"):
               print("Error, file 2 is not fasta extension")
               sys.exit()
           file2 = open(file2)
           text = text + "\n" + file2.read()
        except FileNotFoundError or IndexError or ValueError:
           print("Error , file 2 not found")
           sys.exit()
           
        for Onefile in args:
           try:
               if(Onefile[-6:] != ".fasta"):
                   print("Error, an extra file is not fasta extension")
                   sys.exit()
               temp = open(Onefile)
               text = text + "\n" + temp.read()
           except FileNotFoundError or IndexError or ValueError:
              print("Error , extra file  not found")
              sys.exit()

        if(outputfile!=None):
           sys.stdout = open(outputfile, "w")

        print(text)

#Number Elleven
#convert_to_fasta 
#Usage: convert_to_fasta file 
#Options and arguments: file: specifies the path of a genbank file 
#Description :This command converts the input genbank file with multiple records onto a 
#fasta formatted file. The output is to be written in a different output fasta file. 
#Sample run: convert_to_fasta “ls-orchid.gbk” 
def convert_to_fasta(file1):

    count = 0 

    try:
        if(file1[-3:] != ".gb"):
           print("Error, wrong file extension ,it is not a gene bank file")
           sys.exit()

    except ValueError or IndexError:
        print("Error, wrong file extension ,it is not a gene bank file")
        sys.exit()

    try:
        with open(file1) as input_handle, open( str(file1[0:-3] + '.fasta'), "w" ) as output_handle: 
           sequences = SeqIO.parse(input_handle, "genbank") 
           count = SeqIO.write(sequences, output_handle, "fasta") 
    except FileNotFoundError:
        print("Error, no such file found")
        sys.exit()

    return ("Converted %i records" % count)

#Number Twellve
#Main Function
#Calls other functions and parses input using getopt library by handling all errors possible
def main():

    # first argument is the filename with an index of 0
    # so, we want to start with an index value of 1
    argv = sys.argv[1:]

    #implemented functions available in this program to be called
    commands_List = ["gc","transcribe","reverse_complement","calc_nbases","is_valid","filter_nbases","seq_alignment","seq_alignment_files","online_alignment","merge_fasta","convert_to_fasta"]
    command =""

    try:
        command = argv[0]

    except IndexError:
        print("Error, you need to enter the command and its parameters, nothing found")
        sys.exit()

    if(command not in commands_List):
        print("Error, unknown command entered or wrong spelling, check format entered")
        sys.exit()

    if(command == "gc"):

         try:
            opts, args = getopt.getopt(argv[1:],"s:",["sequence="])

         except getopt.GetoptError as err:
            print(err)  # will print something like "option -a not recognized"
            sys.exit()
         
         if (opts == []):
             if(args == []):
                 print("Error, function gc needs a single argument as input and nothing was found")
                 sys.exit()
             elif(len(args) > 1):
                 print("Error, function gc needs a single argument as input and many was given")
                 sys.exit()
             else:
                print(gc(args[0]))
         else:
             for opt, arg in opts:
                if(len(args) > 0):
                    print("Error, function gc needs a single argument as input and many was given")
                    sys.exit()
                else:
                    if opt in ["-s", "--sequence"] :
                       seq = arg
                       print(gc(seq))

    elif(command == "transcribe"):

         try:
            opts, args = getopt.getopt(argv[1:],"s:",["sequence="])

         except getopt.GetoptError as err:
            print(err)  # will print something like "option -a not recognized"
            sys.exit()
         
         if (opts == []):
             if(args == []):
                 print("Error, function transcribe needs a single argument as input and nothing was found")
                 sys.exit()
             elif(len(args) > 1):
                 print("Error, function transcribe needs a single argument as input and many was given")
                 sys.exit()
             else:
                print(TRANSCRIBE(args[0]))
         else:
             for opt, arg in opts:
                if(len(args) > 0):
                    print("Error, function transcribe needs a single argument as input and many was given")
                    sys.exit()
                else:
                    if opt in ["-s", "--sequence"] :
                       seq = arg
                       print(TRANSCRIBE(seq))

    elif(command == "reverse_complement"):

         try:
            opts, args = getopt.getopt(argv[1:],"s:",["sequence="])

         except getopt.GetoptError as err:
            print(err)  # will print something like "option -a not recognized"
            sys.exit()
         
         if (opts == []):
             if(args == []):
                 print("Error, function reverse_complement needs a single argument as input and nothing was found")
                 sys.exit()
             elif(len(args) > 1):
                 print("Error, function reverse_complement needs a single argument as input and many was given")
                 sys.exit()
             else:
                print(REVERSE_COMPLEMENT(args[0]))
         else:
             for opt, arg in opts:
                if(len(args) > 0):
                    print("Error, function reverse_complement needs a single argument as input and many was given")
                    sys.exit()
                else:
                    if opt in ["-s", "--sequence"] :
                       seq = arg
                       print(REVERSE_COMPLEMENT(seq))

    elif(command == "calc_nbases"):

         try:
            opts, args = getopt.getopt(argv[1:],"s:",["sequence="])

         except getopt.GetoptError as err:
            print(err)  # will print something like "option -a not recognized"
            sys.exit()
         
         if (opts == []):
             if(args == []):
                 print("Error, function calc_nbases needs a single argument as input and nothing was found")
                 sys.exit()
             elif(len(args) > 1):
                 print("Error, function calc_nbases needs a single argument as input and many was given")
                 sys.exit()
             else:
                print(calc_nbases(args[0]))
         else:
             for opt, arg in opts:
                if(len(args) > 0):
                    print("Error, function calc_nbases needs a single argument as input and many was given")
                    sys.exit()
                else:
                    if opt in ["-s", "--sequence"] :
                       seq = arg
                       print(calc_nbases(seq))

    elif(command == "is_valid"):

         try:
            opts, args = getopt.getopt(argv[1:],"s:t:",["sequence=","type="])

         except getopt.GetoptError as err:
            print(err)  # will print something like "option -a not recognized"
            sys.exit()
        
         if (len(opts) + len(args) != 2):
                print("Error, you must use options -s and -t or long options --sequence and --type or enter the 2 arguments parameters as this function has 2 arguments needed with values")
                sys.exit()

         seq = type = " "    
         for opt, arg in opts:
                if opt in ["-s", "--sequence"] :
                    seq = arg
                elif opt in ["-t","--type"]:
                    type = arg

         if(seq !=" " and type != " "):
               if(type in ["dna","rna","protein"] or type in ["DNA","RNA","PROTEIN"]):
                   print('Validity is= {}'.format(is_valid(seq,type)))
               else:
                   print("Error, second argument must make sense and be (rna or dna or protein)")
                   sys.exit()                   
                   
         if(args!=[]):  
       
             if(seq == " "):
                seq = args[0]

             if (type == " "):
                try:
                   type = args[1]
                except IndexError:
                   type = args[0]

             if(type in ["dna","rna","protein"] or type in ["DNA","RNA","PROTEIN"]):
                 print('Validity is= {}'.format(is_valid(seq,type)))
             else:
                 print("Error, second argument must make sense and be (rna or dna or protein)")
                 sys.exit()                   
        
    elif(command == "filter_nbases"):

         try:
            opts, args = getopt.getopt(argv[1:],"s:",["sequence="])

         except getopt.GetoptError as err:
            print(err)  # will print something like "option -a not recognized"
            sys.exit()
         
         if (opts == []):
             if(args == []):
                 print("Error, function filter_nbases needs a single argument as input and nothing was found")
                 sys.exit()
             elif(len(args) > 1):
                 print("Error, function filter_nbases needs a single argument as input and many was given")
                 sys.exit()
             else:
                print(filter_nbases (args[0]))
         else:
             for opt, arg in opts:
                if(len(args) > 0):
                    print("Error, function gc needs a single argument as input and many was given")
                    sys.exit()
                else:
                    if opt in ["-s", "--sequence"] :
                       seq = arg
                       print(filter_nbases(seq))

    elif(command == "seq_alignment"):

         try:
            opts, args = getopt.getopt(argv[1:],"1:2:o:",["first_sequence=","second_sequence=","outputfile="])

         except getopt.GetoptError as err:
            print(err)  # will print something like "option -a not recognized"
            sys.exit()
        
         if (len(opts) + len(args) < 2) or (len(opts) + len(args) > 3 and len(opts) > 1) or (len(opts) == 0 and ('-1' in args or '-2' in args or '-o' in args)):
                print("Error, you must use options -1 and -2 or long options --first_sequence and --second_sequence or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()
   
         seq1 = seq2 = outputfile = None

         if(len(args)==4):
             if(args[0][0] != '-' and args[0][0] !='--'):
                 seq1 = args[0]
             else:
                print("Error, you must use options -1 and -2 or long options --first_sequence and --second_sequence or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()
             if(args[1][0] != '-' and args[1][0] !='--'):
                 seq2 = args[1]
             else:
                print("Error, you must use options -1 and -2 or long options --first_sequence and --second_sequence or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()
             if(args[2] == '-o' and args[3][0] != '-' and args[3][0] !='--'):
                 outputfile = args[3]
                 seq_alignment(seq1,seq2,outputfile)
                 sys.exit()
             else:
                print("Error, you must use options -1 and -2 or long options --first_sequence and --second_sequence or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()

         for opt, arg in opts:
                if opt in ["-1", "--first_sequence"] :
                    seq1 = arg
                elif opt in ["-2","--second_sequence"]:
                    seq2 = arg
                elif opt in ["-o","--outputfile"] :
                    outputfile = arg 
         
         flag0 = flag1 =0
         if(args!=[]):         
             if(seq1 == None):
                seq1 = args[0]
                flag0 = 1 
             if (seq2 == None):
                try:
                   if(flag0 == 1 ) :
                      seq2 = args[1]
                      flag1 =1
                   else:
                       seq2=args[0]
                       flag0 = 1
                except IndexError:
                   seq2 = args[0]
                   flag0 =1
             if (outputfile == None):
                try:
                   outputfile = args[2]
                except IndexError:
                   try:
                      if(flag1 ==0 ):
                         outputfile = args[1]    
                   except IndexError:
                       if(flag0 == 0):
                          outputfile = args[0]
         seq_alignment(seq1,seq2,outputfile)

    elif(command == "seq_alignment_files"):

         try:
            opts, args = getopt.getopt(argv[1:],"1:2:o:",["first_file=","second_file=","outputfile="])

         except getopt.GetoptError as err:
            print(err)  # will print something like "option -a not recognized"
            sys.exit()
        
         if (len(opts) + len(args) < 2) or (len(opts) + len(args) > 3 and len(opts) > 1 ):
                print("Error, you must use options -1 and -2 or long options --first_file and --second_file or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()
   
         seq1 = seq2 = outputfile = None

         if(len(args)==4):
             if(args[0][0] != '-' and args[0][0] !='--'):
                 seq1 = args[0]
             else:
                print("Error, you must use options -1 and -2 or long options --first_file and --second_file or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()
             if(args[1][0] != '-' and args[1][0] !='--'):
                 seq2 = args[1]
             else:
                print("Error, you must use options -1 and -2 or long options --first_file and --second_file or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()
             if(args[2] == '-o' and args[3][0] != '-' and args[3][0] !='--'):
                 outputfile = args[3]
                 seq_alignment_files(seq1,seq2,outputfile)
                 sys.exit()
             else:
                print("Error, you must use options -1 and -2 or long options --first_file and --second_file or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()

         for opt, arg in opts:
                if opt in ["-1", "--first_file"] :
                    seq1 = arg
                elif opt in ["-2","--second_file"]:
                    seq2 = arg
                elif opt in ["-o","--outputfile"] :
                    outputfile = arg 
                    if((seq2 ==None or seq1==None) and len(opts)!=3):
                       print("Error, you must use options -1 and -2 or long options --first_file and --second_file or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                       sys.exit()

         
         flag0 = flag1 =0
         if(args!=[]):         
             if(seq1 == None):
                seq1 = args[0]
                flag0 = 1 
             if (seq2 == None):
                try:
                   if(flag0 == 1) :
                      seq2 = args[1]
                      flag1 =1
                   else:
                       seq2=args[0]
                       flag0 = 0
                except IndexError:
                   seq2 = args[0]
                   flag0 =1
             if (outputfile == None):
                try:
                   outputfile = args[2]
                except IndexError:
                   try:
                      if(flag1 ==0 ):
                         outputfile = args[1]    
                   except IndexError:
                       if(flag0 == 0):
                          outputfile = args[0]
         seq_alignment_files(seq1,seq2,outputfile)

    elif(command == "online_alignment"):

         try:
            opts, args = getopt.getopt(argv[1:],"s:o:",["sequence=","outputfile="])

         except getopt.GetoptError as err:
            print(err)  # will print something like "option -a not recognized"
            sys.exit()
        
         if (len(opts) + len(args) < 1) or (len(opts) + len(args) > 2 and len(opts) > 1 ):
                print("Error, you must use options -s and -o or long options --sequence and --outputfile or enter the one parameter as this function has 1 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()
   
         seq = outputfile = None

         if(len(args)==3):
             if(args[0][0] != '-' and args[0][0] !='--'):
                 seq = args[0]
             else:
                print("Error, you must use options -1 and -2 or long options --first_file and --second_file or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()
             if(args[1] == '-o' and args[2][0] != '-' and args[2][0] !='--'):
                 outputfile = args[2]
                 online_alignment(seq,outputfile)
                 sys.exit()
             else:
                print("Error, you must use options -s and -o or long options --sequence and --outputfile or enter the one parameter as this function has 1 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()

         for opt, arg in opts:
                if opt in ["-s", "--sequence"] :
                    seq = arg
                elif opt in ["-o","--outputfile"] :
                    outputfile = arg 
                    if(seq == None and len(opts)!=2):
                       print("Error, you must use options -s and -o or long options --sequence and --outputfile or enter the one parameter as this function has 1 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                       sys.exit()

         flag0 = 0
         if(args!=[]):         
             if(seq == None):
                seq = args[0]
                flag0 = 1 

             if (outputfile == None):
                try:
                   outputfile = args[1]
                except IndexError:
                      if(flag0 ==0 ):
                         outputfile = args[0]    
                       
         online_alignment(seq,outputfile)

    elif(command == "merge_fasta"):

         try:
            opts, args = getopt.getopt(argv[1:],"o:",["outputfile="])

         except getopt.GetoptError as err:
            print(err)  # will print something like "option -a not recognized"
            sys.exit()
        
         if len(args) < 2 or (len(args) < 4 and len(opts)==0 and "-o" in args) :
                print("Error, you must use option -o or long option --outputfile or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()
   
         outputfile = None
         file1 = file2 =""
         for opt, arg in opts:
                if opt in ["-o","--outputfile"] :
                    outputfile = arg 
                    file1 = args[0]
                    file2 = args[1]
                    try :
                       myTuple = (i for i in args[2:])
                       merge_fasta(file1,file2,outputfile,*myTuple)
                       sys.exit()
                    except:
                       pass

         if(len(opts) == 0):

             index = None
             try:
                index = args.index('-o',0,len(args)+1)
             except ValueError:
                index = 2000
             try:
               outputfile = args[index+1]
             except IndexError:
                if(index != 2000):
                    print("Error, -o must be followed by a file name")
                    sys.exit()
                else:
                    outputfile = None
             try:
                file1 = args[0]
                file2 =args[1]
             except IndexError:
                print("Error, you must use option -o or long option --outputfile or enter the two parameters as this function has 2 arguments needed with values minimum and maximum add -o or --outputfile or enter the parameter only not more")
                sys.exit()
             try:
                 if(outputfile !=None):
                    myTuple = (i for i in args[2:-2])
                    merge_fasta(file1,file2,outputfile,*myTuple)
                    sys.exit()
                 else:
                    myTuple = (i for i in args[2:])
                    merge_fasta(file1,file2,outputfile,*myTuple)
                    sys.exit()
             except:
                 pass

    elif(command == "convert_to_fasta"):

         try:
            opts, args = getopt.getopt(argv[1:],"f:",["file="])

         except getopt.GetoptError as err:
            print(err)  # will print something like "option -a not recognized"
            sys.exit()
         
         if (opts == []):
             if(args == []):
                 print("Error, function convert_to_fasta needs a single argument as input and nothing was found")
                 sys.exit()
             elif(len(args) > 1):
                 print("Error, function convert_to_fasta needs a single argument as input and many was given")
                 sys.exit()
             else:
                print(convert_to_fasta(args[0]))
         else:
             for opt, arg in opts:
                if(len(args) > 0):
                    print("Error, function convert_to_fasta needs a single argument as input and many was given")
                    sys.exit()
                else:
                    if opt in ["-f", "--file"] :
                       seq = arg
                       print(convert_to_fasta(seq))

    else:
        print("Error in input commands and parameters")
        sys.exit()

#Here we start exetcuting the program
main()