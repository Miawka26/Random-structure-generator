import random
import os
import argparse
import sys
import subprocess
import re
from ushuffle import shuffle, Shuffler 

def random_sequences(gc):
    """The function randomizes without returning nucleobases from the list containing the enumerated amount of Guanine, Cytosine, Adenine, Uracil based on GC-content.
    
    Args:
       gc (int): Percentage of gc that the user specified on the input.
    
    Returns:
        y (list): Returns the drawn nucleobase.
    """
    nucleobases = []

    Guanine = gc * 0.1
    Cytosine = gc * 0.1
    Adenine = (10 - Guanine)
    Uracil = (10 - Cytosine)
    Thymine= (10 - Cytosine)
    
    Guanine1 = 'G ' * int(Guanine)
    nucleobases.append(Guanine1.split())
    Cytosine1 = 'C ' * int(Cytosine)
    nucleobases.append(Cytosine1.split())
    Adenine1 = 'A ' * int(Adenine)
    nucleobases.append(Adenine1.split())
    Uracil1 = 'U ' * int(Uracil)
    nucleobases.append(Uracil1.split())

    flat_list = []
    for row in nucleobases:
        flat_list.extend(row)
    y = random.choices(flat_list) 
    return y

def designating_a_place_with_range(length_range):
    """The function takes a range as string, extracts numbers from it, and then randomizes one number from the range to be the length.
    
    Args:
        length_range (str): Length range that the user specified on the input.
    
    Returns:
        z(int): The length that the program drew from the range.
    """
    range1 = length_range.split("-")
    beginning_of_range = int(range1[0])
    end_of_range = int(range1[1])
    z = random.randint(beginning_of_range, end_of_range)
    return z


def checking_fasta(path):
    """Function that checks the correctness of the entered fasta file.
    
    Args:
        path (str): The path to the fasta file entered by the user as input.

    Returns:
        nucleotide(str): Variable containing information about the nucleotide in the sequence - thymine or uracil.
    """
    with open(path) as file:
        if os.stat(path).st_size==0:
            sys.stderr.write("The file '{}' is empty.".format(path))
            sys.exit(0)  
        first_line = file.readline()
        if not first_line.startswith('>'):
            sys.stderr.write("The file '{}' should start with the identifier: \'>...\'".format(path))
            sys.exit(0)

    with open(path) as file:
        z=0
        for line in file:
            if line.startswith('>'):
                z=z+1
            if not line.startswith(('>','.',')','(')):
                    z=z+1
            if line.startswith(('.',')','(')):
                z=z+1
                if z!=3:
                    sys.stderr.write('Scheme: ID, sequence, dot-bracket was not preserved in the file.')
                    sys.exit(0)
                z=0
    
    with open(path) as file:
        for line in file:
            if line.startswith('(' or ')' or '.'):
                for a in line:
                    if a not in ['.','(',')','\n']:             
                        sys.stderr.write("'{}' file format error. Invalid character in dotbracket '{}'\n".format(path, a))
                        sys.exit(0)

    with open(path) as file:
        for line in file:
          if not line.startswith(('>','.',')','(')):
                    for a in line:
                        if a not in ['G','C','A','U','\n', 'T']:         
                                sys.stderr.write("'{}' file format error. Invalid character '{}' in sequence '{}' \n".format(path, a, line))
                                sys.exit(0)
    Thymine=0
    Uracil=0
    with open(path) as file:
        for line in file:
          if not line.startswith(('>','.',')','(')):
                    for a in line:
                        if a in ['U']: 
                        	    Uracil=Uracil+1
                        else:
                                Uracil=Uracil+0
                        if a in ['T']:         
                                Thymine=Thymine+1
                        else:
                                Thymine=Thymine+0
        if Thymine>0 and Uracil>0:
            sys.stderr.write('The sequence cannot contain thymine and uracil at the same time.')
            sys.exit(0)
        if Thymine>0:
            nucleotide='T'
        else:
        	nucleotide='U'
    return nucleotide

def reading_fasta_files(path):
    """Function that reads a fasta file and adds its content to a dictionary consisting of: identifier: sequence, dot bracket.
    
    Args:
        path (str): The path to the fasta file entered by the user as input.
    
    Returns:
        seq_struc_dot(dict): Dictionary consisting of: id: sequence, dot bracket.
    """
    with open(path) as file:
        seq_struc_dot = {}
        for line in file:
            if line.startswith('>'):
                identifiers = line.strip()
            if not line.startswith(('>','.',')','(')):
                fasta_sequences = line.strip()
            if line.startswith(('.',')','(')):
                dot_bracket = line.strip()
                seq_struc_dot[identifiers] = fasta_sequences, dot_bracket
    return seq_struc_dot


def loading_an_additional_file(path1):
    """The function reads an additional fasta file containing the sequences to be shuffled and writes the sequences from the file to the variable.
    
    Args:
        path1 (str): Path to an additional fasta file containing sequences to be shuffled entered by the user at the input.
    
    Returns:
        sequences_to_shuffle (str): Variable containing the sequences to be shuffled.
        shuffle_nucleotide (str): Variable containing information about the nucleotide in the sequence to be shuffled - thymine or uracil.
    """
    sequences_to_shuffle1=[]
    with open(path1) as file:
        if os.stat(path1).st_size==0:
            sys.stderr.write("The file '{}' is empty.".format(path1))
            sys.exit(0)  

    Thymine=0
    Uracil=0
    with open(path1) as file:
        for line in file:
          if not line.startswith(('>','.',')','(')):
                    for a in line:
                        if a in ['U']: 
                        	    Uracil=Uracil+1
                        else:
                                Uracil=Uracil+0
                        if a in ['T']:         
                                Thymine=Thymine+1
                        else:
                                Thymine=Thymine+0
        if Thymine>0 and Uracil>0:
            sys.stderr.write('The sequence cannot contain thymine and uracil at the same time.')
            sys.exit(0)
        if Thymine>0:
            shuffle_nucleotide='T'
        else:
        	shuffle_nucleotide='U'
 

    with open(path1) as file:
        for line in file:
            if not line.startswith('>'):
                if not line.startswith('(' or ')' or '.'):
                    sequences_to_shuffle1.append(line.rstrip('\n'))
    sequences_to_shuffle=(sequences_to_shuffle1[0])
    return sequences_to_shuffle, shuffle_nucleotide



def creating_random_sequences(sequence_length):
    """A function that creates a random sequence of a specified length.
    
    Args:
        sequence_length (int): Variable containing the length of the sequence.
    
    Returns:
        lst(list): List containing sequences.
    """   
    lst = []
    for _ in range(sequence_length):
        lst.append(random_sequences(GC_content)[0])
    return lst


def creating_random_sequences__RANGE(length_range):
    """A function that creates a random sequence of a specified length.
    
    Args:
       length_range (str): Length range that the user specified on the input.
    
    Returns:
        lst(list): List containing sequences.
    """
    lst = []
    sequence_length= designating_a_place_with_range(length_range)
    for _ in range(sequence_length):
        lst.append(random_sequences(GC_content)[0])
    return lst


def number_of_sequences(number, sequence_length):
    """A function that generates the number of sequences given by the user.
    
    Args:
        number (int): Number of sequences given by the user as input.
        sequence_length (int): Variable containing the length of the sequence.
    
    Returns:
        sequences(list): A list containing the given number of sequences.
    """
    sequences = []
    for _ in range(number):
        sequences.append(creating_random_sequences(sequence_length))
    return sequences


def number_of_sequences__RANGE(number, length_range):
    """A function that generates the number of sequences given by the user.
    
    Args:
       number (int): Number of sequences given by the user as input.
       length_range (str): Length range that the user specified on the input.
    
    Returns:
        sequences(list): A list containing the given number of sequences.
    """
    sequences = []
    for _ in range(number):
        sequences.append(creating_random_sequences__RANGE(length_range))
    return sequences


def sequence_shuffling(number, sequences_to_shuffle, shuffle_option):
    """A function that generates the number of sequences given by the user.
    
    Args:
       number (int): Number of sequences given by the user as input.
       sequences_to_shuffle (str): Variable containing the sequences to be shuffled.
       shuffle_option (int): A variable containing the users choice of how the sequences should be shuffled.

    Returns:
       sequences(list): Variable containing the shuffled sequence.
    """
    sequences = []
    for _ in range(number):
        if shuffle_option == 2:
            sequence_shuffled1=shuffle(sequences_to_shuffle.encode('utf-8'), 2).decode("utf-8")
        else:
            sequence_shuffled1=shuffle(sequences_to_shuffle.encode('utf-8'), 3).decode("utf-8")
        sequences.append(list(sequence_shuffled1))
    return sequences



def writing_to_txt_file(sequence, output_file_name, Thymine_or_uracil):
    """Function that writes sequences to a .txt file.
    
    Args:
        sequences (list): A list containing the given number of sequences.
        output_file_name (str): Variable containing the given file name entered by the user on input.
        Thymine_or_uracil (str): Choosing what is in the sequence - Thymine or Uracil.
    """
    nucleobases1 = ""  
    filename_txt=output_file_name+'.txt'
    with open(filename_txt, 'w') as txt_file:
        for nucleobase in sequence:
            if Thymine_or_uracil == 'U':
                txt_file.write((nucleobases1.join(nucleobase)).replace('T', 'U') + '\n')
            else:
                txt_file.write((nucleobases1.join(nucleobase)).replace('U', 'T') + '\n')


def random_position_in_sequence(sequence_length, fragment_length):
    """A function that calculates the positions at which a fragment can be inserted, and then randomly picks one of them.
    
    Args:
        sequence_length (int): Variable containing the length of the sequence.
        fragment_length (int): Variable containing the length of the fragment to be inserted.
    
    Returns:
        random_place(int): Variable containing random position on which the fragment will be inserted.
    """
    random1=sequence_length-fragment_length
    random_place=random.randint(0,random1)
    return random_place


def inserting_known_sequences(sequences, number, insertion_position):
        """The function that inserts user-specified fragments in the fasta file at the appropriate positions.
        
        Args:
            sequences (list): A list containing the given number of sequences.
            number (int): Number of sequences given by the user as input.
            insertion_position (int): Variable containing the fragment insertion position specified by the user as input.
        
        Returns:
            random_to_bed(list) - Variable containing the positions on which the fragments were inserted.
            fragment_length_to_bed(list) - Variable containing the length of inserted fragments.
            number_of_fragments(int) - Variable containing the number of inserted fragments.
            inserted(list)- Complete sequences with fragments inserted.
        """
        random_to_bed=[]
        fragment_length_to_bed=[]
        only_sequences_from_fasta1=[]
        fragment_length=[]
        inserted=[]
        number_of_fragments= calculating_the_number_of_fragments_in_a_fasta(seq_struc_dot)
        for fasta_sequences in seq_struc_dot.values():
            only_sequences_from_fasta1.append(list(fasta_sequences[0]))
        for i in range(number):
            sequence_length=len(sequences[i])
            if insertion_position!=None:
                random_place=insertion_position
            else:
                if i<number_of_fragments:
                    fragment_length=len(only_sequences_from_fasta1[i])
                    if fragment_length>sequence_length:
                        sys.stderr.write("Sequence '{}' cannot be inserted. Total length ({}) is shorter than fragment length ({}).".format(''.join(only_sequences_from_fasta1[i]),sequence_length,fragment_length))
                        sys.exit(0)
                    random_place=random_position_in_sequence(sequence_length,fragment_length)
                if i>=number_of_fragments:
                    random_place = None
 
            if i<number_of_fragments:
                random_sequences_from_fasta=("".join(((sequences[i][:random_place]+only_sequences_from_fasta1[i])+sequences[i][random_place:])[:sequence_length]))
                random_to_bed.append(random_place)
                inserted.append(random_sequences_from_fasta)
                fragment_length_to_bed.append(len(only_sequences_from_fasta1[i]))
 
            if i>=number_of_fragments:
                inserted.append("".join(sequences[i]))
        return random_to_bed, fragment_length_to_bed, number_of_fragments, inserted


def calculating_the_number_of_fragments_in_a_fasta(seq_struc_dot):  
    """A function that calculates the number of fragments in fasta.

    Args:
        seq_struc_dot(dict): Dictionary consisting of: id: sequence, dot bracket.

    Returns:
        number_of_fragments1(int): Variable containing the number of inserted fragments.
    """
    number_of_fragments=str(seq_struc_dot.keys())
    number_of_fragments1=number_of_fragments.count('>') 
    return number_of_fragments1


def calculating_the_total_number_of_fragments(percent,seq_struc_dot):
    """Calculating the number of sequences based on the percentage of known fragments of the total sequence.
    
    Args:
        percent (int): The percentage of known fragments relative to all sequences.
        seq_struc_dot(dict): Dictionary consisting of: id: sequence, dot bracket.
    
    Returns:
        amount_of_all_fragments(int): Total number of sequences
    """
    number_of_fragments2=calculating_the_number_of_fragments_in_a_fasta(seq_struc_dot)
    amount_of_all_fragments=int((number_of_fragments2*100)/percent)
    return amount_of_all_fragments
    
   
def BED_format(number, random_to_bed, fragment_length_to_bed, number_of_fragments, insertion_position):
    """The function saves data to a file in the .bed format.
    
    Args:
        number (int): Number of sequences given by the user as input.
        random_to_bed (list) - Variable containing the positions on which the fragments were inserted.
        fragment_length_to_bed (list): Variable containing the length of inserted fragments.
        number_of_fragments(int) - Variable containing the number of inserted fragments.
        insertion_position (int): Variable containing the fragment insertion position specified by the user as input.
    
    Returns:
        BED1(list): Variable containing the content of the file with the extension .bed.
    """
    keys=[]
    BED1=[]
    keys=list(seq_struc_dot.keys())
    for i in range(number):
        if i<number_of_fragments:
            if insertion_position is None:
                BED='{}{}   {}  {}  {}  0   +'.format('rand_',i,random_to_bed[i],random_to_bed[i]+fragment_length_to_bed[i],keys[i].lstrip('>'))
                BED1.append(BED)
            else:
                end_position=(insertion_position+fragment_length_to_bed[i])
                BED='{}{}   {}  {}  {}  0   +'.format('rand_',i,insertion_position,end_position,keys[i].lstrip('>'))
                BED1.append(BED)           
    return BED1


def ID_record(number, random_to_bed, fragment_length_to_bed, number_of_fragments):
    """A function that creates fasta identifiers based on the given data.
    
    Args:
        number (int): Number of sequences given by the user as input.
        random_to_bed (list): Variable containing the positions on which the fragments were inserted.
        fragment_length_to_bed (list): Variable containing the length of inserted fragments.
        number_of_fragments(int): Variable containing the number of inserted fragments.
    
    Returns:
        ID(list): List of fasta identifiers based on the provided data.
    """
    keys=[]
    ID=[]
    keys=list(seq_struc_dot.keys())
    for i in range(number):
        if i<number_of_fragments:
            ID1='>{}{}|{}-{}|{}'.format('rand_', i,random_to_bed[i],random_to_bed[i]+fragment_length_to_bed[i],keys[i].lstrip('>'))
            ID.append(ID1)
        else:
            ID1='>{}{}'.format('rand_',i)
            ID.append(ID1)
    return ID

def ID_record_without_file(number):
    """Function that creates fasta identifiers with sequence number only.
    
    Args:
        number (int): Number of sequences given by the user as input.
    
    Returns:
         ID(list): List of fasta identifiers with sequence number only.
    """
    ID=[]
    for i in range(number):
        ID1='{}{}{}'.format('>','rand_',i)
        ID.append(ID1)
    return ID


def writing_to_BED_file(BED, output_file_name):
    """The function saves data to a file in the .bed format.
    
    Args:
        BED (list):  Variable containing the content of the file with the extension .bed.
        output_file_name (str): Variable containing the given file name entered by the user on input.
    """
    file_name=output_file_name+'.bed'
    with open(file_name, 'w') as bed_file:
        for a in BED:
            bed_file.write(a + '\n')


def writing_to_FASTA_file(sequence, output_file_name,ID,number, Thymine_or_uracil):
    """The function saves data to a file in the .fasta format.
    
    Args:
        sequence (list): The sequence to be written to the fasta file.
        output_file_name (str): Variable containing the given file name entered by the user on input.
        ID(list): List of fasta identifiers with sequence number only.
        number (int): Number of sequences given by the user as input.
        Thymine_or_uracil (str): Choosing what is in the sequence - Thymine or Uracil.
    
    Returns:
        sequence4(list): Variable containing the content of the file with the extension .fasta.
    """
    sequence3=''
    sequence4=sequence
    fasta_file_name=(output_file_name+'_sequences.fasta')
    with open(fasta_file_name, 'w') as file_with_seq:
        number_of_seq=len(sequence)
        for i in range(number):
            sequence2=sequence[i]
            file_with_seq.write(str(ID[i]) + '\n')
            if Thymine_or_uracil == 'U':
                file_with_seq.write((sequence3.join(sequence2)).replace('T', 'U') + '\n')
            else:
            	file_with_seq.write((sequence3.join(sequence2)).replace('U','T') + '\n')
    return sequence4



def RNAfold(output_file_name):
    """A function that predicts the secondary structure.
    
    Args:
        output_file_name (str): Variable containing the given file name entered by the user on input.
    
    Returns:
        console_output(str): Variable containing console output.
    """
    fasta_file_name=(output_file_name+'_sequences.fasta')
    console_output = subprocess.run(["RNAfold", "--noPS", fasta_file_name], stdout=subprocess.PIPE).stdout.decode("utf-8")
    return console_output


def RNAfold_record (output_file_name, Thymine_or_uracil):
    """A function that records the intended secondary structure.
    
    Args:
        output_file_name (str): Variable containing the given file name entered by the user on input.
        Thymine_or_uracil (str): Choosing what is in the sequence - Thymine or Uracil.
    """
    RNAfold1=RNAfold(output_file_name)
    RNAfold2=re.sub(r' \(.*?\)', '', RNAfold1)
    file_name=output_file_name+'.fasta'
    with open(file_name, 'w') as fasta_file: 
        if Thymine_or_uracil == 'U':
            fasta_file.write(RNAfold2.replace('T', 'U') + '\n')
        else:
            fasta_file.write(RNAfold2.replace('U', 'T') + '\n')
       



def RNAfold_known_fragments(seq_struc_dot):
    """A function that modifies the dot bracket of known fragments.
    
    Args:
        seq_struc_dot(dict): Dictionary consisting of: id: sequence, dot bracket.
    
    Returns:
        dot_bracket1(list): Modified dot-bracket of known fragments.
    """
    dot_bracket1=[]
    dictionary_values=list(seq_struc_dot.values())
    for m in dictionary_values:
        dot_bracket1.append(m[1].replace (".", "x"))
    return dot_bracket1


def inserting_into_RNAfold(sequences,number,insertion_position,random_to_bed):
        """Function generating a list of modified dot-bracket and dot-bracket random sequences.
        
        Args:
            sequences (list): Sequence list.
            number (int): Number of sequences given by the user as input.
            insertion_position (int): Variable containing the fragment insertion position specified by the user as input.
            random_to_bed (list): Variable containing the positions on which the fragments were inserted.
        
        Returns:
            inserted1(list): A list containing the modified dot-bracket and dot-bracket random sequences.
        """
        fragment_length=[]
        inserted1=[]
        known_fragments=[]
        number_of_fragments=calculating_the_number_of_fragments_in_a_fasta(seq_struc_dot)
        known_fragments=RNAfold_known_fragments(seq_struc_dot)
        for i in range(number):
            sequences[i]=([u.replace ("A", ".").replace("C", ".").replace("U", ".").replace("G", ".").replace("T", ".") for u in sequences[i]])
            sequence_length=len(sequences[i])
            if insertion_position:
                random_place=insertion_position
            else:
                if i<number_of_fragments:
                    fragment_length=len(known_fragments[i])
                    if fragment_length>sequence_length:
                        sys.exit(0)
                    random_place=random_to_bed[i]
                if i>=number_of_fragments:
                    random_place = None
 
            if i<number_of_fragments:
                random_sequences_from_fasta="".join(((sequences[i][:random_place]+list(known_fragments[i]))+sequences[i][random_place:])[:sequence_length])
                inserted1.append(random_sequences_from_fasta)
 
            if i>=number_of_fragments:
                inserted1.append("".join(sequences[i]))
        return inserted1


def writing_RNAfold2(sequence, inserted1, output_file_name, ID,number, Thymine_or_uracil):
    """Write id, sequence and list containing modified dot-bracket.
    
    Args:
        sequences (list): Sequence list.
        inserted1(list): A list containing the modified dot-bracket and dot-bracket random sequences.
        output_file_name (str): Variable containing the given file name entered by the user on input.
        ID(list): List of fasta identifiers with sequence number only.
        number (int): Number of sequences given by the user as input.
        Thymine_or_uracil (str): Choosing what is in the sequence - Thymine or Uracil.
    """
    sequence3=''
    inserted2=''
    fasta_file_name=(output_file_name+'RNAfold2.fasta')
    with open(fasta_file_name, 'w') as rnafold_file:
        number_of_seq=len(sequence)
        for i in range(number):
            sequence2=sequence[i]
            inserted=inserted1[i]
            rnafold_file.write(str(ID[i]) + '\n')
            if Thymine_or_uracil == 'U':
                rnafold_file.write((sequence3.join(sequence2)).replace('T', 'U') + '\n')
            else:
                rnafold_file.write((sequence3.join(sequence2)).replace('U', 'T') + '\n')
            rnafold_file.write(inserted2.join(inserted)+'\n')


def RNAfold2_from_file(output_file_name):
    """A function predicting secondary structure from a file with id, sequence and list containing modified dot-bracket.
    
    Args:
        output_file_name (str): Variable containing the given file name entered by the user on input.
    
    Returns:
        console_output(str): Variable containing console output.
    """
    file_name=(output_file_name+'RNAfold2.fasta')
    console_output1 = subprocess.run(["RNAfold", "-C", "--noPS", "--enforceConstrain", file_name], stdout=subprocess.PIPE).stdout.decode("utf-8")
    return console_output1

def RNAfold_record2 (output_file_name, Thymine_or_uracil):
    """A function that records the intended secondary structure.
    
    Args:
        output_file_name (str): Variable containing the given file name entered by the user on input.
        Thymine_or_uracil (str): Choosing what is in the sequence - Thymine or Uracil.
    
    Returns:
        RNAfold4(list): Variable containing the content of the file with the extension .fasta.
    """
    RNAfold3=RNAfold2_from_file(output_file_name)
    RNAfold4=re.sub(r' \(.*?\)', '', RNAfold3)
    file_name=output_file_name+'.fasta'
    with open(file_name, 'w') as fasta_file2:
        if Thymine_or_uracil == 'U':
            fasta_file2.write(RNAfold4.replace('T', 'U') + '\n')
        else:
            fasta_file2.write(RNAfold4.replace('U', 'T') + '\n')
    return RNAfold4

def RNAfold4_split(RNAfold4):
    """Function splitting the output file into lists containing identifiers, sequences and dot-bracket sequences with known fragments and those without known fragments.
    
    Args:
        RNAfold4(list): Variable containing the content of the file with the extension .fasta.
    
    Returns:
        known1(list): Variable containing the identifier, sequence, and dot-bracket of the sequence with inserted fragments.
        unknown1(list): Variable containing the identifier, sequence, and dot-bracket of the sequence without inserted fragments.
    """
    known=[]
    known1=[]
    unknown=[]
    unknown1=[]
    RNAfold4=RNAfold4.split('\n')
    for line in RNAfold4:
            known.append(line)
            unknown.append(line)
    number_of_lines=len(known)

    for i in range(number_of_lines):
        if '|' in known[i]:
            known1.append(known[i]+'\n')
            known1.append(known[i+1]+'\n')
            known1.append(known[i+2]+'\n')

    for m in range(number_of_lines):
        if not '|' in unknown[m] and 'rand' in unknown[m]:
            unknown1.append(unknown[m]+'\n')
            unknown1.append(unknown[m+1]+'\n')
            unknown1.append(unknown[m+2]+'\n')
    return known1, unknown1


def known_rnafold_record(RNAfold4, output_file_name, Thymine_or_uracil):
    """A function that writes a list containing identifiers, sequences and dot-bracket sequences with known fragments.
    
    Args:
        RNAfold4(list): Variable containing the content of the file with the extension .fasta.
        output_file_name (str): Variable containing the given file name entered by the user on input.
        Thymine_or_uracil (str): Choosing what is in the sequence - Thymine or Uracil.
    """
    known1, unknown1=RNAfold4_split(RNAfold4)
    fasta_file_name=(output_file_name+'_known.fasta')
    with open(fasta_file_name, 'w') as file_known:
        for i in range(len(known1)):
            if Thymine_or_uracil == 'U':
                file_known.write(known1[i].replace('T', 'U'))
            else:
                file_known.write(known1[i].replace('U', 'T'))

def unknown_rnafold_record(RNAfold4, output_file_name,number_of_seq, Thymine_or_uracil):
    """Function that writes a list containing identifiers, sequences, and dot-bracket sequences without known fragments.
    
    Args:
        RNAfold4(list): Variable containing the content of the file with the extension .fasta.
        output_file_name (str): Variable containing the given file name entered by the user on input.
        number_of_seq (int): Number of sequences given by the user as input.
        Thymine_or_uracil (str): Choosing what is in the sequence - Thymine or Uracil.
    """
    known1, unknown1=RNAfold4_split(RNAfold4)
    number_of_known_seq=(len(known1)/3)
    if number_of_seq<=number_of_known_seq:
        sys.stderr.write("The number of sequences is less or equal than the number of fragments. Unable to create random sequence file.")
    else:
        fasta_file_name=(output_file_name+'_random.fasta')
        with open(fasta_file_name, 'w') as file_rand:
            number_of_seq=len(unknown1)
            for i in range(number_of_seq):
                if Thymine_or_uracil == 'U':
                    file_rand.write(unknown1[i].replace('T', 'U'))
                else:
                    file_rand.write(unknown1[i].replace('U', 'T'))


def FASTA_RECORD (sequence, seq_struc_dot, output_file_name, ID,number, Thymine_or_uracil):
    """Function that writes id, sequences and dot-bracket to the fasta file. (without predicting structure)
    
    Args:
        sequences (list): Sequence list.
        seq_struc_dot(dict): Dictionary consisting of: id: sequence, dot bracket.
        output_file_name (str): Variable containing the given file name entered by the user on input.
        ID(list): List of fasta identifiers with sequence number only.
        number (int): Number of sequences given by the user as input.
        Thymine_or_uracil (str): Choosing what is in the sequence - Thymine or Uracil.

    Returns:
    	fasta_file2 (str): Variable containing the contents of the fasta file.
    """
    sequence3=''
    inserted2=''
    fasta_file2=''


    number_of_dot_bracket=len(seq_struc_dot.keys())
    values=list(seq_struc_dot.values())
    fasta_file_name=(output_file_name+'.fasta')
    with open(fasta_file_name, 'w') as fasta_file:
        number_of_seq=len(sequence)
        for i in range(number):
            sequence2=sequence[i]
            if (i<number_of_dot_bracket):
                seq_struc_dot=values[i][1]
            else:
            	seq_struc_dot=''
            fasta_file.write(str(ID[i]) + '\n')
            fasta_file2=fasta_file2+str(ID[i]) + '\n'
            if Thymine_or_uracil == 'U':
                fasta_file.write((sequence3.join(sequence2)).replace('T', 'U') + '\n')
                fasta_file2=fasta_file2+((sequence3.join(sequence2)).replace('T', 'U') + '\n')
            else:
                fasta_file.write((sequence3.join(sequence2)).replace('U', 'T') + '\n')
                fasta_file2=fasta_file2+((sequence3.join(sequence2)).replace('U', 'T') + '\n')
            fasta_file.write(inserted2.join(seq_struc_dot)+'\n')
            fasta_file2=fasta_file2+(inserted2.join(seq_struc_dot)+'\n')
    return fasta_file2


def FASTA_RECORD1 (sequence, output_file_name, ID,number, Thymine_or_uracil):
    """Function that writes id and sequences to the fasta file. (without predicting structure)
    
    Args:
        sequences (list): Sequence list.
        output_file_name (str): Variable containing the given file name entered by the user on input.
        ID(list): List of fasta identifiers with sequence number only.
        number (int): Number of sequences given by the user as input.
        Thymine_or_uracil (str): Choosing what is in the sequence - Thymine or Uracil.
        
    Returns:
    	fasta_file2 (str): Variable containing the contents of the fasta file.
    """
    sequence3=''
    inserted2=''
    fasta_file2=''

    fasta_file_name=(output_file_name+'.fasta')
    with open(fasta_file_name, 'w') as fasta_file:
        number_of_seq=len(sequence)
        for i in range(number):
            sequence2=sequence[i]
            fasta_file.write(str(ID[i]) + '\n')
            fasta_file2=fasta_file2+str(ID[i]) + '\n'
            if Thymine_or_uracil == 'U':
                fasta_file.write((sequence3.join(sequence2)).replace('T', 'U') + '\n')
                fasta_file2=fasta_file2+((sequence3.join(sequence2)).replace('T', 'U') + '\n')
            else:
                fasta_file.write((sequence3.join(sequence2)).replace('U', 'T') + '\n')
                fasta_file2=fasta_file2+((sequence3.join(sequence2)).replace('U', 'T') + '\n')
    return fasta_file2



def known_FASTA_RECORD(fasta_file, output_file_name, Thymine_or_uracil):
    """A function that writes a list containing identifiers, sequences and dot-bracket sequences with known fragments. (without predicting structure)
    
    Args:
        fasta_file(str): Variable containing the content of the file with the extension .fasta.
        output_file_name (str): Variable containing the given file name entered by the user on input.
        Thymine_or_uracil (str): Choosing what is in the sequence - Thymine or Uracil.
    """
    known1, unknown1=RNAfold4_split(fasta_file)
    fasta_file_name=(output_file_name+'_known.fasta')
    with open(fasta_file_name, 'w') as file_known:
        for i in range(len(known1)):
            if Thymine_or_uracil == 'U':
                file_known.write(known1[i].replace('T', 'U'))
            else:
                file_known.write(known1[i].replace('U', 'T'))



def unknown_FASTA_RECORD(fasta_file, output_file_name,number_of_seq, Thymine_or_uracil):
    """Function that writes a list containing identifiers, sequences, and dot-bracket sequences without known fragments. (without predicting structure)
    
    Args:
        fasta_file(str): Variable containing the content of the file with the extension .fasta.
        output_file_name (str): Variable containing the given file name entered by the user on input.
        number_of_seq (int): Number of sequences given by the user as input.
        Thymine_or_uracil (str): Choosing what is in the sequence - Thymine or Uracil.
    """
    known1, unknown1=RNAfold4_split(fasta_file)
    number_of_known_seq=(len(known1)/3)
    if number_of_seq<=number_of_known_seq:
        sys.stderr.write("The number of sequences is less or equal than the number of fragments. Unable to create random sequence file.")
    else:
        fasta_file_name=(output_file_name+'_random.fasta')
        with open(fasta_file_name, 'w') as file_rand:
            number_of_seq=len(unknown1)
            for i in range(number_of_seq):
                if Thymine_or_uracil == 'U':
                    file_rand.write(unknown1[i].replace('T', 'U'))
                else:
                    file_rand.write(unknown1[i].replace('U', 'T'))




parser = argparse.ArgumentParser(prog='PROG', description='PROJEKT')
group = parser.add_mutually_exclusive_group()
group.add_argument('-e', dest='length_range', type=str, help='The range including the length of the random sequences.')
group.add_argument('-l', dest='length', type=int, help='The length of the generated random sequences.')
group2 = parser.add_mutually_exclusive_group()
group2.add_argument('-k', dest='percentage_of_known_fragments', type=int, help='Percentage of known fragments in series.')
group2.add_argument('-n', dest='number', type=int, help='Number of random sequences generated.')
parser.add_argument('-p', dest='path', type=str, help='Fasta file path.')
parser.add_argument('-r', dest='reshuffle', type=int, choices=[2, 3], help='Reshuffle for dinucleotides(2) or trinucleotides(3)?')
parser.add_argument('-g', dest='gc', type=int, help='The percent of GC nucleotides in a random sequence.')
parser.add_argument('-o', dest='output_file_name', type=str, help='The name of the output file.')
parser.add_argument('-i', dest='insertion_position', type=int, help='Position in which known fragment is inserted into a random sequence.')
parser.add_argument('-a', dest='additional_fasta_file', action='store_true', help='Additional fasta file.')
parser.add_argument('-s', dest='separate', action='store_true', help='Separate fasta files with sequences and structures for those fully random and those with given fragments.')
parser.add_argument('-b', dest='file_with_sequences_to_shuffle',  type=str, help='Fasta file path with sequences to shuffle')
parser.add_argument('-t', dest='Thymine_or_uracil',  type=str, help='Should the sequences in the output files contain thymine (T) or Uracil (U)?')
parser.add_argument('-w', dest='without_predicting',  action='store_true', help='Option without predicting the structure')

args = parser.parse_args()

if args.path and not os.path.exists(args.path):
    sys.stderr.write("Could not open fasta file '{}': No such file\n".format(args.path))
    sys.exit(0)

if not args.without_predicting:
#################################################################################################################### with file
	if args.path and os.path.exists(args.path):

	    if args.length_range or args.length:
	        if args.length_range:
	            length_range=args.length_range
	        if args.length:
	            sequence_length=args.length
	    else:
	        sequence_length=200


	    nucleotide=checking_fasta(args.path)
	    if args.Thymine_or_uracil:
		    Thymine_or_uracil=args.Thymine_or_uracil
	    else:
		    Thymine_or_uracil= nucleotide


	    seq_struc_dot = reading_fasta_files(args.path)

	    if args.number:
	        number=args.number
	    elif args.percentage_of_known_fragments:
	        number=calculating_the_total_number_of_fragments(args.percentage_of_known_fragments,seq_struc_dot)
	    else:
	        number=10

	    GC_content = args.gc if args.gc else 50

	    if args.file_with_sequences_to_shuffle:
	        sequences_to_shuffle, shuffle_nucleotide =loading_an_additional_file(args.file_with_sequences_to_shuffle)
	        if args.gc:
	            sys.stderr.write("The gc percentage parameter is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.length:
	           sys.stderr.write("The length is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.length_range:
	           sys.stderr.write("The range including the length of the random sequences is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.reshuffle:
	            sequences1= sequence_shuffling(number, sequences_to_shuffle, args.reshuffle)
	        else:
	            args.reshuffle == 2
	            sequences1= sequence_shuffling(number, sequences_to_shuffle, args.reshuffle)



	    if not args.file_with_sequences_to_shuffle:
	        if not args.length_range:
	           sequences1= number_of_sequences(number, sequence_length)
	        else:
	           sequences1 = number_of_sequences__RANGE(number, length_range)



	    if args.reshuffle:
	        if not args.file_with_sequences_to_shuffle:
	            sys.stderr.write('If you specify the way you want to shuffle sequences, you must supply that sequence.')
	            sys.exit(0)

	    random_to_bed, fragment_length_to_bed, number_of_fragments, inserted=inserting_known_sequences(sequences1,number,args.insertion_position)
	    BED1 = BED_format(number, random_to_bed, fragment_length_to_bed, number_of_fragments, args.insertion_position)
	    ID = ID_record(number, random_to_bed, fragment_length_to_bed, number_of_fragments)
	    inserted1=inserting_into_RNAfold(sequences1,number,args.insertion_position, random_to_bed)



	    if args.output_file_name:
	        output_file_name1 = args.output_file_name
	        writing_to_txt_file1= writing_to_txt_file(inserted, output_file_name1, Thymine_or_uracil)
	        BED_record= writing_to_BED_file(BED1, output_file_name1)
	        sequence4= writing_to_FASTA_file(inserted, output_file_name1,ID,number, Thymine_or_uracil)
	        writing_RNAfold2=writing_RNAfold2(inserted, inserted1, output_file_name1,ID,number,Thymine_or_uracil)
	        RNAfold4=RNAfold_record2(output_file_name1,Thymine_or_uracil)
	        if args.separate:
	            unknown_rnafold_record=unknown_rnafold_record(RNAfold4, output_file_name1,number, Thymine_or_uracil)
	            known_rnafold_record=known_rnafold_record(RNAfold4, output_file_name1, Thymine_or_uracil)
	            os.remove(output_file_name1+'.fasta')
	        os.remove(output_file_name1+'.txt')
	        os.remove(output_file_name1+'RNAfold2.fasta')
	        if not args.additional_fasta_file:
	            os.remove(output_file_name1+'_sequences.fasta')

	       

	    else:
	        output_file_name1 = 'name'
	        writing_to_txt_file1 = writing_to_txt_file(inserted, output_file_name1, Thymine_or_uracil)
	        BED_record= writing_to_BED_file(BED1, output_file_name1)
	        sequence4=writing_to_FASTA_file(inserted, output_file_name1,ID,number, Thymine_or_uracil)
	        writing_RNAfold2=writing_RNAfold2(inserted, inserted1, output_file_name1,ID,number,Thymine_or_uracil)
	        RNAfold4=RNAfold_record2(output_file_name1,Thymine_or_uracil)
	        if args.separate:
	            unknown_rnafold_record=unknown_rnafold_record(RNAfold4, output_file_name1,number, Thymine_or_uracil)
	            known_rnafold_record=known_rnafold_record(RNAfold4, output_file_name1, Thymine_or_uracil)
	            os.remove(output_file_name1+'.fasta')
	        os.remove(output_file_name1+'.txt')
	        os.remove(output_file_name1+'RNAfold2.fasta')
	        if not args.additional_fasta_file:
	            os.remove(output_file_name1+'_sequences.fasta')

	 
	############################################################################################################################ without file
	if args.path is None:

	    if args.length_range or args.length:
	        if args.length_range:
	            length_range=args.length_range
	        if args.length:
	            sequence_length=args.length
	    else:
	        sequence_length=200


	    if args.separate:
	        sys.stderr.write('You cannot create separate fasta files with sequences and structures for those fully random and those with given fragments without giving a known fragment file.')
	        sys.exit(0)

	    if args.number:
	        number=args.number
	    elif args.percentage_of_known_fragments:
	        sys.stderr.write("You cannot give percentage of known fragments without giving a known fragment file.")
	        sys.exit(0)
	    else:
	        number=10


	    GC_content = args.gc if args.gc else 50

	    shuffle_nucleotide = None
	    if args.file_with_sequences_to_shuffle:
	        sequences_to_shuffle, shuffle_nucleotide=loading_an_additional_file(args.file_with_sequences_to_shuffle)
	        if args.gc:
	            sys.stderr.write("The gc percentage parameter is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.length:
	           sys.stderr.write("The length is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.length_range:
	           sys.stderr.write("The range including the length of the random sequences is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.reshuffle:
	            sequences1=sequence_shuffling(number, sequences_to_shuffle, args.reshuffle)

	        if not args.reshuffle:
	            args.reshuffle == 2
	            sequences1=sequence_shuffling(number, sequences_to_shuffle, args.reshuffle)


	    if args.Thymine_or_uracil:
	        Thymine_or_uracil = args.Thymine_or_uracil
	    elif shuffle_nucleotide is not None:
	        Thymine_or_uracil = shuffle_nucleotide
	    else:
	        Thymine_or_uracil = 'U'



	    if not args.file_with_sequences_to_shuffle:
	        if not args.length_range:
	           sequences1 = number_of_sequences(number, sequence_length)
	        else:
	           sequences1 = number_of_sequences__RANGE(number, length_range)


	    if args.reshuffle:
	        if not args.file_with_sequences_to_shuffle:
	            sys.stderr.write('If you specify the way you want to shuffle sequences, you must supply that sequence.')
	            sys.exit(0)


	    ID= ID_record_without_file(number)

	    if args.insertion_position:
	        sys.stderr.write("You cannot give an insertion position for known fragments without giving a file with known fragments.")
	        sys.exit(0)

	    if args.output_file_name:
	        output_file_name1 = args.output_file_name
	        writing_to_txt_file1 = writing_to_txt_file(sequences1, output_file_name1, Thymine_or_uracil)
	        sequence4=writing_to_FASTA_file(sequences1, output_file_name1,ID,number, Thymine_or_uracil)
	        RNAfold_record=RNAfold_record(output_file_name1, Thymine_or_uracil)
	        os.remove(output_file_name1+'.txt')
	        if not args.additional_fasta_file:
	            os.remove(output_file_name1+'_sequences.fasta')



	    else:
	        output_file_name1 = 'name'
	        writing_to_txt_file1 = writing_to_txt_file(sequences1, output_file_name1, Thymine_or_uracil)
	        sequence4=writing_to_FASTA_file(sequences1, output_file_name1,ID,number, Thymine_or_uracil)
	        RNAfold_record=RNAfold_record(output_file_name1,Thymine_or_uracil)
	        os.remove(output_file_name1+'.txt')
	        if not args.additional_fasta_file:
	            os.remove(output_file_name1+'_sequences.fasta')

















############################################################################################################################
############################################################################################################################ without predicting structure

if args.without_predicting:
     
#################################################################################################################### with file
	if args.path and os.path.exists(args.path):
	    if args.length_range or args.length:
	        if args.length_range:
	            length_range=args.length_range
	        if args.length:
	            sequence_length=args.length
	    else:
	        sequence_length=200


	    nucleotide=checking_fasta(args.path)
	    if args.Thymine_or_uracil:
		    Thymine_or_uracil=args.Thymine_or_uracil
	    else:
		    Thymine_or_uracil= nucleotide


	    seq_struc_dot = reading_fasta_files(args.path)

	    if args.number:
	        number=args.number
	    elif args.percentage_of_known_fragments:
	        number=calculating_the_total_number_of_fragments(args.percentage_of_known_fragments,seq_struc_dot)
	    else:
	        number=10

	    GC_content = args.gc if args.gc else 50

	    if args.file_with_sequences_to_shuffle:
	        sequences_to_shuffle, shuffle_nucleotide =loading_an_additional_file(args.file_with_sequences_to_shuffle)
	        if args.gc:
	            sys.stderr.write("The gc percentage parameter is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.length:
	           sys.stderr.write("The length is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.length_range:
	           sys.stderr.write("The range including the length of the random sequences is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.reshuffle:
	            sequences1= sequence_shuffling(number, sequences_to_shuffle, args.reshuffle)
	        else:
	            args.reshuffle == 2
	            sequences1= sequence_shuffling(number, sequences_to_shuffle, args.reshuffle)



	    if not args.file_with_sequences_to_shuffle:
	        if not args.length_range:
	           sequences1= number_of_sequences(number, sequence_length)
	        else:
	           sequences1 = number_of_sequences__RANGE(number, length_range)



	    if args.reshuffle:
	        if not args.file_with_sequences_to_shuffle:
	            sys.stderr.write('If you specify the way you want to shuffle sequences, you must supply that sequence.')
	            sys.exit(0)

	    random_to_bed, fragment_length_to_bed, number_of_fragments, inserted=inserting_known_sequences(sequences1,number,args.insertion_position)
	    BED1 = BED_format(number, random_to_bed, fragment_length_to_bed, number_of_fragments, args.insertion_position)
	    ID = ID_record(number, random_to_bed, fragment_length_to_bed, number_of_fragments)
	    inserted1=inserting_into_RNAfold(sequences1,number,args.insertion_position, random_to_bed)



	    if args.output_file_name:
	        output_file_name1 = args.output_file_name
	        writing_to_txt_file1= writing_to_txt_file(inserted, output_file_name1, Thymine_or_uracil)
	        BED_record= writing_to_BED_file(BED1, output_file_name1)
	        sequence4= writing_to_FASTA_file(inserted, output_file_name1,ID,number, Thymine_or_uracil)
	        fasta_file=FASTA_RECORD(inserted, seq_struc_dot, output_file_name1,ID,number,Thymine_or_uracil)
	        if args.separate:
	            os.remove(output_file_name1+'.fasta')
	            unknown_FASTA_RECORD=unknown_FASTA_RECORD(fasta_file, output_file_name1,number, Thymine_or_uracil)
	            known_FASTA_RECORD=known_FASTA_RECORD(fasta_file, output_file_name1, Thymine_or_uracil)
	        os.remove(output_file_name1+'.txt')
	        if not args.additional_fasta_file:
	           os.remove(output_file_name1+'_sequences.fasta')

	       

	    else:
	        output_file_name1 = 'name'
	        writing_to_txt_file1 = writing_to_txt_file(inserted, output_file_name1, Thymine_or_uracil)
	        BED_record= writing_to_BED_file(BED1, output_file_name1)
	        sequence4=writing_to_FASTA_file(inserted, output_file_name1,ID,number, Thymine_or_uracil)
	        fasta_file=FASTA_RECORD(inserted, seq_struc_dot, output_file_name1,ID,number,Thymine_or_uracil)
	        if args.separate:
	            os.remove(output_file_name1+'.fasta')
	            unknown_FASTA_RECORD=unknown_FASTA_RECORD(fasta_file, output_file_name1,number, Thymine_or_uracil)
	            known_FASTA_RECORD=known_FASTA_RECORD(fasta_file, output_file_name1, Thymine_or_uracil)
	        os.remove(output_file_name1+'.txt')
	        if not args.additional_fasta_file:
	            os.remove(output_file_name1+'_sequences.fasta')

	 
	############################################################################################################################ without file
	if args.path is None:

	    if args.length_range or args.length:
	        if args.length_range:
	            length_range=args.length_range
	        if args.length:
	            sequence_length=args.length
	    else:
	        sequence_length=200


	    if args.separate:
	        sys.stderr.write('You cannot create separate fasta files with sequences and structures for those fully random and those with given fragments without giving a known fragment file.')
	        sys.exit(0)

	    if args.number:
	        number=args.number
	    elif args.percentage_of_known_fragments:
	        sys.stderr.write("You cannot give percentage of known fragments without giving a known fragment file.")
	        sys.exit(0)
	    else:
	        number=10


	    GC_content = args.gc if args.gc else 50

	    shuffle_nucleotide = None
	    if args.file_with_sequences_to_shuffle:
	        sequences_to_shuffle, shuffle_nucleotide=loading_an_additional_file(args.file_with_sequences_to_shuffle)
	        if args.gc:
	            sys.stderr.write("The gc percentage parameter is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.length:
	           sys.stderr.write("The length is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.length_range:
	           sys.stderr.write("The range including the length of the random sequences is not available in sequence shuffle mode. The specified value was omitted.")

	        if args.reshuffle:
	            sequences1=sequence_shuffling(number, sequences_to_shuffle, args.reshuffle)

	        if not args.reshuffle:
	            args.reshuffle == 2
	            sequences1=sequence_shuffling(number, sequences_to_shuffle, args.reshuffle)


	    if args.Thymine_or_uracil:
	        Thymine_or_uracil = args.Thymine_or_uracil
	    elif shuffle_nucleotide is not None:
	        Thymine_or_uracil = shuffle_nucleotide
	    else:
	        Thymine_or_uracil = 'U'



	    if not args.file_with_sequences_to_shuffle:
	        if not args.length_range:
	           sequences1 = number_of_sequences(number, sequence_length)
	        else:
	           sequences1 = number_of_sequences__RANGE(number, length_range)


	    if args.reshuffle:
	        if not args.file_with_sequences_to_shuffle:
	            sys.stderr.write('If you specify the way you want to shuffle sequences, you must supply that sequence.')
	            sys.exit(0)


	    ID= ID_record_without_file(number)

	    if args.insertion_position:
	        sys.stderr.write("You cannot give an insertion position for known fragments without giving a file with known fragments.")
	        sys.exit(0)

	    if args.output_file_name:
	        output_file_name1 = args.output_file_name
	        writing_to_txt_file1 = writing_to_txt_file(sequences1, output_file_name1, Thymine_or_uracil)
	        sequence4=writing_to_FASTA_file(sequences1, output_file_name1,ID,number, Thymine_or_uracil)
	        fasta_file=FASTA_RECORD1(sequences1, output_file_name1,ID,number,Thymine_or_uracil)
	        os.remove(output_file_name1+'.txt')
	        if not args.additional_fasta_file:
	            os.remove(output_file_name1+'_sequences.fasta')



	    else:
	        output_file_name1 = 'name'
	        writing_to_txt_file1 = writing_to_txt_file(sequences1, output_file_name1, Thymine_or_uracil)
	        sequence4=writing_to_FASTA_file(sequences1, output_file_name1,ID,number, Thymine_or_uracil)
	        fasta_file=FASTA_RECORD1(sequences1, output_file_name1,ID,number,Thymine_or_uracil)
	        os.remove(output_file_name1+'.txt')
	        if not args.additional_fasta_file:
	            os.remove(output_file_name1+'_sequences.fasta')

