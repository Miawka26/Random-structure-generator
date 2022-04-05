# Random structure generator 
## Requirements:
Required python version: Python>=3.5
Libraries (library name/command):
* Random (import random),
* Os (import os),
* Argparse (import argparse),
* Sys (import sys),
* Subprocess (import subproces),
* Re (import re),
* ushuffle -> Shuffle, shuffler (from ushuffle import shuffle, Shuffler).

You will also need the RNA structural analysis package, Vienna RNA. It can be downloaded from https://www.tbi.univie.ac.at/RNA/ where you just need to select your operating system version and then follow the installation instructions.

## Program options:
* Length of generated fragments:
    - -e (optional) - Range in which sequence lengths are to be (string)
   - -l (optional) - Sequence length (int)
* Number of random sequences generated:
   - -k (optional) - Percentage to be sequences with known fragments of the total set of sequences (int)
    - -n (optional) - Number of generated sequences (int)
*  -p (optional) - Path to the file with the extension '.fasta' containing the known fragments (string)
*  -r (optional) - The option to choose whether the shuffle will be for dinucleotides or trinucleotides (int)
*  -g (optional) - This option allows you to specify the percentage of gc in the sequence (int)
*  -o (optional) - Output file name (string)
*  -i (optional) - Position where the known fragment should be inserted (int)
*  -a (optional) - Option to create an additional fasta file (store_true)
*  -s (optional) - The option to get two separate files with sequences and structures: fully random and those with the specified fragments (store_true)
*  -b (optional) - Path to the file with the extension '.fasta' containing the sequence to be reshuffled (string)


## Input files - description of formats:
1. A non-empty file with extension '.fasta' containing known fragments to be inserted into a random or reshuffled sequence. It should have at least one: identifier (starting with '>'), RNA sequence (containing 'A', 'U', 'C', 'G') and dot-bracket (containing '(', ')', '.'). 
2. A non-empty file with the extension '.fasta' containing the RNA sequence (containing 'A', 'U', 'C', 'G') to be shuffled.

## Output files - description of formats:
1. name*.fasta - File with extension '.fasta' containing identifiers, sequences and structures. 
2. name*.bed - File with extension '.bed' containing coordinates of inserted fragments. 
3. name*_sequences.fasta - File with extension '.fasta' containing only sequences with identifiers. 
4. name*_known.fasta, name*_random.fasta - Files with extension '.fasta' containing identifiers, sequences and structures: fully random/random sequences and those with specified fragments. 
*By default, the program saves under the name "name", but the user has the possibility to give his own name.

## Program modes:
### 1. Create random sequences. 
The user can create random sequences without entering any additional files.
He can then specify: what should be the number (-n) of random sequences, length (-e/-l) of the sequence, the percentage of GC in the sequence (-g), a declaration to create an additional fasta file (-a) and the file name (-o) under which the resulting files should be saved. If you enter other options, the program will display an error or another program mode will be executed.


#### Workflow
1. Draw without returning nucleic bases from a list containing calculated amounts of guanine, cytosine, adenine, uracil based on GC content.
2. Adding in a loop (of sequence length) the drawn nucleic bases to the list. 
3. Creating the appropriate number of sequences based on the default value or given by the user. 
4. Predicting the structure with the RNAfold program. 
5. Save the sequence and predicted structure to the file 'name.fasta' or to a file with the extension '.fasta' specified by the user. 
6. If the user chooses to create an additional fasta file, it will create a file 'name_sequences.fasta' (or instead of name - the name given on the input by the name given by the user on input) containing only sequences with identifiers.

### 2. Inserting known fragments into a random sequence. 
The user can insert known fragments into a random sequence. He must then specify a file with the extension '.fasta' containing the known fragments (-p). He then also has the option to enter: what should be the number (-n) of sequences or the percentage that sequences (of the given length) with known fragments should be in relation to the whole set of sequences (-k), the chosen length (-e/-l) of the sequence, the percentage of GC in the sequence (-g), the position at which the known fragment should be inserted (-i), a declaration to create an additional fasta file (-a), a declaration to create two separate files with sequences and structures: fully random and those with the specified fragments (-s), and the file name (-o) under which the resulting files are to be saved. If you enter any other options the program will display an error or a different program mode will be executed.

#### Workflow
1. Load a file with known fragments. 
2. Check the loaded fasta file to ensure that it meets the following conditions:
    * is not empty, 
    * contains an identifier (starting with '>'), 
    * contains an RNA sequence (containing 'A', 'U', 'C', 'G') and no invalid character, 
    * contains a dot-bracket (containing '(', ')', '.') and there is no invalid character, 
    * the order: identifier, sequence, dot-bracket is preserved. 
3. Writing values from a file to a dictionary containing the schema: identifier: sequence, dot bracket. 
4. Drawing without returning nucleic bases from a list containing calculated amounts of guanine, cytosine, adenine, uracil based on GC content. 
5. Adding in a loop (of sequence length) the drawn nucleic bases to the list. 
6. Creating the appropriate number of sequences based on a default value or specified by the user. 
7. Inserting known fragments into random sequences at a random position or at a position specified by the user.
8. Creating a bed format. 
9. Preparing the structure for prediction with the RNAfold program (in random sequences, nucleotides replace '.' . On the other hand in the structures of known fragments change '.' to 'x', the rest of the structure remains unchanged). 
10. Predict the structure with the RNAfold program. 
11. Save the sequence and the predicted structure to the file 'name.fasta' or to a file with a user specified name with '.fasta' extension.
12. Saving the coordinates of the inserted fragments in the file name.bed or to the file with the extension '.bed' specified by the user.
13. If the user chooses to create an additional fasta file, a 'name_sequences.fasta' file will be created (or instead of name - the name given on input by the user) containing only sequences with IDs. 
14. If the user chooses to create two separate files with sequences and structures: fully random and those with specified fragments, the files 'name_known.fasta' and 'name_random.fasta' will be created (or instead of name - the name given on input by the user).


### 3. Shuffling of the given fragment. 
The user has the option to shuffle the given sequence. He then has to specify a file with the extension '.fasta' containing the sequences to be reshuffled (-b). He can also specify the number (-n) of sequences, a declaration of creating an additional fasta file (-a), a declaration of how the sequence is to be reshuffled (-r), and the name file (-o) under which the resulting files are to be saved. Entering any other option will cause the program to display an error or a different program mode will be executed.

#### Workflow
1. Load a sequence file for shuffling. 
2. Checking whether the specified file with the sequence to be shuffled is empty. 
3. Shuffling of the sequence, taking into account the user's choice, whether it should be for dinucleotides or trinucleotides. 
4. Predict the structure with the RNAfold program. 
5. Save the sequence and predicted structure to the file 'name.fasta' or to a file with the extension '.fasta' specified by the user. 
6. If the user chooses to create an additional fasta file, it will create a file 'name_sequences.fasta' (or instead of name - the name given on the input by the user input) containing only sequences with identifiers.

### 4. Inserting known fragments into a shuffled sequence
User can also specify the number (-n) of sequences or the percentage of sequences with known fragments in relation to the whole set of sequences (-k), the position at which a known fragment is to be inserted (-i), a declaration to create an additional fasta file (-a), a declaration to create two separate files with sequences and structures: fully random and the one with the specified fragments (-s), a declaration on how to reshuffled sequence (-r) and the file name (-o) under which the files to be saved. If other options are specified, an error is generated or a different program mode is executed.



#### Workflow
1. Load a fasta file with fragments. 
2. Checking the loaded fasta file to ensure that it meets the following conditions:
    * is not empty, 
    * contains an identifier (starting with '>'), 
    * contains an RNA sequence (containing 'A', 'U', 'C', 'G') and no invalid character, 
    * contains a dot-bracket (containing '(', ')', '.') and no invalid character is present,
    * the order: identifier, sequence, dot-bracket is preserved. 
3. Writing values from a file to a dictionary containing the schema: identifier: sequence, dot bracket. 
4. Loading the sequence file for shuffling. 
5. Checking that the specified file with the sequence to be shuffled is not empty. 
6. Shuffle the sequence. 
7. Inserting known fragments into random sequences at random position or at position specified by the user. 
8. Creating a bed format. 
9. Preparation of the structure to be predicted with the RNAfold program (in the reshuffled sequences nucleotides replace '.' . In the structures of known fragments '.' are changed to 'x', the rest of the structure remains unchanged). 
10. Predict the structure with the RNAfold program. 
11. Save the sequence and the predicted structure to the file 'name.fasta' or to a file with a user defined name with extension '.fasta'. 
12. Saving the coordinates of the inserted fragments in the file name.bed or to the file with the extension '.bed' specified by the user. 
13. If the user chooses to create an additional fasta file, the following file will be created 'name_sequences.fasta' is created (or instead of name - the name given on input by the user) containing only sequences with identifiers. 
14. If the user chooses to create two separate files with sequences and structures: fully random and those with specified fragments the files 'name_known.fasta' and 'name_random.fasta' will be created (or instead of name - the name given on input by the user).

### Choosing a method to create random sequences
When choosing a method, the following options were considered: randomization without returning, drawing with return and drawing with the help of the Nullseq program. Due to the very similar results of the obtained percentage compared to the desired GC percentage, a returnless random draw was chosen to create random sequences. This method also did not produce any error. The draw using the Nullseq program produced an error when the GC percentage was set to values of 10%, 30%, 90%.

<img src="images/Test_methods_for_creating_random_sequences.png">
<img src="images/Test_methods_for_creating_random_sequences_table.png">


### Shuffle validation
In order to test the correctness of the reshuffle, it was checked how the dinucleotides/trinucleotides in the sequence to be reshuffled were distributed versus the starting sequences (the average of 100 sequences was used here). The results obtained were then compared and graphs were created from them. 
This test was performed for the non-overlapping dinucleotides/trinucleotides by taking the sequences and splitting them from the beginning into dinucleotides/trinucleotides, as well as for the overlapping dinucleotides/trinucleotides (splitting the sequences from position 0 and 1). The results show that the sequence shuffling function works correctly, as the frequencies in the sequence to be shuffled as well as the outgoing sequences are preserved.


<img src="images/non_overlapping_dinucleotides.png">
<img src="images/non_overlapping_trinucleotides.png">
<img src="images/overlapping_dinucleotides.png">
<img src="images/overlapping_trinucleotides.png">
