


Friday, 3/8/2024:

#ToDo encapsulate the methods into a class to make the method calls less verbose

Monday, 3/4/2024

Propective Algorithm:

1. Get the first 40 nucleotides of the sequence.
2. get the minimum and maximum melting temperatures of each length, 20 - 40
3. Shift the 40 long window over by 1 nucleotide?

Sunday, 3/3/2024

I have the following function: generate_melting_temp_dict_from_sequence()

My question is, why do I want to generate a melting temp dict from a sequence in the first place.

So basically,

for every overlap end:
   for every overlap start:
      calculate min and max melting temperature?
      

What do I want?

For every given temperature and starting? overlap, what are all the possible ending overlaps.

Thursday:

1. Compute_temp_dict
2. calculate min_temp_dict and max_temp_dict for every position for lengths 20 - 40
3. filter global_temp_dict using min and max temp_dict for every position.
4. Now we have allowable overlaps for every position and temperature

One limitation: Overlap sequences are generated that have CG content that is well below or well above the 
standard cg amount range.

THere are many temperature dictionaries, and I need to come up with a standard naming convention so that I do
not confuse them all.

1. global_temperature_cg_overlap_dict: temp dict for all deisred overlaps and all desired cg_content
2. positional_min_temperature_cg_overlap_dict:
3. positional_max_temperature_cg_overlap_dict
4. positional_filtered_temperature_cg_overlap_dict:
5. temperature_cg_overlap_dict_array: Array of positional_filtered_temperature_cg_overlap_dicts



Then, there will be an array of positional_filtered_temperature_dict, each position will have one dictionary. 

Tuesday, 2/6/2024

I finished function that determines minimum and maximum cg_amount for any position yesterday. However, I would rather that it \



Monday, 2/5/2024

Yesterday, I successfully tested the compute_melting_temp_dict function. Today, I will make the next function in the sequence.
This function takes in as input the melting_temp_dict, as well as the minimum and maximum melting temperature binary sequences.
The function will filter the values from the melting_temp_dict that are within bounds of the minimum and maximum allowable
cg_content. 

However, this function is a bit too complex. Really, the binary format of the minimum and maximum dna sequences is very
unwieldy. I would like to refactor those functions to instead output a list of numbers that correspond to the total number
of cg found so far.

For example. If the function got "GTCCGAT" as input. It would give [4, 3, 3, 2, 1, 0, 0] as output. This makes the minimum 
and maximum cg_amount for any given position more readily available, as it does not have to be calculated from a binary 
dna sequence at the time the information is needed.

Sunday, 2/4/2024

Today, I will try testing the dna melting temp dict.

Thursday, 2/1/2024

Melting temperature as a function of overlap length and cg_content with overlap on the x-axis and
cg_content on the y-axis.

The first 40 nucleotides of the reverse translation DNA from the ferrodoxin_oxidoreduxase.fasta file were used. <br>
low_cg_dna: "0010101000100000101000001001001000000000" <br>
high_cg_dna: "0010111010110011111011011011011010011010" <br>
      40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20  <br>
1  [  __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ ]<br>
2  [  __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ ]<br>
3  [  __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ ]<br>
4  [  __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ 38 37 ]<br>
5  [  __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ 42 41 40 39 ]<br>
6  [  __ __ __ __ __ __ __ __ __ __ __ 48 48 47 46 46 45 44 43 42 41 ]<br>
7  [  __ __ __ __ __ __ __ 51 51 51 50 50 49 49 48 47 47 46 45 44 43 ]<br>
8  [  __ __ __ __ __ 53 53 53 52 52 52 51 51 50 50 49 48 48 47 46 45 ]<br>
9  [  __ __ __ 55 55 55 54 54 54 53 53 53 52 52 51 51 50 49 49 48 47 ]<br>
10 [  57 57 56 56 56 56 55 55 55 55 54 54 54 53 53 52 52 51 51 50 49 ]<br>
11 [  58 58 57 57 57 57 57 56 56 56 56 55 55 55 54 54 53 53 52 52 51 ]<br>
12 [  59 59 59 58 58 58 58 58 57 57 57 57 56 56 56 56 55 55 54 54 53 ]<br>
13 [  60 60 60 60 59 59 59 59 59 59 58 58 58 58 57 57 57 57 56 __ __ ]<br>
14 [  61 61 61 61 61 60 60 60 60 60 60 60 59 59 59 59 59 58 __ __ __ ]<br>
15 [  62 62 62 62 62 62 62 61 61 61 61 61 61 61 61 60 60 __ __ __ __ ]<br>
16 [  63 63 63 63 63 63 63 63 63 63 62 62 62 62 62 62 __ __ __ __ __ ]<br>
17 [  64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 __ __ __ __ __ __ ]<br>
18 [  65 65 65 65 65 65 65 65 65 65 65 65 __ __ __ __ __ __ __ __ __ ]<br>
19 [  66 66 66 66 66 66 66 66 66 67 67 __ __ __ __ __ __ __ __ __ __ ]<br>
20 [  67 67 67 67 67 67 68 68 68 __ __ __ __ __ __ __ __ __ __ __ __ ]<br>
21 [  68 68 68 68 69 69 69 __ __ __ __ __ __ __ __ __ __ __ __ __ __ ]<br>
22 [  69 69 69 69 70 70 __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ ]<br>
23 [  70 70 70 71 71 __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ ]<br>
24 [  71 71 72 __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ __ ]<br>


Here is the plan of attack:

melting_temp_dict[melting_temperature] = [(overlap1, cg_content1), (overlap2, cg_congent2), ... , (overlap-n, cg_conent-n)]

use the melting temperature dict along with minimuma nd maximum cg_overlap information to determine overlap_sets for each


Edit 3/10/2024:
The problem with this plan of attack is that not all temperatures will be covered by minimum and maximum cg_overlap information.
What if, instead, I calculate an overlap dict, that stores (melting_temp, cg_content) tuples, in order from lowest cg_content
to highest cg_content. Then, I can filter this array, dict, only taking the entries in the range (min_cg_content, max_cg_content)
for each overlap length.

completed: function: generate_global_overlap_dict(min_overlap, max_overlap, min_cg, max_cg) -> global_overlap_dict[overlap] = [(min_cg, temp_1), (min_cg + 1, temp_2), ... , (max_cg, temp_n)]
function: generate_overlap_dict_from_sequence(binary_sequence, min_overlap, max_overlap) -> sequence_overlap_dict[overlap] = (cg_content)
function: generate_overlap_dict_array(binary_sequence, ....)
min_temp_overlap_dict = generate_overlap_dict_from_sequence(min_binary_sequence, min_overlap, max_overlap)
max_temp_overlap_dict = generate_overlap_dict_from_sequence(max_binary_sequence, min_overlap, max_overlap)
Friday, 1/26/2024

Hello, long time no update. I have began implementing the overlap calculator. Over the last three days, I have gone
down the rabbit hole for calculating DNA melting temperature. There are two primary ways that I could use, There is
the basic approach, which is 64.9 + 41 * (cg_amount - 16.4) / dna_length. Then, there is the nearest neighbor approach,
that calculates the temperature based on neighboring pairs. There is also the salt adjusted model, but that model only
becomes more accurate for longer overlap lengths than are going to be used in this program. For now, I will use the basic
approach. I will talk about the limitations of the approach, and what would be needed to adjust this algorithm to the 
nearest neighbor temperature model. First off, I would have to write the new melting temperature calculator. That should
not be too bad, and I think would only take a few days. But then how do you determine the maximum and minimum overlap range
for a given temperature and DNA index. Now, since temperature is entirely dependent on cg_amount, I just precompute 
the DNA sequence with the minimal and maximal cg_amount, and plug those sequences into the overlap calculator to determine
the minimum and maximum overlap for every temperature and starting index. However, I do not know how to calculate such a
minimum and maximum sequence for the nearest neighbor model. I believe such a solution would require dynamic programming.
I think I would be interested in tackling it once the rest of this program is complete.

For now, I will stick with the basic algorithm. The problem I am running into now, is as follows. I have calculated the
minimum and maximum cg_content DNA strings, and have used them to compute the minimum and maximum allowable overlap matrix.
At the moment, however, there are some position temperature coordinates that have no overlap length that melts within
the desired temperature range. What I want to do in these cases is to repeatedly increment (or decrement) the cg count
until I find a .... Apologies for stopping at the middle of a thought. I will run some tests to see if I can get a better
understanding of the problem, and then get back to you.



Tuesday, 1/23/2024

Over Thursday, Friday, and Saturday I made a lot of progress. I now can fully convert a protein-sequence to a
dna-sequence and vis-versa. It is now time to begin implementing the actual algorithm portion of the program.
I am very nervous for this part, and we will see how it goes.



Thursday, 1/18/2024

Daily Debrief:
Yesterday, I completed the codon_apportionment function. Today, I will hopefully complete 
the rest of the protein_sequence_to_dna_sequence function

Wednesday, 1/17/2024

Daily Debrief:
Yesterday, I finished refactoring data_parser.py. Today, I am working on the conversion functions between protein
sequence and DNA sequence. First, I need to write a function that combines protein_data and multi_amino data into 
one data structure. i.e.: ({Y}, {A}, {'K', 'X'} ... {'I'}, {'G', 'N', 'P', 'S'})

Update: Before I write any further code, I am trying to wrap my head around the algorithm I used last time around, 
in my java implementation. Here are my findings:

1. altAminoToDecodons(multi_amino) &rarr; multi_codon_data
<br>&emsp; converts list of multi_aminos into a list of multi_codons
2. proteinSequenceToDna(protein_sequence) &rarr; dna_sequence
3. proteinSequenceToBinary(protein_sequence, mode : (min or max)) &rarr; dna_binary_sequence
<br>&emsp; output is the corresponding DNA sequence with CG content either maximised or minimised. The sequence is in
binary with a 0 : (A or T) and 1 : (C or G)
4. overlapCalculator(minTempSequence, minTemp, maxTemp, maxLen, mode: (min or max)) &rarr; overlap_length_array
5. Go through the DNA sequence, using the cutting information to actually perform the cut.

Tuesday, 1/16/2024

Today, my primary task is to document my code better. Once that is done, I may start working on implementing the protein to
DNA conversion function

Update: I have finished refactoring data_parser. For now, I am happy with the refactoring and level of documentation.

These are some of the tasks ahead of me:
Converting protein to oligonucleotide.
Temperature function


Monday , 1/15/2024

Daily Debrief:
Yesterday, I completed the round trip of amino set to codon set. I finished the refactoring, and 
both convert_from_amino_set_to_codon_set and convert_from_codon_set_to_amino_set have been completed. A roundtrip 
test starting from an amino set, and going to a codon test has passed successfully. Perhaps I would want a round trip
test that starts from a codon set and goes to an amino set, however I will leave this for now.

Update 1: I have made some progress on documentation, I am approximately 1/3 the way through. At that point, I got sidetracked
by wanting to set up an API using Sphinx. This is not a simple matter, and I will leave this task for a future time.

Sunday, 1/14/2024

Daily Debrief:
Over the weekend, I have been working on completing the amino_acid nucleotide conversion functions. Hopefully, I will be
completed with them today. This will require a bit of refactoring, to decide exactly which funcitons should be doing
which jobs. Then, I suppose, I shall move on to the problem of converting an entire amino-acid chain into
DNA, taking into account codon frequencies.

Brainstorming Session:
So, here is the current flow of data regarding amino_acids and codons.
1. We start off with an amino acid set: {'N', 'I', 'T', 'S'}. 
2. We then convert it into a binary string using letter_to_binary: '00100000010000011000'
3. We then pass the binary string into the amino_acid to codon dictionary through the convert_from_multi_amino_to_multi_codon
function to return a codon set: {'CTG', 'GAC', 'TGS'}
   1. Right now this function is pointless, as it is literally a single line long. Any time I want to convert from binary to 
codon_set, I could literally just pass the binary number through the dictionary, instead of passing the binary number
and the dictionary through a pointless function that has the same exact functionality, but with extra steps.

To summarize: 
1. letter_to_binary({'N', 'I', 'T', 'S'}) &rarr; '00100000010000011000'
2. convert_from_multi_amino_to_multi_codon('00100000010000011000') &rarr; {'CTG', 'GAC', 'TGS'}

I will combine 

Now let's document the flow of data going the other way.

1. We convert one or more multi_codons into amino_acid sets using convert_from_multi_codon_to_multi_amino
   1. We start off with one or more multi_codons that are converted into a set of corresponding codon sets using 
the helper functionL convert_from_multi_codon_to_single_codon_set
   2. We then convert the combined codon set into an amino_acid set using the helper function: 
convert_from_codon_set_to_amino_acid_set


The only function we have not used is binary_to_letter. Is there any reason why we would need a function that converts
a binary representation into a letter representation? I do not know. For testing purposes, I think that we should have
a test that tests every single combination of amino acids in letter format, going on a round trip. I believe that the binary
representation is only ever an intermediate value, and we never actually want to look at the data in binary format. Perhaps
there may be some situation where having the protein sequence in binary format would result in a more efficient algorithm, 
but if that is the case, then I think that I could just add another function then.

So really, I need only two main functions.
1. convert_from_amino_set_to_codon_set(amino_set, codon_dict) &rarr; codon_set
2. convert_from_codon_set_to_amino_set(codon_set) &rarr; amino_set

Function 1 will be a combination of letter_to_binary() and convert_from_multi_amino_to_multi_codon
Function 2 may have additional helper methods, as it is now, or it may be one long method, depending on what will be more legible




Friday, 1/12/2024

On Sunday, I came to the conclusion that tackling the multi-amino to multi-decodon function is very important.

Update: I have not done any coding, really, and have spent time researching. I have the following problem:
It will be easy to confirm the correctness of the conversion of binary <---> string representation. But how do I confirm
the codons to binary? I think I am confused and don't understand the problem entirely, and am having a hard time putting
it into words.

I need to go from binary to codons, and codons to binary. I need to ensure that the codons are actually coding for the
amino acids I think they are coding for, and that there isn't a faulty conversion function resulting in the test passing
even if the codon does not actually code for said amino acids. I believe that a fault in the conversion from binary to 
string function will be caught, as converting from codons -> amino acid is a function that I will implement, but whose
mechanics exist beyond the scope of the program, Even if I implement the conversion from codon -> amino acid correctly,
the round trip test will not pass if the binary -> amino acid conversion function is a faulty conversion

binary -> string -> codons
binary: 01100110000110110000	 
acid: RNQEKMPS
codons: CMK KGB TDT


Sunday, 1/7/2024

Daily Debrief: I think I can finally make a project outline, for what tasks still need to be completed. 

Update 1: Brainstorming Session

The following is a list of the tasks that still need to be completed.
1. Amino Data function: perhaps have a file and data-structure that contains important information for every amino acid.
Information to be stored could potentially be full name, 3 letter abbreviation, 1 letter abbreviation, codons, as well 
as minimum and maximum number of C/Gs used. The file could be as simple as a csv file, and the data structure would be 
a dictionary
2. Conversion functions from going from DNA to Protein and vis-versa. Perhaps for DNA to Protein, I could use the inbuilt
bio-python library.
3. DNA melting temperature calculator
4. Uncertain number of further tasks associated with actually performing the calculation.

In regard to the Amino_acid data, I am unsure if this data is already stored in biopython.
I would need 
1. the 1 letter abbreviation to convert from DNA to protein, but I believe that biopython already has such a function
2. THe min and max CG content is currently needed

There are 2 use cases for amino acid data. Use case 1 is going from protein to DNA. Use case 2 is going from DNA to protein.

Use Case 1: Going from Protein to DNA involves converting the entire protein sequence plus multi-amino information into
a DNA sequence with degenerate codons. This really only involves the codon-frequency information as well as the decodon-set
information. It is still up in the air on how to ensure that the resulting DNA contains the right codon frequencies for 
the given organism, but that is not too big of a deal. All I really have to do is go back to version 1 of the project, and
check to see if the solution I used then is satisfactory. The real aspect that I need to tackle is getting from a multi-amino
set to a multi-decodon set. In fact, this is what I will be working on right now.


1/6/2024

Update 1: I completed all 4 data extraction and storage functions. All 4 unit tests have passed successfully. Time to
move on to the next stage of the project. : )

Daily Debrief: Yesterday, I ran into a bit of a roadblock in regard to finding good quality data on organism codon frequencies. 
I finally found a source that I would feel comfortable citing, as it gets its data from RefSeq. I think I will be able 
to finish file input tonight.

1/5/2024

Note: I got the human codon frequencies data from the following link https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs
which contains 123938 CDSs (80831647 codons), taken from RefSeq.

Daily Debrief: Yesterday, I was able to implement the extraction and storage of the multi-amino acid data. I then moved 
on to dealing with the multi-amino data files. I completed that as well, however, all the storage functions return a
huge string with no tabs, or new line characters. Therefore, the output is not very human-readable. THis is because it makes
string comparison easier. If in the future,  I want to be able to make sense of the output, the function will have to be 
rewritten, or a new function will have to be made. Another note is that the two functions handling the multi-amino files
take like 6 seconds to run consecutively. That is quite long, but it is a one time cost that can be run in the background,
so it is not too big of a deal. For now, all that I have left for file i/o is the organism codon-frequency data. Hopefully
that will be completed today.

1/4/2024

Daily Debrief: Yesterday, I completed the implementation of the store_protein_sequence and extract_protein_sequence. I 
realized that I would have to extract and store the multi-amino data in a separate function. I realized that these are 2
separate tasks and so require 2 separate functions. Implementing the extraction and storage of the multi-amino data is 
today's goal. We will see how that goes.

1/3/2024

Good morning. Today, I will continue the progress I have made from yesterday. I would like to complete step 4 today, and get working
on step 5.


Update 2: 

I have now finished writing the unit test for all the sequence information given by the user. The Unit Test successfully
fails. It is now time to implement all the functions needed to get the Unit Test working. That means I finally need to 
decde on the data structure of the data once it has been created. I know that I had indicated that I do not want to use
OOP in this project. I believe that this was a wrong decision. The classes can be relatively basic, but I think that they
are necessary.




Update1: 

Right now, I am working on setting up unit testing for the function(s) dealing with parsing protein sequence
information. As noted on 12/27/2023, I have decided to take a round trip approach. I am re-encoding the data

back into a file. Then, I will compare the resulting file with the original to ensure no loss of data. The trouble I am 
currently running into is decinding where to put the convert_back_to_file function, that will return the data to
its original file format. Here is what I know so far. The file needs to be a temporary file. 

Option 1:
convert_back_to_file is a function in the data_parser class.

Option 2:
convert_back_to_file is a helper function in the test_data_parser file.

I think that option 2 is the better option, as there is no reason to call the function except for testing purposes. 
Perhaps I am mistaken, however, and there will be a case in the future. I am perfectly fine with the function being
housed in the data_parser class. The real problem, is that I need to make sure that the temporary file is deleted
after it has been used. But it needs to be analysed by the unit test. So The file really needs to be created by
the unit test and then passed to the convert_back_to_file function. But I don't like the idea of the return object being
passed to the function as a parameter. Is this a bad idea? More research is needed.


1/2/2023

Happy New Year! I just finished my Application to RICE university, and now have time to dedicate fully towards
this project. Let us see how much progress we can make in the next three days!

Things that are on the to-do list
1. Make a list of all the files that need to be read to
2. decide naming conventions for the data that comes from each of these files
3. set up the function shells for:
   1. parsing all the data from the files and turning them into data objects
   2. transferring the data back from data objects into data stored in files
   3. unit tests to confirm that round trip results in no lost information
4. Code up the unit tests ensuring that they all fail
5. code all the individual functions until the unit tests pass

Note: Step 5 will be problematic, as it will require me to make difficult decisions on how to store the data
Look at notes from 12/26/2023 to see some of my initial brainstorming on the matter.

Step 1: a list of all data sources:
1. protein sequence information (provided by the user of the program)
2. codon frequency information: the frequency of each codon in the DNA of a given organism
3. degenerate codon information: the file containing the solutions to optimally coding for any given codon-set
codon set.
4. Technically not data from a file, but all the various parameters set by the user of the program
such as temperature and length ranges.

I will only worry about 1, 2 and 3 for now.

Step 2: Naming conventions
How about:
1. **protein_Sequence**
2. this is harder. Do we store every single organism's codon frequency as a different object? I guess that for simplicity,
we will keep all the organismal codon frequency information in their file form, until the user submits their data
which will include what organism they want to optimize for. Then, the program will parse said's organisms data, and put it
into a data object called: **organism_Codon_Frequencies**
3. For the degenerate codon information, does this information change based on what organism is being selected? I do not
know. I do not believe so, but I will have to ask Dr. Papamichail when I talk to him next. I will continue as if the information
is universal. The information can be decoded at the start of the program, even before the user submits their request. Name:
**degenerate_Codons_Key**

Step 3:
12/27/2023

Today, I want to set up unit testing for this project. It will make coding much easier long term. However,
this leads me down the path of how I want to organize my code.  

Idea: I could do round-trip verification as a unit test for the various file readers


12/26/2023

I have been away from this project for over a month, and I am having a hard time getting back into the swing
of things. The following is a short list of some tasks that I have in front of me:

Right now, I need to decide the format of the fasta files. I suppose I would like the ability to change the input
mode in the future without having to rewrite a bunch of other extraneous code, as I had to do last time. In that case,
the input method doesn't even really matter. Its just making smart decisions when it comes to data storage that is
the real issue.

For simplicity, there will be only one fasta file as input.

Update: I have spent about a half hour thinking about the problem. Hopefully, I will be using BioPython
to streamline some of the data storage. However, BioPython does not have a way to store sequences with degenerate
amino acids at certain positions. Perhaps, I can modify some of the classes to shoehorn in this functionality. As a 
side note, using biopython means that to some extent, my code will be opbect oriented, despite what I outlined on 11/13/2023.
I don't see this as particularly problematic.

11/13/2023

I have made a few decisions today in regard to the overall design that I would like to document.
1. No Database: There is not much writing of data, and the bulk data i.e. the organism data only really needs to be read
2. Not OOP: The format of this application is a glorified script, and as such, OOP is not the correct paredigm to use.

