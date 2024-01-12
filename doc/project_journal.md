

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
   2. transfering the data back from data objects into data stored in files
   3. unit tests to confirm that round trip results in no lossed information
4. Code up the unit tests ensuring that they all fail
5. code all the individual functions until the unit tests pass

Note: Step 5 will be problematic, as it will require me to make difficult descisions on how to how to store the data
Look at notes from 12/26/2023 to see some of my initial brainstorming on the matter.

Step 1: a list of all data sources:
1. protein sequence information (provided by the user of the program)
2. codon frequency information: the frequency of each codon in the DNA of a given organism
3. degenerate codon information: the file containing the solutions to optimally coding for any given codon set
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

