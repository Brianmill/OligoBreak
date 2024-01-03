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

