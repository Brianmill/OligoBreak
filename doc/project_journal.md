
12/27/2023

Today, I want to set up unit testing for this project. It will make coding much easier long term. However,
this leads me down the path of how I want to organize my code.  

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

