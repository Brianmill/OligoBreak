# OligoBreak
Program to break a DNA sequence into oligonucleotides that minimizes the manufacturing cost while staying compatable with Gibson Assembly

# Data Sources

primer.dat comes from dnaMATE, a DNA calculator created by the Melo lab in 2005. According to the website the server is hosted on, 
primer.dat "contains a total of three hundred thousands (300.000) oligonucleotide sequences in a size range between 16 
and 30 bases. For each oligonucleotide length, 20.000 sequences are available, which are distributed homogeneously in 10
classes of CG-content percentage. Thus, for each particular length and CG-content class, a total of 2.000 oligonucleotide
sequences are available. For each oligonucleotide sequence, several features and properties were extracted and/or 
calculated, and stored in this database." primer.dat was used to verify the dna melting temperature calculator used in this program

The paper outlining the dnaMATE server:

Citation: Panjkovich A, Norambuena T, Melo F. dnaMATE: a consensus melting temperature prediction server for short DNA 
sequences. Nucleic Acids Res. 2005 Jul 1;33(Web Server issue):W570-2. doi: 10.1093/nar/gki379. PMID: 15980538; PMCID: 
PMC1160140.

The actual server can be found at :http://melolab.org/dnaMATE/ with the data being found at http://melolab.org/dnaMATE/suppl-material-rd.html
These links were working as of January 25th, 2024
