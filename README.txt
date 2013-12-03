################################################
The included code is liscenced under GPL, see attached 
Copyright 2013 Robert Cope cope.robert.c@gmail.com
################################################

This program creates simulated populations, designed to mimic the characteristics of wild mammals with long lifespans, displaying random and inter-generational mating.
Population size and sampling proportion can be varied, along with intended genetic variability. Output includes inferred maturity/size class of individuals.
This similation was designed for the testing of 'PR-genie' software developed for:

Cope RC, Lanyon JM, Seddon JM, Pollett PK (In Press) Development and testing of a genetic marker based pedigree reconstruction system ''PR-genie'' incorporating size-class data. Molecular Ecology Resources.


################################################
Instructions:

1) Compile the code via command line, e.g.
gcc genData.h genData.c generate.c -o generatePopulation

2) Run the program. It takes 4 arguments:
- A (2-digit integer) tuning parmeter that determines the genetic variability in the population. Smaller numbers will produce higher heterozygosity, e.g. 10 gives an average heterozygosity >0.80, 50 gives ~0.4, 80 gives ~0.35.
- The proportion (i.e., percentage) of the population to be returned in the sample (integer 1-100). 
- The size of the population. Integer.
- The file to which this should all be output, e.g. text.txt
e.g.:
./generatePopulation 50 80 500 text.txt

Lines of the output file look like this:
101585 (101075 101333) 1 165 2 1 163 1 4.4 0.1 5.10 11.12 4.3 5.12 10.10 9.8 6.3 4.6 3.2 10.10 4.7 1.1 2.4 0.2 2.4 1.3 7.7 2.2 3.2 0.9 1.5 4.2 0.3 
The format is:
a (b c) d e f g h i u.v w.x y.z ....
a is the id of the individual,
b and c are the ids of its parents.
d is a dummy variable not used in this instance.
e is the timestep in which the individual was sampled.
f is its size-class at sampling. 0 for adult, 1 for subadult, 2 for calf.
g is its sex. 1 for male, 0 for female.
h is the year of birth of the individual (for reference; not used).
i is a dummy variable also not used here.
u.v w.x y.z are the alleles at each locus, i.e., a pair is one locus with two alleles.

3) You can use the python script parents-to-xml.py to convert to xml. It takes one parameter - the name of the file to be converted, without the .txt extension.


