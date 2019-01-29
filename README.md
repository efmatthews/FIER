# README #

### The FIER Repository ###

* This repository is used for the storage of and collaboration on the FIER (Fission Induced Electromagnetic Response) code.

* For more details on FIER, please see our website: https://bang.berkeley.edu/fier

* Use C++ 2011 or higher. 

* The FIER code is stored in the FIER directory.

* For a detailed walkthrough of building and running FIER, compile (with LaTeX) and read FIER\_Manual\main.pdf. A makefile is included to build the LaTeX manual. 

### Set Up ###
* All code is in the FIER directory.

* Because the program was written in C++, make sure you have access to a C++ compilier, g++ on unix systems is most commonly used.

* Download FIER. To run the code open the terminal and enter the bash command: `make`

* FIER takes an input deck argument. Running `make` will call the default (but editable) input deck: FIER/deck.txt

* No other configuration is necessary.


### Input Deck Structure  ###
The input deck has a specific structure. An example can be found in input\_deck.txt  
Each column in the table represents a line. Include all lines unless marked optional.

---
SYNTAX   |   Definition
-----------|:-------
MODE:X | X = MONTECARLO N run Monte Carlo uncertainty analysis with N trials, X = SINGLE no Monte Carlo uncertainty analysis.  
X DECAY PREDICTION      |X = ON, decay prediction on X = OFF, decay prediction off  
"/path/to/isotopes.csv" |half lives file input with .csv format (see isotopes2.csv for example)  
"/path/to/decays.csv"   |Decay mode file input with .csv format   (see decays2.csv for example)  
"/path/to/gammas.csv"   |Gamma intensity file input with .csv format   (see gammas_5per.csv for example)  
"/path/to/yields.csv"   |Fission Yields file with .csv format        (see U235_500.csv for example)  
"/output/chains.csv"    |Decay chains output file. Note: for any output file replace with NONE to suppress creation  
"/output/stems.csv"     |Decay stems output file 
"/output/pops.csv"      |Decay populations output file 
"/output/Ys.csv"        |Gamma energies output file
"/output/errlog.txt"    |Error log output file
SELECT:X                |Allows to specify the species one is looking for. X = ON to specify, X = OFF to not.  
Z1,A1,I1                |Z = atomic number, A = atmoic mass, I = isomeric number. Use a new line for each species sought  
...  |
Zn,An,In  |
INITIALIZE              |Populations of initial species (required key word)
Z1,A1,I1,POP1           |Optional initial populations
...  |
Zn,An,In,POPn |  
IRRADIATION             |Signifies start of irradiation profile (required key word)
ATI1,FR1                |ATI = absolute time interval (seconds), FR = fissions/second
...                     |Add one per irradiation step desired
ATIn,FRn |
POPULATIONS             |Signifies start of population calculations (required key word) 
TIME                    |Add absolute time to record populations at desired times. TIME = time (seconds w.r.t. t = 0)  |
COUNTS                  |Signifies start of counting scheme (required key word)
METHOD:ANALYTICAL       | T1,T2 start and end times of counting bin
T1,T2                   |
.... |
Tn1,Tn2  |

### Input Files needed by the Input Deck ###

* The input libraries isotopes.csv, decays.csv, gammas.csv can be found pre-parsed in the FIER\input\_data subdirectory.  

* The yields files must be supplied by the user but are taken from England and Rider (see manual).  

These can be found for individual species by going to https://www.nndc.bnl.gov/sigma/index.jsp > selecting neutron-induced-fission-yields from  
the sublibrary dropdown > clicking on the element, then isotope desired (isotope on the right bar) > then clicking on the interpreted link next to fission yields.
This produces a library of yields in the format of Z, A, FPS, Yield, Uncertainty, which can then be parsed into a .csv to be in the correct yields format.

### Cite This Work ###

* Please cite the FIER manuscript published in NIM-A if you use FIER in your published work! A BibTeX file for this manuscript can be found in the file FIER_NIMA.bib. 


### Contribution Guidelines ###

* Experimental work is done on the FIER branch. Production code is hosted in the Master branch. Once FIER-branch work has been tested to satisfaction, merge FIER into Master.



### Contacts ###

* Project Lead and Author:
Eric Matthews - efmatthews@berkeley.edu

* Principle Investigator:
Dr. Bethany Goldblum, bethany@nuc.berkeley.edu

* Other Contributors:
Matthew Shinner - matthew.shinner@berkeley.edu