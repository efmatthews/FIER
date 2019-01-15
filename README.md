# README #

### The FIER Repository ###

* This repository is used for the storage of and collaboration on the FIER (Fission Induced Electromagnetic Response) code.

* Use C++ 2011 or higher. 

* Use the C++ FIER code, stored in the FIER directory.

* For a detailed walkthrough of building and running FIER, compile (with LaTeX) and read FIER\_Manual\main.pdf.

### Set Up ###
* All code is in the FIER\_C++ directory.

* Because the program was written in C++, make sure you have access to a C++ compilier, g++ on unix systems is the most commonly used one.

* Download FIER (FIER\_v1.cpp). To run the code open the terminal and enter the bash command: $make

* FIER takes an input deck argument. Running make will call the default (but editable) FIER/deck.txt

* No other configuration is necessary at this time.


### Input Deck Structure  ###
The input deck has a specific structure. An example can be found in input\_deck.txt  
Each column in the table represents a line. Include all lines unless marked optional.

---
SYNTAX   |   Definition
-----------|:-------
MODE:X | X = MONTECARLO N run montecarlo analysis with N trials, X = SINGLE no montecarlo analysis.  
X DECAY PREDICTION      |X = ON, decay prediction on X = OFF, decay prediction off  
"/path/to/isotopes.csv" |half lives file input with .csv format (see isotopes2.csv for example)  
"/path/to/decays.csv"   |Decay mode file input with .csv format   (see decays2.csv for example)  
"/path/to/gammas.csv"   |Gamma intensity file input with .csv format   (see gammas_5per.csv for example)  
"/path/to/yields.csv"   |Fission Yields file with .csv format        (see U235_500.csv for example)  
"/output/chains.csv"    |Decay chains //For any output replace with NONE to supress creation  
"/output/stems.csv"     |Decay stems  
"/output/pops.csv"      |Decay populations  
"/output/Ys.csv"        |Gamma energies  
"/output/errlog.txt"    |Error log  
SELECT:X                |Allows to specify the species one is looking for. X = ON to specify, X = OFF to not.  
Z1,A1,I1                |Z = atomic number, A = atmoic mass. Use a new line for each species sought  
...  |
Zn,An,In  |
INITIALIZE              |Populations of initial species
Z1,A1,I1,POP1           |Optional
...  |
Zn,An,In,POPn |  
IRRADIATION             |signifies start of radiation profile
ATI1,FR1                |ATI = absolute time interval (seconds), FR = fissions/second (optional)
...                     |Add one per irradiation step desired
ATIn,FRn |
POPULATIONS             |Non optional line  
TIME                    |Add absolute time to record populations at desired times. TIME = time (seconds w.r.t. t = 0)  |
COUNTS                  |
METHOD:ANALYTICAL       | T1,T2 start and end times to integrate over
T1,T2                   |Optional
.... |
Tn1,Tn2  |

### Input Files needed by the Input Deck ###

* The input libraries isotopes2.csv, decays2.csv, gammas\_endf.csv can be found pre-parsed in the FIER\input\_data subdirectory.  

* The yields files must be supplied by the user but are taken from England and Rider (see manual).  

These can be found for individual species by going to https://www.nndc.bnl.gov/sigma/index.jsp > selecting neutron-induced-fission-yields from  
the sublibrary dropdown > clickign on the element, then isotope desired (isotope on the right bar) > then clicking on the interpreted link next to fission yields.
This produces a library of yields in the format of Z, A, FPS, Yield, Uncertainty, which can then be parsed into a .csv to be in the correct yields format.


### Contribution Guidelines ###

* Experimental work is done on the FIER branch. Production code is hosted in the Master branch. Once FIER-branch work has been tested to satisfaction, merge FIER into Master.



### Contacts ###

* Project Lead:
Eric Matthews - efmatthews@berkeley.edu

* Admins:
Dr. Bethany Goldblum, bethany@nuc.berkeley.edu

* Other contributors:
Matthew Shinner - matthew.shinner@berkeley.edu