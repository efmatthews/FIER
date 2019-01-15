CXX = g++
CFLAGS = -std=c++11 -g
DECK = deck.txt
TESTDECK = testing/testdeck.txt
OBJS = species_data.o helper_functions.o chains_data.o product_data.o monte_carlo.o

all: fier.exe run clean

%.o: %.cpp makefile
	$(CXX) -c $(CFLAGS) -fPIC $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

fier.exe: $(OBJS)
	$(CXX) $(CFLAGS) -o fier.exe main.cpp $(OBJS)

run:
	./fier.exe $(DECK)

test: fier.exe runtest cleantestoutput cleantest clean
	
runtestdeck:
	rm -f testing/*.csv
	./fier.exe $(TESTDECK)
	
buildtest:
	$(CXX) $(CFLAGS) -o test testing/test.cpp
	
runtest:
	python testing/test.py

cleantestoutput:
	rm -f testing/output/*.csv
	rm -f testing/output/*.txt

cleantest:
	rm -f testing/test.exe
	rm -f *.o
	rm -f test.exe
	
clean:
	rm -f *.o *.d fier.exe
	rm -rf fier.exe.dSYM

docs:
	rm -f index.html
	doxygen Doxyfile
	ln -s ./html/index.html ./index.html

docsw:
	rm -f index.html
	doxygen Doxyfile
	cmd /c "mklink .\index.html .\html\index.html"
	
.PHONY: clean cleantest

montecarlo: format fier.exe shell join clean

format:
	rm -f threads/*.txt
	rm -f threads/*.exe
	rm -f threads/gammas/*.csv
	rm -f threads/populations/*.csv
	rm -f threads/stdevs/*.csv
	$(CXX) $(CFLAGS) -o threads/format.exe threads/format.cpp
	./threads/format.exe
	rm -f threads/format.exe

shell:
	chmod u+x threads/threadscript.sh
	threads/threadscript.sh
	
join:
	$(CXX) $(CFLAGS) -o threads/join threads/join.cpp
	./threads/join
	rm -f threads/join
	
help:
	@echo -e  "\n---------------- FIER MAKE HELP -----------------\nTyping make without a target builds and \
	runs FIER with \"deck.txt\"\
	\n\nAvaliable Targets\n\trun: runs but does NOT compile FIER with deck.txt \
	as an input.\n\tfier.exe: compiles but does NOT run FIER.\n\ttest: runs the unit test \
	that checks for a correct build.\
	\n\tmontecarlo: executes FIER in parallel for easier monte carlo.\
	\n\tdocs: creates Doxygen documentation for unix users.\
	\n\tdocsw: creates Doxygen documentation for windows users.\
	\n\thelp: Prints this message.\
	\n\n NOTE: you need g++ and c++11 to run this software!"
	

-include $(wildcard *.d)

