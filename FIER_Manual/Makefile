#Variables
#----------------------------------------------------------------------------------------
manuscript = main
latexopt = -file-line-error -halt-on-error -shell-escape
#----------------------------------------------------------------------------------------



#Main - LaTeX build the manuscript
#----------------------------------------------------------------------------------------
$(manuscript).pdf: clean
	pdflatex $(latexopt) $(manuscript).tex
	pdflatex $(latexopt) $(manuscript).tex
	pdflatex $(latexopt) $(manuscript).tex
#----------------------------------------------------------------------------------------

#Clean the repo
#----------------------------------------------------------------------------------------
clean :
	rm -f *.aux *.log *.bbl *.lof *.lot *.blg *.out *.toc *.run.xml *.bcf *.spl
	rm -f $(manuscript).pdf
#----------------------------------------------------------------------------------------