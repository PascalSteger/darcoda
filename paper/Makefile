# makefile for latex files
# modified from a Makefile by H. Katzgraber 1999

# choose compiler 
CC = pdflatex

# source files ( no .tex extension!!!)
SRC = main
BIB = $(SRC).bib

all: pdf
# plain latex compilation, produces a .dvi
tex:	bib
	$(CC) $(SRC)
	$(CC) $(SRC)

# first generate .aux, then compile library file
bib:	$(BIB)
	$(CC) $(SRC)
	bibtex $(SRC)

# make a ps file
ps:	tex
#	dvips -o $(SRC).ps -t a4 $(SRC).dvi

# make a pdf file
pdf:	tex
	cp main.pdf Steger_2013_MCMC_spherical.pdf
#	dvips -Ppdf -G0  -o $(SRC).ps $(SRC).dvi
#	ps2pdf -dPDFsettings=/prepress $(SRC).ps

# rm the tex crappy files and other twiggle files
clean:
	rm -rf *.aux *.dvi *.log *.toc *.lof *.lot *.blg *.bbl *.end *~

# rm ps and pdf files as well
tidy: clean
	rm -rf $(SRC).ps $(SRC).pdf
