
# ----------------------------------------------------------------------
# Makefile for processing latex/bibtex files, with commands for
# generating files which can be submitted electronically (uses the
# nat2jour.pl script).  Edit the first section below as necessary.
# The following commands operate on a file FOO.tex which uses the
# BIB file for its bibliography:
# 
#   make FOO.bbl        (generate bibliography file)
#   make FOO.dvi        (generate dvi file)
#   make FOO.ps         (generate postscript file)
#   make FOO.pr         (print postscript)
#   make clean          (get rid of some tex-generated files)
#   make dviclean       (get rid of dvi, useful for forcing re-make)
#   make realclean      (get rid of all generated files)
#
# FOO can be replaced by FOO-JSUF (see below) to make submitted 
# versions.  Make sure TEXINPUTS, BSTINPUTS, and BIBINPUTS 
# environment variables are set correctly if the necessary .sty,
# .bst, and .bib files are not in the current directory.
#
# 1998/08/11  Jonathan Baker <jbaker@astro.berkeley.edu>
#
# $Log: Makefile,v $
# Revision 1.2  1998/08/30 21:06:13  jbaker
# Added RCS keyword.
#
# ----------------------------------------------------------------------

# Commands
LATEX = latex
BIBTEX = bibtex -verbose
DVIPS = dvips
LPR = lp
NAT2JOUR = nat2jour.pl
RM = /bin/rm -f

# Flags
DVIPS_FLAGS =
LPR_FLAGS = -c -d p544 -o duplex

# Names of the default files to process, WITHOUT .tex extensions
TARGS = example1

# Names of bibliography database files WITH full path names (if not
# in .) and WITHOUT .bib extensions.  Note if these are not in the 
# current directory, you must set the BIBINPUTS environment variable 
# in order for BIBTEX to find them.
BIB = example1

# Suffix added to files generated for journal submission
JSUF = -apj
#JSUF = -mn

# ----------------------------------------------------------------------
# Special function targets

# Clear the default implicit rules
.SUFFIXES:

# Don't let GNU make delete these if they are intermediates
.SECONDARY : $(TARGS:=.bbl) $(TARGS:=$(JSUF).bbl) $(TARGS:=$(JSUF).tex)

# ----------------------------------------------------------------------
# Default rule: make all the dvi's
# This is invoked when `make' is called with no arguments

default : $(TARGS:=.dvi) $(TARGS:=$(JSUF).dvi)

# ----------------------------------------------------------------------
# Implicit rules
# These are invoked by `make FOO.X' if `X' is one of the suffixes and
# there is a file `FOO.tex'

# Build files for journal submission, must come before the rule for
# %.bbl to override it (GNU make)
%$(JSUF) %$(JSUF).tex %$(JSUF).bbl : %.bbl
	$(NAT2JOUR) $* $*$(JSUF)

%.bbl : %.tex $(BIB:=.bib)             # bibliography
	$(LATEX) $*
	$(BIBTEX) $*

%.dvi : %.bbl                          # dvi files
	$(LATEX) $*
	$(LATEX) $*

%.ps : %.dvi                           # postscript files
	$(DVIPS) $(DVIPS_FLAGS) -o $@ $<

%.pr : %.dvi                           # print
	$(DVIPS) $(DVIPS_FLAGS) $< | $(LPR) $(LPR_FLAGS)

# ----------------------------------------------------------------------
# Clean up

clean :                                # log and auxiliary files
	-$(RM) *.aux *.log *.blg

dviclean :
	-$(RM) *.dvi

bibclean : clean                       # above + bibliography files
	-$(RM) *.bbl

realclean : clean dviclean bibclean    # all of the above + ps and journal
	-$(RM) *.ps *$(JSUF).*
