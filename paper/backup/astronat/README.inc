%nat----------------------------------------------------------------------
%nat  This file uses the natbib style package along with a modified
%nat  natbib bibliography style (.bst) file.  This allows us to use
%nat  BIBTeX with AAS/MNRAS-like styles.  References should be
%nat  added to the .bib file specified in the \bibliography{} command.
%nat  These need not be limited to the references actually cited; the
%nat  package automatically picks out the right ones.
%nat
%nat  The first time this file is processed, you will need to run
%nat  `latex' followed by `bibtex' followed by two more iterations
%nat  of `latex'.  This must be done whenever bibliographic references 
%nat  need to be updated; otherwise you only need the usual `latex'.
%nat  
%nat  To generate .tex and .bbl files for submission, use `nat2jour.pl'.
%nat  This will produces a .tex file in which all references have been
%nat  inlined and dependencies on natbib have been removed.  The
%nat  appearance of the bibliography may change since natbib redefines
%nat  it.  (This step no longer necessary with AASTeX version 5.)
%nat
%nat  Here are examples of natbib citation commands which may be used:
%nat    \citet{key}              ==>>  Jones et al. (1990)
%nat    \citep{key}              ==>>  (Jones et al. 1990)
%nat    \citep{key1,key2,...}    ==>>  (Jones et al. 1990; Smith 1989; ...)
%nat                                or (Jones et al. 1990, 1991; ...)
%nat                                or (Jones et al. 1990a,b; ...)
%nat    \citep*{key}             ==>>  (Jones, Baker, & Williams 1990)
%nat    \citep[chap. 2]{key}     ==>>  (Jones et al., 1990, chap. 2)
%nat    \citep[e.g.,][]{key}     ==>>  (e.g., Jones et al., 1990)
%nat    \citep[see][p. 34]{key}  ==>>  (see Jones et al., 1990, p. 34)
%nat    \citealt{key}            ==>>  Jones et al., 1990
%nat    \citealt*{key}           ==>>  Jones, Baker, & Williams, 1990
%nat    \citealp{key}            ==>>  Jones et al. 1990
%nat    \citealp*{key}           ==>>  Jones, Baker, & Williams 1990
%nat    \citeauthor{key}         ==>>  Jones et al.
%nat    \citeauthor*{key}        ==>>  Jones, Baker, & Williams
%nat    \citeyear{key}           ==>>  1990
%nat    \citeyearpar{key}        ==>>  (1990)
%nat    \citetext{priv. comm.}   ==>>  (priv. comm.)
%nat
%nat $Name: release-1-4 $
%nat $Id: README.inc,v 1.3 1999/09/01 02:58:23 jbaker Exp $
%nat----------------------------------------------------------------------
