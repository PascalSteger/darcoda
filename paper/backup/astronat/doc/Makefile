
# $Id: Makefile,v 1.4 1999/09/01 03:33:59 jbaker Exp $

MAKEINFO = makeinfo
TEXI2DVI = texi2dvi
TEXI2HTML = texi2html
DVIPS = dvips

TEXI2HTML_FLAGS =  -verbose -split_chapter

HTMLDIR = ./html
TEXDIR = ./tex

SRC = astronat
TEXI = $(SRC).texi


all : info text html # ps

info : $(TEXI)
	$(MAKEINFO) --verbose $(TEXI)

text : $(TEXI)
	$(MAKEINFO) --no-headers $(TEXI) > $(SRC).txt

html : $(TEXI)
	if [ ! -d $(HTMLDIR) ]; then mkdir $(HTMLDIR); fi
	cd $(HTMLDIR); $(TEXI2HTML) $(TEXI2HTML_FLAGS) ../$(TEXI); \
	if [ -f index.html ]; then $(RM) index.html; fi; \
	ln -s astronat_toc.html index.html

ps : $(TEXI)
	if [ ! -d $(TEXDIR) ]; then mkdir $(TEXDIR); fi
	cd $(TEXDIR); $(TEXI2DVI) ../$(TEXI); $(DVIPS) -o $(SRC).ps $(SRC); \
	mv $(SRC).ps $(SRC).dvi ..

clean :
	-$(RM) $(SRC).info $(SRC).txt $(SRC).dvi $(SRC).ps
	-$(RM) -r $(HTMLDIR) $(TEXDIR)
