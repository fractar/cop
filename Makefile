# LaTeX compiler
TEX_TO_PDF=pdflatex

# Index builder
MI=makeindex

# Index style option
MIFLAGS=-s

# Tex, style, db
SRC=compterendu.tex

# Directories
DIR=`$(PWD)`

# PDF name
PDF_NAME="compterendu"

PDF_READER=evince

.PHONY:

all: index view

view: $(PDF_NAME)
	$(PDF_READER) "$<.pdf" &

$(PDF_NAME):
	$(TEX_TO_PDF) --jobname=$(PDF_NAME) $(SRC)	

index: $(PDF_NAME)
	$(TEX_TO_PDF) --jobname=$(PDF_NAME) $(SRC)

clean:
	rm -f *.acn
	rm -f *.acr
	rm -f *.alg
	rm -f *.aux
	rm -f *.glg
	rm -f *.glo
	rm -f *.gls
	rm -f *.idx
	rm -f *.ilg
	rm -f *.ild
	rm -f *.ind
	rm -f *.ist
	rm -f *.log
	rm -f *.lof
	rm -f *.lol
	rm -f *.lot
	rm -f *.out
	rm -f *.toc
	rm -f latex_files/*.aux
	rm -f latex_files/others/*.aux

mrproper: clean
	rm -f $(PDF_NAME).pdf

distclean: mrproper
	rm -f *~
	rm -f .log
	rm -f *.bak
	rm -f *.bck
	rm -f .depend 
	rm -f Makefile.rules
	rm -f $(PDF_NAME).tar

dist: distclean
	tar -cvf /tmp/$(PDF_NAME).tar ../$(DIR)
	mv /tmp/$(PDF_NAME).tar .
