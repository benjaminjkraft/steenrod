SOURCES=$(wildcard *.tex)
TARGETS=$(SOURCES:.tex=.pdf)

all: $(TARGETS)

%.pdf: %.tex %.bib
	pdflatex $<
	biber $(<:.tex=)
	pdflatex $<

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

clean:
	rm -f *.aux *.bbl *.bcf *.blg *.log *.out *.run.xml *._aux *._log *.nav *.toc *.synctex.gz *.tex.project.vim
