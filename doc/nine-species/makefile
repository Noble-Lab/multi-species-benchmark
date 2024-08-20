all: nine-species.pdf \
	cover-letter.pdf \
	response.pdf

%.pdf: %.tex
	xelatex $*.tex
	if ( grep -q bibliographystyle $*.tex ) ; then \
		bibtex $* ; \
		pdflatex $*.tex ; \
	fi
	xelatex $*.tex

%.pdf: %.fig
	fig2dev -L pdftex -p dummyarg $*.fig $*.pdf

clean:
	rm -f *.dvi *.bbl *.aux *.blg *.log *.ps *.out

