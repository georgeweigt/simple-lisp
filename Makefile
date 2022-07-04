.PHONY: default clean

default:
	make lisp
	make man.pdf

lisp: lisp.c
	gcc -O0 -o lisp lisp.c

man.pdf: man.tex
	pdflatex man.tex
	rm -f man.aux man.log

clean:
	rm -f lisp man.pdf
