.PHONY: clean

lisp: lisp.c
	gcc -O0 -o lisp lisp.c -lm

clean:
	rm -f lisp
