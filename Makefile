.PHONY: clean

lisp: lisp.c
	gcc -O0 -o lisp lisp.c

clean:
	rm -f lisp
