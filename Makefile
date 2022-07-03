.PHONY: clean

simple-lisp: simple-lisp.c
	gcc -O0 -o simple-lisp simple-lisp.c -lm

clean:
	rm -f simple-lisp
