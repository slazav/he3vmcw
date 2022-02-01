# grad functions and tests
CC=g++

all: grad_test grad_speed
	./grad_test

grad.c grad_test.c: grad.h check.h
grad_test: grad_test.o grad.o check.o
grad_speed: grad_speed.o grad.o
legg2_test: check.o

convert_prof: convert_prof.o grad.o
convert_prof.o: grad.h

clean:
	rm -f *.o grad_test grad_speed

