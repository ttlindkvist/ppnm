CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONY: all clean

all: output.out
	cat $<
		
output.out: main
	./$< > $@
	
main: main.o matrix.o jacobi_diag.o

clean:
	$(RM) main *.o *.out 