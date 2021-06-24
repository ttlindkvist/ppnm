CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONY: all clean

all: output.out
	cat $<

output.out: main
	./$< > output.out
	
main: main.o matrix.o SVD_two_jaco.o

clean:
	$(RM) main *.o *.out *.png