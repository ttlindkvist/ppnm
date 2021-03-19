CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS =  -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all clean

all: output.out
	cat $<

output.out: main
	./main > output.out

main: main.o matrix.o matrix.h

clean:
	$(RM) main *.o *.out