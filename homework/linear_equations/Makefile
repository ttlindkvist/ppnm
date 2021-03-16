CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all clean

all: result.out
	cat $<

result.out: main
	./$< > result.out

main: main.o

clean:
	$(RM) main *.o *.out