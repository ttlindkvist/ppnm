CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

all: out.txt
	cat out.txt

out.txt: main
	./$< > out.txt

main: main.c integration.c integration.h

clean:
	$(RM) main *.o out.txt *.out