CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONY: all clean

all: out.txt test.png
	cat out.txt
	
out.txt: main
	./$< > $@

test.png: data.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "";\
	set key top left;\
	set tics out;\
	set xlabel "";\
	set ylabel "";\
	set xrange [0:20];\
	set yrange [-1.8:1.8];\
	plot \
	 "$<" using 1:2 with line color red title "test"\
	' | pyxplot

main: main.o ode.o ode.h matrix.c matrix.h

clean:
	$(RM) main *.o *.out