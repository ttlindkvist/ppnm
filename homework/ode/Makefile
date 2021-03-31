CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONY: all clean

all: SHM.png
	
SHM.out: main
	./$<

SHM.png: SHM.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "";\
	set key top left;\
	set tics out;\
	set xlabel "t";\
	set ylabel "x";\
	set xrange [0:20];\
	set yrange [-1.8:1.8];\
	plot \
	 "$<" using 1:2 with line color black title "position",\
	 "$<" using 1:3 with line color red title "velocity"\
	' | pyxplot

main: main.o ode.o ode.h matrix.c matrix.h

clean:
	$(RM) main *.o *.out *.png