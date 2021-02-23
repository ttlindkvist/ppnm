CFLAGS = -O -Wall $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: all clean

all: out.txt erf.pyxplot.png gamma.pyxplot.png
	cat $<

out.txt erf.dat gamma.dat: main
	./$< > out.txt

erf.pyxplot.png: erf.dat Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf by integration"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	' | pyxplot
gamma.pyxplot.png: gamma.dat Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set yrange[-5:7];\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "gamma by integration"\
	,"$<" using 1:3 with line title "gsl-sf-gamma"\
	' | pyxplot

main: main.o
main.o: main.c

clean:
	$(RM) main *.o *.dat