CFLAGS = -O -Wall $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: all clean

all: out.txt erf.pyxplot.png bessel.pyxplot.png
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
bessel.pyxplot.png: bessel.dat Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set yrange[-0.5:1.1];\
	set title "Bessel functions of the first kind";\
	plot \
	 "$<" using 1:2 with line title "$$ J_0(x) $$"\
	,"$<" using 1:3 with line title "$$ J_1(x) $$"\
	,"$<" using 1:4 with line title "$$ J_2(x) $$"\
	' | pyxplot

main: main.o
main.o: main.c

clean:
	$(RM) main *.o *.dat *.png