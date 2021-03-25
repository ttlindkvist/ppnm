CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONY: all clean

all: output.out eigenfunc.png
	cat $<

eigenfunc.png: eigenfunc.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Eigenfunctions to particle in a box";\
	set key above;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [0:1];\
	set yrange [-1.5:1.5];\
	plot \
	 "$<" index 0 using 1:2 with line color black title "n=0",\
	 "$<" index 1 using 1:2 with line color red title "n=1",\
	 "$<" index 2 using 1:2 with line color blue title "n=2",\
	 "$<" index 3 using 1:2 with line lt 2 lw 2 color black title "n=0 - analytic",\
	 "$<" index 4 using 1:2 with line lt 2 lw 2 color red title "n=1 - analytic",\
	 "$<" index 5 using 1:2 with line lt 2 lw 2 color blue title "n=2 - analytic"\
	' | pyxplot

output.out eigenfunc.out: main
	./$< > output.out
	
main: main.o matrix.o jacobi_diag.o

clean:
	$(RM) main *.o *.out *.png