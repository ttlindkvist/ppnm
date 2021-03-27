CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONY: all clean

all: output.out eigenfunc.png jacobi.timing.png gsl.timing.png
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
jacobi.timing.png: timing.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Jacobi diagonalization";\
	set key top left;\
	set tics out;\
	set xlabel "Matrix size N";\
	set ylabel "Time [ms]";\
	set xrange [10:200];\
	set yrange [0:800];\
	a=1;\
	f(x)=a*x**3;\
	fit f(x) "$<" using 1:2 via a;\
	b=1;\
	g(x)=b*x**3;\
	fit g(x) "$<" using 1:4 via b;\
	c=1;\
	h(x)=c*x**3;\
	fit h(x) "$<" using 1:5 via c;\
	plot \
	 "$<" using 1:2 with points color red title "Unchanging diag as convergence criteria",\
	 "$<" using 1:5 with points color black title "Only updating upper triangle on sweep",\
	 "$<" using 1:4 with points color black title "Off-diag sum as convergence criteria",\
	 f(x) with line lt 2 color black title "c1*N**3",\
	 g(x) with line lt 2 color black title "c2*N**3",\
	 h(x) with line lt 2 color black title "c3*N**3"\
	' | pyxplot
gsl.timing.png: timing.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Timing";\
	set key top left;\
	set tics out;\
	set xlabel "Matrix size N";\
	set ylabel "Time [ms]";\
	set xrange [10:200];\
	set yrange [0:20];\
	b=1;\
	g(x)=b*x**3;\
	fit g(x) "$<" using 1:3 via b;\
	plot \
	 "$<" using 1:3 with points color black title "gsl implementation",\
	 g(x) with line lt 2 color black title "c2*N**3"\
	' | pyxplot

output.out eigenfunc.out timing.out: main
	./$< > output.out
	
main: main.o matrix.o jacobi_diag.o

clean:
	$(RM) main *.o *.out *.png