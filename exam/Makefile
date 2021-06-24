CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONY: all clean

all: output.out timing.png
	cat $<

timing.png: timing.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "SVD Timing (on NxN matrix)";\
	set key top left;\
	set tics out;\
	set xlabel "N";\
	set ylabel "Time [ms]";\
	set xrange [10:500];\
	set yrange [0:1000];\
	a=1;\
	b=2;\
	f(x)=a*x**b;\
	fit f(x) "$<" using 1:2 via a,b;\
	plot \
	 "$<" using 1:2 with line color red title "My two sided jacobi SVD",\
	 f(x) with line lt 2 color black title "$$\\alpha N^{%.2f}$$"%(b)\
	' | pyxplot
	
output.out timing.out: main
	./$< > output.out
	
main: main.o matrix.o SVD_two_jaco.o

clean:
	$(RM) main *.o *.out *.png