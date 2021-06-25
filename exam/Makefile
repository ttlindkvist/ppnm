CFLAGS = -Ofast -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONY: all clean

all: output.out timing.png
	cat $<

report.pdf: report.tex timing.png
	pdflatex $<
	pdflatex $<

timing.png: timing.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "SVD Timing (on NxN matrix)";\
	set key top left;\
	set tics out;\
	set xlabel "N";\
	set ylabel "Time [ms]";\
	set xrange [0:260];\
	set yrange [0:7000];\
	a=1;\
	b=2;\
	f(x)=a*x**b;\
	fit f(x) "$<" using 1:2 via a,b;\
	c=1;\
	d=2;\
	g(x)=c*x**d;\
	fit g(x) "$<" using 1:3 via c,d;\
	plot \
	 "$<" using 1:2 with points color red title "My two-sided jacobi SVD",\
	 "$<" using 1:3 with points color black title "GSL one-sided jacobi SVD",\
	 f(x) with line lt 2 color black title "$$\\alpha N^{%.2f}$$"%(b),\
	 g(x) with line lt 2 color black title "$$\\beta N^{%.2f}$$"%(d)\
	' | pyxplot
	
output.out: main
	./$< test > output.out

timing.out: main
	./$< time

main: main.o matrix.o SVD_two_jaco.o

clean:
	$(RM) main *.o *.out *.png *.pdf *.log *.aux