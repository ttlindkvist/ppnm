CXX = g++
CXXFLAGS = -Wall -std=c++17 -fopenmp -lgomp
LDLIBS = -lstdc++ -lm

all: out.txt errorplot.png
	cat out.txt

errorplot.png: error.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Error - area of a circle";\
	set key bottom left;\
	set tics out;\
	set xlabel "N";\
	set ylabel "Error";\
	set xrange [10:2e5];\
	set yrange [1e-7:1e-1];\
	set logscale y;\
	a=1;b=-0.5;\
	f(x)=a*x**b;\
	fit f(x) withouterrors "$<" using 1:2 via a,b;\
	c=1;d=-1;\
	g(x)=c*x**d;\
	fit g(x) withouterrors "$<" using 1:3 via c,d;\
	plot \
	 "$<" using 1:2 with line color red title "pseudo",\
	 "$<" using 1:3 with line color black title "quasi",\
	 f(x) with line lt 2 color black title "$$%.2fn^{%.2f}$$"%(a,b),\
	 g(x) with line lt 2 color black title "$$%.2fn^{%.2f}$$"%(c,d)\
	' | pyxplot

out.txt error.out: main
	./$< > out.txt

main: main.cpp montecarlo.cpp montecarlo.hpp

clean:
	find . -type f -executable -delete
	$(RM) *.o out.* *.out