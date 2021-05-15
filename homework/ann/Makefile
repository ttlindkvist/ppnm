CXX = g++
CXXFLAGS = -Wall -std=c++17
CXXFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lstdc++
LDLIBS += $(shell gsl-config --libs)

all: out.txt plotA.png plotB.png
	cat out.txt

out.txt plot.out: main
	./$< > out.txt

main: main.cpp ann.cpp ann.hpp minimization.cpp minimization.hpp integration.c integration.h

plotA.png: plot.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "ANN fit";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [-1:1];\
	set yrange [-1:1];\
	plot \
	 "$<" index 0 using 1:2 with points color black title "points",\
	 "$<" index 1 using 1:2 with line color red title "ANN fit"\
	' | pyxplot
plotB.png: plot.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "ANN fit";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [-1:1];\
	set yrange [-4.5:5.5];\
	plot \
	 "$<" index 2 using 1:2 with line lt 2 color black lw 2 title "ANN derivative",\
	 "$<" index 2 using 1:4 with line color red title "Actual derivative",\
	 "$<" index 2 using 1:3 with line lt 2 color black lw 2 title "ANN integral",\
	 "$<" index 2 using 1:5 with line color green title "Numeric integral"\
	' | pyxplot

clean:
	find . -type f -executable -delete
	$(RM) *.o out.* *.out *.png