CXX = g++
CXXFLAGS = -Wall -std=c++17 -Wno-sign-compare
CXXFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lstdc++
LDLIBS += $(shell gsl-config --libs)

all: out.txt wavefunc.png convergence.png
	cat out.txt

wavefunc.png: wavefunc.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Lowest bound state of H";\
	set key top right;\
	set tics out;\
	set xlabel "y";\
	set ylabel "x";\
	set xrange [0:8];\
	set yrange [0:0.4];\
	plot \
	 "$<" using 1:2 with line color black title "Calculated",\
	 "$<" using 1:4 with line color red title "Exact"\
	' | pyxplot
	
convergence.png: convergence.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Convergence";\
	set key top right;\
	set tics out;\
	set xlabel "$$R_{max}$$";\
	set ylabel "Energy";\
	set xrange [0.9:10];\
	set yrange [-0.6:1];\
	plot \
	 "$<" using 1:2 with line color black title "Old boundary condition",\
	 "$<" using 1:3 with line color red title "New boundary condition",\
	 "$<" using 1:4 with line lt 2 color black title "$$\\epsilon$$"\
	' | pyxplot
out.txt error.out wavefunc.out convergence.out: main
	./$< > out.txt

main: main.cpp roots.cpp roots.hpp matrix.c matrix.h function.cpp function.hpp ode.h ode.c dyn_matrix.c dyn_matrix.h

clean:
	find . -type f -executable -delete
	$(RM) *.o out.* *.out