CXX = g++
CXXFLAGS = -Wall -std=c++2a
CXXFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lstdc++
LDLIBS += $(shell gsl-config --libs)

all: out.txt interp.png
	cat out.txt

out.txt: main
	./$< > out.txt

main: main.cpp ann.cpp ann.hpp minimization.cpp minimization.hpp

interp.png: interp.out Makefile
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

clean:
	find . -type f -executable -delete
	$(RM) *.o out.* *.out *.png