CXX = g++
CXXFLAGS = -Wall -std=c++17
CXXFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lstdc++
LDLIBS += $(shell gsl-config --libs)

all: out.txt higgs_fit.png
	cat out.txt

out.txt higgs_fit.out: main
	./$< > out.txt

higgs_fit.png: higgs_fit.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Higgs fit";\
	set key top right;\
	set tics out;\
	set xlabel "E (GeV)";\
	set ylabel "signal";\
	set xrange [100:160];\
	set yrange [-5:10];\
	plot \
	 "$<" using 1:2 with line color black title "Fit",\
	 "higgs.in" using 1:2:3 with errorbars color red title "Data"\
	' | pyxplot

main: main.cpp minimization.cpp minimization.hpp

clean:
	find . -type f -executable -delete
	$(RM) *.o out.* *.out