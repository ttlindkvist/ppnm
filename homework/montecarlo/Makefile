CXX = g++
CXXFLAGS = -Wall -std=c++17 -fopenmp -lgomp
LDLIBS = -lstdc++ -lm

all: out.txt
	cat out.txt

out.txt: main
	./$< > $@

main: main.cpp montecarlo.cpp montecarlo.hpp

clean:
	find . -type f -executable -delete
	$(RM) *.o out.* *.out