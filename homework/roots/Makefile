CXX = g++
CXXFLAGS = -Wall -std=c++17 -Wno-sign-compare
CXXFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lstdc++
LDLIBS += $(shell gsl-config --libs)

all: out.txt 
	cat out.txt

out.txt error.out: main
	./$< > out.txt

main: main.cpp roots.cpp roots.hpp matrix.c matrix.h function.cpp function.hpp ode.h ode.c dyn_matrix.c dyn_matrix.h

clean:
	find . -type f -executable -delete
	$(RM) *.o out.* *.out