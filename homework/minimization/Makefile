CXX = g++
CXXFLAGS = -Wall -std=c++17
CXXFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lstdc++
LDLIBS += $(shell gsl-config --libs)

all: out.txt 
	cat out.txt

out.txt : main
	./$< > out.txt

main: main.cpp minimization.cpp minimization.hpp

clean:
	find . -type f -executable -delete
	$(RM) *.o out.* *.out