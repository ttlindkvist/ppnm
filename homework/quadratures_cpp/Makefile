CXX = g++
CXXFLAGS = -Wall -std=c++17
LDLIBS = -lstdc++

all: out.txt
	cat out.txt

out.txt: main
	./$< > $@

main: main.cpp integration.cpp integration.hpp

clean:
	find . -type f -executable -delete
	$(RM) *.o out.* *.out