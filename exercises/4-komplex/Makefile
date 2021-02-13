CFLAGS = -Wall -std=gnu99 -O
LDLIBS = -lm
.PHONEY: all clean
all: out.txt Makefile
	cat $<
out.txt: main
	./$< > $@
main: main.o komplex.o
main.o komplex.o: komplex.h
clean:
	$(RM) main *.o out.txt 
