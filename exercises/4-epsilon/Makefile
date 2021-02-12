CFLAGS = -Wall -std=gnu11 -fwrapv -O 
LDLIBS = -lm
.PHONEY: all clean
all: out.txt Makefile
	cat $<
out.txt: main
	./$< > $@
main: main.o
main.o:
clean:
	$(RM) main *.o out.txt 