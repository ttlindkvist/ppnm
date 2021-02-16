CFLAGS = -Wall -std=gnu99 -O
LDLIBS = -lm
.PHONEY: all clean
all: out.txt input.txt Makefile
	cat $<
out.txt: main
	./$< CMD $(shell cat input.txt) > $@
	./$< STDIN < input.txt >> $@
	./$< FILE input.txt $@

	
main: main.o
main.o: main.c
clean:
	$(RM) main *.o out.txt 