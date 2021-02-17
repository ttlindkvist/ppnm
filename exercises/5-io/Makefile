CFLAGS = -Wall -std=gnu11 -O
LDLIBS = -lm
.PHONEY: all clean
all: out.txt input.txt Makefile
	cat $<
out.txt: main
	./$< CMD $(shell cat input.txt) > $@
	./$< STDIN < input.txt >> $@
	./$< FILE input.txt $@ #appends output to existing out.txt 
	
main: main.o
main.o: main.c
clean:
	$(RM) main *.o out.txt 