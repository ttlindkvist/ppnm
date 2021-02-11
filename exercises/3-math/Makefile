CC = gcc
CFLAGS = -O -std=gnu11 -Wall
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: math
	./math > out.txt

math: math.o
	$(CC) -o math math.o $(LDLIBS)

math.o: math.c
	$(CC) $(CFLAGS) -c math.c

.PHONEY: clean
clean:
	rm -rf math.o ./math out.txt

.PHONEY: test
test:
	echo $(LDLIBS)
	echo $(CC)
