CC = gcc
CFLAGS = -O -std=gnu11
LDLIBS = -lm

out.txt: hello
	./hello > out.txt

hello: hello.o
	$(CC) -o hello hello.o $(LDLIBS)

hello.o: hello.c
	$(CC) $(CFLAGS) -c hello.c

clean:
	rm -f hello.o hello out.txt
