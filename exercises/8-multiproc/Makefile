CFLAGS += -O -Wall 
LDLIBS += -lm -pthread -lpthread

.PHONEY: all clean

all: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

main: main.o
main.o: main.c

clean:
	$(RM) out.txt *.o main
