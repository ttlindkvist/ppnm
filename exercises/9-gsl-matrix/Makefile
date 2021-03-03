CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all clean

all: out.txt
	cat out.txt

out.txt: main
	./main > out.txt

clean:
	$(RM) *.o main *.txt