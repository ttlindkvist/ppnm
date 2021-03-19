CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all clean

all: result.out timing.png
	cat $<

timing.png: timing.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Timing";\
	set key top left;\
	set tics out;\
	set xlabel "N";\
	set ylabel "Time [ms]";\
	set xrange [0:500];\
	set yrange [0:300];\
	plot \
	 "$<" using 1:2 with line color red title "GS QR decomp",\
	 "$<" using 1:3 with line color black title "GSL HH QR decomp"\
	' | pyxplot
	
result.out timing.out: main
	@echo "This may take a while..."
	./$< > result.out

main: main.o matrix.o matrix.h

clean:
	$(RM) main *.o *.out *.png