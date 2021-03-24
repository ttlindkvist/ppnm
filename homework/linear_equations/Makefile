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
	set xrange [200:800];\
	set yrange [0:1000];\
	a=1;\
	f(x)=a*x**3;\
	fit f(x) "$<" using 1:2 via a;\
	b=1;\
	g(x)=b*x**3;\
	fit g(x) "$<" using 1:3 via b;\
	plot \
	 "$<" using 1:2 with line color red title "GS QR decomp",\
	 "$<" using 1:3 with line color black title "gsl-linalg-QR-decomp",\
	 f(x) with line lt 2 color black title "c1*N**3",\
	 g(x) with line lt 2 color black title "c2*N**3"\
	' | pyxplot
	
result.out timing.out: main
	@echo "This may take a while..."
	./$< > result.out

main: main.o matrix.o matrix.h

clean:
	$(RM) main *.o *.out *.png