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
	set yrange [0:700];\
	a=1;\
	b=1;\
	f(x)=a*x**b;\
	fit f(x) "$<" using 1:2 via a,b;\
	c=1;\
	d=1;\
	g(x)=c*x**d;\
	fit g(x) "$<" using 1:3 via c,d;\
	plot \
	 "$<" using 1:2 with points color red title "GS QR decomp",\
	 "$<" using 1:3 with points color black title "gsl-linalg-QR-decomp",\
	 f(x) with line lt 2 color black title "$$\\alpha N^{%.2f}$$"%(b),\
	 g(x) with line lt 2 color black title "$$\\beta N^{%.2f}$$"%(d)\
	' | pyxplot
	
result.out timing.out: main
	@echo "This may take a while..."
	./$< > result.out

main: main.o matrix.o matrix.h

clean:
	$(RM) main *.o *.out *.png