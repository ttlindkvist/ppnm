CFLAGS = -O -Wall -Werror
LDLIBS = -lm

.PHONEY: all clean

all: plot1.png plot2.png report.pdf

report.pdf: report.tex plot1.png plot2.png
	pdflatex $<
	pdflatex $<

plot1.png: out.dat Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set logscale y;\
	set title "Exponential function";\
	plot \
	 "$<" using 1:2 with line color red title "Our exp(x)",\
	 "$<" using 1:3 with line color black title "exp(x) from math.h"\
	' | pyxplot
plot2.png: out.dat Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "";\
	set xrange [0:100];\
	plot \
	 "$<" using 1:4 with line color red title "relative error"\
	' | pyxplot
out.dat: main
	./$<

main: main.o
main.o: main.c

clean:
	$(RM) main *.dat *.png *.o *.pdf