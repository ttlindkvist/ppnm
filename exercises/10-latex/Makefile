CFLAGS = -O -Wall -Werror
LDLIBS = -lm

.PHONEY: all clean

all: report.pdf

report.pdf: report.tex plot1.png plot2.png plot3.png plot4.png plot5.png
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
	set xrange [0:100];\
	plot \
	 "$<" using 1:4 with line color red title "relative error"\
	' | pyxplot
plot3.png: out.dat Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set logscale y;\
	set xrange [-100:0];\
	plot \
	 "$<" using 1:4 with line color red title "without negative",\
	 "$<" using 1:5 with line color black title "with negative"\
	' | pyxplot
plot4.png: out.dat Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set logscale y;\
	set xrange [0:100];\
	plot \
	 "$<" using 1:4 with line lw 3 color red title "convoluted",\
	 "$<" using 1:6 with line color black title "backwards ordinary sum",\
	 "$<" using 1:7 with line color green title "forwards ordinary sum"\
	' | pyxplot
	
plot5.png: largeN.dat Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [0:100];\
	set yrange [0:2e-13];\
	set arrow 1 from 8,0 to 8,1 with nohead linetype 2 color black;\
	set arrow 2 from 16,0 to 16,1 with nohead linetype 2 color black;\
	set arrow 3 from 32,0 to 32,1 with nohead linetype 2 color black;\
	set arrow 4 from 64,0 to 64,1 with nohead linetype 2 color black;\
	set arrow 5 from 4,0 to 4,1 with nohead linetype 2 color black;\
	plot \
	 "$<" using 1:2 with line color red title "relative error"\
	' | pyxplot
	
out.dat largeN.dat: main
	./$<

main: main.o
main.o: main.c

clean:
	$(RM) main *.dat *.png *.o *.pdf *.log *.aux