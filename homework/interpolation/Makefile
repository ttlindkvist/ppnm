CFLAGS = -O -Wall -Werror
LDLIBS = -lm

.PHONEY: all clean

all: linspline.png quadspline.png

linear.out quad.out: main *.in
	./main
	
main: main.o lspline.o qspline.o

linspline.png: linear.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [0:10];\
	set yrange [-4:20];\
	plot \
	 "$<" using 1:2 with line color red title "Linear spline through points",\
	 "$<" using 1:3 with line color black title "integrated spline",\
	 "lpoints.in" every ::1 using 1:2 with points title "points"\
	' | pyxplot
quadspline.png: quad.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [0:6];\
	set yrange [0:26];\
	plot \
	 "$<" using 1:2 with line color red title "Quad spline through points",\
	 "$<" using 1:3 with line color black title "derivative",\
	 "qpoints.in" every ::1 using 1:2 with points title "points"\
	' | pyxplot
clean:
	$(RM) *.out *.o main *.png