CFLAGS = -O -Wall -Werror
LDLIBS = -lm

.PHONEY: all clean

all: linspline.png

linear.out: main points.in
	./main
	
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
	 "points.in" every ::1 using 1:2 with points title "points"\
	' | pyxplot
	
clean:
	$(RM) *.out *.o main *.png