CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

FIGDIR = figs
OBJDIR = obj

.PHONEY: all clean

all: linear.png quadratic.png cubic.png
	mkdir -p $(FIGDIR)
	mv $^ $(FIGDIR)/

linear.out quad.out cubic.out: main *.in
	./main
	
main: main.o lspline.o qspline.o cspline.o

linear.png: linear.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Linear spline through points";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [0:10];\
	set yrange [-4:20];\
	plot \
	 "$<" using 1:2 with line color red title "Linear spline through points",\
	 "$<" using 1:3 with line color black title "integrated spline",\
	 "$<" using 1:4 with line color red title "gsl linear spline",\
	 "$<" using 1:5 with line color black title "gsl integrated",\
	 "lpoints.in" every ::1 using 1:2 with points title "points"\
	' | pyxplot
quadratic.png: quad.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Quadratic spline through points";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [0:6];\
	set yrange [-6:14];\
	plot \
	 "$<" using 1:2 with line color red title "Quad spline through points",\
	 "$<" using 1:3 with line color black title "derivative",\
	 "$<" using 1:4 with line color green title "integrated",\
	 "qpoints.in" every ::1 using 1:2 with points title "points"\
	' | pyxplot
	
cubic.png: cubic.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Cubic spline through points";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [0:20];\
	set yrange [-3:30];\
	plot \
	 "$<" index 0 using 1:2 with points title "",\
	 "$<" index 1 using 1:2 with line color black title "$$ f(x) = x + 10\sin^5(2x/\pi) $$",\
	 "$<" index 1 using 1:3 with line color red title "Cubic spline"\
	' | pyxplot
	
clean:
	$(RM) *.out *.o main $(FIGDIR)/*.png