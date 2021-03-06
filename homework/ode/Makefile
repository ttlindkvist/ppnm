CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONY: all clean dyn_test

all: SHM.png SIR.png output.out 3body.png 4body.png
	cat output.out


dyn_test: dyn_matrix_test
	./$<
	
SHM.out SIR.out output.out nbody.out: main
	./$< 1> output.out

SHM.png: SHM.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "";\
	set key top left;\
	set tics out;\
	set xlabel "t";\
	set ylabel "x";\
	set xrange [0:20];\
	set yrange [-1.2:1.8];\
	plot \
	 "$<" using 1:2 with line color black title "position",\
	 "$<" using 1:3 with line color red title "velocity"\
	' | pyxplot
	
SIR.png: SIR.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "";\
	set key outside;\
	set tics out;\
	set xlabel "t (days)";\
	set ylabel "Population (\%) - (Tc, Tr)";\
	set xrange [0:100];\
	set yrange [0:1];\
	plot \
	 "$<" index 0 using 1:2 with line color black title "S (1, 14)",\
	 "$<" index 0 using 1:3 with line color red title "I",\
	 "$<" index 0 using 1:4 with line color green title "R",\
	 "$<" index 1 using 1:2 with line lt 2 color black title "S (2, 14)",\
	 "$<" index 1 using 1:3 with line lt 2 color red title "I",\
	 "$<" index 1 using 1:4 with line lt 2 color green title "R",\
	 "$<" index 2 using 1:2 with line lt 3 color black title "S (4, 14)",\
	 "$<" index 2 using 1:3 with line lt 3 color red title "I",\
	 "$<" index 2 using 1:4 with line lt 3 color green title "R"\
	' | pyxplot
	
3body.png: nbody.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "";\
	set key top left;\
	set tics out;\
	set xlabel "y";\
	set ylabel "x";\
	set xrange [-2:2];\
	set yrange [-1:1];\
	plot \
	 "$<" index 0 using 1:2 with line color black title "m1",\
	 "$<" index 0 using 4:5 with line color red title "m2",\
	 "$<" index 0 using 7:8 with line color green title "m3"\
	' | pyxplot
	
4body.png: nbody.out Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Almost atomic symbol";\
	set key top left;\
	set tics out;\
	set xlabel "y";\
	set ylabel "x";\
	set xrange [-3:3];\
	set yrange [-1.5:1.5];\
	plot \
	 "$<" index 1 using 1:2 with line color black title "m1",\
	 "$<" index 1 using 4:5 with line color red title "m2",\
	 "$<" index 1 using 7:8 with line color green title "m3",\
	 "$<" index 1 using 10:11 with line color blue title "m4"\
	' | pyxplot
main: main.o ode.c ode.h matrix.c matrix.h dyn_matrix.c dyn_matrix.h
dyn_matrix_test: dyn_matrix_test.o dyn_matrix.c dyn_matrix.h

clean:
	$(RM) main *.o *.out *.png dyn_matrix_test