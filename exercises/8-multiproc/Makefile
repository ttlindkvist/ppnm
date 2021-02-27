CFLAGS = -O -Wall -std=gnu11
CFLAGS += -pthread
CFLAGS += -fopenmp -lgomp
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm -pthread
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all clean

all: out.txt time.txt convergence.png
	cat out.txt
	cat time.txt

out.txt time.txt: pthreads_mc_pi_new single_thread openMP_mc_pi Makefile
	time ./single_thread 10000000 1> out.txt 2> time.txt 
	time ./pthreads_mc_pi_new 4 2500000 1>> out.txt 2>> time.txt
	time ./openMP_mc_pi 4 2500000 1>> out.txt 2>> time.txt

convergence.png: convergence.dat Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom left;\
	set tics out;\
	set logscale x y;\
	set xlabel "N";\
	set ylabel "relative error";\
	set title "Convergence of Monte Carlo";\
	plot \
	 "$<" using 1:2 with points pointsize 0.4 title "rand-r",\
	 "$<" using 1:3 with points pointsize 0.4 color red title "Low discrepancy",\
	 "$<" using 1:4 with line color black title "$$\\frac{1}{\\sqrt{N}}$$"\
	' | pyxplot

convergence.dat: convergence
	./convergence 4

pthreads_mc_pi_new: pthreads_mc_pi_new.c
single_thread: single_thread.c
openMP_mc_pi: openMP_mc_pi.c

clean:
	$(RM) *.txt *.o pthreads_mc_pi_new single_thread openMP_mc_pi convergence *.png *.dat