CFLAGS = -O -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
LDLIBS =  -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all clean

all: output.out fit_plot.png
	cat $<

fit_plot.png: fit.out data.in Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Best fit";\
	set key top right;\
	set tics out;\
	set xlabel "Time [days]";\
	set ylabel "Relative activity";\
	set xrange [0:16];\
	set yrange [9:143];\
	plot \
	 "$<" index 0 using 1:2 with line color blue title "Best fit",\
	 "$<" index 1 using 1:2 with line color green title "$$c+\delta c$$",\
	 "$<" index 2 using 1:2 with line color red title "$$c-\delta c$$",\
	 "data.in" every ::1 using 1:2:($$2/20.0) with errorbars title "Data"\
	' | pyxplot

output.out fit.out: main data.in
	./main > output.out

main: main.o matrix.o matrix.h

clean:
	$(RM) main *.o *.out