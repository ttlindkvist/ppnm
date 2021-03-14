#CFLAGS = `gsl-config --cflags`
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)
.PHONEY: all clean

all: erf.gnuplot.png erf.pyxplot.png gamma.pyxplot.png gamma.gnuplot.png lgamma.pyxplot.png lgamma.gnuplot.png cgamma.gnuplot.png

erf.pyxplot.png: erf.dat Makefile tabulated.erf.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	,"tabulated.erf.txt" using 1:2 with points title "tabulated data"\
	' | tee erf.plot.gpi | pyxplot

gamma.pyxplot.png: gamma.dat Makefile tabulated.gamma.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [-5:5];\
	set yrange [-5:26];\
	set title "Gamma function";\
	plot \
	"$<" using 1:2 with line title "tgamma from math.h",\
	"$<" using 1:3 with line title "gsl-sf-gamma",\
	"$<" using 1:4 with line title "homemade Gamma",\
	"tabulated.gamma.txt" using 1:2 with points title "tabulated data"\
	' | tee gamma.plot.gpi | pyxplot

lgamma.pyxplot.png: lgamma.dat Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set xrange [0:10];\
	set yrange [-1:15];\
	set title "log of the gamma function";\
	plot \
	"$<" using 1:2 with line title "lgamma from math.h",\
	"$<" using 1:3 with line title "gsl-sf-lgamma",\
	"$<" using 1:4 with line title "homemade LGamma"\
	' | tee lgamma.plot.gpi | pyxplot

cgamma.gnuplot.png: cgamma.gpi cgamma.dat 
	gnuplot $<

cgamma.gpi: Makefile
	>$@
	echo 'set term png size 1280,960 \' >> $@
	echo 'font arial 28 background rgb "white"' >> $@
	echo 'set out "cgamma.gnuplot.png"'  >> $@
	echo 'set view ,310' >>$@
	echo 'set title "|Γ(z)|"' >>$@
	echo 'set hidden3d' >>$@
	echo 'set pm3d' >>$@
	echo 'set grid' >>$@
	echo 'set xlabel "Re(z)"' >>$@
	echo 'set ylabel "Im(z)"' >>$@
	echo '#set zlabel "|Gamma(z)|" rotate parallel' >>$@
	echo 'set contour both' >>$@
	echo 'set cntrlabel onecolor' >>$@
	echo 'set linetype 1 linecolor rgb "black"' >>$@
	echo 'zmax=6' >>$@
	echo 'set cntrparam levels incremental 0.0,0.25,zmax' >>$@
	echo 'set palette rgb 33,13,10' >>$@
	echo 'cut(z)=z>zmax?zmax:z' >>$@
	echo 'splot [][][:zmax]\' >>$@
	echo ' "cgamma.dat" using 1:2:(cut($$3)) with lines notitle' >> $@

erf.plot.gpi: erf.pyxplot.png
gamma.plot.gpi: gamma.pyxplot.png
lgamma.plot.gpi: lgamma.pyxplot.png


erf.gnuplot.png: erf.plot.gpi
	cat $< | sed 's/pyxplot/gnuplot/' | gnuplot

gamma.gnuplot.png: gamma.plot.gpi
	cat $< | sed 's/pyxplot/gnuplot/' | gnuplot

lgamma.gnuplot.png: lgamma.plot.gpi
	cat $< | sed 's/pyxplot/gnuplot/' | gnuplot

erf.dat gamma.dat lgamma.dat: main
	./$<

main: main.o

clean:
	$(RM) *.o main *.png *.dat *.gpi 
test:
	echo $(CFLAGS)
	echo $(LDLIBS)