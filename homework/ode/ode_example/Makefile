ODE=../ode
CFLAGS = -I$(ODE)
LDLIBS = -lm

all:Out.path.png
Out.path.png:out.path.txt Makefile
	echo '\
set term png ;\
set out "$@" ;\
set xlabel "x" ;\
set ylabel "y" ;\
set key out ;\
set title "Figure-8 solution to 3-body gravitational problem" ;\
plot \
 "$<" using 2:3 with lines title "body 1"\
,"$<" using 4:5 with lines title "body 2"\
,"$<" using 6:7 with lines title "body 3"\
	'|pyxplot

out.path.txt:nbody
	./$< 2> $@

nbody:nbody.o ode.o
nbody.o: ode.h
ode.o: ode.c ode.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<

clean:
	find . -type f -executable -delete
	$(RM) *.o [Oo]ut*
