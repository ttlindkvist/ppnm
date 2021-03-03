# CFLAGS = -o -w -std=gnu99
CFLAGS = -O -Wall -std=gnu99
# CFLAGS += 'gsl-config --cflags'
CFLAGS += `gsl-config --cflags`

# LDLIBS += $(gsl-config --libs)
LDLIBS += $(shell gsl-config --libs)

# all: output ; cat $@
all: output ; cat $<
# output: main : $@ > $<
output: main ; ./$< > $@

# clean: ; $RM main output
clean: ; $(RM) main output
