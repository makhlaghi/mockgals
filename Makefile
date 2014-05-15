src=./src

objects=main.o mock.o stats.o attaavv.o raddist.o convolve.o   \
	forqsort.o rand.o integtwod.o profiles.o sll.o         \
	arraymanip.o fitsarrayvv.o pix.o ui.o

vpath %.h $(src)
vpath %.c $(src)

CC=gcc
CFLAGS=-Wall -W -O -pedantic -I$(src)
LDLIBS=-lcfitsio -lfftw3f -pthread -lgsl -lgslcblas -lm

mockgals: $(objects) 
	@$(CC) -o mockgals $(LDFLAGS) $(objects) $(LDLIBS) 
	@rm *.o
#	./mockgals -mne

.SILENT: $(objects)

.PHONY: install

install:
	cp ./mockgals /usr/local/bin/
