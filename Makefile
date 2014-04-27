src=./src

objects=main.o mock.o stats.o attaavv.o raddist.o freqdomain.o \
        forqsort.o arraymanip.o fitsarrayvv.o pix.o ui.o

vpath %.h $(src)
vpath %.c $(src)

CC=gcc
CFLAGS=-Wall -W -O -pedantic -I$(src)
LDFLAGS=-L/usr/lib64/atlas/
LDLIBS=-lcfitsio -lfftw3 -pthread -lcblas -lgsl -lm

mockgals: $(objects) 
	rm *.fits *.txt ./mockgals
	@$(CC) -o mockgals $(LDFLAGS) $(objects) $(LDLIBS) 
	@rm *.o
	./mockgals -mnp

.SILENT: $(objects)
