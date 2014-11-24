## Process this file with automake to produce Makefile.in
##
## Hand written file: used as input into Autotools.
##
## This file will be called by the Makefile.am in the top directory,
## so this file does not need to be called directly by automake to
## make the full MockGals package.
##
## This is the automake instructions to build MockGals. Using an input
## catalog of mock galaxies, MockGals will build an image of the
## galaxies in a FITS image. It is free software released under GNU
## GPLv3+.
##
## Copyright (C) 2014 Mohammad Akhlaghi
## Tohoku University Astronomical Institute, Sendai, Japan.
## http://astr.tohoku.ac.jp/~akhlaghi/

AM_CFLAGS = -Wall -O3
##AM_CFLAGS = -Wall -g -O0
AM_LDFLAGS = -lcfitsio -lfftw3f -pthread -lgsl -lgslcblas -lm

bin_PROGRAMS = mockgals
mockgals_SOURCES = arraymanip.c arraymanip.h attaavv.c attaavv.h	\
convolve.c convolve.h fitsarrayvv.c fitsarrayvv.h forqsort.c		\
forqsort.h integtwod.c integtwod.h macrofunctions.h main.c mock.c	\
mock.h pix.c pix.h profiles.c profiles.h raddist.c raddist.h rand.c	\
rand.h sll.c sll.h stats.c stats.h ui.c ui.h argpparser.h

## Necessary commands for `make check`:
TESTS = basic.sh

CLEANFILES = testinfo.fits testinfo_cat.txt 
