mockgals
=========

Make mock astronomical profiles (galaxy, star, ...) in a FITS file.

Description:
------------

`mockgals` is a tool to make a pre-defined or random set of mock
galaxies in a simple `C` array or a FITS image.  The unique aspect of
`mockgals` is integration of the central regions of the mock profiles
and the PSFs.  If asked to make random profiles (no options are
input), it will make random Sersic profiles (of galaxies) and convolve
them with a Moffat PSF. As an option, `mockgals` will also produce
histograms of the no noised and noised images.  An example random set
of mock galaxies with 50 mock galaxies placed randomly in the image
can be seen below. Left, image prior to convolution with with PSF.
Middle, after convolving with a FWHM of 5 pixels PSF. Right: After
adding noise.

<img src="https://raw.github.com/makhlaghi/mockgals/master/doc/mockgals-figures/s1_noconv.jpg" />
<img src="https://raw.github.com/makhlaghi/mockgals/master/doc/mockgals-figures/s1_conv.jpg" />
<img src="https://raw.github.com/makhlaghi/mockgals/master/doc/mockgals-figures/s1_noised.jpg" />
    
The noise for every pixel is a random value taken from a Gaussian
distribution with sigma=sqrt(sky+pixel value) and mean=sky+pixel value.


Prerequisits:
------------
`mockgals` requires several packages to be installed on your
machine, installing them is easy and straightforward. 

- [GNU Scientific Library](http://www.gnu.org/software/gsl/).
  For integration and random number generation.
- [FFTW](http://www.fftw.org/) 3.
  For convolution
- [cfitsio](http://heasarc.nasa.gov/fitsio/fitsio.html)
  For reading and writing FITS files.
- [ds9](http://ds9.si.edu/site/Home.html)
  To view the resulting Multiextension FITS files.

Installing and running:
------------

To install this `mockgals`, after downloading or cloning it, all you
have to do is to run `make` in the downloaded directory.  To run it
with default (45 random Sersic profiles and 5 random stars) you just
have to run: `./mockgals`. Some command line options can be given so
you can customize the output, to learn them, run `./mockgals -h`.  A
full list of all the options will be provided, nealy all the operation
of `mockgals` can be defined by these input options and their
arguments. In the future long arguments and an configuration file will
also be provided for user customization. The [POSIX argument syntax 
conventions](http://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html#Argument-Syntax) apply.

The output is a possibly multi extension FITS file (depending on what
you ask for in the options).


###Input file

You can define an ASCII file as input into the program.  In this case
the `mockgals` will read the number and properties of the mock
galaxies from this table. In short, it has to have the same number of
columns as the `mock.txt` file generated for random mock galaxies with
the same column definitions:

1. ID (Won't be used!)
2. 0: Sersic, 1: Moffat, 2: Gaussian, 3: Point source
3. X position (FITS standard, not C)
4. Y position (FITS standard, not C)
5. Sersic n, Moffat beta or Gaussian sigma.
6. Sersic r_e, Moffat FWHM or won't be used in Gaussian 
7. Position angle.
8. Axis ratio.
9. Signal to noise (average profile flux to sqrt(sky))
10. Total flux (Won't be used)

For the cases that "won't be used" you can just define a zero in the
input file. The last column will be filled once `mockgals` is finished.

###Histogram

If you want a histogram of the image before and after noise is added
you can do so with the `-t`, `-c` and `-d` arguments to the
program. By default the first is zero, this means that no histogram
will be made. But if you give it any positive value, it will make a
histogram in the range of the values you give to `-c` and `-d`, with
the number of bins between the two being the value of `-t`.

###Viewing Multi extension FITS files:

I recommend `ds9` to view your FITS files, ordinarily opening a FITS
file in `ds9` will only open the first extension.  In order to view
all extensions you can open `ds9`, then: `file` -> `open other` ->
`Open Multi Ext Cube...` or `Open Multi Ext Multi Frames` and then
choose the file.  You can alternatively open `ds9` on the command line
with the option: `-medatacube`, for example: `ds9 -medatacube
mock.fits`.  This second method is the most convenient. On a linux
system you can even set this command to be run when you double click
only any fits file, if it has several extensions you can go through them
and if it doesn't, it will open just like normal (this is how I have 
personally configured my linux system).

Future updates:
------------
0. Make a nice documentation.
1. Add more functionality:
   - Positioning of stars.
2. Work on a better user experience.
3. Fix any bugs I have not found yet!

Comments and suggestions:
----------------------------------------

I hope `mockgals` will be useful for you. If you find any problems in
this program please contact me so I can correct them. I would also be
very glad to hear any suggestions or comments you might have, thank
you.

makhlaghi@gmail.com 

akhlaghi@astr.tohoku.ac.jp

http://astr.tohoku.ac.jp/~akhlaghi/

----------------------------------------
Copyright:
----------------------------------------
Copyright (C) 2014 Mohammad Akhlaghi

Tohoku University Astronomical Institute

http://astr.tohoku.ac.jp/~akhlaghi/

`mockgals` is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

`mockgals` is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with `mockgals`.  If not, see <http://www.gnu.org/licenses/>.
