mockgals
=========

Make mock galaxy(s) in a C array or FITS image and add noise to it.

About this tool:
------------
`mockgals` is a tool to make a pre-defined or random
set of mock galaxies in a simple `C` array or a FITS image.
As an option, this program will also produce histograms of
the no noised and noised images.
The unique aspect of `mockgals` is integration of the central 
regions of the mock profiles and the PSFs. 
If asked to make random profiles, it will make random Sersic 
profiles (of galaxies) and convolve them with a Moffat PSF 
(beta=3, FWHM=3 pixels). 
An example random set of mock galaxies with 50 mock galaxies
placed randomly in the image can be seen below:

<img src="https://raw.github.com/makhlaghi/mockgals/master/jpgs/nonoise.jpg"
    width=300 />
<img src="https://raw.github.com/makhlaghi/mockgals/master/jpgs/withnoise.jpeg"
    width=300 />
    
The noise for every pixel is a random value taken from a 
Gaussian distribution with sigma=sqrt(sky+pixel value) and
mean of zero (which you can think of as sky subtracted). 


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
To install this `mockgals`, after downloading or cloning it,
all you have to do is to run `make` in the downloaded directory.
To run it with default (50 random Sersic profiles) you just have
to run: `./mockgals`. A FITS image titled `mock.fits` will be created
along with a text file listing all the parameters of the mock profiles.

###Set number of random galaxies
In case you want a certain number of random mock galaxies, 
then simply add that number as an option to `./mockgals`, for example
for 100 mock galaxies you can run: `./mockgals 100`.

The random values are taken from a uniform distribution, with the range 
of parameters that are defined in the function `setprflprms()` that is 
defined in `./src/mock.c`. There are three headers in `mock.h` that will 
facilitate particular cases, their names are fairly descriptive:

- `ONLYONEPROFILE`: if `1`, then only one profile with the parameters
  set in `setprflprms()` will be created.
- `ONGRID`: if `1` then 25 mock profiles are created on a grid.
  Their parameters can be set in `setprflprms()`.


###Input file
You can alternatively define an ASCII file as input into the program.
In this case the `mockgals` will read the number and properties of 
the mock galaxies from this table. In short, it has to have the same
number of columns as the `mock.txt` file generated for random mock
galaxies with the same column definitions:

1. ID (Won't be used!)
2. 0: Sersic, 1: Moffat, 2: Gaussian
3. X position (FITS standard, not C)
4. Y position (FITS standard, not C)
5. Sersic n, Moffat beta or Gaussian sigma.
6. Sersic r_e, Moffat FWHM or won't be used in Gaussian 
7. Position angle.
8. Axis ratio.
9. Signal to noise (average profile flux to sqrt(sky))
10. Total flux (Won't be used)

For the cases that "won't be used" you can just define a zero 
in the input file.

###Histogram
If you activate the header macro: `MOCKHIST`, then another output
text file will be created, a 3 column data file showing the left
histogram bin value in the first column, the number of pixels in 
that bin in the no noised image on the second column and the number
of pixels in that bin for the noised image. You can use this
as input into any bar plotting program you would like to generate
a histogram.

###Setting other parameters:
Currently the only way to change other paramters, e.g., the PSF 
used to convolve the image or the output image size, is to modify 
the `./scr/main.c` file. The necessary parameters are nicely named
to define their use.

###Viewing Multi extension FITS files:
I recommend `ds9` to view your FITS files, ordinarily opening a
FITS file in `ds9` will only open the first extension.
In order to view all extensions you can open `ds9`, then:
`file` -> `open other` -> `Open Multi Ext Cube...` or 
`Open Multi Ext Multi Frames` and then choose the file.
You can alternatively open `ds9` on the command line with
the option: `-medatacube`, for example: `ds9 -medatacube mock.fits`.
This second method is the most convenient.

Future updates:
------------
1. Work on a better user experience.
2. Fix any bugs I have not found yet!

Comments and suggestions:
----------------------------------------
I hope `mockgals` will be useful for you. If you find any problems in 
this program please contact me so I can correct them. I would also be 
very glad to hear any suggestions or comments you might have, thank you.

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
