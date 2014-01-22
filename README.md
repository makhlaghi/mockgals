mockgals
=========

Make mock galaxy(s) in a C array and add noise to it.

About this tool:
------------
`mockgals` is a tool to make a pre-defined or random
set of mock galaxies in a simple `C` array. The unique
aspect of `mockgals` is integration of the central regions
of the mock profiles and the PSFs. If asked to make random
profiles, it will make random Sersic profiles (of galaxies) 
and convolve them with a Moffat PSF (beta=3, FWHM=3 pixels). 
An example random set of mock galaxies with 50 mock galaxies
placed randomly in the image can be seen below:

<img src="https://raw.github.com/makhlaghi/mockgals/master/jpgs/nonoise.jpg"
    width=300 />
<img src="https://raw.github.com/makhlaghi/mockgals/master/jpgs/withnoise.jpeg"
    width=300 />

