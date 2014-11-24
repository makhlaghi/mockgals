#! /bin/bash
#
# Tests (examples) of using MockGals:
#
#
#
# General information:
# ====================
#
# testinfo.txt is a short catalog of 5 stars and 50 galaxies randomly
# positioned. We will run MockGals to make the 50 profiles mentioned
# in this catalog.
#
# Since MockGals is not yet installed on the system, it is running
# from this directory, hence the `./mockgals`. After installation,
# run it with "mockgals" in any directory.
#
# After installation, see `mockgals --help` for all options.
#
#
# This test:
# -----------
#
# This is the most basic usage of MockGals. Two files will be output:
# a FITS image and a text image with the 50 profiles placed in the
# noisy image. If you want to see how the image would be prior to
# convolution or prior to adding noise, add the `--viewnoconv` and
# `--viewconv` options.

./mockgals $srcdir/testinfo.txt
