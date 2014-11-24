# Process this file with autoconf to produce a configure script.
#
# Hand written file: used as input into GNU Autotools (autoconf).
#
# This is the autoconf instructions to build MockGals. MockGals is a
# tool to make mock galaxies and stars in a FITS image from an input
# catalog. It is free software released under GNU GPLv3+.
#
# Copyright (C) 2014 Mohammad Akhlaghi
# Tohoku University Astronomical Institute, Sendai, Japan.
# http://astr.tohoku.ac.jp/~akhlaghi/

# Definitions:
AC_PREREQ([2.69])
AC_INIT([MockGals], [0.2], [makhlaghi@gmail.com], [mockgals],
        [http://www.astr.tohoku.ac.jp/~akhlaghi/mockgals/index.html])
AM_INIT_AUTOMAKE([-Wall -Werror gnu]) 
AC_CONFIG_SRCDIR([src/ui.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
: ${CFLAGS=""}
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [main], [], [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([pthread], [main], [], [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([cfitsio], [main], [], [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([fftw3f], [main], [], [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([gslcblas], [main], [], [AC_MSG_ERROR([Cannot continue.])], [-lm])
AC_CHECK_LIB([gsl], [main], [],[AC_MSG_ERROR([Cannot continue.])], [-lgslcblas -lm])

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h sys/time.h unistd.h argp.h],
                 [], [AC_MSG_ERROR([Cannot continue.])])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([gettimeofday pow sqrt strtol])

# Tell autoconf what to work on:
AC_CONFIG_FILES([Makefile
                 src/Makefile])

#Read the configuration date.
AC_MSG_CHECKING([for date of compilation])
AC_DEFINE_UNQUOTED([CONFIGDATE], ["`date +"%d %B %Y"`"], 
                   [Date of compilation])
AC_DEFINE_UNQUOTED([CONFIGTIME], ["`date +"%H:%M"`"], 
                   [Time of compilation])
AC_MSG_RESULT([done])

#Check if `defaults.ac` exists:
AC_MSG_CHECKING([if \`defaults.ac\` exists])
if [[ -f "$srcdir/defaults.ac" ]]; then
  source $srcdir/defaults.ac
else
  printf "\n$srcdir/defaults.ac does not exist!"; exit 1; fi
AC_MSG_RESULT([done])

#Check to see if the number of threads is defined in `default.ac`, if
#it is, use that value, if not, read the number of threads from the system.
if [[ -n "$numthreads" ]]; then
    AC_MSG_CHECKING([number of threads in \`defaults.ac\`])
    AC_DEFINE_UNQUOTED([DP_NUMTHREADS_T], ["$numthreads"], [Default parameter])
    AC_DEFINE_UNQUOTED([DP_NUMTHREADS_V], [$numthreads], [Default parameter])
    AC_MSG_RESULT([$numthreads])	  
else
    AC_MSG_CHECKING([for number of available threads to this system])
    AC_DEFINE_UNQUOTED([DP_NUMTHREADS_T], ["`nproc`"], 
                       [Number of threads available to system, text])
    AC_DEFINE_UNQUOTED([DP_NUMTHREADS_V], [`nproc`], 
                       [Number of threads available to system, value])
    AC_MSG_RESULT([$(nproc)])
fi

#Function to check if input is a number:
checkifnumber()
{
  if ! [[[ $1 =~ ^-?[0-9]+([.][0-9]+)?$ ]]] ; then
    printf "\n\n\"$1\" in \'defaults.ac\` is not a number.\n"
    printf "Note: Replace numbers written like \".1\" with \"0.1\".\n\n"
    exit 1
  fi
}

#Read the rest of the default values and check them:
AC_MSG_CHECKING([default values in \`defaults.ac\`])

AC_DEFINE_UNQUOTED([DP_PSF], ["$psf"], [Default parameter])

if [[ -n "$truncation" ]]; then
  checkifnumber $truncation
  if [[ `echo "$truncation>0" | bc` == 1 ]]; then
    AC_DEFINE_UNQUOTED([DP_TRUNCATION_T], ["$truncation"], [Default parameter])
    AC_DEFINE_UNQUOTED([DP_TRUNCATION_V], [$truncation], [Default parameter])
  else
    printf "\n\nIn \`defaults.ac\`, truncation must be >0.\n\n"; exit 1; fi
else printf "\n\ntruncation not in \`defaults.ac\` or not set.\n\n"; exit 1; fi

if [[ -n "$tolerance" ]]; then
  checkifnumber $tolerance
  if [[ `echo "$tolerance>=0 || $tolerance<=1" | bc` == 1 ]]; then
    AC_DEFINE_UNQUOTED([DP_TOLERANCE_T], ["$tolerance"], [Default parameter])
    AC_DEFINE_UNQUOTED([DP_TOLERANCE_V], [$tolerance], [Default parameter])
  else
    printf "\n\nIn \`defaults.ac\`, tolerance must be >=0 or <=1.\n\n"; exit 1; fi
else printf "\n\ntolerance not in \`defaults.ac\` or not set.\n\n"; exit 1; fi

if [[ -n "$psffunction" ]]; then
  if [[ $psffunction = "moffat" ]]; then
    AC_DEFINE_UNQUOTED([DP_PSFFUNCTION_T], ["$psffunction"], [Default parameter])  
    AC_DEFINE_UNQUOTED([DP_PSFFUNCTION_V], [1], [Default parameter])
  elif [[ $psffunction = "gaussian" ]]; then
    AC_DEFINE_UNQUOTED([DP_PSFFUNCTION_T], ["$psffunction"], [Default parameter])    
    AC_DEFINE_UNQUOTED([DP_PSFFUNCTION_V], [2], [Default parameter])
  else
    printf "\n\nIn \`defaults.ac\`, psffunction must be \"moffat\" or \"gaussian\".\n\n"
    exit 1; fi
else printf "\n\npsffunction not in \`defaults.ac\` or not set.\n\n"; exit 1; fi

if [[ -n "$psffwhm" ]]; then
  checkifnumber $psffwhm
  if [[ `echo "$psffwhm>0" | bc` == 1 ]]; then
    AC_DEFINE_UNQUOTED([DP_PSFFWHM_T], ["$psffwhm"], [Default parameter])
    AC_DEFINE_UNQUOTED([DP_PSFFWHM_V], [$psffwhm], [Default parameter])
  else
    printf "\n\nIn \`defaults.ac\`, psffwhm must be >0.\n\n"; exit 1; fi
else printf "\n\npsffwhm not in \`defaults.ac\` or not set.\n\n"; exit 1; fi

if [[ -n "$moffatbeta" ]]; then
  checkifnumber $moffatbeta
  if [[ `echo "$moffatbeta>0" | bc` == 1 ]]; then
    AC_DEFINE_UNQUOTED([DP_MOFFATBETA_T], ["$moffatbeta"], [Default parameter])
    AC_DEFINE_UNQUOTED([DP_MOFFATBETA_V], [$moffatbeta], [Default parameter])
  else
    printf "\n\nIn \`defaults.ac\`, moffatbeta must be >0.\n\n"; exit 1; fi
else printf "\n\nmoffatbeta not in \`defaults.ac\` or not set.\n\n"; exit 1; fi

if [[ -n "$psftrunc" ]]; then
  checkifnumber $psftrunc
  if [[ `echo "$psftrunc>0" | bc` == 1 ]]; then
    AC_DEFINE_UNQUOTED([DP_PSFTRUNC_T], ["$psftrunc"], [Default parameter])
    AC_DEFINE_UNQUOTED([DP_PSFTRUNC_V], [$psftrunc], [Default parameter])
  else
    printf "\n\nIn \`defaults.ac\`, psftrunc must be >0.\n\n"; exit 1; fi
else printf "\n\npsftrunc not in \`defaults.ac\` or not set.\n\n"; exit 1; fi

if [[ -n "$naxis1" ]]; then
  checkifnumber $naxis1
  if [[ `echo "$naxis1>0" | bc` == 1 ]]; then
    AC_DEFINE_UNQUOTED([DP_NAXIS1_T], ["$naxis1"], [Default parameter])
    AC_DEFINE_UNQUOTED([DP_NAXIS1_V], [$naxis1], [Default parameter])
  else
    printf "\n\nIn \`defaults.ac\`, naxis1 must be >0.\n\n"; exit 1; fi
else printf "\n\nnaxis1 not in \`defaults.ac\` or not set.\n\n"; exit 1; fi

if [[ -n "$naxis2" ]]; then
  checkifnumber $naxis2
  if [[ `echo "$naxis2>0" | bc` == 1 ]]; then
    AC_DEFINE_UNQUOTED([DP_NAXIS2_T], ["$naxis2"], [Default parameter])
    AC_DEFINE_UNQUOTED([DP_NAXIS2_V], [$naxis2], [Default parameter])
  else
    printf "\n\nIn \`defaults.ac\`, naxis2 must be >0.\n\n"; exit 1; fi
else printf "\n\nnaxis2 not in \`defaults.ac\` or not set.\n\n"; exit 1; fi

if [[ -n "$background" ]]; then
  checkifnumber $background
  if [[ `echo "$background>0" | bc` == 1 ]]; then
    AC_DEFINE_UNQUOTED([DP_BACKGROUND_T], ["$background"], [Default parameter])
    AC_DEFINE_UNQUOTED([DP_BACKGROUND_V], [$background], [Default parameter])
  else
    printf "\n\nIn \`defaults.ac\`, background must be >0.\n\n"; exit 1; fi
else printf "\n\nbackground not in \`defaults.ac\` or not set.\n\n"; exit 1; fi

if [[ -n "$zeropoint" ]]; then
  checkifnumber $zeropoint
  AC_DEFINE_UNQUOTED([DP_ZEROPOINT_T], ["$zeropoint"], [Default parameter])
  AC_DEFINE_UNQUOTED([DP_ZEROPOINT_V], [$zeropoint], [Default parameter])
else printf "\n\nzeropoint not in \`defaults.ac\` or not set.\n\n"; exit 1; fi

AC_MSG_RESULT([done])

AC_OUTPUT

# Print a bye-bye message.
echo
echo "==================================================="
echo "AC_PACKAGE_STRING is now configured:"
echo "  - Your defaults read from \`"defaults.ac"\`."
echo "  - Compilation configurations ready."
echo
echo "You are now ready to compile AC_PACKAGE_STRING".
echo "Simply run the commands: \`make\` then \`make install\`"
echo "You can run \`make check\` before installing."
echo "==================================================="
echo
