/*********************************************************************
MockGals - Make mock galaxies and stars from a catalog.

Copyright (C) 2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

MockGals is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MockGals is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MockGals. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef ARGPPARSER_H
#define ARGPPARSER_H

#include <argp.h>

#include "mock.h"
#include "ui.h"			/* Needs mock.h. */

/* Definition parameters for the argp: */
const char *argp_program_version=PACKAGE_STRING
  "\nCopyright (C) 2014 Mohammad Akhlaghi.\n"
  "License GPLv3+: GNU GPL version 3 or later "
  "<http://gnu.org/licenses/gpl.html>\n"
  "This is free software: you are free to change "
  "and redistribute it.\n"
  "There is NO WARRANTY, to the extent permitted by law.";





const char *argp_program_bug_address=PACKAGE_BUGREPORT;





static char args_doc[] = "[OutFITSimage.fits] catalog";





const char doc[] = 
  /* Before the list of options: */
  "\n"PACKAGE_STRING" -- Make mock galaxies and stars from a catalog. \n"
  "Configured for this machine on "CONFIGDATE", "CONFIGTIME".\n\n"
  "See \"Running "PACKAGE_NAME"\" in the official documentation for more "
  "detailed explanation. Official documentation can be found by running:"
  "`info "PACKAGE"`.\n"
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* The options are classified into these categories:
   1. Operating mode like quiet, help and version.
   2. Input (image, mask and kernel name, extensions and ...)
   3. Meshs and threads.
   4. Detection 

   Available letters for short options:
   a c d e g h i j k m n u s v
   A E F G H I J L M P Q R S T U W X Y Z

   Number keys used: <=502.

   Options with keys (second structure element) larger than 500 do not
   have a short version.
 */
static struct argp_option options[] =
  {
    /* Such cases are group headers */
    {
      0, 0, 0, 0,  		/* These have to be zero for groups. */
      "Operating modes:", 	/* Explanation for the group. */
      -1			/* Group ID. */
    },
    {
      "quiet",		      /* Long name for this option.         */
      'q',		      /* Short name or key for this option. */
      0,		      /* Informative type of value it gets. */
      OPTION_ARG_OPTIONAL,    /* Flags for this option. */
      "Only report errors.",
      -1		      /* Option group ID. */
    },
    {
      "numthreads",
      'N',
      "INT",
      0,
      "["DP_NUMTHREADS_T"] The number of CPU threads to use.",
      -1
    },
    {
      "keepnamedir",
      'K',
      0,
      0,
      "Keep the output name directory information",
      -1
    },
    {
      "dontdelete",
      'D',
      0,
      0,
      "Don't delete an existing output, exit.",
      -1
    },



    {
      0, 0, 0, 0,
      "Galaxy profiles:",
      1		
    },
    {
      "truncation",
      't',
      "FLT",
      0,
      "["DP_TRUNCATION_T"] Truncation distance, multiple of radius.",
      1
    },
    {
      "tolerance",
      'l',
      "FLT",
      0,
      "["DP_TOLERANCE_T"] When to switch to less accurate method.",
      1
    },





    {
      0, 0, 0, 0,
      "Point Spread function",
      2
    },
    {
      "psf",
      'p',
      "STR",
      0,
      "["DP_PSF"] File name of input PSF FITS image.",
      2
    },
    {
      "psffunction",
      'f',
      "STR",
      0,
      "["DP_PSFFUNCTION_T"] PSF function: `moffat` or `gaussian`.",
      2
    },
    {
      "psffwhm",
      'w',
      "FLT",
      0,
      "["DP_PSFFWHM_T"] FWHM of PSF in units of pixels.",
      2
    },
    {
      "moffatbeta",
      'B',
      "FLT",
      0,
      "["DP_MOFFATBETA_T"] Moffat function's beta value.",
      2
    },    
    {
      "psftrunc",
      'r',
      "FLT",
      0,
      "["DP_PSFTRUNC_T"] PSF truncation in units of FWHM/2.",
      2
    },
    {
      "savepsf",
      500,
      0,
      0,
      "Save the PSF used in `psf.fits`.",
      2
    },
    {
      "onlypsf",
      501,
      0,
      0,
      "Only make the PSF, no galaxies.",
      2
    },
    {
      "sepsf",
      502,
      0,
      0,
      "Save PSF for input into SExtractor.",
      2
    },




    {
      0, 0, 0, 0,
      "Outputs:",
      3
    },
    {
      "naxis1",
      'x',
      "INT",
      0,
      "["DP_NAXIS1_T"] Number of pixels along first FITS axis.",
      3
    },    
    {
      "naxis2",
      'y',
      "INT",
      0,
      "["DP_NAXIS2_T"] Number of pixels along second FITS axis.",
      3
    },
    {
      "output",
      'o',
      "STR",
      0,
      "Output name in any format",
      3
    },
    {
      "background",
      'b',
      "FLT",
      0,
      "["DP_BACKGROUND_T"] Image background (amplitude of noise).",
      3
    },
    {
      "zeropoint",
      'z',
      "FLT",
      0,
      "["DP_ZEROPOINT_T"] Zeropoint magnitude of the image.",
      3
    },
    {
      "viewnoconv",
      'O',
      0,
      0,
      "Keep the unconvolved image as an extension.",
      2
    },
    {
      "viewconv",
      'C',
      0,
      0,
      "Keep the convolved image as an extension.",
      2
    },

    


    {0}
  };





/* Parse a single option: */
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{ 
  /* Save the arguments structure: */
  struct mockparams *p = state->input;
  
  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then `arg`
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if(arg && arg[0]=='=')
    argp_error(state, "Incorrect use of the equal sign (`=`). For short "
	       "options, `=` should not be used and for long options, "
	       "there should be no space between the option, equal sign "
	       "and value.");
 
  switch(key)
    {
    /* Profile parameters:  */
    case 't':
      floatl0(arg, &p->truncation, "truncation", key);
      break;
    case 'l':
      floatel0(arg, &p->tolerance, "tolerance", key);
      break;



      

    /* PSF: */
    case 'p':
      p->up.psfname=arg;
      break;
    case 'f':
      if(strcmp(arg, "moffat")==0)
	p->psffunction=1;
      else if(strcmp(arg, "gaussian"))
	p->psffunction=2;
      else
	argp_error(state, "The value of the `--psffunction` (`-f`) option "
		   "should be either `moffat` or `gaussian`.");
      break;
    case 'w':
      floatl0(arg, &p->psf_p1, "psffwhm", key);
      break;
    case 'B':
      floatl0(arg, &p->psf_p2, "moffatbeta", key);
      break;
    case 'r':
      floatl0(arg, &p->psf_t, "psftrunc", key);
      break;
    case 500:
      p->outpsfname="junk";
      break;
    case 501:
      p->onlypsf=1;
      break;
    case 502:
      p->sepsfname="junk";
      break;


      

    /* Output parameters */
    case 'x':
      sizetlzero(arg, &p->s1, "naxis1", key);
      break;
    case 'y':
      sizetlzero(arg, &p->s0, "naxis2", key);
      break;	
    case 'o':
      if(p->up.outname)
	argp_error(state, "Only one output name can be given.");
      else
	p->up.outname=arg;
      break;
    case 'b':
      floatel0(arg, &p->background, "background", key);
      break;
    case 'z':
      anyfloat(arg, &p->zeropoint, "zeropoint", key);
      break;
    case 'C':
      p->viewconv=1;
      break;
    case 'O':
      p->viewnoconv=1;
      break;

      
      

    /* Operating modes: */
    case 'q':
      p->verb=0;
      break;
    case 'N':
      fprintf(stderr, PACKAGE": Warning:`--numthreads` (`-N`) is "
	      "not yet activated, one thread will be used.\n");
      break;
    case 'K':
      p->up.removenamedir=0;
      break;
    case 'D':
      p->up.dontdelete=1;
      break;



      
      /* Read the non-option argument: */
    case ARGP_KEY_ARG:
      /* Since there are only two arguments, state->arg_num should
	 never be more than 2. Note that it starts from zero.*/
      if(state->arg_num >= 2)
	argp_error(state, "Too many arguments! Arguments can either be an "
		   "optional `FITSimage.fits` (output name) or a mandatory "
		   "catalog (text file with any extension).");

      /* See what type of input value it is and put it in. */
      if( strcmp(&arg[strlen(arg)-5], ".fits")==0 )
	{
	  if(p->up.outname)
	    argp_error(state, "Only one output name can be given.");
	  else
	    p->up.outname=arg;
	}
      else
        p->up.infoname=arg;
      break;
      
      /* Make sure an argument is given: */
    case ARGP_KEY_END:
      if(state->arg_num==0)
	argp_error(state, "No argument given!");
      if(p->up.infoname==NULL)
	argp_error(state, "No input catalog specified!");
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}





/* Basic structure defining the whole argument reading process. */
static struct argp argp = {options, parse_opt, args_doc, doc};


#endif
