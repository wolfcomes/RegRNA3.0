/* esl_config.h.  Generated from esl_config.h.in by configure.  */
/* esl_config.h.in  [input to configure]
 * 
 * System-dependent configuration of Easel, by autoconf.
 * 
 * This file should be included in all Easel .c files before
 * anything else, because it may set #define's that control
 * behaviour of system includes and system libraries. An example
 * is large file support.
 * 
 * SVN $Id$
 * SRE, Fri Mar  3 08:03:32 2006 [St. Louis]
 */
#ifndef ESL_CONFIG_INCLUDED
#define ESL_CONFIG_INCLUDED

/* Version info.
 */
#define EASEL_VERSION "0.1"
#define EASEL_DATE "March 2007"
#define EASEL_COPYRIGHT "Copyright (C) 2004-2007 HHMI Janelia Farm Research Campus"
#define EASEL_LICENSE "Freely distributed under the Janelia Software License."

/* Large file support
 * Must precede any header file inclusion.
 */
#define _FILE_OFFSET_BITS 64
/* #undef _LARGE_FILES */
#define _LARGEFILE_SOURCE 1

/* Debugging verbosity (0=none;3=most verbose)
 */
/* #undef eslDEBUGLEVEL */

/* System headers
 */
#define HAVE_UNISTD_H 1
#define HAVE_STDINT_H 1
#define HAVE_INTTYPES_H 1

/* Types
 */
/* #undef WORDS_BIGENDIAN */
/* #undef uint8_t */
/* #undef uint16_t */
/* #undef uint32_t */
/* #undef uint64_t */
/* #undef off_t */

/* Optional packages
 */
/* #undef HAVE_LIBGSL */

/* Functions
 */
#define HAVE_MKSTEMP 1
#define HAVE_POPEN 1
#define HAVE_STRCASECMP 1
#define HAVE_TIMES 1
#define HAVE_FSEEKO 1

/*****************************************************************
 * Available augmentations.
 * 
 * If you grab a single module from Easel to use it by itself,
 * leave all these #undef'd; you have no augmentations.
 * 
 * If you grab additional Easel .c files, you can enable any
 * augmentations they provide to other modules by #defining the
 * modules you have below. Alternatively, you can -D them on
 * the compile line, as in cc -DeslAUGMENT_SSI -DeslAUGMENT_MSA.
 * 
 * If you compile and install the complete Easel library, all of these
 * get #defined automatically by ./configure, plus the eslLIBRARY flag
 * which means the full library with all augmentations is
 * available. So, if you steal files from an installed library, just
 * set these all back to #undef (depending on which files you have).
 *****************************************************************/

#define eslAUGMENT_ALPHABET 1
#define eslAUGMENT_DMATRIX 1
#define eslAUGMENT_FILEPARSER 1
#define eslAUGMENT_GEV 1
#define eslAUGMENT_GUMBEL 1
#define eslAUGMENT_HISTOGRAM 1
#define eslAUGMENT_KEYHASH 1
#define eslAUGMENT_MINIMIZER 1
#define eslAUGMENT_MSA 1
#define eslAUGMENT_RANDOM 1
#define eslAUGMENT_SSI 1
#define eslAUGMENT_STATS 1
#define eslLIBRARY 1

#endif /*ESL_CONFIG_INCLUDED*/
/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/
