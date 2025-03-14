/* src/config.h.  Generated from config.h.in by configure.  */
/* @configure_input@
 * DO NOT EDIT config.h!! 
 * config.h is generated from config.h.in by autoconf.
 * 
 * Configurable compile-time constants in INFERNAL.
 * 
 * Because this header may configure the behavior of system headers
 * (for example, LFS support), it must be included before any other
 * header file.
 * 
 * SRE, Sun Jun  3 20:22:38 2001 [St. Louis]
 * SVN $Id: config.h.in 1934 2007-04-12 15:00:23Z nawrockie $
 */
#ifndef CONFIGH_INCLUDED
#define CONFIGH_INCLUDED

/****************************************************************
 * This first section can be edited manually before compilation
 ****************************************************************/

/* RAMLIMIT (in MB) defines how much memory we're
 * allowed to expend on alignment algorithms without
 * switching to more efficient memory forms - e.g.
 * in smallcyk.c
 */
#ifndef RAMLIMIT
#define RAMLIMIT 0
#endif                                           

/* SRE_CONLEVEL will prob move to squid somewhere.
 *  Set to 1 to activate contract checking, during debugging.
 */                          
#define SRE_CONLEVEL 1
#if (SRE_CONLEVEL >= 1)
#include <assert.h>
#endif

/*****************************************************************
 * Everything else that follows is configured automatically 
 * by the ./configure script. DO NOT EDIT.
 *****************************************************************/

/* Version info - set once for whole package in configure.ac
 */
#define PACKAGE_NAME "Infernal"
#define PACKAGE_VERSION "0.81"
#define PACKAGE_DATE "May 2007"
#define PACKAGE_COPYRIGHT "Copyright (C) 2001-2007 HHMI Janelia Farm"
#define PACKAGE_LICENSE "Freely distributed under the GNU General Public License (GPL)"

/* Information about location of alloca()
 * Used by rigfilters/cm2hmm-1.0/MiscExceptions.cpp
 * This function is known to have portability issues (including
 * variable locations in system headers, and broken implementations
 * on certain platforms) and may be problematic.  Hopefully autoconf
 * will prevent major issues.
 */
#define HAVE_ALLOCA_H 1

/* --enable-debugging=x  debugging diagnostics (development versions only)
 */
#ifndef DEBUGLEVEL
/* #undef DEBUGLEVEL */
#endif

/* --enable-lfs          Large File Summit (LFS) support for >2Gb files
 */
/* #undef _LARGEFILE_SOURCE */
/* #undef _LARGEFILE64_SOURCE */
/* #undef _FILE_OFFSET_BITS */

/* --enable-mpi            MPI parallelization
 */
/* #undef USE_MPI */

#endif /* CONFIGH_INCLUDED */

