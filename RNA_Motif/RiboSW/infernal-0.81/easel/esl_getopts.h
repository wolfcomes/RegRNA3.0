/* esl_getopts.h
 * 
 * Command line, config file, and environment variable
 * configuration of an application. Extends standard
 * UNIX/POSIX/GNU getopt().
 * 
 * SVN $Id: esl_getopts.h 144 2007-01-03 22:35:37Z eddys $
 * SRE, Thu Jan 13 08:38:28 2005 [St. Louis]
 */
#ifndef ESL_GETOPTS_INCLUDED
#define ESL_GETOPTS_INCLUDED

/* Object: ESL_OPTIONS
 * 
 * The application main.c defines an array of <ESL_OPTIONS> structures to
 * define what configuration options are used. The array is 
 * terminated by a structure containing { NULL, NULL, NULL, 0, NULL,
 * NULL, NULL, NULL} (or more simply, just 0 in all 8 fields.)
 */
/*::cexcerpt::options_object::begin::*/
typedef struct {
  char *name;           /* either short "-a" or long "--foo" style               */
  int   type;           /* arg type, for type checking: (eslARG_INT, etc.)       */
  char *defval;         /* default setting, or NULL ("default" is a C keyword)   */
  char *envvar;         /* associated environ var ("BLASTDB"), or NULL           */
  char *range;          /* for range checking arg: ("0<=x<=1", etc.)             */
  char *toggle_opts;    /* comma-sep'd optlist: turn these off if this opt is on */
  char *required_opts;  /* comma-sep'd optlist: these must also be set           */
  char *incompat_opts;  /* comma-sep'd optlist: these must not be set            */
  char *help;           /* help/usage string                                     */
  int   docgrouptag;    /* integer tag for documentation groups                  */
} ESL_OPTIONS;
/*::cexcerpt::options_object::end::*/

/* Argument types: the "type" variable in <ESL_OPTIONS>.
 */
#define eslARG_NONE      0	/* option takes no argument (so, is boolean)   */
#define eslARG_INT       1	/* arg convertable by atoi()               <n> */
#define eslARG_REAL      2	/* arg convertable by atof()               <x> */
#define eslARG_CHAR      3	/* arg is a single character               <c> */
#define eslARG_STRING    4	/* unchecked arg type                      <s> */
#define eslARG_INFILE    5      /* input file - same as string, shown as   <f> */
#define eslARG_OUTFILE   6      /* output file - same as string, shown as  <f> */



/* Object: ESL_GETOPTS
 * 
 * An <ESL_GETOPTS> object is created to parse configuration
 * from command line options, config file(s), and environment
 * variables.
 */
typedef struct {
  ESL_OPTIONS *opt;       /* array of app-defined options              */
  int          nopts;     /* number of options                         */

  int    argc;		  /* argc from command line                    */
  char **argv;		  /* argv from command line                    */
  char  *usage;		  /* command-line usage                        */
  int    optind;	  /* where we are in argc                      */
  int    argi;		  /* what command line arg we're on            */
  int    nfiles;	  /* # of cfgfiles that have been processed    */

  char **val;		  /* config'ed val for each option (as string) */
  int   *setby;		  /* array [0..nopts-1] for who set option i   */
  int   *valloc;          /* 0, or length of alloc for val[i]          */

  char  *optstring;	  /* internal: ptr into string of 1-char opts in argv[] */
} ESL_GETOPTS;


/* Possible values of the <setby> variable in ESL_GETOPTS.
 * Additionally, values of >3 also indicate a config file, in order 
 * of _ProcessConfigFile() calls (that is, setby=3 is the first 
 * config file, setby=4 is the second, etc.).
 */
#define eslARG_SETBY_DEFAULT  0
#define eslARG_SETBY_CMDLINE  1
#define eslARG_SETBY_ENV      2
#define eslARG_SETBY_CFGFILE  3


/* The visible API.
 */
extern ESL_GETOPTS *esl_getopts_Create(ESL_OPTIONS *opt, char *usage);
extern void         esl_getopts_Destroy(ESL_GETOPTS *g);
extern void         esl_getopts_Dump(FILE *ofp, ESL_GETOPTS *g);

extern int esl_opt_ProcessConfigfile(ESL_GETOPTS *g, char *filename, FILE *fp);
extern int esl_opt_ProcessEnvironment(ESL_GETOPTS *g);
extern int esl_opt_ProcessCmdline(ESL_GETOPTS *g, int argc, char **argv);
extern int esl_opt_VerifyConfig(ESL_GETOPTS *g);

extern int   esl_opt_IsSet(ESL_GETOPTS *g, char *optname);
extern int   esl_opt_GetBooleanOption(ESL_GETOPTS *g, char *optname, int *ret_state);
extern int   esl_opt_GetIntegerOption(ESL_GETOPTS *g, char *optname, int *ret_n);
extern int   esl_opt_GetFloatOption(ESL_GETOPTS *g, char *optname, float *ret_x);
extern int   esl_opt_GetDoubleOption(ESL_GETOPTS *g, char *optname, double *ret_x);
extern int   esl_opt_GetCharOption(ESL_GETOPTS *g, char *optname, char *ret_c);
extern int   esl_opt_GetStringOption(ESL_GETOPTS *g, char *optname, char **ret_s);
extern char *esl_opt_GetCmdlineArg(ESL_GETOPTS *g, int type, char *range);
#define esl_opt_ArgNumber(g)  ((g)->argc - (g)->optind)

extern int esl_opt_DisplayHelp(FILE *ofp, ESL_GETOPTS *go, int docgroup, int indent, int textwidth);

#endif /* ESL_GETOPTS_INCLUDED */
/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/
