/* esl_msa.h
 * Multiple sequence alignment file i/o.
 * 
 * SVN $Id: esl_msa.h 131 2006-11-14 19:07:57Z eddys $
 * SRE, Wed Jan 19 19:16:28 2005
 */
#ifndef eslMSA_INCLUDED
#define eslMSA_INCLUDED
#include <stdio.h>

#ifdef eslAUGMENT_KEYHASH
#include <esl_keyhash.h>
#endif
#ifdef eslAUGMENT_ALPHABET
#include <esl_alphabet.h>
#endif
#ifdef eslAUGMENT_SSI
#include <esl_ssi.h>
#endif

/* The following constants define the Pfam/Rfam cutoff set we propagate
 * from Stockholm format msa's into HMMER and Infernal models.
 */
/*::cexcerpt::msa_cutoffs::begin::*/
#define eslMSA_TC1     0
#define eslMSA_TC2     1
#define eslMSA_GA1     2
#define eslMSA_GA2     3
#define eslMSA_NC1     4
#define eslMSA_NC2     5
#define eslMSA_NCUTS   6
/*::cexcerpt::msa_cutoffs::end::*/

/* Object: ESL_MSA
 * 
 * A multiple sequence alignment.
 */
typedef struct {
  /* Mandatory information associated with the alignment.
   * (The important stuff.)
   */
  /*::cexcerpt::msa_mandatory::begin::*/
  char  **aseq;                /* alignment itself, [0..nseq-1][0..alen-1] */
  char  **sqname;              /* sequence names, [0..nseq-1][]            */
  double *wgt;                 /* sequence weights [0..nseq-1]             */
  int     alen;                /* length of alignment (columns)            */
  int     nseq;                /* number of seqs in alignment              */
  int     flags;               /* flags for what info has been set         */
  /*::cexcerpt::msa_mandatory::end::*/

#ifdef eslAUGMENT_ALPHABET
  /* When augmented w/ digital alphabets, we can store pre-digitized data in
   * ax[][], instead of the text info in aseq[][].
   */
  ESL_ALPHABET  *abc;		/* reference ptr to alphabet            */
  ESL_DSQ      **ax;		/* digitized aseqs [0..nseq-1][1..alen] */
#endif

  /* Optional information that we understand, and that we might have.
   * (The occasionally useful stuff.)
   */
  /*::cexcerpt::msa_optional::begin::*/
  char  *name;                  /* name of alignment, or NULL               */
  char  *desc;                  /* description of alignment, or NULL        */
  char  *acc;                   /* accession of alignment, or NULL          */
  char  *au;                    /* "author" information, or NULL            */
  char  *ss_cons;               /* consensus secondary structure, or NULL   */
  char  *sa_cons;               /* consensus surface accessibility, or NULL */
  char  *rf;                    /* reference coordinate system, or NULL     */
  char **sqacc;                 /* accession numbers for sequences i        */
  char **sqdesc;                /* description lines for sequences i        */
  char **ss;                    /* per-seq secondary structures, or NULL    */
  char **sa;                    /* per-seq surface accessibilities, or NULL */
  float  cutoff[eslMSA_NCUTS];  /* NC/TC/GA cutoffs propagated to Pfam/Rfam */
  int    cutset[eslMSA_NCUTS];  /* TRUE if a cutoff is set; else FALSE      */
  /*::cexcerpt::msa_optional::end::*/

  /* Info needed for maintenance of the data structure 
   * (internal stuff.)
   */
  int    sqalloc;		/* # seqs currently allocated for   */
  int   *sqlen;                 /* individual seq lengths during parsing    */
  int   *sslen;                 /* individual ss lengths during parsing     */
  int   *salen;                 /* individual sa lengths during parsing     */
  int    lastidx;		/* last index we saw; use for guessing next */

  /* Optional information, especially Stockholm markup.
   * (The stuff we don't understand, but we can regurgitate.)
   *
   * That is, we know what type of information it is, but it's
   * either (interpreted as) free-text comment, or it's Stockholm 
   * markup with unfamiliar tags.
   */
  char  **comment;              /* free text comments, or NULL      */
  int     ncomment;		/* number of comment lines          */
  int     alloc_ncomment;	/* number of comment lines alloc'ed */

  char  **gf_tag;               /* markup tags for unparsed #=GF lines  */
  char  **gf;                   /* annotations for unparsed #=GF lines  */
  int     ngf;			/* number of unparsed #=GF lines        */
  int     alloc_ngf;		/* number of gf lines alloc'ed          */

  char  **gs_tag;               /* markup tags for unparsed #=GS lines     */
  char ***gs;                   /* [0..ngs-1][0..nseq-1][free text] markup */
  int     ngs;                  /* number of #=GS tag types                */
  
  char  **gc_tag;               /* markup tags for unparsed #=GC lines  */
  char  **gc;                   /* [0..ngc-1][0..alen-1] markup         */
  int     ngc;                  /* number of #=GC tag types             */

  char  **gr_tag;               /* markup tags for unparsed #=GR lines   */
  char ***gr;                   /* [0..ngr][0..nseq-1][0..alen-1] markup */
  int     ngr;			/* number of #=GR tag types              */

  /* Optional augmentation w/ keyhashes. 
   * This can significantly speed up parsing of large alignments
   * with many (>1,000) sequences.
   */
#ifdef eslAUGMENT_KEYHASH 
  ESL_KEYHASH  *index;	        /* name ->seqidx hash table */
  ESL_KEYHASH  *gs_idx;         /* hash of #=GS tag types   */
  ESL_KEYHASH  *gc_idx;         /* hash of #=GC tag types   */
  ESL_KEYHASH  *gr_idx;         /* hash of #=GR tag types   */
#endif /*eslAUGMENT_KEYHASH*/
} ESL_MSA;

/* Flags for msa->flags
 */
#define eslMSA_HASWGTS (1 << 0)  /* 1 if wgts were set, 0 if default 1.0's */
#define eslMSA_DIGITAL (1 << 1)  /* if ax[][] is used instead of aseq[][]  */

                                     
/* Object: ESL_MSAFILE
 * 
 * Defines an alignment file that we open for reading.
 */
typedef struct {
  FILE *f;                      /* open file pointer                         */
  char *fname;			/* name of file. used for diagnostic output  */
  int   linenumber;		/* what line are we on in the file           */
  char  errbuf[eslERRBUFSIZE];  /* buffer for holding parse error info       */

  char *buf;			/* buffer for line input w/ sre_fgets()      */
  int   buflen;			/* current allocated length for buf          */

  int   do_gzip;		/* TRUE if f is "gzip -dc |" (will pclose(f))*/
  int   do_stdin;		/* TRUE if f is stdin (won't close f)        */
  int   format;			/* format of alignment file we're reading    */

  int   do_digital;		/* TRUE to digitize seqs directly into ax    */
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc;		/* AUGMENTATION (alphabet): digitized input  */
#endif
#ifdef eslAUGMENT_SSI		/* AUGMENTATION: SSI indexing of an MSA db   */
  ESL_SSI *ssi;		        /* open SSI index file; or NULL, if none.    */
#endif
} ESL_MSAFILE;


/* Alignment file format codes.
 * Must coexist with sqio.c/squid.h unaligned file format codes.
 * Rules:
 *     - 0 is an unknown/unassigned format 
 *     - <=100 reserved for unaligned formats
 *     - >100 reserved for aligned formats
 */
#define eslMSAFILE_UNKNOWN   0	  /* unknown format                              */
#define eslMSAFILE_STOCKHOLM 101  /* Stockholm format, interleaved               */
#define eslMSAFILE_PFAM      102  /* Pfam/Rfam one-line-per-seq Stockholm format */


/* Declarations of the API
 */
/* 1. The ESL_MSA object */
extern ESL_MSA *esl_msa_Create(int nseq, int alen);
extern ESL_MSA *esl_msa_CreateFromString(char *s, int fmt);
extern void     esl_msa_Destroy(ESL_MSA *msa);
extern int      esl_msa_Expand(ESL_MSA *msa);

/* 2. The ESL_MSAFILE object */
extern int  esl_msafile_Open(char *filename, int format, char *env, 
			     ESL_MSAFILE **ret_msafp);
extern void esl_msafile_Close(ESL_MSAFILE *afp);

/* 3. Digitized MSA's (ALPHABET augmentation required) */
#ifdef eslAUGMENT_ALPHABET
extern ESL_MSA *esl_msa_CreateDigital(ESL_ALPHABET *abc, int nseq, int alen);
extern int      esl_msa_Digitize(ESL_ALPHABET *abc, ESL_MSA *msa);
extern int      esl_msa_Textize(ESL_MSA *msa);
extern int      esl_msafile_OpenDigital(ESL_ALPHABET *abc, char *filename, 
					int format, char *env, ESL_MSAFILE **ret_msafp);
#endif

/* 4. General i/o API, all alignment formats */
extern int esl_msa_Read(ESL_MSAFILE *afp, ESL_MSA **ret_msa);
extern int esl_msa_Write(FILE *fp, ESL_MSA *msa, int fmt);
extern int esl_msa_GuessFileFormat(ESL_MSAFILE *afp);

/* 5. Miscellaneous functions for manipulating MSAs */
extern int esl_msa_SequenceSubset(ESL_MSA *msa, int *useme, ESL_MSA **ret_new);
extern int esl_msa_MinimGaps(ESL_MSA *msa, char *gaps);
extern int esl_msa_NoGaps(ESL_MSA *msa, char *gaps);
extern int esl_msa_SymConvert(ESL_MSA *msa, char *oldsyms, char *newsyms);

#endif /*eslMSA_INCLUDED*/

/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/
