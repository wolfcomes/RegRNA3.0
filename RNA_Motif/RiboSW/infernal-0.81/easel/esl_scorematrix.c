/* Routines for manipulating sequence alignment score matrices,
 * such as the BLOSUM and PAM matrices.
 * 
 * Contents:
 *   1. The ESL_SCOREMATRIX object.
 *   2. Reading/writing score matrices.
 *   3. Interpreting score matrices probabilistically.
 *   4. Utility programs.
 *   5. Unit tests.
 *   6. Test driver.
 *   7. Example program.
 *   8. License and copyright.
 * 
 * SRE, Mon Apr  2 08:25:05 2007 [Janelia]
 * SVN $Id$
 */

#include <esl_config.h>

#include <string.h>
#include <math.h>

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_dmatrix.h>
#include <esl_fileparser.h>
#include <esl_rootfinder.h>
#include <esl_ratematrix.h>
#include <esl_scorematrix.h>

/*****************************************************************
 * 1. The ESL_SCOREMATRIX object
 *****************************************************************/

/* Function:  esl_scorematrix_Create()
 * Incept:    SRE, Mon Apr  2 08:38:10 2007 [Janelia]
 *
 * Purpose:   Allocates a score matrix for alphabet <abc>, initializes
 *            all scores to zero.
 *
 * Args:      abc   - pointer to digital alphabet 
 *
 * Returns:   a pointer to the new object.
 *
 * Throws:    <NULL> on allocation failure.
 */
ESL_SCOREMATRIX *
esl_scorematrix_Create(const ESL_ALPHABET *abc)
{
  int status;
  int i;
  ESL_SCOREMATRIX *S = NULL;

  ESL_ALLOC(S, sizeof(ESL_SCOREMATRIX));
  S->s          = NULL;
  S->K          = abc->K;
  S->Kp         = abc->Kp;
  S->isval      = NULL;
  S->abc_r      = abc;
  S->nc         = 0;
  S->outorder   = NULL;
  S->has_stop   = FALSE;
  S->stopsc     = 0;
  S->stopstopsc = 0;

  ESL_ALLOC(S->s, sizeof(int *) * abc->Kp);
  for (i = 0; i < abc->Kp; i++) S->s[i] = NULL;
  ESL_ALLOC(S->isval, sizeof(char) * abc->Kp);
  for (i = 0; i < abc->Kp; i++) S->isval[i] = FALSE;
  ESL_ALLOC(S->outorder, sizeof(char) * (abc->Kp+1)); /* maximum col/row count in output = Kp + stop character * */
  S->outorder[0] = '\0';		/* init to empty string. */

  ESL_ALLOC(S->s[0], sizeof(int) * abc->Kp * abc->Kp);
  for (i = 1; i < abc->Kp; i++) S->s[i] = S->s[0] + abc->Kp * i;

  for (i = 0; i < abc->Kp*abc->Kp; i++) S->s[0][i] = 0;
  return S;

 ERROR:
  esl_scorematrix_Destroy(S);
  return NULL;
}

/* Function:  esl_scorematrix_SetBLOSUM62
 * Incept:    SRE, Tue Apr  3 13:22:03 2007 [Janelia]
 *
 * Purpose:   Set the 20x20 canonical residue scores in an 
 *            allocated amino acid score matrix <S> to BLOSUM62
 *            scores \citep{Henikoff92}.
 *
 * Returns:   <eslOK> on success, and the scores in <S> are set.
 */
int
esl_scorematrix_SetBLOSUM62(ESL_SCOREMATRIX *S)
{
  int x,y;
  static int blosum62[28][28] = {
    /*  A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y    -    B    J    Z    O    U    X    ~  */
    {   4,   0,  -2,  -1,  -2,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,  -1,  -1,   1,   0,   0,  -3,  -2,   0,  -2,   0,  -1,   0,   0,   0,   0,  },
    {   0,   9,  -3,  -4,  -2,  -3,  -3,  -1,  -3,  -1,  -1,  -3,  -3,  -3,  -3,  -1,  -1,  -1,  -2,  -2,   0,  -3,   0,  -3,   0,   0,  -2,   0,  },
    {  -2,  -3,   6,   2,  -3,  -1,  -1,  -3,  -1,  -4,  -3,   1,  -1,   0,  -2,   0,  -1,  -3,  -4,  -3,   0,   4,   0,   1,   0,   0,  -1,   0,  },
    {  -1,  -4,   2,   5,  -3,  -2,   0,  -3,   1,  -3,  -2,   0,  -1,   2,   0,   0,  -1,  -2,  -3,  -2,   0,   1,   0,   4,   0,   0,  -1,   0,  },
    {  -2,  -2,  -3,  -3,   6,  -3,  -1,   0,  -3,   0,   0,  -3,  -4,  -3,  -3,  -2,  -2,  -1,   1,   3,   0,  -3,   0,  -3,   0,   0,  -1,   0,  },
    {   0,  -3,  -1,  -2,  -3,   6,  -2,  -4,  -2,  -4,  -3,   0,  -2,  -2,  -2,   0,  -2,  -3,  -2,  -3,   0,  -1,   0,  -2,   0,   0,  -1,   0,  },
    {  -2,  -3,  -1,   0,  -1,  -2,   8,  -3,  -1,  -3,  -2,   1,  -2,   0,   0,  -1,  -2,  -3,  -2,   2,   0,   0,   0,   0,   0,   0,  -1,   0,  },
    {  -1,  -1,  -3,  -3,   0,  -4,  -3,   4,  -3,   2,   1,  -3,  -3,  -3,  -3,  -2,  -1,   3,  -3,  -1,   0,  -3,   0,  -3,   0,   0,  -1,   0,  },
    {  -1,  -3,  -1,   1,  -3,  -2,  -1,  -3,   5,  -2,  -1,   0,  -1,   1,   2,   0,  -1,  -2,  -3,  -2,   0,   0,   0,   1,   0,   0,  -1,   0,  },
    {  -1,  -1,  -4,  -3,   0,  -4,  -3,   2,  -2,   4,   2,  -3,  -3,  -2,  -2,  -2,  -1,   1,  -2,  -1,   0,  -4,   0,  -3,   0,   0,  -1,   0,  },
    {  -1,  -1,  -3,  -2,   0,  -3,  -2,   1,  -1,   2,   5,  -2,  -2,   0,  -1,  -1,  -1,   1,  -1,  -1,   0,  -3,   0,  -1,   0,   0,  -1,   0,  },
    {  -2,  -3,   1,   0,  -3,   0,   1,  -3,   0,  -3,  -2,   6,  -2,   0,   0,   1,   0,  -3,  -4,  -2,   0,   3,   0,   0,   0,   0,  -1,   0,  },
    {  -1,  -3,  -1,  -1,  -4,  -2,  -2,  -3,  -1,  -3,  -2,  -2,   7,  -1,  -2,  -1,  -1,  -2,  -4,  -3,   0,  -2,   0,  -1,   0,   0,  -2,   0,  },
    {  -1,  -3,   0,   2,  -3,  -2,   0,  -3,   1,  -2,   0,   0,  -1,   5,   1,   0,  -1,  -2,  -2,  -1,   0,   0,   0,   3,   0,   0,  -1,   0,  },
    {  -1,  -3,  -2,   0,  -3,  -2,   0,  -3,   2,  -2,  -1,   0,  -2,   1,   5,  -1,  -1,  -3,  -3,  -2,   0,  -1,   0,   0,   0,   0,  -1,   0,  },
    {   1,  -1,   0,   0,  -2,   0,  -1,  -2,   0,  -2,  -1,   1,  -1,   0,  -1,   4,   1,  -2,  -3,  -2,   0,   0,   0,   0,   0,   0,   0,   0,  },
    {   0,  -1,  -1,  -1,  -2,  -2,  -2,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,   1,   5,   0,  -2,  -2,   0,  -1,   0,  -1,   0,   0,   0,   0,  },
    {   0,  -1,  -3,  -2,  -1,  -3,  -3,   3,  -2,   1,   1,  -3,  -2,  -2,  -3,  -2,   0,   4,  -3,  -1,   0,  -3,   0,  -2,   0,   0,  -1,   0,  },
    {  -3,  -2,  -4,  -3,   1,  -2,  -2,  -3,  -3,  -2,  -1,  -4,  -4,  -2,  -3,  -3,  -2,  -3,  11,   2,   0,  -4,   0,  -3,   0,   0,  -2,   0,  },
    {  -2,  -2,  -3,  -2,   3,  -3,   2,  -1,  -2,  -1,  -1,  -2,  -3,  -1,  -2,  -2,  -2,  -1,   2,   7,   0,  -3,   0,  -2,   0,   0,  -1,   0,  },
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  },
    {  -2,  -3,   4,   1,  -3,  -1,   0,  -3,   0,  -4,  -3,   3,  -2,   0,  -1,   0,  -1,  -3,  -4,  -3,   0,   4,   0,   1,   0,   0,  -1,   0,  },
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  },
    {  -1,  -3,   1,   4,  -3,  -2,   0,  -3,   1,  -3,  -1,   0,  -1,   3,   0,   0,  -1,  -2,  -3,  -2,   0,   1,   0,   4,   0,   0,  -1,   0,  },
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  },
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  },
    {   0,  -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -2,  -1,  -1,   0,   0,  -1,  -2,  -1,   0,  -1,   0,  -1,   0,   0,  -1,   0,  },
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  },
  };
  /* The BLOSUM62 background frequencies are the actual frequencies used to create
   * the matrix in 1992. */
  /*                           A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y */
  /* double blosum62f[20] = { 0.074, 0.025, 0.054, 0.054, 0.047, 0.074, 0.026, 0.068, 0.058, 0.099, 0.025, 0.045, 0.039, 0.034, 0.052, 0.057, 0.051, 0.073, 0.013, 0.032 };
   */

  for (x = 0;           x < S->K;  x++)      S->isval[x] = TRUE;
  for (x = S->abc_r->K; x < S->Kp; x++)      S->isval[x] = FALSE;
  x = esl_abc_DigitizeSymbol(S->abc_r, 'B'); S->isval[x] = TRUE;
  x = esl_abc_DigitizeSymbol(S->abc_r, 'Z'); S->isval[x] = TRUE;
  x = esl_abc_DigitizeSymbol(S->abc_r, 'X'); S->isval[x] = TRUE;
    
  for (x = 0; x < S->Kp; x++)
    for (y = 0; y < S->Kp; y++)
      S->s[x][y] = blosum62[x][y];

  /* Bookkeeping necessary to be able to reproduce BLOSUM62 output format exactly, if we need to Write() */
  strcpy(S->outorder, "ARNDCQEGHILKMFPSTWYVBZX*");
  S->nc         = strlen(S->outorder);
  S->has_stop   = TRUE;
  S->stopsc     = -4;
  S->stopstopsc = 1;

  return eslOK;
}


/* Function:  esl_scorematrix_SetWAG()
 * Synopsis:  Parameterize a score matrix from the WAG evolutionary model.           
 * Incept:    SRE, Thu Apr 12 13:23:28 2007 [Janelia]
 *
 * Purpose:   Parameterize an amino acid score matrix <S> using the WAG
 *            rate matrix \citep{WhelanGoldman01} as the underlying
 *            evolutionary model, at a distance of <t>
 *            substitutions/site, with scale factor <lambda>.
 *
 * Args:      S      - score matrix to set parameters in. Must be created for
 *                     an amino acid alphabet.
 *            lambda - scale factor for scores     
 *            t      - distance to exponentiate WAG to, in substitutions/site         
 *                 
 * Returns:   <eslOK> on success, and the 20x20 residue scores in <S> are set.
 *
 * Throws:    <eslEINVAL> if <S> isn't an allocated amino acid score matrix.
 *            <eslEMEM> on allocation failure.
 */
int
esl_scorematrix_SetWAG(ESL_SCOREMATRIX *S, const double lambda, const double t)
{
  int status;
  int i,j;
  ESL_DMATRIX *Q = NULL;
  ESL_DMATRIX *P = NULL;
  static double wagpi[20];

  if (S->K != 20) ESL_EXCEPTION(eslEINVAL, "Must be using an amino acid alphabet (K=20) to make WAG-based matrices");

  if (( Q = esl_dmatrix_Create(20, 20)) == NULL)  goto ERROR;
  if (( P = esl_dmatrix_Create(20, 20)) == NULL)  goto ERROR;
  if ( esl_composition_WAG(wagpi)       != eslOK) goto ERROR;
  if ( esl_rmx_SetWAG(Q, wagpi)         != eslOK) goto ERROR;
  if ( esl_dmx_Exp(Q, t, P)             != eslOK) goto ERROR;

  for (i = 0; i < 20; i++) 
    for (j = 0; j < 20; j++)
      P->mx[i][j] *= wagpi[i];	/* P_ij = P(j|i) pi_i */
  
  esl_scorematrix_SetFromProbs(S, lambda, P, wagpi, wagpi);

  esl_dmatrix_Destroy(Q);
  esl_dmatrix_Destroy(P);
  return eslOK;

 ERROR:
  if (Q != NULL) esl_dmatrix_Destroy(Q);
  if (Q != NULL) esl_dmatrix_Destroy(P);
  return status;
}


/* Function:  esl_scorematrix_SetFromProbs()
 * Synopsis:  Set new score matrix scores from target and background probabilities.
 * Incept:    SRE, Wed Apr 11 17:37:45 2007 [Janelia]
 *
 * Purpose:   Sets the scores in a new score matrix <S> from target joint
 *            probabilities in <P>, query background probabilities <fi>, and 
 *            target background probabilities <fj>, with scale factor <lambda>:
 *                 $s_{ij} = \frac{1}{\lambda} \frac{p_{ij}}{f_i f_j}$.
 *                 
 *            Size of everything must match the canonical alphabet
 *            size in <S>. That is, <S->abc->K> is the canonical
 *            alphabet size of <S>; <P> must contain $K times K$
 *            probabilities $P_{ij}$, and <fi>,<fj> must be vectors of
 *            K probabilities. All probabilities must be nonzero.
 *            
 * Args:      S      - score matrix to set scores in
 *            lambda - scale factor     
 *            P      - matrix of joint probabilities P_ij (KxK)
 *            fi     - query background probabilities (0..K-1)
 *            fj     - target background probabilities 
 *
 * Returns:   <eslOK> on success, and <S> contains the calculated score matrix.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_scorematrix_SetFromProbs(ESL_SCOREMATRIX *S, const double lambda, const ESL_DMATRIX *P, const double *fi, const double *fj)
{
  int    i,j;
  double sc;
  
  for (i = 0; i < S->abc_r->K; i++)
    for (j = 0; j < S->abc_r->K; j++)
      {
	sc = log(P->mx[i][j] / (fi[i] * fj[j])) / lambda;
	S->s[i][j] = (int) (sc + (sc>0 ? 0.5 : -0.5)); /* that's rounding to the nearest integer */
      }

  for (i = 0; i < S->abc_r->K; i++)
    S->isval[i] = TRUE;
  S->nc = S->abc_r->K;

  strncpy(S->outorder, S->abc_r->sym, S->abc_r->K);
  S->outorder[S->nc] = '\0';
  return eslOK;
}





/* Function:  esl_scorematrix_Compare()
 * Incept:    SRE, Tue Apr  3 14:17:12 2007 [Janelia]
 *
 * Purpose:   Compares two score matrices; returns <eslOK> if they 
 *            are identical, <eslFAIL> if they differ.
 */
int
esl_scorematrix_Compare(const ESL_SCOREMATRIX *S1, const ESL_SCOREMATRIX *S2)
{
  int a,b;

  if (strcmp(S1->outorder, S2->outorder) != 0) return eslFAIL;
  if (S1->nc         != S2->nc)                return eslFAIL;
  if (S1->has_stop   != S2->has_stop)          return eslFAIL;
  if (S1->stopsc     != S2->stopsc)            return eslFAIL;
  if (S1->stopstopsc != S2->stopstopsc)        return eslFAIL;
  
  for (a = 0; a < S1->nc; a++)
    if (S1->isval[a] != S2->isval[a])          return eslFAIL;
  
  for (a = 0; a < S1->Kp; a++)
    for (b = 0; b < S1->Kp; b++)
      if (S1->s[a][b] != S2->s[a][b])          return eslFAIL;
  return eslOK;
}


/* Function:  esl_scorematrix_Max()
 * Synopsis:  Returns maximum value in score matrix.
 * Incept:    SRE, Thu Apr 12 18:04:35 2007 [Janelia]
 *
 * Purpose:   Returns the maximum value in score matrix <S>.
 */
int
esl_scorematrix_Max(const ESL_SCOREMATRIX *S)
{
  int i,j;
  int max = S->s[0][0];

  for (i = 0; i < S->K; i++)
    for (j = 0; j < S->K; j++)
      if (S->s[i][j] > max) max = S->s[i][j];
  return max;
}

/* Function:  esl_scorematrix_Min()
 * Synopsis:  Returns minimum value in score matrix.
 * Incept:    SRE, Thu Apr 12 18:06:50 2007 [Janelia]
 *
 * Purpose:   Returns the minimum value in score matrix <S>.
 */
int
esl_scorematrix_Min(const ESL_SCOREMATRIX *S)
{
  int i,j;
  int min = S->s[0][0];

  for (i = 0; i < S->K; i++)
    for (j = 0; j < S->K; j++)
      if (S->s[i][j] < min) min = S->s[i][j];
  return min;
}



/* Function:  esl_scorematrix_Destroy()
 * Incept:    SRE, Mon Apr  2 08:46:44 2007 [Janelia]
 *
 * Purpose:   Frees a score matrix.
 *
 * Returns:   (void).
 */
void
esl_scorematrix_Destroy(ESL_SCOREMATRIX *S)
{
  if (S == NULL) return;
  if (S->s != NULL) {
    if (S->s[0] != NULL) free(S->s[0]);
    free(S->s);
  }
  if (S->isval    != NULL) free(S->isval);
  if (S->outorder != NULL) free(S->outorder);
  free(S);
  return;
}


/*****************************************************************
 * 2. Reading/writing score matrices.
 *****************************************************************/

/* Function:  esl_scorematrix_Read()
 * Incept:    SRE, Mon Apr  2 08:26:40 2007 [Janelia]
 *
 * Purpose:   Given a pointer <efp> to an open file parser for a file
 *            containing a score matrix (such as a PAM or BLOSUM
 *            matrix), parse the file and create a new score matrix
 *            object. The scores are expected to be for the alphabet
 *            <abc>. 
 *            
 *            The score matrix file is in the format that BLAST or
 *            FASTA use. The first line is a header contains N
 *            single-letter codes for the residues. Each of N
 *            subsequent rows optionally contains a residue row label
 *            (in the same order as the columns), followed by N
 *            residue scores.  (Older matrix files do not contain the
 *            leading row label; newer ones do.) The residues may
 *            appear in any order. They must minimally include the
 *            canonical K residues (K=4 for DNA, K=20 for protein),
 *            and may also contain none, some, or all degeneracy
 *            codes. Any other residue code that is not in the Easel
 *            digital alphabet (including, in particular, the '*' code
 *            for a stop codon) is ignored by the parser.
 *
 * Returns:   <eslOK> on success, and <ret_S> points to a newly allocated 
 *            score matrix. 
 *
 *            Returns <eslEFORMAT> on parsing error; in which case, <ret_S> is
 *            returned <NULL>.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_scorematrix_Read(ESL_FILEPARSER *efp, ESL_ALPHABET *abc, ESL_SCOREMATRIX **ret_S)
{
  int status;
  ESL_SCOREMATRIX *S     = NULL;
  int             *map   = NULL; /* maps col/row index to digital alphabet x */
  char            *tok;
  int              toklen;
  int              c, x;
  int              row,col;

  /* Allocate the matrix
   */
  if ((S = esl_scorematrix_Create(abc)) == NULL) { status = eslEMEM; goto ERROR; }

  /* Make sure we've got the comment character set properly in the fileparser.
   * Score matrices use #.
   */
  esl_fileparser_SetCommentChar(efp, '#');

  /* Look for the first non-blank, non-comment line in the file.  That line
   * gives us the single-letter codes in the order that the file's using.
   */
  if ((status = esl_fileparser_NextLine(efp)) != eslOK) ESL_XFAIL(eslEFORMAT, efp->errbuf, "file appears to be empty");

  /* Read the characters: count them and store them in order in label[0..nc-1].
   * nc cannot exceed Kp+1 in our expected alphabet (+1, for the stop character *)
   */
  S->nc = 0;
  while ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) == eslOK)
    {
      if (S->nc >= abc->Kp+1) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Header contains more residues than expected for alphabet");
      if (toklen != 1)        ESL_XFAIL(eslEFORMAT, efp->errbuf, "Header can only contain single-char labels; %s is invalid", tok);
      S->outorder[S->nc++] = *tok;
    }
  if (status != eslEOL) ESL_XFAIL(status, efp->errbuf, "Unexpected failure of esl_fileparser_GetTokenOnLine()");
  S->outorder[S->nc] = '\0';	/* NUL terminate */
  
  /* Verify that these labels for the score matrix seem plausible, given our alphabet.
   * This sets S->isval array: which residues we have scores for.
   * It also sets the map[] array, which maps coord in label[] to x in alphabet.
   * It's possible to see a residue in the score matrix that's not in the alphabet (main example is '*', a stop codon)
   */
  ESL_ALLOC(map, sizeof(int) * S->nc);
  for (c = 0; c < S->nc; c++)
    {
      if (esl_abc_CIsValid(abc, S->outorder[c])) 
	{  
	  x = esl_abc_DigitizeSymbol(abc, S->outorder[c]);
	  map[c] = x;
	  S->isval[x] = TRUE;
	}
      else if (S->outorder[c] == '*')
	{
	  S->has_stop = TRUE;
	  map[c] = -1;
	}
      else
	ESL_XFAIL(eslEFORMAT, efp->errbuf, "Don't know how to deal with residue %c in matrix file", S->outorder[c]);
    }
  for (x = 0; x < abc->K; x++)
    if (! S->isval[x]) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Expected to see a column for residue %c", abc->sym[x]);


  /* Read nc rows, one at a time;
   * on each row, read nc+1 or nc tokens, of which nc are scores (may lead with a label or not)
   */
  for (row = 0; row < S->nc; row++)
    {
      if ((status = esl_fileparser_NextLine(efp)) != eslOK) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Unexpectedly ran out of lines in file");
      for (col = 0; col < S->nc; col++)
	{
	  if ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Unexpectedly ran out of fields on line");
	  if (col == 0 && *tok == S->outorder[row]) { col--; continue; } /* skip leading label */

	  if (map[row] >= 0 && map[col] >= 0)
	    S->s[map[row]][map[col]] = atoi(tok);
	  else if (map[row] == -1 && map[col] == -1) /* stop/stop alignment */
	    S->stopstopsc = atoi(tok);
	  else 
	    S->stopsc = atoi(tok); /* this'll reset the stop score 2*nc-1 times, wastefully, and assuming they're all identical */
	}
      if ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslEOL)  ESL_XFAIL(eslEFORMAT, efp->errbuf, "Too many fields on line");
    }
  if ((status = esl_fileparser_NextLine(efp)) != eslEOF) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Too many lines in file");
  

  free(map);
  *ret_S = S;
  return eslOK;

 ERROR:
  esl_scorematrix_Destroy(S);
  if (map != NULL) free(map);
  *ret_S = NULL;
  return status;
}

/* Function:  esl_scorematrix_Write()
 * Incept:    SRE, Tue Apr  3 13:55:10 2007 [Janelia]
 *
 * Purpose:   Writes a score matrix <S> to an open stream <fp>, in 
 *            format compatible with BLAST, FASTA, and other common
 *            sequence alignment software.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_scorematrix_Write(FILE *fp, const ESL_SCOREMATRIX *S)
{
  int a,b;			
  int x,y;
  int nc = 0;
  
  /* Total paranoia: we have two redundant ways to determine the
   * number of residues in this matrix, and they should match:
   * S->nc, or the sum of the isval[] flags + has_stop.
   */
  if (S->has_stop) nc++;
  for (x = 0; x < S->Kp; x++)
    if (S->isval[x]) nc++;
  if (nc != S->nc) ESL_EXCEPTION(eslEINVAL, "nc's don't match. matrix is corrupt");

  /* The header line, with column labels for residues */
  fprintf(fp, "  ");
  for (a = 0; a < nc; a++) fprintf(fp, "  %c ", S->outorder[a]);
  fprintf(fp, "\n");
  
  /* The data. Watch out for those pesky *'s, which aren't in the Easel digital alphabet (yet)
   */
  for (a = 0; a < nc; a++)
    {
      fprintf(fp, "%c ", S->outorder[a]);
      for (b = 0; b < nc; b++)
	{
	  if (S->outorder[a] != '*' && S->outorder[b] != '*') 
	    {
	      x = esl_abc_DigitizeSymbol(S->abc_r, S->outorder[a]);
	      y = esl_abc_DigitizeSymbol(S->abc_r, S->outorder[b]);
	      fprintf(fp, "%3d ", S->s[x][y]);
	    } 
	  else if (S->outorder[a] != '*' || S->outorder[b] != '*')
	    fprintf(fp, "%3d ", S->stopsc);
	  else
	    fprintf(fp, "%3d ", S->stopstopsc);
	}
      fprintf(fp, "\n");
    }
  
  return eslOK;
}

/*****************************************************************
 * 3. Interpreting score matrices probabilistically.
 *****************************************************************/ 

/* Function:  esl_scorematrix_ObtainPij()
 * Synopsis:  Obtain $P_{ij}$ from score matrix with known $\lambda$ and background $f$'s.
 * Incept:    SRE, Thu Apr 12 17:46:20 2007 [Janelia]
 *
 * Purpose:   Given a score matrix <S> with known <lambda> and known
 *            query and target background frequencies <fi> and <fj>, calculate
 *            joint probabilities <P>.
 *            
 *            Caller provides square <P> matrix allocated for <S->K>
 *            $\times$ <S->K> joint residue probabilities $p_{ij}$. 
 *
 * Returns:   <eslOK> on success, and the joint probabilities are in <P>.
 */
int
esl_scorematrix_ObtainPij(const ESL_SCOREMATRIX *S, const double *fi, const double *fj, const double lambda, ESL_DMATRIX *P)
{
  int i,j;

  for (i = 0; i < S->K; i++)
    for (j = 0; j < S->K; j++)
      P->mx[i][j] = fi[i] * fj[j] * exp(lambda * (double) S->s[i][j]);
  return eslOK;
}



struct lambda_params {
  const double *fi;
  const double *fj;
  const ESL_SCOREMATRIX *S;
};

static int
lambda_fdf(double lambda, void *params, double *ret_fx, double *ret_dfx)
{
  struct lambda_params *p = (struct lambda_params *) params;
  int    i,j;
  double tmp;
  
  *ret_fx  = 0.;
  *ret_dfx = 0.;
  for (i = 0; i < p->S->K; i++)
    for (j = 0; j < p->S->K; j++)
      {
	tmp      = p->fi[i] * p->fj[j] * exp(lambda * (double) p->S->s[i][j]);
	*ret_fx  += tmp;
	*ret_dfx += tmp * (double) p->S->s[i][j];
      }
  *ret_fx -= 1.0;
  return eslOK;
}

/* Function:  esl_scorematrix_SolveLambda()
 * Synopsis:  Find $\lambda$ for score matrix, given background.
 * Incept:    SRE, Thu Apr 12 18:09:41 2007 [Janelia]
 *
 * Purpose:   Given valid score matrix <S> and query and target background
 *            probabilities <fi>, <fj>, solve for $\lambda$. 
 *
 * Args:      
 *
 * Returns:   <eslOK> on success, and $\lambda$ is returned in <ret_lambda>;
 *            additionally, if caller provides a non-<NULL> <P> allocated
 *            for joint probabilities $p_{ij}$, the joint probabilities are
 *            calculated and left there.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_scorematrix_SolveLambda(const ESL_SCOREMATRIX *S, const double *fi, const double *fj, ESL_DMATRIX *P, double *ret_lambda)
{
  int status;
  ESL_ROOTFINDER *R = NULL;
  struct lambda_params p;
  double lambda_guess;
  double fx, dfx;
  p.fi = fi;
  p.fj = fj;
  p.S  = S;

  /* It's important that we come at the root from the far side, where
   * f(lambda) is positive; else we may identify the root we don't want
   * at lambda=0.
   */
  lambda_guess = 1. / (double) esl_scorematrix_Max(S);
  for (; lambda_guess < 50.; lambda_guess *= 2.0) {
    lambda_fdf(lambda_guess, &p, &fx, &dfx);
    if (fx > 0) break;
  }
  if (fx <= 0) ESL_EXCEPTION(eslEINVAL, "Failed to bracket root for solving lambda");

  /* Find lambda by Newton/Raphson */
  if ((  R = esl_rootfinder_CreateFDF(lambda_fdf, &p) ) == NULL) { status = eslEMEM; goto ERROR; }
  if (( status = esl_root_NewtonRaphson(R, lambda_guess, ret_lambda))  != eslOK) goto ERROR;
  
  if (P != NULL) {
    if ((status = esl_scorematrix_ObtainPij(S, fi, fj, *ret_lambda, P)) != eslOK) goto ERROR;
  }
  esl_rootfinder_Destroy(R);
  return eslOK;

 ERROR:
  if (R != NULL) esl_rootfinder_Destroy(R);
  *ret_lambda = 0.;
  return status;
}  




/* This section is an implementation of one of the ideas in
 * Yu and Altschul, PNAS 100:15688, 2003 [YuAltschul03]:
 * given a valid score matrix, calculate its probabilistic
 * basis (P_ij, f_i, f_j, and lambda), on the assumption that
 * the background probabilities are the marginals of P_ij.
 */

struct yualtschul_params {
  ESL_DMATRIX *S;   /* pointer to the KxK score matrix w/ values cast to doubles */		
  ESL_DMATRIX *M;   /* not a param per se: alloc'ed storage for M matrix provided to the objective function */
  ESL_DMATRIX *Y;   /* likewise, alloc'ed storage for Y (M^-1) matrix provided to obj function */
};

/* yualtschul_func()
 *
 * This is the objective function we try to find a root of. 
 * Its prototype is dictated by the esl_rootfinder API.
 */
static int
yualtschul_func(double lambda, void *params, double *ret_fx)
{
  int status;
  struct yualtschul_params *p = (struct yualtschul_params *) params;
  ESL_DMATRIX  *S = p->S;
  ESL_DMATRIX  *M = p->M;
  ESL_DMATRIX  *Y = p->Y;
  int i,j;

  /* the M matrix has entries M_ij = e^{lambda * s_ij} */
  for (i = 0; i < S->n; i++)
    for (j = 0; j < S->n; j++)
      M->mx[i][j] = exp(lambda * S->mx[i][j]);

  /* the Y matrix is the inverse of M */
  if ((status = esl_dmx_Invert(M, Y)) != eslOK) return status;

  /* We're trying to find the root of \sum_ij Y_ij - 1 = 0 */
  *ret_fx = esl_dmx_Sum(Y) - 1.;
  return eslOK;
}

/* yualtschul_engine()
 *
 * This function backcalculates the probabilistic basis for a score
 * matrix S, when S is a double-precision matrix. Providing this
 * as a separate "engine" and writing esl_scorematrix_ReverseEngineer()
 * as a wrapper around it allows us to separately test inaccuracy
 * due to numerical performance of our linear algebra, versus 
 * inaccuracy due to integer roundoff in integer scoring matrices.
 * 
 * It is not uncommon for this to fail, when S is derived from
 * integer scores. Because the scores may have been provided by the
 * user, and this may be our first chance to detect the "user error"
 * of an invalid matrix, this engine returns <eslENORESULT> as a normal error
 * if it can't reach a valid solution.
 */
static int 
yualtschul_engine(ESL_DMATRIX *S, ESL_DMATRIX *P, double *fi, double *fj, double *ret_lambda)
{
  int status;
  ESL_ROOTFINDER *R = NULL;
  struct yualtschul_params p;
  double lambda;
  double xl, xr;
  double fx;
  int    i,j;

  /* Set up a bisection method to find lambda */
  p.S = S;
  p.M = p.Y = NULL;
  if ((p.M = esl_dmatrix_Create(S->n, S->n))           == NULL) { status = eslEMEM; goto ERROR; }
  if ((p.Y = esl_dmatrix_Create(S->n, S->n))           == NULL) { status = eslEMEM; goto ERROR; }
  if ((R = esl_rootfinder_Create(yualtschul_func, &p)) == NULL) { status = eslEMEM; goto ERROR; }

  /* Need a reasonable initial guess for lambda; if we use extreme
   * lambda guesses, we'll introduce numeric instability in the
   * objective function, and may even blow up the values of e^{\lambda
   * s_ij} in the M matrix. Appears to be safe to start with lambda on
   * the order of 2/max(s_ij).
   */
  xr = 1. / esl_dmx_Max(S);
  
  /* Identify suitable brackets on lambda. */
  for (xl = xr; xl > 1e-10; xl /= 1.6) {
    if ((status = yualtschul_func(xl, &p, &fx))  != eslOK) goto ERROR;
    if (fx > 0.) break;
  }
  if (fx <= 0.) { status = eslENORESULT; goto ERROR; }

  for (; xr < 100.; xr *= 1.6) {
    if ((status = yualtschul_func(xr, &p, &fx))  != eslOK) goto ERROR;
    if (fx < 0.) break;
  }
  if (fx >= 0.) { status = eslENORESULT; goto ERROR; }

  /* Find lambda by bisection */
  if (esl_root_Bisection(R, xl, xr, &lambda) != eslOK)     goto ERROR;

  /* Find fi, fj from Y: fi are column sums, fj are row sums */
  for (i = 0; i < S->n; i++) {
    fi[i] = 0.;
    for (j = 0; j < S->n; j++) fi[i] += p.Y->mx[j][i];
  }
  for (j = 0; j < S->n; j++) {
    fj[j] = 0.;
    for (i = 0; i < S->n; i++) fj[j] += p.Y->mx[j][i];
  }

  /* Find p_ij */
  for (i = 0; i < S->n; i++) 
    for (j = 0; j < S->n; j++)
      P->mx[i][j] = fi[i] * fj[j] * p.M->mx[i][j];

  *ret_lambda = lambda;
  esl_dmatrix_Destroy(p.M);
  esl_dmatrix_Destroy(p.Y);
  esl_rootfinder_Destroy(R);
  return eslOK;

 ERROR:
  if (p.M != NULL) esl_dmatrix_Destroy(p.M);
  if (p.Y != NULL) esl_dmatrix_Destroy(p.Y);
  if (R   != NULL) esl_rootfinder_Destroy(R);
  return status;
}

/* Function:  esl_scorematrix_ReverseEngineer()
 * Synopsis:  Calculate the probabilistic basis of a score matrix.
 * Incept:    SRE, Wed Apr 11 07:56:44 2007 [Janelia]
 *
 * Purpose:   Reverse engineering of a score matrix: given a "valid"
 *            substitution matrix <S>, obtain implied joint
 *            probabilities <P>, query composition <fi>, target
 *            composition <fj>, and scale <lambda>, by assuming that
 *            <fi> and <fj> are the appropriate marginals of <P>.
 *
 *            This implements the algorithm described by [YuAltschul03].
 *            
 *            Caller provides the allocated space for a $K \times K$
 *            matrix <P>, and $K$-vectors <fi> and <fj>, where
 *            $K$ is the base alphabet size (<S->K>). 
 *            
 *            This algorithm works fine in principle, but when it is
 *            applied to rounded integer scores with small dynamic
 *            range (the typical situation for score matrices) it may
 *            fail due to integer roundoff error. It works best for
 *            score matrices built using small values of $\lambda$. Yu
 *            and Altschul use $\lambda = 0.00635$ for BLOSUM62, which
 *            amounts to scaling default BLOSUM62 up 50-fold. It
 *            happens that default BLOSUM62 (which was created with
 *            lambda = 0.3466, half-bits) can be successfully reverse
 *            engineered (albeit with some loss of accuracy;
 *            calculated lambda is 0.3240) but other common matrices
 *            may fail.
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success, and <P>, <fi>, <fj>, and <ret_lambda>
 *            contain the results.
 *            
 *            <eslENORESULT> if the algorithm fails to determine a
 *            valid solution.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      J1/35.
 */
int
esl_scorematrix_ReverseEngineer(const ESL_SCOREMATRIX *S, ESL_DMATRIX *P, double *fi, double *fj, double *ret_lambda)
{
  int status;
  int i,j;
  ESL_DMATRIX    *Sd  = NULL;

  if (( Sd = esl_dmatrix_Create(S->K, S->K))  == NULL)  goto ERROR;

  /* Construct a double-precision dmatrix from S.
   * I've tried integrating over the rounding uncertainty by
   * averaging over trials with values jittered by +/- 0.5,
   * but it doesn't appear to help much, if at all.
   */
  for (i = 0; i < S->K; i++) 
    for (j = 0; j < S->K; j++)
      Sd->mx[i][j] = (double) S->s[i][j];

  /* Reverse engineer the doubles */
  if ((status = yualtschul_engine(Sd, P, fi, fj, ret_lambda)) != eslOK) goto ERROR;
      
  esl_dmatrix_Destroy(Sd);
  return eslOK;

 ERROR:
  if (Sd  != NULL) esl_dmatrix_Destroy(Sd);
  return status;
}








/*****************************************************************
 * 4. Utilities
 *****************************************************************/ 

/* Reformat a score matrix file, canonical residues only, into
 * Easel internal digital alphabet order, suitable for making 
 * a static data structure.
 */
#ifdef eslSCOREMATRIX_UTILITY1
/* 
    gcc -g -Wall -o utility -I. -L. -DeslSCOREMATRIX_UTILITY1 esl_scorematrix.c -leasel -lm
    ./utility BLOSUM62
*/
#include <easel.h>
#include <esl_alphabet.h>
#include <esl_scorematrix.h>
#include <esl_fileparser.h>

int
main(int argc, char **argv)
{
  char *infile = argv[1];
  ESL_ALPHABET    *abc;
  ESL_FILEPARSER  *efp;
  ESL_SCOREMATRIX *S;
  int x,y;

  abc = esl_alphabet_Create(eslAMINO);

  if (esl_fileparser_Open(infile, &efp)  != eslOK) esl_fatal("Failed to open %s\n", infile);
  if (esl_scorematrix_Read(efp, abc, &S) != eslOK) esl_fatal("parse failed: %s", efp->errbuf);

  for (x = 0; x < abc->Kp; x++) {
    printf("{ ");
    for (y = 0; y < abc->Kp; y++)
      printf("%3d, ", S->s[x][y]);
    printf(" },\n");
  }
  
  esl_scorematrix_Destroy(S);
  esl_fileparser_Close(efp);
  esl_alphabet_Destroy(abc);
  return eslOK;
}
#endif /*eslSCOREMATRIX_UTILITY1*/







/*****************************************************************
 * 5. Unit tests.
 *****************************************************************/

#ifdef eslSCOREMATRIX_TESTDRIVE
#include <esl_dirichlet.h>

static void
utest_ReadWrite(ESL_ALPHABET *abc, ESL_SCOREMATRIX *S)
{
  char tmpfile[16]     = "esltmpXXXXXX";
  FILE            *fp  = NULL;
  ESL_SCOREMATRIX *S2  = NULL;
  ESL_FILEPARSER  *efp = NULL;
  
  if (esl_tmpfile_named(tmpfile, &fp)       != eslOK) esl_fatal("failed to open tmp file");
  if (esl_scorematrix_Write(fp, S)          != eslOK) esl_fatal("failed to write test matrix");
  fclose(fp);

  if (esl_fileparser_Open(tmpfile, &efp)    != eslOK) esl_fatal("failed to open tmpfile containing BLOSUM62 matrix");
  if (esl_scorematrix_Read(efp, abc, &S2)   != eslOK) esl_fatal("failed to read tmpfile containing BLOSUM62 matrix");
  if (esl_scorematrix_Compare(S, S2)        != eslOK) esl_fatal("the two test matrices aren't identical");
  
  remove(tmpfile); 
  esl_fileparser_Close(efp);
  esl_scorematrix_Destroy(S2);
  return;
}


static void
utest_SolveLambda(ESL_ALPHABET *abc, ESL_SCOREMATRIX *S0, ESL_DMATRIX *P0, double *wagpi, double lambda0)
{
  char *msg = "SolveLambda() unit test failed";
  ESL_DMATRIX     *P    = NULL;
  double           lambda;

  if ((P   = esl_dmatrix_Create(S0->K, S0->K))== NULL)  esl_fatal(msg);

  if (esl_scorematrix_SolveLambda(S0, wagpi, wagpi, P, &lambda) != eslOK) esl_fatal(msg);

  if (esl_DCompare(lambda0, lambda, 1e-3)     != eslOK) esl_fatal("lambda is wrong");
  if (esl_DCompare(esl_dmx_Sum(P), 1.0, 1e-9) != eslOK) esl_fatal("P doesn't sum to 1");
  if (esl_dmatrix_Compare(P0, P, 1e-2)        != eslOK) esl_fatal("P is wrong");

  esl_dmatrix_Destroy(P);
  return;
}
 

/* The scores->pij reverse engineering engine works with scores in doubles,
 * so we can separate effects of rounding to integers in standard
 * score matrices.
 */
static void 
utest_yualtschul(ESL_DMATRIX *P0, double *wagpi)
{
  char *msg = "reverse engineering engine test failed";
  ESL_DMATRIX     *S   = NULL;	/* original score matrix, in double form, not rounded to ints (calculated from P, fi, fj) */
  ESL_DMATRIX     *P   = NULL;	/* backcalculated P_ij joint probabilities */
  double          *fi  = NULL;	/* backcalculated f_i query composition */
  double          *fj  = NULL;	/* backcalculated f'_j target composition */
  double           lambda0;	/* true lambda */
  double           lambda;	/* backcalculated lambda */
  int              i,j;

  /* Allocations */
  if (( S  = esl_dmatrix_Create(20, 20))     == NULL)  esl_fatal(msg);
  if (( P  = esl_dmatrix_Create(20, 20))     == NULL)  esl_fatal(msg);
  if ((fi  = malloc(sizeof(double) * 20))    == NULL)  esl_fatal(msg);
  if ((fj  = malloc(sizeof(double) * 20))    == NULL)  esl_fatal(msg);

  /* Make a WAG-based score matrix in double-precision, without rounding to integers */
  lambda0 = 0.3;
  for (i = 0; i < 20; i++) 
    for (j = 0; j < 20; j++)
      S->mx[i][j] = log(P0->mx[i][j] / (wagpi[i] * wagpi[j])) / lambda0;

  /* Reverse engineer it in double precision */
  if ( yualtschul_engine(S, P, fi, fj, &lambda) != eslOK) esl_fatal("reverse engineering engine failed");

  /* Validate the solution (expect more accuracy from this than from integer scores) */
  if (esl_DCompare(lambda0, lambda, 1e-4)      != eslOK) esl_fatal("failed to get right lambda");
  if (esl_DCompare(esl_dmx_Sum(P),  1.0, 1e-6) != eslOK) esl_fatal("reconstructed P doesn't sum to 1");  
  if (esl_dmatrix_Compare(P, P0, 1e-3)         != eslOK) esl_fatal("failed to recover correct P_ij");
  for (i = 0; i < 20; i++) 
    {
      if (esl_DCompare(fi[i],    fj[i],  1e-6) != eslOK) esl_fatal("background fi, fj not the same");
      if (esl_DCompare(wagpi[i], fi[i],  1e-3) != eslOK) esl_fatal("failed to reconstruct WAG backgrounds");  
    }

  free(fj);
  free(fi);
  esl_dmatrix_Destroy(S);
  esl_dmatrix_Destroy(P);
  return;
}


static void
utest_ReverseEngineer(ESL_SCOREMATRIX *S0, ESL_DMATRIX *P0, double *wagpi, double lambda0)
{
  ESL_DMATRIX     *P  = NULL;
  double          *fi = NULL;
  double          *fj = NULL;
  double           lambda;	/* reconstructed lambda */

  if ((P   = esl_dmatrix_Create(S0->K, S0->K)) == NULL)  esl_fatal("P allocation failed");
  if ((fi  = malloc(sizeof(double) * S0->K))   == NULL)  esl_fatal("fi allocation failed");
  if ((fj  = malloc(sizeof(double) * S0->K))   == NULL)  esl_fatal("fj allocation failed");

  if (esl_scorematrix_ReverseEngineer(S0, P, fi, fj, &lambda) != eslOK) esl_fatal("reverse engineering failed");

  /* Validate the solution, gingerly (we expect significant error due to integer roundoff) */
  if (esl_DCompare(lambda0, lambda, 0.01)       != eslOK) esl_fatal("failed to get right lambda");
  if (esl_DCompare(esl_dmx_Sum(P),  1.0, 0.001) != eslOK) esl_fatal("reconstructed P doesn't sum to 1");
  if (esl_dmatrix_Compare(P, P0, 0.1)           != eslOK) esl_fatal("reconstructed P is wrong");

  free(fj);
  free(fi);
  esl_dmatrix_Destroy(P);
  return;
}


#endif /*eslSCOREMATRIX_TESTDRIVE*/


/*****************************************************************
 * 6. Test driver.
 *****************************************************************/
/* 
    gcc -g -Wall -I. -L. -o test -DeslSCOREMATRIX_TESTDRIVE esl_scorematrix.c -leasel -lm
    ./test
*/
#ifdef eslSCOREMATRIX_TESTDRIVE
#include <easel.h>
#include <esl_scorematrix.h>

int 
main(int argc, char **argv)
{
  ESL_ALPHABET    *abc = NULL;	/* amino acid alphabet */
  ESL_SCOREMATRIX *BL62= NULL;	/* BLOSUM62 matrix */
  ESL_SCOREMATRIX *S0  = NULL;	/* original score matrix (calculated from P, fi, fj) */
  ESL_DMATRIX     *P0  = NULL;	/* original P_ij joint probabilities */
  ESL_DMATRIX     *Q   = NULL;	/* WAG rate matrix */
  double           lambda0;	/* true lambda used to construct S */
  double           t;
  int              i,j;
  static double    wagpi[20];

  /* Allocations */
  if ((abc = esl_alphabet_Create(eslAMINO))      == NULL)  esl_fatal("allocation of alphabet failed");
  if ((BL62= esl_scorematrix_Create(abc))        == NULL)  esl_fatal("allocation of BLOSUM62 failed");
  if ((S0  = esl_scorematrix_Create(abc))        == NULL)  esl_fatal("allocation of scorematrix failed");
  if ((P0  = esl_dmatrix_Create(abc->K, abc->K)) == NULL)  esl_fatal("P allocation failed");
  if ((Q   = esl_dmatrix_Create(abc->K, abc->K)) == NULL)  esl_fatal("Q allocation failed");

  /* Make a BLOSUM matrix */
  if ( esl_scorematrix_SetBLOSUM62(BL62) != eslOK) esl_fatal("failed to set a BLOSUM matrix");

  /* Make a WAG-based score matrix with small lambda. */
  lambda0 = 0.00635;
  t    = 2.0;
  esl_scorematrix_SetWAG(S0, lambda0, t);
  esl_composition_WAG(wagpi);

  /* Redo some calculations to get the known probabilistic basis of that S */
  if ( esl_rmx_SetWAG(Q, wagpi)  != eslOK) esl_fatal("failed to set WAG");
  if ( esl_dmx_Exp(Q, t, P0)     != eslOK) esl_fatal("failed to exponentiate WAG");
  for (i = 0; i < 20; i++) 
    for (j = 0; j < 20; j++)
      P0->mx[i][j] *= wagpi[i];	/* P_ij = P(j|i) pi_i */

  /* The unit test battery
   */
  utest_ReadWrite(abc, BL62);
  utest_ReadWrite(abc, S0);
  utest_SolveLambda(abc, S0, P0, wagpi, lambda0);
  utest_yualtschul(P0, wagpi);
  utest_ReverseEngineer(S0, P0, wagpi, lambda0); 

  esl_dmatrix_Destroy(Q);
  esl_dmatrix_Destroy(P0);
  esl_scorematrix_Destroy(BL62);
  esl_scorematrix_Destroy(S0);
  esl_alphabet_Destroy(abc);

  return 0;
}
#endif /*eslSCOREMATRIX_TESTDRIVE*/

/*****************************************************************
 * 7. Example program
 *****************************************************************/
/* 
    gcc -g -Wall -I. -L. -o example -DeslSCOREMATRIX_EXAMPLE esl_scorematrix.c -leasel -lm
    ./example <score matrix file>
*/
#ifdef eslSCOREMATRIX_EXAMPLE
#include <easel.h>
#include <esl_alphabet.h>
#include <esl_fileparser.h>
#include <esl_dmatrix.h>
#include <esl_scorematrix.h>

int main(int argc, char **argv)
{
  char            *scorefile = argv[1];
  ESL_ALPHABET    *abc       = esl_alphabet_Create(eslAMINO);
  ESL_FILEPARSER  *efp       = NULL;
  ESL_SCOREMATRIX *S         = NULL;
  ESL_DMATRIX     *P         = esl_dmatrix_Create(abc->K, abc->K);
  double          *fi        = malloc(sizeof(double) * abc->K);
  double          *fj        = malloc(sizeof(double) * abc->K);
  double           lambda;
  
  /* Input an amino acid score matrix from a file. */
  if ( esl_fileparser_Open(scorefile, &efp) != eslOK) esl_fatal("failed to open score file %s", scorefile);
  if ( esl_scorematrix_Read(efp, abc, &S)   != eslOK) esl_fatal("failed to read matrix from %s", scorefile);
  esl_fileparser_Close(efp);

  /* Reverse engineer it to get implicit probabilistic model. */
  if ( esl_scorematrix_ReverseEngineer(S, P, fi, fj, &lambda) != eslOK) esl_fatal("reverse engineering failed");

  /* Print lambda, and the joint probabilities. */
  printf("Lambda is %g\n\n", lambda);
  printf("Implicit joint probabilities are:\n");
  esl_dmatrix_Dump(stdout, P, abc->sym, abc->sym);

  free(fi);
  free(fj);
  esl_dmatrix_Destroy(P);
  esl_scorematrix_Destroy(S);
  esl_alphabet_Destroy(abc);
  return 0;
}
#endif /*eslSCOREMATRIX_EXAMPLE*/


/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/ 
