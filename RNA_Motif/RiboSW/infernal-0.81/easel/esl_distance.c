/* Pairwise identities, distances, and distance matrices.
 *
 * Contents:
 *    1. Pairwise distances for aligned text sequences.
 *    2. Pairwise distances for aligned digital seqs.     [alphabet]
 *    3. Distance matrices for aligned text sequences.    [dmatrix]
 *    4. Distance matrices for aligned digital sequences. [alphabet,dmatrix]
 *    5. Private (static) functions.
 *    6. Unit tests.
 *    7. Test driver.
 *    8. Example.
 *    9. Copyright notice and license.
 *    
 * SVN $Id: esl_distance.c 153 2007-02-18 03:46:56Z eddys $
 * SRE, Mon Apr 17 20:05:43 2006 [St. Louis]
 */
#include <esl_config.h>

#include <ctype.h>
#include <string.h>
#include <math.h>

#include <easel.h>
#ifdef eslAUGMENT_ALPHABET
#include <esl_alphabet.h>
#endif
#ifdef eslAUGMENT_DMATRIX
#include <esl_dmatrix.h>
#endif

/* Forward declaration of our static functions.
 */
static int jukescantor(int n1, int n2, int alphabet_size, double *ret_distance, double *ret_variance);





/*****************************************************************
 * 1. Pairwise distances for aligned text sequences.
 *****************************************************************/

/* Function:  esl_dst_CPairId()
 * Synopsis:  Pairwise identity of two aligned text strings.
 * Incept:    SRE, Mon Apr 17 20:06:07 2006 [St. Louis]
 *
 * Purpose:   Calculates pairwise fractional identity between two
 *            aligned character strings <asq1> and <asq2>. 
 *            Return this distance in <ret_pid>; return the
 *            number of identities counted in <ret_nid>; and
 *            return the denominator <MIN(len1,len2)> in
 *            <ret_n>.
 *            
 *            Alphabetic symbols <[a-zA-Z]> are compared
 *            case-insensitively for identity. Any nonalphabetic
 *            character is assumed to be a gap symbol.
 *            
 *            This simple comparison rule is unaware of synonyms and
 *            degeneracies in biological alphabets.  For a more
 *            sophisticated and biosequence-aware comparison, use
 *            digitized sequences and the <esl_dst_XPairId()> function
 *            instead.
 *
 * Args:      asq1         - aligned character string 1
 *            asq2         - aligned character string 2
 *            ret_pid      - optRETURN: pairwise identity, 0<=x<=1
 *            ret_nid      - optRETURN: # of identities
 *            ret_n        - optRETURN: denominator MIN(len1,len2)
 *
 * Returns:   <eslOK> on success. <ret_pid>, <ret_nid>, <ret_n>
 *            contain the answers (for whichever were passed non-NULL). 
 *
 * Throws:    <eslEINVAL> if the strings are different lengths
 *            (not aligned).
 */
int
esl_dst_CPairId(char *asq1, char *asq2, 
		double *ret_pid, int *ret_nid, int *ret_n)
{
  int     status;
  int     idents;               /* total identical positions  */
  int     len1, len2;           /* lengths of seqs            */
  int     i;                    /* position in aligned seqs   */

  idents = len1 = len2 = 0;
  for (i = 0; asq1[i] != '\0' && asq2[i] != '\0'; i++) 
    {
      if (isalpha(asq1[i])) len1++;
      if (isalpha(asq2[i])) len2++;
      if (isalpha(asq1[i]) && isalpha(asq2[i])
	  && toupper(asq1[i]) == toupper(asq2[i])) 
	idents++;
    }
  if (asq1[i] != '\0' || asq2[i] != '\0') 
    ESL_XEXCEPTION(eslEINVAL, "strings not same length, not aligned");

  if (ret_pid  != NULL)  *ret_pid = ( len1==0 ? 0. : (double) idents / (double) ESL_MIN(len1,len2));
  if (ret_nid  != NULL)  *ret_nid = idents;
  if (ret_n    != NULL)  *ret_n   = len1;
  return eslOK;

 ERROR:
  if (ret_pid  != NULL)  *ret_pid = 0.;
  if (ret_nid  != NULL)  *ret_nid = 0;
  if (ret_n    != NULL)  *ret_n   = 0;
  return status;
}


/* Function:  esl_dst_CJukesCantor()
 * Synopsis:  Jukes-Cantor distance for two aligned strings.
 * Incept:    SRE, Tue Apr 18 14:00:37 2006 [St. Louis]
 *
 * Purpose:   Calculate the generalized Jukes-Cantor distance between
 *            two aligned character strings <as1> and <as2>, in
 *            substitutions/site, for an alphabet of <K> residues
 *            (<K=4> for nucleic acid, <K=20> for proteins). The
 *            maximum likelihood estimate for the distance is
 *            optionally returned in <ret_distance>. The large-sample
 *            variance for the distance estimate is
 *            returned in <ret_variance>.
 *            
 *            Alphabetic symbols <[a-zA-Z]> are compared
 *            case-insensitively to count the number of identities
 *            (<n1>) and mismatches (<n2>>). Any nonalphabetic
 *            character is assumed to be a gap symbol, and aligned
 *            columns containing gap symbols are ignored.  The
 *            fractional difference <D> used to calculate the
 *            Jukes/Cantor distance is <n2/n1+n2>.
 *            
 * Args:      K            - size of the alphabet (4 or 20)
 *            as1          - 1st aligned seq, 0..L-1, \0-terminated
 *            as2          - 2nd aligned seq, 0..L-1, \0-terminated 
 *            ret_distance - RETURN: ML estimate of distance d
 *            ret_variance - RETURN: large-sample variance of d
 *
 * Returns:   <eslOK> on success.
 * 
 *            Infinite distances are possible, in which case distance
 *            and variance are both <HUGE_VAL>. Caller has to deal
 *            with this case as it sees fit, perhaps by enforcing
 *            an arbitrary maximum distance.
 *
 * Throws:    <eslEINVAL> if the two strings aren't the same length (and
 *            thus can't have been properly aligned).
 *            <eslEDIVZERO> if no aligned residues were counted.
 *            On either failure, distance and variance are both returned
 *            as <HUGE_VAL>.
 */
int
esl_dst_CJukesCantor(int K, char *as1, char *as2, 
		     double *ret_distance, double *ret_variance)
{
  int     status;
  int     n1, n2;               /* number of observed identities, substitutions */
  int     i;                    /* position in aligned seqs   */

  /* 1. Count identities, mismatches.
   */
  n1 = n2 = 0;
  for (i = 0; as1[i] != '\0' && as2[i] != '\0'; i++) 
    {
      if (isalpha(as1[i]) && isalpha(as2[i]))
	{
	  if (toupper(as1[i]) == toupper(as2[i])) n1++; else n2++;
	}
    }
  if (as1[i] != '\0' || as2[i] != '\0') 
    ESL_XEXCEPTION(eslEINVAL, "strings not same length, not aligned");
  
  return jukescantor(n1, n2, K, ret_distance, ret_variance); /* can throw eslEDIVZERO */

 ERROR:
  if (ret_distance != NULL)  *ret_distance = HUGE_VAL;
  if (ret_variance != NULL)  *ret_variance = HUGE_VAL;
  return status;
}
/*------- end, pairwise distances for aligned text seqs ---------*/





/*****************************************************************
 * 2. Pairwise distances for aligned digitized sequences. [alphabet]
 *****************************************************************/
#ifdef eslAUGMENT_ALPHABET

/* Function:  esl_dst_XPairId()
 * Synopsis:  Pairwise identity of two aligned digital seqs.
 * Incept:    SRE, Tue Apr 18 09:24:05 2006 [St. Louis]
 *
 * Purpose:   Digital version of <esl_dst_PairId()>: <adsq1> and
 *            <adsq2> are digitized aligned sequences, in alphabet
 *            <abc>. Otherwise, same as <esl_dst_PairId()>.
 *            
 * Args:      abc          - digital alphabet in use
 *            ax1          - aligned digital seq 1
 *            ax2          - aligned digital seq 2
 *            ret_pid      - RETURN: pairwise identity, 0<=x<=1
 *            ret_nid      - RETURN: # of identities
 *            ret_n        - RETURN: denominator MIN(len1,len2)
 *
 * Returns:   <eslOK> on success. <ret_distance>, <ret_nid>, <ret_n>
 *            contain the answers.
 *
 * Throws:    <eslEINVAL> if the strings are different lengths (not aligned).
 */
int
esl_dst_XPairId(ESL_ALPHABET *abc, ESL_DSQ *ax1, ESL_DSQ *ax2, 
		double *ret_distance, int *ret_nid, int *ret_n)
{
  int     status;
  int     idents;               /* total identical positions  */
  int     len1, len2;           /* lengths of seqs            */
  int     i;                    /* position in aligned seqs   */

  idents = len1 = len2 = 0;
  for (i = 1; ax1[i] != eslDSQ_SENTINEL && ax2[i] != eslDSQ_SENTINEL; i++) 
    {
      if (esl_abc_XIsCanonical(abc, ax1[i])) len1++;
      if (esl_abc_XIsCanonical(abc, ax2[i])) len2++;

      if (esl_abc_XIsCanonical(abc, ax1[i]) && esl_abc_XIsCanonical(abc, ax2[i])
	  && ax1[i] == ax2[i])
	idents++;
    }
  if (len2 < len1) len1 = len2;

  if (ax1[i] != eslDSQ_SENTINEL || ax2[i] != eslDSQ_SENTINEL) 
    ESL_XEXCEPTION(eslEINVAL, "strings not same length, not aligned");

  if (ret_distance != NULL)  *ret_distance = ( len1==0 ? 0. : (double) idents / (double) len1 );
  if (ret_nid      != NULL)  *ret_nid      = idents;
  if (ret_n        != NULL)  *ret_n        = len1;
  return eslOK;

 ERROR:
  if (ret_distance != NULL)  *ret_distance = 0.;
  if (ret_nid      != NULL)  *ret_nid      = 0;
  if (ret_n        != NULL)  *ret_n        = 0;
  return status;
}


/* Function:  esl_dst_XJukesCantor()
 * Synopsis:  Jukes-Cantor distance for two aligned digitized seqs.
 * Incept:    SRE, Tue Apr 18 15:26:51 2006 [St. Louis]
 *
 * Purpose:   Calculate the generalized Jukes-Cantor distance between two
 *            aligned digital strings <ax> and <ay>, in substitutions/site, 
 *            using alphabet <abc> to evaluate identities and differences.
 *            The maximum likelihood estimate for the distance is returned in
 *            <ret_distance>. The large-sample variance for the distance
 *            estimate is optionally returned in <ret_variance>.
 *            
 *            Identical to <esl_dst_CJukesCantor()>, except that it takes
 *            digital sequences instead of character strings.
 *
 * Args:      abc          - bioalphabet to use for comparisons
 *            ax           - 1st digital aligned seq
 *            ay           - 2nd digital aligned seq
 *            ret_distance - RETURN: ML estimate of distance d
 *            ret_variance - RETURN: large-sample variance of d
 *
 * Returns:   <eslOK> on success. As in <esl_dst_CJukesCantor()>, the
 *            distance and variance may be infinite, in which case they
 *            are returned as <HUGE_VAL>.
 *
 * Throws:    <eslEINVAL> if the two strings aren't the same length (and
 *            thus can't have been properly aligned).
 *            <eslEDIVZERO> if no aligned residues were counted.
 *            On either failure, the distance and variance are set
 *            to <HUGE_VAL>.
 */
int
esl_dst_XJukesCantor(ESL_ALPHABET *abc, ESL_DSQ *ax, ESL_DSQ *ay, 
		     double *ret_distance, double *ret_variance)
{
  int     status;
  int     n1, n2;               /* number of observed identities, substitutions */
  int     i;                    /* position in aligned seqs   */

  n1 = n2 = 0;
  for (i = 1; ax[i] != eslDSQ_SENTINEL && ay[i] != eslDSQ_SENTINEL; i++) 
    {
      if (esl_abc_XIsCanonical(abc, ax[i]) && esl_abc_XIsCanonical(abc, ay[i]))
	{
	  if (ax[i] == ay[i]) n1++;
	  else                n2++;
	}
    }
  if (ax[i] != eslDSQ_SENTINEL || ay[i] != eslDSQ_SENTINEL) 
    ESL_XEXCEPTION(eslEINVAL, "strings not same length, not aligned");
  
  return jukescantor(n1, n2, abc->K, ret_distance, ret_variance);

 ERROR:
  if (ret_distance != NULL)  *ret_distance = HUGE_VAL;
  if (ret_variance != NULL)  *ret_variance = HUGE_VAL;
  return status;
}

#endif /*eslAUGMENT_ALPHABET*/
/*---------- end pairwise distances, digital seqs --------------*/




/*****************************************************************
 * 3. Distance matrices for aligned text sequences.
 *****************************************************************/
#ifdef eslAUGMENT_DMATRIX

/* Function:  esl_dst_CPairIdMx()
 * Synopsis:  NxN identity matrix for N aligned text sequences.
 * Incept:    SRE, Thu Apr 27 08:46:08 2006 [New York]
 *
 * Purpose:   Given a multiple sequence alignment <as>, consisting
 *            of <N> aligned character strings; calculate
 *            a symmetric fractional pairwise identity matrix by $N(N-1)/2$
 *            calls to <esl_dst_CPairId()>, and return it in 
 *            <ret_D>.
 *
 * Args:      as      - aligned seqs (all same length), [0..N-1]
 *            N       - # of aligned sequences
 *            ret_S   - RETURN: symmetric fractional identity matrix
 *
 * Returns:   <eslOK> on success, and <ret_S> contains the fractional
 *            identity matrix. Caller free's <S> with
 *            <esl_dmatrix_Destroy()>.
 *
 * Throws:    <eslEINVAL> if a seq has a different
 *            length than others. On failure, <ret_D> is returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_CPairIdMx(char **as, int N, ESL_DMATRIX **ret_S)
{
  ESL_DMATRIX *S = NULL;
  int status;
  int i,j;

  if (( S = esl_dmatrix_Create(N,N) ) == NULL) goto ERROR;
  
  for (i = 0; i < N; i++)
    {
      S->mx[i][i] = 1.;
      for (j = i+1; j < N; j++)
	{
	  status = esl_dst_CPairId(as[i], as[j], &(S->mx[i][j]), NULL, NULL);
	  if (status != eslOK)
	    ESL_XEXCEPTION(status, "Pairwise identity calculation failed at seqs %d,%d\n", i,j);
	  S->mx[j][i] =  S->mx[i][j];
	}
    }
  if (ret_S != NULL) *ret_S = S; else esl_dmatrix_Destroy(S);
  return eslOK;

 ERROR:
  if (S     != NULL)  esl_dmatrix_Destroy(S);
  if (ret_S != NULL) *ret_S = NULL;
  return status;
}


/* Function:  esl_dst_CDiffMx()
 * Synopsis:  NxN difference matrix for N aligned text sequences.
 * Incept:    SRE, Fri Apr 28 06:27:20 2006 [New York]
 *
 * Purpose:   Same as <esl_dst_CPairIdMx()>, but calculates
 *            the fractional difference <d=1-s> instead of the
 *            fractional identity <s> for each pair.
 *
 * Args:      as      - aligned seqs (all same length), [0..N-1]
 *            N       - # of aligned sequences
 *            ret_D   - RETURN: symmetric fractional difference matrix
 *
 * Returns:   <eslOK> on success, and <ret_D> contains the
 *            fractional difference matrix. Caller free's <D> with 
 *            <esl_dmatrix_Destroy()>.
 *
 * Throws:    <eslEINVAL> if any seq has a different
 *            length than others. On failure, <ret_D> is returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_CDiffMx(char **as, int N, ESL_DMATRIX **ret_D)
{
  ESL_DMATRIX *D = NULL;
  int status;
  int i,j;

  status = esl_dst_CPairIdMx(as, N, &D);
  if (status != eslOK) goto ERROR;

  for (i = 0; i < N; i++)
    {
      D->mx[i][i] = 0.;
      for (j = i+1; j < N; j++) 
	{
	  D->mx[i][j] = 1. - D->mx[i][j];
	  D->mx[j][i] = D->mx[i][j];
	}
    }

  if (ret_D != NULL) *ret_D = D; else esl_dmatrix_Destroy(D);
  return eslOK;

 ERROR:
  if (D     != NULL)  esl_dmatrix_Destroy(D);
  if (ret_D != NULL) *ret_D = NULL;
  return status;

}

/* Function:  esl_dst_CJukesCantorMx()
 * Synopsis:  NxN Jukes/Cantor distance matrix for N aligned text seqs.
 * Incept:    SRE, Tue Apr 18 16:00:16 2006 [St. Louis]
 *
 * Purpose:   Given a multiple sequence alignment <aseq>, consisting of
 *            <nseq> aligned character sequences in an alphabet of
 *            <K> letters (usually 4 for DNA, 20 for protein);
 *            calculate a symmetric Jukes/Cantor pairwise distance
 *            matrix for all sequence pairs, and return the distance
 *            matrix in <ret_D>, and a symmetric matrix of the
 *            large-sample variances for those ML distance estimates
 *            in <ret_V>.
 *            
 *            Infinite distances (and variances) are possible; they
 *            are represented as <HUGE_VAL> in <D> and <V>. Caller must
 *            be prepared to deal with them as appropriate.
 *
 * Args:      K      - size of the alphabet (usually 4 or 20)
 *            aseq   - aligned sequences [0.nseq-1][0..L-1]
 *            nseq   - number of aseqs
 *            ret_D  - RETURN: [0..nseq-1]x[0..nseq-1] symmetric distance mx
 *            ret_V  - RETURN: matrix of variances.
 *
 * Returns:   <eslOK> on success. <D> and <V> contain the
 *            distance matrix (and variances); caller frees these with
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslEINVAL> if any pair of sequences have differing lengths
 *            (and thus cannot have been properly aligned). 
 *            <eslEDIVZERO> if some pair of sequences had no aligned
 *            residues. On failure, <D> and <V> are both returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_CJukesCantorMx(int K, char **aseq, int nseq, 
		       ESL_DMATRIX **ret_D, ESL_DMATRIX **ret_V)
{
  int          status;
  ESL_DMATRIX *D = NULL;
  ESL_DMATRIX *V = NULL;
  int          i,j;

  if (( D = esl_dmatrix_Create(nseq, nseq) ) == NULL) goto ERROR;
  if (( V = esl_dmatrix_Create(nseq, nseq) ) == NULL) goto ERROR;

  for (i = 0; i < nseq; i++)
    {
      D->mx[i][i] = 0.;
      V->mx[i][i] = 0.;
      for (j = i+1; j < nseq; j++)
	{
	  status = esl_dst_CJukesCantor(K, aseq[i], aseq[j], 
					&(D->mx[i][j]), &(V->mx[i][j]));
	  if (status != eslOK) 
	    ESL_XEXCEPTION(status, "J/C calculation failed at seqs %d,%d", i,j);

	  D->mx[j][i] = D->mx[i][j];
	  V->mx[j][i] = V->mx[i][j];
	}
    }
  if (ret_D != NULL) *ret_D = D;  else esl_dmatrix_Destroy(D);
  if (ret_V != NULL) *ret_V = V;  else esl_dmatrix_Destroy(V);
  return eslOK;

 ERROR:
  if (D     != NULL) esl_dmatrix_Destroy(D);
  if (V     != NULL) esl_dmatrix_Destroy(V);
  if (ret_D != NULL) *ret_D = NULL;
  if (ret_V != NULL) *ret_V = NULL;
  return status;
}

#endif /*eslAUGMENT_DMATRIX*/
/*----------- end, distance matrices for aligned text seqs ---------*/




/*****************************************************************
 * 4. Distance matrices for aligned digital sequences.
 *****************************************************************/
#if defined(eslAUGMENT_ALPHABET) && defined(eslAUGMENT_DMATRIX)


/* Function:  esl_dst_XPairIdMx()
 * Synopsis:  NxN identity matrix for N aligned digital seqs.
 * Incept:    SRE, Thu Apr 27 09:08:11 2006 [New York]
 *
 * Purpose:   Given a digitized multiple sequence alignment <ax>, consisting
 *            of <N> aligned digital sequences in alphabet <abc>; calculate
 *            a symmetric pairwise fractional identity matrix by $N(N-1)/2$
 *            calls to <esl_dst_XPairId()>, and return it in <ret_S>.
 *            
 * Args:      abc   - digital alphabet in use
 *            ax    - aligned dsq's, [0..N-1][1..alen]                  
 *            N     - number of aligned sequences
 *            ret_S - RETURN: NxN matrix of fractional identities
 *
 * Returns:   <eslOK> on success, and <ret_S> contains the distance
 *            matrix. Caller is obligated to free <S> with 
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslEINVAL> if a seq has a different
 *            length than others. On failure, <ret_S> is returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_XPairIdMx(ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_S)
{
  int status;
  ESL_DMATRIX *S = NULL;
  int i,j;

  if (( S = esl_dmatrix_Create(N,N) ) == NULL) goto ERROR;
  
  for (i = 0; i < N; i++)
    {
      S->mx[i][i] = 1.;
      for (j = i+1; j < N; j++)
	{
	  status = esl_dst_XPairId(abc, ax[i], ax[j], &(S->mx[i][j]), NULL, NULL);
	  if (status != eslOK)
	    ESL_XEXCEPTION(status, "Pairwise identity calculation failed at seqs %d,%d\n", i,j);
	  S->mx[j][i] =  S->mx[i][j];
	}
    }
  if (ret_S != NULL) *ret_S = S; else esl_dmatrix_Destroy(S);
  return eslOK;

 ERROR:
  if (S     != NULL)  esl_dmatrix_Destroy(S);
  if (ret_S != NULL) *ret_S = NULL;
  return status;
}


/* Function:  esl_dst_XDiffMx()
 * Synopsis:  NxN difference matrix for N aligned digital seqs.         
 * Incept:    SRE, Fri Apr 28 06:37:29 2006 [New York]
 *
 * Purpose:   Same as <esl_dst_XPairIdMx()>, but calculates fractional
 *            difference <1-s> instead of fractional identity <s> for
 *            each pair.
 *
 * Args:      abc   - digital alphabet in use
 *            ax    - aligned dsq's, [0..N-1][1..alen]                  
 *            N     - number of aligned sequences
 *            ret_D - RETURN: NxN matrix of fractional differences
 *            
 * Returns:   <eslOK> on success, and <ret_D> contains the difference
 *            matrix; caller is obligated to free <D> with 
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslEINVAL> if a seq has a different
 *            length than others. On failure, <ret_D> is returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_XDiffMx(ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_D)
{
  int status;
  ESL_DMATRIX *D = NULL;
  int i,j;

  status = esl_dst_XPairIdMx(abc, ax, N, &D);
  if (status != eslOK) goto ERROR;

  for (i = 0; i < N; i++)
    {
      D->mx[i][i] = 0.;
      for (j = i+1; j < N; j++) 
	{
	  D->mx[i][j] = 1. - D->mx[i][j];
	  D->mx[j][i] = D->mx[i][j];
	}
    }
  if (ret_D != NULL) *ret_D = D; else esl_dmatrix_Destroy(D);
  return eslOK;

 ERROR:
  if (D     != NULL)  esl_dmatrix_Destroy(D);
  if (ret_D != NULL) *ret_D = NULL;
  return status;
}

/* Function:  esl_dst_XJukesCantorMx()
 * Synopsis:  NxN Jukes/Cantor distance matrix for N aligned digital seqs.
 * Incept:    SRE, Thu Apr 27 08:38:08 2006 [New York City]
 *
 * Purpose:   Given a digitized multiple sequence alignment <ax>,
 *            consisting of <nseq> aligned digital sequences in
 *            bioalphabet <abc>, calculate a symmetric Jukes/Cantor
 *            pairwise distance matrix for all sequence pairs;
 *            return the distance matrix in <ret_D> and 
 *            a matrix of the large-sample variances for those ML distance
 *            estimates in <ret_V>.
 *            
 *            Infinite distances (and variances) are possible. They
 *            are represented as <HUGE_VAL> in <D> and <V>. Caller must
 *            be prepared to deal with them as appropriate.
 *
 * Args:      abc    - bioalphabet for <aseq>
 *            ax     - aligned digital sequences [0.nseq-1][1..L]
 *            nseq   - number of aseqs
 *            ret_D  - RETURN: [0..nseq-1]x[0..nseq-1] symmetric distance mx
 *            ret_V  - RETURN: matrix of variances.
 *
 * Returns:   <eslOK> on success. <D> (and optionally <V>) contain the
 *            distance matrix (and variances). Caller frees these with
 *            <esl_dmatrix_Destroy()>. 
 *
 * Throws:    <eslEINVAL> if any pair of sequences have differing lengths
 *            (and thus cannot have been properly aligned). 
 *            <eslEDIVZERO> if some pair of sequences had no aligned
 *            residues. On failure, <D> and <V> are both returned <NULL>
 *            and state of inputs is unchanged.
 */
int
esl_dst_XJukesCantorMx(ESL_ALPHABET *abc, ESL_DSQ **ax, int nseq, 
		       ESL_DMATRIX **ret_D, ESL_DMATRIX **ret_V)
{
  ESL_DMATRIX *D = NULL;
  ESL_DMATRIX *V = NULL;
  int          status;
  int          i,j;

  if (( D = esl_dmatrix_Create(nseq, nseq) ) == NULL) goto ERROR;
  if (( V = esl_dmatrix_Create(nseq, nseq) ) == NULL) goto ERROR;

  for (i = 0; i < nseq; i++)
    {
      D->mx[i][i] = 0.;
      V->mx[i][i] = 0.;
      for (j = i+1; j < nseq; j++)
	{
	  status = esl_dst_XJukesCantor(abc, ax[i], ax[j], 
					&(D->mx[i][j]), &(V->mx[i][j]));
	  if (status != eslOK) 
	    ESL_XEXCEPTION(status, "J/C calculation failed at digital aseqs %d,%d", i,j);

	  D->mx[j][i] = D->mx[i][j];
	  V->mx[j][i] = V->mx[i][j];
	}
    }
  if (ret_D != NULL) *ret_D = D;  else esl_dmatrix_Destroy(D);
  if (ret_V != NULL) *ret_V = V;  else esl_dmatrix_Destroy(V);
  return eslOK;

 ERROR:
  if (D     != NULL) esl_dmatrix_Destroy(D);
  if (V     != NULL) esl_dmatrix_Destroy(V);
  if (ret_D != NULL) *ret_D = NULL;
  if (ret_V != NULL) *ret_V = NULL;
  return status;
}
#endif /*eslAUGMENT_ALPHABET && eslAUGMENT_DMATRIX*/
/*------- end, distance matrices for digital alignments ---------*/



/*****************************************************************
 * 5. Private (static) functions
 *****************************************************************/

/* jukescantor()
 * 
 * The generalized Jukes/Cantor distance calculation.
 * Given <n1> identities and <n2> differences, for a
 * base alphabet size of <alphabet_size> (4 or 20);
 * calculate J/C distance in substitutions/site and
 * return it in <ret_distance>; calculate large-sample
 * variance and return it in <ret_variance>.
 *
 * Returns <eslEDIVZERO> if there are no data (<n1+n2=0>).
 */
static int
jukescantor(int n1, int n2, int alphabet_size, double *ret_distance, double *ret_variance)
{
  int    status;
  double D, K, N;
  double x;
  double distance, variance;

  ESL_DASSERT1( (n1 >= 0) );
  ESL_DASSERT1( (n2 >= 0) );
  ESL_DASSERT1( (alphabet_size >= 0) );

  if (n1+n2 == 0) { status = eslEDIVZERO; goto ERROR; }

  K = (double) alphabet_size;
  D = (double) n2 / (double) (n1+n2);
  N = (double) (n1+n2);

  x = 1. - D * K/(K-1.);
  if (x <= 0.) 
    {
      distance = HUGE_VAL;
      variance = HUGE_VAL;
    }
  else
    {
      distance =   -log(x) * K/(K-1);
      variance =  exp( 2.*K*distance/(K-1) ) * D * (1.-D) / N;
    }
  if (ret_distance != NULL)  *ret_distance = distance;
  if (ret_variance != NULL)  *ret_variance = variance;
  return eslOK;

 ERROR:
  if (ret_distance != NULL)  *ret_distance = HUGE_VAL;
  if (ret_variance != NULL)  *ret_variance = HUGE_VAL;
  return status;
}
/*--------------- end of private functions ----------------------*/


/*****************************************************************
 * 6. Unit tests.
 *****************************************************************/ 
#ifdef eslDISTANCE_TESTDRIVE

/* Each unit test is given an alignment with certain known
 * properties:
 *    seqs 0,1 are identical
 *    seqs 0,2 are completely different
 *    seqs 3..N are random
 * The alignment may contain gaps, so don't assume that the 
 * # of compared residues == alignment length. The alignment
 * contains only canonical residues, because one of our tests
 * is that C and X functions give the same results.
 */
static int 
utest_CPairId(char **as, int N)
{
  double pid;
  int    nid;
  int    nres;	
  int    L;
  int    i,j;

  /* Self comparison gives identity = 1. */
  L = strlen(as[0]);
  if (esl_dst_CPairId(as[0], as[0], &pid, &nid, &nres) != eslOK) abort();
  if (pid  != 1.0 || nid != L || nres > L) abort();

  /* So does 0,1 comparison  */
  if (esl_dst_CPairId(as[0], as[1], &pid, &nid, &nres) != eslOK) abort();
  if (pid  != 1.0 || nid != L || nres > L) abort();

  /* 0,2 comparison gives 0.0, 0 */
  if (esl_dst_CPairId(as[0], as[2], &pid, &nid, &nres) != eslOK) abort();
  if (pid  != 0.0 || nid != 0 || nres > L) abort();
  
  /* remaining comparisons shouldn't fail */
  for (i = 3; i < N; i++)
    for (j = i; j < N; j++)
      {
	if (esl_dst_CPairId(as[i], as[j], &pid, &nid, &nres) != eslOK) abort();
	if (pid < 0. || pid > 1. || nid < 0 || nid > L || nres > L)    abort();
      }

  /* API should accept NULL for return values */
  if (esl_dst_CPairId(as[0], as[0], NULL, NULL, NULL) != eslOK) abort();  
  return eslOK;
}

static int 
utest_CJukesCantor(int K, char **as, int N)
{
  double d, V;
  int    i,j;

  /* Self comparison gives distance = 0. */
  if (esl_dst_CJukesCantor(K, as[0], as[0], &d, &V) != eslOK) abort();
  if (d != 0.0) abort();

  /* So does 0,1 comparison  */
  if (esl_dst_CJukesCantor(K, as[0], as[1], &d, &V) != eslOK) abort();
  if (d != 0.0) abort();

  /* 0,2 comparison gives infinite distance (HUGE_VAL) */
  if (esl_dst_CJukesCantor(K, as[0], as[2], &d, &V) != eslOK) abort();
  if (d != HUGE_VAL) abort();
  
  /* remaining comparisons shouldn't fail */
  for (i = 3; i < N; i++)
    for (j = i; j < N; j++)
      if (esl_dst_CJukesCantor(K, as[i], as[j], &d, &V) != eslOK) abort();

  /* API should accept NULL for return values */
  if (esl_dst_CJukesCantor(K, as[0], as[0], NULL, NULL) != eslOK) abort();  
  return eslOK;
}

#ifdef eslAUGMENT_ALPHABET
static int 
utest_XPairId(ESL_ALPHABET *abc, char **as, ESL_DSQ **ax, int N)
{
  double pid, pid2;
  int    nid, nid2;
  int    nres, nres2;	
  int    dL, L;
  int    i,j;

  /* Self comparison gives identity = 1. */
  dL = esl_abc_dsqlen(ax[0]);
  L  = strlen(as[0]); 
  if (dL != L) abort();
  if (esl_dst_XPairId(abc, ax[0], ax[0], &pid, &nid, &nres) != eslOK) abort();
  if (pid  != 1.0 || nid != L || nres > dL) abort();

  /* So does 0,1 comparison  */
  if (esl_dst_XPairId(abc, ax[0], ax[1], &pid, &nid, &nres) != eslOK) abort();
  if (pid  != 1.0 || nid != L || nres > L) abort();

  /* 0,2 comparison gives 0.0, 0 */
  if (esl_dst_XPairId(abc, ax[0], ax[2], &pid, &nid, &nres) != eslOK) abort();
  if (pid  != 0.0 || nid != 0 || nres > L) abort();
  
  /* remaining comparisons shouldn't fail, and should be identical to text mode */
  for (i = 3; i < N; i++)
    for (j = i; j < N; j++)
      {
	if (esl_dst_XPairId(abc, ax[i], ax[j], &pid, &nid, &nres) != eslOK) abort();
	if (esl_dst_CPairId(as[i], as[j], &pid2, &nid2, &nres2)   != eslOK) abort();
	if (pid < 0. || pid > 1. || nid < 0 || nid > L || nres > L)         abort();
	if (pid != pid2 || nid != nid2 || nres != nres2)                    abort();
      }

  /* API should accept NULL for return values */
  if (esl_dst_XPairId(abc, ax[0], ax[0], NULL, NULL, NULL) != eslOK) abort();  
  return eslOK;

}

static int 
utest_XJukesCantor(ESL_ALPHABET *abc, char **as, ESL_DSQ **ax, int N)
{
  double d, V;
  int    i,j;

  /* Self comparison gives distance = 0. */
  if (esl_dst_XJukesCantor(abc, ax[0], ax[0], &d, &V) != eslOK) abort();
  if (d != 0.0) abort();

  /* So does 0,1 comparison  */
  if (esl_dst_XJukesCantor(abc, ax[0], ax[1], &d, &V) != eslOK) abort();
  if (d != 0.0) abort();

  /* 0,2 comparison gives infinite distance (HUGE_VAL) */
  if (esl_dst_XJukesCantor(abc, ax[0], ax[2], &d, &V) != eslOK) abort();
  if (d != HUGE_VAL) abort();
  
  /* remaining comparisons shouldn't fail */
  for (i = 3; i < N; i++)
    for (j = i; j < N; j++)
      if (esl_dst_XJukesCantor(abc, ax[i], ax[j], &d, &V) != eslOK) abort();

  /* API should accept NULL for return values */
  if (esl_dst_XJukesCantor(abc, ax[0], ax[0], NULL, NULL) != eslOK) abort();  
  return eslOK;

}
#endif /*eslAUGMENT_ALPHABET*/


#ifdef eslAUGMENT_DMATRIX
static int 
utest_CPairIdMx(char **as, int N)
{
  ESL_DMATRIX *S;
  int          i,j;
  double       pid;

  if (esl_dst_CPairIdMx(as, N, &S) != eslOK) abort();

  for (i = 0; i < N; i++) 
    if (S->mx[i][i] != 1.0) abort();

  pid = 0.;
  for (i = 3; i < N; i++)
    for (j = i+1; j < N; j++)
      pid += S->mx[i][j];
  pid /= (double) ((N-3) * (N-4) / 2); /* first 3 don't count */
  if (pid < 0.15 || pid > 0.35) abort(); /* should be 0.25 */

  esl_dmatrix_Destroy(S);
  return eslOK;
}

static int 
utest_CDiffMx(char **as, int N)
{
  ESL_DMATRIX *D;
  int          i,j;
  double       diff;

  if (esl_dst_CDiffMx(as, N, &D) != eslOK) abort();

  for (i = 0; i < N; i++) 
    if (D->mx[i][i] != 0.0) abort();

  diff = 0.;
  for (i = 3; i < N; i++)
    for (j = i+1; j < N; j++)
      diff += D->mx[i][j];
  diff /= (double) ((N-3) * (N-4) / 2);	/* first 3 don't count */
  if (diff < 0.65 || diff > 0.85) abort(); /* should be 0.75 */

  esl_dmatrix_Destroy(D);
  return eslOK;
}

static int 
utest_CJukesCantorMx(int K, char **as, int N)
{
  ESL_DMATRIX *D, *V;
  /* just a crash test */
  if (esl_dst_CJukesCantorMx(K, as, N, &D, &V) != eslOK) abort();
  esl_dmatrix_Destroy(D);
  esl_dmatrix_Destroy(V);
  return eslOK;
}
#endif /*eslAUGMENT_DMATRIX*/

#if defined (eslAUGMENT_ALPHABET) && defined (eslAUGMENT_DMATRIX)
 
static int 
utest_XPairIdMx(ESL_ALPHABET *abc, char **as, ESL_DSQ **ax, int N)
{
  ESL_DMATRIX *S, *S2;
  int i, j;

  if (esl_dst_XPairIdMx(abc, ax, N, &S) != eslOK) abort();
  if (esl_dst_CPairIdMx(as, N, &S2)     != eslOK) abort();

  for (i = 0; i < N; i++) 
    for (j = i; j < N; j++)
      if (fabs(S->mx[i][j] - S2->mx[j][i]) > 0.01) abort();

  esl_dmatrix_Destroy(S);
  esl_dmatrix_Destroy(S2);
  return eslOK;
}

static int 
utest_XDiffMx(ESL_ALPHABET *abc, char **as, ESL_DSQ **ax, int N)
{
  ESL_DMATRIX *D, *D2;
  int i, j;

  if (esl_dst_XDiffMx(abc, ax, N, &D) != eslOK) abort();
  if (esl_dst_CDiffMx(as, N, &D2)     != eslOK) abort();

  for (i = 0; i < N; i++) 
    for (j = i; j < N; j++)
      if (fabs(D->mx[i][j] - D2->mx[j][i]) > 0.01) abort();

  esl_dmatrix_Destroy(D);
  esl_dmatrix_Destroy(D2);
  return eslOK;
}

static int 
utest_XJukesCantorMx(ESL_ALPHABET *abc, char **as, ESL_DSQ **ax, int N)
{
  ESL_DMATRIX *D, *D2, *V, *V2;
  int i, j;

  if (esl_dst_XJukesCantorMx(abc, ax, N, &D, &V)      != eslOK) abort();
  if (esl_dst_CJukesCantorMx(abc->K, as, N, &D2, &V2) != eslOK) abort();

  for (i = 0; i < N; i++) 
    for (j = i; j < N; j++)
      {
	if (fabs(D->mx[i][j] - D2->mx[j][i]) > 0.01) abort();
	if (fabs(V->mx[i][j] - V2->mx[j][i]) > 0.01) abort();
      }

  esl_dmatrix_Destroy(D);
  esl_dmatrix_Destroy(D2);
  esl_dmatrix_Destroy(V);
  esl_dmatrix_Destroy(V2);
  return eslOK;
}
#endif /*eslAUGMENT_ALPHABET && eslAUGMENT_DMATRIX*/

/*------------------ end of unit tests --------------------------*/
#endif /* eslDISTANCE_TESTDRIVE */



/*****************************************************************
 * Test driver.
 *****************************************************************/ 

/* 
   gcc -g -Wall -o testdriver -I. -DeslDISTANCE_TESTDRIVE esl_distance.c esl_getopts.c esl_random.c easel.c -lm

   gcc -g -Wall -o testdriver -I. -DeslDISTANCE_TESTDRIVE\
     -DeslAUGMENT_ALPHABET \
     esl_distance.c esl_getopts.c esl_random.c \
     esl_alphabet.c easel.c -lm

   gcc -g -Wall -o testdriver -I. -DeslDISTANCE_TESTDRIVE\
     -DeslAUGMENT_DMATRIX \
     esl_distance.c esl_getopts.c esl_random.c \
     esl_dmatrix.c easel.c -lm

   gcc -g -Wall -o testdriver -I. -DeslDISTANCE_TESTDRIVE\
     -DeslALPHABET -DeslAUGMENT_DMATRIX \
     esl_distance.c esl_getopts.c esl_random.c \
     esl_alphabet.c esl_dmatrix.c easel.c -lm

   gcc -g -Wall -o testdriver -I. -L. -DeslDISTANCE_TESTDRIVE esl_distance.c -leasel -lm
 */
#ifdef eslDISTANCE_TESTDRIVE
#include <easel.h>
#include <esl_getopts.h>
#include <esl_random.h>
#include <esl_distance.h>
#ifdef eslAUGMENT_ALPHABET
#include <esl_alphabet.h>
#endif
#ifdef eslAUGMENT_DMATRIX
#include <esl_dmatrix.h>
#endif

static ESL_OPTIONS options[] = {
  /* name        type       def   env  range toggles reqs incomp help                       docgroup*/
  { "-h",     eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",            0},
  { "-N",     eslARG_INT,    "10", NULL,"n>3", NULL, NULL, NULL, "number of iid seqs in alignment",0},
  { "-L",     eslARG_INT,    "50", NULL,"n>0", NULL, NULL, NULL, "length of seqs in alignment",    0},
  { "--seed", eslARG_INT,    "42", NULL,"n>0", NULL, NULL, NULL, "random # seed",                  0},
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[] = "Usage: ./testdrive-distance [-options]";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS   *go = NULL;
  ESL_RANDOMNESS *r = NULL;
  char  **as = NULL;		/* aligned character seqs (random, iid) */
  int     N,L;			/* # of seqs, and their aligned lengths */
  int show_help;
  int seed;
  int i,j;
  int status;
  double p[4];			/* ACGT probabilities */
#ifdef eslAUGMENT_ALPHABET
  ESL_DSQ      **ax = NULL;		/* digitized alignment                  */
  ESL_ALPHABET *abc = NULL;
#endif

  /* Process command line
   */
  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);
  esl_opt_GetBooleanOption(go, "-h", &show_help);
  if (show_help) {
    puts(usage); 
    puts("\n  where options are:");
    esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=all docgroups; 2=indentation; 80=width */
    return 0;
  }
  esl_opt_GetIntegerOption(go, "-L",     &L);
  esl_opt_GetIntegerOption(go, "-N",     &N);
  esl_opt_GetIntegerOption(go, "--seed", &seed);
  if (esl_opt_ArgNumber(go) != 0) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return 1;
  }
  esl_getopts_Destroy(go);

  /* Create a random DNA alignment;
   * force it to obey the conventions of the unit tests:
   *   0,1 are identical
   *   0,2 are completely dissimilar
   */
  r   = esl_randomness_Create(seed);
  for (i = 0; i < 4; i++) p[i] = 0.25;
  ESL_ALLOC(as, sizeof(char *) * N);
  for (i = 0; i < N; i++) 
    ESL_ALLOC(as[i], sizeof(char) * (L+1));
  esl_rnd_IID(r, "ACGT", p, 4, L, as[0]);
  strcpy(as[1], as[0]);
  esl_rnd_IID(r, "ACGT", p, 4, L, as[2]);
  for (j = 0; j < L; j++)
    while (as[2][j] == as[0][j])
      as[2][j] = "ACGT"[esl_rnd_Choose(r, 4)];
  for (i = 3; i < N; i++)
    esl_rnd_IID(r, "ACGT", p, 4, L, as[i]);

#ifdef eslAUGMENT_ALPHABET
  abc = esl_alphabet_Create(eslDNA);
  ESL_ALLOC(ax, sizeof(ESL_DSQ *) * N);
  for (i = 0; i < N; i++) 
    esl_abc_CreateDsq(abc, as[i], &(ax[i]));
#endif /*eslAUGMENT_ALPHABET*/


  /* Unit tests
   */
  if (utest_CPairId(as, N)               != eslOK) return eslFAIL;
  if (utest_CJukesCantor(4, as, N)       != eslOK) return eslFAIL;

#ifdef eslAUGMENT_ALPHABET
  if (utest_XPairId(abc, as, ax, N)      != eslOK) return eslFAIL;
  if (utest_XJukesCantor(abc, as, ax, N) != eslOK) return eslFAIL;
#endif /*eslAUGMENT_ALPHABET*/

#ifdef eslAUGMENT_DMATRIX
  if (utest_CPairIdMx(as, N)             != eslOK) return eslFAIL;
  if (utest_CDiffMx(as, N)               != eslOK) return eslFAIL;
  if (utest_CJukesCantorMx(4, as, N)     != eslOK) return eslFAIL;
#endif /* eslAUGMENT_DMATRIX*/

#if defined (eslAUGMENT_ALPHABET) && defined (eslAUGMENT_DMATRIX)
  if (utest_XPairIdMx(abc, as, ax, N)       != eslOK) return eslFAIL;
  if (utest_XDiffMx(abc, as, ax, N)         != eslOK) return eslFAIL;
  if (utest_XJukesCantorMx(abc, as, ax, N)  != eslOK) return eslFAIL;
#endif

  esl_randomness_Destroy(r);
  esl_Free2D((void **) as, N);
#ifdef eslAUGMENT_ALPHABET
  esl_alphabet_Destroy(abc);
  esl_Free2D((void **) ax, N);
#endif
  return eslOK;

 ERROR:
  return eslFAIL;
}
#endif /*eslDISTANCE_TESTDRIVE*/




/*****************************************************************
 * 8. Example.
 *****************************************************************/ 

#ifdef eslDISTANCE_EXAMPLE
/*::cexcerpt::distance_example::begin::*/
/* gcc -g -Wall -o example -I. -DeslDISTANCE_EXAMPLE esl_distance.c\
       esl_dmatrix.c esl_msa.c easel.c -lm
   ./example <msa file>
 */
#include <easel.h>
#include <esl_distance.h>
#include <esl_dmatrix.h>
#include <esl_msa.h>

int main(int argc, char **argv)
{
  ESL_MSAFILE *afp; 
  ESL_MSA     *msa;
  ESL_DMATRIX *P;
  int          status;
  int          i,j;
  double       min, avg, max;

  esl_msafile_Open(argv[1], eslMSAFILE_UNKNOWN, NULL, &afp);
  esl_msa_Read(afp, &msa);

  esl_dst_CPairIdMx(msa->aseq, msa->nseq, &P);

  min = 1.0;
  max = 0.0;
  avg = 0.0;
  for (i = 0; i < msa->nseq; i++)
    for (j = i+1; j < msa->nseq; j++)
      {
	avg += P->mx[i][j];
	if (P->mx[i][j] < min) min = P->mx[i][j];
	if (P->mx[i][j] > max) max = P->mx[i][j];
      }
  avg /= (double) (msa->nseq * (msa->nseq-1) / 2);

  printf("Average pairwise %% id:  %.1f%%\n", avg * 100.);
  printf("Minimum pairwise %% id:  %.1f%%\n", min * 100.);
  printf("Maximum pairwise %% id:  %.1f%%\n", max * 100.);

  esl_dmatrix_Destroy(P);
  esl_msa_Destroy(msa);
  esl_msafile_Close(afp);
  return 0;
}
/*::cexcerpt::distance_example::end::*/
#endif /*eslDISTANCE_EXAMPLE*/


 /*****************************************************************
  *    This copyrighted source code is freely distributed 
  *    under the terms of the GNU General Public License. See
  *    the files COPYRIGHT and LICENSE for details.
  *****************************************************************/


