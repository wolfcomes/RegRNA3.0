/* esl_vectorops.c
 * Operations on vectors of floats or doubles.
 * 
 * Can operate on vectors of doubles, floats, or integers - appropriate
 * routine is prefixed with D, F, or I. For example, esl_vec_DSet() is
 * the Set routine for a vector of doubles; esl_vec_ISet() is for integers.
 * 
 * esl_vec_{D,F,I}Set()         - set all items in vector to value.
 * esl_vec_{D,F,I}Scale()       - multiply all items in vector by scale
 * esl_vec_{D,F,I}Increment()   - add a scalar to all items in vector
 * esl_vec_{D,F,I}Sum()         - return sum of values in vector
 * esl_vec_{D,F,I}Add()         - add vec2 to vec1.
 * esl_vec_{D,F,I}Copy()        - set vec1 to be same as vec2. 
 * esl_vec_{D,F,I}Compare()     - compare vec1 to vec2 for equality.
 * esl_vec_{D,F,I}Dot()         - return dot product of two vectors.
 * esl_vec_{D,F,I}Max()         - return value of maximum element in vector
 * esl_vec_{D,F,I}Min()         - return value of minimum element in vector 
 * esl_vec_{D,F,I}ArgMax()      - return index of maximum element in vector
 * esl_vec_{D,F,I}ArgMin()      - return index of minimum element in vector
 * esl_vec_{D,F,I}SortIncreasing() - sort from smallest to largest
 * esl_vec_{D,F,I}SortDecreasing() - sort from largest to smallest
 * esl_vec_{D,F,I}Dump()        - dump vector to an output stream
 *
 * esl_vec_D2F()                - copy double vector to float vector
 * esl_vec_F2D()                - copy float vector to double vector
 * esl_vec_I2F()                - copy integer vector to float vector
 * esl_vec_I2D()                - copy integer vector to double vector
 * 
 * esl_vec_{D,F}Norm()          - normalize a probability vector of length n.
 * esl_vec_{D,F}Log()           - convert to log probabilities 
 * esl_vec_{D,F}Entropy()       - return Shannon entropy of probability vector, in bits
 *
 * esl_vec_{D,F}Exp()           - convert log p's back to probabilities
 * esl_vec_{D,F}LogSum()        - given vector of log p's; return log of summed p's.
 * esl_vec_{D,F}LogNorm()       - given vec of unnormalized log p's; normalize, make it a p vector
 *                        
 * esl_vec_{D,F}Validate()      - validate a probability vector
 * esl_vec_{D,F}LogValidate()   - validate a log probability vector
 *                        
 * SRE, Tue Oct  1 15:23:25 2002 [St. Louis]
 * SVN $Id: esl_vectorops.c 160 2007-03-23 13:02:04Z eddys $
 */                      

#include <esl_config.h>

#include <math.h>
#include <float.h>

#include <easel.h>
#include <esl_vectorops.h>

/* Function:  esl_vec_DSet()
 *            
 * Purpose:   Sets all <n> items in <vec> to <value>.
 *                        
 *            <esl_vec_FSet()> and <esl_vec_ISet()> do the same,
 *            for float and integer vectors.
 */
void
esl_vec_DSet(double *vec, int n, double value)
{
  int x; 
  for (x = 0; x < n; x++) vec[x] = value;
}
void
esl_vec_FSet(float *vec, int n, float value)
{
  int x; 
  for (x = 0; x < n; x++) vec[x] = value;
}
void
esl_vec_ISet(int *vec, int n, int value)
{
  int x; 
  for (x = 0; x < n; x++) vec[x] = value;
}


/* Function:  esl_vec_DScale()
 *            
 * Purpose:   Multiplies all <n> items in <vec> by <scale>.
 *            
 *            <esl_vec_FScale()> and <esl_vec_IScale()> do the same,
 *            for float and integer vectors.
 *            
 *            Essentially the same as BLAS1's xSCAL().
 */
void
esl_vec_DScale(double *vec, int n, double scale)
{
  int x;
  for (x = 0; x < n; x++) vec[x] *= scale;
}
void
esl_vec_FScale(float *vec, int n, float scale)
{
  int x;
  for (x = 0; x < n; x++) vec[x] *= scale;
}
void
esl_vec_IScale(int *vec, int n, int scale)
{
  int x;
  for (x = 0; x < n; x++) vec[x] *= scale;
}


/* Function:  esl_vec_DIncrement()
 * Incept:    SRE, Mon Mar 21 11:56:57 2005 [St. Louis]
 *
 * Purpose:   Adds scalar <x> to all items in the <n>-vector <v>.
 * 
 *            <esl_vec_FIncrement()> and <esl_vec_IIncrement()> do the
 *            same, for float and integer vectors.
 */
void
esl_vec_DIncrement(double *v, int n, double x)
{
  int i;
  for (i = 0; i < n; i++) v[i] += x;
}
void
esl_vec_FIncrement(float *v, int n, float x)
{
  int i;
  for (i = 0; i < n; i++) v[i] += x;
}
void
esl_vec_IIncrement(int *v, int n, int x)
{
  int i;
  for (i = 0; i < n; i++) v[i] += x;
}



/* Function:  esl_vec_DSum()
 *            
 * Purpose:   Returns the scalar sum of the <n> items in <vec>.
 *            
 *            <esl_vec_FSum()> and <esl_vec_ISum()> do the same,
 *            but for float and integer vectors.
 */
double 
esl_vec_DSum(double *vec, int n)
{
  double sum = 0.;
  int    x;
  for (x = 0; x < n; x++) sum += vec[x];
  return sum;
}
float 
esl_vec_FSum(float *vec, int n)
{
  float sum = 0.;
  int   x;
  for (x = 0; x < n; x++) sum += vec[x];
  return sum;
}
int
esl_vec_ISum(int *vec, int n)
{
  int sum = 0;
  int   x;
  for (x = 0; x < n; x++) sum += vec[x];
  return sum;
}


/* Function:  esl_vec_DAdd()
 *
 * Purpose:   Vector addition. Adds <vec2> to <vec1>, leaving
 *            result in <vec1>. (<vec2> is unchanged.). 
 *            Both vectors are of size <n>.
 *            
 *            <esl_vec_FAdd()> and <esl_vec_IAdd()> do the same,
 *            for float and integer vectors.
 */
void
esl_vec_DAdd(double *vec1, double *vec2, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] += vec2[x];
}
void
esl_vec_FAdd(float *vec1, float *vec2, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] += vec2[x];
}
void
esl_vec_IAdd(int *vec1, int *vec2, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] += vec2[x];
}


/* Function: esl_vec_DAddScaled()
 * 
 * Purpose:  Scales <vec2> by scalar <a>, and adds that
 *           to <vec1>. Both vectors are of size <n>. 
 *           
 *            <esl_vec_FAddScaled()> and <esl_vec_IAddScaled()> do the same,
 *            for float and integer vectors.
 *            
 *            Essentially the same as BLAS1 xAXPY().
 */
void
esl_vec_DAddScaled(double *vec1, double *vec2, double a, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] += vec2[x] * a;
}
void
esl_vec_FAddScaled(float *vec1, float *vec2, float a, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] += vec2[x] * a;
}
void
esl_vec_IAddScaled(int *vec1, int *vec2, int a, int n)
{
  int x;
  for (x = 0; x < n; x++) vec1[x] += vec2[x] * a;
}



/* Function:  esl_vec_DCopy()
 *
 * Purpose:   Copies <src> to <dest>. <src> is
 *            unchanged. Both vectors are of size <n>.
 *            
 *            <esl_vec_FCopy()> and <esl_vec_ICopy()> do the same,
 *            for float and integer vectors.
 *            
 *            Essentially the same as BLAS1 xCOPY().
 */
void
esl_vec_DCopy(const double *src, const int n, double *dest)
{
  int x;
  for (x = 0; x < n; x++) dest[x] = src[x];
}
void
esl_vec_FCopy(const float *src, const int n, float *dest)
{
  int x;
  for (x = 0; x < n; x++) dest[x] = src[x];
}
void
esl_vec_ICopy(const int *src, const int n, int *dest)
{
  int x;
  for (x = 0; x < n; x++) dest[x] = src[x];
}


/* Function:  esl_vec_DCompare()
 * Incept:    SRE, Mon Nov  6 10:20:28 2006 [Janelia]
 *
 * Purpose:   Compare <vec1> to <vec2> for equality, by
 *            comparing each cognate element pair. Both vectors 
 *            are of size <n>. Equality of elements is
 *            defined by being $\leq$ fractional tolerance <tol> 
 *            for floating point comparisons, and strict equality
 *            for integer comparisons. Return <eslOK>
 *            if the vectors are equal, and <eslFAIL> if not.
 *
 *            <esl_vec_FCompare()> and <esl_vec_ICompare()> do the same,
 *            for float and integer vectors.
 */
int
esl_vec_DCompare(double *vec1, double *vec2, int n, double tol)
{
  int i;
  for (i = 0; i < n; i++) if (esl_DCompare(vec1[i], vec2[i], tol) == eslFAIL) return eslFAIL;
  return eslOK;
}
int
esl_vec_FCompare(float *vec1, float *vec2, int n, float tol)
{
  int i;
  for (i = 0; i < n; i++) if (esl_DCompare(vec1[i], vec2[i], tol) == eslFAIL) return eslFAIL;
  return eslOK;
}
int
esl_vec_ICompare(int *vec1, int *vec2, int n)
{
  int i;
  for (i = 0; i < n; i++) if (vec1[i] != vec2[i]) return eslFAIL;
  return eslOK;
}



/* Function:  esl_vec_DSwap()
 *
 * Purpose:   Swaps <vec2> and <vec1>. 
 *            Both vectors are of size <n>.
 *            
 *            <esl_vec_FSwap()> and <esl_vec_ISwap()> do the same,
 *            for float and integer vectors.
 *            
 *            Essentially the same as BLAS1 xSWAP().
 */
void
esl_vec_DSwap(double *vec1, double *vec2, int n)
{
  int    x;
  double tmp;

  for (x = 0; x < n; x++) 
    { tmp = vec1[x]; vec1[x] = vec2[x]; vec2[x] = tmp; }
}
void
esl_vec_FSwap(float *vec1, float *vec2, int n)
{
  int   x;
  float tmp;

  for (x = 0; x < n; x++) 
    { tmp = vec1[x]; vec1[x] = vec2[x]; vec2[x] = tmp; }
}
void
esl_vec_ISwap(int *vec1, int *vec2, int n)
{
  int    x;
  int tmp;

  for (x = 0; x < n; x++) 
    { tmp = vec1[x]; vec1[x] = vec2[x]; vec2[x] = tmp; }
}




/* Function:  esl_vec_DDot()
 *
 * Purpose:   Returns the scalar dot product <vec1> $\cdot$ <vec2>.
 *            Both vectors are of size <n>.
 *            
 *            <esl_vec_FDot()> and <esl_vec_IDot()> do the same,
 *            for float and integer vectors.
 */
double
esl_vec_DDot(double *vec1, double *vec2, int n)
{
  double result = 0.;
  int x;
  for (x = 0; x < n; x++) result += vec1[x] * vec2[x];
  return result;
}
float
esl_vec_FDot(float *vec1, float *vec2, int n)
{
  float result = 0.;
  int x;
  for (x = 0; x < n; x++) result += vec1[x] * vec2[x];
  return result;
}
int
esl_vec_IDot(int *vec1, int *vec2, int n)
{
  int result = 0;
  int x;
  for (x = 0; x < n; x++) result += vec1[x] * vec2[x];
  return result;
}



/* Function:  esl_vec_DMax()
 *
 * Purpose:   Returns the maximum value of the <n> values
 *            in <vec>.
 *            
 *            <esl_vec_FMax()> and <esl_vec_IMax()> do the same,
 *            for float and integer vectors.
 */
double
esl_vec_DMax(double *vec, int n)
{
  int i;
  double best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] > best) best = vec[i];
  return best;
}
float
esl_vec_FMax(float *vec, int n)
{
  int   i;
  float best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] > best) best = vec[i];
  return best;
}
int
esl_vec_IMax(int *vec, int n)
{
  int   i;
  int   best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] > best) best = vec[i];
  return best;
}


/* Function:  esl_vec_DMin()
 *
 * Purpose:   Returns the minimum value of the <n> values
 *            in <vec>.
 *            
 *            <esl_vec_FMin()> and <esl_vec_IMin()> do the same,
 *            for float and integer vectors.
 */
double
esl_vec_DMin(double *vec, int n)
{
  int i;
  double best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] < best) best = vec[i];
  return best;
}
float
esl_vec_FMin(float *vec, int n)
{
  int   i;
  float best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] < best) best = vec[i];
  return best;
}
int
esl_vec_IMin(int *vec, int n)
{
  int   i;
  int   best;

  best = vec[0];
  for (i = 1; i < n; i++)
    if (vec[i] < best) best = vec[i];
  return best;
}


/* Function:  esl_vec_DArgMax()
 *
 * Purpose:   Returns the index of the maximum value in the <n> values
 *            in <vec>.
 *            
 *            <esl_vec_FArgMax()> and <esl_vec_IArgMax()> do the same,
 *            for float and integer vectors.
 */
int
esl_vec_DArgMax(double *vec, int n)
{
  int i;
  int best = 0;

  for (i = 1; i < n; i++)
    if (vec[i] > vec[best]) best = i;
  return best;
}
int
esl_vec_FArgMax(float *vec, int n)
{
  int i;
  int best = 0;

  for (i = 1; i < n; i++)
    if (vec[i] > vec[best]) best = i;
  return best;
}
int
esl_vec_IArgMax(int *vec, int n)
{
  int i;
  int best = 0;

  for (i = 1; i < n; i++)
    if (vec[i] > vec[best]) best = i;
  return best;
}


/* Function:  esl_vec_DArgMin()
 *
 * Purpose:   Returns the index of the minimum value in the <n> values
 *            in <vec>.
 *            
 *            <esl_vec_FArgMin()> and <esl_vec_IArgMin()> do the same,
 *            for float and integer vectors.
 */
int
esl_vec_DArgMin(double *vec, int n)
{
  int i;
  int best = 0;
  for (i = 1; i < n; i++)
    if (vec[i] < vec[best]) best = i;
  return best;
}
int
esl_vec_FArgMin(float *vec, int n)
{
  int   i;
  int   best = 0;

  for (i = 1; i < n; i++)
    if (vec[i] < vec[best]) best = i;
  return best;
}
int
esl_vec_IArgMin(int *vec, int n)
{
  int   i;
  int   best = 0;

  for (i = 1; i < n; i++)
    if (vec[i] < vec[best]) best = i;
  return best;
}


/* some static functions to pass to qsort() that the 
 * upcoming Sort() functions will call
 */
static int
qsort_DIncreasing(const void *xp1, const void *xp2)
{
  double x1 = * (double *) xp1;
  double x2 = * (double *) xp2; 
  if (x1 < x2) return -1;
  if (x1 > x2) return 1;
  return 0;
}
static int
qsort_FIncreasing(const void *xp1, const void *xp2)
{
  float x1 = * (float *) xp1;
  float x2 = * (float *) xp2; 
  if (x1 < x2) return -1;
  if (x1 > x2) return 1;
  return 0;
}
static int
qsort_IIncreasing(const void *xp1, const void *xp2)
{
  int x1 = * (int *) xp1;
  int x2 = * (int *) xp2; 
  if (x1 < x2) return -1;
  if (x1 > x2) return 1;
  return 0;
}
static int
qsort_DDecreasing(const void *xp1, const void *xp2)
{
  double x1 = * (double *) xp1;
  double x2 = * (double *) xp2; 
  if (x1 > x2) return -1;
  if (x1 < x2) return 1;
  return 0;
}
static int
qsort_FDecreasing(const void *xp1, const void *xp2)
{
  float x1 = * (float *) xp1;
  float x2 = * (float *) xp2; 
  if (x1 > x2) return -1;
  if (x1 < x2) return 1;
  return 0;
}
static int
qsort_IDecreasing(const void *xp1, const void *xp2)
{
  int x1 = * (int *) xp1;
  int x2 = * (int *) xp2; 
  if (x1 > x2) return -1;
  if (x1 < x2) return 1;
  return 0;
}

/* Function:  esl_vec_DSortIncreasing()
 * Incept:    SRE, Wed Aug 17 10:44:31 2005 [St. Louis]
 *
 * Purpose:   Sorts <vec> in place, from smallest to largest value.
 *            (That is, <vec[0]> is the minimum and <vec[n-1]> is
 *            the maximum.)
 *            
 *            <esl_vec_FSortIncreasing()> and <esl_vec_ISortIncreasing()>
 *            do the same, for float and integer vectors.
 */
void
esl_vec_DSortIncreasing(double *vec, int n)
{
  qsort((void *) vec, n, sizeof(double), qsort_DIncreasing);
}
void
esl_vec_FSortIncreasing(float *vec, int n)
{
  qsort((void *) vec, n, sizeof(float), qsort_FIncreasing);
}
void
esl_vec_ISortIncreasing(int *vec, int n)
{
  qsort((void *) vec, n, sizeof(int), qsort_IIncreasing);
}

/* Function:  esl_vec_DSortDecreasing()
 * Incept:    SRE, Wed Aug 17 10:44:31 2005 [St. Louis]
 *
 * Purpose:   Sorts <vec> in place, from largest to smallest value.
 *            (That is, <vec[0]> is the maximum and <vec[n-1]> is
 *            the minimum.)
 *            
 *            <esl_vec_FSortDecreasing()> and <esl_vec_ISortDecreasing()>
 *            do the same, for float and integer vectors.
 */
void
esl_vec_DSortDecreasing(double *vec, int n)
{
  qsort((void *) vec, n, sizeof(double), qsort_DDecreasing);
}
void
esl_vec_FSortDecreasing(float *vec, int n)
{
  qsort((void *) vec, n, sizeof(float), qsort_FDecreasing);
}
void
esl_vec_ISortDecreasing(int *vec, int n)
{
  qsort((void *) vec, n, sizeof(int), qsort_IDecreasing);
}


/* Function:  esl_vec_DDump()
 * Incept:    ER, Thu Jul 21 12:54:56 CDT 2005 [St. Louis]
 *
 * Purpose:   Given a vector, dump it to stream <ofp>.
 * 
 *            If <label> is non-NULL, they represent
 *            single-character labels to put on the vector. 
 *            (For example, these might be a sequence alphabet).
 *            Numbers 1..n is used if <label> is NULL.
 *
 * Args:      ofp   -  output file pointer; stdout, for example.
 *            v     -  vector to dump.
 *            label -  optional: NULL, or character labels
 *
 * Returns:   <eslOK> on success.
 */
int
esl_vec_DDump(FILE *ofp, double *v, int n, char *label)
{
  int a;

  fprintf(ofp, "     ");
  if (label != NULL) 
    for (a = 0; a < n; a++) fprintf(ofp, "       %c ", label[a]);
  else
    for (a = 0; a < n; a++) fprintf(ofp, "%8d ", a+1);
  fprintf(ofp, "\n");
  
  fprintf(ofp, "      ");
  for (a = 0; a < n; a++) fprintf(ofp, "%8.8f ", v[a]);
  fprintf(ofp, "\n");

  return eslOK;
}
int
esl_vec_FDump(FILE *ofp, float *v, int n, char *label)
{
  int a;

  fprintf(ofp, "     ");
  if (label != NULL) 
    for (a = 0; a < n; a++) fprintf(ofp, "       %c ", label[a]);
  else
    for (a = 0; a < n; a++) fprintf(ofp, "%8d ", a+1);
  fprintf(ofp, "\n");
  
  fprintf(ofp, "      ");
  for (a = 0; a < n; a++) fprintf(ofp, "%8.4f ", v[a]);
  fprintf(ofp, "\n");

  return eslOK;
}
int
esl_vec_IDump(FILE *ofp, int *v, int n, char *label)
{
  int a;

  fprintf(ofp, "     ");
  if (label != NULL) 
    for (a = 0; a < n; a++) fprintf(ofp, "       %c ", label[a]);
  else
    for (a = 0; a < n; a++) fprintf(ofp, "%8d ", a+1);
  fprintf(ofp, "\n");
  
  fprintf(ofp, "      ");
  for (a = 0; a < n; a++) fprintf(ofp, "%8d ", v[a]);
  fprintf(ofp, "\n");

  return eslOK;
}

/* Function:  esl_vec_D2F()
 * Incept:    SRE, Thu Mar 30 09:04:17 2006 [St. Louis]
 *
 * Purpose:   Copy a double vector <src> to a float vector <dst>. Caller
 *            provides space in the float vector that is at
 *            least <n>.
 *            
 *            Similarly, <esl_vec_F2D()> converts float to double; 
 *            <esl_vec_I2D()> converts integer to double; 
 *            <esl_vec_I2F()> converts integer to float.
 */
void
esl_vec_D2F(double *src, int n, float *dst)
{
  int i;
  for (i = 0; i < n; i++) dst[i] = src[i];
}
void
esl_vec_F2D(float *src, int n, double *dst)
{
  int i;
  for (i = 0; i < n; i++) dst[i] = src[i];
}
void
esl_vec_I2F(int *src, int n, float *dst)
{
  int i;
  for (i = 0; i < n; i++) dst[i] = src[i];
}
void
esl_vec_I2D(int *src, int n, double *dst)
{
  int i;
  for (i = 0; i < n; i++) dst[i] = src[i];
}




/* Function:  esl_vec_DNorm()
 *
 * Purpose:   Normalizes a probability vector <vec>,
 *            such that $\sum_{i=1}{n} \mathrm{vec}_i = 1.0$.
 *            
 *            <esl_vec_FNorm()> does the same, for a probability vector
 *            of floats.
 */
void
esl_vec_DNorm(double *vec, int n)
{
  int    x;
  double sum;

  sum = esl_vec_DSum(vec, n);
  if (sum != 0.0) for (x = 0; x < n; x++) vec[x] /= sum;
  else            for (x = 0; x < n; x++) vec[x] = 1. / (double) n;
}
void
esl_vec_FNorm(float *vec, int n)
{
  int    x;
  float  sum;

  sum = esl_vec_FSum(vec, n);
  if (sum != 0.0) for (x = 0; x < n; x++) vec[x] /= sum;
  else            for (x = 0; x < n; x++) vec[x] = 1. / (float) n;
}


/* Function:  esl_vec_DLog()
 *
 * Purpose:   Converts a probability vector <vec> to a log
 *            probability vector: takes the log of each of the <n> 
 *            values in the vector.
 *
 *            <esl_vec_FLog()> does the same, for a probability vector
 *            of floats.
 */
void
esl_vec_DLog(double *vec, int n)
{
  int x;
  for (x = 0; x < n; x++) 
    if (vec[x] > 0.) vec[x] = log(vec[x]);
    else vec[x] = -DBL_MAX;
}
void
esl_vec_FLog(float *vec, int n)
{
  int x;
  for (x = 0; x < n; x++) 
    if (vec[x] > 0.) vec[x] = log(vec[x]);
    else vec[x] = -FLT_MAX;
}


/* Function:  esl_vec_DEntropy()
 *
 * Purpose:   Returns the Shannon entropy of a probability vector <p>,
 *            in bits ($\log_2$).
 *
 *            <esl_vec_FEntropy()> does the same, for a probability vector
 *            of floats.
 */
double
esl_vec_DEntropy(double *p, int n)
{
  int    i;
  double entropy;
 
  entropy = 0.;
  for(i = 0; i < n; i++)
    if (p[i] > 0.) entropy += p[i] * log(p[i]);
  return(-1.44269504 * entropy); /* converts to bits */
}
float
esl_vec_FEntropy(float *p, int n)
{
  int    i;
  float  entropy;

  entropy = 0.;
  for(i = 0; i < n; i++)
    if (p[i] > 0.) entropy += p[i] * log(p[i]);
  return(-1.44269504 * entropy); /* converts to bits */
}


/* Function:  esl_vec_DExp()
 *
 * Purpose:   Converts a log probability vector <vec> back to a 
 *            probability vector: exponentiates each of the <n> 
 *            values in the vector.
 *
 *            <esl_vec_FExp()> does the same, for a log probability vector
 *            of floats.
 */
void
esl_vec_DExp(double *vec, int n)
{
  int x;
  for (x = 0; x < n; x++) vec[x] = exp(vec[x]);
  esl_vec_DNorm(vec, n);
}
void
esl_vec_FExp(float *vec, int n)
{
  int x;
  for (x = 0; x < n; x++) vec[x] = exp(vec[x]);
  esl_vec_FNorm(vec, n);
}

/* Function:  esl_vec_DLogSum()
 *
 * Purpose:   <vec> is a log probability vector; return the log of the scalar sum
 *            of the probabilities in <vec>. That is, the <n> elements in <vec>
 *            are log probabilities, but the summation is done in probability
 *            space, by exponentiating each of the <n> values in the vector,
 *            summing, and returning the log of the sum. 
 *            
 *            That is: return $\log \sum_i e^{v_i}$.
 *
 *            The trick is to do this without numerical underflow or overflow.
 *
 *            <esl_vec_FLogSum()> does the same, for a log probability vector
 *            of floats.
 */
double
esl_vec_DLogSum(double *vec, int n)
{
  int x;
  double max, sum;
  
  max = esl_vec_DMax(vec, n);
  if (max == eslINFINITY) return eslINFINITY; /* avoid inf-inf below! */
  sum = 0.0;
  for (x = 0; x < n; x++)
    if (vec[x] > max - 50.)
      sum += exp(vec[x] - max);
  sum = log(sum) + max;
  return sum;
}
float
esl_vec_FLogSum(float *vec, int n)
{
  int x;
  float max, sum;
  
  max = esl_vec_FMax(vec, n);
  sum = 0.0;
  for (x = 0; x < n; x++)
    if (vec[x] > max - 50.)
      sum += exp(vec[x] - max);
  sum = log(sum) + max;
  return sum;
}


/* Function:  esl_vec_DLogNorm()
 * Incept:    SRE, Thu Apr  7 17:45:39 2005 [St. Louis]
 *
 * Purpose:   Given an unnormalized log probability vector <vec>   
 *            of length <n>, normalize it and make it a 
 *            probability vector. 
 *            
 *            <esl_vec_FLogNorm()> does the same, but for a vector
 *            of floats instead of doubles.
 *
 * Returns:   (void); <vec> is changed in place.
 */
void
esl_vec_DLogNorm(double *vec, int n)
{
  double denom;
  
  denom = esl_vec_DLogSum(vec, n);
  esl_vec_DIncrement(vec, n, -1.*denom);
  esl_vec_DExp(vec, n);
}
void
esl_vec_FLogNorm(float *vec, int n)
{
  float denom;
  
  denom = esl_vec_FLogSum(vec, n);
  esl_vec_FIncrement(vec, n, -1.*denom);
  esl_vec_FExp(vec, n);
}

/* Function:  esl_vec_DValidate()
 * Incept:    ER, Tue Dec  5 09:38:54 EST 2006 [janelia]
 *
 * Purpose:   Validate a probability vector <vec> of length <n>.
 *            Each element has to be between 0 and 1, and
 *            the sum of all elements has to be 1.
 *
 * Args:      v      - p vector to validate.
 *            n      - dimensionality of v
 *            tol    - convergence criterion applied to sum of v
 *            errbuf - NULL, or a failure message buffer allocated
 *                     for at least <eslERRBUFSIZE> chars. 
 *
 * Returns:   <eslOK> on success, or <eslFAIL> on validation failure.
 *            Upon failure, if caller provided a non-<NULL> <errbuf>,
 *            an informative message is left there.
 */
int
esl_vec_DValidate(double *vec, int n, double tol, char *errbuf)
{
  int    status;
  int    x;
  double sum = 0.;

  if (errbuf) *errbuf = 0;
  if (n == 0) return eslOK;

  for (x = 0; x < n; x++) {
    if (vec[x] < 0.0 || vec[x] > 1.0)
      ESL_XFAIL(eslFAIL, errbuf, "value %d is not a probability between 0..1", x);
    sum += vec[x];
  }

  if (fabs(sum - 1.0) > tol) 
    ESL_XFAIL(eslFAIL, errbuf, "vector does not sum to 1.0");
  return eslOK;

 ERROR:
  return status;
}
int
esl_vec_FValidate(float *vec, int n, float tol, char *errbuf)
{
  int   status;
  int   x;
  float sum = 0.;

  if (errbuf) *errbuf = 0;
  if (n == 0) return eslOK;

  for (x = 0; x < n; x++) {
    if (vec[x] < 0.0 || vec[x] > 1.0)
      ESL_XFAIL(eslFAIL, errbuf, "value %d is not a probability between 0..1", x);
    sum += vec[x];
  }

  if (fabs(sum - 1.0) > tol) 
    ESL_XFAIL(eslFAIL, errbuf, "vector does not sum to 1.0");
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_vec_DLogValidate()
 * Incept:    ER,  Tue Dec  5 09:46:51 EST 2006 [janelia]
 *
 * Purpose:   Validate a log probability vector <vec> of length <n>.
 *            The exp of each element has to be between 0 and 1, and
 *            the sum of all elements has to be 1.
 *
 * Args:      v      - log p vector to validate.
 *            n      - dimensionality of v
 *            tol    - convergence criterion applied to sum of exp v
 *            errbuf - NULL, or a failure message buffer allocated
 *                     for at least p7_ERRBUFSIZE chars. 
 *
 * Returns:   <eslOK> on success, or <eslFAIL> on failure; upon failure,
 *            if caller provided a non-<NULL> <errbuf>, an informative
 *            message is left there.
 *            
 * Throws:    <eslEMEM> on allocation failure.           
 */
int
esl_vec_DLogValidate(double *vec, int n, double tol, char *errbuf)
{
  int     status;
  double *expvec = NULL;

  if (errbuf) *errbuf = 0;
  if (n == 0) return eslOK;

  ESL_ALLOC(expvec, sizeof(double)*n);
  esl_vec_DCopy(vec, n, expvec);
  esl_vec_DExp(expvec, n); 
  if ((status = esl_vec_DValidate(expvec, n, tol, errbuf)) != eslOK) goto ERROR;
  free(expvec);
  return eslOK;

 ERROR:
  if (expvec != NULL) free(expvec);
  return status;
}
int
esl_vec_FLogValidate(float *vec, int n, float tol, char *errbuf)
{
  int     status;
  float  *expvec = NULL;

  if (errbuf) *errbuf = 0;
  if (n == 0) return eslOK;

  ESL_ALLOC(expvec, sizeof(float)*n);
  esl_vec_FCopy(vec, n, expvec);
  esl_vec_FExp(expvec, n); 
  if ((status = esl_vec_FValidate(expvec, n, tol, errbuf)) != eslOK) goto ERROR;
  free(expvec);
  return eslOK;

 ERROR:
  if (expvec != NULL) free(expvec);
  return eslOK;
}

/*****************************************************************  
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/
