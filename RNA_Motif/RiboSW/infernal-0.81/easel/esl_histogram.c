/* esl_histogram.c
 * Collecting and displaying histograms.
 *
 * SRE, Fri Jul  1 13:21:45 2005 [St. Louis]
 * SVN $Id: esl_histogram.c 155 2007-03-04 17:32:24Z eddys $
 */
#include <esl_config.h>

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <easel.h>
#include <esl_stats.h>
#include <esl_histogram.h>
#include <esl_vectorops.h>

static int esl_histogram_sort(ESL_HISTOGRAM *h);


/*****************************************************************
 * Creating/destroying histograms and collecting data.
 *****************************************************************/

/* Function:  esl_histogram_Create()
 * Incept:    SRE, Fri Jul  1 13:40:26 2005 [St. Louis]
 *
 * Purpose:   Creates and returns a new histogram object, initially
 *            allocated to count scores $>$ <xmin> and $<=$ <xmax> into
 *            bins of width <w>. Thus, a total of <xmax>-<xmin>/<w> bins
 *            are initially created. 
 *            
 *            The lower bound <xmin> and the width <w> permanently
 *            determine the offset and width of the binning, but not
 *            the range.  For example, <esl_histogram_Create(-100,
 *            100, 0.5)> would init the object to collect scores into
 *            400 bins $[-100< x \leq -99.5],[-99.5 < x \leq
 *            -99.0]...[99.5 <x \leq 100.0]$.  Aside from this, the
 *            range specified by the bounds <xmin> and <xmax> only
 *            needs to be an initial guess. The histogram object will
 *            reallocate itself dynamically as needed to accommodate
 *            scores that exceed current bounds.
 *
 *            You can be sloppy about <xmax>; it does not have to
 *            exactly match a bin upper bound. The initial allocation
 *            is for all full-width bins with upper bounds $\leq
 *            xmax$.
 *
 *            <esl_histogram_Create()> creates a simplified histogram
 *            object that collates only the "display" histogram. For
 *            a more complex object that also keeps the raw data samples,
 *            better suited for fitting distributions and goodness-of-fit
 *            testing, use <esl_histogram_CreateFull()>.
 *  
 * Args:      xmin - caller guesses that minimum score will be > xmin
 *            xmax - caller guesses that max score will be <= xmax
 *            w    - size of bins (1.0, for example)
 *            
 * Returns:   ptr to new <ESL_HISTOGRAM> object, which caller is responsible
 *            for free'ing with <esl_histogram_Destroy()>.
 *
 * Throws:    <NULL> on allocation failure.
 */
ESL_HISTOGRAM *
esl_histogram_Create(double xmin, double xmax, double w)
{
  ESL_HISTOGRAM *h = NULL;
  int status;

  ESL_ALLOC(h, sizeof(ESL_HISTOGRAM));

  h->xmin      =  DBL_MAX;	/* xmin/xmax are the observed min/max */
  h->xmax      = -DBL_MAX;
  h->n         = 0;
  h->obs       = NULL;		/* briefly... */
  h->bmin      = xmin;		/* bmin/bmax are the allocated bounds */
  h->bmax      = xmax;
  h->nb        = (int)((xmax-xmin)/w);
  h->imin      = h->nb;
  h->imax      = -1;
  h->w         = w;

  h->x         = NULL;
  h->nalloc    = 0;

  h->phi       = 0.;
  h->cmin      = h->imin;	/* sentinel: no observed data yet */
  h->z         = 0;
  h->Nc        = 0;
  h->No        = 0;

  h->expect    = NULL;		/* 'til a Set*() call */
  h->emin      = -1;            /* sentinel: no expected counts yet */
  h->tailbase  = 0.;		/* unused unless is_tailfit TRUE */
  h->tailmass  = 1.0;		/* <= 1.0 if is_tailfit TRUE */

  h->is_full       = FALSE;
  h->is_done       = FALSE;
  h->is_sorted     = FALSE;
  h->is_tailfit    = FALSE;
  h->is_rounded    = FALSE;
  h->dataset_is    = COMPLETE;

  ESL_ALLOC(h->obs, sizeof(int) * h->nb);
  esl_vec_ISet(h->obs, h->nb, 0);
  return h;

 ERROR:
  esl_histogram_Destroy(h);
  return NULL;
}

/* Function:  esl_histogram_CreateFull()
 * Incept:    SRE, Tue Jul 26 13:19:27 2005 [St. Louis]
 *
 * Purpose:   Alternative form of <esl_histogram_Create()> that 
 *            creates a more complex histogram that will contain not just the 
 *            display histogram, but also keeps track of all
 *            the raw sample values. Having a complete vector of raw
 *            samples improves distribution-fitting and goodness-of-fit 
 *            tests, but will consume more memory. 
 */
ESL_HISTOGRAM *
esl_histogram_CreateFull(double xmin, double xmax, double w)
{
  int status;
  ESL_HISTOGRAM *h = esl_histogram_Create(xmin, xmax, w);
  if (h == NULL) return NULL;

  h->n      = 0;		/* make sure */
  h->nalloc = 128;		/* arbitrary initial allocation size */
  ESL_ALLOC(h->x, sizeof(double) * h->nalloc);
  h->is_full = TRUE;
  return h;

 ERROR:
  esl_histogram_Destroy(h);
  return NULL;
}


/* Function:  esl_histogram_Destroy()
 * Incept:    SRE, Sat Jul  2 19:41:17 2005 [St. Louis]
 *
 * Purpose:   Frees an <ESL_HISTOGRAM> object <h>.
 */
void
esl_histogram_Destroy(ESL_HISTOGRAM *h)
{
  if (h ==  NULL) return;
  if (h->x      != NULL) free(h->x);
  if (h->obs    != NULL) free(h->obs); 
  if (h->expect != NULL) free(h->expect);
  free(h);
  return;
}

/* Function:  esl_histogram_Add()
 * Incept:    SRE, Sat Jul  2 19:41:45 2005 [St. Louis]
 *
 * Purpose:   Adds score <x> to a histogram <h>.
 *           
 *            The histogram will be automatically reallocated as
 *            needed if the score is smaller or larger than the
 *            current allocated bounds.  
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 *            Throws <eslEINVAL> for cases where something has been done
 *            to the histogram that requires it to be 'finished', and
 *            adding more data is prohibited; for example, 
 *            if tail or censoring information has already been set.
 *            On either failure, initial state of <h> is preserved.
 */
int
esl_histogram_Add(ESL_HISTOGRAM *h, double x)
{
  int   status;
  void *tmp;
  int b;			/* what bin we're in                       */
  int nnew;			/* # of new bins created by a reallocation */

  /* Censoring info must only be set on a finished histogram;
   * don't allow caller to add data after configuration has been declared
   */
  if (h->is_done)
    ESL_EXCEPTION(eslEINVAL, "can't add more data to this histogram");

  /* If we're a full histogram, check whether we need to reallocate
   * the full data vector.
   */
  if (h->is_full && h->nalloc == h->n) 
    {
      ESL_RALLOC(h->x, tmp, sizeof(double) * h->nalloc * 2);
      h->nalloc *= 2;
    }

  /* Which bin will we want to put x into?
   */
  b = esl_histogram_Score2Bin(h,x);

  /* Make sure we have that bin. Realloc as needed.
   * If that reallocation succeeds, we can no longer fail;
   * so we can change the state of h.
   */
  if (b < 0)    /* Reallocate below? */
    {				
      nnew = -b*2;	/* overallocate by 2x */
      ESL_RALLOC(h->obs, tmp, sizeof(int) * (nnew+ h->nb));
      
      memmove(h->obs+nnew, h->obs, sizeof(int) * h->nb);
      h->nb    += nnew;
      b        += nnew;
      h->bmin  -= nnew*h->w;
      h->imin  += nnew;
      h->cmin  += nnew;
      if (h->imax > -1) h->imax += nnew;
      esl_vec_ISet(h->obs, nnew, 0);
    }
  else if (b >= h->nb)  /* Reallocate above? */
    {
      nnew = (b-h->nb+1) * 2; /* 2x overalloc */
      ESL_RALLOC(h->obs, tmp, sizeof(int) * (nnew+ h->nb));
      esl_vec_ISet(h->obs+h->nb, nnew, 0);
      if (h->imin == h->nb) { /* boundary condition of no data yet*/
	h->imin+=nnew; 
	h->cmin+=nnew;
      }
      h->bmax  += nnew*h->w;
      h->nb    += nnew;
    }

  /* If we're a full histogram, then we keep the raw x value,
   * reallocating as needed.
   */
  if (h->is_full)  h->x[h->n] = x;
  h->is_sorted = FALSE;		/* not any more! */

  /* Bump the bin counter, and all the data sample counters.
   */
  h->obs[b]++;
  h->n++;
  h->Nc++;
  h->No++;

  if (b > h->imax) h->imax = b;
  if (b < h->imin) { h->imin = b; h->cmin = b; }
  if (x > h->xmax) h->xmax = x;
  if (x < h->xmin) h->xmin = x;
  return eslOK;

 ERROR:
  return status;
}
  

/* esl_histogram_sort()
 * Incept:    SRE, Thu Aug 18 10:45:46 2005 [St. Louis]
 *
 * Purpose:   Sort the raw scores in a full histogram, from smallest to
 *            largest. Has no effect on a normal histogram, or on a full
 *            histogram that is already sorted.
 *
 * Returns:   <eslOK> on success.
 *            Upon return, <h->x[h->n-1]> is the high score, <h->x[0]> is the 
 *            low score. 
 */
int
esl_histogram_sort(ESL_HISTOGRAM *h)
{
  if (h->is_sorted) return eslOK; /* already sorted, don't do anything */
  if (! h->is_full) return eslOK; /* nothing to sort */
  
  esl_vec_DSortIncreasing(h->x, h->n);
  h->is_sorted = TRUE;
  return eslOK;
}

/*****************************************************************
 * Declarations about the binned data before parameter fitting
 *****************************************************************/ 

/* Function:  esl_histogram_DeclareCensoring()
 * Incept:    SRE, Tue Aug 23 10:00:14 2005 [St. Louis]
 *
 * Purpose:   Declare that the dataset collected in <h> is known to be a
 *            censored distribution, where <z> samples were unobserved because
 *            they had values $\leq$ some threshold <phi> ($\phi$).
 *            
 *            No more data can be added to the histogram with <_Add()>
 *            after censoring information has been set.
 *            
 *            This function is for "true" censored datasets, where
 *            the histogram truly contains no observed points
 *            $x \leq \phi$, and the number that were censored is known
 *            to be <z>. 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if you try to set <phi> to a value that is
 *            greater than the minimum <x> stored in the histogram.
 */
int
esl_histogram_DeclareCensoring(ESL_HISTOGRAM *h, int z, double phi)
{
  if (phi > h->xmin) ESL_EXCEPTION(eslEINVAL, "no uncensored x can be <= phi");

  h->phi         = phi;
  h->cmin        = h->imin;
  h->z           = z;
  h->Nc          = h->n + z;
  h->No          = h->n;
  h->dataset_is  = TRUE_CENSORED;
  h->is_done     = TRUE;
  return eslOK;
}

/* Function:  esl_histogram_DeclareRounding()
 * Incept:    SRE, Tue Jan 31 13:52:10 2006 [St. Louis]
 *
 * Purpose:   Declare that the data sample values in the histogram <h>
 *            are rounded off. Ideally, your bins in <h> should exactly 
 *            match the rounding procedure. This raises a flag that
 *            binned parameter fitting routines will use when they set
 *            an origin, using the lower bound of the bin instead of
 *            the lowest raw value in the bin.
 */
int
esl_histogram_DeclareRounding(ESL_HISTOGRAM *h)
{
  h->is_rounded = TRUE;
  return eslOK;
}


/* Function:  esl_histogram_SetTail()
 * Incept:    SRE, Tue Aug 23 09:01:10 2005 [St. Louis]
 *
 * Purpose:   Suggest a threshold <phi> to split a histogram <h>
 *            into "unobserved" data (values $\leq \phi$) and "observed" 
 *            data (values $> \phi$). 
 *
 *            The suggested <phi> is revised downwards to a $\phi$ at the next 
 *            bin lower bound, because operations on binned data in <h>
 *            need to know unambiguously whether all the data in a given bin
 *            will be counted as observed or unobserved. 
 *
 *            The probability mass that is in the resulting right tail
 *            is optionally returned in <ret_newmass>. You need to know
 *            this number if you're fitting a distribution solely to the
 *            tail (an exponential tail, for example).
 *
 *            Any data point $x_i \leq \phi$ is then considered to be
 *            in a censored (unobserved) region for purposes of parameter
 *            fitting, calculating expected binned counts,
 *            and binned goodness-of-fit tests. 
 *            
 *            No more data can be added to the histogram after
 *            censoring information has been set.
 *            
 *            This function defines a "virtual" left-censoring: the
 *            histogram actually contains complete data, but appropriate
 *            flags are set to demarcate the "observed" data in the right
 *            tail.
 *            
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_SetTail(ESL_HISTOGRAM *h, double phi, double *ret_newmass)
{
  int b;

  /* Usually, put true phi at the next bin lower bound, but
   * watch for a special case where phi is already exactly equal to a 
   * bin upper bound.
   */
  h->cmin = esl_histogram_Score2Bin(h,phi);
  if (phi == esl_histogram_Bin2UBound(h,h->cmin)) h->phi = phi;
  else   h->phi  = esl_histogram_Bin2LBound(h, h->cmin);

  h->z    = 0;
  for (b = h->imin; b < h->cmin; b++)
    h->z += h->obs[b];
  h->Nc         = h->n;		/* (redundant) */
  h->No         = h->n - h->z;
  h->dataset_is = VIRTUAL_CENSORED;
  h->is_done    = TRUE;
  if (ret_newmass != NULL) *ret_newmass = (double) h->No / (double) h->Nc;
  return eslOK;
}

/* Function:  esl_histogram_SetTailByMass()
 * Incept:    SRE, Tue Aug 23 08:10:39 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h> (with or without raw data samples),
 *            find a cutoff score that at least fraction <pmass> of the samples
 *            exceed. This threshold is stored internally in the histogram
 *            as <h->phi>. The number of "virtually censored" samples (to the 
 *            left, with scores $\leq \phi$) is stored internally in <h->z>.
 *            
 *            The identified cutoff score must be a lower bound for some bin
 *            (bins can't be partially censored). The censored mass
 *            will thus usually be a bit greater than <pmass>, as the
 *            routine will find the highest satisfactory <h->phi>. The
 *            narrower the bin widths, the more accurately the routine
 *            will be able to satisfy the requested <frac>. The actual
 *            probability mass in the right tail is optionally returned
 *            in <ret_newmass>. You need to know this number if you're 
 *            fitting a distribution solely to the tail (an exponential tail,
 *            for example). It is safe for <ret_newmass> to point at 
 *            <pmass>, in which case the suggested <pmass> will be overwritten
 *            with the actual mass upon return.
 *
 *            This function defines that the binned data will be
 *            fitted either as a tail, or as a (virtually) left-censored dataset.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_SetTailByMass(ESL_HISTOGRAM *h, double pmass, double *ret_newmass)
{
  int b;
  int sum = 0;
	    
  for (b = h->imax; b >= h->imin; b--)
    {
      sum += h->obs[b];
      if (sum >= (pmass * (double)h->n)) break;
    }

  h->phi         = esl_histogram_Bin2LBound(h,b);
  h->z           = h->n - sum;
  h->cmin        = b;
  h->Nc          = h->n;	/* (redundant) */
  h->No          = h->n - h->z;
  h->dataset_is  = VIRTUAL_CENSORED;
  h->is_done     = TRUE;
  if (ret_newmass != NULL) *ret_newmass = (double) h->No / (double) h->Nc;
  return eslOK;
}



/*****************************************************************
 * Routines for accessing data samples in a full histogram.
 *****************************************************************/

/* Function:  esl_histogram_GetRank()
 * Incept:    SRE, Thu Jul 28 08:39:52 2005 [St. Louis]
 *
 * Purpose:   Retrieve the <rank>'th highest score from a 
 *            full histogram <h>. <rank> is <1..n>, for
 *            <n> total samples in the histogram; return it through
 *            <ret_x>.
 *            
 *            If the raw scores aren't sorted, they are sorted
 *            first (an $N \log N$ operation).
 *            
 *            This can be called at any time, even during data
 *            collection, to see the current <rank>'th highest score.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if the histogram is display-only,
 *            or if <rank> isn't in the range 1..n.
 */
int
esl_histogram_GetRank(ESL_HISTOGRAM *h, int rank, double *ret_x)
{
  if (! h->is_full) 
    ESL_EXCEPTION(eslEINVAL, 
	      "esl_histogram_GetRank() needs a full histogram");
  if (rank > h->n)
    ESL_EXCEPTION(eslEINVAL, 
	      "no such rank: not that many scores in the histogram");
  if (rank < 1)
    ESL_EXCEPTION(eslEINVAL, "histogram rank must be a value from 1..n");

  esl_histogram_sort(h);	/* make sure */
  *ret_x = h->x[h->n - rank + 1];
  return eslOK;
}

/* Function:  esl_histogram_GetData()
 * Incept:    SRE, Fri Jan 27 07:57:21 2006 [St. Louis]
 *
 * Purpose:   Retrieve the raw data values from the histogram <h>.
 *            Return them in the vector <ret_x>, and the number
 *            of values in <ret_n>. The values are indexed <[0..n-1]>,
 *            from smallest to largest (<x[n-1]> is the high score).
 *            
 *            <ret_x> is a pointer to internal memory in the histogram <h>.
 *            The histogram <h> is still responsible for that storage;
 *            its memory will be free'd when you call
 *            <esl_histogram_Destroy()>.
 *            
 *            You can only call this after you have finished collecting
 *            all the data. Subsequent calls to <esl_histogram_Add()>
 *            will fail.
 *            
 * Internal note:
 *            The prohibition against adding more data (by raising
 *            the h->is_done flag) is because we're passing a pointer
 *            to internal data storage back to the caller. Subsequent
 *            calls to Add() will modify that memory -- in the worst case,
 *            if Add() has to reallocate that storage, completely invalidating
 *            the pointer that the caller has a copy of. We want to make
 *            sure that the <ret_x> pointer stays valid.
 *            
 * Args:      h     - histogram to retrieve data values from
 *            ret_x - RETURN: pointer to the data samples, [0..n-1] 
 *            ret_n - RETURN: number of data samples
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if the histogram <h> is not a full histogram.
 */
int
esl_histogram_GetData(ESL_HISTOGRAM *h, double **ret_x, int *ret_n)
{
  if (! h->is_full) ESL_EXCEPTION(eslEINVAL, "not a full histogram");
  esl_histogram_sort(h);

  *ret_x = h->x;
  *ret_n = h->n;

  h->is_done = TRUE;
  return eslOK;
}


/* Function:  esl_histogram_GetTail()
 * Incept:    SRE, Fri Jan 27 07:56:38 2006 [St. Louis]
 *
 * Purpose:   Given a full histogram <h>, retrieve all data values 
 *            above the threshold <phi> in the right (high scoring) 
 *            tail, as a ptr <ret_x> to an array of <ret_n> values 
 *            indexed <[0..n-1]> from lowest to highest score. 
 *            Optionally, it also returns the number of values in 
 *            rest of the histogram in <ret_z>;
 *            this number is useful if you are going to fit
 *            the tail as a left-censored distribution.
 *            
 *            The test is strictly greater than <phi>, not greater
 *            than or equal to.
 *            
 *            <ret_x> is a pointer to internal memory in the histogram <h>.
 *            The histogram <h> is still responsible for that storage;
 *            its memory will be free'd when you call 
 *            <esl_histogram_Destroy()>.
 *            
 *            You can only call this after you have finished collecting
 *            all the data. Subsequent calls to <esl_histogram_Add()>
 *            will fail.             
 *            
 * Args:      h     - histogram to retrieve the tail from
 *            phi   - threshold: tail is all scores > phi
 *            ret_x - optRETURN: ptr to vector of data values [0..n-1]
 *            ret_n - optRETURN: number of data values in tail
 *            ret_z - optRETURN: number of data values not in tail.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if the histogram is not a full histogram.
 */
int
esl_histogram_GetTail(ESL_HISTOGRAM *h, double phi, 
		      double **ret_x, int *ret_n, int *ret_z)
{
  int hi, lo, mid;

  if (! h->is_full) ESL_EXCEPTION(eslEINVAL, "not a full histogram");
  esl_histogram_sort(h);

  if      (h->n         == 0)   mid = h->n;  /* we'll return NULL, 0, n */  
  else if (h->x[0]       > phi) mid = 0;     /* we'll return x, n, 0    */
  else if (h->x[h->n-1] <= phi) mid = h->n;  /* we'll return NULL, 0, n */
  else /* binary search, faster than a brute force scan */
    {
      lo = 0;
      hi = h->n-1; /* know hi>0, because above took care of n=0 and n=1 cases */
      while (1) {
	mid = (lo + hi + 1) / 2;  /* +1 makes mid round up, mid=0 impossible */
	if      (h->x[mid]  <= phi) lo = mid; /* we're too far left  */
	else if (h->x[mid-1] > phi) hi = mid; /* we're too far right */
	else break;		              /* ta-da! */
      }
    }

  if (ret_x != NULL) *ret_x = h->x + mid;
  if (ret_n != NULL) *ret_n = h->n - mid;
  if (ret_z != NULL) *ret_z = mid;
  h->is_done = TRUE;
  return eslOK;
}


/* Function:  esl_histogram_GetTailByMass()
 * Incept:    SRE, Sun Jan 29 17:56:37 2006 [St. Louis]
 *
 * Purpose:   Given a full histogram <h>, retrieve the data values in
 *            the right (high scoring) tail, as a pointer <ret_x>
 *            to an array of <ret_n> values indexed <[0..n-1]> from
 *            lowest to highest score. The tail is defined by a
 *            given mass fraction threshold <pmass>; the mass in the returned
 *            tail is $\leq$ this threshold. <pmass> is a probability,
 *            so it must be $\geq 0$ and $\leq 1$.
 *            
 *            Optionally, the number of values in the rest of the
 *            histogram can be returned in <ret_z>. This is useful
 *            if you are going to fit the tail as a left-censored
 *            distribution.
 *            
 *            <ret_x> is a pointer to internal memory in <h>. 
 *            The histogram <h> remains responsible for its storage,
 *            which will be free'd when you call <esl_histogram_Destroy()>.
 *            As a consequence, you can only call 
 *            <esl_histogram_GetTailByMass()> after you have finished
 *            collecting data. Subsequent calls to <esl_histogram_Add()>
 *            will fail.
 *
 * Args:      h     - histogram to retrieve the tail from
 *            pmass - fractional mass threshold; tail contains <= pmass
 *            ret_x - optRETURN: ptr to vector of data values [0..n-1]
 *            ret_n - optRETURN: number of data values in tail x
 *            ret_z - optRETURN: number of data values not in tail
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if the histogram is not a full histogram, 
 *            or <pmass> is not a probability.
 */
int
esl_histogram_GetTailByMass(ESL_HISTOGRAM *h, double pmass,
			    double **ret_x, int *ret_n, int *ret_z)
{
  int n;

  if (! h->is_full) 
    ESL_EXCEPTION(eslEINVAL, "not a full histogram");
  if (pmass < 0. || pmass > 1.) 
    ESL_EXCEPTION(eslEINVAL, "pmass not a probability");

  esl_histogram_sort(h);

  n = (int) ((float) h->n * pmass); /* rounds down, guaranteeing <= pmass */

  if (ret_x != NULL) *ret_x = h->x + (h->n - n);
  if (ret_n != NULL) *ret_n = n;
  if (ret_z != NULL) *ret_z = h->n - n;
  h->is_done = TRUE;
  return eslOK;
}






/*****************************************************************
 * Setting expected counts
 *****************************************************************/ 

/* Function:  esl_histogram_SetExpect()
 * Incept:    SRE, Wed Aug 17 17:36:58 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h> containing some number of empirically
 *            observed binned counts, and a pointer to a function <(*cdf)()>
 *            that describes the expected cumulative distribution function 
 *            (CDF) for the complete data, conditional on some parameters 
 *            <params>; calculate the expected counts in each bin of the 
 *            histogram, and hold that information internally in the structure.
 *            
 *            The caller provides a function <(*cdf)()> that calculates
 *            the CDF via a generic interface, taking only two
 *            arguments: a quantile <x> and a void pointer to whatever
 *            parameters it needs, which it will cast and interpret.
 *            The <params> void pointer to the given parameters is
 *            just passed along to the generic <(*cdf)()> function. The
 *            caller will probably implement this <(*cdf)()> function as
 *            a wrapper around its real CDF function that takes
 *            explicit (non-void-pointer) arguments.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; state of <h> is preserved.
 */
int
esl_histogram_SetExpect(ESL_HISTOGRAM *h, 
			double (*cdf)(double x, void *params), void *params)
{
  int    status;
  int    i;
  double ai,bi;			/* ai < x <= bi : lower,upper bounds in bin */

  if (h->expect == NULL) 
    ESL_ALLOC(h->expect, sizeof(double) * h->nb);

  for (i = 0; i < h->nb; i++)
    {
      ai = esl_histogram_Bin2LBound(h, i);
      bi = esl_histogram_Bin2UBound(h, i);
      h->expect[i] = h->Nc * ( (*cdf)(bi, params) - (*cdf)(ai, params) );

      if (h->emin != -1 && h->expect[i] > 0.) h->emin = i;
    }

  h->is_done = TRUE;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_histogram_SetExpectedTail()
 * Incept:    SRE, Mon Jan 30 08:57:57 2006 [St. Louis]
 *
 * Purpose:   Given a histogram <h>, and a pointer to a generic function
 *            <(*cdf)()> that describes the expected cumulative
 *            distribution function for the right (high-scoring) tail
 *            starting at <base_val> (all expected <x> $>$ <base_val>) and
 *            containing a fraction <pmass> of the complete data
 *            distribution (<pmass> $\geq 0$ and $\leq 1$);
 *            set the expected binned counts for all complete bins
 *            $\geq$ <base_val>. 
 *            
 *            If <base_val> falls within a bin, that bin is considered
 *            to be incomplete, and the next higher bin is the starting
 *            point. 
 *           
 * Args:      h          - finished histogram
 *            base_val   - threshold for the tail: all expected x > base_val
 *            pmass      - fractional mass in the tail: 0 <= pmass <= 1
 *            cdf        - generic-interface CDF function describing the tail
 *            params     - void pointer to parameters for (*cdf)()
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on memory allocation failure.
 */
int
esl_histogram_SetExpectedTail(ESL_HISTOGRAM *h, double base_val, double pmass,
			      double (*cdf)(double x, void *params), 
			      void *params)
{
  int status;
  int b;
  double ai, bi;

  if (h->expect == NULL)  ESL_ALLOC(h->expect, sizeof(double) * h->nb);

  h->emin = 1 + esl_histogram_Score2Bin(h, base_val);
  esl_vec_DSet(h->expect, h->emin, 0.);

  for (b = h->emin; b < h->nb; b++)
    {
      ai = esl_histogram_Bin2LBound(h, b);
      bi = esl_histogram_Bin2UBound(h, b);
      h->expect[b] = pmass * (double) h->Nc * 
	             ( (*cdf)(bi, params) - (*cdf)(ai, params) );
    }
  
  h->tailbase   = base_val;
  h->tailmass   = pmass;
  h->is_tailfit = TRUE;
  h->is_done    = TRUE;
  return eslOK;

 ERROR:
  return status;
}




/*****************************************************************
 * Output and display of binned data.
 *****************************************************************/ 

/* Function:  esl_histogram_Print() 
 * Incept:    SRE, Sat Jul  2 16:03:37 2005 [St. Louis]
 *
 * Purpose:   Print a "prettified" display histogram <h> to a file 
 *            pointer <fp>.
 *            Deliberately a look-and-feel clone of Bill Pearson's 
 *            excellent FASTA output.
 *            
 *            Also displays expected binned counts, if they've been
 *            set.
 *            
 *            Display will only work well if the bin width (w) is 0.1 or more,
 *            because the score labels are only shown to one decimal point.
 * 
 * Args:      fp     - open file to print to (stdout works)
 *            h      - histogram to print
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_Print(FILE *fp, ESL_HISTOGRAM *h)
{
  int    i;
  double x;
  int maxbar;
  int imode;
  int units;
  int num;
  char buffer[81];		/* output line buffer */
  int  pos;			/* position in output line buffer */
  int  ilowbound, lowcount;	/* cutoffs on the low side  */
  int  ihighbound, highcount;	/* cutoffs on the high side */
  int  emptybins = 3;

  /* Find out how we'll scale the histogram.  We have 58 characters to
   * play with on a standard 80-column terminal display: leading "%6.1f
   * %6d %6d|" occupies 21 chars.  Save the peak position, we'll use
   * it later.
   */
  maxbar = 0;
  imode  = 0;
  for (i = 0; i < h->nb; i++)
    if (h->obs[i] > maxbar) 
      {
	maxbar  = h->obs[i];     /* max height    */
	imode   = i;
      }

  /* Truncate histogram display on both sides, ad hoc fashion.
   * Start from the peak; then move out until we see <emptybins> empty bins,
   * and stop.
   */
  for (num = 0, ihighbound = imode; ihighbound < h->imax; ihighbound++)
    {
      if (h->obs[ihighbound] > 0) { num = 0; continue; } /* reset */
      if (++num == emptybins)     { break;             } /* stop  */
    }
  for (num = 0, ilowbound = imode; ilowbound > h->imin; ilowbound--)
    {
      if (h->obs[ilowbound] > 0)  { num = 0; continue; } /* reset */
      if (++num == emptybins)     { break;             } /* stop  */
    }

		/* collect counts outside of bounds */
  for (lowcount = 0, i = h->imin; i < ilowbound; i++)
    lowcount += h->obs[i];
  for (highcount = 0, i = h->imax; i > ihighbound; i--)
    highcount += h->obs[i];

		/* maxbar might need to be raised now; then set our units  */
  if (lowcount  > maxbar) maxbar = lowcount;
  if (highcount > maxbar) maxbar = highcount;
  units = ((maxbar-1)/ 58) + 1;

  /* Print the histogram
   */
  fprintf(fp, "%6s %6s %6s  (one = represents %d sequences)\n", 
	  "score", "obs", "exp", units);
  fprintf(fp, "%6s %6s %6s\n", "-----", "---", "---");
  buffer[80] = '\0';
  buffer[79] = '\n';
  for (i = h->imin; i <= h->imax; i++)
    {
      memset(buffer, ' ', 79 * sizeof(char));
      x = i*h->w + h->bmin;

      /* Deal with special cases at edges
       */
      if      (i < ilowbound)  continue;
      else if (i > ihighbound) continue;
      else if (i == ilowbound && i != h->imin) 
	{
	  sprintf(buffer, "<%5.1f %6d %6s|", x+h->w, lowcount, "-");
	  if (lowcount > 0) {
	    num = 1+(lowcount-1) / units;
	    for (pos = 21; num > 0; num--)  buffer[pos++] = '=';
	  }
	  fputs(buffer, fp);
	  continue;
	}
      else if (i == ihighbound && i != h->imax)
	{
	  sprintf(buffer, ">%5.1f %6d %6s|", x, highcount, "-");
	  if (highcount > 0) {
	    num = 1+(highcount-1) / units;
	    for (pos = 21; num > 0; num--)  buffer[pos++] = '=';
	  }
	  fputs(buffer, fp);
	  continue;
	}

      /* Deal with most cases
       */
      if (h->expect != NULL) 
	sprintf(buffer, "%6.1f %6d %6d|", 
		x, h->obs[i], (int) h->expect[i]);
      else
	sprintf(buffer, "%6.1f %6d %6s|", x, h->obs[i], "-");
      buffer[21] = ' ';		/* sprintf writes a null char */

      /* Mark the histogram bar for observed hits
       */ 
      if (h->obs[i] > 0) {
	num = 1 + (h->obs[i]-1) / units;
	for (pos = 21; num > 0; num--)  buffer[pos++] = '=';
      }
	  
      /* Mark the theoretically expected value
       * (The test > 0. also suffices to remove any censored region.)
       */
      if (h->expect != NULL && h->expect[i] > 0.)
	{
	  pos = 21 + (int)(h->expect[i]-1) / units;
	  if (pos >= 78) pos = 78; /* be careful of buffer bounds */
	  buffer[pos] = '*';
	}

      /* Print the line
       */
      fputs(buffer, fp);
    }

  return eslOK;
}
  
/* Function:  esl_histogram_Plot()
 * Incept:    SRE, Mon Jan 30 11:09:01 2006 [St. Louis]
 *
 * Purpose:   Print observed (and expected, if set) binned counts
 *            in a histogram <h> to open file pointer <fp>
 *            in xmgrace XY input file format.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_Plot(FILE *fp, ESL_HISTOGRAM *h)
{
  int    i;
  double x;

  /* First data set is the observed histogram
   */
  for (i = h->imin; i <= h->imax; i++)
    if (h->obs[i] > 0)
      {
	x = esl_histogram_Bin2LBound(h,i);
	fprintf(fp, "%f %d\n", x, h->obs[i]);
      }
  fprintf(fp, "&\n");

  /* Second data set is the theoretical (expected) histogram
   */
  if (h->expect != NULL)
    {
      for (i = 0; i < h->nb; i++)
	if (h->expect[i] > 0.)	/* >0 suffices to remove censored region */
	  {
	    x = esl_histogram_Bin2LBound(h,i);
	    fprintf(fp, "%.2f %g\n", x, h->expect[i]);
	  }
      fprintf(fp, "&\n");
    }
  return eslOK;
}

/* Function:  esl_histogram_PlotSurvival()
 * Incept:    SRE, Mon Jan 30 11:11:05 2006 [St. Louis]
 *
 * Purpose:   Given a histogram <h>, output the observed (and
 *            expected, if available) survival function $P(X>x)$
 *            to file pointer <fp> in xmgrace XY input file format.
 *            
 *            One point is plotted per bin, so the narrower the
 *            bin width, the more smooth and accurate the resulting
 *            plots will be.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_PlotSurvival(FILE *fp, ESL_HISTOGRAM *h)
{
  int i;
  double c;
  double ai;
  
  /* The observed binned counts:
   */
  c = 0.;
  for (i = h->imax; i >= h->imin; i--)
    {
      if (h->obs[i] > 0) {
	c   += h->obs[i];
	ai = esl_histogram_Bin2LBound(h, i);
	fprintf(fp, "%f\t%g\n", ai, c / (double) h->Nc);
      }
    }
  fprintf(fp, "&\n");

  /* The expected binned counts:
   */
  if (h->expect != NULL) 
    {
      c = 0.;
      for (i = h->nb-1; i >= 0; i--)
	{
	  if (h->expect[i] > 0.) {
	    c += h->expect[i];
	    ai = esl_histogram_Bin2LBound(h, i);
	    fprintf(fp, "%f\t%g\n", ai, c / (double) h->Nc);
	  }
	}
      fprintf(fp, "&\n");
    }
  return eslOK;
}

/* Function:  esl_histogram_PlotQQ()
 * Incept:    SRE, Sat Aug 20 14:15:01 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h> containing an empirically observed
 *            distribution, and a pointer to a function <(*invcdf)()>
 *            for an expected inverse cumulative distribution
 *            function conditional on some parameters <params>;
 *            output a Q-Q plot in xmgrace XY format to file <fp>.
 *            
 *            Same domain limits as goodness-of-fit testing: output
 *            is restricted to overlap between observed data (excluding
 *            any censored data) and expected data (which may be limited
 *            if only a tail was fit).
 *
 * Returns:   <eslOK> on success.
 */
int
esl_histogram_PlotQQ(FILE *fp, ESL_HISTOGRAM *h, 
		     double (*invcdf)(double x, void *params), void *params)
{
  int    i;
  double cdf;
  double bi;
  int    bbase;
  double sum;

  /* on censored data, start counting observed cdf at z, not 0
   */
  if (h->dataset_is == TRUE_CENSORED || h->dataset_is == VIRTUAL_CENSORED)
    sum = h->z; 
  else
    sum = 0.;

  /* Determine smallest bin included in goodness of fit eval
   */
  bbase = h->cmin;
  if (h->is_tailfit && h->emin > bbase) bbase = h->emin;
  for (i = h->cmin; i < bbase; i++) sum += (double) h->obs[i];
  
  /* The q-q plot:
   */
  for (i = bbase; i < h->imax; i++) /* avoid last bin where upper cdf=1.0 */
    {
      sum += (double) h->obs[i];
      cdf = sum / (double) h->Nc;

      if (h->is_tailfit) cdf = (cdf + h->tailmass - 1.) / (h->tailmass);

      bi = esl_histogram_Bin2UBound(h, i);
      fprintf(fp, "%f\t%f\n", bi, (*invcdf)(cdf, params));
    }
  fprintf(fp, "&\n");

  /* Plot a 45-degree expected QQ line:
   */
  bi = esl_histogram_Bin2LBound(h, bbase);
  fprintf(fp, "%f\t%f\n", bi,  bi);
  fprintf(fp, "%f\t%f\n", h->xmax, h->xmax);
  fprintf(fp, "&\n");

  return eslOK;
}



/* Function:  esl_histogram_Goodness()
 * Incept:    SRE, Wed Aug 17 12:46:05 2005 [St. Louis]
 *
 * Purpose:   Given a histogram <h> with observed and expected counts,
 *            where, for the expected counts, <nfitted> ($\geq 0$)
 *            parameters were fitted (and thus should be subtracted
 *            from the degrees of freedom);
 *            Perform a G-test and/or a $\chi^2$ test for goodness of 
 *            fit between observed and expected, and optionally return
 *            the number of bins the data were sorted into
 *            (<ret_bins>), the G statistic and its probability (<ret_G> and
 *            <ret_Gp>), and the $\chi^2$ statistic and its probability
 *            (<ret_X2> and <ret_X2p>). 
 *            
 *            If a goodness-of-fit probability is less than some threshold
 *            (usually taken to be 0.01 or 0.05), that is considered to
 *            be evidence that the observed data are unlikely to be consistent
 *            with the tested distribution.
 *            
 *            The two tests should give similar
 *            probabilities. However, both tests are sensitive to
 *            arbitrary choices in how the data are binned, and
 *            neither seems to be on an entirely sound theoretical footing.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if expected counts have not been set in
 *            the histogram; <eslERANGE> or <eslENOHALT> on different internal
 *            errors that can arise in calculating the probabilities;
 *            <eslEMEM> on internal allocation failure.
 */
int
esl_histogram_Goodness(ESL_HISTOGRAM *h, 
		       int nfitted, int *ret_nbins,
		       double *ret_G,  double *ret_Gp,
		       double *ret_X2, double *ret_X2p)
{
  int     *obs  = NULL;		/* observed in bin i, [0..nb-1]   */
  double  *exp  = NULL;		/* expected in bin i, [0..nb-1]   */
  double  *topx = NULL;		/* all values in bin i <= topx[i] */
  int      nb;			/* # of re-bins                   */
  int      minc;		/* minimum target # of counts/bin */
  int      i,b;
  double   G, Gp;
  double   X2, X2p;
  double   tmp;
  int      status;
  int      bbase;
  int      hmax;
  int      nobs;
  double   nexp;

  if (h->expect == NULL) ESL_EXCEPTION(eslEINVAL, "no expected counts in that histogram");

  /* Determine the smallest histogram bin included in 
   * the goodness of fit evaluation.
   */
  bbase = h->cmin;		
  if (h->is_tailfit && h->emin > bbase) bbase = h->emin;
  
  /* How many observed total counts are in the evaluated range,
   * and what is the maximum in any given histogram bin?
   */
  nobs = 0;
  hmax = 0;
  for (i = bbase; i <= h->imax; i++)
    {
      nobs += h->obs[i];
      if (h->obs[i] > hmax) hmax = h->obs[i];
    }

  /* Figure out how many eval bins we'd like to have, then allocate
   * for re-binning.
   * Number of bins for goodness-of-fit tests like G and X^2 
   * is crucial but arbitrary, unfortunately. Some literature suggests
   * using 2*n^{0.4}, which gives:
   *        n    nbins     #/bin
   *    -----    ------   ------
   *     1000      31       32
   *    10000      79      127
   *   100000     200      500
   *  1000000     502     1992
   *  
   * The most important thing seems to be to get the # of counts
   * in each bin to be roughly equal.
   */
  nb   = 2* (int) pow((double) nobs, 0.4); /* "desired" nb. */
  minc = 1 + nobs / (2*nb);	/* arbitrarily set min = 1/2 of the target # */
  ESL_ALLOC(obs,  sizeof(int)    * (nb*2+1)); /* final nb must be <= 2*nb+1 */
  ESL_ALLOC(exp,  sizeof(double) * (nb*2+1));
  ESL_ALLOC(topx, sizeof(double) * (nb*2+1));

  /* Determine the observed counts in each bin: that is, partition 
   * the <sum> in the evaluated region.
   * Sweep left to right on the histogram bins,
   * collecting sum of counts, dropping the sum into the next re-bin 
   * whenever we have more than <minc> counts.
   */
  nobs = nexp = 0;
  for (i = 0, b = bbase; b <= h->imax; b++) 
    {
      nobs += h->obs[b];
      nexp += h->expect[b];

      /* if we have enough counts, drop into bin i: */
      if (nobs >= minc && nexp >= minc) {
	ESL_DASSERT1( (i < (nb*2+1)) );
	obs[i]  = nobs;
	exp[i]  = nexp;
	topx[i] = esl_histogram_Bin2UBound(h,b);
	nobs = nexp = 0;
	i++;
      }
    }
  obs[i-1]  += nobs;		/* add the right tail to final bin */
  exp[i-1]  += nexp;
  topx[i-1]  = esl_histogram_Bin2UBound(h, h->imax);
  nb         = i;		/* nb is now actual # of bins, not target */

  /* Calculate the X^2 statistic: \sum (obs_i - exp_i)^2 / exp_i */
  X2 = 0.;
  for (i = 0; i < nb; i++)
    {
      tmp = obs[i] - exp[i];
      X2 += tmp*tmp / exp[i];
    }
  /* X^2 is distributed approximately chi^2. */
  if (nb-nfitted >= 0 && X2 != eslINFINITY)
    {
      status = esl_stats_ChiSquaredTest(nb-nfitted, X2, &X2p);
      if (status != eslOK) return status;
    }
  else X2p = 0.;

  /* The G test assumes that #exp=#obs (the X^2 test didn't).
   * If that's not true, renormalize to make it so. 
   */
  nobs = nexp = 0;
  for (i = 0; i < nb; i++) 
    {
      nobs += obs[i];
      nexp += exp[i];
    }
  for (i = 0; i < nb; i++)
    exp[i] = exp[i] * (double) nobs / nexp;
  
  /* Calculate the G statistic: 2 * LLR  */
  G = 0.;
  for (i = 0; i < nb; i++)
    G += (double) obs[i] * log ((double) obs[i] / exp[i]);
  G *= 2;
  
  /* G is distributed approximately as \chi^2.
   * -1 is because total #obs=#exp (which is must be)
   */
  ESL_DASSERT1( (G >= 0.));
  if (nb-nfitted-1 >= 0 && G != eslINFINITY)
    {
      status = esl_stats_ChiSquaredTest(nb-nfitted-1, G, &Gp);
      if (status != eslOK) return status;
    }
  else Gp = 0.;

  if (ret_nbins != NULL) *ret_nbins = nb;
  if (ret_G     != NULL) *ret_G     = G;
  if (ret_Gp    != NULL) *ret_Gp    = Gp;
  if (ret_X2    != NULL) *ret_X2    = X2;
  if (ret_X2p   != NULL) *ret_X2p   = X2p;
  free(obs);
  free(exp);
  free(topx);
  return eslOK;

 ERROR:
  if (ret_nbins != NULL) *ret_nbins = 0;
  if (ret_G     != NULL) *ret_G     = 0.;
  if (ret_Gp    != NULL) *ret_Gp    = 0.;
  if (ret_X2    != NULL) *ret_X2    = 0.;
  if (ret_X2p   != NULL) *ret_X2p   = 0.;
  if (obs  != NULL) free(obs);
  if (exp  != NULL) free(exp);
  if (topx != NULL) free(topx);
  return status;
}



/*****************************************************************
 * Five example main()'s for five use cases:
 *    - complete data, fit to complete Gumbel
 *    - complete data, high scores fit as censored Gumbel
 *    - complete data, high scores fit to exponential tail
 *    - censored data, fit as censored Gumbel
 *    - complete data, binned, high scores fit to exponential tail
 *
 * (These same five cases are tested by ./test -t1 through ./test -t5.)
 *****************************************************************/
/* Case 1. Complete data fit to complete Gumbel.
 * compile: gcc -I. -L. -o example -DeslHISTOGRAM_EXAMPLE1 esl_histogram.c -leasel -lm
 * run:     ./example 
 */
#ifdef eslHISTOGRAM_EXAMPLE1
/*::cexcerpt::histogram_example1::begin::*/
#include <easel.h>
#include <esl_random.h>
#include <esl_histogram.h>
#include <esl_gumbel.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r  = esl_randomness_CreateTimeseeded();
  ESL_HISTOGRAM  *h  = esl_histogram_CreateFull(-100, 100, 0.2);
  int     nsamples    = 10000;
  double  mu          = 10.0;
  double  lambda      = 0.8;
  double  params[2];
  int     i;
  double  x;
  double *xv;
  int     n;
  double  G, Gp, X2, X2p;

  for (i = 0; i < nsamples; i++) {
    x = esl_gumbel_Sample(r, mu, lambda);
    esl_histogram_Add(h, x);
  }

  esl_histogram_GetData(h, &xv, &n);
  esl_gumbel_FitComplete(xv, n, &mu, &lambda);

  params[0] = mu;
  params[1] = lambda;
  esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, &params);

  esl_histogram_Print(stdout, h);
  esl_histogram_Goodness(h, 0, NULL, &G, &Gp, &X2, &X2p);
  printf("G   = %f  p = %f\n", G, Gp);
  printf("X^2 = %f  p = %f\n", X2, X2p);

  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::histogram_example1::end::*/
#endif /*eslHISTOGRAM_EXAMPLE1*/



/* Case 2. complete data, high scores fit as censored Gumbel 
 * compile: gcc -I. -L. -o example -DeslHISTOGRAM_EXAMPLE2 esl_histogram.c -leasel -lm
 * run:     ./example 
 */
#ifdef eslHISTOGRAM_EXAMPLE2
/*::cexcerpt::histogram_example2::begin::*/
#include <easel.h>
#include <esl_random.h>
#include <esl_histogram.h>
#include <esl_gumbel.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r  = esl_randomness_CreateTimeseeded();
  ESL_HISTOGRAM  *h  = esl_histogram_CreateFull(-100, 100, 0.2);
  int     nsamples    = 10000;
  double  mu          = 10.0;
  double  lambda      = 0.8;
  double  params[2];
  int     i;
  double  x;
  double *xv;
  int     n, z;
  double  G, Gp, X2, X2p;

  for (i = 0; i < nsamples; i++) {
    x = esl_gumbel_Sample(r, mu, lambda);
    esl_histogram_Add(h, x);
  }

  esl_histogram_GetTailByMass(h, 0.5, &xv, &n, &z); /* fit to right 50% */
  esl_gumbel_FitCensored(xv, n, z, xv[0], &mu, &lambda);

  params[0] = mu;
  params[1] = lambda;
  esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, &params);

  esl_histogram_Print(stdout, h);
  esl_histogram_Goodness(h, 0, NULL, &G, &Gp, &X2, &X2p);
  printf("G   = %f  p = %f\n", G, Gp);
  printf("X^2 = %f  p = %f\n", X2, X2p);

  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::histogram_example2::end::*/
#endif /*eslHISTOGRAM_EXAMPLE2*/


/* Case 3. complete data, high scores fit to exponential tail
 * compile: gcc -I. -L. -o example -DeslHISTOGRAM_EXAMPLE3 esl_histogram.c -leasel -lm
 * run:     ./example 
 */
#ifdef eslHISTOGRAM_EXAMPLE3
/*::cexcerpt::histogram_example3::begin::*/
#include <easel.h>
#include <esl_random.h>
#include <esl_histogram.h>
#include <esl_gumbel.h>
#include <esl_exponential.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r  = esl_randomness_CreateTimeseeded();
  ESL_HISTOGRAM  *h  = esl_histogram_CreateFull(-100, 100, 0.2);
  int     nsamples    = 10000;
  double  mu          = 10.0;
  double  lambda      = 0.8;
  double  params[2];
  int     i;
  double  x;
  double *xv;
  int     n;
  double  G, Gp, X2, X2p;

  for (i = 0; i < nsamples; i++) {
    x = esl_gumbel_Sample(r, mu, lambda);
    esl_histogram_Add(h, x);
  }

  esl_histogram_GetTailByMass(h, 0.1, &xv, &n, NULL); /* fit to 10% tail */
  esl_exp_FitComplete(xv, n, &mu, &lambda);

  params[0] = mu;
  params[1] = lambda;
  esl_histogram_SetExpectedTail(h, mu, 0.1, &esl_exp_generic_cdf, &params);

  esl_histogram_Print(stdout, h);
  esl_histogram_Goodness(h, 0, NULL, &G, &Gp, &X2, &X2p);
  printf("G   = %f  p = %f\n", G, Gp);
  printf("X^2 = %f  p = %f\n", X2, X2p);

  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::histogram_example3::end::*/
#endif /*eslHISTOGRAM_EXAMPLE3*/

/* Case 4. censored data, high scores fit as a censored Gumbel tail
 * compile: 
     gcc -I. -L. -o example -DeslHISTOGRAM_EXAMPLE4 esl_histogram.c -leasel -lm
 * run:     ./example 
 */
#ifdef eslHISTOGRAM_EXAMPLE4
/*::cexcerpt::histogram_example4::begin::*/
#include <easel.h>
#include <esl_random.h>
#include <esl_histogram.h>
#include <esl_gumbel.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r  = esl_randomness_CreateTimeseeded();
  ESL_HISTOGRAM  *h  = esl_histogram_CreateFull(-100, 100, 0.2);
  int     nsamples    = 10000;
  double  mu          = 10.0;
  double  lambda      = 0.8;
  double  phi         = 9.0;
  double  params[2];
  int     i;
  double  x;
  double *xv;
  int     n, z;
  double  G, Gp, X2, X2p;

  z = 0;
  for (i = 0; i < nsamples; i++) {
    x = esl_gumbel_Sample(r, mu, lambda);
    if (x > phi) esl_histogram_Add(h, x);
    else         z++;
  }

  esl_histogram_GetData(h, &xv, &n);
  esl_gumbel_FitCensored(xv, n, z, phi, &mu, &lambda);

  params[0] = mu;
  params[1] = lambda;
  esl_histogram_SetExpect(h, &esl_gumbel_generic_cdf, &params);

  esl_histogram_Print(stdout, h);
  esl_histogram_Goodness(h, 0, NULL, &G, &Gp, &X2, &X2p);
  printf("G   = %f  p = %f\n", G, Gp);
  printf("X^2 = %f  p = %f\n", X2, X2p);

  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::histogram_example4::end::*/
#endif /*eslHISTOGRAM_EXAMPLE4*/

/* Case 5. complete data, binned high scores fit to exponential tail
 * compile:
     gcc -I. -L. -o example -DeslHISTOGRAM_EXAMPLE5 esl_histogram.c -leasel -lm
 * run:     ./example 
 */
#ifdef eslHISTOGRAM_EXAMPLE5
/*::cexcerpt::histogram_example5::begin::*/
#include <easel.h>
#include <esl_random.h>
#include <esl_histogram.h>
#include <esl_gumbel.h>
#include <esl_exponential.h>

int
main(int argc, char **argv)
{
  ESL_RANDOMNESS *r  = esl_randomness_CreateTimeseeded();
  ESL_HISTOGRAM  *h  = esl_histogram_Create(-100, 100, 1.0);
  int     nsamples    = 10000;
  double  mu          = 10.0;
  double  lambda      = 0.8;
  double  params[2];
  int     i;
  double  x;
  double  actual_mass;
  double  G, Gp, X2, X2p;

  for (i = 0; i < nsamples; i++) {
    x = esl_gumbel_Sample(r, mu, lambda);
    x = ceil(x);      /* crudely simulate an x of limited precision */
    esl_histogram_Add(h, x);
  }

  esl_histogram_SetTailByMass(h, 0.1, &actual_mass);
  esl_histogram_DeclareRounding(h);
  esl_exp_FitCompleteBinned(h, &mu, &lambda);

  params[0] = mu;
  params[1] = lambda;
  esl_histogram_SetExpectedTail(h, mu, actual_mass, &esl_exp_generic_cdf, &params);

  esl_histogram_Print(stdout, h);
  esl_histogram_Goodness(h, 0, NULL, &G, &Gp, &X2, &X2p);
  printf("G   = %f  p = %f\n", G, Gp);
  printf("X^2 = %f  p = %f\n", X2, X2p);

  esl_histogram_Destroy(h);
  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::histogram_example5::end::*/
#endif /*eslHISTOGRAM_EXAMPLE5*/



/*****************************************************************
 * Test driver code
 *****************************************************************/
#ifdef eslHISTOGRAM_TESTDRIVE
/* compile: 
 *   gcc -g -Wall -I. -L. -o test -DeslHISTOGRAM_TESTDRIVE esl_histogram.c -leasel -lm
 * run:     
 *   ./test -t1; ./test -t2; ./test -t3; ./test -t4; ./test -t5
 *   
 * Some suggestions for manual testing:
 *   ./test -t1 -j1 -v --surv test.xy; xmgrace test.xy          
 *        examine survivor plot fit, for -t1 
 *        do -t2 thru -t5 too
 *
 *   ./test -t1 --j1 -v -qq test.xy; xmgrace test.xy          
 *        examine QQ plot fit, for -t1 
 *        do -t2 thru -t5 too
 *        
 *   ./test -t1 -v > foo
 *   grep "^Estimated" foo | awk '{print $9}' | sort -g > test.xy
 *        Look for straight line fit to G-test p values.
 *        sub $9->$13 for chi-squared
 *        sub Estimated -> Parametric for the parametric fits
 */
#include <stdio.h>
#include <stdlib.h>

#include <easel.h>
#include <esl_stats.h>
#include <esl_gumbel.h>
#include <esl_exponential.h>
#include <esl_random.h>
#include <esl_getopts.h>

static ESL_OPTIONS options[] = {
  /* name         type      default   env_var   range   toggles     reqs   incompat */
  { "-j",       eslARG_INT,   "100",  NULL,     "n>0",     NULL,  NULL,   NULL },
  { "-m",       eslARG_INT,     "0",  NULL,    "n>=0",     NULL,  NULL,   NULL },
  { "-n",       eslARG_INT, "10000",  NULL,     "n>0",     NULL,  NULL,   NULL },
  { "-t",       eslARG_INT,     "1",  NULL, "1<=n<=5",     NULL,  NULL,   NULL },
  { "-v",       eslARG_NONE,  FALSE,  NULL,      NULL,     NULL,  NULL,   NULL },
  { "--ascii",  eslARG_STRING, NULL,  NULL,      NULL,     NULL,  NULL,   NULL },
  { "--cmass",  eslARG_REAL,  "0.7",  NULL, "0<=x<=1",     NULL,  NULL,   NULL },
  { "--lambda", eslARG_REAL,  "0.8",  NULL,     "x>0",     NULL,  NULL,   NULL },
  { "--mu",     eslARG_REAL, "10.0",  NULL,      NULL,     NULL,  NULL,   NULL },
  { "--phi",    eslARG_REAL, "10.0",  NULL,      NULL,     NULL,  NULL,   NULL },
  { "--plot",   eslARG_STRING, NULL,  NULL,      NULL,     NULL,  NULL,   NULL },
  { "--qq",     eslARG_STRING, NULL,  NULL,      NULL,     NULL,  NULL,   NULL },
  { "--surv",   eslARG_STRING, NULL,  NULL,      NULL,     NULL,  NULL,   NULL },
  { "--tail",   eslARG_REAL,  "0.1",  NULL, "0<=x<=1",     NULL,  NULL,   NULL },
  { 0,0,0,0,0,0,0,0 },
};

static char usage[] = "\
Usage: ./test [-options]\n\
where options are:\n\
\n\
  -t <n>  : choice of test scenario [-t1 is default]:\n\
     -t1    - complete data, fit to complete Gumbel\n\
     -t2    - complete data, high scores fit as censored Gumbel\n\
     -t3    - complete data, high scores fit to exponential tail\n\
     -t4    - censored data, fit as censored Gumbel\n\
     -t5    - complete data, binned, high scores fit to exponential tail\n\
  -j <n>  : number of trials [100]\n\
  -n <n>  : number of training set samples [10000]\n\
  -m <n>  : number of independent test set samples [default: use training set]\n\
\n\
Output options: (best to use -j1 when saving files)\n\
  -v          : report verbose output [default is silent success/failure]\n\
  --ascii <f> : output ASCII histogram to file <f>\n\
  --plot <f>  : output histogram (densities) to xmgrace file <f>\n\
  --qq <f>    : output Q-Q goodness of fit plot to xmgrace file <f>\n\
  --surv <f>  : output survival (P(X>x)) plots to xmgrace file <f>\n\
\n\
  --mu <x>      : set Gumbel mu parameter to <x> [10.0]\n\
  --lambda <x>  : set Gumbel lambda parameter to <x> [0.8]\n\
  --phi <x>     : set data censoring threshold to <x> [10.0]\n\
  --cmass <x>   : set virtual censoring mass to <x> [0.7]\n\
  --tail <x>    : set tail mass for censored/tail fitting to <x> [0.1]\n\
";


static int
binmacro_test(void)
{
  ESL_HISTOGRAM *h = esl_histogram_Create(-100, 100, 1.0);
  double trialx[3]  = { -42.42, 0, 42.42 };
  double x, ai, bi;  
  int    i,b;

  /* test bin<->score conversion macros.
   */
  for (i = 0; i < 3; i++)
    {
      x  = trialx[i];
      b  = esl_histogram_Score2Bin(h, x);
      ai = esl_histogram_Bin2LBound(h, b);
      bi = esl_histogram_Bin2UBound(h, b);
      if (x <= ai || x > bi) {
	fprintf(stderr,
		"failed: (ai=%.1f) <= (x=%.2f) < (bi=%.1f) in bin %d, bin macro test\n",
		ai, x, bi, b);
	esl_histogram_Destroy(h);
	return 0;
      }
    }
  esl_histogram_Destroy(h);
  return 1;
}
int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go;
  ESL_RANDOMNESS *r;
  ESL_HISTOGRAM  *h;
  ESL_HISTOGRAM  *h1;
  double          p[2];		/* parametric mu, lambda */
  double          ep[2];	/* estimated mu, lambda  */
  double          avg_ep[2];	/* average estimated mu, lambda over many trials */
  int             ntrials, trial;
  int             ntrain, ntest;
  int             test_type;
  enum { COLLECT_COMPLETE, COLLECT_CENSORED }   cstrategy;
  enum { FIT_BINNED, FIT_SAMPLES }              bstrategy;
  enum { FIT_COMPLETE, FIT_CENSORED, FIT_TAIL}  fstrategy;
  double          phi;		/* censoring threshold   */
  int             z;
  double          cmass;
  double          tailmass, save_tailmass;
  int             nfitted;
  int             nbins;
  double          G, Gp, X2, X2p, minGp, minX2p;
  int             verbose;
  FILE           *outfp;
  char           *ascfile, *plotfile, *survfile, *qqfile;
  int     i;
  double  x;
  double *xv;
  int     n;

  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_GetIntegerOption(go, "-t",       &test_type);
  esl_opt_GetIntegerOption(go, "-j",       &ntrials);
  esl_opt_GetIntegerOption(go, "-n",       &ntrain);
  esl_opt_GetIntegerOption(go, "-m",       &ntest);
  esl_opt_GetBooleanOption(go, "-v",       &verbose);
  esl_opt_GetDoubleOption (go, "--cmass",  &cmass);
  esl_opt_GetDoubleOption (go, "--lambda", &(p[1]));
  esl_opt_GetDoubleOption (go, "--mu",     &(p[0]));
  esl_opt_GetDoubleOption (go, "--phi",    &phi);
  esl_opt_GetDoubleOption (go, "--tail",   &save_tailmass);
  esl_opt_GetStringOption (go, "--ascii",  &ascfile);
  esl_opt_GetStringOption (go, "--plot",   &plotfile);
  esl_opt_GetStringOption (go, "--qq",     &qqfile);
  esl_opt_GetStringOption (go, "--surv",   &survfile);
  esl_getopts_Destroy(go);

  r         = esl_randomness_Create(42);
  avg_ep[0] = 0.;
  avg_ep[1] = 0.;
  minGp     = 1.;
  minX2p    = 1.;
  tailmass  = save_tailmass;

  if (test_type == 1)
    {
      cstrategy = COLLECT_COMPLETE;
      bstrategy = FIT_SAMPLES;
      fstrategy = FIT_COMPLETE;
    }
  else if (test_type == 2)
    {
      cstrategy = COLLECT_COMPLETE;
      bstrategy = FIT_SAMPLES;
      fstrategy = FIT_CENSORED;
    }
  else if (test_type == 3)
    {
      cstrategy = COLLECT_COMPLETE;
      bstrategy = FIT_SAMPLES;
      fstrategy = FIT_TAIL;
    }
  else if (test_type == 4)
    {
      cstrategy = COLLECT_CENSORED;
      bstrategy = FIT_SAMPLES;
      fstrategy = FIT_CENSORED;
    }
  else if (test_type == 5)
    {
      cstrategy = COLLECT_COMPLETE;
      bstrategy = FIT_BINNED;
      fstrategy = FIT_TAIL;
    }


  for (trial = 0; trial < ntrials; trial++)
    {
      /* Collection of the training data in <h>.
       * Data set can either be complete, true censored, or virtual censored.
       */
      h = esl_histogram_CreateFull(-100, 100, 0.1);
      z = 0;
      for (i = 0; i < ntrain; i++) {
	x = esl_gumbel_Sample(r, p[0], p[1]);
	if (cstrategy != COLLECT_CENSORED || x > phi)
	  esl_histogram_Add(h, x);
	else
	  z++;
      }
      if (cstrategy == COLLECT_CENSORED)
	esl_histogram_DeclareCensoring(h, z, phi);

      /* Parameter fitting.
       * We test for four of twelve possible combinations of
       * collection strategy, binned vs. raw data, and complete,
       * censored, vs. tail fitting.
       *   1. complete Gumbel data, raw, fit to a Gumbel.
       *   2. complete Gumbel data, raw, tail fit as a censored Gumbel
       *   3. complete Gumbel data, raw, tail fit to an exponential tail
       *   4. censored Gumbel data, raw, censored fit to a Gumbel
       *   5  complete Gumbel data, binned, fit to an exponential tail.
       */
      if (cstrategy == COLLECT_COMPLETE &&
	  bstrategy == FIT_SAMPLES &&
	  fstrategy == FIT_COMPLETE)
	{
	  esl_histogram_GetData(h, &xv, &n);
	  esl_gumbel_FitComplete(xv, n, &(ep[0]), &ep[1]);
	}
      else if (cstrategy == COLLECT_COMPLETE &&
	       bstrategy == FIT_SAMPLES &&
	       fstrategy == FIT_CENSORED)
	{
	  esl_histogram_GetTailByMass(h, cmass, &xv, &n, &z);
	  esl_gumbel_FitCensored(xv, n, z, xv[0], &(ep[0]), &ep[1]);
	}
      else if (cstrategy == COLLECT_COMPLETE &&
	       bstrategy == FIT_SAMPLES &&
	       fstrategy == FIT_TAIL)
	{
	  esl_histogram_GetTailByMass(h, tailmass, &xv, &n, &z);
	  esl_exp_FitComplete(xv, n, &(ep[0]), &ep[1]);
	}
      else if (cstrategy == COLLECT_CENSORED &&
	       bstrategy == FIT_SAMPLES &&
	       fstrategy == FIT_CENSORED)
	{
	  esl_histogram_GetData(h, &xv, &n);
	  esl_gumbel_FitCensored(xv, n, h->z, h->phi, &(ep[0]), &ep[1]);
	}
      else if (cstrategy == COLLECT_COMPLETE &&
	       bstrategy == FIT_BINNED &&
	       fstrategy == FIT_TAIL)
	{
	  tailmass = save_tailmass; /* reset to original for each trial. */
	  esl_histogram_SetTailByMass(h, tailmass, &tailmass);
	  esl_exp_FitCompleteBinned(h, &(ep[0]), &ep[1]);
	}
      else
	ESL_EXCEPTION(eslEINVAL, "not a scenario we currently test");

      /* Keep track of average estimated mu, lambda
       * for automated testing purposes.
       */
      avg_ep[0] += ep[0] / (double) ntrials;
      avg_ep[1] += ep[1] / (double) ntrials;

      /* Test data can either be the same as the training data,
       * or a new test set.
       */
      if (ntest > 0)
	{
	  h1 = esl_histogram_CreateFull(-100.05, 100.05, 0.2);
	  z = 0;
	  for (i = 0; i < ntest; i++) {
	    x = esl_gumbel_Sample(r, p[0], p[1]);
	    if (cstrategy != COLLECT_CENSORED || x > phi)
	      esl_histogram_Add(h1, x);
	    else
	      z++;
	  }
	  if (cstrategy == COLLECT_CENSORED)
	    esl_histogram_DeclareCensoring(h, z, phi);
	}
      else h1 = h;
      

      /* Set expected binned counts in the test data, h1:
       */
      if (fstrategy == FIT_TAIL)
	esl_histogram_SetExpectedTail(h1, ep[0], tailmass, 
				      &esl_exp_generic_cdf, ep);
      else
	esl_histogram_SetExpect(h1, &esl_gumbel_generic_cdf, ep);

  
      /* Evaluate goodness-of-fit
       */
      nfitted =  (ntest == 0)? 2 : 0;
      esl_histogram_Goodness(h1, nfitted, &nbins, &G, &Gp, &X2, &X2p);

      /* Track minimum goodness of fit probs, for automated testing
       */
      if (Gp  < minGp)  minGp  = Gp;
      if (X2p < minX2p) minX2p = X2p;

      if (verbose)
	printf("Estimated:  %6.2f %6.4f nb %4d G %g\tGp %g\tX2 %g\tX2p %g\n",
	       ep[0], ep[1], nbins, G, Gp, X2, X2p);

      /* Output files, if requested.
       * (Best if ntrials=1. Will overwrite previous trials.)
       */
      if (ascfile != NULL)
	{
	  outfp = fopen(ascfile, "w");
	  esl_histogram_Print(outfp, h1);
	  fclose(outfp);
	}
      if (plotfile != NULL)
	{
	  outfp = fopen(plotfile, "w");
	  esl_histogram_Plot(outfp,  h1);
	  fclose(outfp);
	}
      if (survfile != NULL)  
	{
	  outfp = fopen(survfile, "w");
	  esl_histogram_PlotSurvival(outfp,  h1);
	  fclose(outfp);
	}
      if (qqfile != NULL)
	{
	  outfp = fopen(qqfile, "w");
	  if (fstrategy == FIT_TAIL)
	    esl_histogram_PlotQQ(outfp, h1, &esl_exp_generic_invcdf, ep);
	  else
	    esl_histogram_PlotQQ(outfp, h1, &esl_gumbel_generic_invcdf, ep);
	  fclose(outfp);
	}

      esl_histogram_Destroy(h);
      if (ntest > 0) esl_histogram_Destroy(h1);
    }

  /* Trap badness in an automated test.
   */
  if (fstrategy != FIT_TAIL && fabs(avg_ep[0] - p[0]) > 0.1)
    ESL_EXCEPTION(eslFAIL, "Something awry with Gumbel mu fit");
  if (fabs(avg_ep[1] - p[1]) > 0.1)
    ESL_EXCEPTION(eslFAIL, "Something awry with lambda fit");
 if (minGp < 1. / (1000. * ntrials))
    ESL_EXCEPTION(eslFAIL, "Something awry with G-test");
  if (minX2p < 1. / (1000. * ntrials))
    ESL_EXCEPTION(eslFAIL, "Something awry with chi squared test");

  /* Smaller final tests
   */
  if (! binmacro_test()) exit(1);
  
  esl_randomness_Destroy(r);
  return 0;
}
#endif /*eslHISTOGRAM_TESTDRIVE*/
