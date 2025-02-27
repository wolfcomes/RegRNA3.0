#ifndef VIENNA_RNA_PACKAGE_LFOLD_H
#define VIENNA_RNA_PACKAGE_LFOLD_H

/**
 *  @file Lfold.h
 *  @ingroup  local_fold
 *  @brief    Functions for locally optimal MFE structure prediction
 */

#ifdef VRNA_WARN_DEPRECATED
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

#include <ViennaRNA/mfe.h>

/**
 *  @brief Local MFE prediction using a sliding window approach (simplified interface)
 *
 *  This simplified interface to vrna_mfe_window() computes the MFE and locally
 *  optimal secondary structure using default options. Structures are predicted
 *  using a sliding window approach, where base pairs may not span outside the
 *  window. Memory required for dynamic programming (DP) matrices will
 *  be allocated and free'd on-the-fly. Hence, after return of this function, the recursively filled
 *  matrices are not available any more for any post-processing.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *  you require other conditions than specified by the default model details, use vrna_mfe_window(),
 *  and the data structure #vrna_fold_compound_t instead.
 *
 *  @ingroup local_mfe_fold
 *
 *  @see  vrna_mfe_window(), vrna_Lfoldz(), vrna_mfe_window_zscore(), vrna_fold_compound(),
 *        #vrna_fold_compound_t
 *
 *  @param  string      The nucleic acid sequence
 *  @param  window_size The window size for locally optimal structures
 *  @param  file        The output file handle where predictions are written to (if NULL, output is written to stdout)
 */
float
vrna_Lfold(const char *string,
           int        window_size,
           FILE       *file);


float
vrna_Lfold_cb(const char                *string,
              int                       window_size,
              vrna_mfe_window_callback  *cb,
              void                      *data);


#ifdef VRNA_WITH_SVM
/**
 *  @brief Local MFE prediction using a sliding window approach with z-score cut-off (simplified interface)
 *
 *  This simplified interface to vrna_mfe_window_zscore() computes the MFE and locally
 *  optimal secondary structure using default options. Structures are predicted
 *  using a sliding window approach, where base pairs may not span outside the
 *  window. Memory required for dynamic programming (DP) matrices will
 *  be allocated and free'd on-the-fly. Hence, after return of this function, the recursively filled
 *  matrices are not available any more for any post-processing.
 *  This function is the z-score version of vrna_Lfold(), i.e.
 *  only predictions above a certain z-score cut-off value are
 *  printed.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *  you require other conditions than specified by the default model details, use vrna_mfe_window(),
 *  and the data structure #vrna_fold_compound_t instead.
 *
 *  @ingroup local_mfe_fold
 *
 *  @see  vrna_mfe_window_zscore(), vrna_Lfold(), vrna_mfe_window(), vrna_fold_compound(),
 *        #vrna_fold_compound_t
 *
 *  @param  string      The nucleic acid sequence
 *  @param  window_size The window size for locally optimal structures
 *  @param  min_z       The minimal z-score for a predicted structure to appear in the output
 *  @param  file        The output file handle where predictions are written to (if NULL, output is written to stdout)
 */
float
vrna_Lfoldz(const char  *string,
            int         window_size,
            double      min_z,
            FILE        *file);


float
vrna_Lfoldz_cb(const char                       *string,
               int                              window_size,
               double                           min_z,
               vrna_mfe_window_zscore_callback  *cb,
               void                             *data);


#endif


/**
 *  @addtogroup local_consensus_fold
 *  @{
 *
 *  @}
 */

float vrna_aliLfold(const char  **AS,
                    int         maxdist,
                    FILE        *fp);


float vrna_aliLfold_cb(const char               **AS,
                       int                      maxdist,
                       vrna_mfe_window_callback *cb,
                       void                     *data);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief The local analog to fold().
 *
 *  Computes the minimum free energy structure including only base pairs
 *  with a span smaller than 'maxdist'
 *
 *  @ingroup local_mfe_fold
 *
 *  @deprecated Use vrna_mfe_window() instead!
 */
DEPRECATED(float Lfold(const char *string,
                       char       *structure,
                       int        maxdist),
"Use vrna_Lfold() or vrna_Lfold_cb() instead");

/**
 *  @brief
 *
 *  @ingroup local_mfe_fold
 *
 *  @deprecated Use vrna_mfe_window_zscore() instead!
 */
DEPRECATED(float Lfoldz(const char  *string,
                        char        *structure,
                        int         maxdist,
                        int         zsc,
                        double      min_z),
"Use vrna_Lfoldz() or vrna_Lfoldz_cb() instead");

/**
 *  @brief
 *
 *  @ingroup local_consensus_fold
 *
 *  @param strings
 *  @param structure
 *  @param maxdist
 *  @return
 */
DEPRECATED(float aliLfold(const char **AS,
               char       *structure,
               int        maxdist),
"Use vrna_aliLfold() or vrna_aliLfold_cb() instead");


DEPRECATED(float aliLfold_cb(const char                **AS,
                  int                       maxdist,
                  vrna_mfe_window_callback  *cb,
                  void                      *data),
"Use vrna_aliLfold() or vrna_aliLfold_cb() instead");


#endif

#endif
