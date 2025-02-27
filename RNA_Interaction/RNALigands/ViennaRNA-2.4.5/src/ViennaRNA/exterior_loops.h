#ifndef VIENNA_RNA_PACKAGE_EXTERIOR_LOOPS_H
#define VIENNA_RNA_PACKAGE_EXTERIOR_LOOPS_H

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params.h>

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

/**
 *  @file     exterior_loops.h
 *  @ingroup  loops
 *  @brief    Energy evaluation of exterior loops for MFE and partition function calculations
 */

/**
 *  @addtogroup   loops
 *  @{
 *
 */

/**
 *  @brief  Auxiliary helper arrays for fast exterior loop computations
 *
 *  @see vrna_exp_E_ext_fast_init(), vrna_exp_E_ext_fast_rotate(),
 *  vrna_exp_E_ext_fast_free(), vrna_exp_E_ext_fast()
 */
typedef struct {
  FLT_OR_DBL  *qq;
  FLT_OR_DBL  *qq1;

  int         qqu_size;
  FLT_OR_DBL  **qqu;
} vrna_mx_pf_aux_el_t;


/**
 *  @brief  Evaluate a stem branching off the exterior loop
 *
 *  Given a base pair @f$(i,j)@f$ encoded by @em type, compute the energy contribution
 *  including dangling-end/terminal-mismatch contributions. Instead of returning the
 *  energy contribution per-se, this function returns the corresponding Boltzmann factor.
 *  If either of the adjacent nucleotides @f$(i - 1)@f$ and @f$(j+1)@f$ must not
 *  contribute stacking energy, the corresponding encoding must be @f$-1@f$.
 *
 *  @see vrna_E_exp_stem()
 *
 *  @param  type  The base pair encoding
 *  @param  n5d   The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1)
 *  @param  n3d   The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1)
 *  @param  p     The pre-computed energy parameters
 *  @return       The energy contribution of the introduced exterior-loop stem
 */
int
vrna_E_ext_stem(unsigned int  type,
                int           n5d,
                int           n3d,
                vrna_param_t  *p);


/**
 *  @brief  Evaluate a stem branching off the exterior loop (Boltzmann factor version)
 *
 *  Given a base pair @f$(i,j)@f$ encoded by @em type, compute the energy contribution
 *  including dangling-end/terminal-mismatch contributions. Instead of returning the
 *  energy contribution per-se, this function returns the corresponding Boltzmann factor.
 *  If either of the adjacent nucleotides @f$(i - 1)@f$ and @f$(j+1)@f$ must not
 *  contribute stacking energy, the corresponding encoding must be @f$-1@f$.
 *
 *  @see vrna_E_ext_stem()
 *
 *  @param  type  The base pair encoding
 *  @param  n5d   The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1)
 *  @param  n3d   The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1)
 *  @param  p     The pre-computed energy parameters (Boltzmann factor version)
 *  @return The Boltzmann weighted energy contribution of the introduced exterior-loop stem
 */
FLT_OR_DBL
vrna_exp_E_ext_stem(unsigned int      type,
                    int               n5d,
                    int               n3d,
                    vrna_exp_param_t  *p);


int
vrna_E_ext_loop(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j);


int
vrna_E_ext_loop_5(vrna_fold_compound_t *vc);


int
vrna_E_ext_loop_3(vrna_fold_compound_t  *fc,
                  int                   i);


int
vrna_BT_ext_loop_f5(vrna_fold_compound_t  *vc,
                    int                   *k,
                    int                   *i,
                    int                   *j,
                    vrna_bp_stack_t       *bp_stack,
                    int                   *stack_count);


int
vrna_BT_ext_loop_f3(vrna_fold_compound_t  *vc,
                    int                   *k,
                    int                   maxdist,
                    int                   *i,
                    int                   *j,
                    vrna_bp_stack_t       *bp_stack,
                    int                   *stack_count);


int
vrna_BT_ext_loop_f3_pp(vrna_fold_compound_t *fc,
                       int                  *i,
                       int                  maxdist);


vrna_mx_pf_aux_el_t *
vrna_exp_E_ext_fast_init(vrna_fold_compound_t *vc);


void
vrna_exp_E_ext_fast_rotate(vrna_fold_compound_t *vc,
                           vrna_mx_pf_aux_el_t  *aux_mx);


void
vrna_exp_E_ext_fast_free(vrna_fold_compound_t *vc,
                         vrna_mx_pf_aux_el_t  *aux_mx);


FLT_OR_DBL
vrna_exp_E_ext_fast(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    vrna_mx_pf_aux_el_t   *aux_mx);


void
vrna_exp_E_ext_fast_update(vrna_fold_compound_t *fc,
                           int                  j,
                           vrna_mx_pf_aux_el_t  *aux_mx);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @brief Compute the energy contribution of a stem branching off a loop-region
 *
 *  This function computes the energy contribution of a stem that branches off
 *  a loop region. This can be the case in multiloops, when a stem branching off
 *  increases the degree of the loop but also <I>immediately interior base pairs</I>
 *  of an exterior loop contribute free energy.
 *  To switch the behavior of the function according to the evaluation of a multiloop-
 *  or exterior-loop-stem, you pass the flag 'extLoop'.
 *  The returned energy contribution consists of a TerminalAU penalty if the pair type
 *  is greater than 2, dangling end contributions of mismatching nucleotides adjacent to
 *  the stem if only one of the si1, sj1 parameters is greater than 0 and mismatch energies
 *  if both mismatching nucleotides are positive values.
 *  Thus, to avoid incorporating dangling end or mismatch energies just pass a negative number,
 *  e.g. -1 to the mismatch argument.
 *
 *  This is an illustration of how the energy contribution is assembled:
 *  <PRE>
 *        3'  5'
 *        |   |
 *        X - Y
 *  5'-si1     sj1-3'
 *  </PRE>
 *
 *  Here, (X,Y) is the base pair that closes the stem that branches off a loop region.
 *  The nucleotides si1 and sj1 are the 5'- and 3'- mismatches, respectively. If the base pair
 *  type of (X,Y) is greater than 2 (i.e. an A-U or G-U pair, the TerminalAU penalty will be
 *  included in the energy contribution returned. If si1 and sj1 are both nonnegative numbers,
 *  mismatch energies will also be included. If one of si1 or sj1 is a negative value, only
 *  5' or 3' dangling end contributions are taken into account. To prohibit any of these mismatch
 *  contributions to be incorporated, just pass a negative number to both, si1 and sj1.
 *  In case the argument extLoop is 0, the returned energy contribution also includes
 *  the <I>internal-loop-penalty</I> of a multiloop stem with closing pair type.
 *
 *  @see    E_MLstem()
 *  @see    E_ExtLoop()
 *  @note   This function is threadsafe
 *
 *  @param  type    The pair type of the first base pair un the stem
 *  @param  si1     The 5'-mismatching nucleotide
 *  @param  sj1     The 3'-mismatching nucleotide
 *  @param  extLoop A flag that indicates whether the contribution reflects the one of an exterior loop or not
 *  @param  P       The data structure containing scaled energy parameters
 *  @return         The Free energy of the branch off the loop in dcal/mol
 *
 */
DEPRECATED(int E_Stem(int           type,
                      int           si1,
                      int           sj1,
                      int           extLoop,
                      vrna_param_t  *P),
           "This function is obsolete");


DEPRECATED(int E_ExtLoop(int          type,
                         int          si1,
                         int          sj1,
                         vrna_param_t *P),
           "Use vrna_E_ext_stem() instead");


/**
 *  This is the partition function variant of @ref E_ExtLoop()
 *  @deprecated Use vrna_exp_E_ext_stem() instead!
 *
 *  @see E_ExtLoop()
 *  @return The Boltzmann weighted energy contribution of the introduced exterior-loop stem
 */
DEPRECATED(FLT_OR_DBL exp_E_ExtLoop(int               type,
                                    int               si1,
                                    int               sj1,
                                    vrna_exp_param_t  *P),
           "Use vrna_exp_E_ext_stem() instead");


/**
 *  <H2>Compute the Boltzmann weighted energy contribution of a stem branching off a loop-region</H2>
 *  This is the partition function variant of @ref E_Stem()
 *  @see E_Stem()
 *  @note This function is threadsafe
 *
 *  @return The Boltzmann weighted energy contribution of the branch off the loop
 */
DEPRECATED(FLT_OR_DBL exp_E_Stem(int              type,
                                 int              si1,
                                 int              sj1,
                                 int              extLoop,
                                 vrna_exp_param_t *P),
           "This function is obsolete");


#endif

/**
 * @}
 */


#endif
