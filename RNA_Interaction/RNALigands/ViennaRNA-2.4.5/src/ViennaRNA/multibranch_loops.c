#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/exterior_loops.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/multibranch_loops.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "multibranch_loops.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
E_mb_loop_fast(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               int                  *dmli1,
               int                  *dmli2);


PRIVATE int
E_mb_loop_fast_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j,
                      int                   *dmli1,
                      int                   *dmli2);


PRIVATE int
E_mb_loop_fast_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j,
                           int                  *dmli1,
                           int                  *dmli2);


PRIVATE int
E_mb_loop_fast_comparative_window(vrna_fold_compound_t  *vc,
                                  int                   i,
                                  int                   j,
                                  int                   *dmli1,
                                  int                   *dmli2);


PRIVATE int
E_ml_stems_fast(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j,
                int                   *fmi,
                int                   *dmli);


PRIVATE int
E_ml_stems_fast_comparative(vrna_fold_compound_t  *vc,
                            int                   i,
                            int                   j,
                            int                   *fmi,
                            int                   *dmli);


PRIVATE int
E_ml_stems_fast_comparative_window(vrna_fold_compound_t *vc,
                                   int                  i,
                                   int                  j,
                                   int                  *fmi,
                                   int                  *dmli);


PRIVATE int
E_ml_stems_fast_window(vrna_fold_compound_t *vc,
                       int                  i,
                       int                  j,
                       int                  *fmi,
                       int                  *dmli);


PRIVATE int
extend_fm_3p(int                  i,
             int                  j,
             int                  *fm,
             vrna_fold_compound_t *vc);


PRIVATE int
extend_fm_3p_comparative(int                  i,
                         int                  j,
                         int                  *fm,
                         vrna_fold_compound_t *vc);


PRIVATE int
extend_fm_3p_window(int                   i,
                    int                   j,
                    int                   **fm,
                    vrna_fold_compound_t  *vc);


PRIVATE int E_mb_loop_stack(vrna_fold_compound_t  *vc,
                            int                   i,
                            int                   j);


PRIVATE int E_mb_loop_stack_window(vrna_fold_compound_t *vc,
                                   int                  i,
                                   int                  j);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_E_mb_loop_fast(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    int                   *dmli1,
                    int                   *dmli2)
{
  int e = INF;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_mb_loop_fast_window(vc, i, j, dmli1, dmli2);
        else
          e = E_mb_loop_fast(vc, i, j, dmli1, dmli2);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_mb_loop_fast_comparative_window(vc, i, j, dmli1, dmli2);
        else
          e = E_mb_loop_fast_comparative(vc, i, j, dmli1, dmli2);

        break;
    }
  }

  return e;
}


PUBLIC int
vrna_E_ml_stems_fast(vrna_fold_compound_t *vc,
                     int                  i,
                     int                  j,
                     int                  *fmi,
                     int                  *dmli)
{
  int e = INF;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_ml_stems_fast_window(vc, i, j, fmi, dmli);
        else
          e = E_ml_stems_fast(vc, i, j, fmi, dmli);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_ml_stems_fast_comparative_window(vc, i, j, fmi, dmli);
        else
          e = E_ml_stems_fast_comparative(vc, i, j, fmi, dmli);

        break;
    }
  }

  return e;
}


PUBLIC int
vrna_E_mb_loop_stack(vrna_fold_compound_t *vc,
                     int                  i,
                     int                  j)
{
  int e = INF;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_mb_loop_stack_window(vc, i, j);
        else
          e = E_mb_loop_stack(vc, i, j);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        break;
    }
  }

  return e;
}


PRIVATE int
E_mb_loop_fast_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j,
                           int                  *dmli1,
                           int                  *dmli2)
{
  short                     **S, **S5, **S3;
  int                       *indx, e, decomp, s, n_seq, dangle_model, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq         = vc->n_seq;
  indx          = vc->jindx;
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  scs           = vc->scs;
  dangle_model  = md->dangles;
  e             = INF;
  evaluate      = prepare_hc_default(vc, &hc_dat_local);

  /* multi-loop decomposition ------------------------*/
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    decomp = dmli1[j - 1];

    S   = vc->S;
    S5  = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
    S3  = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */

    if (dangle_model) {
      for (s = 0; s < n_seq; s++) {
        tt      = vrna_get_ptype_md(S[s][j], S[s][i], md);
        decomp  += E_MLstem(tt, S5[s][j], S3[s][i], P);
      }
    } else {
      for (s = 0; s < n_seq; s++) {
        tt      = vrna_get_ptype_md(S[s][j], S[s][i], md);
        decomp  += E_MLstem(tt, -1, -1, P);
      }
    }

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->energy_bp)
            decomp += scs[s]->energy_bp[indx[j] + i];
      }
    }

    e = decomp + n_seq * P->MLclosing;
  }

  return e;
}


PRIVATE int
E_mb_loop_fast_comparative_window(vrna_fold_compound_t  *vc,
                                  int                   i,
                                  int                   j,
                                  int                   *dmli1,
                                  int                   *dmli2)
{
  short                     **S, **S5, **S3;
  int                       e, decomp, s, n_seq, dangle_model, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq         = vc->n_seq;
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  scs           = vc->scs;
  dangle_model  = md->dangles;
  e             = INF;
  evaluate      = prepare_hc_default_window(vc, &hc_dat_local);

  /* multi-loop decomposition ------------------------*/
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    decomp = dmli1[j - 1 - (i + 1)];

    S   = vc->S;
    S5  = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
    S3  = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */

    if (dangle_model) {
      for (s = 0; s < n_seq; s++) {
        tt      = vrna_get_ptype_md(S[s][j], S[s][i], md);
        decomp  += E_MLstem(tt, S5[s][j], S3[s][i], P);
      }
    } else {
      for (s = 0; s < n_seq; s++) {
        tt      = vrna_get_ptype_md(S[s][j], S[s][i], md);
        decomp  += E_MLstem(tt, -1, -1, P);
      }
    }

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->f)
            decomp += scs[s]->f(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, scs[s]->data);
      }
    }

    e = decomp + n_seq * P->MLclosing;
  }

  return e;
}


PRIVATE int
E_mb_loop_fast(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               int                  *dmli1,
               int                  *dmli2)
{
  short                     S_i1, S_j1, *S, *S2;
  unsigned int              strands, *sn, *so, *ss, *se;
  int                       decomp, en, e, *indx, *fc, ij, dangle_model, tt, noGUclosure;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  S             = vc->sequence_encoding;
  S2            = vc->sequence_encoding2;
  indx          = vc->jindx;
  strands       = vc->strands;
  sn            = vc->strand_number;
  so            = vc->strand_order;
  ss            = vc->strand_start;
  se            = vc->strand_end;
  hc            = vc->hc;
  sc            = vc->sc;
  fc            = vc->matrices->fc;
  P             = vc->params;
  md            = &(P->model_details);
  noGUclosure   = md->noGUclosure;
  dangle_model  = md->dangles;

  /* init values */
  e       = INF;
  decomp  = INF;

  evaluate = prepare_hc_default(vc, &hc_dat_local);

  ij  = indx[j] + i;
  tt  = vrna_get_ptype_md(S2[j], S2[i], md);

  if (noGUclosure && ((tt == 3) || (tt == 4)))
    return e;

  if (strands == 1) {
    S_i1  = S[i + 1];
    S_j1  = S[j - 1];
  } else {
    S_i1  = (sn[i] == sn[i + 1]) ? S[i + 1] : -1;
    S_j1  = (sn[j - 1] == sn[j]) ? S[j - 1] : -1;
  }

  if ((S_i1 >= 0) && (S_j1 >= 0)) {
    /* regular multi branch loop */
    /* new closing pair (i,j) with mb part [i+1,j-1] */
    if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      decomp = dmli1[j - 1];

      if (decomp != INF) {
        switch (dangle_model) {
          /* no dangles */
          case 0:
            decomp += E_MLstem(tt, -1, -1, P);
            if (sc) {
              if (sc->energy_bp)
                decomp += sc->energy_bp[ij];

              if (sc->f)
                decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
            }

            break;

          /* double dangles */
          case 2:
            decomp += E_MLstem(tt, S_j1, S_i1, P);
            if (sc) {
              if (sc->energy_bp)
                decomp += sc->energy_bp[ij];

              if (sc->f)
                decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
            }

            break;

          /* normal dangles, aka dangles = 1 || 3 */
          default:
            decomp += E_MLstem(tt, -1, -1, P);
            if (sc) {
              if (sc->energy_bp)
                decomp += sc->energy_bp[ij];

              if (sc->f)
                decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
            }

            break;
        }
      }
    }

    if (dangle_model % 2) {
      /* dangles == 1 || dangles == 3 */
      /* new closing pair (i,j) with mb part [i+2,j-1] */
      if (evaluate(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if (dmli2[j - 1] != INF) {
          en = dmli2[j - 1] +
               E_MLstem(tt, -1, S_i1, P) +
               P->MLbase;

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[i + 1][1];

            if (sc->energy_bp)
              en += sc->energy_bp[ij];

            if (sc->f)
              en += sc->f(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
          }

          decomp = MIN2(decomp, en);
        }
      }

      /* new closing pair (i,j) with mb part [i+2.j-2] */
      if (evaluate(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if (dmli2[j - 2] != INF) {
          en = dmli2[j - 2] +
               E_MLstem(tt, S_j1, S_i1, P) +
               2 * P->MLbase;

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[i + 1][1] +
                    sc->energy_up[j - 1][1];

            if (sc->energy_bp)
              en += sc->energy_bp[ij];

            if (sc->f)
              en += sc->f(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, sc->data);
          }

          decomp = MIN2(decomp, en);
        }
      }

      /* new closing pair (i,j) with mb part [i+1, j-2] */
      if (evaluate(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if (dmli1[j - 2] != INF) {
          en = dmli1[j - 2] +
               E_MLstem(tt, S_j1, -1, P) +
               P->MLbase;
          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[j - 1][1];

            if (sc->energy_bp)
              en += sc->energy_bp[ij];

            if (sc->f)
              en += sc->f(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, sc->data);
          }

          decomp = MIN2(decomp, en);
        }
      }
    } /* end if dangles % 2 */

    if (decomp != INF)
      e = decomp + P->MLclosing;
  } /* end regular multibranch loop */

  if (sn[i] != sn[j]) {
    /* multibrach like cofold structure with cut somewhere between i and j */
    if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      if ((fc[i + 1] != INF) && (fc[j - 1] != INF)) {
        decomp = fc[i + 1] +
                 fc[j - 1];
        switch (dangle_model) {
          case 0:
            decomp += E_ExtLoop(tt, -1, -1, P);
            break;

          case 2:
            decomp += E_ExtLoop(tt, S_j1, S_i1, P);
            break;

          default:
            decomp += E_ExtLoop(tt, -1, -1, P);
            break;
        }
      }
    }

    if (dangle_model % 2) {
      /* dangles == 1 || dangles == 3 */
      if (evaluate(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if ((fc[i + 2] != INF) && (fc[j - 1] != INF)) {
          en = fc[i + 2] +
               fc[j - 1] +
               E_ExtLoop(tt, -1, S_i1, P);
          decomp = MIN2(decomp, en);
        }
      }

      if (evaluate(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if ((fc[i + 1] != INF) && (fc[j - 2] != INF)) {
          en = fc[i + 1] +
               fc[j - 2] +
               E_ExtLoop(tt, S_j1, -1, P);
          decomp = MIN2(decomp, en);
        }
      }

      if (evaluate(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if ((fc[i + 2] != INF) && (fc[j - 2] != INF)) {
          en = fc[i + 2] +
               fc[j - 2] +
               E_ExtLoop(tt, S_j1, S_i1, P);
          decomp = MIN2(decomp, en);
        }
      }
    }

    e = MIN2(e, decomp);
  }

  return e;
}


PRIVATE int
E_mb_loop_fast_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j,
                      int                   *dmli1,
                      int                   *dmli2)
{
  short                     S_i1, S_j1, *S, *S2;
  int                       decomp, en, e, dangle_model, tt, noGUclosure;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  S             = vc->sequence_encoding;
  S2            = vc->sequence_encoding2;
  hc            = vc->hc;
  sc            = vc->sc;
  P             = vc->params;
  md            = &(P->model_details);
  noGUclosure   = md->noGUclosure;
  dangle_model  = md->dangles;

  /* init values */
  e       = INF;
  decomp  = INF;

  evaluate = prepare_hc_default_window(vc, &hc_dat_local);

  S_i1  = S[i + 1];
  S_j1  = S[j - 1];

  tt = vrna_get_ptype_md(S2[j], S2[i], md);

  if (noGUclosure && ((tt == 3) || (tt == 4)))
    return e;

  /* new closing pair (i,j) with mb part [i+1,j-1] */
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    decomp = dmli1[j - 1 - (i + 1)];

    if (decomp < INF) {
      switch (dangle_model) {
        /* no dangles */
        case 0:
          decomp += E_MLstem(tt, -1, -1, P);
          break;
        /* double dangles */
        case 2:
          decomp += E_MLstem(tt, S_j1, S_i1, P);
          break;
        /* normal dangles, aka dangle_model = 1 */
        default:
          decomp += E_MLstem(tt, -1, -1, P);
          break;
      }
      if (sc) {
        if (sc->energy_bp_local)
          decomp += sc->energy_bp_local[i][j - i];

        if (sc->f)
          decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
      }
    }
  }

  if (dangle_model % 2) {
    /* dangles == 1 || dangles == 3 */
    /* new closing pair (i,j) with mb part [i+2,j-1] */
    if (evaluate(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      if (dmli2[j - 1 - (i + 2)] != INF) {
        en = dmli2[j - 1 - (i + 2)] +
             E_MLstem(tt, -1, S_i1, P) +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            decomp += sc->energy_up[i + 1][1];

          if (sc->energy_bp_local)
            decomp += sc->energy_bp_local[i][j - i];

          if (sc->f)
            en += sc->f(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        decomp = MIN2(decomp, en);
      }
    }

    /* new closing pair (i,j) with mb part [i+2.j-2] */
    if (evaluate(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      if (dmli2[j - 2 - (i + 2)] != INF) {
        en = dmli2[j - 2 - (i + 2)] +
             E_MLstem(tt, S_j1, S_i1, P) +
             2 * P->MLbase;

        if (sc) {
          if (sc->energy_up)
            decomp += sc->energy_up[i + 1][1] +
                      sc->energy_up[j - 1][1];

          if (sc->energy_bp_local)
            decomp += sc->energy_bp_local[i][j - i];

          if (sc->f)
            en += sc->f(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        decomp = MIN2(decomp, en);
      }
    }

    /* new closing pair (i,j) with mb part [i+1, j-2] */
    if (evaluate(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      if (dmli1[j - 2 - (i + 1)] != INF) {
        en = dmli1[j - 2 - (i + 1)] +
             E_MLstem(tt, S_j1, -1, P) +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            decomp += sc->energy_up[j - 1][1];

          if (sc->energy_bp_local)
            decomp += sc->energy_bp_local[i][j - i];

          if (sc->f)
            en += sc->f(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        decomp = MIN2(decomp, en);
      }
    }
  } /* end if dangles % 2 */

  if (decomp != INF)
    e = decomp + P->MLclosing;

  return e;
}


PRIVATE int
E_mb_loop_stack(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j)
{
  char                      *ptype;
  int                       *c, *fML, e, decomp, en, i1k, k1j1, ij, k, *indx, turn,
                            type, type_2, *rtype;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  indx  = vc->jindx;
  hc    = vc->hc;
  P     = vc->params;
  md    = &(P->model_details);
  turn  = md->min_loop_size;
  rtype = &(md->rtype[0]);
  c     = vc->matrices->c;
  fML   = vc->matrices->fML;
  ptype = vc->ptype;
  ij    = indx[j] + i;
  type  = vrna_get_ptype(ij, ptype);
  sc    = vc->sc;
  e     = INF;

  evaluate = prepare_hc_default(vc, &hc_dat_local);

  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    decomp  = INF;
    k1j1    = indx[j - 1] + i + 2 + turn + 1;
    for (k = i + 2 + turn; k < j - 2 - turn; k++, k1j1++) {
      i1k = indx[k] + i + 1;

      if (evaluate(i, j, i + 1, k, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
        type_2 = rtype[vrna_get_ptype(i1k, ptype)];

        en = c[i1k] +
             P->stack[type][type_2] +
             fML[k1j1];

        if (sc)
          if (sc->f)
            en += sc->f(i, j, i + 1, k, VRNA_DECOMP_ML_COAXIAL, sc->data);

        decomp = MIN2(decomp, en);
      }

      if (evaluate(i, j, k + 1, j - 1, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
        type_2 = rtype[vrna_get_ptype(k1j1, ptype)];

        en = c[k1j1] +
             P->stack[type][type_2] +
             fML[i1k];

        if (sc)
          if (sc->f)
            en += sc->f(i, j, k + 1, j - 1, VRNA_DECOMP_ML_COAXIAL, sc->data);

        decomp = MIN2(decomp, en);
      }
    }
    /* no TermAU penalty if coax stack */
    decomp += 2 * P->MLintern[1] + P->MLclosing;
    if (sc) {
      if (sc->energy_bp)
        decomp += sc->energy_bp[ij];

      if (sc->f)
        decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
    }

    e = decomp;
  }

  return e;
}


PRIVATE int
E_mb_loop_stack_window(vrna_fold_compound_t *vc,
                       int                  i,
                       int                  j)
{
  char                      **ptype;
  int                       **c, **fML, e, decomp, en, k, turn, *rtype, type, type_2;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  hc        = vc->hc;
  P         = vc->params;
  md        = &(P->model_details);
  turn      = md->min_loop_size;
  rtype     = &(md->rtype[0]);
  c         = vc->matrices->c_local;
  fML       = vc->matrices->fML_local;
  sc        = vc->sc;
  e         = INF;
  evaluate  = prepare_hc_default_window(vc, &hc_dat_local);

  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    ptype = vc->ptype_local;
    type  = vrna_get_ptype_window(i, j, ptype);

    decomp = INF;
    for (k = i + 2 + turn; k < j - 2 - turn; k++) {
      if (evaluate(i, j, i + 1, k, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
        type_2 = rtype[vrna_get_ptype_window(i + 1, k, ptype)];

        en = c[i + 1][k - i - 1] +
             P->stack[type][type_2] +
             fML[k + 1][j - 1 - k - 1];

        if (sc)
          if (sc->f)
            en += sc->f(i, j, i + 1, k, VRNA_DECOMP_ML_COAXIAL, sc->data);

        decomp = MIN2(decomp, en);
      }

      if (evaluate(i, j, k + 1, j - 1, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
        type_2 = rtype[vrna_get_ptype_window(k + 1, j - 1, ptype)];

        en = c[k + 1][j - 1 - k - 1] +
             P->stack[type][type_2] +
             fML[i + 1][k - i - 1];

        if (sc)
          if (sc->f)
            en += sc->f(i, j, k + 1, j - 1, VRNA_DECOMP_ML_COAXIAL, sc->data);

        decomp = MIN2(decomp, en);
      }
    }
    /* no TermAU penalty if coax stack */
    decomp += 2 * P->MLintern[1] + P->MLclosing;
    if (sc) {
      if (sc->energy_bp_local)
        decomp += sc->energy_bp_local[i][j - i];

      if (sc->f)
        decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
    }

    e = decomp;
  }

  return e;
}


PUBLIC int
E_ml_rightmost_stem(int                   i,
                    int                   j,
                    vrna_fold_compound_t  *vc)
{
  if ((vc) && (vc->matrices) && (vc->matrices->fM1)) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        return extend_fm_3p(i, j, vc->matrices->fM1, vc);

      case VRNA_FC_TYPE_COMPARATIVE:
        return extend_fm_3p_comparative(i, j, vc->matrices->fM1, vc);
    }
  }

  return INF;
}


/*
 * compose a multibranch loop part fm[i:j]
 * by either c[i,j]/ggg[i,j] or fm[i:j-1]
 *
 * This function can be used for fM and fM1
 */
PRIVATE int
extend_fm_3p(int                  i,
             int                  j,
             int                  *fm,
             vrna_fold_compound_t *vc)
{
  short                     *S;
  unsigned int              *sn;
  int                       en, length, *indx, *c, *ggg, ij, type,
                            dangle_model, with_gquad, e, u, k, cnt, with_ud;
  vrna_param_t              *P;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  P             = vc->params;
  length        = vc->length;
  S             = vc->sequence_encoding;
  indx          = vc->jindx;
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  c             = vc->matrices->c;
  ggg           = vc->matrices->ggg;
  ij            = indx[j] + i;
  type          = vrna_get_ptype(ij, vc->ptype);
  dangle_model  = P->model_details.dangles;
  with_gquad    = P->model_details.gquad;
  domains_up    = vc->domains_up;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  e             = INF;
  evaluate      = prepare_hc_default(vc, &hc_dat_local);

  if (sn[i - 1] == sn[i]) {
    if (sn[j] == sn[j + 1]) {
      if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        e = c[ij];
        if (e != INF) {
          switch (dangle_model) {
            case 2:
              e += E_MLstem(type, (i == 1) ? S[length] : S[i - 1], S[j + 1], P);
              break;

            default:
              e += E_MLstem(type, -1, -1, P);
              break;
          }
          if (sc)
            if (sc->f)
              e += sc->f(i, j, i, j, VRNA_DECOMP_ML_STEM, sc->data);
        }
      }

      if (with_gquad) {
        if (sn[i] == sn[j]) {
          en  = ggg[ij] + E_MLstem(0, -1, -1, P);
          e   = MIN2(e, en);
        }
      }
    }

    if (sn[j - 1] == sn[j]) {
      if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        if (fm[indx[j - 1] + i] != INF) {
          en = fm[indx[j - 1] + i] +
               P->MLbase;

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[j][1];

            if (sc->f)
              en += sc->f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
          }

          e = MIN2(e, en);
        }
      }
    }

    if (with_ud) {
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        k = j - u + 1;
        if ((k > i) && (sn[j - u] == sn[j])) {
          if (evaluate(i, j, i, k - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
            if (fm[indx[k - 1] + i] != INF) {
              en = domains_up->energy_cb(vc,
                                         k,
                                         j,
                                         VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                         domains_up->data);
              if (en != INF) {
                en += fm[indx[k - 1] + i] +
                      u * P->MLbase;

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[k][u];

                  if (sc->f)
                    en += sc->f(i, j, i, k - 1, VRNA_DECOMP_ML_ML, sc->data);
                }

                e = MIN2(e, en);
              }
            }
          }
        }
      }
    }
  }

  return e;
}


PRIVATE int
extend_fm_3p_comparative(int                  i,
                         int                  j,
                         int                  *fm,
                         vrna_fold_compound_t *vc)
{
  short                     **S, **S5, **S3;
  unsigned int              **a2s, *sn;
  int                       en, length, *indx, *c, *ggg, ij, type, n_seq, s,
                            dangle_model, with_gquad, e, u, k, cnt;
  vrna_md_t                 *md;
  vrna_param_t              *P;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq         = vc->n_seq;
  P             = vc->params;
  md            = &(P->model_details);
  length        = vc->length;
  a2s           = vc->a2s;
  S             = vc->S;
  S5            = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
  S3            = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */
  indx          = vc->jindx;
  sn            = vc->strand_number;
  hc            = vc->hc;
  scs           = vc->scs;
  c             = vc->matrices->c;
  ggg           = vc->matrices->ggg;
  ij            = indx[j] + i;
  dangle_model  = md->dangles;
  with_gquad    = md->gquad;
  e             = INF;
  evaluate      = prepare_hc_default(vc, &hc_dat_local);

  if (sn[i - 1] == sn[i]) {
    if (sn[j] == sn[j + 1]) {
      if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        en = c[ij];
        if (en != INF) {
          switch (dangle_model) {
            case 2:
              for (s = 0; s < n_seq; s++) {
                type  = vrna_get_ptype_md(S[s][i], S[s][j], md);
                en    += E_MLstem(type, S5[s][i], S3[s][j], P);
              }
              break;

            default:
              for (s = 0; s < n_seq; s++) {
                type  = vrna_get_ptype_md(S[s][i], S[s][j], md);
                en    += E_MLstem(type, -1, -1, P);
              }
              break;
          }

          if (scs) {
            for (s = 0; s < n_seq; s++)
              if (scs[s] && scs[s]->f)
                en += scs[s]->f(i, j, i, j, VRNA_DECOMP_ML_STEM, scs[s]->data);
          }

          e = MIN2(e, en);
        }
      }

      if (with_gquad) {
        if (sn[i] == sn[j]) {
          en = ggg[ij] +
               n_seq * E_MLstem(0, -1, -1, P);
          e = MIN2(e, en);
        }
      }
    }

    if (sn[j - 1] == sn[j]) {
      if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        if (fm[indx[j - 1] + i] != INF) {
          en = fm[indx[j - 1] + i] +
               n_seq * P->MLbase;

          if (scs) {
            for (s = 0; s < n_seq; s++) {
              if (scs[s]) {
                if (scs[s]->energy_up)
                  en += scs[s]->energy_up[a2s[s][j]][1];

                if (scs[s]->f)
                  en += scs[s]->f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, scs[s]->data);
              }
            }
          }

          e = MIN2(e, en);
        }
      }
    }
  }

  return e;
}


#ifdef VRNA_WITH_SSE_IMPLEMENTATION
/* SSE modular decomposition -------------------------------*/
#include <emmintrin.h>
#include <smmintrin.h>

//http://stackoverflow.com/questions/9877700/getting-max-value-in-a-m128i-vector-with-sse
int
horizontal_min_Vec4i(__m128i x)
{
  __m128i min1  = _mm_shuffle_epi32(x, _MM_SHUFFLE(0, 0, 3, 2));
  __m128i min2  = _mm_min_epi32(x, min1);
  __m128i min3  = _mm_shuffle_epi32(min2, _MM_SHUFFLE(0, 0, 0, 1));
  __m128i min4  = _mm_min_epi32(min2, min3);

  return _mm_cvtsi128_si32(min4);
}


PRIVATE int
modular_decomposition(const int i,
                      const int ij,
                      const int j,
                      const int turn,
                      const int *fmi,
                      const int *fm)
{
  int       k       = i + turn + 1;
  int       k1j     = ij + turn + 2; //indx[j] + i + 1; //indx[j] + i + turn + 2;
  const int stop    = j - 2 - turn;
  int       decomp  = INF;
  {
    const int end = 1 + stop - k;
    int       i;
    __m128i   inf = _mm_set1_epi32(INF);

    for (i = 0; i < end - 3; i += 4) {
      __m128i   a = _mm_loadu_si128((__m128i *)&fmi[k + i]);
      __m128i   b = _mm_loadu_si128((__m128i *)&fm[k1j + i]);
      __m128i   c = _mm_add_epi32(a, b);
      /* deactivate this part if you are sure to not use any hard constraints */
#if 1
      __m128i   mask1 = _mm_cmplt_epi32(a, inf);
      __m128i   mask2 = _mm_cmplt_epi32(b, inf);
      __m128i   res   = _mm_or_si128(_mm_and_si128(mask1, c),
                                     _mm_andnot_si128(mask1, a));

      res = _mm_or_si128(_mm_and_si128(mask2, res),
                         _mm_andnot_si128(mask2, b));
      const int en = horizontal_min_Vec4i(res);
#else
      const int en = horizontal_min_Vec4i(c);
#endif
      decomp = MIN2(decomp, en);
    }
    for (; i < end; i++) {
      if ((fmi[k + i] != INF) && (fm[k1j + i] != INF)) {
        const int en = fmi[k + i] + fm[k1j + i];
        decomp = MIN2(decomp, en);
      }
    }
  }

  return decomp;
}


/* End SSE modular decomposition -------------------------------*/
#endif


PRIVATE int
extend_fm_3p_window(int                   i,
                    int                   j,
                    int                   **fm,
                    vrna_fold_compound_t  *vc)
{
  short                     *S;
  unsigned int              *sn;
  int                       en, length, **c, **ggg, type,
                            dangle_model, with_gquad, e, u, k, cnt, with_ud;
  vrna_param_t              *P;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  P             = vc->params;
  length        = vc->length;
  S             = vc->sequence_encoding;
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  c             = vc->matrices->c_local;
  ggg           = vc->matrices->ggg_local;
  dangle_model  = P->model_details.dangles;
  with_gquad    = P->model_details.gquad;
  domains_up    = vc->domains_up;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  e             = INF;


  type = vrna_get_ptype_window(i, j, vc->ptype_local);

  evaluate = prepare_hc_default_window(vc, &hc_dat_local);

  if (sn[i - 1] == sn[i]) {
    if (sn[j] == sn[j + 1]) {
      if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        e = c[i][j - i];
        if (e != INF) {
          switch (dangle_model) {
            case 2:
              e += E_MLstem(type, (i == 1) ? S[length] : S[i - 1], S[j + 1], P);
              break;

            default:
              e += E_MLstem(type, -1, -1, P);
              break;
          }
          if (sc)
            if (sc->f)
              e += sc->f(i, j, i, j, VRNA_DECOMP_ML_STEM, sc->data);
        }
      }

      if (with_gquad) {
        if (sn[i] == sn[j]) {
          en = ggg[i][j - i] +
               E_MLstem(0, -1, -1, P);
          e = MIN2(e, en);
        }
      }
    }

    if (sn[j - 1] == sn[j]) {
      if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        if (fm[i][j - 1 - i] != INF) {
          en = fm[i][j - 1 - i] +
               P->MLbase;

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[j][1];

            if (sc->f)
              en += sc->f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
          }

          e = MIN2(e, en);
        }
      }
    }

    if (with_ud) {
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        k = j - u + 1;
        if ((k > i) && (sn[j - u] == sn[j])) {
          if (evaluate(i, j, i, k - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
            if (fm[i][k - 1 - i] != INF) {
              en = domains_up->energy_cb(vc,
                                         k,
                                         j,
                                         VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                         domains_up->data);
              if (en != INF) {
                en += fm[i][k - 1 - i] +
                      u * P->MLbase;

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[k][u];

                  if (sc->f)
                    en += sc->f(i, j, i, k - 1, VRNA_DECOMP_ML_ML, sc->data);
                }

                e = MIN2(e, en);
              }
            }
          }
        }
      }
    }
  }

  return e;
}


PRIVATE int
E_ml_stems_fast(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j,
                int                   *fmi,
                int                   *dmli)
{
  char                      *ptype;
  short                     *S;
  unsigned int              strands, *sn, *so, *ss, *se;
  int                       k, en, decomp, mm5, mm3, type_2, k1j, stop, length, *indx,
                            *c, *fm, ij, dangle_model, turn, type, *rtype, circular, e, u,
                            cnt, with_ud;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_param_t              *P;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = (int)vc->length;
  ptype         = vc->ptype;
  S             = vc->sequence_encoding;
  indx          = vc->jindx;
  strands       = vc->strands;
  sn            = vc->strand_number;
  so            = vc->strand_order;
  ss            = vc->strand_start;
  se            = vc->strand_end;
  hc            = vc->hc;
  sc            = vc->sc;
  c             = vc->matrices->c;
  fm            = vc->matrices->fML;
  P             = vc->params;
  ij            = indx[j] + i;
  dangle_model  = P->model_details.dangles;
  turn          = P->model_details.min_loop_size;
  type          = vrna_get_ptype(ij, ptype);
  rtype         = &(P->model_details.rtype[0]);
  circular      = P->model_details.circ;
  domains_up    = vc->domains_up;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  e             = INF;
  evaluate      = prepare_hc_default(vc, &hc_dat_local);

  /*
   *  extension with one unpaired nucleotide at the right (3' site)
   *  or full branch of (i,j)
   */
  e = extend_fm_3p(i, j, fm, vc);

  /*
   *  extension with one unpaired nucleotide at 5' site
   *  and all other variants which are needed for odd
   *  dangle models
   */
  if (sn[i - 1] == sn[i]) {
    if (sn[i] == sn[i + 1]) {
      if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        if (fm[ij + 1] != INF) {
          en = fm[ij + 1] +
               P->MLbase;

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[i][1];

            if (sc->f)
              en += sc->f(i, j, i + 1, j, VRNA_DECOMP_ML_ML, sc->data);
          }

          e = MIN2(e, en);
        }
      }
    }

    /* extension with bound ligand on 5'site */
    if (with_ud) {
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        k = i + u - 1;
        if ((k < j) && (sn[i] == sn[k + 1])) {
          if (evaluate(i, j, k + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
            if (fm[ij + u] != INF) {
              en = domains_up->energy_cb(vc,
                                         i,
                                         k,
                                         VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                         domains_up->data);
              if (en != INF) {
                en += fm[ij + u] +
                      u * P->MLbase;

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[i][u];

                  if (sc->f)
                    en += sc->f(i, j, k + 1, j, VRNA_DECOMP_ML_ML, sc->data);
                }

                e = MIN2(e, en);
              }
            }
          }
        }
      }
    }

    if (dangle_model % 2) {
      /* dangle_model = 1 || 3 */

      mm5 = ((i > 1) || circular) ? S[i] : -1;
      mm3 = ((j < length) || circular) ? S[j] : -1;

      if (sn[i] == sn[i + 1]) {
        if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
          if (c[ij + 1] != INF) {
            type = vrna_get_ptype(ij + 1, ptype);

            en = c[ij + 1] +
                 E_MLstem(type, mm5, -1, P) +
                 P->MLbase;

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[i][1];

              if (sc->f)
                en += sc->f(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, sc->data);
            }

            e = MIN2(e, en);
          }
        }
      }

      if (sn[j - 1] == sn[j]) {
        if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
          if (c[indx[j - 1] + i] != INF) {
            type = vrna_get_ptype(indx[j - 1] + i, ptype);

            en = c[indx[j - 1] + i] +
                 E_MLstem(type, -1, mm3, P) +
                 P->MLbase;

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[j][1];

              if (sc->f)
                en += sc->f(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, sc->data);
            }

            e = MIN2(e, en);
          }
        }
      }

      if ((sn[j - 1] == sn[j]) && (sn[i] == sn[i + 1])) {
        if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
          if (c[indx[j - 1] + i + 1] != INF) {
            type = vrna_get_ptype(indx[j - 1] + i + 1, ptype);

            en = c[indx[j - 1] + i + 1] +
                 E_MLstem(type, mm5, mm3, P) +
                 2 * P->MLbase;

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[j][1] + sc->energy_up[i][1];

              if (sc->f)
                en += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, sc->data);
            }

            e = MIN2(e, en);
          }
        }
      }
    } /* end special cases for dangles == 1 || dangles == 3 */
  }

  /* modular decomposition -------------------------------*/
  k1j   = indx[j] + i + turn + 2;
  stop  = (strands > 1) ? (se[0]) : (j - 2 - turn);

  /* duplicated code is faster than conditions in loop */
  if (hc->f) {
    if (sc && sc->f) {
      for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k] + fm[k1j];
          en      += sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          decomp  = MIN2(decomp, en);
        }
      }
      k++;
      k1j++;
      for (; k <= j - 2 - turn; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k] + fm[k1j];
          en      += sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          decomp  = MIN2(decomp, en);
        }
      }
    } else {
      for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k] + fm[k1j];
          decomp  = MIN2(decomp, en);
        }
      }
      k++;
      k1j++;
      for (; k <= j - 2 - turn; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k] + fm[k1j];
          decomp  = MIN2(decomp, en);
        }
      }
    }
  } else {
    if (sc && sc->f) {
      for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF)) {
          en      = fmi[k] + fm[k1j];
          en      += sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          decomp  = MIN2(decomp, en);
        }
      }
      k++;
      k1j++;
      for (; k <= j - 2 - turn; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF)) {
          en      = fmi[k] + fm[k1j];
          en      += sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          decomp  = MIN2(decomp, en);
        }
      }
    } else {
#ifdef VRNA_WITH_SSE_IMPLEMENTATION

      /* modular decomposition -------------------------------*/

      decomp = modular_decomposition(i, ij, j, turn, fmi, vc->matrices->fML);
      /* end modular decomposition -------------------------------*/

#else
      for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF)) {
          en      = fmi[k] + fm[k1j];
          decomp  = MIN2(decomp, en);
        }
      }
      k++;
      k1j++;
      for (; k <= j - 2 - turn; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF)) {
          en      = fmi[k] + fm[k1j];
          decomp  = MIN2(decomp, en);
        }
      }
#endif
    }
  }

  dmli[j] = decomp;               /* store for use in fast ML decompositon */
  e       = MIN2(e, decomp);

  /* coaxial stacking */
  if (dangle_model == 3) {
    /* additional ML decomposition as two coaxially stacked helices */
    int ik;
    k1j = indx[j] + i + turn + 2;
    for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++) {
      ik = indx[k] + i;
      if (evaluate(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        type    = rtype[vrna_get_ptype(ik, ptype)];
        type_2  = rtype[vrna_get_ptype(k1j, ptype)];

        en = c[ik] +
             c[k1j] +
             P->stack[type][type_2];

        if (sc)
          if (sc->f)
            en += sc->f(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, sc->data);

        decomp = MIN2(decomp, en);
      }
    }
    k++;
    k1j++;
    for (; k <= j - 2 - turn; k++, k1j++) {
      ik = indx[k] + i;
      if (evaluate(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        type    = rtype[vrna_get_ptype(ik, ptype)];
        type_2  = rtype[vrna_get_ptype(k1j, ptype)];

        en = c[ik] +
             c[k1j] +
             P->stack[type][type_2];

        if (sc)
          if (sc->f)
            en += sc->f(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL, sc->data);

        decomp = MIN2(decomp, en);
      }
    }

    decomp += 2 * P->MLintern[1];        /* no TermAU penalty if coax stack */
#if 0
    /*
     * This is needed for Y shaped ML loops with coax stacking of
     * interior pairts, but backtracking will fail if activated
     */
    DMLi[j] = MIN2(DMLi[j], decomp);
    DMLi[j] = MIN2(DMLi[j], DMLi[j - 1] + P->MLbase);
    DMLi[j] = MIN2(DMLi[j], DMLi1[j] + P->MLbase);
    new_fML = MIN2(new_fML, DMLi[j]);
#endif
    e = MIN2(e, decomp);
  }

  fmi[j] = e;

  return e;
}


PRIVATE int
E_ml_stems_fast_window(vrna_fold_compound_t *vc,
                       int                  i,
                       int                  j,
                       int                  *fmi,
                       int                  *dmli)
{
  char                      **ptype, type, type_2, tt;
  short                     *S1;
  int                       dangle_model, **c, **fML, e, decomp, en, en2, k, turn, *rtype;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  e             = INF;
  ptype         = vc->ptype_local;
  S1            = vc->sequence_encoding;
  P             = vc->params;
  md            = &(P->model_details);
  c             = vc->matrices->c_local;
  fML           = vc->matrices->fML_local;
  hc            = vc->hc;
  sc            = vc->sc;
  type          = vrna_get_ptype_window(i, j, ptype);
  turn          = md->min_loop_size;
  rtype         = &(md->rtype[0]);
  dangle_model  = md->dangles;
  evaluate      = prepare_hc_default_window(vc, &hc_dat_local);

  /*
   *  extension with one unpaired nucleotide at the right (3' site)
   *  or full branch of (i,j)
   */
  e = extend_fm_3p_window(i, j, fML, vc);

  /*
   *  extension with one unpaired nucleotide at 5' site
   *  and all other variants which are needed for odd
   *  dangle models
   */
  if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    decomp = fML[i + 1][j - i - 1] +
             P->MLbase;

    if (sc) {
      if (sc->energy_up)
        decomp += sc->energy_up[i][1];

      if (sc->f)
        decomp += sc->f(i, j, i + 1, j, VRNA_DECOMP_ML_ML, sc->data);
    }

    e = MIN2(e, decomp);
  }

  if (dangle_model % 2) {
    /* i+1,j */
    if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
      tt = vrna_get_ptype_window(i + 1, j, ptype);

      decomp = c[i + 1][j - i - 1] +
               E_MLstem(tt, S1[i], -1, P) +
               P->MLbase;

      if (sc) {
        if (sc->energy_up)
          decomp += sc->energy_up[i][1];

        if (sc->f)
          decomp += sc->f(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, sc->data);
      }

      e = MIN2(e, decomp);
    }

    /* i, j-1 */
    if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
      tt = vrna_get_ptype_window(i, j - 1, ptype);

      decomp = c[i][j - 1 - i] +
               E_MLstem(tt, -1, S1[j], P) +
               P->MLbase;

      if (sc) {
        if (sc->energy_up)
          decomp += sc->energy_up[j][1];

        if (sc->f)
          decomp += sc->f(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, sc->data);
      }

      e = MIN2(e, decomp);
    }

    /* i+1,j-1 */
    if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
      tt = vrna_get_ptype_window(i + 1, j - 1, ptype);

      decomp = c[i + 1][j - 1 - i - 1] +
               E_MLstem(tt, S1[i], S1[j], P) +
               2 * P->MLbase;

      if (sc) {
        if (sc->energy_up)
          decomp += sc->energy_up[i][1] +
                    sc->energy_up[j][1];

        if (sc->f)
          decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, sc->data);
      }

      e = MIN2(e, decomp);
    }
  }

  /* modular decomposition -------------------------------*/
  if (sc && sc->f) {
    if (hc->f) {
      for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++)
        if ((fmi[k - i] != INF) && (fML[k + 1][j - k - 1] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k - i] + fML[k + 1][j - k - 1];
          en      += sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          decomp  = MIN2(decomp, en);
        }
    } else {
      for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++) {
        en  = fmi[k - i];
        en2 = fML[k + 1][j - k - 1];
        if ((en != INF) && (en2 != INF))
          en += en2 + sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

        decomp = MIN2(decomp, en);
      }
    }
  } else {
    if (hc->f) {
      for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++)
        if ((fmi[k - i] != INF) && (fML[k + 1][j - k - 1] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k - i] + fML[k + 1][j - k - 1];
          decomp  = MIN2(decomp, en);
        }
    } else {
      for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++) {
        en  = fmi[k - i];
        en2 = fML[k + 1][j - k - 1];
        if ((en != INF) && (en2 != INF))
          decomp = MIN2(decomp, en + en2);
      }
    }
  }

  dmli[j - i] = decomp;               /* store for use in ML decompositon */
  e           = MIN2(e, decomp);

  /* coaxial stacking */
  if (dangle_model == 3) {
    /* additional ML decomposition as two coaxially stacked helices */
    for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++) {
      if (evaluate(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        type    = rtype[vrna_get_ptype_window(i, k, ptype)];
        type_2  = rtype[vrna_get_ptype_window(k + 1, j, ptype)];

        en = c[i][k - i] +
             c[k + 1][j - k - 1] +
             P->stack[type][type_2];

        if (sc)
          if (sc->f)
            en += sc->f(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, sc->data);

        decomp = MIN2(decomp, en);
      }
    }

    decomp += 2 * P->MLintern[1];          /* no TermAU penalty if coax stack */
#if 0
    /* This is needed for Y shaped ML loops with coax stacking of
    * interior pairts, but backtracking will fail if activated */
    dmli[j - i] = MIN2(dmli[j - i], decomp);
    dmli[j - i] = MIN2(dmli[j - i], dmli[j - 1 - i] + P->MLbase);
    dmli[j - i] = MIN2(dmli[j - i], dmli1[j - (i + 1)] + P->MLbase);
    e           = MIN2(e, dmli[j - i]);
#endif
    e = MIN2(e, decomp);
  }

  fmi[j - i] = e;

  return e;
}


PRIVATE int
E_ml_stems_fast_comparative(vrna_fold_compound_t  *vc,
                            int                   i,
                            int                   j,
                            int                   *fmi,
                            int                   *dmli)
{
  unsigned char             *hard_constraints;
  short                     **S, **S5, **S3;
  unsigned int              **a2s;
  int                       e, energy, *c, *fML, *ggg, ij, *indx, s, n_seq, k,
                            dangle_model, decomp, turn, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_mx_mfe_t             *matrices;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq             = vc->n_seq;
  matrices          = vc->matrices;
  P                 = vc->params;
  md                = &(P->model_details);
  c                 = matrices->c;
  fML               = matrices->fML;
  ggg               = matrices->ggg;
  indx              = vc->jindx;
  hc                = vc->hc;
  scs               = vc->scs;
  hard_constraints  = hc->matrix;
  dangle_model      = md->dangles;
  turn              = md->min_loop_size;
  a2s               = vc->a2s;
  ij                = indx[j] + i;
  e                 = INF;
  evaluate          = prepare_hc_default(vc, &hc_dat_local);

  /*
   *  extension with one unpaired nucleotide at the right (3' site)
   *  or full branch of (i,j)
   */
  e = extend_fm_3p_comparative(i, j, fML, vc);

  /*
   *  extension with one unpaired nucleotide at 5' site
   *  and all other variants which are needed for odd
   *  dangle models
   */
  if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    if (fML[ij + 1] != INF) {
      energy = fML[ij + 1] + n_seq * P->MLbase;
      if (scs) {
        for (s = 0; s < n_seq; s++) {
          if (scs[s]) {
            if (scs[s]->energy_up)
              energy += scs[s]->energy_up[a2s[s][i]][1];

            if (scs[s]->f)
              energy += scs[s]->f(i, j, i + 1, j, VRNA_DECOMP_ML_ML, scs[s]->data);
          }
        }
      }

      e = MIN2(e, energy);
    }
  }

  /* modular decomposition -------------------------------*/
#ifdef VRNA_WITH_SSE_IMPLEMENTATION
  decomp = modular_decomposition(i, ij, j, turn, fmi, vc->matrices->fML);
#else
  decomp = INF;
  if (hc->f) {
    for (k = i + 1 + turn; k <= j - 2 - turn; k++) {
      if (hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        energy  = fmi[k] + fML[indx[j] + k + 1];
        decomp  = (energy < decomp) ? energy : decomp;
      }
    }
  } else {
    for (k = i + 1 + turn; k <= j - 2 - turn; k++) {
      energy  = fmi[k] + fML[indx[j] + k + 1];
      decomp  = (energy < decomp) ? energy : decomp;
    }
  }

#endif
  dmli[j] = decomp; /* store for later use in ML decompositon */

  e = MIN2(e, decomp);

  fmi[j] = e; /* store for later use in ML decompositon */

  return e;
}


PRIVATE int
E_ml_stems_fast_comparative_window(vrna_fold_compound_t *vc,
                                   int                  i,
                                   int                  j,
                                   int                  *fmi,
                                   int                  *dmli)
{
  unsigned char             **hard_constraints;
  short                     **S, **S5, **S3;
  unsigned int              **a2s;
  int                       e, energy, **c, **fML, **ggg, s, n_seq, k,
                            dangle_model, decomp, turn, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_mx_mfe_t             *matrices;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq             = vc->n_seq;
  matrices          = vc->matrices;
  P                 = vc->params;
  md                = &(P->model_details);
  c                 = matrices->c_local;
  fML               = matrices->fML_local;
  ggg               = matrices->ggg_local;
  hc                = vc->hc;
  scs               = vc->scs;
  hard_constraints  = hc->matrix_local;
  dangle_model      = md->dangles;
  turn              = md->min_loop_size;
  a2s               = vc->a2s;
  e                 = INF;
  evaluate          = prepare_hc_default_window(vc, &hc_dat_local);

  if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    energy = fML[i + 1][j - (i + 1)] + n_seq * P->MLbase;
    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->energy_up)
            energy += scs[s]->energy_up[a2s[s][i]][1];
      }
    }

    e = MIN2(e, energy);
  }

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    energy = fML[i][j - 1 - i] + n_seq * P->MLbase;
    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->energy_up)
            energy += scs[s]->energy_up[a2s[s][j]][1];
      }
    }

    e = MIN2(e, energy);
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
    energy = c[i][j - i];

    S   = vc->S;
    S5  = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
    S3  = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */

    if (dangle_model) {
      for (s = 0; s < n_seq; s++) {
        tt      = vrna_get_ptype_md(S[s][i], S[s][j], md);
        energy  += E_MLstem(tt, S5[s][i], S3[s][j], P);
      }
    } else {
      for (s = 0; s < n_seq; s++) {
        tt      = vrna_get_ptype_md(S[s][i], S[s][j], md);
        energy  += E_MLstem(tt, -1, -1, P);
      }
    }

    e = MIN2(e, energy);
  }

  if (md->gquad) {
    decomp = ggg[i][j - i] +
             n_seq * E_MLstem(0, -1, -1, P);
    e = MIN2(e, decomp);
  }

  /* modular decomposition -------------------------------*/
  decomp = INF;
  if (hc->f) {
    for (k = i + 1 + turn; k <= j - 2 - turn; k++) {
      if (hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        energy  = fmi[k - i] + fML[k + 1][j - (k + 1)];
        decomp  = (decomp > energy) ? energy : decomp;
      }
    }
  } else {
    for (k = i + 1 + turn; k <= j - 2 - turn; k++) {
      energy  = fmi[k - i] + fML[k + 1][j - (k + 1)];
      decomp  = (decomp > energy) ? energy : decomp;
    }
  }

  dmli[j - i] = decomp; /* store for later use in ML decompositon */

  e = MIN2(e, decomp);

  fmi[j - i] = e; /* store for later use in ML decompositon */

  return e;
}
