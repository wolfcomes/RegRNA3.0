#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/exterior_loops.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/hairpin_loops.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "hairpin_loops.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
eval_hp_loop_fake(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

/**
 *  @brief  Evaluate the free energy of a hairpin loop
 *          and consider possible hard constraints
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 */
PUBLIC int
vrna_E_hp_loop(vrna_fold_compound_t *vc,
               int                  i,
               int                  j)
{
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  if (vc->hc->type == VRNA_HC_WINDOW)
    evaluate = prepare_hc_default_window(vc, &hc_dat_local);
  else
    evaluate = prepare_hc_default(vc, &hc_dat_local);

  if ((i > 0) && (j > 0)) {
    /* is this base pair allowed to close a hairpin (like) loop ? */
    if (evaluate(i, j, i, j, VRNA_DECOMP_PAIR_HP, &hc_dat_local)) {
      if (j > i)  /* linear case */
        return vrna_eval_hp_loop(vc, i, j);
      else        /* circular case */
        return vrna_eval_ext_hp_loop(vc, j, i);
    }
  }

  return INF;
}


/**
 *  @brief  Evaluate the free energy of an exterior hairpin loop
 *          and consider possible hard constraints
 */
PUBLIC int
vrna_E_ext_hp_loop(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j)
{
  return vrna_E_hp_loop(vc, j, i);
}


/**
 *  @brief Evaluate free energy of an exterior hairpin loop
 *
 *  @ingroup eval
 *
 */
PUBLIC int
vrna_eval_ext_hp_loop(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j)
{
  char          **Ss, loopseq[10];
  unsigned int  **a2s;
  short         *S, *S2, **SS, **S5, **S3;
  int           u1, u2, e, s, type, n_seq, length, noGUclosure;
  vrna_param_t  *P;
  vrna_sc_t     *sc, **scs;
  vrna_md_t     *md;

  length      = vc->length;
  P           = vc->params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  e           = INF;

  u1  = length - j;
  u2  = i - 1;

  if ((u1 + u2) < 3)
    return e;

  switch (vc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S           = vc->sequence_encoding;
      S2          = vc->sequence_encoding2;
      sc          = vc->sc;
      type        = vrna_get_ptype_md(S2[j], S2[i], md);
      loopseq[0]  = '\0';

      if (noGUclosure && ((type == 3) || (type == 4)))
        break;

      /* maximum special hp loop size: 6 */
      if ((u1 + u2) < 7) {
        strcpy(loopseq, vc->sequence + j - 1);
        strncat(loopseq, vc->sequence, i);
      }

      e = E_Hairpin(u1 + u2, type, S[j + 1], S[i - 1], loopseq, P);

      if (sc) {
        if (sc->energy_up)
          e += sc->energy_up[j + 1][u1] +
               sc->energy_up[1][u2];

        if (sc->f)
          e += sc->f(j, i, j, i, VRNA_DECOMP_PAIR_HP, sc->data);
      }

      break;

    /* sequence alignments */
    case  VRNA_FC_TYPE_COMPARATIVE:
      SS    = vc->S;
      S5    = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
      S3    = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
      Ss    = vc->Ss;
      a2s   = vc->a2s;
      scs   = vc->scs;
      n_seq = vc->n_seq;
      e     = 0;

      for (s = 0; s < n_seq; s++) {
        char loopseq[10];
        u1          = a2s[s][length] - a2s[s][j];
        u2          = a2s[s][i - 1];
        loopseq[0]  = '\0';

        if ((u1 + u2) < 7) {
          strcpy(loopseq, Ss[s] + a2s[s][j] - 1);
          strncat(loopseq, Ss[s], a2s[s][i]);
        }

        if ((u1 + u2) < 3) {
          e += 600;
        } else {
          type  = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
          e     += E_Hairpin(u1 + u2, type, S3[s][j], S5[s][i], loopseq, P);
        }
      }
      if (scs) {
        for (s = 0; s < n_seq; s++) {
          if (scs[s]) {
            if (scs[s]->energy_up)
              e += scs[s]->energy_up[1][a2s[s][i - 1]] +
                   scs[s]->energy_up[a2s[s][j + 1]][a2s[s][length] - a2s[s][j]];

            if (scs[s]->f) {
              e += scs[s]->f(a2s[s][j],
                             a2s[s][i],
                             a2s[s][j],
                             a2s[s][i],
                             VRNA_DECOMP_PAIR_HP,
                             scs[s]->data);
            }
          }
        }
      }

      break;

    /* nothing */
    default:
      break;
  }

  return e;
}


/**
 *  @brief Evaluate free energy of a hairpin loop
 *
 *  @ingroup eval
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @param  vc  The #vrna_fold_compound_t for the particular energy evaluation
 *  @param  i   5'-position of the base pair
 *  @param  j   3'-position of the base pair
 *  @returns    Free energy of the hairpin loop closed by @f$ (i,j) @f$ in deka-kal/mol
 */
PUBLIC int
vrna_eval_hp_loop(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j)
{
  char          **Ss;
  unsigned int  **a2s;
  short         *S, *S2, **SS, **S5, **S3;
  unsigned int  *sn;
  int           u, e, s, ij, type, *idx, n_seq, en, noGUclosure;
  vrna_param_t  *P;
  vrna_sc_t     *sc, **scs;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;

  idx         = vc->jindx;
  P           = vc->params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  sn          = vc->strand_number;
  domains_up  = vc->domains_up;
  e           = INF;

  if (sn[j] != sn[i])
    return eval_hp_loop_fake(vc, i, j);

  /* regular hairpin loop */
  switch (vc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      S2    = vc->sequence_encoding2;
      sc    = vc->sc;
      u     = j - i - 1;
      type  = vrna_get_ptype_md(S2[i], S2[j], md);

      if (noGUclosure && ((type == 3) || (type == 4)))
        break;

      e = E_Hairpin(u, type, S[i + 1], S[j - 1], vc->sequence + i - 1, P);

      /* add soft constraints */
      if (sc) {
        if (sc->energy_up)
          e += sc->energy_up[i + 1][u];

        switch (sc->type) {
          case VRNA_SC_DEFAULT:
            if (sc->energy_bp) {
              ij  = idx[j] + i;
              e   += sc->energy_bp[ij];
            }

            break;

          case VRNA_SC_WINDOW:
            if (sc->energy_bp_local)
              e += sc->energy_bp_local[i][j - i];

            break;
        }

        if (sc->f)
          e += sc->f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
      }

      /* consider possible ligand binding */
      if (domains_up && domains_up->energy_cb) {
        en = domains_up->energy_cb(vc,
                                   i + 1, j - 1,
                                   VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                   domains_up->data);
        if (en != INF)
          en += e;

        e = MIN2(e, en);
      }

      break;

    /* sequence alignments */
    case  VRNA_FC_TYPE_COMPARATIVE:
      SS    = vc->S;
      S5    = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
      S3    = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
      Ss    = vc->Ss;
      a2s   = vc->a2s;
      scs   = vc->scs;
      n_seq = vc->n_seq;

      for (e = s = 0; s < n_seq; s++) {
        u = a2s[s][j - 1] - a2s[s][i];
        if (u < 3) {
          e += 600;                          /* ??? really 600 ??? */
        } else {
          type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
          e     += E_Hairpin(u, type, S3[s][i], S5[s][j], Ss[s] + (a2s[s][i - 1]), P);
        }
      }

      if (scs) {
        for (s = 0; s < n_seq; s++) {
          if (scs[s]) {
            u = a2s[s][j - 1] - a2s[s][i];

            if (scs[s]->energy_up)
              e += scs[s]->energy_up[a2s[s][i + 1]][u];

            switch (scs[s]->type) {
              case VRNA_SC_DEFAULT:
                if (scs[s]->energy_bp) {
                  ij  = idx[j] + i;
                  e   += scs[s]->energy_bp[ij];
                }

                break;

              case VRNA_SC_WINDOW:
                if (scs[s]->energy_bp_local)
                  e += scs[s]->energy_bp_local[i][j - i];

                break;
            }

            if (scs[s]->f) {
              e += scs[s]->f(a2s[s][i],
                             a2s[s][j],
                             a2s[s][i],
                             a2s[s][j],
                             VRNA_DECOMP_PAIR_HP,
                             scs[s]->data);
            }
          }
        }
      }

      break;

    /* nothing */
    default:
      break;
  }

  return e;
}


PRIVATE int
eval_hp_loop_fake(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j)
{
  short         *S, *S2;
  unsigned int  *sn;
  int           u, e, ij, type, *idx, en, noGUclosure;
  vrna_param_t  *P;
  vrna_sc_t     *sc;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;

  idx         = vc->jindx;
  P           = vc->params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  sn          = vc->strand_number;
  domains_up  = vc->domains_up;
  e           = INF;

  switch (vc->type) {
    /* single sequences and cofolding hybrids */
    case  VRNA_FC_TYPE_SINGLE:
      S     = vc->sequence_encoding;
      S2    = vc->sequence_encoding2;
      sc    = vc->sc;
      u     = j - i - 1;
      ij    = idx[j] + i;
      type  = vrna_get_ptype_md(S2[j], S2[i], md);

      if (noGUclosure && ((type == 3) || (type == 4)))
        break;

      /* hairpin-like exterior loop (for cofolding) */
      short si, sj;

      si  = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;
      sj  = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;

      if (md->dangles)
        e = E_ExtLoop(type, sj, si, P);
      else
        e = E_ExtLoop(type, -1, -1, P);

      /* add soft constraints */
      if (sc) {
        if (sc->energy_up)
          e += sc->energy_up[i + 1][u];

        if (sc->energy_bp)
          e += sc->energy_bp[ij];

        if (sc->f)
          e += sc->f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
      }

      /* consider possible ligand binding */
      if (domains_up && domains_up->energy_cb) {
        en = domains_up->energy_cb(vc,
                                   i + 1, j - 1,
                                   VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                   domains_up->data);
        if (en != INF)
          en += e;

        e = MIN2(e, en);
      }

      break;

    /* nothing */
    default:
      break;
  }

  return e;
}
