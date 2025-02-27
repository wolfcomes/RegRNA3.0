
/*=============================================================================
 mscores.c                        A.Lambert le 29/11/01  revu le 30/11/01

 Fonctions operant des calculs de scores d'un masque de 'pattern'.
 Le calcul des scores depuis des configurations limitees aux elements brins et
 helices retenus par le masque.
 Il peut s'ensuivre une forte reduction du nombre de configurations.


 cc -O2 -Wall -c mscores.c -I ../include ;

 ar -r ../lib/librnaIV.a mscores.o ;

=============================================================================*/

#include "rnaIV.h"

extern double LOG_ZERO;

void  GetMaskProfiles(Trset *trset, Mask *mask);
void  SetMaskScoresTab(Mask *mask, int len, int level);
void  FreeMaskScoresTab(Mask *mask);
void  GetMaskScoresTab(char *seq, int len, Mask *mask, int level);
float GetMaskScore(Mask *mask, int tab_index, int level);

/*=============================================================================
 GetMaskProfiles(): Cree l'ensemble des profils d'un masque de pattern.
                    Si un seul masque est recherche il n'est pas necessaire de
                    calculer tous les profils du pattern.
                    Afin d'eviter les allocations multiples, 'GetStProfile' et
                    'GetHlxProfile' ne creent les profils que si les pointeurs
                    'St->Profile' et 'Hlx->Profile' ont la valeur NULL dans le
                    fichier 'profs.c'.
=============================================================================*/

void GetMaskProfiles(Trset *trset, Mask *mask)
{
  int i, *m;

  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
  {
      GetStProfile(trset, mask->pattern->std + *m);
  }
  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
  {
      GetHlxProfile(trset, mask->pattern->hlx + *m);
  }
  return;
}
/*=============================================================================
 SetMaskScoresTab(): Cree les tableaux de scores des helices et brins d'un mas-
                     que de pattern.
                     'len' represente la longueur d'une sequence a explorer.
                     pour chaque helice ou brin la longueur du tableau cree
                     sera 'len - hlx.max_len + 1' ou 'len - std.max_len + 1'.
                     Si 'level' est egal a 0 les tableaux 'Scores' sont crees
                     sinon ce sera les tableaux 'ScoresBis'.
                     Un controle est effectue pour ne pas creer plusieurs fois
                     le tableau des brins et helices presents dans plusieurs
                     masques.
=============================================================================*/

void SetMaskScoresTab(Mask *mask, int len, int level)
{
  int   i, *m;

  if (level == 0)
  {
      for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
      {
          if (mask->pattern->std[*m].Scores == NULL )
              mask->pattern->std[*m].Scores =
                                fMat(mask->pattern->std[*m].max_gaps + 1,
                                     len - mask->pattern->std[*m].max_len + 1);
      }
      for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
      {
          if (mask->pattern->hlx[*m].Scores == NULL)          
              mask->pattern->hlx[*m].Scores =
                                fMat(mask->pattern->hlx[*m].max_gaps + 1,
                                     len - mask->pattern->hlx[*m].max_len + 1);
      }
  }
  else
  {
      for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
      {
          if (mask->pattern->std[*m].ScoresBis == NULL)  
              mask->pattern->std[*m].ScoresBis =
                                fMat(mask->pattern->std[*m].max_gaps + 1,
                                     len - mask->pattern->std[*m].max_len + 1);
      }
      for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
      {
          if (mask->pattern->hlx[*m].ScoresBis == NULL)  
              mask->pattern->hlx[*m].ScoresBis =
                                fMat(mask->pattern->hlx[*m].max_gaps + 1,
                                     len - mask->pattern->hlx[*m].max_len + 1);
      }
  }

  return;
}
/*=============================================================================
 FreeMaskScoresTab(): Libere les tableaux de scores d'un masque de pattern.
=============================================================================*/

void FreeMaskScoresTab(Mask *mask)
{
  int i, *m;

  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
  {
      if (mask->pattern->std[*m].Scores != NULL) {
          FreefMat(mask->pattern->std[*m].Scores);
          mask->pattern->std[*m].Scores = NULL;
      }
  }
  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
  {
      if (mask->pattern->hlx[*m].Scores != NULL) {
          FreefMat(mask->pattern->hlx[*m].Scores);
          mask->pattern->hlx[*m].Scores = NULL;
      }
  }
  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
  {
      if (mask->pattern->std[*m].ScoresBis != NULL) {
          FreefMat(mask->pattern->std[*m].ScoresBis);
          mask->pattern->std[*m].ScoresBis = NULL;
      }
  }
  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
  {
      if (mask->pattern->hlx[*m].ScoresBis != NULL) {
          FreefMat(mask->pattern->hlx[*m].ScoresBis);
          mask->pattern->hlx[*m].ScoresBis = NULL;
      }
  }

  return;
}
/*=============================================================================
 GetMaskScoresTab(): Calcule les tableaux de scores d'un masque de pattern.
                     'len' represente la longueur d'une sequence a explorer
                     pour chaque helice ou brin la longueur du tableau cree
                     sera 'len - hlx.max_len' ou 'len - std.max_len'.
                     Si 'level' est egal a 0 les tableaux 'Scores' sont calcules
                     sinon ce sera les tableaux 'ScoresBis'.
=============================================================================*/

void GetMaskScoresTab(char *seq, int len, Mask *mask, int level)
{
  int i, *m;

  if (level == 0)
  {
      for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
      {
          GetStScoresTab(seq, len - mask->pattern->std[*m].max_len + 1,
                         mask->pattern->std + *m);
      }
      for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
      {
          GetHlxScoresTab(seq, len - mask->pattern->hlx[*m].max_len + 1,
                          mask->pattern->hlx + *m);
      }
  }
  else
  {
      for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
      {
          GetStScoresTabBis(seq, len - mask->pattern->std[*m].max_len + 1,
                            mask->pattern->std + *m);
      }
      for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
      {
          GetHlxScoresTabBis(seq, len - mask->pattern->hlx[*m].max_len + 1,
                             mask->pattern->hlx + *m);
      }
  }
  return;
}
/*=============================================================================
 GetMaskScore(): Calcule depuis l'abscisse 'tab_index' des tableaux de scores
                 precalcules le score maximal parmi les configurations d'un
                 masque de pattern.
                 En retour 'mask->cfgindex' pointe l'index, dans le tableau des
                 configurations, de celle dont le score est retenu.
                 Le score maximal est aussi retourne (pour commodite).
                 Si 'level' est egal a 0 les tableaux 'Scores' sont calcules,
                 sinon ce sera les tableaux 'ScoresBis'.
=============================================================================*/

float GetMaskScore(Mask *mask, int tab_index, int level)
{
  int   i, s, h, I, J, *m;
  float score, score_max;

  extern double LOG_ZERO;

  score_max = mask->pattern->max_len * LOG_ZERO;
  mask->cfgindex = 0;

  if (level == 0)
  {
      for (i = 0; i < mask->ncfg; i++)      /* boucle sur les configurations */
      {
          score = 0.0;

           for (s = 0, m = mask->stindex; s < mask->nst; s++, m++)  /* brins */
           {
               I = (int) mask->cfg[i].stgaps[s];
               J = tab_index + (int) mask->cfg[i].stbgn[s];

               score += mask->pattern->std[*m].Scores[I][J];
           }
           for (h = 0, m = mask->hxindex; h < mask->nhx; h++, m++)/* helices */
           {
               I = (int) mask->cfg[i].hxgaps[h];
               J = tab_index + (int) mask->cfg[i].hxbgn[h];

               score += mask->pattern->hlx[*m].Scores[I][J];
           }
           if (score > score_max) {
               score_max = score;
               mask->cfgindex = i;
           }
       }
  }
  else
  {
      for (i = 0; i < mask->ncfg; i++)      /* boucle sur les configurations */
      {
          score = 0.0;

           for (s = 0, m = mask->stindex; s < mask->nst; s++, m++)  /* brins */
           {
               I = (int) mask->cfg[i].stgaps[s];
               J = tab_index + (int) mask->cfg[i].stbgn[s];

               score += mask->pattern->std[*m].ScoresBis[I][J];
           }
           for (h = 0, m = mask->hxindex; h < mask->nhx; h++, m++)/* helices */
           {
               I = (int) mask->cfg[i].hxgaps[h];
               J = tab_index + (int) mask->cfg[i].hxbgn[h];

               score += mask->pattern->hlx[*m].ScoresBis[I][J];
           }
           if (score > score_max) {
               score_max = score;
               mask->cfgindex = i;
           }
       }
  }
  mask->score = score_max;                /* enregistrement du score maximal */

  return score_max;
}
/*===========================================================================*/
