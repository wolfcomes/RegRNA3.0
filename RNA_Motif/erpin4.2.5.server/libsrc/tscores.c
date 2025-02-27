
/*=============================================================================
 tscores.c                        A.Lambert  le 02/12/01 revu le 05/12/01

 Ce code regroupe les fonctions procedant au calcul et a la statistique des sco-
 res de patterns et masques des sequences d'une "base d'entrainement".
 Les calculs sont faits sur les sequences de l'alignement multiple sans suppri-
 mer les gaps.

 cc -O2 -Wall -c tscores.c -I../include ;

 ar -r ../lib/librnaIV.a tscores.o ;

=============================================================================*/

#include "rnaIV.h"

extern double LOG_ZERO;

double GetHlxTScore(char *seq, Helix *Hlx);
double GetStTScore(char *seq, Strand *St);
double GetPatternTScore(char *seq, Pattern *pattern);
double GetMaskTScore(char *seq, Mask *mask);

double *GetPatternTScores(Trset *trset, Pattern *pattern);
double *GetMaskTScores(Trset *trset, Mask *mask);

int    cmpdbl(const void *x1, const void *x2);
void   fPrintfCaptures(double *ts_scores, int nscores, FILE *txt);
double ConvertRatio(int percent, double *ts_scores, int nscores);
void   fPrintfScoresStat(double *scores, int nscores, FILE *txt);

void   fPrintfPatternTStat(Trset *trset, Pattern *pattern, FILE *txt);
void   fPrintfMaskTStat(Trset *trset, Mask *mask, FILE *txt);
double GetPatternThreshold(int percent, Trset *trset, Pattern *pattern);
double GetMaskThreshold(int percent, Trset *trset, Mask *mask);
void   GetMasksThresholds(int *percent, Trset *trset, Mask *mask, int nmask);

/*=============================================================================
 GetHlxTScore(): Calcule, depuis la position pointee par 'seq', le score de
                 l'helice pointee par 'Hlx' dans une base d'entrainement,
                 c'est a dire avec l'eloignement maximal des deux brins de
                 l'helice, tel qu'il se presente dans l'alignement multiple.
                 Variante de 'GetHlxScores' de 'scores.c'.
=============================================================================*/

double GetHlxTScore(char *seq, Helix *Hlx)
{
  int     i, j;
  char    *right;
  double  score;
  extern  short **NtHlxCode;               /* declare dans 'libsrc/ntcode.c' */

  right = seq + Hlx->max_len - 1;

  for (j = 0, score = 0.0; j < Hlx->helix_len; j++, right--)
  {
      i = NtHlxCode[(int) seq[j]][(int) *right];
      score += Hlx->Profile[i][j];
  }
  return score;
}
/*=============================================================================
 GetStTScore(): Calcule, depuis la position pointee par 'seq', le score du brin
                pointe par 'St' dans une base d'entrainement (avec gaps), c'est
	        a dire avec la longueur maximale du brin, tel qu'il se presente
                dans l'alignement multiple.
=============================================================================*/

double GetStTScore(char *seq, Strand *St)
{
  int    i, j;
  double score;
  extern short *NtStCode;                  /* declare dans 'libsrc/ntcode.c' */

  for (j = 0, score = 0.0; j < St->max_len; j++)
  {
      i = NtStCode[(int) seq[j]];
      score += St->Profile[i][j];
  }
  return score;
}
/*=============================================================================
 GetPatternTScore(): Calcule, depuis la position pointee par 'seq', le score du
                     pattern pointe par 'pattern'.
=============================================================================*/

double GetPatternTScore(char *seq, Pattern *pattern)
{
  int     i;
  double  score;

  score = 0.0;

  for (i = 0; i < pattern->nst; i++)                 /* boucle sur les brins */
  {
      score += GetStTScore(seq + pattern->std[i].db_bgn, pattern->std + i);
  }
  for (i = 0; i < pattern->nhx; i++)               /* boucle sur les helices */
  {
      score += GetHlxTScore(seq + pattern->hlx[i].db_bgn1, pattern->hlx + i);
  }
   
  return score;
}
/*=============================================================================
 GetPatternTScores(): Calcule, pour l'ensemble des sequences de la base d'en-
                      trainement pointee par 'trset', le score du pattern pointe
                      par 'pattern'.
                      Retourne un pointeur sur le tableau des scores calcules
                      reordonnes dans le sens ascendant pour un traitement sta-
                      tistique ulterieur.
=============================================================================*/

double *GetPatternTScores(Trset *trset, Pattern *pattern)
{
  int    i;
  double *scores;

  scores = (double *) malloc(trset->nseq * sizeof(double));

  for (i = 0; i < trset->nseq; i++)
      scores[i] = GetPatternTScore(trset->data[i], pattern);

  qsort(scores, trset->nseq, sizeof(double), &cmpdbl);      /* tri ascendant */
 
  return scores;
}
/*=============================================================================
 GetMaskTScore(): Calcule, depuis la position pointee par 'seq', le score du
                  masque pointe par 'mask'.
=============================================================================*/

double GetMaskTScore(char *seq, Mask *mask)
{
  int     i, *m;
  double  score;

  score = 0.0;

  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
  {
      score += GetStTScore(seq + mask->pattern->std[*m].db_bgn,
                           mask->pattern->std + *m);
  }
  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
  {
      score += GetHlxTScore(seq + mask->pattern->hlx[*m].db_bgn1,
                            mask->pattern->hlx + *m);
  }

  return score;
}
/*=============================================================================
 GetMaskTScores(): Calcule, pour l'ensemble des sequences de la base d'entrai-
                   nement pointee par 'trset', le score du masque pointe par
                   'mask'.
                   Retourne un pointeur sur le tableau des scores calcules
                   reordonnes dans le sens ascendant pour un traitement statis-
                   tique ulterieur.
=============================================================================*/

double *GetMaskTScores(Trset *trset, Mask *mask)
{
  int    i;
  double *scores;

  scores = (double *) malloc(trset->nseq * sizeof(double));

  for (i = 0; i < trset->nseq; i++)
      scores[i] = GetMaskTScore(trset->data[i], mask);

  qsort(scores, trset->nseq, sizeof(double), &cmpdbl);      /* tri ascendant */

  return scores;
}
/*=============================================================================
 cmpdbl(): Retourne -1 ou 1 selon que la difference (x1 - x2) est <= 0.0 ou > 0.
           Cette fonction sera utilisee pour trier 'npts' points d'un tableau 
           de doubles dans le sens ascendant par la fonction 'qsort' de la bi-
           bliotheque standard de C:  
           qsort(T, npts, sizeof(double), &cmpdbl);
=============================================================================*/

int cmpdbl(const void *x1, const void *x2)
{
  return (*(double *)x1 - *(double *)x2 <= 0.0 ? -1 : 1);
}
/*=============================================================================
 fPrintfCaptures(): Etant donnes un tableau de scores trie dans l'ordre ascen-
                    dant et un tableau de pourcentages de scores saisis, affi-
                    che pour chaque pourcentage le 'score-seuil' correspondant.
=============================================================================*/

void fPrintfCaptures(double *ts_scores, int nscores, FILE *txt)
{
  int  i, j, ratios[10] = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10};

  fprintf(txt, "Cutoff and ratios of retained sequences:\n");

  for (j = 0, i = 0; j < nscores && i < 10; j++)
  {
      if ((int) rint(100 * (nscores - j) / (double) nscores) <= ratios[i])
      {
          fprintf(txt, "%.2f\t%d%%\n", ts_scores[j!=0 ? j-1 : 0], ratios[i]);
          i++;
      }
  }
  fprintf(txt, "\n");

  return;
}
/*=============================================================================
 ConvertRatio(): Retourne, le seuil correspondant au pourcentage 'percent' de
                 valeurs a saisir dans le tableau pointe par 'ts_scores' de
                 'nscores' elements, lequel est suppose trie dans l'ordre ascen-
                 dant.
=============================================================================*/

double ConvertRatio(int percent, double *ts_scores, int nscores)
{
  int    j;
  double threshold = ts_scores[nscores - 1];      /* initialise au % minimal */

  if (percent > 100) percent = 100;          /* intervalle ramene a [0, 100] */
  else
  if (percent < 0) percent = 0;

  for (j = 0; j < nscores; j++)
      if ((int) rint(100 * (nscores - j) / (double) nscores) <= percent)
      {
          threshold = ts_scores[j != 0 ? j - 1 : 0];
          break;
      }
  return threshold - 1.e-3;                   /* un peu reduit pour securite */
}
/*=============================================================================
 fPrintfScoresStat(): Calcule les valeurs extremes, moyenne et ecart type d'un
                      tableau de scores trie dans l'ordre ascendant.
                      Les resultats sont diriges sur le fichier 'txt'.
=============================================================================*/

void fPrintfScoresStat(double *scores, int nscores, FILE *txt)
{
  int    i;
  double min, max, mean, meansqr, stdev;

  min = scores[0];
  max = scores[nscores - 1];

  mean = meansqr = 0.0;
  
  for (i = 0; i < nscores; i++)
  {
      mean += scores[i];
      meansqr += scores[i] * scores[i];
  }
  mean /= nscores;
  meansqr /= nscores;

  stdev = sqrt(meansqr - mean * mean);

  fprintf(txt, "Statistics of scores:\n");
  fprintf(txt, "min:  %.2f\n", min);
  fprintf(txt, "max:  %.2f\n", max);
  fprintf(txt, "mean: %.2f\n", mean);
  fprintf(txt, "std:  %.2f\n\n", stdev);

  return;
}
/*=============================================================================
 fPrintfPatternTStat(): Procede au calcul des scores d'un pattern et a leur sta-
                        tistique, les resultats sont diriges sur 'txt'.
=============================================================================*/

void fPrintfPatternTStat(Trset *trset, Pattern *pattern, FILE *txt)
{
  double *scores;

  scores = GetPatternTScores(trset, pattern);

  fprintf(txt, "\n");

  fprintf(txt, "===========================================================\n");
  fprintf(txt, " REGION (%s) STATISTICS:\n", pattern->id);
  fprintf(txt, "===========================================================\n");

  fPrintfCaptures(scores, trset->nseq, txt);
  fPrintfScoresStat(scores, trset->nseq, txt);

  free(scores);

  return;
}
/*=============================================================================
 fPrintfMaskTStat(): Procede au calcul des scores d'un masque et a leur statis-
                     tique, les resultats sont diriges sur le fichier 'txt'.
=============================================================================*/

void fPrintfMaskTStat(Trset *trset, Mask *mask, FILE *txt)
{
  double *scores;
  char   *str = (char *) malloc(80);
  int    i;

  scores = GetMaskTScores(trset, mask);

  sprintf(str, "{ ");
  for (i = 0; i < mask->nargs; i++) {
      sprintf(str, "%s%d ", str, mask->args[i]);
  }
  sprintf(str, "%s%c", str, '}');

  fprintf(txt, "===========================================================\n");
  fprintf(txt, " REGION %s STATISTICS:\n", str);
  fprintf(txt, "===========================================================\n");

  fPrintfCaptures(scores, trset->nseq, txt);
  fPrintfScoresStat(scores, trset->nseq, txt);

  free(scores);
  free(str);

  return;
}
/*=============================================================================
 GetPatternThreshold(): Retourne la valeur numerique du seuil correspondant au
                        pourcentage de captures 'percent' de la structure poin-
                        tee par 'pattern' dans la base d'entrainement pointee
                        par 'trset'.
=============================================================================*/

double GetPatternThreshold(int percent, Trset *trset, Pattern *pattern)
{
  double *scores, threshold;

  scores = GetPatternTScores(trset, pattern);
  threshold = ConvertRatio(percent, scores, trset->nseq);

  free(scores);

  return threshold;
}
/*=============================================================================
 GetMaskThreshold(): Retourne la valeur numerique du seuil correspondant au
                     pourcentage de captures 'percent' de la structure pointee
                     par 'mask' dans la base d'entrainement pointee par 'trset'.
                     Cette valeur est aussi enregistree dans 'mask->threshold'.
=============================================================================*/

double GetMaskThreshold(int percent, Trset *trset, Mask *mask)
{
  double *scores;

  scores = GetMaskTScores(trset, mask);
  mask->threshold = ConvertRatio(percent, scores, trset->nseq);

  free(scores);

  return mask->threshold;
}
/*=============================================================================
 GetMasksThresholds(): Identique a 'GetMaskThreshold' mais pour un tableau de
                       'nmask' masques, les pourcentages sont pointes en entree
                       par un tableau d'entiers 'percent'.
=============================================================================*/

void GetMasksThresholds(int *percent, Trset *trset, Mask *mask, int nmask)
{
  int    i;
  double *scores;

  for (i = 0; i < nmask; i++)
  {
      scores = GetMaskTScores(trset, mask + i);
      mask[i].threshold = ConvertRatio(percent[i], scores, trset->nseq);

      free(scores);
  }

  return;
}
/*===========================================================================*/
