
/*=============================================================================
 mhisto.c                           A.Lambert le 24/03/03 revu le 15/04/04

 Ce code contient des fonctions destinees a operer la statistique des scores
 d'un masque sur une sequence aleatoire.
 Les histogrammes des scores des helices et brins sont d'abord calcules, puis
 celui concernant un masque.
 Les sequences aleatoires sont supposees a bases successives independantes, les
 frequences des bases sont parametrees.

 cc -O2 -Wall -c mhisto.c -I../include ;

 ar -rs ../lib/librnaIV.a mhisto.o ;

=============================================================================*/

#include "rnaIV.h"

void   GetSSStat(Mask *mask, double *pfsc, double *min, double *max,
                 double *mean, double *std);
void   GetTScoresStat(double *scores, int nscores,
                      double *min, double *max , double *mean , double *std);
void   GetMaskTStat(Trset *trset, Mask *mask,
                    double *min, double *max , double *mean , double *std);

Histo  GetHlxHist(Mask *mask, double dx, double *Prfsc);
Histo  GetStHist(Mask *mask, double dx, double *Prfsc, int spercol);
Histo  GetStsHist(Mask *mask, double dx, double *Prfsc, int spercol);
Histo  GetMaskHist(Mask *mask, double dx, double *Prfsc, int spercol);

/*=============================================================================
 GetSSStat(): En retour pointe les elements statistiques de la structure secon-
              daire (helices) du masque pointe par 'mask'.
=============================================================================*/

void GetSSStat(Mask *mask, double *pfsc, double *min, double *max,
               double *mean, double *std)
{
  double  **SSprof;
  int     width;

  SSprof = CatHlxProfiles(mask, &width);
  GetSSProfileStat(SSprof, width, pfsc, min, max, mean, std);
  FreedMat(SSprof);

  return;
}
/*=============================================================================
 GetTScoresStat(): calcule les elements statistiques (moyenne, dev.std, valeurs
                   extremes) des scores contenus dans le tableau pointe par
                   'scores', lequel a ete trie dans l'ordre ascendant (voir la
                   fonction 'GetMaskTScores' dans 'tscores.c').
                   les resultats sont pointes par 'min', 'max', 'mean', 'std'.
=============================================================================*/
#define EPSILON (1.e-8)

void GetTScoresStat(double *scores, int nscores,
                    double *min, double *max , double *mean , double *std)
{
  double meansqr, v;
  int    i;

  *min = scores[0];              /* 'scores' est trie dans l'ordre croissant */
  *max = scores[nscores - 1];

  *mean = meansqr = 0.0;

  for (i = 0; i < nscores; i++) {
      *mean += scores[i];
      meansqr += scores[i] * scores[i];
  }
  *mean /= nscores;
  meansqr /= nscores;

  v = meansqr - *mean * *mean;
  *std = (v < EPSILON ? 0.0 : sqrt(v));

  return;
}
#undef EPSILON

/*=============================================================================
 GetMaskTStat(): Procede au calcul des scores d'un masque et a leur statistique
                 les resultats sont, en retour, pointes par les 4 derniers
                 arguments 'min', 'max', 'mean', 'std'.
=============================================================================*/

void GetMaskTStat(Trset *trset, Mask *mask,
                  double *min, double *max , double *mean , double *std)
{
  double *scores;

  scores = GetMaskTScores(trset, mask);
  GetTScoresStat(scores, trset->nseq, min, max , mean , std);
  free(scores);
  return;
}
/*=============================================================================
 GetHlxHist(): Retourne l'histogramme, dont l'integrale est normalisee a 1, des
               scores des helices sur une sequence aleatoire du masque pointe
               par 'mask'.
               'dx' est la largeur des intervalles de cet histogramme.
               En retour 'Prfsc' pointe la probabilite des scores finis.
=============================================================================*/

Histo GetHlxHist(Mask *mask, double dx, double *Prfsc)
{
  Histo  hist;
  double **SSProf;
  int    width;

  hist.vals = NULL;         /*valeur retourne si le nombre d'helices est nul */
  if (mask->nhx == 0) return hist;

  SSProf = CatHlxProfiles(mask, &width);
  hist = GetSSHist(SSProf, width, dx, Prfsc);

  FreedMat(SSProf);
  return hist;
}
/*=============================================================================
 GetStHist(): Retourne l'histogramme des scores des brins (avec et sans gaps)
              sur une sequence aleatoire du masque pointe par 'mask.
              - Si le masque ne contient pas de brins le champ 'vals' de la
              structure retournee est mise a 'NULL'.
              Le nombre d'echantillons examines, dans le cas ou le brin contient
	      des gaps, croit comme le carre du nombre de colonnes du profil:
	      'spercol' fixe le nombre d'echantillons pour la premiere colonne.
              'dx' est la mesure d'un element de cet histogramme.
=============================================================================*/

#define MAX_WIDTH 12

Histo GetStHist(Mask *mask, double dx, double *Prfsc, int spercol)
{
  Strand  *St;
  Histo   hist1, hist2;
  int     i, *m, samples, L;
  double  Probfsc;

  *Prfsc = 1.0;

  hist1.vals = NULL;         /* valeur retournee si le nbre de brins est nul */
  if (mask->nst == 0) return hist1;

  St = mask->pattern->std + mask->stindex[0];          /* 1er brin rencontre */

  if (St->max_gaps == 0)                                   /* brin sans gaps */
      hist1 = GetStNoGHist(St, dx, &Probfsc);
  else
  {
      L = St->max_len <= MAX_WIDTH ? St->max_len : MAX_WIDTH;   /* borne sup */
      samples = spercol * L * L;
      hist1 = GetAlnHist(St, dx, &Probfsc, samples);
  }
  *Prfsc *= Probfsc;

                                    /* examen des brins suivants s'il y en a */

  for (i = 1, m = mask->stindex + 1; i < mask->nst; i++, m++)
  {
      St = mask->pattern->std + *m;

      if (St->max_gaps == 0)                               /* brin sans gaps */
          hist2 = GetStNoGHist(St, dx, &Probfsc);
      else
      {
          L = St->max_len <= MAX_WIDTH ? St->max_len : MAX_WIDTH;
          samples = spercol * L * L;
          hist2 = GetAlnHist(St, dx, &Probfsc, samples);
      }

      *Prfsc *= Probfsc;
      ConvHist(&hist1, hist2);               /* convolution des histogrammes */
      free((double *) hist2.vals);
  }
  NormalizeHist(hist1);                        /* normalisation avant sortie */

  return hist1;
}
#undef MAX_WIDTH

/*=============================================================================
 GetStsHist(): Retourne l'histogramme des scores des brins (avec et sans gaps)
               sur une sequence aleatoire du masque pointe par 'mask.
               - Si le masque ne contient pas de brins le champ 'vals' de la
               structure retournee est mise a 'NULL'.
               Le nombre d'echantillons examines, dans le cas ou le brin con-
	       tient des gaps, croit comme le carre du nombre de colonnes du
	       profil:
	      'spercol' fixe le nombre d'echantillons pour la premiere colonne.
              'dx' est la mesure d'un element de cet histogramme.
=============================================================================*/

#define MAX_WIDTH 12

Histo GetStsHist(Mask *mask, double dx, double *Prfsc, int spercol)
{
  Strand  *St;
  Histo   hist1, hist2;
  int     i, *m, samples, L, width;
  double  **StsProf, Probfsc;

  *Prfsc = 1.0;

  hist1.vals = NULL;         /* valeur retournee si le nbre de brins est nul */
  if (mask->nst == 0) return hist1;

  StsProf = CatStNoGProfs(mask, &width);

  if (StsProf != NULL)                      /* il y a des brins sans gaps .. */
  {
      hist1 = StsNoGHist(StsProf, width, dx, Prfsc);
      FreedMat(StsProf);

      for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
      {
          St = mask->pattern->std + *m;

	  if (St->max_gaps != 0)
	  {
              L = St->max_len <= MAX_WIDTH ? St->max_len : MAX_WIDTH;
              samples = spercol * L * L;
              hist2 = GetAlnHist(St, dx, &Probfsc, samples);

              *Prfsc *= Probfsc;
              ConvHist(&hist1, hist2);       /* convolution des histogrammes */
              free((double *) hist2.vals);
	  }
      }
  }
  else                                /* il n'y a que des brins avec gaps .. */
  {
      St = mask->pattern->std + mask->stindex[0];      /* 1er brin rencontre */

      L = St->max_len <= MAX_WIDTH ? St->max_len : MAX_WIDTH;   /* borne sup */
      samples = spercol * L * L;
      hist1 = GetAlnHist(St, dx, &Probfsc, samples);
      *Prfsc *= Probfsc;
                          /* examen des brins avec gaps suivants s'il y en a */
      if (mask->nst > 1)
          for (i = 1, m = mask->stindex + 1; i < mask->nst; i++, m++)
          {
              St = mask->pattern->std + *m;

              L = St->max_len <= MAX_WIDTH ? St->max_len : MAX_WIDTH;
              samples = spercol * L * L;
              hist2 = GetAlnHist(St, dx, &Probfsc, samples);

              *Prfsc *= Probfsc;
              ConvHist(&hist1, hist2);       /* convolution des histogrammes */
              free((double *) hist2.vals);
          }
  }

  NormalizeHist(hist1);                        /* normalisation avant sortie */

  return hist1;
}
#undef MAX_WIDTH

/*=============================================================================
 GetMaskHist(): Retourne l'histogramme des scores du masque pointe par 'mask'
                sur un echantillonnage aleatoire.
                Si l'histogramme n'a pas pu etre cree le champ 'hist.vals' re-
                tourne est NULL.
                Pour les brins (avec ou sans gaps) le nombre d'echantillons est
                proportionnel au carre du nombre de colonnes du profil: 'spercol'
                fixe le nombre d'echantillons pour la 1ere colonne.
                'dx' est la mesure d'un element de cet histogramme.
                En retour 'Prfsc' pointe la probabilite des scores finis.
                Les profils du masque sont supposes avoir ete crees par la
                fonction 'GetMaskProfiles' de 'mscores.c'.
=============================================================================*/

Histo GetMaskHist(Mask *mask, double dx, double *Prfsc, int spercol)
{
  Histo   hist1, hist2;
  double  Probfsc;

  *Prfsc = 1.0;
  hist1.vals = NULL;

  if (mask->nhx != 0)
  {
      hist1 = GetHlxHist(mask, dx, &Probfsc);
      *Prfsc *= Probfsc;

      if (mask->nst != 0) {
          hist2 = GetStsHist(mask, dx, &Probfsc, spercol);
          *Prfsc *= Probfsc;
          ConvHist(&hist1, hist2);           /* convolution des histogrammes */
          NormalizeHist(hist1);
          free((double *) hist2.vals);
      }
  }
  else
  if (mask->nst != 0) {
      hist1 = GetStsHist(mask, dx, &Probfsc, spercol);
      *Prfsc *= Probfsc;
  }

  return hist1;
}
/*===========================================================================*/
