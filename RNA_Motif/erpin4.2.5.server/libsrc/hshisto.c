
/*=============================================================================
 hshisto.c                                      A.Lambert le 15/04/04
 
 Ce code contient les fonctions permettant de construire les histogrammes des
 scores des helices et brins mesures sur des sequences aleatoires.

 cc -O2 -Wall -c hshisto.c -I../include ;

 ar -rs ../lib/librnaIV.a hshisto.o ;

=============================================================================*/

#include "rnaIV.h"

#define EPSILON (1.e-8)

double **CatStProfiles(Mask *mask, int *width);
double **CatStNoGProfs(Mask *mask, int *width);
double **CatHlxProfiles(Mask *mask, int *width);

/* brins sans gaps */

Histo  StsNoGHist(double **Prof, int width, double dh, double *Prfsc);
Histo  GetStNoGHist(Strand *St, double dh, double *Prfsc);

/* brins avec gaps */

double **GetWFScoresProb(Strand *St, double *bkgfreqs);
float  *RndAlnScores(Strand *St, int samples, double *Prfsc);
Histo  GetAlnHist(Strand *St, double dx, double *Pfs, int samples);

/* helices */

void GetSSProfileStat(double **SSprof, int width, double *Prfsc,
                      double *min, double *max, double *mean, double *std);
Histo GetSSHist(double **SSProf, int width, double dh, double *Prfsc);

/*=============================================================================
 CatStProfiles(): Cree et pointe en retour un tableau ou sont concatenes les
                  profils de l'ensemble des brins du masque pointe par l'argu-
		  ment 'mask'.
		  En retour 'width' pointe la largeur du tableau dont la hau-
		  teur est 'APLHA_LEN + 1'.
=============================================================================*/

double **CatStProfiles(Mask *mask, int *width)
{
  double  **Stprof;
  Strand  *std;
  int     i, j, k, col, *m;

  *width = 0;                               /* mesure de la largeur du profil */
  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
  {
      std = mask->pattern->std + *m;
      *width += std->max_len;
  }

  Stprof = dMat(ALPHA_LEN + 1, *width);               /* creation du tableau */

  col = 0;
  for (k = 0, m = mask->stindex; k < mask->nst; k++, m++)
  {
      std = mask->pattern->std + *m;

      for (j = 0; j < std->max_len; j++)
      {
          for (i = 0; i <= ALPHA_LEN; i++)
              Stprof[i][col] = std->Profile[i][j];
	  col++;
      }
  }
  return Stprof;
}
/*=============================================================================
 CatStNoGProfs(): Cree et pointe en retour un tableau ou sont concatenes les
                  profils des brins sans gaps du masque pointe par l'argument
        	  'mask'.
		  En retour 'width' pointe la largeur du tableau dont la hau-
		  teur est 'APLHA_LEN'.
		  Si il n'y a pas de brins sans gaps retourne NULL.
=============================================================================*/

double **CatStNoGProfs(Mask *mask, int *width)
{
  double  **Stprof;
  Strand  *std;
  int     i, j, k, col, *m;

  *width = 0;                               /* mesure de la largeur du profil */
  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
  {
      std = mask->pattern->std + *m;
      if (std->max_gaps == 0)
          *width += std->max_len;
  }
  if (*width == 0) return NULL;            /* il n'y a pas de brins sans gaps */

  Stprof = dMat(ALPHA_LEN, *width);                    /* creation du tableau */

  col = 0;
  for (k = 0, m = mask->stindex; k < mask->nst; k++, m++)
  {
      std = mask->pattern->std + *m;
      if (std->max_gaps == 0)
          for (j = 0; j < std->max_len; j++)
          {
              for (i = 0; i < ALPHA_LEN; i++)
                  Stprof[i][col] = std->Profile[i][j];
	      col++;
          }
  }
  return Stprof;
}
/*=============================================================================
 CatHlxProfiles(): Retourne un pointeur sur le tableau bidimensionnel reunissant
                   les profils de l'ensemble des helices du masque pointe par
                   'mask'.
                   Les dimensions de ce tableau sont 'width x SQR_ALPHA_LEN' ou
                   'width', qui est pointe en retour, est la somme des largeurs
                   des profils des helices concernees.
                   C'est la longueur de la "structure secondaire generique".
                   Les lignes des profils au dela de SQR_ALPHA_LEN sont ignorees.
=============================================================================*/

double **CatHlxProfiles(Mask *mask, int *width)
{
  double  **Hlxprof;
  Helix   *hlx;
  int     i, j, k, col, *m;

  *width = 0;
  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
  {
      hlx = mask->pattern->hlx + *m;
      *width += hlx->helix_len;                  /* longueur totale de la ss */
  }

  Hlxprof = dMat(SQR_ALPHA_LEN, *width);              /* creation du tableau */

  col = 0;                                              /* copie des profils */
  for (k = 0, m = mask->hxindex; k < mask->nhx; k++, m++)
  {
      hlx = mask->pattern->hlx + *m;

      for (j = 0; j < hlx->helix_len; j++)
      {
          for (i = 0; i < SQR_ALPHA_LEN; i++)
              Hlxprof[i][col] = hlx->Profile[i][j];
          col++;
      }
  }
  return Hlxprof;
}

/* ------------------------------- BRINS SANS GAPS --------------------------*/

/*=============================================================================
 StsNoGHist(): Retourne l'histogramme des scores d'un groupe de brin sans gaps,
               depuis le profil de largeur 'width' pointe par 'Prof'.
               Les scores sont supposes etre evalues sur des sequences de
	       composition aleatoire.
	       'Prfsc' pointe en sortie la probabilite des scores finis.
	       L'histogramme est obtenu en effectuant le produit de convolu-
	       tion des histogrammes des colonnes du profil (distribution de
               la somme des variables aleatoires associees aux colonnes).
=============================================================================*/

Histo StsNoGHist(double **Prof, int width, double dh, double *Prfsc)
{
  extern double *LogDataFreqs;
  const  double minusInf = 0.8 * LOG_ZERO;
  Histo  GetNHistW(double *, double *, int, double);

  Histo  hist1, hist2;
  double *scores, *weights, *bkgfreqs, probfsc;
  int    i, j, nfsc;

  *Prfsc = 1.0;
  scores = padd(0.0, ALPHA_LEN);
  weights = padd(0.0, ALPHA_LEN);
  bkgfreqs = (double *) malloc(ALPHA_LEN * sizeof(double));

  for (i = 0; i < ALPHA_LEN; i++)  bkgfreqs[i] = exp(LogDataFreqs[i]);

                                               /* premiere colonne du profil */

  for (i = nfsc = 0, probfsc = 0.0; i < ALPHA_LEN; i++)  {
      if (Prof[i][0] > minusInf) {
          probfsc += bkgfreqs[i];
	  weights[nfsc] = bkgfreqs[i];
          scores[nfsc] = Prof[i][0];
	  nfsc++;
      }
  }
  *Prfsc *= probfsc;
  hist1 = GetNHistW(scores, weights, nfsc, dh);

  for (j = 1; j < width; j++)           /* boucle sur les colonnes suivantes */
  {
      for (i = nfsc = 0, probfsc = 0.0; i < ALPHA_LEN; i++)  {
          if (Prof[i][j] > minusInf) {
              probfsc += bkgfreqs[i];
	      weights[nfsc] = bkgfreqs[i];
              scores[nfsc] = Prof[i][j];
	      nfsc++;
          }
      }
      *Prfsc *= probfsc;
      hist2 = GetNHistW(scores, weights, nfsc, dh);

      ConvHist(&hist1, hist2);
      NormalizeHist(hist1);
      free((double *) hist2.vals);
  }

  free(scores);
  free(weights);
  free(bkgfreqs);

  return hist1;
}
/*=============================================================================
 GetStNoGHist(): Retourne l'histogramme des scores d'un brin sans gaps,
                 Les scores sont supposes etre evalues sur des sequences de
		 composition aleatoire.
	         'Prfsc' pointe en sortie la probabilite des scores finis.
		 L'histogramme est obtenu en effectuant le produit de convolu-
		 tion des histogrammes des colonnes du profil (distribution de
	         la somme des variables aleatoires associees aux colonnes).
		 Variante de la fonction precedente 'StsNoGHist'.
=============================================================================*/

Histo  GetStNoGHist(Strand *St, double dh, double *Prfsc)
{
  extern double *LogDataFreqs;
  const  double minusInf = 0.8 * LOG_ZERO;
  Histo  GetNHistW(double *, double *, int, double);

  Histo  hist1, hist2;
  double **Prof, *scores, *weights, *bkgfreqs, probfsc;
  int    i, j, nfsc, width;

  Prof = St->Profile;
  width = St->max_len;

  *Prfsc = 1.0;
  scores = padd(0.0, ALPHA_LEN);
  weights = padd(0.0, ALPHA_LEN);
  bkgfreqs = (double *) malloc(ALPHA_LEN * sizeof(double));

  for (i = 0; i < ALPHA_LEN; i++)  bkgfreqs[i] = exp(LogDataFreqs[i]);

                                               /* premiere colonne du profil */

  for (i = nfsc = 0, probfsc = 0.0; i < ALPHA_LEN; i++)  {
      if (Prof[i][0] > minusInf) {
          probfsc += bkgfreqs[i];
	  weights[nfsc] = bkgfreqs[i];
          scores[nfsc] = Prof[i][0];
	  nfsc++;
      }
  }
  *Prfsc *= probfsc;
  hist1 = GetNHistW(scores, weights, nfsc, dh);

  for (j = 1; j < width; j++)           /* boucle sur les colonnes du profil */
  {
      for (i = nfsc = 0, probfsc = 0.0; i < ALPHA_LEN; i++)  {
          if (Prof[i][j] > minusInf) {
              probfsc += bkgfreqs[i];
	      weights[nfsc] = bkgfreqs[i];
              scores[nfsc] = Prof[i][j];
	      nfsc++;
          }
      }
      *Prfsc *= probfsc;
      hist2 = GetNHistW(scores, weights, nfsc, dh);

      ConvHist(&hist1, hist2);
      NormalizeHist(hist1);
      free((double *) hist2.vals);
  }

  free(scores);
  free(weights);
  free(bkgfreqs);

  return hist1;
}
/* -------------------------- BRINS CONTENANT DES GAPS ----------------------*/

/*=============================================================================
 GetWFScoresProb(): Pour le brin avec gaps pointe par 'St' et des frequences de
                    bases pointees par 'bkgfreqs', cette fonction retourne un
                    pointeur sur un tableau donnant pour chaque colonne la dis-
                    tribution des bases possibles conduisant a un score fini.
                    La 5eme ligne de ce tableau enregistre la probabilite de
                    score fini.
                    La 6eme ligne contient des 1 et 0 suivant que la colonne
                    correspondante dans la base d'entrainement contient des
                    gaps ou non.
=============================================================================*/

double **GetWFScoresProb(Strand *St, double *bkgfreqs)
{
  const   double  minusInf = 0.8 * LOG_ZERO;
  double  **fprobs, sum;
  int     i, j, height, width;

  height = ALPHA_LEN + 2;
  width = St->max_len;

  fprobs = dMat(height, width);
  FilldMat(fprobs, height, width, 0.0);

  for (j = 0; j < width; j++)           /* boucle sur les colonnes du profil */
  {
      sum = 0.0;
      for (i = 0; i < ALPHA_LEN; i++)
      {
          if (St->Profile[i][j] > minusInf) {
              fprobs[i][j] = bkgfreqs[i];
              sum += bkgfreqs[i];
          }
      }
      fprobs[ALPHA_LEN][j] = sum;                     /* prob. de score fini */
                                             /* prob. parmi les scores finis */
      for (i = 0; i < ALPHA_LEN; i++) fprobs[i][j] /= sum;

      if (St->Profile[ALPHA_LEN][j] > minusInf)  /* des gaps dans la colonne */
          fprobs[ALPHA_LEN + 1][j] = 1.0;
  }

  return fprobs;
}
/*=============================================================================
 RndAlnScores(): Pour le brin avec gaps pointe par 'St' cette fonction calcule
                 un echantillonnage de 'samples' brins de scores finis et re-
                 tourne un pointeur sur le tableau des scores cree.
                 En retour 'Prfsc' pointe la probabilite des scores finis.
                 Note: Pour chaque brin aleatoire la fonction saisit un score
                 fini, qui est calcule sur une configuration aleatoire de gaps.
                 Le nombre de gaps est distribue uniformement sur l'echantillon.
=============================================================================*/

                                  /* tirage uniforme d'un entier dans [0, n[ */
#define RND(n)   (int)((double)(n) * rand() / RMAXp1)
#define PERMUTS  4

float *RndAlnScores(Strand *St, int samples, double *Prfsc)
{
  extern double *LogDataFreqs;
  const double RMAXp1 = RAND_MAX + 1.0;

  int     *gapcfgs;   /* tableau des indices situant les gaps du brin aligne */
  float   *scores;
  double  *freqs, **fprobs, prfs, sum;
  char    *seq;
  int     i, j, k, p, r, tmp, width, maxgaps, ngaps, gtotal;

  width = St->max_len;
  maxgaps = St->max_gaps;

  seq = (char *) malloc(width + 1);              /* alloc. et initialisation */
  for (j = 0; j < width; j++)  seq[j] = RndNucl(0.25, 0.25, 0.25);
  seq[width] = '\0';

  scores = (float *) malloc(samples * sizeof(float));
  freqs = (double *) malloc(ALPHA_LEN * sizeof(double));
  for (j = 0; j < ALPHA_LEN; j++)  freqs[j] = exp(LogDataFreqs[j]);

  fprobs = GetWFScoresProb(St, freqs);

  for (j = gtotal = 0; j < width; j++)       /* 'gtotal' est >= St->max_gaps */
      if ((int)fprobs[ALPHA_LEN + 1][j] == 1) gtotal++;

          /* 'gapcfgs' contient les indices des colonnes du profil avec gaps */

  gapcfgs = (int *) malloc(gtotal * sizeof(int));
  for (j = k = 0; j < width; j++)
      if ((int)fprobs[ALPHA_LEN + 1][j] == 1) gapcfgs[k++] = j;

         /* --------- Creation de l'echantillon de brins avec gaps --------- */

  for (p = 0, sum = 0.0; p < samples; p++)
  {
      prfs = 1.0;
      ngaps = RND(maxgaps + 1);         /* tirage uniforme dans [0, maxgaps] */

      if (ngaps != 0)
      {
          for (k = 0; k < PERMUTS; k++)   /* permute plusieurs fois les elts */
              for (j = 0; j < gtotal; j++)
	      {
		  r = RND(gtotal);           /* un entier dans [0, gtotal-1] */
                  tmp = gapcfgs[j];
                  gapcfgs[j] = gapcfgs[r];
                  gapcfgs[r] = tmp;
              }
          SortAsc(gapcfgs, ngaps);  /* tri ascend. des 'ngaps' premiers elts */
      }

      for (i = j = k = 0; j < width; j++)
      {
          if (j == gapcfgs[k] && k < ngaps)
	      k++;                          /* position d'insertion d'un gap */
	  else {
              seq[i++] = RndNucl(fprobs[0][j], fprobs[1][j], fprobs[2][j]);
              prfs *= fprobs[ALPHA_LEN][j];
	  }
      }

      St->len = width;
      AlignSProfile(seq, St);       /* alignement du brin cree sur le profil */
                              /* saisie du score correspondant au nb de gaps */
      scores[p] = (float) St->Align[width - ngaps][width];
      sum += prfs;
  }
  *Prfsc = sum / samples;

  free(gapcfgs);
  free(seq);
  free(freqs);
  FreedMat(fprobs);

  return scores;
}
#undef RND
#undef PERMUTS

/*=============================================================================
 GetAlnHist(): Retourne l'histogramme des scores du brin pointe par 'St' conte-
               nant des gaps sur un echantillonnage de 'samples' elements alea-
               toires obtenus a l'issue d'un alignement.
               'dx' est la largeur d'un element de cet histogramme.
               En retour 'Pfs' pointe la probabilite des scores finis.
=============================================================================*/

Histo GetAlnHist(Strand *St, double dx, double *Pfs, int samples)
{
  Histo   hist;
  float   *scores;

  scores = RndAlnScores(St, samples, Pfs);
  hist = GetNHist(scores, samples, dx);             /* histogramme normalise */
  free(scores);

  return hist;
}

/* ---------------------------------- HELICES -------------------------------*/

/*=============================================================================
 GetSSProfileStat(): calcule les elements statistiques du profil pointe par
                     'SSprofile' (genere par 'GetSSProfile') de largeur 'width'.
                     'pfsc' pointe en sortie la probabilite d'obtenir un score
                     fini (calcule sur les elements superieurs a 'logzero' dans
                     les profils initiaux).
=============================================================================*/

void GetSSProfileStat(double **SSprof, int width, double *Prfsc,
                      double *min, double *max, double *mean, double *std)
{
  extern double *LogDataFreqs;

  const  double minusInf = 0.8 * LOG_ZERO;
  int    i, j;
  double scmin, scmax, scmean, scmean2, scvar, pf, *freqs;

  *min = *max = 0.0;

  for (j = 0; j < width; j++)      /* calcul des scores extremes: min et max */
  {
      scmin =  1.e+6;
      scmax = -1.e+6;

      for (i = 0; i < SQR_ALPHA_LEN; i++)
      {
          if (SSprof[i][j] > minusInf)
          {
              if (SSprof[i][j] < scmin) scmin = SSprof[i][j];
              if (SSprof[i][j] > scmax) scmax = SSprof[i][j];
          }
      }
      *min += scmin;
      *max += scmax;
  }
  *min -= 0.01;            /* au cas ou les valeurs extremes seraient egales */
  *max += 0.01;

  freqs = (double *) malloc(ALPHA_LEN * sizeof(double));
  for (i = 0; i < ALPHA_LEN; i++)  freqs[i] = exp(LogDataFreqs[i]);

#ifdef DEBUG
{
  double sum = 0.0;
  for (i = 0; i < ALPHA_LEN; i++)  sum += freqs[i];
  fprintf(stderr, "\nFreq_ATGC_Sum = %.3f\n\n", sum);        /* verification */
}
#endif

  *mean = scvar = 0.0;
  *Prfsc = 1.0;

  for (j = 0; j < width; j++)        /* calcul de la moyenne et l'ecart type */
  {
      scmean = scmean2 = pf = 0.0;

      for (i = 0; i < SQR_ALPHA_LEN; i++)
      {
          double P = freqs[i / ALPHA_LEN] * freqs[i % ALPHA_LEN];

          if (SSprof[i][j] > minusInf)
          {
              pf += P;
              scmean  += P * SSprof[i][j];
              scmean2 += P * SSprof[i][j] * SSprof[i][j];
          }
      }
      scmean  /= pf;                    /* pf est toujours au moins egal a 1 */
      scmean2 /= pf;
                                   /* independance statistique des colonnes: */
      scvar += scmean2 - scmean * scmean;             /* somme des variances */
      *mean += scmean;                                 /* somme des moyennes */
      *Prfsc *= pf;                                      /* produit des prob. */
  }
  *std = (scvar < EPSILON ? 0.0 : sqrt(scvar));         /* il faudra traiter */
                                /* separement le cas ou l'ecart type est nul */
  free(freqs);

  return;
}
/*=============================================================================
 GetSSHist(): Retourne l'histogramme des scores d'un masque d'helices dont le
              profil (concatenation des profils d'helices) est pointe par
	      'SSProf' de largeur 'width'.
	      Cet histogramme est obtenu par le produit de convolution des his-
	      togrammes des colonnes du profil (sauf exclusions).
              'dh' indique l'intervalle elementaire de l'axe de l'histogramme
              et 'Prfsc' pointe en sortie la probabilite des scores finis.
=============================================================================*/

Histo GetSSHist(double **SSProf, int width, double dh, double *Prfsc)
{
  extern double *LogDataFreqs;
  const  double minusInf = 0.8 * LOG_ZERO;
  Histo  GetNHistW(double *, double *, int, double);

  Histo  hist1, hist2;
  double *scores, *weights, *bkgfreqs, probfsc;
  int    i, j, nfsc;

  *Prfsc = 1.0;
  scores = padd(0.0, SQR_ALPHA_LEN);
  weights = padd(0.0, SQR_ALPHA_LEN);

  bkgfreqs = (double *) malloc(SQR_ALPHA_LEN * sizeof(double));
  for (i = 0; i < SQR_ALPHA_LEN; i++)
      bkgfreqs[i] = exp(LogDataFreqs[i / ALPHA_LEN]) *
                    exp(LogDataFreqs[i % ALPHA_LEN]);

                                  /* 1er histogramme: 1ere colonne du profil */

  for (i = nfsc = 0, probfsc = 0.0; i < SQR_ALPHA_LEN; i++)  {
      if (SSProf[i][0] > minusInf) {
          probfsc += bkgfreqs[i];
	  weights[nfsc] = bkgfreqs[i];
          scores[nfsc] = SSProf[i][0];
	  nfsc++;
      }
  }
  *Prfsc *= probfsc;
  hist1 = GetNHistW(scores, weights, nfsc, dh);

  for (j = 1; j < width; j++)    /* col. suivantes: convol. des histogrammes */
  {
      for (i = nfsc = 0, probfsc = 0.0; i < SQR_ALPHA_LEN; i++)  {
          if (SSProf[i][j] > minusInf) {
              probfsc += bkgfreqs[i];
	      weights[nfsc] = bkgfreqs[i];
              scores[nfsc] = SSProf[i][j];
	      nfsc++;
          }
      }
      *Prfsc *= probfsc;
      hist2 = GetNHistW(scores, weights, nfsc, dh);

      ConvHist(&hist1, hist2);
      NormalizeHist(hist1);
      free((double *) hist2.vals);
  }
  free(scores);
  free(weights);
  free(bkgfreqs);

  return hist1;
}
/*===========================================================================*/

#undef EPSILON
