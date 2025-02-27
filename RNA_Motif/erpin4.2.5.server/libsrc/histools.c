
/*=============================================================================
 histools.c                        A.Lambert le 22/03/03 revu le 13/02/04

 Ce fichier regroupe des fonctions destinees a l'etude d'histogrammes.

 cc -O2 -Wall -c histools.c -I../include ;

 ar -rs ../lib/librnaIV.a histools.o ;

=============================================================================*/

#include "rnaIV.h"

#define EPSILON (1.e-8)

Histo  SetupHist(double min, double max, double dh, int samples);
Histo  GetNHist(float *scores, int samples, double dx);
Histo  GetNHistW(double *scores, double *weights, int samples, double dx);
void   NormalizeHist(Histo hist);
double *SetHistAxis(Histo hist);
void   PrHistInfo(FILE *txt, char *cmt, Histo hist);
void   PrHistDat(FILE *txt, Histo hist);
void   GetHistStat(Histo hist, double *mean, double *std);
double Gauss(double x, double mean, double std);
int    RandNucl(double p0, double p1, double p2);
char   RndNucl(double p0, double p1, double p2);
void   SortAsc(int *T, int n);
int    ResetHistBins(Histo *hist, int bins);

/*=============================================================================
 SetupHist(): Determine les elements geometriques d'une structure 'Histo' de-
              puis la donnee des parametres 'min', 'max' (valeurs extremes con-
              siderees), 'dh' (mesure de l'intervalle elementaire des valeurs).
              Le champ 'bins' de la structure 'Histo' aura toujours une valeur
              superieure a 5.
              Le tableau pointe par le champ 'vals' est cree et mis a zero.
              Le nombre d'echantillons 'samples' est passe pour completer
              l'initialisation.
=============================================================================*/

Histo SetupHist(double min, double max, double dh, int samples)
{
  Histo  hist;
  const  double Margin = 2.5;

  if (min > max) {
      fprintf(stderr, "SetupHist: arg#1 must be <= arg#2, exit..\n");
      exit(1);
  }
  if (dh <= 0.0) {
      fprintf(stderr, "SetupHist: arg#3 must be > 0.0, exit..\n");
      exit(1);
  }

  hist.hmin = min - Margin * dh;     /* l'histogramme aura au moins 5 points */
  hist.hmax = max + Margin * dh;
                                                        /* division arrondie */
  hist.bins = 1 + (int) ceil((hist.hmax - hist.hmin) / dh);
  hist.hmax = hist.hmin + (double)(hist.bins - 1) * dh;      /* reajustement */
  hist.samples = samples;

  if ((hist.vals = (double *) calloc(hist.bins, sizeof(double))) == NULL) {
      fprintf(stderr, "SetupHist: allocation failure, exit..\n");
      exit(1);
  }

  return hist;
}
/*=============================================================================
 GetNHist(): Retourne dans une structure 'Histo' les elements d'un histogramme
             reunis depuis 'samples' donnees du tableau pointee par 'scores'.
             Le tableau cree deborde (en general legerement) par rapport aux
             valeurs enregistrees:
             'dx' est la mesure d'une division de l'histogramme d'ou decoule le
             nombre d'intervalles (champ 'hist.bins').
             L'integrale de l'histogramme (pas la somme !) est normalisee a 1.
=============================================================================*/

Histo GetNHist(float *scores, int samples, double dx)
{
  Histo  hist;
  Map    map;
  double min, max, C;
  int    i;

  min = max = scores[0];

  for (i = 1; i < samples; i++) {
      if (scores[i] < min) min = scores[i];
      else
      if (scores[i] > max) max = scores[i];
  }

  hist = SetupHist(min, max, dx, samples);

  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins - 1);      /* etalonnage */

  for (i = 0; i < hist.samples; i++)  hist.vals[ GetX(map, scores[i]) ]++;

  C = (double) hist.samples / map.ratio;
  for (i = 0; i < hist.bins; i++)   hist.vals[i] /= C;  /* normalisation a 1 */
                                                           /* de l'integrale */
  return hist;
}
/*=============================================================================
 GetNHistW(): Retourne dans une structure 'Histo' les elements d'un histogramme
             reunis depuis 'samples' donnees du tableau pointe par 'scores' et
	     dont les poids associes sont donnes par le tableau pointe par
	     'weights'.
             Le tableau cree deborde (en general legerement) par rapport aux
             valeurs enregistrees:
             'dx' est la mesure d'une division de l'histogramme d'ou decoule le
             nombre d'intervalles (champ 'hist.bins').
             L'integrale de l'histogramme (pas la somme !) est normalisee a 1.
             Variante de 'GetHisto' de 'histo.c'.
=============================================================================*/

Histo GetNHistW(double *scores, double *weights, int samples, double dx)
{
  Histo  hist;
  Map    map;
  double min, max, C, total;
  int    i;

  min = max = scores[0];

  for (i = 1; i < samples; i++) {
      if (scores[i] < min) min = scores[i];
      else
      if (scores[i] > max) max = scores[i];
  }

  hist = SetupHist(min, max, dx, samples);

  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins - 1);      /* etalonnage */

  for (i = 0, total = 0.0; i < hist.samples; i++)  {
      hist.vals[ GetX(map, scores[i]) ] += weights[i];
      total += weights[i];
  }

  C = total / map.ratio;
  for (i = 0; i < hist.bins; i++)   hist.vals[i] /= C;  /* normalisation a 1 */
                                                           /* de l'integrale */
  return hist;
}
/*=============================================================================
 NormalizeHist(): Procede a la normalisation de l'histogramme 'hist'.
                  L'integrale (pas la somme !) de l'histogramme est amenee a 1.
=============================================================================*/

void NormalizeHist(Histo hist)
{
  Map    map;
  double sum;
  int    i;

  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins - 1);

  for (i = 0, sum = 0.0; i < hist.bins; i++)  sum += hist.vals[i];

  sum /= map.ratio;

  for (i = 0; i < hist.bins; i++) hist.vals[i] /= sum;

  return;
}
/*=============================================================================
 SetHistAxis(): Retourne un pointeur sur un tableau constituant l'axe des or-
                donnees de l'histogramme contenu dans 'hist'.
=============================================================================*/

double *SetHistAxis(Histo hist)
{
  Map    map;
  double *axis;
  int    i;

  if ((axis = (double *) malloc(hist.bins * sizeof(double))) == NULL) {
      fprintf(stderr, "SetHistAxis: allocation failure, exit..\n");
      exit(1);
  }
  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins-1);

  for (i = 0; i < hist.bins; i++) axis[i] = Getx(map, i);

  return axis;
}
/*=============================================================================
 PrHistInfo():
=============================================================================*/

void PrHistInfo(FILE *txt, char *cmt, Histo hist)
{
  fprintf(txt, "%s\n", cmt);
  fprintf(txt, "bins: %d\n", hist.bins);
  fprintf(txt, "hmin: %.3f\n", hist.hmin);
  fprintf(txt, "hmax: %.3f\n\n", hist.hmax);
  return;
}
/*=============================================================================
 PrHistDat():
=============================================================================*/

void PrHistDat(FILE *txt, Histo hist)
{
  Map  map;
  int  i;

  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins-1);

  for (i = 0; i < hist.bins; i++)
      fprintf(txt, "%.4e %.4e\n", Getx(map, i), hist.vals[i]);
  return;
}
/*=============================================================================
 GetHistStat(): Calcule la moyenne et l'ecart type de l'histogramme 'hist'.
                L'histogramme est suppose normalise (integrale = 1).
=============================================================================*/

void GetHistStat(Histo hist, double *mean, double *std)
{
  int    i;
  double score, prob, v;
  Map    map;

  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins - 1);

  for (i = 0, *mean = *std = 0.0; i < hist.bins; i++)
  {
      score = Getx(map, i);
      prob = hist.vals[i] / map.ratio;
      *mean += prob * score;
      *std  += prob * score * score;
  }
  v = *std - *mean * *mean;
  *std = (v < EPSILON ? 0.0 : sqrt(v));

  return;
}
/*=============================================================================
 Gauss(): Retourne la valeur au point 'x' de la gaussienne normalisee a 1,
          parametree par 'mean' et 'std' (moyenne et deviation standard).
=============================================================================*/

double Gauss(double x, double mean, double std)
{
  double C = 1.0 / sqrt(2.0 * M_PI) / std;

  x = (x - mean) / std;
  return (C * exp( - x * x / 2.0));
}
/*=============================================================================
 RandNucl(): retourne 0,1,2 ou 3 avec les probabilites passees en arguments.
            Cette version est adaptee au contexte de la fonction 'RndNoGScores'
=============================================================================*/

int RandNucl(double p0, double p1, double p2)
{
  double x, rnd = rand() / (RAND_MAX + 1.0);

  x = p0;
  if (rnd < x) return _A_;
  x += p1;
  if (rnd < x) return _T_;
  x += p2;
  if (rnd < x) return _G_;
  return _C_;
}
/*=============================================================================
 RndNucl(): retourne A,T,G ou C avec les probabilites passees en arguments.
            Cette version est adaptee au contexte de la fonction 'RndAlnScores'
=============================================================================*/

char RndNucl(double p0, double p1, double p2)
{
  double x, rnd = rand() / (RAND_MAX + 1.0);

  x = p0;
  if (rnd < x) return 'A';
  x += p1;
  if (rnd < x) return 'T';
  x += p2;
  if (rnd < x) return 'G';
  return 'C';
}
/*=============================================================================
 SortAsc(): Procede au tri ascendant, par insertion directe, des 'n' premiers
            elements du tableau d'entiers pointe par 'T'.
=============================================================================*/

void SortAsc(int *T, int n)
{
  int i, j, a;

  for (i = 1; i < n; i++)
  {
      a = T[i];
      for (j = i - 1; j >= 0 && a < T[j]; j--)  T[j + 1] = T[j];
      T[j + 1] = a;
  }
}
/*=============================================================================
 ResetHistBins(): Modifie l'intervalle elementaire de l'histogramme pointe par
                  'hist' pour le reduire a 'bins' points et recalcule la dist-
	          ribution en conservant sa normalisation (moyenne locale).
                  Si 'bins' est <= hist->bins (le nombre de points de l'histo-
		  gramme passe en argument) la fonction retourne 0 sinon 1.
		  Cette fonction est destinee a permettre diverses visualisa-
		  tions d'un histogramme.
=============================================================================*/

int ResetHistBins(Histo *hist, int bins)
{
  int    i;
  double ratio, *newvals;
  Map    map, newmap;

  if (bins >= hist->bins) return 0;     /* augmentation du nb de pts refusee */

  ratio = (double) hist->bins / (double) bins;
  map = MapSetup(hist->hmin, hist->hmax, 0, hist->bins - 1);
  newmap = MapSetup(hist->hmin, hist->hmax, 0, bins - 1);

  if ((newvals = (double *) calloc(bins, sizeof(double))) == NULL) {
      fprintf(stderr, "ResetHistBins: allocation failure, exit..\n");
      exit(1);
  }

  for (i = 0; i < hist->bins; i++) {
      double x = Getx(map, i);
      newvals[GetX(newmap, x)] += hist->vals[i];
  }
  free(hist->vals);
  hist->vals = newvals;
  hist->bins = bins;
  for (i = 0; i < hist->bins; i++) hist->vals[i] /= ratio;

  return 1;
}
/*===========================================================================*/

#undef EPSILON
