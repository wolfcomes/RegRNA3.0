
/*=============================================================================
 cdf.c                              A.Lambert le 12/04/03 revu le 05/05/03

 "cumulative distribution function"

 cc -O2 -Wall -c cdf.c -I../include ;
 ar -rs ../lib/librnaIV.a cdf.o ;

=============================================================================*/

#include "rnaIV.h"

Histo  GetHistCdf(Histo hist);
double Interpol(Histo cdf, double x);

/*=============================================================================
 GetHistCdf(): Retourne une structure 'Histo' donnant la "cdf" de l'histogramme
               'hist'.
               Si l'integrale de 'hist' est normalisee a 1 chaque element de 
               'cdf.vals' donne l'aire de la partie de 'hist' situee a droite
               de l'element considere.
=============================================================================*/

Histo GetHistCdf(Histo hist)
{
  Map    map;
  Histo  cdf;
  int    i;
  double dh;

  if (hist.bins < 2) {
      fprintf(stderr, "GetHistCdf: too few elements in array, exit..\n");
      exit(0);
  }

  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins - 1);
  dh = 1.0 / map.ratio;          /* largeur d'un intervalle de l'histogramme */

  cdf.hmin = hist.hmin;
  cdf.hmax = hist.hmax;
  cdf.bins = hist.bins;
  if ((cdf.vals = (double *) malloc(cdf.bins * sizeof(double))) == NULL) {
      fprintf(stderr, "GetHistCdf: allocation failure, exit..\n");
      exit(1);
  }

  cdf.vals[cdf.bins-1] = hist.vals[cdf.bins-1];

  for (i = cdf.bins - 2; i >=0; i--)
      cdf.vals[i] = cdf.vals[i+1] + hist.vals[i];

  for (i = 0; i < cdf.bins; i++)  cdf.vals[i] *= dh;

  return cdf;
}
/*=============================================================================
 Interpol(): Retourne la valeur obtenue par interpolation lineaire des elements
             du tableau contenu dans'cdf' correspondant a l'abscisse 'x'.
             
             Si l'histogramme dont est deduit 'cdf' a une integrale normalisee
             a 1 et represente la probabilite associee a une variable aleatoire
             S cette fonction retourne la probabilite que S soit superieur a
             l'argument 'x'.
             Une interpolation lineaire est utilisee.
             L'interpolation est prolongee inferieurement et superieurement par
             les valeurs initiale et finale du tableau pointe par 'cdf.vals'.
=============================================================================*/

double Interpol(Histo cdf, double x)
{
  Map    map;
  int    Xo;
  double dh, xo, E;

  if (x < cdf.hmin) return cdf.vals[0];                     /* prolongements */
  if (x > cdf.hmax) return cdf.vals[cdf.bins - 1];  /* si x est hors limites */

  map = MapSetup(cdf.hmin, cdf.hmax, 0, cdf.bins - 1);
  dh = 1.0 / map.ratio;          /* largeur d'un intervalle de l'histogramme */

  Xo = GetXFloor(map, x);
  xo = Getx(map, Xo);

  E = cdf.vals[Xo] + (x - xo) * (cdf.vals[Xo+1] - cdf.vals[Xo]) / dh;

  return E;
}
/*===========================================================================*/
