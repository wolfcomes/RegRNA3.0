
/*=============================================================================
 lfit.c                                           A.Lambert le 01/04/03

 Methode des moindres carres, fit lineaire sur un tableau de points.

 cc -O2 -Wall -c lfit.c -I../include ;
 ar -rs ../lib/librnaIV.a lfit.o ;

=============================================================================*/

#include "rnaIV.h"

/* structure servant a stocker les coord. d'un point du plan

typedef struct {
   double x, y;
} Point;

*/

void fit(Point *Pts, int ndata, double *a, double *b,
         double *siga, double *sigb, double *chi2);
void lfit(Point *Pts, int ndata, double *a, double *b);

/*=============================================================================
 fit(): Effectue un ajustage de 'ndata' points x,y avec la droite y = a + b*x
        utilisant la methode des moindres carres (determination de a et b par
        minimisation de chi2).
        En retour les parametres de la droite sont pointes par 'a', 'b', les
        deviations standards par 'siga', 'sigb'. 'chi2' est aussi retourne.
=============================================================================*/

void fit(Point *Pts, int ndata, double *a, double *b,
         double *siga, double *sigb, double *chi2)
{
  int   i;
  double meanx, sumx, sumy, t, sumt2, sigdat, dy;

  *b = 0.0;

  for (i = 0, sumx = 0.0, sumy = 0.0; i < ndata; i++) {
     sumx += Pts[i].x;
     sumy += Pts[i].y;
  }

  meanx = sumx / (double)  ndata;
  for (i = 0, sumt2 = 0.0; i < ndata; i++) {
     t = Pts[i].x - meanx;
     sumt2 += t * t;
     *b += t * Pts[i].y;
  }

  *b /= sumt2;
  *a = (sumy - sumx*(*b)) / ndata;
  *siga = sqrt(( 1.0 + sumx*sumx / (ndata*sumt2) ) / ndata);
  *sigb = sqrt(1.0/sumt2);

  for (i = 0, *chi2 = 0.0; i < ndata; i++) {
     dy = Pts[i].y - (*a) - (*b) * Pts[i].x;
     *chi2 += dy*dy;
  }

  sigdat = sqrt((*chi2) / (ndata - 2));
  *siga *= sigdat;
  *sigb *= sigdat;

  return;
}
/*=============================================================================
 lfit(): Interface simplifiee de 'fit'.
         En retour 'a' et 'b' pointent les coefficients de la droite 'a + b*x'.
=============================================================================*/

void lfit(Point *Pts, int ndata, double *a, double *b)
{
  double siga, sigb, chi2;

  fit(Pts, ndata, a, b, &siga, &sigb, &chi2);
  return;
}
/*===========================================================================*/
