
/*=============================================================================
 epnstat.c                          A.Lambert le 17/04/03 revu le 14/05/04

 Ce programme permet le traitement de l'histogramme, lu dans un fichier, des
 scores superieurs a un seuil issu de l'option '-hist' de 'erpin'.
 Un des buts est de comparer ces resultats "experimentaux" a ceux produits par
 'mstat', particulierement a la courbe des E-values vs scores.

 Plutot que l'histogramme cumule (CDF) qui est produit par defaut, on peut, en
 utilisant l'option '-hist', coder les donnees de l'histogramme original (PDF)
 sous une forme directement utilisable par 'gnuplot', 'octave' ou 'matlab'.
 L'option '-bins' permet de modifier le nombre d'intervalles elementaires de
 l'histogramme.

 Les resultats sont diriges sur la sortie standard.

 cc -Wall -O3 -o epnstat epnstat.c -I../include -L../lib -lrnaIV -lm ;
 strip epnstat ;
 chmod 755 epnstat ;
 mv epnstat ../bin ;

 Usage:   epnstat <filename>       fichier de l'histogramme produit par 'erpin'
                  [-hist]         donnees non cumulees de histogramme original
		  [-bins <bins>]                  fixe le nombre d'intervalles

 Exemple: epnstat data/epn.trna.hist.dat > data/epnstat.trna.dat

 Exemple complet concernant sno-CD-archae.epn:
 ---------------------------------------------
 cd ev.pub/sno ;
 ../../bin/epnstat epn.sno.hist.dat > epnstat.sno.dat ;
 ../../bin/mstat ~/devc/projets/bio/data/trsetsIII/sno-CD-archae.epn -2,2 \
                 -mask 5 7 9 -logzero -100 -cutoff 18.43 -Mbs 300 -evals ;
 mv mlogev.dat mlogev.sno.dat ;
 mv mevals.dat mevals.sno.dat ;

 gnuplot ;
 set title "Erpin & mstat: sno-CD-archae.epn, 300Mb data" ;
 set grid ;
 set nokey ;
 set xlabel "region {-2,2 - {5 7 9}} scores" ;
 set ylabel "Log(E-value)" ;
 set xrange [18.43:26] ;
 set size square ;
 logscale 'epnstat.sno.dat' , 'mevals.sno.dat' u 1:2 with lines ;

 set terminal postscript portrait ;
 set output 'epnstat.sno.ps' ;
 replot ;

 Exemple complet concernant trna-typeI.epn:
 ------------------------------------------
 cd ev.pub/trna6.8.13 ;
 ../../bin/epnstat epn.trna.hist.dat > epnstat.trna.dat ;
 ../../bin/mstat ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -6,8 \
                 -umask 6 8 13 -cutoff 15 -logzero -50 -Mbs 300 -evals ;
 mv mlogev.dat mlogev.trna.dat ;
 mv mevals.dat mevals.trna.dat ;

 gnuplot ;
 set title "Erpin & mstat: trna-typeI.epn, 300Mb data" ;
 set grid ;
 set nokey ;
 set xlabel "region {6 8 13} scores" ;
 set ylabel "Log(E-value)" ;
 set xrange [15:21] ;
 set size square ;
 logscale 'epnstat.trna.dat' , 'mevals.trna.dat' u 1:2 with lines ;

 set terminal postscript portrait ;
 set output 'epnstat.trna.ps' ;
 replot ;

=============================================================================*/

#include "rnaIV.h"

#define SMALL    (1.e-20)
#define PDF 0
#define CDF 1

Histo GetEpnCdf(Histo hist);

/*=============================================================================
 main():
=============================================================================*/

int main(int argc, char *argv[])
{
  Histo  hist, cdf;
  Map    map;
  int    i, bins, mode = CDF;

  if (argc < 2) {
      fprintf(stderr, "Usage: epnstat <filename> [-hist][-bins <bins>]\n");
      exit(1);
  }
  hist = ReadHisto(argv[1]);   /* histogramme brut des detections de 'erpin' */
  bins = hist.bins;                                     /* valeur par defaut */

  for (i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-hist") == 0)  mode = PDF;
      if (strcmp(argv[i], "-bins") == 0)  bins = atoi(argv[++i]);
  }

  if (bins < hist.bins)
      ResetHistBins(&hist, bins);

  if (mode == PDF)                       /* simple recodage de l'histogramme */
  {
      map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins - 1);

      for (i = 0; i < hist.bins; i++) {
          printf("%.3e %.4e\n", Getx(map, i), hist.vals[i]);
      }
      free(hist.vals);
      return 0;
  }
                                                  /* histogramme cumule: CDF */
  cdf = GetEpnCdf(hist);
  free(hist.vals);

  map = MapSetup(cdf.hmin, cdf.hmax, 0, cdf.bins - 1);

  for (i = 0; i < cdf.bins; i++)
  {
      double x, y;
      x = Getx(map, i);
      y = cdf.vals[i];              /* 'cdf.vals' est un tableau decroissant */
      if (y < SMALL) break;            /* eviter ulterieurement les log de 0 */
      printf("%.3e %.4e\n", x, y);
  }

  free(cdf.vals);
  return 0;
}
#undef SMALL
/*=============================================================================
 GetEpnCdf(): Dresse a partir de l'histogramme brut 'hist' celui des detections
              cumulees conduisant a la 'E-value', nombre moyen des detections
              de score superieur au seuil donne par l'abscisse.
              Variante de 'GetHistCdf' de 'cdf.c'.
=============================================================================*/

Histo GetEpnCdf(Histo hist)
{
  Map    map;
  Histo  cdf;
  int    i;

  if (hist.bins < 2) {
      fprintf(stderr, "GetEpnCdf: too few elements in array, exit..\n");
      exit(1);
  }

  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins - 1);

  cdf.hmin = hist.hmin;
  cdf.hmax = hist.hmax;
  cdf.bins = hist.bins;
  if ((cdf.vals = (double *) malloc(cdf.bins * sizeof(double))) == NULL) {
      fprintf(stderr, "GetEpnCdf: allocation failure, exit..\n");
      exit(1);
  }

  cdf.vals[cdf.bins-1] = hist.vals[cdf.bins-1];

  for (i = cdf.bins - 2; i >= 0; i--)
      cdf.vals[i] = cdf.vals[i+1] + hist.vals[i];

  return cdf;
}
/*===========================================================================*/
