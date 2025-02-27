
/*=============================================================================
 dhisto.c                  A.Lambert le 24/02/03    revu le 11/04/03

 fonctions permettant la saisie de l'histogramme non normalise des scores des
 detections d'un masque a l'issue d'une operation de recherche par 'erpin',
 (on y cumule simplement les detections successives).

 La fonction 'PrintScoresHisto' est executee a la sortie de 'erpin', le calcul
 de l'histogramme est effectue seulement dans le cas de l'option '-hist'.

 cc -O2 -Wall -c dhisto.c -I../include ;
 ar -rs ../lib/librnaIV.a dhisto.o ;

=============================================================================*/

#include "rnaIV.h"

/*
voir "rnaIV.h"

typedef struct {
  double  min, max, ratio;
  int     Pixmin, Pixmax;
  } Map;

typedef struct {
  double  hmin, hmax;
  int     bins, samples;
  double  *vals;
  } Histo;
*/

void   InitDetectsHisto(Context *ctxt);
void   AddToDetectsHisto(double x);
void   PrintHisto(Histo hist, char *filename);
Histo  ReadHisto(char *filename);
void   PrintScoresHisto(Context *ctxt);

/*=============================================================================
 InitDetectsHisto(): Initialise l'histogramme des detections gerees par 'ctxt'.
                     fonction appelee par 'InitOutputContext' de 'msearch.c'
=============================================================================*/

void InitDetectsHisto(Context *ctxt)
{
  extern Histo MainMaskDetects;
  extern Histo MainMaskEvals;
  extern Map   MainMaskMap;
  
  MainMaskDetects = SetupHist(ctxt->mask->threshold, MainMaskEvals.hmax,
                              DELTA_H, 0);
  MainMaskMap = MapSetup(MainMaskDetects.hmin, MainMaskDetects.hmax,
                         0, MainMaskDetects.bins - 1);
  return;
}
/*=============================================================================
 AddToDetectsHisto(): Ajoute a l'histogramme 'MainMaskDetects' qui est suppose
                      correctement initialise ainsi que 'MainMaskMap' la valeur
                      detecte 'x'.
                      fonction appelee par 'PrintOutput' de 'outputs.c'
=============================================================================*/

void AddToDetectsHisto(double x)
{
  extern Histo MainMaskDetects;
  extern Map   MainMaskMap;
  int index;

  if (x > MainMaskDetects.hmin && x < MainMaskDetects.hmax)
  {
      index = GetX(MainMaskMap, x);
      MainMaskDetects.vals[index]++;
      MainMaskDetects.samples++;
  }
  return;
}
/*=============================================================================
 PrintHisto(): Enregistre le contenu de l'histogramme pointe par 'hist'.
=============================================================================*/

void PrintHisto(Histo hist, char *filename)
{
  FILE  *txt;
  int   i;

  txt = fopen(filename, "w");

  fprintf(txt, " %.4f  %.4f  %d  %d\n",
          hist.hmin, hist.hmax, hist.bins, hist.samples);

  for (i = 0; i < hist.bins; i++)  fprintf(txt, " %d\n", (int) hist.vals[i]);

  fclose(txt);
  return;
}
/*=============================================================================
 ReadHisto(): Lit dans un fichier les donnees d'un histogramme et les retourne
              dans une structure 'Histo'.
=============================================================================*/

Histo ReadHisto(char *filename)
{
  Histo  hist;
  FILE   *txt;
  int    i;

  if ((txt = fopen(filename, "r")) == NULL) {
      fprintf(stderr, "ReadHist: file '%s' not found, exit..\n", filename);
      exit(1);
  }

  fscanf(txt, "%lf %lf %d %d",                      /* 1ere ligne du fichier */
         &hist.hmin, &hist.hmax, &hist.bins, &hist.samples);

  hist.vals = (double *) malloc(hist.bins * sizeof(double));

  for (i = 0; i < hist.bins; i++)
      fscanf(txt, "%lf", &hist.vals[i]);

  fclose(txt);

  return hist;
}
/*=============================================================================
 PrintScoresHisto(): Interface des fonctions precedentes.
                     L'histogramme des scores n'est dresse que dans le cas ou
                     l'option '-hist' (erpin) a ete entree.
                     On s'assure avant la creation de l'histogramme de la pre-
                     sence d'un nombre suffisant de points 'MIN_DETECTS'.
=============================================================================*/

void PrintScoresHisto(Context *ctxt)
{
  extern Histo MainMaskDetects;

  if (ctxt->hist == ON)
  {
      if (MainMaskDetects.samples > MIN_DETECTS)
      {
          PrintHisto(MainMaskDetects, HISTO_FNAME);
          free(MainMaskDetects.vals);
          fprintf(ctxt->outfile,
                  "writing histogram of hits to '%s' file..\n", HISTO_FNAME);
      }
      else
          fprintf(ctxt->outfile,
                  "too few hits for writing histogram..\n");
  }
  return;
}
/*===========================================================================*/
