
/*=============================================================================
 thresholds.c                        A.Lambert le 27/12/01

 Fonctions destinees a lire les arguments 'cutoff' entres a un programme sous
 forme de nombres servant de seuils ou sous forme de pourcentages de captures
 dans une base d'entrainement.

 cc -O2 -Wall -c thresholds.c -I../include ;

 ar -r ../lib/librnaIV.a thresholds.o ;

=============================================================================*/

#include "rnaIV.h"

/*
 rappel de 'rnaIV.h'

typedef struct {
      int    input;                           type d'entree: VALUE ou PERCENT
      int    percent;                 pourcentage de captures dans un 'trset'
      double val;                   valeur du seuil entre ou apres convertion
  }
  Threshold;
*/

Threshold *GetThresholdsArgs(int argc, char *argv[], int nmask);
void  GetThresholds(Threshold *threshold, Trset *trset, Mask *mask, int nmask);
void  fPrintfThresholds(Mask *mask, int nmask, FILE *txt);

/*=============================================================================
 GetThresholdsArgs(): Saisie des arguments 'cutoff' dans les arguments de
                      'main', lesquels seront eventuellement convertis ulteri-
                      eurement.
                      Les seuils sont entres soit par leur valeur, soit par un
                      pourcentage (entier) de captures, dans ce cas la chaine
                      de l'argument se termine par le caractere '%'.
                      Prealablement a la lecture, on affecte a chaque seuil
                      une valeur par defaut egale a 100%.
                      L'option '-cutoff' est suivi d'un nombre de seuils egal
                      ou inferieur au nombre de masques 'nmask' qui doit etre
                      prealablement determine, sinon ou en cas d'argument mal
                      forme, la valeur par defaut (100%) est conservee.
                      L'ordre des seuils lus suit celui des masques.
                      Note: lire les commentaires dans le code.
=============================================================================*/

Threshold *GetThresholdsArgs(int argc, char *argv[], int nmask)
{
  int       i, j, len;
  char      status, *str = (char *) malloc(12);
  Threshold *threshold;
  double    x;

  threshold = (Threshold *) malloc(nmask * sizeof(Threshold));

  for (i = 0; i < nmask; i++) {                        /* valeurs par defaut */
      threshold[i].percent = 100;
      threshold[i].input = PERCENT;
  } 

  for (i = 0; i < argc-1; i++)      /* argc-1: ignore "-cutoff" s'il termine */
  {                                                  /* la ligne de commande */
      if (strcmp(argv[i], "-cutoff") == 0)
      {
          for (j = 0; j < nmask && i < argc-1; j++)   /* argc-1: incrementer */
          {                                       /* l'index 'i' sans danger */
              strcpy(str, argv[++i]);
              len = strlen(str);
              if (str[len - 1] != '%')
              {
                  x = StrToD(str, &status);
                  if (status == 1) {              /* enregistre la valeur si */
                      threshold[j].val = x;         /* elle est correctement */
                      threshold[j].input = VALUE;                /* detectee */
                  }
                  else break;  /* arg. mal forme ou option suivante atteinte */
              }
              else {
                  str[len - 1] = '\0';             /* le car. '%' est ecrase */
                  x = StrToD(str, &status);
                  if (status == 1) {            /* enregistre le pourcentage */
                      threshold[j].percent = (int) x;
                      threshold[j].input = PERCENT;
                  }
              }
          }
          break;                           /* une seule option "-cutoff" lue */
      }
  }
  free(str);

  return threshold;
}
/*=============================================================================
 GetThresholds(): Initialise les champs 'threshold' d'un tableau de masques
                  apres convertion eventuelle des seuils entres sous forme de
                  pourcentages de captures dans un 'trset'.
=============================================================================*/

void GetThresholds(Threshold *threshold, Trset *trset, Mask *mask, int nmask)
{
  int   i;

  for (i = 0; i < nmask; i++)
  {
      if (threshold[i].input == PERCENT)
          GetMaskThreshold(threshold[i].percent, trset, mask + i);
      else
          mask[i].threshold = threshold[i].val;
  }

  return;
}
/*=============================================================================
 fPrintfThresholds(): Affiche les valeurs des seuils associes a un tableau de
                      'nmasks' masques.
=============================================================================*/

void fPrintfThresholds(Mask *mask, int nmask, FILE *txt)
{
  int i;

  fprintf(txt, "Cutoff:\t\t");

  for (i = 0; i < nmask; i++)
      fprintf(txt, "%.2f  ", mask[i].threshold);

  fprintf(txt, "\n\n");

  return;
}
/*===========================================================================*/
