
/*=============================================================================
 bkgstat.c                                     A.Lambert le 07/04/03

 Complements sur la statistique du "fond" de sequences aleatoires cibles.
 La fonction 'ChMaskStat' qui concerne un masque isole ne devrait pas etre
 utilisee par 'erpin' qui travaille plutot avec une region entiere (pattern).

 cc -O2 -Wall -c bkgstat.c -I../include ;
 ar -rs ../lib/librnaIV.a bkgstat.o ;

=============================================================================*/

#include "rnaIV.h"

void ReadBkgFreqs(double *freqs, int offset, char *argv[]);
void ResetBkgFreqs(double *freqs);
void ChMaskStat(Mask *mask);

/*=============================================================================
 ReadBkgFreqs(): Lit dans la liste pointee par 'argv' des arguments du program-
                 me, a partir de l'indice 'offset+1' les frequences de "fond"
                 des bases, qui seront stockees dans le tableau 'freqs'.
                 Seulement 3 nombres sont saisis, le 4eme en sera deduit.
                 Les nombres saisis sont des pourcentages, divises ensuite par
                 100 pour normalisation a 1.
=============================================================================*/

void ReadBkgFreqs(double *freqs, int offset, char *argv[])
{
  int    i;
  double sum;

  for (i = 0, sum = 0.0; i < ALPHA_LEN - 1; i++)
  {
      freqs[i] = atof(argv[++offset])/100.0;
      sum += freqs[i];
  }
  if (sum > 1.) {
      fprintf(stderr, "invalid arguments for 'bkg' option, exit..\n");
      exit(1);
  }
  freqs[ALPHA_LEN - 1] = 1.0 - sum;

  return;
}
/*=============================================================================
 ResetBkgFreqs(): reinitialise les tableaux fixant la statistique du 'fond',
                  a l'aide du tableau pointe par 'freqs'. Les bases successi-
                  sont supposees statitiquement independantes.
=============================================================================*/

void ResetBkgFreqs(double *freqs)
{
  extern double *NewLogDataFreqs;
  int i;

  for (i = 0; i < ALPHA_LEN; i++)  NewLogDataFreqs[i] = log(freqs[i]);

  return;
}
/*=============================================================================
 ChMaskStat(): Interface, pour la structure pointee par 'mask', des fonctions
               relatives a la modification de statistique du 'fond' pour ses
               elements constitutifs: brins et helices.
               L'operation s'acheve par la permutation des tableaux contenant
               les donnees statistiques (en vue d'une nouvelle modification).
=============================================================================*/

void ChMaskStat(Mask *mask)
{
  int  i, *m;

  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++) {
      ChStStat(mask->pattern->std + *m);
  }
  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++) {
      ChHlxStat(mask->pattern->hlx + *m);
  }
  SwapStatTables();

  return;
}
/*===========================================================================*/
