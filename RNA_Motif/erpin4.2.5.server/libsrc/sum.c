
/*=============================================================================
 sum.c                                          A.Lambert le 15/04/04
 
 Fonctions concernant la lecture et l'utilisation d'une matrice de substitution.
 
 cc -O2 -Wall -c sum.c -I../include ;
 
 ar -rs ../lib/librnaIV.a sum.o ;

=============================================================================*/

#include "rnaIV.h"

void ReadSUM(char *fname, double hpcw, double spcw, Trset *trset);
void SumImg(double **Sbstm, int dim, double **prof, int len, double weight);

/*=============================================================================
 ReadSUM(): lit dans le fichier identifie par 'fname' le contenu de 2 matrices
            de substitution correspondant a des helices et brins de cotes res-
            pectifs SQR_ALPHA_LEN et ALPHA_LEN, et les charge dans 2 tableaux
            pointes par 'trset->hlxsum' et 'trset->stsum'.
	    Les poids 'hpcw' et 'spcw' des pseudo-comptes sont enregistres.
=============================================================================*/

void ReadSUM(char *fname, double hpcw, double spcw, Trset *trset)
{
  int  i, j;
  FILE *txt;

  if ((txt = fopen(fname, "r")) == NULL)  {
      fprintf(stderr, "ReadSUM: file '%s' not found, exit...\n", fname);
      exit(1);
  }
  trset->pcflag = ON;
  trset->hpcw = hpcw;
  trset->spcw = spcw;

  trset->hlxsum = dMat(SQR_ALPHA_LEN, SQR_ALPHA_LEN);          /* allocation */
  trset->stsum  = dMat(ALPHA_LEN, ALPHA_LEN);

  for (i = 0; i < SQR_ALPHA_LEN; i++)        /* matrice relative aux helices */
      for (j = 0; j < SQR_ALPHA_LEN; j++)
          fscanf(txt, "%lf", trset->hlxsum[i] + j);

  for (i = 0; i < ALPHA_LEN; i++)              /* matrice relative aux brins */
      for (j = 0; j < ALPHA_LEN; j++)
          fscanf(txt, "%lf", trset->stsum[i] + j);

  fclose(txt);
}
/*=============================================================================
 SumImg(): Calcule l'image d'un profil d'helice ou de brin, suivant la dimen-
           sion 'dim' de la matrice de substitution pointee par 'Sbstm',
	   constitue de 'len' colonnes.
           Le parametre 'weight' donne le poids (entre 0 et 1) des pseudo-comp-
	   tes.
=============================================================================*/

void SumImg(double **Sbstm, int dim, double **prof, int len, double weight)
{
  int     i, j, k;
  double  *tmp;

  if (weight < 0 || weight > 1) {
      fprintf(stderr, "SumImg: arg #4 must be > 0 and < 1, exit..\n");
      exit(1);
  }

  tmp = (double *) malloc(dim * sizeof(double));

  for (j = 0; j < len; j++)             /* boucle sur les colonnes du profil */
  {
      FilldTab(tmp, dim, 0.0);
      for (i = 0; i < dim; i++)
          for (k = 0; k < dim; k++)  tmp[i] += Sbstm[i][k] * prof[k][j];

      for (i = 0; i < dim; i++)
          prof[i][j] = (1. - weight) * prof[i][j] + weight * tmp[i];
  }
  free(tmp);
  return;
}
/*===========================================================================*/
