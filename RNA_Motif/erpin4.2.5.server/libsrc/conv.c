
/*=============================================================================
 conv.c                                            A.Lambert le 20/03/03

 Ce code contient les fonctions necessaires au calcul du produit de convolution
 de 2 vecteurs.
 On n'utilise pas ici la fft, procedant a un calcul direct.
 calcul direct du produit de convolution de 2 vecteurs.

 cc -O2 -Wall -c conv.c -I../include ;
 ar -rs ../lib/librnaIV.a conv.o ;

=============================================================================*/

#include "rnaIV.h"

double *convlv(double *a, double *b, int na, int nb, int *ncnv);
void   conv(double **a, int *na, double *b, int nb);
void   ConvHist(Histo *h1, Histo h2);

/*=============================================================================
 convlv(): Calcule de facon directe le produit de convolution de 2 tableaux de
           type 'double' pointes par 'a' et 'b' de longueurs respectives 'na'
           et 'nb'.
           La fonction retourne un pointeur sur un tableau dans lequel les ele-
           ents du produit de convolution 'a*b' sont enregistres, la taille de
           ce tableau est 'na + nb -1', 'ncnv' pointe en sortie cette valeur.
=============================================================================*/

double *convlv(double *a, double *b, int na, int nb, int *ncnv)
{
  double *cnv, x;
  int    j, k, s, min, max;

  *ncnv = na + nb - 1;        /* taille du tableau du produit de convolution */

  if ((cnv = (double *) malloc(*ncnv * sizeof(double))) == NULL) {
      fprintf(stderr, "convlv: allocation failure, exit..\n");
      exit(1);
  }

  for (k = 0; k < *ncnv; k++)
  {
      min = (s = k+1-nb) > 0 ? s : 0;
      max = (s = na-1) < k ? s : k;
      x = 0.0;
      for (j = min; j <= max; j++)  x += a[j] * b[k-j];
      cnv[k] = x;
  }
  return cnv;
}
/*=============================================================================
 conv(): interface de 'convlv' adapte a des appels iteres.
         Le resultat est pointe par 'a' dont la taille est aussi modifiee,
         le contenu pointe par 'a' en entree est donc ecrase.
         En entree 'na' pointe la taille du 1er tableau, et en sortie celle du
         produit de convolution calcule.
=============================================================================*/

void conv(double **a, int *na, double *b, int nb)
{
  double  *cnv, *p1, *p2;
  int     i, ncnv;

  cnv = convlv(*a, b, *na, nb, &ncnv);

  free(*a);

  if ((*a = (double *) malloc(ncnv * sizeof(double))) == NULL) {
      fprintf(stderr, "conv: allocation failure, exit..\n");
      exit(1);
  }
                                                 /* copie de 'cnv' dans '*a' */
  for (i = 0, p1 = *a, p2 = cnv; i < ncnv; i++, p1++, p2++) *p1 = *p2;

  free(cnv);
  *na = ncnv;
  return;
}
/*=============================================================================
 ConvHist(): Procede au calcul du produit de convolution de deux structures
             'Histo': 'h1' pointee par 'h1' et 'h2' passee par valeur.
             Le resultat est en retour pointe par 'h1'.
             Cette fonction sera utilise pour calculer le produit de convolu-
             tion de 2 fonctions discretisees avec un pas commun.
=============================================================================*/

void ConvHist(Histo *h1, Histo h2)
{

  conv(&(h1->vals), &(h1->bins), h2.vals, h2.bins);

  h1->hmin += h2.hmin;
  h1->hmax += h2.hmax;

  return;
}
/*===========================================================================*/
