
/*=============================================================================
 stsum.c                             A.Lambert le 22/03/04 revu le 31/12/04

 Ce code regroupe les fonctions procedant a la construction d'une matrice de
 substitution pour les bases impliquees dans les simples brins de la struc-
 ture des ARN regroupees dans un alignement multiple, lequel est utilise comme
 base d'entrainement.

 cc -O2 -Wall -c stsum.c -I../include ;

=============================================================================*/

#include "rnaIV.h"

double **GetWeights3(char **data, int nbstr, int bgn, int len, double *sweights);
void   GetStProfile3(Trset *trset, Strand *St, double *sweights);
void   GetMaskStProfiles3(Trset *trset, Mask *mask, double *sweights);
void   FreeMaskStProfiles(Mask *mask);
double **StSuMat(Trset *trset, Mask *mask, double *sweights, int Vmode);


void PrintdMat(double **m, int nrow, int ncol, FILE *txt);

/*=============================================================================
 PrintdMat(): Enregistre le contenu de la matrice 'm' de doubles dans le
             fichier 'txt'. ('stdout' pour un affichge au shell).
=============================================================================*/

void PrintdMat(double **m, int nrow, int ncol, FILE *txt)
{
  int    i, j;

  for (i = 0; i < nrow; i++)
  {
      for (j = 0; j < ncol; j++)   fprintf(txt, "%6.1f ", m[i][j]);
      fprintf(txt, "\n");
  }
  fprintf(txt, "\n");

  return;
}
/*=============================================================================
 GetWeights3(): Mesure sur les sequences contenues dans un tableau 'data' la
               frequence des caracteres rencontres (A, T, G, C, - ) a partir de
               'bgn' et sur une longueur 'len', les 'N' sont ignores.
               Stocke les resultats dans un profil dont le nombre de lignes est
               egal a la longueur de l'alphabet utilise plus 1 pour la frequence
               des gaps.
               Les resultats sont retournes dans un tableau de doubles.
               Remarque:
               Cette fonction ne prend pas de structure 'Strand' comme argument
               car elle est aussi destinee a etre utilisee sur des fragments de
               sequences non structures.
	       'sweights' pointe le tableau de ponderation des sequences, il
	       est suppose alloue et correctement initialise.
=============================================================================*/

double **GetWeights3(char **data, int nbstr, int bgn, int len, double *sweights)
{
  double **profile;
  int    i, j, height;

  height = ALPHA_LEN + 1;                 /* hauteur du profil: nb de lignes */

  profile = dMat(height, len);
  FilldMat(profile, height, len, 0.0);

  for (j = 0; j < len; j++)           /* boucle sur les colonnes successives */
  {
      for (i = 0; i < nbstr; i++)                   /* boucle sur les lignes */
          switch (data[i][j + bgn]) {
              case 'A': profile[_A_][j] += sweights[i]; break;
              case 'T': profile[_T_][j] += sweights[i]; break;
              case 'G': profile[_G_][j] += sweights[i]; break;
              case 'C': profile[_C_][j] += sweights[i]; break;
              case '-': profile[_X_][j] += sweights[i]; break;
              default : break;
          }
  }
  return profile;
}
/*=============================================================================
 GetStProfile3(): Interface de 'GetWeights3' prenant pour argument les adresses
                 de 'trset' et d'une structure 'Strand'.
                 Le profil n'est cree que si St->Profile == NULL.
		 'sweights' pointe le tableau de ponderation des sequences.
=============================================================================*/

void GetStProfile3(Trset *trset, Strand *St, double *sweights)
{
  if (St->Profile == NULL)
      St->Profile = GetWeights3(trset->data, trset->nseq,
                                St->db_bgn, St->max_len, sweights);
  else {
    fprintf(stderr, "GetStProfile3: allocation failure, exit..\n");
    exit(1);
  }
  return;
}
/*=============================================================================
 GetMaskStProfiles3(): Cree l'ensemble des profils des brins d'un masque de
                       pattern (limite au calcul du decompte des couples).
=============================================================================*/

void GetMaskStProfiles3(Trset *trset, Mask *mask, double *sweights)
{
  int i, *m;

  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
  {
      GetStProfile3(trset, mask->pattern->std + *m, sweights);
  }
  return;
}
/*=============================================================================
 FreeMaskStProfiles(): Libere la memoire allouee aux profils d'un Mask.
=============================================================================*/

void FreeMaskStProfiles(Mask *mask)
{
  int    i, *m;
  Strand *Std;

  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
  {
      Std = mask->pattern->std + *m;
      if (Std->Profile != NULL) {
          FreedMat(Std->Profile);
          Std->Profile = NULL;
      }
  }
  return;
}
/*=============================================================================
 StSuMat(): Calcule la matrice de substitution des brins du masque pointe par
            'mask', et retourne un pointeur sur ce tableau carre de cote
	    ALPHA_LEN.
	    'sweights' pointe le tableau de ponderation des sequences.
	    Si 'Vmode' vaut 1 (ou ON) le profil des brins est affiche.
	    La matrice produite est normalisee de sorte que la somme des ele-
	    ments de chaque colonne soit egale a 1.
=============================================================================*/

double **StSuMat(Trset *trset, Mask *mask, double *sweights, int Vmode)
{
  double **Stprof, **SbstMat, total;
  int    i, j, k, width;

  GetMaskStProfiles3(trset, mask, sweights);
  Stprof = CatStProfiles(mask, &width);

  if (Vmode == ON) PrintdMat(Stprof, ALPHA_LEN+1, width, stdout);

  printf("Stds: %d x %d\n", trset->nseq, width);

  FreeMaskStProfiles(mask);

       /* --- Creation et calcul de la matrice de substitution des brins --- */

  SbstMat = dMat(ALPHA_LEN, ALPHA_LEN);
  FilldMat(SbstMat, ALPHA_LEN, ALPHA_LEN, 0.0);

  for (i = 0; i < ALPHA_LEN; i++)
      for (j = 0; j < ALPHA_LEN; j++)
          for (k = 0; k < width; k++)
              SbstMat[i][j] += Stprof[i][k] * Stprof[j][k];

  for (j = 0; j < ALPHA_LEN; j++)
  {
      total = 0.0;                 /* somme des elts de la col. j de SbstMat */
      for (i = 0; i < ALPHA_LEN; i++) total += SbstMat[i][j];
                                        /* et normalisation a 1 de la col. j */
      for (i = 0; i < ALPHA_LEN; i++) SbstMat[i][j] /= total;
  }

  FreedMat(Stprof);

  return SbstMat;
}
/*===========================================================================*/
