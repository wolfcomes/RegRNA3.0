
/*=============================================================================
 hlxsum.c                           A.Lambert le 19/03/04 revu le 31/12/04

 Ce code regroupe les fonctions procedant a la construction d'une matrice de
 substitution pour les paires de bases impliquees dans les helices de la struc-
 ture des ARN regroupees dans un alignement multiple, lequel est utilise comme
 base d'entrainement.

 cc -O2 -Wall -c hlxsum.c -I../include ;

=============================================================================*/

#include "rnaIV.h"

double **GetCorrels3(char **data, int nbstr,
                     int bgn1, int bgn2, int len, double *sweights);
void   GetHlxProfile3(Trset *trset, Helix *Hlx, double *sweights);
void   GetMaskHlxProfiles3(Trset *trset, Mask *mask, double *sweights);
void   FreeMaskHlxProfiles(Mask *mask);
double **HlxSuMat(Trset *trset, Mask *mask, double *sweights, int Vmode);

void PrintdMat(double **m, int nrow, int ncol, FILE *txt);

/*=============================================================================
 GetCorrels3(): calcule le tableau des correlations entre 2 brins d'egale lon-
               gueur 'len' debutant respectivement a 'bgn1' et 'bgn2', struc-
	       tures comme les 2 parties d'une helice, sur un alignement de
	       sequences pointe par 'data'.
               Si la longueur de l'helice est 'len' le tableau est un rectangle
               de 'len' colonnes et SQR_ALPHA_LEN (soit 16) lignes.
	       Les couples contenant un caracteres autre que 'ATGC' (N, gap..)
               sont ignores.
=============================================================================*/

double **GetCorrels3(char **data, int nbstr,
                     int bgn1, int bgn2, int len, double *sweights)
{
  int     j, n, i1 = 0, i2 = 0, findN, height, end2;
  char    *h1, *h2;
  double  **hprofile, weight;

  end2 = bgn2 + len - 1;
  height = SQR_ALPHA_LEN;

  hprofile = dMat(height, len);                        /* creation du profil */
  FilldMat(hprofile, height, len, 0.0);

  for (n = 0; n < nbstr; n++)  /* boucle sur les seq. de la base des donnees */
  {
      h1 = data[n] + bgn1;           /* adresse du 1er caractere du 1er brin */
      h2 = data[n] + end2;      /* adresse du dernier caractere du 2eme brin */
      weight = sweights[n];               /* poids de la sequence consideree */

      for (j = 0; j < len; j++)
      {
          findN = 0;
          switch(h1[j]) {
              case 'A': i1 = ALPHA_LENxA; break;
              case 'T': i1 = ALPHA_LENxT; break;
              case 'G': i1 = ALPHA_LENxG; break;
              case 'C': i1 = ALPHA_LENxC; break;
              default : findN += 1;       break;        /* caractere inconnu */
          }
	  if (findN == 0)
              switch(h2[-j]) {
                  case 'A': i2 = _A_;   break;
                  case 'T': i2 = _T_;   break;
                  case 'G': i2 = _G_;   break;
                  case 'C': i2 = _C_;   break;
                  default : findN += 1; break;          /* caractere inconnu */
              }
          if (findN == 0)   hprofile[i1 + i2][j] += weight;
      }
  }
  return hprofile;
}
/*=============================================================================
 GetHlxProfile3(): Interface de 'GetCorrels3' dont les arguments pointent les
                   structures 'trset', 'Helix' et le tableau de ponderation des
	           sequences.
                   Le profil n'est cree que si Hlx->Profile == NULL.
	           'sweights' pointe le tableau de ponderation des sequences.
=============================================================================*/

void GetHlxProfile3(Trset *trset, Helix *Hlx, double *sweights)
{
  if (Hlx->Profile == NULL)
      Hlx->Profile = GetCorrels3(trset->data, trset->nseq,
                         Hlx->db_bgn1, Hlx->db_bgn2, Hlx->helix_len, sweights);
  else {
    fprintf(stderr, "GetHlxProfile3: allocation failure, exit..\n");
    exit(1);
  }
  return;
}
/*=============================================================================
 GetMaskHlxProfiles3(): Cree l'ensemble des profils des helices d'un masque de
                        pattern (limite au calcul du decompte des couples).
                        Afin d'eviter les allocations multiples, 'GetHlxProfile3'
      		        ne cree le profil que si le pointeur 'Hlx->Profile' a
	                la valeur NULL.
=============================================================================*/

void GetMaskHlxProfiles3(Trset *trset, Mask *mask, double *sweights)
{
  int i, *m;

  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
  {
      GetHlxProfile3(trset, mask->pattern->hlx + *m, sweights);
  }
  return;
}
/*=============================================================================
 FreeMaskHlxProfiles(): Libere la memoire precedemment allouee aux profils
                        d'helices du masque pointe par 'mask'.
			Les pointeurs concernes sont ensuite remis a 'NULL'.
=============================================================================*/

void FreeMaskHlxProfiles(Mask *mask)
{
  int    i, *m;
  Helix  *Hlx;

  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
  {
      Hlx = mask->pattern->hlx + *m;
      if (Hlx->Profile != NULL) {
          FreedMat(Hlx->Profile);
	  Hlx->Profile = NULL;
      }
  }
  return;
}
/*=============================================================================
 HlxSuMat(): Calcule la matrice de substitution des helices du masque pointe
             par 'mask', et retourne un pointeur sur ce tableau carre de cote
	     SQR_ALPHA_LEN.
	     'sweights' pointe le tableau de ponderation des sequences.
	     Si 'Vmode' vaut 1 (ou ON) le profil des brins est affiche.
	     La matrice produite est normalisee de sorte que la somme des ele-
	     ments de chaque colonne soit egale a 1.
=============================================================================*/

double **HlxSuMat(Trset *trset, Mask *mask, double *sweights, int Vmode)
{
  double **SSProf, **SbstMat, total;
  int    i, j, k, width;

  GetMaskHlxProfiles3(trset, mask, sweights);
  SSProf = CatHlxProfiles(mask, &width);

  if (Vmode == ON)  PrintdMat(SSProf, SQR_ALPHA_LEN, width, stdout);

  printf("Hlxs: %d x %d\n", trset->nseq, width);

  FreeMaskHlxProfiles(mask);

     /* --- Creation et calcul de la matrice de substitution des helices --- */

  SbstMat = dMat(SQR_ALPHA_LEN, SQR_ALPHA_LEN);
  FilldMat(SbstMat, SQR_ALPHA_LEN, SQR_ALPHA_LEN, 0.0);

  for (i = 0; i < SQR_ALPHA_LEN; i++)
      for (j = 0; j < SQR_ALPHA_LEN; j++)
          for (k = 0; k < width; k++)
              SbstMat[i][j] += SSProf[i][k] * SSProf[j][k];

  for (j = 0; j < SQR_ALPHA_LEN; j++)
  {
      total = 0.0;                 /* somme des elts de la col. j de SbstMat */
      for (i = 0; i < SQR_ALPHA_LEN; i++) total += SbstMat[i][j];
                                        /* et normalisation a 1 de la col. j */
      for (i = 0; i < SQR_ALPHA_LEN; i++) SbstMat[i][j] /= total;
  }

  FreedMat(SSProf);

  return SbstMat;
}
/*===========================================================================*/
