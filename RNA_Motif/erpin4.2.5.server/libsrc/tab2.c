
/*=============================================================================
 tab2.c                     A.Lambert le 15/04/2000 revu le 06/10/2000

 quelques fonctions de gestion de tableaux bidimensionnels.

 cc -O2 -Wall -c tab2.c -I../include ;

 ar -rs ../lib/librnaIV.a tab2.o ;

=============================================================================*/

#include "rnaIV.h"

char   **cMat(int nrow, int ncol);
short  **sMat(int nrow, int ncol);
int    **iMat(int nrow, int ncol);
double **dMat(int nrow, int ncol);
float  **fMat(int nrow, int ncol);
void   FreecMat(char **m);
void   FreesMat(short **m);
void   FreeiMat(int **m);
void   FreedMat(double **m);
void   FreefMat(float **m);
void   FillcMat(char **m, int nrow, int ncol, char val);
void   FillsMat(short **m, int nrow, int ncol, short val);
void   FilliMat(int **m, int nrow, int ncol, int val);
void   FilldMat(double **m, int nrow, int ncol, double val);
void   FillfMat(float **m, int nrow, int ncol, float val);
void   CopysMat(short **src, int nrow, int ncol, short **target);

/*=============================================================================
 cMat(): cree une matrice rectangulaire de caracteres (char). 
=============================================================================*/

char **cMat(int nrow, int ncol)
{
  int  i;
  char **m;

  if ((m = (char **) malloc(nrow * sizeof(char *))) == NULL) {
      fprintf(stderr, "cMat: allocation failure 1, exit..\n");
      exit(1);
  }
  if ((m[0] = (char *) malloc(nrow * ncol * sizeof(char))) == NULL) {
      fprintf(stderr, "cMat: allocation failure 2, exit..\n");
      exit(1);
  }
  for (i = 1; i < nrow; i++) m[i] = m[i-1] + ncol;

  return m;
}
/*=============================================================================
 sMat(): cree une matrice rectangulaire d'entiers courts (2 octets). 
=============================================================================*/

short **sMat(int nrow, int ncol)
{
  int  i;
  short  **m;

  if ((m = (short **) malloc(nrow * sizeof(short *))) == NULL) {
      fprintf(stderr, "sMat: allocation failure 1, exit..\n");
      exit(1);
  }
  if ((m[0] = (short *) malloc(nrow * ncol * sizeof(short))) == NULL) {
      fprintf(stderr, "sMat: allocation failure 2, exit..\n");
      exit(1);
  }
  for (i = 1; i < nrow; i++) m[i] = m[i-1] + ncol;

  return m;
}
/*=============================================================================
 iMat(): cree une matrice rectangulaire d'entiers (4 octets). 
=============================================================================*/

int **iMat(int nrow, int ncol)
{
  int  i;
  int  **m;

  if ((m = (int **) malloc(nrow * sizeof(int *))) == NULL) {
      fprintf(stderr, "iMat: allocation failure 1, exit..\n");
      exit(1);
  }
  if ((m[0] = (int *) malloc(nrow * ncol * sizeof(int))) == NULL) {
      fprintf(stderr, "iMat: allocation failure 2, exit..\n");
      exit(1);
  }
  for (i = 1; i < nrow; i++) m[i] = m[i-1] + ncol;

  return m;
}
/*=============================================================================
 dMat(): cree une matrice rectangulaire de nombres de type double (8 octets). 
=============================================================================*/

double **dMat(int nrow, int ncol)
{
  int    i;
  double **m;

  if ((m = (double **) malloc(nrow * sizeof(double *))) == NULL) {
      fprintf(stderr, "dMat: allocation failure 1, exit..\n");
      exit(1);
  }
  if ((m[0] = (double *) malloc(nrow * ncol * sizeof(double))) == NULL) {
      fprintf(stderr, "dMat: allocation failure 2, exit..\n");
      exit(1);
  }
  for (i = 1; i < nrow; i++) m[i] = m[i-1] + ncol;

  return m;
}
/*=============================================================================
 fMat(): cree une matrice rectangulaire de nombres de type float (4 octets). 
=============================================================================*/

float **fMat(int nrow, int ncol)
{
  int    i;
  float  **m;

  if ((m = (float **) malloc(nrow * sizeof(float *))) == NULL) {
      fprintf(stderr, "fMat: allocation failure 1, exit..\n");
      exit(1);
  }
  if ((m[0] = (float *) malloc(nrow * ncol * sizeof(float))) == NULL) {
      fprintf(stderr, "fMat: allocation failure 2, exit..\n");
      exit(1);
  }
  for (i = 1; i < nrow; i++) m[i] = m[i-1] + ncol;

  return m;
}
/*=============================================================================
 FreecMat(): restitution de la memoire attribuee a un tableau de caracteres. 
=============================================================================*/

void FreecMat(char **m) 
{
  free((char *) m[0]);  
  free((char *) m);
  return;
}
/*=============================================================================
 FreesMat(): restitution de la memoire attribuee a un tableau d'entiers courts. 
=============================================================================*/

void FreesMat(short **m) 
{
  free((char *) m[0]);  
  free((char *) m);
  return;
}
/*=============================================================================
 FreeiMat(): restitution de la memoire attribuee a un tableau d'entiers. 
=============================================================================*/

void FreeiMat(int **m) 
{
  free((char *) m[0]);  
  free((char *) m);
  return;
}
/*=============================================================================
 FreedMat(): restitution de la memoire attribuee a un tableau de double. 
=============================================================================*/

void FreedMat(double **m) 
{
  free((char *) m[0]);  
  free((char *) m);
  return;
}
/*=============================================================================
 FreefMat(): restitution de la memoire attribuee a un tableau de float. 
=============================================================================*/

void FreefMat(float **m) 
{
  free((char *) m[0]);  
  free((char *) m);
  return;
}
/*=============================================================================
 FilldMat(): initialise les elements d'une matrice de doubles a 'val'.
=============================================================================*/

void FilldMat(double **m, int nrow, int ncol, double val)
{
  int     i;
  double  *ptr = m[0];

  for (i = 0; i < nrow * ncol; i++, ptr++)  *ptr = val;

  return;
}
/*=============================================================================
 FillfMat(): initialise les elements d'une matrice de floats a 'val'.
=============================================================================*/

void FillfMat(float **m, int nrow, int ncol, float val)
{
  int     i;
  float  *ptr = m[0];

  for (i = 0; i < nrow * ncol; i++, ptr++)  *ptr = val;

  return;
}
/*=============================================================================
 FilliMat(): initialise les elements d'une matrice d'entiers a 'val'.
=============================================================================*/

void FilliMat(int **m, int nrow, int ncol, int val)
{
  int  i, *ptr = m[0];

  for (i = 0; i < nrow * ncol; i++, ptr++)  *ptr = val;

  return;
}
/*=============================================================================
 FillsMat(): initialise les elements d'une matrice d'entiers courts a 'val'.
=============================================================================*/

void FillsMat(short **m, int nrow, int ncol, short val)
{
  int   i;
  short *ptr = m[0];

  for (i = 0; i < nrow * ncol; i++, ptr++)  *ptr = val;

  return;
}
/*=============================================================================
 FillcMat(): initialise les elements d'une matrice de type char a 'val'.
=============================================================================*/

void FillcMat(char **m, int nrow, int ncol, char val)
{
  int   i;
  char  *ptr = m[0];

  for (i = 0; i < nrow * ncol; i++, ptr++)  *ptr = val;

  return;
}
/*=============================================================================
 CopysMat(): copie le contenu du tableau pointe par 'src' dans celui pointe par
             'target'. Les dimensions 'nrow' et 'ncol' de la copie sont suppo-
             sees compatibles avec celles des 2 tableaux d'entiers courts.
=============================================================================*/

void CopysMat(short **src, int nrow, int ncol, short **target)
{
  int i, j;

  for (i = 0; i < nrow; i++)
      for (j = 0; j < ncol; j++) target[i][j] = src[i][j];

  return;
}
/*===========================================================================*/
