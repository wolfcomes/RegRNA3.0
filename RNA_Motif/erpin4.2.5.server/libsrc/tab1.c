
/*=============================================================================
 tab1.c                         A.Lambert le 17/09/01   revu le 26/09/01

 Fonctions concernant la gestion de tableaux unidimensionnels d'entiers, impli-
 ques pour la plupart dans le codage de structures secondaires d'ARN.

 La valeur 0 est reservee a la gestion des tableaux, par exemple la detection
 de la fin des donnees.
 Cette contrainte ne s'applique pas a 'FilldTab' et 'FilliTab'.

 cc -O2 -Wall -c tab1.c -I../include ;

 ar -r ../lib/librnaIV.a tab1.o ;

=============================================================================*/

#include "rnaIV.h"

int   DataLen(char *str, int len);
int   TabLen(const int *t);
int   *TabSearch(const int *t, int val);
int   *TabSearchn(int *t, int val, int p);
int   TabRepeats(int *t, int n);

void  FilldTab(double *tab, int len, double val);
void  FilliTab(int *tab, int len, int val);
void  FillStr(char *str, int len, char val);

double *padd(double val, int n);

/*=============================================================================
 DataLen(): Retourne la longueur d'une sequence 'str' de longueur 'len' privee
            des caracteres 'espace' ('\n', ' ' ..) situes dans sa partie finale.
=============================================================================*/

int DataLen(char *str, int len)
{
  char *s = str + (len - 1);     /* s pointe le dernier caractere avant '\0' */
                                                    /* si len == strlen(str) */
  while (isspace(*s)) {                       /* elimine les '\n', espace .. */
      len--;                                      /* situes en fin de lignes */
      s--;
  }
  return len;
}
/*=============================================================================
 TabLen(): Mesure le nombre d'elements d'un tableau d'entiers, dans lequel la
           valeur '0' est reservee pour detecter le dernier element, comme '\0'
           dans le cas d'une chaine de caracteres.
=============================================================================*/

int TabLen(const int *t)
{
  const int *n;

  for (n = t; *n != 0; ++n)   /* nothing */ ;

  return (n - t);
}
/*=============================================================================
 TabSearch(): Cherche la 1ere occurence de la valeur 'val' dans le tableau t[].
              Retourne le pointeur sur la position trouvee, sinon NULL en cas
              d'echec.
=============================================================================*/

int *TabSearch(const int *t, int val)
{
  const int m = val;

  for (; *t != m; ++t)   if (*t == 0) return NULL;

  return ((int *) t);
}
/*=============================================================================
 TabSearchn(): Variante de 'TabSearch' qui cherche l'occurence Nb 'p' (p >= 1)
               de la valeur 'val' dans le tableau t[].
               Retourne le pointeur sur la position trouvee, sinon NULL en cas
               d'echec, par exemple si p <= 0.
=============================================================================*/

int *TabSearchn(int *t, int val, int p)
{
  int i, j, *pi1, *pi2;

  if (p <= 0) return NULL;

  for (i = 0, pi1 = t, j = 0; i < p; i++)
  {
      if ((pi2 = TabSearch(pi1, val)) == NULL) return NULL;
      j += (pi2 - pi1) + 1;
      pi1 = pi2 + 1;
  }
  return (t + j - 1);
}
/*=============================================================================
 TabRepeats(): Compte dans le tableau d'entiers pointe par 't' le nombre de 
               valeurs contigues egales a 'n' et retourne ce nombre.
               Cette fonction est destinee a mesurer la longueur d'un brin dans
               le tableau d'entiers codant la structure secondaire d'un ARN.

               Si 't' est l'adresse du 1er entier de 'nnnnn...' on recuperera
               l'adresse qui le suit immediatement par:
               int *pi;   pi += StrRepeats(t, n);
=============================================================================*/

int TabRepeats(int *t, int n)
{
  int *bgn, *pi;

  if ((pi = TabSearch(t, n)) == NULL)     /* recherche la 1ere occur. de 'n' */
      return 0;
  bgn = pi;

  while (*pi == n) pi++;

  return (pi - bgn);
}
/*=============================================================================
    	ICI LA CONTRAINTE: "0 EST UNE VALEUR RESERVEE" N'EST PAS IMPOSEE.
=============================================================================*/

/*=============================================================================
 FilldTab(): Initialise les 'len' premiers elements du tableau de double pointe
             par 'tab' a la valeur 'val'. La longueur du tableau est supposee
             suffisante.
=============================================================================*/

void  FilldTab(double *tab, int len, double val)
{
  int    i;
  double *ptd;

  for (i = 0, ptd = tab; i < len; i++, ptd++) *ptd = val;

  return;
}
/*=============================================================================
 FilliTab(): Initialise les 'len' premiers elements du tableau d'entiers pointe
             par 'tab' a la valeur 'val'. La longueur du tableau est supposee
             suffisante.
=============================================================================*/

void FilliTab(int *tab, int len, int val)
{
  int    i, *pti;

  for (i = 0, pti = tab; i < len; i++, pti++) *pti = val;

  return;
}
/*=============================================================================
 FillStr(): Initialise les 'len' premiers caracteres du tableau pointe par
            'str' a la valeur 'val'. La longueur du tableau est supposee suf-
            fisante.
=============================================================================*/

void FillStr(char *str, int len, char val)
{
  int    i;
  char   *ptc;

  for (i = 0, ptc = str; i < len; i++, ptc++) *ptc = val;

  return;
}
/*=============================================================================
 padd(): retourne un pointeur sur un tableau de type double de longueur 'n'
         dont les elements sont mis a la valeur 'val'.
=============================================================================*/

double *padd(double val, int n)
{
  double *m, *p;
  int    i;

  if ((m = (double *) malloc(n * sizeof(double))) == NULL) {
      fprintf(stderr, "padd: allocation failure, exit..\n");
      exit(1);
  }
  for (i = 0, p = m; i < n; i++, p++) *p = val;

  return m;
}
/*===========================================================================*/
