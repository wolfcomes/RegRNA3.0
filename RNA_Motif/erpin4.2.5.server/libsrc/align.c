
/*=============================================================================
 align.c                           A.Lambert le 20/09/00 revu le 22/11/01

 Ce code regroupe les fonctions procedant a l'alignement d'une sequence sur un
 profil statistique calcule depuis un alignement multiple.
 - Le profil est obtenu en denombrant pour chaque colonne les lettres A,T,G,C
 et les Gaps, les bases inconnues notees 'N' sont incluses dans le decompte des
 autres lettres ou elles sont reparties avec un poids egal (1/ALPHA_LEN).
 - La sequence a aligner sur le profil a une longueur inferieure ou egale a
 celle du profil, et, lors de l'alignement aucun gap ne sera introduit dans le
 profil.
 - Les lettres 'E' de la sequence a aligner donnent une contribution egale a
 LOG_ZERO: les lettres 'E' placees en fin de sequence en signalent la fin.

 Si L est la longueur du profil, |A| celle de l'alphabet utilise, et Pi(X) la
 probabilite (suivant la distribution deduite du profil) du caractere K en po-
 sition 'i', et Po = 1/|A| on procede comme suit:

 Ex: L = 8, seq = "ATGAC"
 P = P1(A).P2(T).P3(-).P4(-).P5(G).P6(-).P7(A).P8(C) est la probabilite de la
 sequence alignee "AT--G-AC".
 L'alignement consiste a maximiser 'P' en deplacant les '-'.
 Plus precisement on maximise l'expression:

 log(P/Po^L) = log(P) + L.log(|A|)        (x^y represente "x la puissance y")

 L'algorithme est une variante de S.B.Needleman et C.D.Wunsch.

 cc -O2 -Wall -c align.c -I../include ;

 ar -rs ../lib/librnaIV.a align.o ;

=============================================================================*/

#include "rnaIV.h"

extern double LOG_ZERO;

void AlignSProfile(char *seq, Strand *St);
void AlignBack(char *seq, Strand *St);

/*=============================================================================
 AlignSProfile(): Calcule les elements de la matrice d'alignement, cette fonc-
                  tion prend ses arguments dans une structure 'Strand' dont les
	          champs 'len', 'min_len', 'max_len' et 'Profile' sont supposes
	          initialises. (notamment par 'ReadStrand').
                  Les elements du tableau sont stockes dans 'St->Align'.
                  L'alignement s'opere entre les 'St->len' premiers caracteres
         	  de 'seq' et les 'St->max_len' caracteres du profil represente
	          par 'St->Profile'.
                  'St->len' est initialise par la fonction appelante.

                  - On gere les caracteres inconnus 'N' de 'seq' en leur asso-
		  ciant une distribution equiprobable 1/ALPHA_LEN = 0.25.
                  Cela revient a interpreter 'seq' comme un profil (trivial).
=============================================================================*/

void AlignSProfile(char *seq, Strand *St)
{
  int     i, j, jm1, jmax;
  double  tmp, max, u;
  extern  short *NtStCode;                 /* declare dans 'libsrc/ntcode.c' */

  St->Align[0][0] = 0.0;
                                          /* initialisation de la 1ere ligne */
  for (j = 1; j <= St->max_gaps; j++)
      St->Align[0][j] = St->Align[0][j-1] + St->Profile[_X_][j-1];

  for (j = 1; j <= St->len; j++)                   /* calcul de la diagonale */
  {
      jm1 = j - 1;
      u = St->Profile[ NtStCode[(int) seq[jm1]] ][jm1];
      St->Align[j][j] = St->Align[jm1][jm1] + u;
  }

  for (i = 1; i <= St->len; i++)               /* fin du calcul des elements */
  {
      jmax = (i + St->max_gaps < St->max_len ? i + St->max_gaps : St->max_len);

      for (j = i + 1; j <= jmax; j++)
      {
          jm1 = j - 1;
	  u = St->Profile[ NtStCode[(int) seq[i-1]] ][jm1];
          max = St->Align[i-1][jm1] + u;
          tmp = St->Align[i][jm1] + St->Profile[_X_][jm1];

          if (tmp > max) max = tmp;

          St->Align[i][j] = max;
      }
  }
  return;
}

/*=============================================================================
 AlignBack(): extraction de la sequence alignee (en fait d'une solution !),
	ou sont positionnes les gaps, qui, en retour, est pointee par 'St.str'.
=============================================================================*/

void AlignBack(char *seq, Strand *St)
{
  int  i = St->len,
       j = St->max_len, jm1;

  St->score = St->Align[i][j];       /* enregistrement du score d'alignement */

  while (j > i && i > 0)
  {
      jm1 = j - 1;
      if (St->Align[i][jm1] > St->Align[i-1][jm1] &&
      fabs(St->Align[i][j] - St->Align[i][jm1] - St->Profile[_X_][jm1]) < 1.e-5)
          St->str[--j] = '-';
      else
          St->str[--j] = seq[--i];
  }

  if (j == i)  while (i-- > 0)  St->str[i] = seq[i];
  else
  if (i == 0)  while (j-- > 0)  St->str[j] = '-';

  St->str[St->max_len] = '\0';

  return;
}
/*=============================================================================
=============================================================================*/
