
/*=============================================================================
 ntcode.c                                     A.Lambert le 11/01/05

 fonctions destinees a creer des tableaux codant les nucleotides en vue d'opti-
 miser le calcul des scores.

 cc -O2 -Wall -c ntcode.c -I ../include ;

 ar -rs ../lib/librnaIV.a ntcode.o ;

=============================================================================*/

#include "rnaIV.h"


void  SetNtStCode(void);
void  SetNtHlxCode(void);
void  SetNtCodes(void);
void  FreeNtCodes(void);

short *NtStCode = NULL;                            /* codage des nucleotides */
short **NtHlxCode = NULL;                                 /* voir 'ntcode.c' */


/*=============================================================================
 SetNtStCode(): Procede au recodage numerique des nucleotides, aux caracteres
                suivants: A,(T,U),G,C,N et E sont affectees, respectivement,
	        les valeurs 0,1,2,3,5,6 qui indiquent l'ordre dans lequel ils
		seront ranges (a quelle ligne de profil ils sont associes).
                Le code ASCII est utilise pour indexer les elements a coder.
		La ligne relative aux '-' est destinee aux scores du training
		set.
=============================================================================*/

void SetNtStCode(void)
{
  extern short *NtStCode;
  int   i, n = 128;

  if (NtStCode != NULL) return;                /* le tableau est deja alloue */

  NtStCode = (short *) malloc(n * sizeof(short));

  for (i = 0; i < n; i++) NtStCode[i] = 7; /* indice d'1 ligne remplie de 0 */

  NtStCode['A'] = NtStCode['a'] = _A_;                   /* soit la valeur 0 */
  NtStCode['T'] = NtStCode['t'] = _T_;                   /* soit la valeur 1 */
  NtStCode['U'] = NtStCode['u'] = _T_;                   /* soit la valeur 1 */
  NtStCode['G'] = NtStCode['g'] = _G_;                   /* soit la valeur 2 */
  NtStCode['C'] = NtStCode['c'] = _C_;                   /* soit la valeur 3 */
  NtStCode['X'] = NtStCode['-'] = _X_;                   /* soit la valeur 4 */
  NtStCode['N'] = NtStCode['n'] = _N_;                   /* soit la valeur 5 */
  NtStCode['E'] = 6;               /* indice d'une ligne remplie de LOG_ZERO */
}
/*=============================================================================
 SetNtHlxCode(): Cree un tableau carre qui met en correspondance ses elememts
                 avec les paires appartenant aux 2 brins d'une helice. L'ele-
		 ment i,j de ce tableau determine une ligne precise d'un profil
                 d'helice.
=============================================================================*/

void SetNtHlxCode(void)
{
  extern short **NtHlxCode;
  int   n = 128;

  if (NtHlxCode != NULL) return;               /* le tableau est deja alloue */

  NtHlxCode = sMat(n, n);
  FillsMat(NtHlxCode, n, n, 25);  /* indice d'1 ligne du profil remplie de 0 */

  NtHlxCode['A']['A'] = ALPHA_LENxA + _A_;               /* soit la valeur 0 */
  NtHlxCode['A']['T'] = ALPHA_LENxA + _T_;               /* soit la valeur 1 */
  NtHlxCode['A']['G'] = ALPHA_LENxA + _G_;
  NtHlxCode['A']['C'] = ALPHA_LENxA + _C_;

  NtHlxCode['T']['A'] = ALPHA_LENxT + _A_;               /* soit la valeur 4 */
  NtHlxCode['T']['T'] = ALPHA_LENxT + _T_;
  NtHlxCode['T']['G'] = ALPHA_LENxT + _G_;
  NtHlxCode['T']['C'] = ALPHA_LENxT + _C_;

  NtHlxCode['G']['A'] = ALPHA_LENxG + _A_;               /* soit la valeur 8 */
  NtHlxCode['G']['T'] = ALPHA_LENxG + _T_;
  NtHlxCode['G']['G'] = ALPHA_LENxG + _G_;
  NtHlxCode['G']['C'] = ALPHA_LENxG + _C_;

  NtHlxCode['C']['A'] = ALPHA_LENxC + _A_;              /* soit la valeur 12 */
  NtHlxCode['C']['T'] = ALPHA_LENxC + _T_;
  NtHlxCode['C']['G'] = ALPHA_LENxC + _G_;
  NtHlxCode['C']['C'] = ALPHA_LENxC + _C_;


  NtHlxCode['N']['A'] = 16 + _A_;
  NtHlxCode['N']['T'] = 16 + _T_;
  NtHlxCode['N']['G'] = 16 + _G_;
  NtHlxCode['N']['C'] = 16 + _C_;

  NtHlxCode['A']['N'] = 20 + _A_;
  NtHlxCode['T']['N'] = 20 + _T_;
  NtHlxCode['G']['N'] = 20 + _G_;
  NtHlxCode['C']['N'] = 20 + _C_;


  NtHlxCode['E']['A'] = 24;
  NtHlxCode['E']['T'] = 24;
  NtHlxCode['E']['G'] = 24;
  NtHlxCode['E']['C'] = 24;

  NtHlxCode['A']['E'] = 24;
  NtHlxCode['T']['E'] = 24;
  NtHlxCode['G']['E'] = 24;
  NtHlxCode['C']['E'] = 24;

  NtHlxCode['E']['E'] = 24; /*indice d'1 ligne du profil remplie de LOG_ZERO */
  NtHlxCode['N']['N'] = 25;       /* indice d'1 ligne du profil remplie de 0 */
}
/*=============================================================================
 SetNtCodes(): Cree les 2 tableaux de codage.
=============================================================================*/

void SetNtCodes(void)
{
  SetNtStCode();
  SetNtHlxCode();
}
/*=============================================================================
 FreeNtCodes(): Libere les pointeurs des tableaux de codage des nucleotides et
                les reinitialise a la valeur NULL.
=============================================================================*/
void FreeNtCodes(void)
{
  extern short *NtStCode;
  extern short **NtHlxCode;

  if (NtStCode != NULL) {
      free(NtStCode);
      NtStCode = NULL;
  }

  if (NtHlxCode != NULL) {
      FreesMat(NtHlxCode);
      NtHlxCode = NULL;
  }
}
