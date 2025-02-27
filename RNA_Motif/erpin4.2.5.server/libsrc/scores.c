
/*=============================================================================
 scores.c                         A.Lambert le 20/11/01 revu le 12/02/04

 Ce code regroupe les fonctions relatives au calcul des scores obtenus en pre-
 sentant une sequence devant les divers profils calcules a partir d'une base de
 donnees.
 Les profils concernent des brins et des helices.

 Pour l'evaluation des scores 'GetStNoGScore' et 'GetHlxScores' sont les 2
 fonctions principales, auxquelles il faut ajouter 'AlignSProfile' de 'align.c'
 utilisee si des gaps etaient presents lors de la construction du profil.
 Les autres fonctions en sont derivees.

 cc -O2 -Wall -c scores.c -I ../include ;

 ar -r ../lib/librnaIV.a scores.o ;

=============================================================================*/

#include "rnaIV.h"

extern double LOG_ZERO;

double GetStNoGScore(char *seq, Strand *St);
void   GetHlxScores(char *seq, Helix *Hlx, double *scores);
double GetHlxBestScore(char *seq, Helix *Hlx);

void   GetStScores(char *seq, Strand *St, double *scores);
double GetStBestScore(char *seq, Strand *St);

void   GetHlxScoresTab(char *seq, int length, Helix *Hlx);
void   GetStScoresTab(char *seq, int length, Strand *St);

void   GetHlxScoresTabBis(char *seq, int length, Helix *Hlx);
void   GetStScoresTabBis(char *seq, int length, Strand *St);

/*=============================================================================
 GetStNoGScore(): Pour un brin pointe par 'St' ne contenant pas de gaps, cette
	          fonction calcule et retourne le score obtenu par 'seq' sur
                  une longueur egale au nombre de colonnes du profil pointe par
        	  'St'.
	          - On gere les caracteres inconnus 'N' de 'sequence' en leur
	          associant une distribution equiprobable 1 / ALPHA_LEN = 0.25.
                  - Les 'E' explicitement absents des profils contribuent pour
		  LOG_ZERO (leur presence marque la fin d'une sequence).
=============================================================================*/

double GetStNoGScore(char *seq, Strand *St)
{
  int    i, j;
  double score;
  extern short *NtStCode;                  /* declare dans 'libsrc/ntcode.c' */
                                                 /* bouclee sur les colonnes */
  for (j = 0, score = 0.0; j < St->max_len; j++)
  {
      i = NtStCode[(int) seq[j]];
      score += St->Profile[i][j];
  }
  return score;
}
/*=============================================================================
 GetHlxScores(): calcule, depuis la position pointee par 'seq', les scores ob-
	         tenus en presentant devant les caracteres de la chaine le pro-
		 fil de l'helice pointee par 'Hlx'.

                 Le tableau 'scores' (indexe de 0 au nombre max. de gaps) re-
		 coit les valeurs obtenues lorsque la distance entre les deux
		 parties supposees de l'helice varie entre 'Hlx.min_dist' et
	        'Hlx.max_dist', suivant les indications du modele et de la base
		 de donnees utilisee.
=============================================================================*/

void GetHlxScores(char *seq, Helix *Hlx, double *scores)
{
  int     i, j, k;
  char    *h, *right;
  extern  short **NtHlxCode;               /* declare dans 'libsrc/ntcode.c' */

  h = seq + Hlx->min_len - 1;

  for (k = 0; k <= Hlx->max_gaps; k++)
  {
      scores[k] = 0.0;
      right = h + k;

      for (j = 0; j < Hlx->helix_len; j++, right--)
      {
          i = NtHlxCode[(int) seq[j]][(int) *right];
          scores[k] += Hlx->Profile[i][j];
      }
  }
}
/*=============================================================================
 GetHlxBestScore(): Calcule, a la position pointee par 'seq', le meilleur score
	            obtenu en presentant devant les caracteres de la chaine le
	            profil d'une helice pointee par 'Hlx'.
=============================================================================*/

double GetHlxBestScore(char *seq, Helix *Hlx)
{
  int     i, j, k, d;
  char    *h, *right;
  double  score, score_max;
  extern  short **NtHlxCode;               /* declare dans 'libsrc/ntcode.c' */

  score_max = 2.0 * LOG_ZERO * Hlx->helix_len; /* sera tout de suite modifie */
  Hlx->dist = Hlx->max_dist;

  h = seq + Hlx->min_len - 1;

  for (k = d = 0; k <= Hlx->max_gaps; k++)
  {
      score = 0.0;
      right = h + k;

      for (j = 0; j < Hlx->helix_len; j++, right--)
      {
          i = NtHlxCode[(int) seq[j]][(int) *right];
          score += Hlx->Profile[i][j];
      }
      if (score > score_max) {
          score_max = score;
          d = k;                          /* indice associe au score maximal */
      }
  }
  Hlx->score = score_max;                        /* extraction des resultats */
  Hlx->dist  = Hlx->min_dist + d;

  return score_max;
}
/*=============================================================================
 GetStScores(): calcule, a la position pointee par 'seq', les scores
	obtenus en presentant devant les caracteres de la chaine le profil d'un
	brin (structure Strand)..

        Le tableau 'scores' (indexe de 0 au nombre max. de gaps) recoit les
        valeurs obtenues lorsque la distance entre les deux parties supposees
        de l'helice varie entre 'St.min_len' et 'St.max_len', suivant les
        indications du modele et de la base de donnees utilisee (gaps).
        Si St.min_len == St.max_len (pas de gaps) un seul nombre est retoune,
        calcule par 'GetStTransitsScore'.
=============================================================================*/

void GetStScores(char *seq, Strand *St, double *scores)
{
  int  i;

  if (St->max_gaps == 0)
      scores[0] = GetStNoGScore(seq, St);
  else
  {
      St->len = St->max_len;
      AlignSProfile(seq, St);

      for (i = 0; i <= St->max_gaps; i++)
          scores[i] = St->Align[St->min_len + i][St->max_len];
  }
  return;
}
/*=============================================================================
 GetStBestScore(): Calcule, a la position pointee par 'seq', le meilleur score
        obtenu en presentant devant les caracteres de la chaine le profil d'un
        brin (structure Strand).
=============================================================================*/

double GetStBestScore(char *seq, Strand *St)
{
  int    i, len;
  double score, tmp;


  if (St->max_gaps == 0)
  {
      score = GetStNoGScore(seq, St);
      len = St->max_len;
  }
  else
  {
      St->len = St->max_len;
      AlignSProfile(seq, St);

      score = St->Align[St->min_len][St->max_len];
      len = St->min_len;

      for (i = 1; i <= St->max_gaps; i++) {
          if ((tmp = St->Align[St->min_len + i][St->max_len]) > score)
          {
              score = tmp;
              len = St->min_len + i;
          }
      }
  }
  St->len = len;
  St->score = score;

  return score;
}
/*=============================================================================
 GetHlxScoresTab(): Reprend le corps de la fonction 'GetHlxScores' et balaye
                    une sequence sur une longueur 'length' afin de dresser le
	            tableau bidimensionnel 'Hlx.Scores' ('Hlx.max_gaps + 1' li-
	            gnes et 'length' colonnes) des scores d'helices.
                    L'espace memoire necessaire est suppose avoir ete alloue a
	            'Hlx.Scores'.
=============================================================================*/

void GetHlxScoresTab(char *seq, int length, Helix *Hlx)
{
  int     s, i, j, k;
  char    *h, *left, *right;
  extern  short **NtHlxCode;               /* declare dans 'libsrc/ntcode.c' */

  if (Hlx->Scores == NULL) {
      fprintf(stderr, "GetHlxScoresTab: array not allocated, exit..\n");
      exit(1);
  }

  for (s = 0; s < length; s++)
  {
      left = seq + s;          /* pointe le car. le + a gauche de l'ensemble */
      h = left + Hlx->min_len - 1;

      for (k = 0; k <= Hlx->max_gaps; k++)
      {
          double u = 0.0;
          right = h + k;

          for (j = 0; j < Hlx->helix_len; j++, right--)
          {
              i = NtHlxCode[(int) left[j]][(int) *right];
              u += Hlx->Profile[i][j];
          }
          Hlx->Scores[k][s] = (float) u;
      }
  }
}
/*=============================================================================
 GetHlxScoresTabBis(): Variante de 'GetHlxScoresTab' qui enregistre les scores
                       dans le 2eme tableau 'ScoresBis'.
=============================================================================*/

void GetHlxScoresTabBis(char *seq, int length, Helix *Hlx)
{
  int     s, i, j, k;
  char    *h, *left, *right;
  extern  short **NtHlxCode;               /* declare dans 'libsrc/ntcode.c' */

  if (Hlx->ScoresBis == NULL) {
      fprintf(stderr, "GetHlxScoresTabBis: array not allocated, exit..\n");
      exit(1);
  }

  for (s = 0; s < length; s++)
  {
      left = seq + s;          /* pointe le car. le + a gauche de l'ensemble */
      h = left + Hlx->min_len - 1;

      for (k = 0; k <= Hlx->max_gaps; k++)
      {
          double u = 0.0;
          right = h + k;

          for (j = 0; j < Hlx->helix_len; j++, right--)
          {
              i = NtHlxCode[(int) left[j]][(int) *right];
              u += Hlx->Profile[i][j];
          }
          Hlx->ScoresBis[k][s] = (float) u;
      }
  }
}
/*=============================================================================
 GetStScoresTab(): reprend le corps de la fonction 'GetStScores' et balaye
        une sequence sur une longueur 'length' pour dresser le tableau bidimen-
        sionnel 'St->Scores' (St->max_len - St->min_len + 1 lignes, soit:
        St.max_gaps + 1, et 'length' colonnes) des scores de brins.
        L'espace memoire necessaire est suppose avoir ete alloue a 'St->Scores'.
=============================================================================*/

void GetStScoresTab(char *seq, int length, Strand *St)
{
  int  s, i;

  if (St->Scores == NULL) {
      fprintf(stderr, "GetStScoresTab: array not allocated, exit..\n");
      exit(1);
  }
  if (St->max_gaps == 0)
      for (s = 0; s < length; s++)
      {
          St->Scores[0][s] = (float) GetStNoGScore(seq + s, St);
      }
  else
      for (s = 0, St->len = St->max_len; s < length; s++)
      {
          AlignSProfile(seq + s, St);

          for (i = 0; i <= St->max_gaps; i++)
          St->Scores[i][s] = (float) St->Align[St->min_len + i][St->max_len];
      }
  return;
}
/*=============================================================================
 GetStScoresTabBis(): Variante de 'GetStScoresTab' qui enregistre les scores
                      dans le 2eme tableau 'ScoresBis'.
=============================================================================*/

void GetStScoresTabBis(char *seq, int length, Strand *St)
{
  int  s, i;

  if (St->ScoresBis == NULL) {
      fprintf(stderr, "GetStScoresTabBis: array not allocated, exit..\n");
      exit(1);
  }
  if (St->max_gaps == 0)
      for (s = 0; s < length; s++)
      {
          St->ScoresBis[0][s] = (float) GetStNoGScore(seq + s, St);
      }
  else
      for (s = 0, St->len = St->max_len; s < length; s++)
      {
          AlignSProfile(seq + s, St);

          for (i = 0; i <= St->max_gaps; i++)
          St->ScoresBis[i][s] = (float) St->Align[St->min_len + i][St->max_len];
      }
  return;
}
/*===========================================================================*/
