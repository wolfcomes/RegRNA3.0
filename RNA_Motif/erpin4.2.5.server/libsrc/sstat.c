
/*=============================================================================
 sstat.c                      A.Lambert le 16/11/2001 revu le 13/02/04

 Ce code concerne un ensemble de fonctions destinees a faire la statistique des 
 bases dans des sequences, afin de completer ou actualiser le contenu des
 profils statistiques:

 - denombrer les frequences des bases A,T,G,C,

 et creer des sequences possedant une composition aleatoire de A,T,G,C donnee.

 2 tableaux seront utilises:
 - de longueur ALPHA_LEN pour les bases A,T,G,C: LogDataFreqs1 et ..2

 2 pointeurs accompagnent ces tableaux pour organiser les permutations de leur
 role dans un processus iteratif:
 - enregistrement de donnees saisies,
 - memorisation apres usage des donnees.

 double *NewLogDataFreqs, *LogDataFreqs;

 cc -Wall -O2 -c sstat.c -I../include ;

 ar -rs ../lib/librnaIV.a sstat.o ;

=============================================================================*/

#include "rnaIV.h"

void sGetStat(char *seq, int len);
void InitStatTables(void);
void ReInitStatTables(double *freqs);
void SwapStatTables(void);
void ReverseStat(void);
char *MakeRandSeq(int len);
double *Cumuls(double *freqs);
void SetRandSeq(char *seq, int len, double *cumuls);
char *RandSeq(int len, double *freqs);

void ChStStat(Strand *St);
void ChHlxStat(Helix *Hlx);
void ChPatternStat(Pattern *pattern);

extern double *NewLogDataFreqs, *LogDataFreqs;
extern double *LogDataFreqs1, *LogDataFreqs2;

/*=============================================================================
 sGetStat(): opere sur la sequence pointee par 'seq' le decompte des bases
            ATGC (ignore les 'N') et des transitions sur une longueur 'len'
            (qui devra rester <= a la longueur de la sequence lue).
            'sGetStat' sera appliquee a des sequences 'preparees': converties
            en majuscules, 'U' change en 'T'.. mais supporte les erreurs.

            Le tableau 'NewLogDataFreqs', declare en variables globale, est
            reinitialise en debut d'action, avant le traitement d'une sequence.
            Afin d'eviter les problemes associes aux elements nuls (possibles
            dans certaines sequences courtes) l'initialisation des frequences
            est faite a 0.01.
=============================================================================*/

void sGetStat(char *seq, int len)
{
  int     i;
  char    *ptc;
  double  sum;

  extern double *NewLogDataFreqs;

  FilldTab(NewLogDataFreqs, ALPHA_LEN, 0.01);

  for (i = 0, ptc = seq; i < len; i++, ptc++)
  {
      switch(*ptc) {
          case 'A': NewLogDataFreqs[_A_]++;  break;
          case 'T': NewLogDataFreqs[_T_]++;  break;
          case 'G': NewLogDataFreqs[_G_]++;  break;
          case 'C': NewLogDataFreqs[_C_]++;  break;
          default :                       break;
      }
  }
                                                 /* fin du calcul du tableau */

  for (i = 0, sum = 0.0; i < ALPHA_LEN; i++)  sum += NewLogDataFreqs[i];

  for (i = 0; i < ALPHA_LEN; i++)  NewLogDataFreqs[i] /= sum;  /* normalisation */

  for (i = 0; i < ALPHA_LEN; i++)                  /* passage aux logarithme */
      NewLogDataFreqs[i] = log(NewLogDataFreqs[i]);

  return; 
}
/*=============================================================================
 InitStatTables(): Cree et initialise (a 1/ALPHA_LEN) les tableaux de la sta-
                   tistique sur les sequences a explorer ulterieurement. Ces
                   tableaux permettent d'actualiser les profils statistiques,
                   et les pointeurs a gerer les permutations (saisie/memoire).
=============================================================================*/

void InitStatTables(void)
{
  double logpo = log(1.0 / ALPHA_LEN);

  extern double *NewLogDataFreqs, *LogDataFreqs;
  extern double *LogDataFreqs1, *LogDataFreqs2;

  LogDataFreqs1 = (double *) malloc(ALPHA_LEN * sizeof(double));
  LogDataFreqs2 = (double *) malloc(ALPHA_LEN * sizeof(double));

  FilldTab(LogDataFreqs1, ALPHA_LEN, logpo);
  FilldTab(LogDataFreqs2, ALPHA_LEN, logpo);
                                               /* etat initial des pointeurs */
  NewLogDataFreqs = LogDataFreqs1;

  LogDataFreqs = LogDataFreqs2;

  return;
}
/*=============================================================================
 ReInitStatTables(): Reinitialise les valeurs pointees par 'freqs' les tableaux
                     de la statistique sur les sequences a explorer ulterieure-
		     ment. Ces tableaux permettent d'initialiser les profils
        	     statistiques, et les pointeurs a gerer les permutations
		     (saisie/memoire).
=============================================================================*/

void ReInitStatTables(double *freqs)
{
  int    i;
  double small = exp(LOG_ZERO);
  extern double *NewLogDataFreqs, *LogDataFreqs;
  extern double *LogDataFreqs1, *LogDataFreqs2;

  for (i = 0; i < ALPHA_LEN; i++) {
      if (freqs[i] > small)
          LogDataFreqs1[i] = LogDataFreqs2[i] = log(freqs[i]);
      else
          LogDataFreqs1[i] = LogDataFreqs2[i] = LOG_ZERO;
  }
                                               /* etat initial des pointeurs */
  NewLogDataFreqs = LogDataFreqs1;
  LogDataFreqs = LogDataFreqs2;

  return;
}
/*=============================================================================
 SwapStatTables(): Permute les valeurs des pointeurs sur les tableaux des don-
                   nees statistiques (plutot que de permuter le contenu).
=============================================================================*/

void SwapStatTables(void)
{
  double *tmpFreqs;

  tmpFreqs = NewLogDataFreqs;
  NewLogDataFreqs = LogDataFreqs;
  LogDataFreqs = tmpFreqs;

  return;
}
/*=============================================================================
 ReverseStat(): echange, a partir des frequences et des transitons des ATGC en-
                registrees dans 'LogDataFreqs', ces valeurs pour s'adapter au
                au brin complementaire:
                A->T, T->A, G->C et C->G + renversement de l'ordre.
                Les nouvelles valeurs sont enregistrees dans 'NewLogDataFreqs'
                (comme dans la saisie d'une nouvelle sequence).
=============================================================================*/

void ReverseStat(void)
{
  int i, I;

  for (i = 0; i < ALPHA_LEN; i++)       /* tableau des frequences des bases */
  {
      I = (i%2 == 0 ? i+1 : i-1);        /*  A,T,G,C (= 0,1,2,3) -> T,A,C,G */
      NewLogDataFreqs[i] = LogDataFreqs[I];
  }
  return;
}
/*=============================================================================
 MakeRandSeq(): Genere une sequence aleatoire de longueur 'len', dont la compo-
                sition en bases ATGC est fixee par les elements du tableau
	        'NewLogDataFreqs', qui est declare en variable globale dans
	        'env.c'.
                Retourne un pointeur sur le tableau cree.
=============================================================================*/

char *MakeRandSeq(int len)
{
  const double RANDMAXp1 = RAND_MAX + 1.0;
  int    i;
  char   *seq;
  double randval, A_ratio, AT_ratio, ATG_ratio;

  extern double *NewLogDataFreqs;

  if ((seq = (char *) malloc(len + 1)) == NULL) {
      fprintf(stderr, "MakeRandSeq: allocation failure, exit..\n");
      exit(1);
  }
  A_ratio = exp(NewLogDataFreqs[_A_]);
  AT_ratio = A_ratio + exp(NewLogDataFreqs[_T_]);
  ATG_ratio = AT_ratio + exp(NewLogDataFreqs[_G_]);

  for (i = 0; i < len; i++)
  {
      randval = rand()/RANDMAXp1;
      if (randval < A_ratio)   seq[i] = 'A';
      else
      if (randval < AT_ratio)  seq[i] = 'T';
      else
      if (randval < ATG_ratio) seq[i] = 'G';
      else
          seq[i] = 'C';
  }
  seq[len] = '\0';

  return seq;
}
/*=============================================================================
 Cumuls(): Donne les frequences cumulees des bases A, AT et ATG depuis les don-
           nees du tableau pointe par 'freqs'.
=============================================================================*/

double *Cumuls(double *freqs)
{
  double *cumuls = (double *) malloc(ALPHA_LEN * sizeof(double));

  cumuls[0] = freqs[0];
  cumuls[1] = cumuls[0] + freqs[1];
  cumuls[2] = cumuls[1] + freqs[2];

  return cumuls;
}
/*=============================================================================
 SetRandSeq(): Initialise une sequence aleatoire, prealablement creee, pointee
               par 'seq' de longueur 'len', dont la composition en bases ATGC
               est lue dans les elements du tableau 'cumuls', donnant les fre-
               quences cumulees des A, AT, ATG (pour optimisation).
               'seq' est supposee pointer un volume suffisant (pas de controle)
=============================================================================*/

void SetRandSeq(char *seq, int len, double *cumuls)
{
  int    i;
  double rndval;

  for (i = 0; i < len; i++)
  {
      rndval = rand() / (RAND_MAX + 1.0);   /* nbre aleat. unif. dans (0,1) */
      if (rndval < cumuls[0])  seq[i] = 'A';
      else
      if (rndval < cumuls[1])  seq[i] = 'T';
      else
      if (rndval < cumuls[2])  seq[i] = 'G';
      else
          seq[i] = 'C';
  }
  seq[len] = '\0';

  return;
}
/*=============================================================================
 RandSeq(): Genere une sequence aleatoire de longueur 'len', dont la composi-
            tion en bases ATGC est issue des elements du tableau 'freqs'.
            Retourne un pointeur sur la sequence creee.
=============================================================================*/

char *RandSeq(int len, double *freqs)
{
  int    i;
  char   *seq;
  double randval, A_ratio, AT_ratio, ATG_ratio;


  if ((seq = (char *) malloc(len + 1)) == NULL) {
      fprintf(stderr, "RandSeq: allocation failure, exit..\n");
      exit(1);
  }
  A_ratio = freqs[0];                                       /* freq. des 'A' */
  AT_ratio = A_ratio + freqs[1];                           /* freq. des 'AT' */
  ATG_ratio = AT_ratio + freqs[2];                        /* freq. des 'ATG' */

  for (i = 0; i < len; i++)
  {
      randval = rand() / (RAND_MAX + 1.0);
      if (randval < A_ratio)   seq[i] = 'A';
      else
      if (randval < AT_ratio)  seq[i] = 'T';
      else
      if (randval < ATG_ratio) seq[i] = 'G';
      else
          seq[i] = 'C';
  }
  seq[len] = '\0';

  return seq;
}
/*=============================================================================
 ChStStat(): Change dans le profil d'un brin les frequences des ATGC du "fond".
             'NewLogDataFreqs' remplace 'LogDataFreqs' (supposes actualises).
=============================================================================*/

void ChStStat(Strand *St) 
{
  int    i, j;
  double u;

  extern double *LogDataFreqs, *NewLogDataFreqs;

  for (i = 0; i < ALPHA_LEN; i++)
  {
      u = LogDataFreqs[i] - NewLogDataFreqs[i];
      for (j = 0; j < St->max_len; j++)  St->Profile[i][j] += u;
  }
  return;
}
/*=============================================================================
 ChHlxStat(): Change dans le tableau des correlations d'une Helice les frequ-
              ences "attendues" des ATGC.
             'NewLogDataFreqs' remplace 'LogDataFreqs' (supposes actualises).
=============================================================================*/

void ChHlxStat(Helix *Hlx)
{
  int    i, j, I, J;
  double u;

  extern double *LogDataFreqs, *NewLogDataFreqs;

  for (i = 0; i < SQR_ALPHA_LEN; i++)
  {
      I = i / ALPHA_LEN;
      J = i % ALPHA_LEN; 
      u = LogDataFreqs[I] + LogDataFreqs[J] -
          (NewLogDataFreqs[I] + NewLogDataFreqs[J]);
      for (j = 0; j < Hlx->helix_len; j++)  Hlx->Profile[i][j] += u;
  }
  return;
}
/*=============================================================================
 ChPatternStat(): Interface, pour une structure 'Pattern', des fonctions rela-
                  tives a la modification de statistique pour ses elements
                  constitutifs: brins et helices.
                  L'operation s'acheve par la permutation des tableaux conte-
                  nant les donnees statistiques.
=============================================================================*/

void ChPatternStat(Pattern *pattern)
{
  int i;

  for (i = 0; i < pattern->nst; i++)  ChStStat(pattern->std + i);

  for (i = 0; i < pattern->nhx; i++)  ChHlxStat(pattern->hlx + i);

  SwapStatTables();

  return;
}
/*===========================================================================*/
