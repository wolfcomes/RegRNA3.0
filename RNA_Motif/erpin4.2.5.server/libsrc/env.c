
/*=============================================================================
 env.c                           A.Lambert le 16/11/01  revu le 12/02/04

 Ce code contient les declarations de quelques variables fixant l'environnement
 des programmes et des variables globales manipulees par certaines fonctions.
 Sont aussi presentes les fonctions qui les initialisent et les modifient.

 cc -O2 -Wall -c env.c ;

 ar -rs ../lib/librnaIV.a env.o ;

=============================================================================*/


#define DEFAULT_LOG_ZERO      -20.0            /* borne inferieure de log(x) */
#define DEFAULT_SCORE_TAB_LEN  1024          /* long. des tableaux de scores */

double  LOG_ZERO = DEFAULT_LOG_ZERO;           /* borne inferieure de log(x) */
                                 /* tableaux de la statistique des sequences */

double  *NewLogDataFreqs, *LogDataFreqs;            /* pointeurs non alloues */
double  *LogDataFreqs1, *LogDataFreqs2;                 /* pointeurs alloues */

int     ScoreTabLen = DEFAULT_SCORE_TAB_LEN; /* long. des tableaux de scores */

unsigned long long int TotalNucScans = 0;  /* nombre total de bases traitees */

/*=============================================================================
 ChScoreTabLen(): modifie la longueur des tableaux de scores precalcules pour
                  la valeur passee en argument.
=============================================================================*/

void ChScoreTabLen(int val)
{
  ScoreTabLen = val;
  return;
}
/*=============================================================================
 ChLogZero(): modifie la valeur de 'LOG_ZERO' pour celle passee en argument.
=============================================================================*/

void ChLogZero(double val)
{
  LOG_ZERO = val;
  return;
}
/*===========================================================================*/
