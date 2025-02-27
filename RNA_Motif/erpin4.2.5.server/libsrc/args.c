
/*=============================================================================
 args.c                           A.Lambert le 16/11/01

 Ce code contient des fonctions destinees a saisir la plupart des arguments
 entrees aux programmes, notamment 'erpin'.

 cc -O2 -Wall -c args.c -I../include ;

 ar -rs ../lib/librnaIV.a args.o ;

=============================================================================*/

#include "rnaIV.h"

void ReadHelpArgs(int argc, char *argv[], int argcmin, void (*helpfunc)(void));
void GetOutputArgs(int argc, char *argv[], int *style, int *warnings,
                   int *evflag, int *histflag);
void GetDirArg(int argc, char *argv[], int *dir);
void GetStatArg(int argc, char *argv[], int *stat);
void GetSeqArgs(int argc, char *argv[],
                int *seqnb1, int *nseq, int *bgn, int *range);
void GetEnvArgs(int argc, char *argv[], int *chrono);
void GetProcArg(int argc, char *argv[], int *maskproc);
char *GetSumArgs(int argc, char *argv[], int *pcflag, double *hpcw, double *spcw);

/*=============================================================================
 ReadHelpArgs(): Lit dans les arguments d'un programme l'appel a un 'help' et
                 affiche un texte approprie par l'intermediaire de l'argument
                 fonction 'helpfunc'.
                 Controle que 'argc' est superieur au nombre 'argcmin', qui est
                 le nombre minimal d'arguments necessaires..
=============================================================================*/

void ReadHelpArgs(int argc, char *argv[], int argcmin, void (*helpfunc)(void))
{
  int  i, prthelp = 0;
  char *tail;

  for (i = prthelp = 0; i < argc; i++)
  {
      if (argv[i][0] == '-')
          if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0) {
              prthelp = 1;
              break;
          }
  }
  if (argc == 1 || prthelp) {
      helpfunc();
      exit(0);
  }
  if (argc < argcmin + 1) {
      GetPathTail(argv[0], &tail);
      fprintf(stderr, "\n\"%s\" needs at least %d argument%s, exit..\n",
      tail, argcmin, argcmin > 1 ? "s" : "");
      free(tail);
      helpfunc();
      exit(1);
  }
  return;
}
/*=============================================================================
 GetOutputArgs(): Saisie de les arguments de sortie des resultats d'une recher-
                  che:
                  'style': la valeur par defaut est 'LONG'.
                  'warnings: la valeur par defaut est 'OFF'.
                  'evflag': la valeur par defaut est 'ON'.
                  'histflag': la valeur par defaut est 'OFF'.
=============================================================================*/

void GetOutputArgs(int argc, char *argv[], int *style, int *warnings,
                   int *evflag, int *histflag)
{
  int i;

  *style = LONG;
  *warnings = OFF;
  *evflag = ON;
  *histflag = OFF;

  for (i = 0; i < argc; i++)
  {
      if (argv[i][0] == '-')
      {
          if (strcmp(argv[i], "-short") == 0)    *style = SHORT;
          else
          if (strcmp(argv[i], "-mute") == 0)     *style = MUTE;
          else
          if (strcmp(argv[i], "-warnings") == 0) *warnings = ON;
          else
          if (strcmp(argv[i], "-Eoff") == 0)     *evflag = OFF;
          else
          if (strcmp(argv[i], "-hist") == 0)     *histflag = ON;
      }
  }
  if (*histflag == ON)  *evflag = ON;      /* l'histo a besoin de la E-value */

  return;
}
/*=============================================================================
 GetDirArg(): Saisie de l'argument 'dir' par lecture des options '-rev', '-fwd'
              ou '-fwd+rev', soit: REV_CMPL ou FORWARD ou FWD_AND_REV.
              La valeur par default est FWD_AND_REV.
=============================================================================*/

void GetDirArg(int argc, char *argv[], int *dir)
{
  int i;

  *dir = FWD_AND_REV;

  for (i = 0; i < argc; i++)
  {
      if (argv[i][0] == '-')
      {
          if (strcmp(argv[i], "-rev") == 0) *dir = REV_CMPL;
          else
          if (strcmp(argv[i], "-fwd") == 0) *dir = FORWARD;
      }
  }
  return;
}
/*=============================================================================
 GetStatArg(): Saisie de l'argument 'stat' par lecture des options '-locstat',
               '-globstat' ou '-unifstat', soit: LOCAL GLOBAL ou UNIFORM.
               La valeur par defaut est GLOBAL: la statistique porte sur l'en-
               semble des sequences explorees.
=============================================================================*/

void GetStatArg(int argc, char *argv[], int *stat)
{
  int i;

  *stat = GLOBAL;

  for (i = 0; i < argc; i++)
  {
      if (argv[i][0] == '-')
      {
          if (strcmp(argv[i], "-locstat") == 0) *stat = LOCAL;
          else
          if (strcmp(argv[i], "-unifstat") == 0) *stat = UNIFORM;
      }
  }
  return;
}
/*=============================================================================
 GetSeqArgs(): Saisie des parametres de la selection portant sur les sequences.
               Pour l'interface utilisateur le debut d'une sequence est repere
               par 1, en interne il est mis a 0.
               Pour la sortie de resultats un decallage inverse sera fait.
=============================================================================*/

void GetSeqArgs(int argc, char *argv[],
                int *seqnb1, int *nseq, int *bgn, int *range)
{
  int i;

  *seqnb1 = 1;                                         /* valeurs par defaut */
  *nseq = SEQ_MAX_NB;
  *bgn = 1;
  *range = SEQ_MAX_LEN;

  for (i = 0; i < argc-1; i++)   /* argc-1: ignore les options qui terminent */
  {                                                  /* la ligne de commande */
      if (argv[i][0] == '-')
      {
          if (strcmp(argv[i], "-seq1") == 0) *seqnb1 = atoi(argv[++i]);
          else
          if (strcmp(argv[i], "-nseq") == 0) *nseq = atoi(argv[++i]);
          else
          if (strcmp(argv[i], "-bgn") == 0) *bgn = atoi(argv[++i]);
          else
          if (strcmp(argv[i], "-len") == 0) *range = atoi(argv[++i]);
      }
  }
  if (*seqnb1 <= 0) *seqnb1 = 1;                                 /* securite */
  if (GetFileType(argv[2]) == TRSET) (*seqnb1)++;      /* ignore la sequence */
                                                             /* d'annotation */
  (*bgn)--;                          /* reperage "interne" dans une sequence */
  if (*bgn < 0) *bgn = 0;                                        /* securite */

  return;
}
/*=============================================================================
 GetEnvArgs(): Saisie des variables globales 'logzero' et 'tablen' pour une mo-
               dification des valeurs utilisees par le programme, les valeurs
               par defaut sont fixees dans 'env.c'.
               Enfin, saisie de l'option 'chrono' pour mesure du temps CPU, par
               defaut OFF (inactive).
=============================================================================*/

void GetEnvArgs(int argc, char *argv[], int *chrono)
{
  int i;

  *chrono = OFF;

  for (i = 0; i < argc; i++)
  {
      if (argv[i][0] == '-')
      {
          if (strcmp(argv[i], "-logzero") == 0 && i < argc-1)
              ChLogZero(atof(argv[++i]));
          else
          if (strcmp(argv[i], "-tablen") == 0 && i < argc-1)
              ChScoreTabLen(atoi(argv[++i]));
          else
          if (strcmp(argv[i], "-chrono") == 0) *chrono = ON;
      }
  }
  return;
}
/*=============================================================================
 GetProcArg(): Saisie la variable 'maskproc' qui fixe le type gestion, statique
               ou dynamique, des configurations des masques successifs.
=============================================================================*/

void GetProcArg(int argc, char *argv[], int *maskproc)
{
  int i;

  *maskproc = DYNAMIC;

  for (i = 0; i < argc; i++)
  {
      if (argv[i][0] == '-')
          if (strcmp(argv[i], "-smp") == 0) *maskproc = STATIC;
  }
  return;
}
/*=============================================================================
 GetSumArgs(): saisie des arguments relatifs au traitement des pseudo-comptes
               a l'aide de matrices de substitution.
	       Les options permises sont '-sumf' pour le fichier des matrices
	       de substitution a lire, et '-hpcw' et/ou '-spcw' ou '-pcw' pour
               le poids des pseudo-comptes associes, respectivement, aux heli-
	       ces, aux brins, ou a l'ensemble.
	       Les arguments consideres ici sont tous facultatifs.
	       'pcflag': la valeur par defaut est 'OFF'. Elle ne sera modifiee
	       si un nom de fichier est entre, pour la lecture des matrices.
=============================================================================*/

char *GetSumArgs(int argc, char *argv[], int *pcflag, double *hpcw, double *spcw)
{
  int    i;
  char   *fname = NULL;

  *pcflag = OFF;                                       /* valeurs par defaut */
  *hpcw = *spcw = PSEUDO_COUNTS_FACTOR * PSEUDO_COUNTS_USER ;

  for (i = 0; i < argc-1; i++) /* argc-1: ignore les arguments qui terminent */
  {                                                  /* la ligne de commande */
      if (strcmp(argv[i], "-sumf") == 0) {
          *pcflag = ON;
	  fname = argv[++i];
      }
  }
  if (*pcflag == OFF) {
      *hpcw = *spcw = 0.0;
      return NULL;                          /* pas de pseudo-comptes: sortie */
  }
  fTest(fname);        /* teste la presence du fichier dont le nom est entre */

  for (i = 0; i < argc-1; i++)
  {
      if (argv[i][0] == '-')
      {
          if (strcmp(argv[i], "-pcw") == 0) {
	      *hpcw = *spcw = atof(argv[++i]) * PSEUDO_COUNTS_FACTOR ;
	      break;
	  }
          if (strcmp(argv[i], "-hpcw") == 0) {
	      *hpcw = atof(argv[++i]) * PSEUDO_COUNTS_FACTOR ;
	  }
          if (strcmp(argv[i], "-spcw") == 0) {
	      *spcw = atof(argv[++i]) * PSEUDO_COUNTS_FACTOR ;
	  }
      }
  }
                                                                 /* controle */
  if (*hpcw < 0 || *hpcw > 1) {
      fprintf(stderr, "GetSumArgs: pseudo-counts weight is out of range, exit..\n");
      exit(1);
  }
  if (*spcw < 0 || *spcw > 1) {
      fprintf(stderr, "GetSumArgs: pseudo-counts weight is out of range, exit..\n");
      exit(1);
  }

  return fname;             /* pointe le nom du fichier de donnees en sortie */
}
/*===========================================================================*/

