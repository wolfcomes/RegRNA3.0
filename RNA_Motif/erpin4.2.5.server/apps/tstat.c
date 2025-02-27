
/*=============================================================================
 tstat.c                           A.Lambert le 03/12/01  revu le 13/01/05

 Calcule dans un training-set, les scores d'une region et d'un nombre arbitrai-
 re de masques associes a cette region, ainsi que la statistique des scores.
 Ce programme devrait aider au choix des seuils en vue d'un run de erpin.

 L'option [-tbkg|-bkg <%A> <%T> <%G>] fixe les frequences du "fond", res-
 pectivement aux valeurs mesurees dans la base d'entrainement (-tbkg) ou en-
 trees explicitement (-bkg).
 Par defaut une distribution uniforme des bases sera utilisee.

 cc -O3 -Wall -o tstat tstat.c -I../include -L../lib -lrnaIV -lm ;
 strip tstat ;
 chmod 755 tstat ;
 mv tstat ../bin ;

 tstat <trset> <region> <mask>
       [(-mask|-umask|-add) <arg1> <arg2> ..][..]
       [-tbkg|-bkg <%A><%T><%G>]
       [-sumf <fname>]           fichier contenant les matrices de substitution
       [-pcw <pcw>]                       poids des pseudo-comptes, defaut: 0.1
       [-hpcw <hpcw>]                 poids des pseudo-comptes pour les helices
       [-spcw <spcw>]                   poids des pseudo-comptes pour les brins

 tstat ~/devc/projets/bio/data/trsetsII/trna.db  -7,7  -umask 7 -umask 8 \
       -sumf ../sum/SUM.dat -pcw 0.25 ;
 tstat ~/devc/projets/bio/data/trsetsII/trna.db  -7,7  -umask 7 -umask 8 -tbkg ;
 tstat ~/devc/projets/bio/data/trsetsII/trna.db  -7,7  -umask 7 -umask 8 \
       -bkg 28 30 20 ;
 tstat ~/devc/projets/bio/data/trsetsII/trna.db  -7,7  -mask 8 ;
 tstat ~/devc/projets/bio/data/trsetsII/23S301proPK.db  1,8 -umask 2 4 8 ;
 tstat ~/devc/projets/bio/data/trsetsII/mirna1.epn -2,2 \
       -umask 10 13 15 17 -add 8 11 19 -add 6 9 21 -nomask ;

=============================================================================*/

#include "rnaIV.h"

void TrsetHelp(void);

/*===========================================================================*/

int main(int argc, char *argv[])
{
  Trset   *trset;
  Pattern *pattern;
  Mask    *mask;
  double  freqs[4];
  int     i, nmask;

  int       pcflag;          /* variables concernant les pseudo-comptes */
  double    hpcw, spcw;
  char      *sumfname;

  ReadHelpArgs(argc, argv, 3, TrsetHelp);
  fTest(argv[1]);
  sumfname = GetSumArgs(argc, argv, &pcflag, &hpcw, &spcw);

  mask = ReadMasksArgs(argc, argv, &nmask);
  trset = ReadTrset(argv[1], 'S', stdout);
  pattern = ReadPattern(argv[2], trset);
  ParseMasksArgs(mask, nmask, pattern);
  InitStatTables();
  SetNtCodes();        /* tableau de codage des nt pour le calcul des scores */

  if (pcflag) {                                            /* pseudo-comptes */
      fprintf(stdout, "Pseudo-counts weights:  %.2f  %.2f\n", hpcw, spcw);
      ReadSUM(sumfname, hpcw, spcw, trset);
      GetPatternProfilesSM(trset, pattern);                      /* Profiles */
  }
  else  GetPatternProfiles(trset, pattern);

                                      /* selection de la statistique du fond */
  for (i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-tbkg") == 0) {
          fGetStat(argv[1], 1, SEQ_MAX_NB, 'V', stdout);
          ChPatternStat(pattern);
          break;
      }
      if (strcmp(argv[i], "-bkg") == 0) {
          ReadBkgFreqs(freqs, i, argv);
          ResetBkgFreqs(freqs);
          ChPatternStat(pattern);               /* actualisation des profils */
          break;
      }
  }
  fPrintfPatternTStat(trset, pattern, stdout);

  GetMasks(mask, nmask, pattern);

  for (i = 0; i < nmask; i++)
      if (mask[i].mode != NOMASK) fPrintfMaskTStat(trset, mask + i, stdout);

  FreeNtCodes();
  DelTrset(trset);
  DelMasks(mask, nmask);
  DelPattern(pattern);

  exit(0);
}
/*=============================================================================
 TrsetHelp(): help
=============================================================================*/

void TrsetHelp(void)
{
  fprintf(stderr,

"tstat: Training set scores STATistics\n\n"
"Usage:\n"
"tstat [-h]\n"
"tstat <training-set> <region> <mask>\n"
"      [(-mask|-umask|-add) <arg1> <arg2> ..][..]\n"
"      [-tbkg|-bkg <%%A><%%T><%%G>]       background percents, default: 25 ..\n"
"      [-sumf <fname>]           substitution matrix file name, default: none\n"
"      [-pcw <pcw>]                        pseudo-counts weight, default: 0.1\n"
"      [-hpcw <hpcw>]                               pseudo-counts for helices\n"
"      [-spcw <spcw>]                               pseudo-counts for strands\n"
);

  return;
}
/*===========================================================================*/
