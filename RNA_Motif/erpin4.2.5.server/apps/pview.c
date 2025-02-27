
/*=============================================================================
 pview.c                                        A.Lambert le 01/04/04

 Visualisation des profils des helices et des brins d'un masque: ce programme
 permet d'observer l'incidence de l'introduction de pseudo-comptes sur le con-
 tenu des profils.

 cc -O3 -Wall -o pview pview.c -I../include -L../lib -lrnaIV -lm ;
 chmod 755 pview ;
 strip pview ;
 mv pview ../bin ;

 pview <trset> <region> <masque>
       [-logzero <logzero>]                                         defaut: -30
       [-bkg <%A><%T><%G>[<%C>]]    pourcent. des bases du 'fond', defaut: 25..
       [-sumf <fname>]           fichier contenant les matrices de substitution
       [-pcw <pcw>]                       poids des pseudo-comptes, defaut: 0.1
       [-hpcw <hpcw>]                 poids des pseudo-comptes pour les helices
       [-spcw <spcw>]                   poids des pseudo-comptes pour les brins

 pview ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -2,2 -nomask \
       -sumf ../sum/SUM.dat -pcw 0.1 > pview.dat ;
 pview ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -2,2 -nomask > pview0.dat ;

 pview ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -6,8 -nomask \
       -sumf ../sum/SUM.dat -hpcw 0.1 -spcw 0.0 > pview.dat ;
 pview ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -6,8 -nomask > pview0.dat ;

 pview ~/devc/projets/bio/data/trsetsIII/LET-7.epn -2,2 -umask 2 4 6 8 \
       -sumf ../sum/SUM.dat -pcw 0.05 > pview.dat ;

=============================================================================*/

#include "rnaIV.h"

void   PrHProfs(FILE *txt, double **M, int width, const char *fmt);
void   PrSProfs(FILE *txt, double **M, int width, const char *fmt);

/*===========================================================================*/

int main(int argc, char *argv[])
{
  double freqs[4] = {0.25, 0.25, 0.25, 0.25};    /* freqs supposees des ATGC */

  Trset   *trset;
  Pattern *pattern;
  Mask    *mask;
  double  **profs;
  int     i, nmask, width;

  int     pcflag;                 /* variables concernant les pseudo-comptes */
  double  hpcw, spcw;
  char    *sumfname;

  if (argc < 4) {
      fprintf(stderr, "'%s' needs at least 3 arguments, exit..\n", argv[0]);
      exit(1);
  }
  fTest(argv[1]);
  ChLogZero(-30.0);                                     /* valeur par defaut */
  InitStatTables();

  for (i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-logzero") == 0)  ChLogZero(atof(argv[++i]));
      if (strcmp(argv[i], "-bkg") == 0)  {
          ReadBkgFreqs(freqs, i, argv);  /* pourcentages des bases du "fond" */
          ReInitStatTables(freqs);
      }
  }
  fprintf(stdout, "Bkg ATGC ratios: %.2f %.2f %.2f %.2f\n",
                   freqs[0], freqs[1], freqs[2], freqs[3]);

  sumfname = GetSumArgs(argc, argv, &pcflag, &hpcw, &spcw);
  fprintf(stdout, "H&S pseudo-count weights: %.2e %.2e\n", hpcw, spcw);

  mask = ReadMasksArgs(argc, argv, &nmask);
  if (nmask != 1) {
      fprintf(stderr, "'%s' needs 1 (and only 1) mask arg, exit..\n", argv[0]);
      exit(1);
  }

  trset = ReadTrset(argv[1], 'S', stdout);
  pattern = ReadPattern(argv[2], trset);
  ParseMasksArgs(mask, nmask, pattern);
  GetMask(mask, pattern);

  if (pcflag) {                                            /* pseudo-comptes */
      ReadSUM(sumfname, hpcw, spcw, trset);
      GetMaskProfilesSM(trset, mask);                             /* Profiles */
  }
  else GetMaskProfiles(trset, mask);                             /* Profiles */

  if (mask->nhx != 0) {
      profs = CatHlxProfiles(mask, &width);             /* profils d'helices */
      PrHProfs(stdout, profs, width, "%6.2f");
      FreedMat(profs);
  }
  if (mask->nst != 0) {
      profs = CatStProfiles(mask, &width);               /* profils de brins */
      PrSProfs(stdout, profs, width, "%6.2f");
      FreedMat(profs);
  }
  DelTrset(trset);
  DelMasks(mask, nmask);
  DelPattern(pattern);

  return 0;
}
/*=============================================================================
 PrHProfs(): Affiche dans le fichier pointe par 'txt' le tableau 'M' de largeur
             'width', des profils d'helices.
	     Le format d'affichage des elements est pointe par 'fmt'.
=============================================================================*/

void PrHProfs(FILE *txt, double **M, int width, const char *fmt)
{
  int  i, j;
  char AUGC[ALPHA_LEN] = "AUGC";

  fprintf(txt, "\n");

  for (i = 0; i < SQR_ALPHA_LEN; )
  {
      fprintf(txt, " %c%c ", AUGC[i / ALPHA_LEN], AUGC[i % ALPHA_LEN]);
      for (j = 0; j < width; j++) {
          fprintf(txt, " ");          /* un espace pour separer les elements */
          fprintf(txt, fmt, M[i][j]);
      }
      fprintf(txt, "\n");
      if ((++i) % ALPHA_LEN == 0) fprintf(txt, "\n");
  }
  fprintf(txt, "\n");
}
/*=============================================================================
 PrSProfs(): Affiche dans le fichier pointe par 'txt' le tableau 'M' de largeur
             'width', des profils de brins.
	     Le format d'affichage des elements est pointe par 'fmt'.
=============================================================================*/

void PrSProfs(FILE *txt, double **M, int width, const char *fmt)
{
  int  i, j;
  char AUGC[ALPHA_LEN] = "AUGC";

  fprintf(txt, "\n");

  for (i = 0; i < ALPHA_LEN; i++)
  {
      fprintf(txt, " %c ", AUGC[i % ALPHA_LEN]);
      for (j = 0; j < width; j++) {
          fprintf(txt, " ");          /* un espace pour separer les elements */
          fprintf(txt, fmt, M[i][j]);
      }
      fprintf(txt, "\n");
  }

  fprintf(txt, " - ");
  for (j = 0; j < width; j++) {
      fprintf(txt, " ");          /* un espace pour separer les elements */
      fprintf(txt, fmt, M[ALPHA_LEN][j]);
  }
  fprintf(txt, "\n\n");
}
/*===========================================================================*/
