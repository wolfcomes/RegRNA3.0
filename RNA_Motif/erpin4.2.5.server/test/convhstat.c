
/*=============================================================================
 convhstat.c                                       A.Lambert le 23/04/04

 Histogramme des scores d'une config. d'helices par prod. de convol des distri-
 butions associees aux colonnes des profils.

 Comparaison avec l'histogramme des scores d'une configuration des helices d'un
 masque issus de l'exploration d'une sequence de composition aleatoire.
 On peut introduire dans les profils une ponderation par des matrices de subs-
 titution.

 cc -O3 -Wall -o convhstat convhstat.c -I../include -L../lib -lrnaIV -lm ;
 strip convhstat ;
 chmod 755 convhstat ;

 convhstat <trset-name> <region> <hmask>
         [-Mbs <len>]                  nbre (en Mb) de sites visites, defaut: 1
         [-logzero <logzero>]                                       defaut: -30
         [-dh <dh>]          histo, mesure de la division de l'axe, defaut: 0.1
         [-bkg <%A><%T><%G>[<%C>]]  pourcent. des bases du 'fond', defaut: 25..
	 [-bins <bins>]        nombre de points pour la sortie de l'histogramme
         [-sumf <fname>]         fichier contenant les matrices de substitution
         [-pcw <pcw>]                     poids des pseudo-comptes, defaut: 0.1
         [-hpcw <hpcw>]               poids des pseudo-comptes pour les helices
         [-spcw <spcw>]                 poids des pseudo-comptes pour les brins

 convhstat ~/devc/projets/bio/data/trsetsIII/LET-7.epn -6,6 \
           -umask 6 8 -Mbs 20 -sumf ../sum/SUM.dat -dh 0.05 -bins 120 ;
 gnuplot ;
 set grid ;
 set title \
     "LET-7.epn -6,6 -umask 6 8 -Mbs 20 -sumf SUM.dat -pcw 0.1  -dh 0.05" ;
 set size square ;
 set xrange [-140:0] ;
 plot 'convhstat.dat' u 1:2 with lines, 'randshstat.dat' u 1:2 with boxes ;

 set terminal postscript ;
 set output "hstat.ps" ;
 replot ;
 set terminal x11 ;

=============================================================================*/

#include "rnaIV.h"

double HlxScore(char *seq, Helix *Hlx);
void   PrHistData(char *fname, Histo hist, int gflag);
void   PrHProfs(FILE *txt, double **M, int width, const char *fmt);

/*===========================================================================*/

int main(int argc, char *argv[])
{
  int     datavol = 1000000,                   /* volume suppose des donnees */
          chbinsflag = OFF;      /* signal de modif. du nb de pts de l'histo */
  double  dh = 0.1,          /* mesure de la division elementaire des histo. */
          freqs[4] = {0.25, 0.25, 0.25, 0.25};   /* freqs supposees des ATGC */

  Trset   *trset;
  Pattern *pattern;
  Mask    *mask;
  Histo   hist;
  Map     map;
  float   *scores = NULL;
  int     i, j, nmask, width, nfsc, bins = 100;
  double  Prfsc, **SSProf, min, max, mean, std, score;

  int     pcflag;                 /* variables concernant les pseudo-comptes */
  double  hpcw, spcw;
  char    *sumfname;
  char    *seq = NULL;

  if (argc < 4) {
      fprintf(stderr, "'%s' needs at least 3 arguments, exit..\n", argv[0]);
      exit(1);
  }
  fTest(argv[1]);
  ChLogZero(-30.0);                                     /* valeur par defaut */
  InitStatTables();

  for (i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-logzero") == 0)  ChLogZero(atof(argv[++i]));
      if (strcmp(argv[i], "-dh") == 0)       dh = atof(argv[++i]);
      if (strcmp(argv[i], "-Mbs") == 0)      datavol = atoi(argv[++i])*1000000;
      if (strcmp(argv[i], "-bkg") == 0)  {
          ReadBkgFreqs(freqs, i, argv);  /* pourcentages des bases du "fond" */
          ReInitStatTables(freqs);
      }
      if (strcmp(argv[i], "-bins") == 0) {   /* nb de pts de l'histo produit */
          chbinsflag = ON;
	  bins = atoi(argv[++i]);
      }
  }
  fprintf(stdout, "Bkg ATGC ratios: %.2f %.2f %.2f %.2f\n",
                   freqs[0], freqs[1], freqs[2], freqs[3]);

  sumfname = GetSumArgs(argc, argv, &pcflag, &hpcw, &spcw);
  fprintf(stdout, "H&S pseudo-count weights: %.2f %.2f\n", hpcw, spcw);

  mask = ReadMasksArgs(argc, argv, &nmask);
  if (nmask != 1) {
      fprintf(stderr, "'%s' needs 1 (and only 1) mask arg, exit..\n", argv[0]);
      exit(1);
  }

  trset = ReadTrset(argv[1], 'S', stdout);
  pattern = ReadPattern(argv[2], trset);
  ParseMasksArgs(mask, nmask, pattern);
  GetMask(mask, pattern);
  if (mask->nhx == 0) {
      fprintf(stderr, "0 helix found in pattern (%s), exit..\n", argv[2]);
      exit(1);
  }
  SetNtCodes();        /* tableau de codage des nt pour le calcul des scores */

  if (pcflag) {                                            /* pseudo-comptes */
      ReadSUM(sumfname, hpcw, spcw, trset);
      GetMaskProfilesSM(trset, mask);                            /* Profiles */
  }
  else GetMaskProfiles(trset, mask);                             /* Profiles */

                  /* elements statistiques du profil de structure secondaire */

  SSProf = CatHlxProfiles(mask, &width);
  GetSSProfileStat(SSProf, width, &Prfsc, &min, &max, &mean, &std);
  PrHProfs(stdout, SSProf, width, "%6.2f");       /* visualisation du profil */
  printf("Scmin:  %.3f\n", min);
  printf("Scmax:  %.3f\n", max);

                     /* --- Creation de la sequence aleatoire a explorer --- */

  seq = RandSeq(datavol + mask->max_len, freqs);

                                         /* histogramme des scores d'une cfg */
  hist = SetupHist(min, max, dh, datavol);
  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins-1);        /* etalonnage */

  for (i = nfsc = 0; i < datavol; i++)                        /* exploration */
  {
      int *m;
      score = 0.0;
      for (j = 0, m = mask->hxindex; j < mask->nhx; j++, m++)
      {
          score += HlxScore(seq + i, mask->pattern->hlx + *m);
      }
      if (score > min) {
          nfsc++;
          hist.vals[GetX(map, score)]++;
      }
  }
  hist.samples = nfsc;
  NormalizeHist(hist);                   /* normalisation de l'integrale a 1 */

  if (chbinsflag)
      ResetHistBins(&hist, bins);
  PrHistData("randshstat.dat", hist, ON); /* enregistrement de l'histogramme */
  fprintf(stderr, "%d finite scores recorded\n", nfsc);
  free(hist.vals);

                 /* --- histogramme issu du produit de convolution itere --- */

  hist = GetSSHist(SSProf, width, dh, &Prfsc);
  PrHistData("convhstat.dat", hist, OFF); /* enregistrement de l'histogramme */
  fprintf(stderr, "OK\n");

  FreedMat(SSProf);
  free(hist.vals);
  free(seq);
  free(scores);
  DelTrset(trset);
  DelMasks(mask, nmask);
  DelPattern(pattern);
  FreeNtCodes();

  return 0;
}
/*=============================================================================
 HlxScore(): calcule, depuis la position pointee par 'seq', le score obtenu
             en presentant devant les caracteres de la chaine le profil de
             l'helice pointee par 'Hlx', la distance entre les 2 brins de l'he-
             lice est fixee a Hlx->min_dist.
             Variante de 'GetHlxScores' dans 'libsrc/scores.c'.
=============================================================================*/

double HlxScore(char *seq, Helix *Hlx)
{
  int     i, j;
  char    *right;
  double  score;
  extern  short **NtHlxCode;               /* declare dans 'libsrc/ntcode.c' */

  right = seq + Hlx->min_len - 1;

  for (j = 0, score = 0.0; j < Hlx->helix_len; j++, right--)
  {
      i = NtHlxCode[(int) seq[j]][(int) *right];
      score += Hlx->Profile[i][j];
  }
  return score;
}
/*=============================================================================
 PrHistData(): Enregistre dans le fichier identifie par 'fname' les donnees,
               abscisses et ordonnees de l'histogramme 'hist' dont l'integrale
               est supposee etre egal a 1 (pas de controle).
               Si l'argument 'gflag' vaut 1 une gaussienne de meme moyenne et
               ecart type que 'hist' est enregistree en 3eme colonne pour com-
	       paraison.
=============================================================================*/

void PrHistData(char *fname, Histo hist, int gflag)
{
  FILE   *txt;
  Map     map;
  double  mean, std;
  int     i;

  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins-1);

  txt = fopen(fname, "w");

  if (gflag == ON)
  {
      GetHistStat(hist, &mean, &std);

      for (i = 0; i < hist.bins; i++) {
          double x, y;
          x = Getx(map, i);
          y = Gauss(x, mean, std);
          fprintf(txt, "%.4e %.4e %.4e\n", x, hist.vals[i], y);
      }
  }
  else
      for (i = 0; i < hist.bins; i++)
          fprintf(txt, "%.4e %.4e\n", Getx(map, i), hist.vals[i]);

  fclose(txt);
  return;
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
/*===========================================================================*/
