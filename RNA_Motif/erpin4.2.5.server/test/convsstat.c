
/*=============================================================================
 convsstat.c                                       A.Lambert le 23/04/04

 Test sur l'histogramme des scores de brins sans gaps par prod. de convolution
 des distributions associees aux colonnes des profils.

 Comparaison avec l'histogramme des scores d'une configuration des brins d'un
 masque issus de l'exploration d'une sequence de composition aleatoire.
 On peut introduire dans les profils une ponderation par des matrices de subs-
 titution.

 cc -O3 -Wall -o convsstat convsstat.c -I../include -L../lib -lrnaIV -lm ;
 strip convsstat ;
 chmod 755 convsstat ;

 convsstat <trset-name> <region> <hmask>
         [-Mbs <len>]                  nbre (en Mb) de sites visites, defaut: 1
         [-logzero <logzero>]                                       defaut: -30
         [-dh <dh>]          histo, mesure de la division de l'axe, defaut: 0.1
         [-bkg <%A><%T><%G>[<%C>]]  pourcent. des bases du 'fond', defaut: 25..
	 [-bins <bins>]        nombre de points pour la sortie de l'histogramme
         [-sumf <fname>]         fichier contenant les matrices de substitution
         [-pcw <pcw>]                     poids des pseudo-comptes, defaut: 0.1
         [-hpcw <hpcw>]               poids des pseudo-comptes pour les helices
         [-spcw <spcw>]                 poids des pseudo-comptes pour les brins

 convsstat ~/devc/projets/bio/data/trsets.serv/polya_signal.epn 11,7 \
           -umask 11 4 5 6 7 -Mbs 30 -sumf ../sum/SUM.dat ;
 convsstat ~/devc/projets/bio/data/trsets.serv/polya_signal.epn 11,5 \
           -umask 11 -Mbs 10 -sumf ../sum/SUM.dat ;

 convsstat ~/devc/projets/bio/data/trsets.serv/polya_signal.epn 11,5 \
           -umask 11 4 5 -Mbs 20 -sumf ../sum/SUM.dat -dh 0.05 -bins 120;
 gnuplot ;
 set grid ;
 set title \
     "polya_signal.epn 11,5 -umask 11 4 5 -Mbs 20 -sumf SUM.dat -dh 0.05" ;
 set size square ;
 plot 'convsstat.dat' u 1:2 with lines, 'randsstat.dat' u 1:2 with boxes ;

 set terminal postscript ;
 set output "sstat.ps" ;
 replot ;
 set terminal x11 ;

=============================================================================*/

#include "rnaIV.h"

void   PrHistData(char *fname, Histo hist, int gflag);
void   PrHProfs(FILE *txt, double **M, int width, const char *fmt);
double StsNoGScore(char *seq, double **Prof, int width);

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
  int     i, nmask, nfsc, width, bins = 100;
  double  **StsProf, Prfsc, min, max, score;

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
  if (mask->nst == 0) {
      fprintf(stderr, "0 strand found in pattern (%s), exit..\n", argv[2]);
      exit(1);
  }
  SetNtCodes();        /* tableau de codage des nt pour le calcul des scores */

  if (pcflag) {                                            /* pseudo-comptes */
      ReadSUM(sumfname, hpcw, spcw, trset);
      GetMaskProfilesSM(trset, mask);                            /* Profiles */
  }
  else GetMaskProfiles(trset, mask);                             /* Profiles */

	         /* --- histogramme issu du produit de convolution itere --- */

  if ((StsProf = CatStNoGProfs(mask, &width)) == NULL) {
      fprintf(stderr, "no strand without gaps in mask, exit..\n");
      exit(1);
  }
  hist = StsNoGHist(StsProf, width, dh, &Prfsc);
  PrHistData("convsstat.dat", hist, OFF); /* enregistrement de l'histogramme */
  min = hist.hmin;
  max = hist.hmax;
  printf("Scmin:  %.3f\n", min);
  printf("Scmax:  %.3f\n", max);
  free(hist.vals);

                     /* --- Creation de la sequence aleatoire a explorer --- */

  seq = RandSeq(datavol + width, freqs);

                                       /* ----- histogramme des scores ----- */

  hist = SetupHist(min, max, dh, datavol);
  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins-1);        /* etalonnage */

  for (i = nfsc = 0; i < datavol; i++)                        /* exploration */
  {
      score = StsNoGScore(seq + i, StsProf, width);
      if (score > min) {
          nfsc++;
          hist.vals[GetX(map, score)]++;
      }
  }
  hist.samples = nfsc;
  NormalizeHist(hist);                   /* normalisation de l'integrale a 1 */

  if (chbinsflag)
      ResetHistBins(&hist, bins);
  PrHistData("randsstat.dat", hist, ON);  /* enregistrement de l'histogramme */
  fprintf(stderr, "%d finite scores recorded\n", nfsc);
  free(hist.vals);

  fprintf(stderr, "OK\n");

  FreedMat(StsProf);
  free(seq);
  free(scores);
  DelTrset(trset);
  DelMasks(mask, nmask);
  DelPattern(pattern);
  FreeNtCodes();

  return 0;
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

  if (gflag == 1)
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
/*=============================================================================
 StsNoGScore(): Pour un groupe de brins de profil pointe par 'Prof' ne conte-
                nant pas de gaps, cette fonction calcule et retourne le score
        	obtenu par 'seq' sur une longueur egale a 'width', le nombre de
		colonnes du profil.
		Version simplifiee reservee a un usage statistique.
=============================================================================*/

double StsNoGScore(char *seq, double **Prof, int width)
{
  int    i, j;
  double score;
  extern short *NtStCode;                  /* declare dans 'libsrc/ntcode.c' */

  for (j = 0, score = 0.0; j < width; j++)
  {
      i = NtStCode[(int) seq[j]];
      score += Prof[i][j];
  }
  return score;
}
/*===========================================================================*/
