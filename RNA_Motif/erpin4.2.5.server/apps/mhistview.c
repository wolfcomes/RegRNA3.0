
/*=============================================================================
 mhistview.c                                       A.Lambert le 07/04/04

 On calcule dans ce programme l'histogramme des scores par iteration du produit
 de convolution, et l'histogramme des scores de la base d'entrainement.

 On visualise les histogrammes des scores du masque passe en argument, calcules
 sur la base d'entrainement et depuis les profils statistiques simulant des se-
 quences aleatoires.
 On a la possibilite d'introduire des pseudo-comptes par l'intermediaire de ma-
 trices de substitution des brins et helices. Cette operation modifie les pro-
 fils, notamment supprime les (ou des) exclusions.
 La E-value issue du calcul est aussi sortie.

 cc -O3 -Wall -o mhistview mhistview.c -I../include -L../lib -lrnaIV -lm ;
 strip mhistview ;
 chmod 755 mhistview ;
 mv mhistview ../bin ;

 mhistview <trset-name> <region> <mask>
           [-dh <dh>]      histo, mesure de la division sur l'axe, defaut: 0.05
           [-logzero <logzero>]                                     defaut: -30
           [-spc <spc>]               tirages aleatoires: nombre d'echantillons
                                    pour la 1ere colonne de profil, defaut: 300
           [-cutoff <cutoff>]  alternative a la val. par defaut: Tmean - 2*Tstd
           [-bkg <%A><%T><%G>[<%C>]]  pourcent. des bases du 'fond', defaut: 25
	   [-Mbs <Mbs>]       nbre suppose (en Mb) de sites visites, defaut: 10

           [-sumf <fname>]       fichier contenant les matrices de substitution
           [-pcw <pcw>]                   poids des pseudo-comptes, defaut: 0.1
           [-hpcw <hpcw>]             poids des pseudo-comptes pour les helices
           [-spcw <spcw>]               poids des pseudo-comptes pour les brins

 mhistview ~/devc/projets/bio/data/trsetsIII/LET-7.epn -2,2 -umask 2 4 6 8 3 \
           -Mbs 100 -cutoff -39 -sumf ../sum/SUM.dat ;
================================================================
E-value at cutoff -39.0 for 100.0Mb single strand data: 4.63e+02
================================================================
prendre aussi pour cutoff: -30  -22

 gnuplot ;
 set grid ;
 set xrange [-200: 50] ;
 set yrange [0: 0.03] ;
 plot 'Thistview.dat' u 1:2 with boxes, 'Rhistview.dat' u 1:2 with lines ;

=============================================================================*/

#include "rnaIV.h"

float  *MaskTScores(Trset *trset, Mask *mask);
void   PrHistData(char *fname, Histo hist, int gflag);

/*===========================================================================*/

int main(int argc, char *argv[])
{
  int     spercol = SAMPLES1_SINGLE;   /* echant. pour la 1ere col. de prof. */
  double  dh = 0.05,         /* mesure de la division elementaire des histo. */
          datavol = 10.,                 /* volume suppose des donnees en Mb */
          freqs[4] = {0.25, 0.25, 0.25, 0.25};   /* freqs supposees des ATGC */

  Trset   *trset;
  Pattern *pattern;
  Mask    *mask;
  Histo   hist, cdf;
  float   *scores = NULL;
  int     i, nmask;

  int     pcflag;                 /* variables concernant les pseudo-comptes */
  double  hpcw, spcw;
  char    *sumfname;

  double  Tmean, Tstd, Rmean, Rstd, D, Prfsc, cutoff, Prob, cfgs, E, L;


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
      if (strcmp(argv[i], "-spc") == 0)      spercol = atoi(argv[++i]);
      if (strcmp(argv[i], "-Mbs") == 0)      datavol = atof(argv[++i]);
      if (strcmp(argv[i], "-bkg") == 0)  {
          ReadBkgFreqs(freqs, i, argv);  /* pourcentages des bases du "fond" */
          ReInitStatTables(freqs);
      }
  }
  fprintf(stdout, "\nBkg ATGC ratios: %.2f %.2f %.2f %.2f\n",
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
  SetNtCodes();        /* tableau de codage des nt pour le calcul des scores */
  
  if (pcflag) {                                            /* pseudo-comptes */
      ReadSUM(sumfname, hpcw, spcw, trset);
      GetMaskProfilesSM(trset, mask);                            /* Profiles */
  }
  else GetMaskProfiles(trset, mask);                             /* Profiles */

           /* ---- histogramme des scores du masque sur le training set ---- */

  scores = MaskTScores(trset, mask);
  hist = GetNHist(scores, trset->nseq, dh);
  GetHistStat(hist, &Tmean, &Tstd);

  printf("Tmin:  %.2f\n", hist.hmin);
  printf("Tmax:  %.2f\n", hist.hmax);
  printf("Tmean: %.2f\n", Tmean);
  printf("Tstd : %.2f\n", Tstd);

  PrHistData("Thistview.dat", hist, 0);   /* enregistrement de l'histogramme */
  free(hist.vals);

        /* ---- histogramme des scores du masque sur sequence aleatoire ---- */

  hist = GetMaskHist(mask, dh, &Prfsc, spercol);
  GetHistStat(hist, &Rmean, &Rstd);
  D = (Tmean - Rmean)/(Tstd + Rstd);

  printf("Rmin:  %.2f\n", hist.hmin);
  printf("Rmax:  %.2f\n", hist.hmax);
  printf("Rmean: %.2f\n", Rmean);
  printf("Rstd : %.2f\n", Rstd);
  printf("D = (Tmean - Rmean)/(Tstd + Rstd) = %.2f\n\n", D);

  PrHistData("Rhistview.dat", hist, 0);   /* enregistrement de l'histogramme */

                  /* ------ calcul de la E-value pour un cutoff donne ------ */

  cfgs = rint(pow(2.0, mask->log2ncfg));   /* nbre de config. voir 'GetMask' */
  cdf = GetHistCdf(hist);
  cutoff = Tmean - 2.0*Tstd;                /* valeur par defaut de 'cutoff' */

  for (i = 1; i < argc; i++)         /* modification optionnelle de 'cutoff' */
      if (strcmp(argv[i], "-cutoff") == 0)
      {
          double tmp = atof(argv[++i]);
          if (tmp < cdf.hmin || tmp >= cdf.hmax)
              fprintf(stdout,
              "WARNING: 'cutoff' input is out of range, use default..\n\n");
          else cutoff = tmp;
      }

  L = 1.e+6 * datavol;                      /* volume, en bases, des donnees */
  Prob = Prfsc * Interpol(cdf, cutoff);              /* Prob(score > cutoff) */
  E = L * evprob(Prob, (unsigned int) cfgs);                      /* E-value */

  fprintf(stdout, "cfgs = %.2e\n", cfgs);
  fprintf(stdout, "cutoff = %.2f\n", cutoff);
  fprintf(stdout, "Prob(score > -inf) = %.2e\n", Prfsc);
  fprintf(stdout, "================================================================\n");
  fprintf(stdout, "E-value at cutoff %.1f for %.1fMb single strand data: %.2e\n",
                  cutoff, datavol, E);
  fprintf(stdout, "================================================================\n");

                                                                   /* sortie */
  free(hist.vals);
  free(scores);
  DelTrset(trset);
  DelMasks(mask, nmask);
  DelPattern(pattern);

  return 0;
 }
/*=============================================================================
 MaskTScores(): Calcule, pour l'ensemble des sequences de la base d'entraine-
                ment pointee par 'trset', le score du pattern pointe par 'mask'
                Retourne un pointeur sur le tableau des scores calcules.
	        Ceux-ci ne sont pas reordonnes comme dans 'GetMaskTScores'.
=============================================================================*/

float *MaskTScores(Trset *trset, Mask *mask)
{
  int   i;
  float *scores;

  scores = (float *) malloc(trset->nseq * sizeof(float));

  for (i = 0; i < trset->nseq; i++)              /* boucle sur les sequences */
      scores[i] = (float) GetMaskTScore(trset->data[i], mask);

  return scores;
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
/*===========================================================================*/
