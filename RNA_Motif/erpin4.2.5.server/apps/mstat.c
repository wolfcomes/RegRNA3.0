
/*=============================================================================
 mstat.c    --  Mask Statistics  --    A.Lambert le 12/04/03 revu le 04/05/04


 Le programme calcule la E-value d'un motif (un seul masque !): le nombre de
 scores du motif superieurs a un seuil donne attendu dans une sequence aleatoi-
 re de taille donnee.


 Des informations statistiques sont affichees notamment:
 - La E-value des scores > cutoff,
 - Le "score discriminant" (Tmean - Rmean) / (Tstd + Rstd) qui relie les moyen-
   ne et dev. standard des scores du Training-set (T) et scores aleatoires (R).

 - L'option '-hist' permet l'enregistrement dans un fichier de l'histogramme
   de l'ensemble des scores finis du motif, une 3eme colonne donne les ordon-
   nees d'une gaussienne de parametres (mean, std) identiques pour comparaison.
 - L'option '-evals' permet l'enregistrement dans un fichier des valeurs de la
   E-value versus les scores > cutoff sur un intervalle allant jusqu'au score
   maximal de la base d'entrainement.

 - Par defaut la distribution des ATGC du "fond" est uniforme, mais le prog-
   ramme accepte en argument une distribution quelconque (option '-bkg').
 - Le nombre d'echantillons pour la statistique des brins est proportionnel au
   carre du nombre de colonnes, l'option '-spc' (defaut: 300) en fixe le nom-
   bre pour la 1ere colonne.

 make -f apps.mk mstat ;
 mv mstat ../bin ;

 cc -O3 -Wall -o mstat mstat.c -I../include -L../lib -lrnaIV -lm ;
 strip mstat ;
 chmod 755 mstat ;

 mstat <trset> <region> <mask>     3 arguments obligatoires, 1 masque et 1 seul
       [-logzero <logzero>]                                         defaut: -30
       [-dh <dh>]           histo, division elementaire sur l'axe, defaut: 0.05
       [-spc <spc>]              tirages aleatoires: nombre d'echantillons pour
                                        la 1ere colonne de profil, defaut:  300
       [-cutoff <cutoff>]      alternative a la val. par defaut: Tmean - 2*Tstd
       [-Mbs <Mbs>]           nbre suppose (en Mb) de sites visites, defaut: 10
       [-bkg <%A><%T><%G>[<%C>]] pourcentages des bases du 'fond', defaut: 25..
       [-hist]         cree un fichier contenant l'histogramme des scores finis
       [-evals]             cree un fichier des 'E-values vs scores' (> cutoff)

       [-sumf <fname>]           fichier contenant les matrices de substitution
       [-pcw <pcw>]                       poids des pseudo-comptes, defaut: 0.1
       [-hpcw <hpcw>]                 poids des pseudo-comptes pour les helices
       [-spcw <spcw>]                   poids des pseudo-comptes pour les brins

 mstat  ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -8,8 \
        -umask 8 -logzero -50                                          1 helice
 mstat  ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -6,8 \
        -umask 6 8 -logzero -50                                       2 helices
 mstat  ~/devc/projets/bio/data/trsetsII/trna-typeI.epn 5,5 \
        -umask 5 -logzero -50                                   1 brin avec gap
 mstat  ~/devc/projets/bio/data/trsetsII/trna-typeI.epn 9,9 \
        -umask 9 -logzero -50                                   1 brin sans gap

 mstat  ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -4,4 \
        -nomask
 mstat  ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -4,8 \
        -nomask -hist -evals                                                OK
 mstat  ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -2,2 \
        -nomask -logzero -50                                     trna 'complet'
 mstat  ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -2,2 \
        -nomask -logzero -50 -bkg 20 26 22
 mstat  ~/devc/projets/bio/data/trsetsIII/uptake-ire.epn 3,3 \
        -nomask -logzero -50 -spc 1000 -Mbs 100     brin sans gap de longeur 1

 mstat ~/devc/projets/bio/data/trsetsIII/23S301proPK.epn 1,8 -nomask -Mbs 100

 gnuplot ;
 set grid ;
 plot 'mhisto.dat' u 1:2 with boxes, 'mhisto.dat' u 1:3 with lines ;
 plot 'mlogev.dat' u 1:2 with lines, 'mlogev.dat' u 1:3 with lines ;

=============================================================================*/

#include "rnaIV.h"

#define  H_OUTPUT  "mhisto.dat"
#define  E_OUTPUT  "mevals.dat"

void  PrTStat(FILE *txt, double min, double max, double mean, double std);
void  PrRStat(FILE *txt, Histo hist, double mean, double std);
void  PrHist(char *fname, Histo hist);
void  PrEval(char *fname, Histo hist, double min, double max);
int   PrLogEval(char *fname, Histo hist,
                double min, double max, double *A, double *B);

void  PrMaskTStat(FILE *txt, Trset *trset, Mask *mask,
                  double *min, double *max , double *mean , double *std);
void  MstatHelp(void);

/*===========================================================================*/

int main(int argc, char *argv[])
{
                                   /* parametres modifies par arg. de 'main' */
  int     spercol = SAMPLES1_SINGLE;   /* echant. pour la 1ere col. de prof. */
  double  dh = 0.05,         /* mesure de la division elementaire des histo. */
          datavol = 10.,                 /* volume suppose des donnees en Mb */
          freqs[4] = {0.25, 0.25, 0.25, 0.25};   /* freqs supposees des ATGC */

  Trset   *trset;
  Pattern *pattern;
  Mask    *mask;
  Histo   hist, cdf;
  int     i, nmask;
  double  Prfsc,                           /* probabilite des scores > - Inf */
          Tmin, Tmax, Tmean, Tstd, Rmean, Rstd,
          cutoff, Prob, cfgs,
          duration,
          E, L;
  unsigned long bgntime;

  int     pcflag;                 /* variables concernant les pseudo-comptes */
  double  hpcw, spcw;
  char    *sumfname;

  ReadHelpArgs(argc, argv, 3, MstatHelp);
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

  trset = ReadTrset(argv[1], 'S', stderr);
  pattern = ReadPattern(argv[2], trset);
  ParseMasksArgs(mask, nmask, pattern);
  GetMask(mask, pattern);
  SetNtCodes();        /* tableau de codage des nt pour le calcul des scores */

  if (pcflag) {                                            /* pseudo-comptes */
      ReadSUM(sumfname, hpcw, spcw, trset);
      GetMaskProfilesSM(trset, mask);                            /* Profiles */
  }
  else GetMaskProfiles(trset, mask);                             /* Profiles */

  PrMaskTStat(stdout, trset, mask, &Tmin, &Tmax, &Tmean, &Tstd);

  DelTrset(trset);

                    /* ---------- fin de l'utilisation de 'trset' ---------- */

  bgntime = clock();                             /* demarrage du chronometre */

  cfgs = rint(pow(2.0, mask->log2ncfg));   /* nbre de config. voir 'GetMask' */

  hist = GetMaskHist(mask, dh, &Prfsc, spercol);
  GetHistStat(hist, &Rmean, &Rstd);
  PrRStat(stdout, hist, Rmean, Rstd);

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
                       /* ------------- calcul de la E-value --------------- */

  L = 1.e+6 * datavol;                      /* volume, en bases, des donnees */
  Prob = Prfsc * Interpol(cdf, cutoff);              /* Prob(score > cutoff) */
  E = L * evprob(Prob, (unsigned int) cfgs);                      /* E-value */


  if (Tstd + Rstd > 1.e-8)
  {
      double dsc = (Tmean - Rmean) / (Tstd + Rstd);    /* score discriminant */
      fprintf(stdout, "D = (Tm - Rm)/(Ts + Rs) = %.2f\n", dsc);
  }
  fprintf(stdout, "cfgs = %.2e\n", cfgs);
  fprintf(stdout, "cutoff = %.2f\n", cutoff);
  fprintf(stdout, "Prob(score > -inf) = %.2e\n", Prfsc);
  fprintf(stdout, "================================================================\n");
  fprintf(stdout, "E-value at cutoff %.1f for %.1fMb single strand data: %.2e\n",
                  cutoff, datavol, E);
  fprintf(stdout, "================================================================\n");

                     /* ---------- sortie optionnelle des courbes ---------- */

  for (i = 1; i < argc; i++)
  {
      if (strcmp(argv[i], "-hist") == 0)           /* histo des scores finis */
          PrHist(H_OUTPUT, hist);
      if (strcmp(argv[i], "-evals") == 0)     /* logE versus scores > cutoff */
      {
          double min = cutoff, max = hist.hmax ;

          if (Tstd + Rstd > 1.e-8)
          {
              for (i = 0; i < cdf.bins; i++) {
                  Prob = Prfsc*cdf.vals[i];
                  cdf.vals[i] = L * evprob(Prob, (unsigned int) cfgs);
              }
              PrEval(E_OUTPUT, cdf, min, max);
          }
      }
  }
  duration = (double) (clock() - bgntime) / (double) CLOCKS_PER_SEC;
  fprintf(stdout, "\ncputime = %.2f sec.\n", duration);

  free(hist.vals);
  free(cdf.vals);
  DelMasks(mask, nmask);
  DelPattern(pattern);
  return 0;
}
/*=============================================================================
 PrTStat(): Affiche sur 'txt' les elements statistiques du Trset.
            Note: 'Tmean - (3/2)*Tstd' capture la majeure partie des motifs de
            la base d'entrainement.
=============================================================================*/

void PrTStat(FILE *txt, double min, double max, double mean, double std)
{
  fprintf(txt, "Scores statistics from the training set:\n");
  fprintf(txt, "Tmin:   %.2f \n", min);
  fprintf(txt, "Tmax:   %.2f \n", max);
  fprintf(txt, "Tmean:  %.2f \n", mean);
  fprintf(txt, "Tstd:   %.2f \n", std);
  fprintf(txt, "Tmean - 2*Tstd:   %.2f (default value of 'cutoff')\n\n",
          mean - 2.0*std);
  return;
}
/*=============================================================================
 PrRStat(): Affiche sur 'txt' les elements statistiques des scores aleatoires
            concernant l'histogramme 'hist'.
=============================================================================*/

void PrRStat(FILE *txt, Histo hist, double mean, double std)
{
  fprintf(txt, "Scores statistics from random sampling:\n");
  fprintf(txt, "Rmin:   %.2f \n", hist.hmin);
  fprintf(txt, "Rmax:   %.2f \n", hist.hmax);
  fprintf(txt, "Rmean:  %.2f \n", mean);
  fprintf(txt, "Rstd:   %.2f \n\n", std);
  return;
}
/*=============================================================================
 PrHist(): Enregistre dans le fichier identifie par 'fname' les donnees, abs-
           cisses et ordonnees de l'histogramme 'hist' dont l'integrale est
           supposee etre egal a 1 (pas de controle).
           Une gaussienne de meme moyenne et ecart type que 'hist' est enregis-
           tree en 3eme colonne pour comparaison.
=============================================================================*/

void PrHist(char *fname, Histo hist)
{
  FILE   *txt;
  Map     map;
  double  mean, std;
  int     i;

  GetHistStat(hist, &mean, &std);

  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins-1);

  txt = fopen(fname, "w");

  for (i = 0; i < hist.bins; i++) {
      double x, y;
      x = Getx(map, i);
      y = Gauss(x, mean, std);
      fprintf(txt, "%.4e %.4e %.4e\n", x, hist.vals[i], y);
  }
  fclose(txt);

  return;
}
/*=============================================================================
 PrEval(): Enregistre dans le fichier identifie par 'fname' les elements du
           tableau contenu dans 'hist' dont les abscisses sont limitees par
           'min' et 'max' et les abscisses correspondantes.
=============================================================================*/
#define SMALL    (1.e-30)

void PrEval(char *fname, Histo hist, double min, double max)
{
  FILE   *txt;
  Point  *Pts;
  Map    map;
  int    i, X1, X2, len;

  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins - 1);      /* etalonnage */

  X1 = GetX(map, min);
  X2 = GetX(map, max);

  if (X1 >= X2 || X1 < 0 || X2 > hist.bins - 1) {
      fprintf(stderr, "PrEval: invalid arguments 3 and 4, exit..\n");
      exit(1);
  }

  len = X2 - X1 + 1;
  Pts = (Point *) malloc(len * sizeof(Point));

  for (i = 0; i < len; i++)
  {
      double u = hist.vals[X1 + i];
      if (u < SMALL) break;            /* eviter ulterieurement les log de 0 */

      Pts[i].x = Getx(map, X1 + i);
      Pts[i].y = u;
  }
  len = i;                               /* reduit la longueur si necessaire */

  txt = fopen(fname, "w");

  for (i = 0; i < len; i++)
      fprintf(txt, "%.3e %.3e\n", Pts[i].x, Pts[i].y);

  fclose(txt);
  free(Pts);

  return;
}
#undef SMALL
/*=============================================================================
 PrLogEval(): Enregistre dans le fichier identifie par 'fname' les elements du
           logarithme du tableau contenu dans 'hist' dont les abscisses sont
           limitees par 'min' et 'max' et les abscisses correspondantes.
           La droite ajustee a ces donnees par la methode des moindres carres
           est aussi enregistree, ses parametres 'A' et 'B' (B: pente) sont
           pointes en retour.
           Retourne 1 si le fit lineaire a pu etre effectue sinon 0.
=============================================================================*/
#define SMALL    (1.e-30)
#define MIN_LEN  2

int PrLogEval(char *fname, Histo hist, double min, double max,
              double *A, double *B)
{
  FILE   *txt;
  Point  *Pts;
  Map    map;
  int    i, X1, X2, len;

  map = MapSetup(hist.hmin, hist.hmax, 0, hist.bins - 1);      /* etalonnage */

  X1 = GetX(map, min);
  X2 = GetX(map, max);

  if (X1 >= X2 || X1 < 0 || X2 > hist.bins - 1) {
      fprintf(stderr, "PrEval: invalid arguments 3 and 4, exit..\n");
      exit(1);
  }

  len = X2 - X1 + 1;

  Pts = (Point *) malloc(len * sizeof(Point));

  for (i = 0; i < len; i++)
  {
      double u = hist.vals[X1 + i];
      if (u < SMALL) break;               /* sortie pour eviter les log de 0 */

      Pts[i].x = Getx(map, X1 + i);
      Pts[i].y = log(u);
  }
  len = i;                               /* reduit la longueur si necessaire */

  if (len < MIN_LEN) {
      free(Pts);
      return 0;                               /* echec de l'operation de fit */
  }

  lfit(Pts, len, A, B);                 /* coefficients de la droite A + B*x */

  txt = fopen(fname, "w");

  for (i = 0; i < len; i++)
  {
      double  X, Y, Z;

      X = Pts[i].x;
      Y = Pts[i].y;
      Z = *A + *B * X;
      fprintf(txt, "%.3e %.3e %.3e\n", X, Y, Z);
  }
  fclose(txt);

  free(Pts);
  return 1;
}
#undef SMALL
#undef MIN_LEN
/*=============================================================================
 PrMaskTStat(): Procede au calcul des scores d'un masque et a leur statisti-
                que, les resultats sont diriges sur le fichier 'txt'.
                Les seuils correspndant a divers pourcentages de captures sont
                affiches ainsi que les moyenne, ecart type, maximum et minimum
                qui sont aussi pointes en sortie par les derniers arguments.
=============================================================================*/

void PrMaskTStat(FILE *txt, Trset *trset, Mask *mask,
                 double *min, double *max , double *mean , double *std)
{
  double *scores;
  char   *str = (char *) malloc(80);
  int    i;

  scores = GetMaskTScores(trset, mask);

  sprintf(str, "{ ");
  for (i = 0; i < mask->nargs; i++) {
      sprintf(str, "%s%d ", str, mask->args[i]);
  }
  sprintf(str, "%s%c", str, '}');

  fprintf(txt, "\n");
  fprintf(txt, "===========================================================\n");
  fprintf(txt, " REGION %s STATISTICS:\n", str);
  fprintf(txt, "===========================================================\n");

  fPrintfCaptures(scores, trset->nseq, txt);
  GetTScoresStat(scores, trset->nseq, min, max, mean, std);

  PrTStat(txt, *min, *max, *mean, *std);

  free(scores);
  free(str);

  return;
}
/*=============================================================================
 MstatHelp():
=============================================================================*/

void MstatHelp(void)
{
  fprintf(stderr,

"mstat: Mask STATistics\n\n"
"Usage:\n"
"mstat [-h]\n"
"mstat <training-set> <region> <mask>\n"
"      [-bkg <%%A><%%T><%%G>[<%%C>]]        background percents, default: 25,25..\n"
"      [-logzero <logzero>]                                      default: -30\n"
"      [-dh <dh>]           elementary length of histogram axis, default: 0.2\n"
"      [-spc1 <spc1>]            random sampling: number for the first column\n"
"                                                   of profile, default:  300\n"
"      [-cutoff <cutoff>]               score cutoff, default: Tmean - 2*Tstd\n"
"      [-Mbs <Mbs>]             number (in Mb) of involved sites, default: 10\n"
"      [-hist]         stores histogram of finite scores in 'mhisto.dat' file\n"
"      [-evals]              stores 'E-values vs scores' in 'mevals.dat' file\n"
"      [-sumf <fname>]           substitution matrix file name, default: none\n"
"      [-pcw <pcw>]                        pseudo-counts weight, default: 0.1\n"
"      [-hpcw <hpcw>]                        pseudo-counts weight for helices\n"
"      [-spcw <spcw>]                        pseudo-counts weight for strands\n\n"

);

  return;
}
/*===========================================================================*/
