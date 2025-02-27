
/*=============================================================================
 ev.c                                         A.Lambert le 19/04/04

 Ce programme calcule la E-Value du motif ultime parmi ceux entres en arguments
 comme le ferait une ligne de commande equivalente passee a 'erpin'.

 Le volume d'une hypothetique sequence aleatoire ainsi que sa composition peu-
 vent etre passes au programme soit explicitement (option '-bkg'), soit par re-
 ference a une base de donnees (option '-data').
 Les profils construits depuis la base d'entrainement peuvent integrer les don-
 nees de matrices de substitution.

 cc -O3 -Wall -o ev ev.c -I../include -L../lib -lrnaIV -lm ;
 chmod 755 ev ;
 strip ev ;
 mv ev ../bin ;

 ev <trset> <region>
    -nomask|((-mask|-umask|-add) <elt1> <elt2>..)          level1,
    [-nomask|((-mask|-umask|-add) <elt1> <elt2>..)]        level2, defaut: void
    [...]                                                  level..         idem
    [-cutoff <cutoff1> <cutoff2> ..]                               defaut: 100%
    [-logzero <logzero>]               borne inf. imposee aux log., defaut: -20

    [-bkg <%A><%T><%G>[<%C>]]    pourcentages des bases du 'fond', defaut: 25..
    [-data <database>]                mesure des frequences depuis une database
    [-Mbs <Mbs>]              nbre suppose (en Mb) de sites visites, defaut: 10
                                        ou le vol. du fichier associe a '-data'
    [-sumf <fname>]              fichier contenant les matrices de substitution
    [-pcw <pcw>]                          poids des pseudo-comptes, defaut: 0.1
    [-hpcw <hpcw>]                    poids des pseudo-comptes pour les helices
    [-spcw <spcw>]                      poids des pseudo-comptes pour les brins

 ev ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -2,2 \
    -cutoff 100% 100% 90% \
    -umask 6 8 -mask 2 3 -nomask -logzero -30 \
    -data ~/devc/projets/bio/data/sequences/ecoli.fna ;
E-value at cutoff 60.3 for 4.6Mb single strand data = 3.18e-14

 ev ~/devc/projets/bio/data/trsetsIII/LET-7.epn -2,2 \
    -umask 2 4 6 8 3 -cutoff -39 -sumf ../sum/SUM.dat  -Mbs 10 -logzero -30 ;
E-value at cutoff -39.0 for 10.0Mb single strand data = 4.63e+01

 ev ~/devc/projets/bio/data/trsetsIII/LET-7.epn -2,2 \
    -umask 2 4 6 8 3 -cutoff -30 -sumf ../sum/SUM.dat -Mbs 10 -logzero -30 ;
E-value at cutoff -30.0 for 10.0Mb single strand data = 6.38e+00
 erpin ~/devc/projets/bio/data/trsetsIII/LET-7.epn \
       ~/devc/projets/bio/data/sequences/randseq10M.fna -2,2 \
       -umask 2 4 6 8 3 -cutoff -30 -sumf ../sum/SUM.dat -Mbs 10 -logzero -30 \
       -fwd -mute ;
22 config. per site
6 hits

=============================================================================*/

#include "rnaIV.h"

int ReadStatArgs(int argc, char *argv[], double *freqs, double *Mbs);

/*===========================================================================*/

int main(int argc, char *argv[])
{
  double freqs[4] = {0.25, 0.25, 0.25, 0.25};    /* freqs supposees des ATGC */
  double Mbs = 10.;                      /* volume suppose des donnees en Mb */

  Trset     *trset;
  Pattern   *pattern;
  Mask      *mask;
  Threshold *thresholds;
  double    E, datavol;
  int       i, nmask;

  int       pcflag;               /* variables concernant les pseudo-comptes */
  double    hpcw, spcw;
  char      *sumfname;


  if (argc < 4) {
      fprintf(stderr, "'%s' needs at least 3 arguments, exit..\n", argv[0]);
      exit(1);
  }
  fTest(argv[1]);
  ChLogZero(-20.0);                                     /* valeur par defaut */
  InitStatTables();
  if (ReadStatArgs(argc, argv, freqs, &Mbs))  ReInitStatTables(freqs);

  for (i = 3; i < argc; i++) {
      if (strcmp(argv[i], "-logzero") == 0)  ChLogZero(atof(argv[++i]));
      if (strcmp(argv[i], "-Mbs") == 0)      Mbs = atof(argv[++i]);
  }

  fprintf(stdout, "\nBkg ATGC ratios: %.3f %.3f %.3f %.3f\n",
                   freqs[0], freqs[1], freqs[2], freqs[3]);

  sumfname = GetSumArgs(argc, argv, &pcflag, &hpcw, &spcw);
  fprintf(stdout, "H&S pseudo-count weights: %.2f %.2f\n", hpcw, spcw);

  mask = ReadMasksArgs(argc, argv, &nmask);
  thresholds = GetThresholdsArgs(argc, argv, nmask);
  trset = ReadTrset(argv[1], 'S', stdout);
  pattern = ReadPattern(argv[2], trset);
  ParseMasksArgs(mask, nmask, pattern);
  SetNtCodes();        /* tableau de codage des nt pour le calcul des scores */

  if (pcflag) {                                            /* pseudo-comptes */
      ReadSUM(sumfname, hpcw, spcw, trset);
      GetPatternProfilesSM(trset, pattern);                      /* Profiles */
  }
  else GetPatternProfiles(trset, pattern);                       /* Profiles */

  GetMasks(mask, nmask, pattern);                                   /* Masks */
  GetThresholds(thresholds, trset, mask, nmask);               /* Thresholds */
  fPrintfThresholds(mask, nmask, stdout);            /* affichage des seuils */

  datavol = 1.e+6 * Mbs;
  GetEvalues(mask + nmask - 1, datavol, FORWARD, GLOBAL);
  E = Evalue(mask[nmask-1].threshold);

  fprintf(stdout, "================================================================\n");
  fprintf(stdout, "E-value at cutoff %.1f for %.1fMb single strand data: %.2e\n",
                  mask[nmask-1].threshold, Mbs, E);
  fprintf(stdout, "================================================================\n");


  FreeEvalsTab();
  free(thresholds);
  DelMasks(mask, nmask);
  DelPattern(pattern);
  return 0;
}
/*=============================================================================
 ReadStatArgs(): Lit les arguments qui fixent les frequences des ATGC du 'fond'.
                 Les 2 options '-bkg' (entree explicite) et '-data' (entree des
	         frequences mesurees dans un fichier de donnees) s'excluent mu-
		 tuellement.
		 Sous l'option '-data' l'argument 'Mbs' pointe en retour le vo-
		 lume des donnees (en Megabases).
		 La fonction retourne la valeur 1 si l'une des 2 options a ete
		 utilisee, sinon 0.
=============================================================================*/

int ReadStatArgs(int argc, char *argv[], double *freqs, double *Mbs)
{
  int i;

  for (i = 1; i < argc-1; i++)
  {
      if (strcmp(argv[i], "-bkg") == 0)  {
          ReadBkgFreqs(freqs, i, argv);  /* pourcentages des bases du "fond" */
	  return 1;
      }
      if (strcmp(argv[i], "-data") == 0)  {
	  *Mbs = 1.e-6 * fGetFreqs(argv[++i], freqs);
	  return 1;
      }
  }
  return 0;
}
/*===========================================================================*/
