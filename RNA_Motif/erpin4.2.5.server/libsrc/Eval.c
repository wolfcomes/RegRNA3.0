
/*=============================================================================
 Eval.c                          A.Lambert le 04/05/03   revu le 30/04/04

 Code des fonctions qui permettent le calcul de la "E-value" des scores d'un
 'masque' sur une sequence aleatoire par produit de convolution des colonnes
 des profils d'helices.

 cc -O2 -Wall -c Eval.c -I../include ;

 ar -rs ../lib/librnaIV.a Eval.o ;

=============================================================================*/

#include "rnaIV.h"

Histo  GetEvals1(Mask *mask, double datavol, int smplc1);
Histo  GetEvals(Mask *mask, double datavol, int dir, int stat);
void   GetEvalues(Mask *mask, double datavol, int dir, int stat);
double Evalue(double x);
void   FreeEvalsTab(void);
double evprob(double p, unsigned int N);
void  fPrintfEvalue(FILE *fout, double E, double cutoff, int dir, double datavol);

/*=============================================================================
 GetEvals1(): Etant donne le masque pointe par 'mask' suppose correctement ini-
              tialise par 'GetMask', cette fonction retourne le tableau permet-
              tant de determiner le nombre moyen de scores du masque superieurs
              a un seuil donne sur le volume de donnees 'datavol' (exprime en
              nucleotides).
              'smplc1' designe le nombre d'echantillons examines pour la 1ere
              colonne des profils de brins.
=============================================================================*/

Histo GetEvals1(Mask *mask, double datavol, int smplc1)
{
  Histo   hist, cdf;
  double  Prfsc, Prob, cfgs;
  int     i;

  cfgs = rint(pow(2.0, mask->log2ncfg));         /* nombre de configurations */

  hist = GetMaskHist(mask, DELTA_H, &Prfsc, smplc1);
  cdf = GetHistCdf(hist);

  for (i = 0; i < cdf.bins; i++) {
      Prob = Prfsc * cdf.vals[i];
      cdf.vals[i] = datavol * evprob(Prob, (unsigned int) cfgs);
  }
  free(hist.vals);

  return cdf;
}
/*=============================================================================
 GetEvals(): Interface de 'GetEvals1' destinee au traitement des simples et
             doubles brins. Les frequences des bases du "fond" sont actualisees
             avant le traitement du brin complementaire.
             L'argument 'dir' peut prendre les valeurs 'FORWARD' 'REV_CMPL' ou
             FWD_AND_REV' (soit: 0, 1 ou 2).
	     L'argument 'stat' peut prendre les valeurs 'UNIFORM' 'LOCAL' ou
	     'GLOBAL'.
             La fonction suppose que les variables globales 'LogDataFreqs' ont
             ete initialisees pour le brin direct (pas de controle): c'est le
             cas dans 'erpin' ou cette fonction est appelee avant le demarrage
             des operations de recherche.
             Avant la sortie la fonction remet les 'LogDataFreqs' et profils a
             leurs valeurs trouvees en entree.
             Note:
             La E-value correspondant a un cutoff quelconque x sera obtenue
             simplement par:
             E = Interpol(Evals, x);
             rappel: La fonction 'Interpol' (cdf.c) prolonge des 2 cotes l'his-
             togramme passe en argument par ses valeurs extremes
=============================================================================*/

Histo GetEvals(Mask *mask, double datavol, int dir, int stat)
{
  Histo Evals, Evals2;
  Map   map;
  int   i;

  if (dir != FORWARD && dir != REV_CMPL && dir != FWD_AND_REV) {
      fprintf(stderr, "GetEvals: unknown argument 3, exit..\n");
      exit(1);
  }
  if (stat != UNIFORM && stat != LOCAL && stat != GLOBAL) {
      fprintf(stderr, "GetEvals: unknown argument 4, exit..\n");
      exit(1);
  }

  if (stat == UNIFORM)
  {
      Evals = GetEvals1(mask, datavol, SAMPLES1_SINGLE);
      if (dir == FWD_AND_REV)                               /* deux passages */
          for (i = 0; i < Evals.bins; i++) Evals.vals[i] *= 2.0;

      return Evals;
  }
                                        /* 'stat' est egal a LOCAL ou GLOBAL */

  if (dir == FORWARD) {
      Evals = GetEvals1(mask, datavol, SAMPLES1_SINGLE);
  }
  else
  if (dir == REV_CMPL) {
      ReverseStat();        /* passage aux frequences du brin complementaire */
      ChMaskStat(mask);                         /* actualisation des profils */
      Evals = GetEvals1(mask, datavol, SAMPLES1_SINGLE);
      ReverseStat();          /* retour au brin direct (involution: I^2 = 1) */
      ChMaskStat(mask);                       /* retour aux profils initiaux */
  }
  else {                                           /* direction: FWD_AND_REV */
      Evals = GetEvals1(mask, datavol, SAMPLES1_DOUBLE);

      ReverseStat();        /* passage aux frequences du brin complementaire */
      ChMaskStat(mask);                         /* actualisation des profils */
      Evals2 = GetEvals1(mask, datavol, SAMPLES1_DOUBLE);
      ReverseStat();          /* retour au brin direct (involution: I^2 = 1) */
      ChMaskStat(mask);                       /* retour aux profils initiaux */

      map = MapSetup(Evals.hmin, Evals.hmax, 0, Evals.bins - 1);

      for (i = 0; i < Evals.bins; i++) {
          double x = Getx(map, i);
          Evals.vals[i] += Interpol(Evals2, x);                 /* voir Note */
      }
      free(Evals2.vals);
  }

  return Evals;
}

/*=============================================================================
 GetEvalues(): Interface de 'GetEvals' destinee a manipuler silencieusement
               l'histogramme 'MainMaskEvals' declare en variable globale.
	       Version utilisant l'iteration des convolutions.
=============================================================================*/

void GetEvalues(Mask *mask, double datavol, int dir, int stat)
{
  extern Histo MainMaskEvals;

  MainMaskEvals.vals = NULL;            /* initialisation du pointeur 'vals' */
  MainMaskEvals = GetEvals(mask, datavol, dir, stat);
  return;
}
/*=============================================================================
 Evalue(): Retourne pour un cutoff fixe a la valeur 'x', la E-value calculee a
           l'aide de l'histogramme 'MainMaskEvals' declare en variable globale.
=============================================================================*/

double Evalue(double x)
{
  extern Histo MainMaskEvals;
  double E = Interpol(MainMaskEvals, x);
  return E;
}
/*=============================================================================
 FreeEvalsTab(): Libere la memoire allouee a l'histogramme 'MainMaskEvals'
                 declare en variable globale.
=============================================================================*/

void FreeEvalsTab(void)
{
  extern Histo MainMaskEvals;
  if (MainMaskEvals.vals != NULL) free(MainMaskEvals.vals);
  return;
}
/*=============================================================================
 evprob(): "extreme value probability"
           Soient 'N' variables aleatoires 'Xi', i = 1,2..N, independantes et
           identiquement distribuees, p = P(Xi > x) et M = max(Xi), cette fonc-
           tion retourne la probabilite P(M > x) = 1 - (1 - p)^N.

	   p doit etre inferieur ou egal a 1. Si, par suite d'une erreur d'ar-
	   rondi, p > 1, alors il est remis a 1, la fonction retournant 1.
	   N est un entier, si N vaut 1 la fonction retourne p.

           Si N*p >= 1 'evprob' utilise la fonction 'pow' de la bibliotheque
           standard, sinon elle utilise le developpement du binome a l'ordre
           'MAX_ORDER' en p.
           Avec q = -p, P(Xi > x) est donne par:
           A = - { Nq + N(N-1)q^2/2 + N(N-1)(N-2)q^3/2.3 + ... + q^N }
=============================================================================*/
#define MAX_ORDER 30

double evprob(double p, unsigned int N)
{
  unsigned int i;
  double  A, b;

  if (p >= 1.) return 1. ;
  if (N == 1)  return p ;

  if (p*N >= 1.)  return  (1. - pow(1. - p, (double) N)) ;

  p = -p;
  A = 0.0;
  b = 1.0;
  for (i = 1; i <= N  && i <= MAX_ORDER ; i++)
  {
      b *= (double)(N-i+1) * p / (double)i;
      A += b;
  }

  return (-A) ;
}
#undef MAX_ORDER

/*=============================================================================
 fPrintfEvalue(): Affiche sur 'fout' une E-value en precisant, suivant la va-
                  leur de l'argument 'dir': FORWARD ou REV_CMPL ou FWD_AND_REV,
                  si elle a ete obtenue sur un simple ou double brin.
                  'datavol' est le nombre de nucleotides visitees (1 ou 2 fois)
=============================================================================*/

void fPrintfEvalue(FILE *fout, double E, double cutoff, int dir, double datavol)
{
  char   str[12];

  if (datavol < 1000)   sprintf(str, "%db", (int) datavol);
  else
  if (datavol < 1.e+5)  sprintf(str, "%.1fKb", datavol*1.e-3);
  else
  if (datavol < 1.e+8)  sprintf(str, "%.1fMb", datavol*1.e-6);
  else
      sprintf(str, "%.1fGb", datavol*1.e-9);

  fprintf(fout, "E-value at cutoff %.1f for %s %s strand data: %.2e\n\n",
          cutoff, str, dir == FWD_AND_REV ? "double" : "single", E);

  return;
}
/*===========================================================================*/
