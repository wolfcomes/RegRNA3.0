
/*=============================================================================
 erpin.c                            A.Lambert le 23/12/01  revu le 12/01/05

 - recherche de masques a plusieurs niveaux, dans l'ordre des arg. de 'main'.
 - exploration complete de toutes les sequences.
 - exploration "transparente" des bases d'entrainement.
 - les seuils sont entres en valeurs ou en % de captures dans 'trset'.
 - on dispose de divers styles d'output (LONG | SHORT | MUTE).
 - on traite le brin complementaire, en plus ou separement.
 - on utilise la statistique des sequences analysees (poids et transitions).
 - on opere des selections sur les sequences a analyser (position, longueur).
 - on calcule la E-value des detections au dernier niveau,
 - ainsi que l'histogramme des detections a ce meme niveau, option '-hist'.
 - on gere les parametres 'logzero', 'tablen', 'chrono', 'warnings'.
 - les repetitions et overlap sont elimines avant affichage.
 - les masques successifs sont geres de facon statique ou dynamique.
 - un controle des besoins en memoire est prealablement effectue.
 - des matrices de substitution peuvent etre introduites, lesquelles modifient
   les profils des brins et/ou helices.

 make -f apps.mk erpin ;
 mv erpin ../bin ;

 cc -O3 -Wall -o erpin erpin.c Eval2.o -I../include -L../lib -lrnaIV -lm ;
 chmod 755 erpin ;
 strip erpin ;

 erpin <trset> <data> <region>
       -nomask|((-mask|-umask|-add) <elt1> <elt2>..)      level1,
       [-nomask|((-mask|-umask|-add) <elt1> <elt2>..)]    level2, default: void
       [...]                                              level..          idem
       [-cutoff <cutoff1> <cutoff2> ..]                           default: 100%
       [-dmp|-smp]                                                default: -dmp
       [-fwd|-rev|-fwd+rev]                                   default: -fwd+rev
       [-long|-short|-mute]                                      default: -long
       [-warnings]                                                 default: OFF
       [-globstat|-locstat|-unifstat]                        default: -globstat
       [-Eon|-Eoff]                                               default: -Eon
       [-hist]                                                     default: OFF
       [-seq1 <seqnb1>][-nseq <nseq>]               defaut: 1, SEQ_MAX_NB (all)
       [-bgn <seqbgn>][-len <range>]                     defaut: 1, SEQ_MAX_LEN
       [-logzero <logzero>]                                        default: -20
       [-tablen <tablen>]                                         default: 1024
       [-chrono]                                                   default: OFF

       [-sumf <fname>]           fichier contenant les matrices de substitution
       [-pcw <pcw>]                       poids des pseudo-comptes, defaut: 0.1
       [-hpcw <hpcw>]                 poids des pseudo-comptes pour les helices
       [-spcw <spcw>]                   poids des pseudo-comptes pour les brins

 default command line: "erpin5 <trset> <data> <region> <mask>" stands for:

 erpin <trset> <data> <region> <mask> -cutoff 100% -fwd+rev -long -globstat \
       -seq1 1 -nseq 10M -bgn 1 -len 300M -logzero -20. -tablen 1024 -dmp ;


 erpin ~/devc/projets/bio/data/trsetsII/trna-typeI.epn \
       ~/devc/projets/bio/data/sequences/ecoli.fna -2,2 \
       -cutoff 100% 100% 90% \
       -umask 6 8 -mask 2 3 -nomask -chrono -short -logzero -30 ;


 erpin ~/devc/projets/bio/data/trsetsII/trna-typeI.epn \
       ~/devc/projets/bio/data/sequences/Chr_22_01-12-1999.fasta -2,2 \
       -cutoff 100% 100% 90% \
       -umask 6 8 -mask 2 3 -nomask -chrono -short -logzero -30 ;


 erpin ~/devc/projets/bio/data/trsetsII/trna-typeI.epn \
       ~/devc/projets/bio/data/sequences/ecoli.fna -2,2 \
       -cutoff 100% 100% 90% \
       -umask 6 8 -mask 2 3 -nomask -fwd -chrono -short -logzero -30 \
       -sumf ../sum/SUM.dat ;

 erpin ~/devc/projets/bio/data/trsetsII/trna.db \
       ~/devc/projets/bio/data/sequences/ecoli.fna -2,2 \
       -cutoff 100% 100% 100% \
       -umask 7 20 -mask 2 3 -nomask -fwd -chrono -smp -short -logzero -30 ;

 erpin ~/devc/projets/bio/data/trsetsIII/LET-7.epn \
       ~/devc/projets/bio/data/sequences/randseq10M.fna -2,2 \
       -umask 2 4 6 8 3 -cutoff -39 -sumf ../sum/SUM.dat -fwd \
       -chrono -mute -logzero -30 ;

E-value at cutoff -39.0 for 10.0Mb single strand data: 46.2e+01
37 hits

 erpin ~/devc/projets/bio/data/trsetsIII/LET-7.epn \
       ~/devc/projets/bio/data/sequences/randseq10M.fna -6,6 \
       -umask 6 -sumf ../sum/SUM.dat -fwd -chrono -mute -logzero -30 ;

E-value at cutoff 13.1 for 10.0Mb single strand data: 2.93e+01
21 hits

 erpin ~/devc/projets/bio/data/trsetsIII/LET-7.epn \
       ~/devc/projets/bio/data/sequences/randseq10M.fna -6,6 -umask 6 \
       -cutoff 5 -sumf ../sum/SUM.dat -fwd -chrono -mute -logzero -30 ;

E-value at cutoff 5.0 for 10.0Mb single strand data: 2.28e+02
230 hits
=============================================================================*/

#include "rnaIV.h"

void ErpinHelp(void);

/*===========================================================================*/

int main(int argc, char *argv[])
{
  FILE      *txt;
  Sequence  *seq;
  Threshold *thresholds;
  Trset     *trset;
  Pattern   *pattern;
  Mask      *mask;
  Context   *ctxt;
  int       nmask, style, dir, stat, warnings, chrono, maskproc,
            first, nseq, bgn, range, eval, hist,
            total = 0;
  double    datavol;

  int       pcflag;               /* variables concernant les pseudo-comptes */
  double    hpcw, spcw;
  char      *sumfname;

  /*--------------- lecture des arguments -----------------------------------*/

  ReadHelpArgs(argc, argv, 4, ErpinHelp);
  fTest(argv[1]);
  fTest(argv[2]);

  mask = ReadMasksArgs(argc, argv, &nmask);

  thresholds = GetThresholdsArgs(argc, argv, nmask);
  GetEnvArgs(argc, argv, &chrono);
  GetProcArg(argc, argv, &maskproc);
  GetOutputArgs(argc, argv, &style, &warnings, &eval, &hist);
  GetDirArg(argc, argv, &dir);
  GetStatArg(argc, argv, &stat);
  GetSeqArgs(argc, argv, &first, &nseq, &bgn, &range);
  sumfname = GetSumArgs(argc, argv, &pcflag, &hpcw, &spcw);

  /*----------------- Trset, Pattern & Masques, Profils  --------------------*/

  trset = ReadTrset(argv[1], 'S', stdout);                          /* Trset */
  pattern = ReadPattern(argv[3], trset);                          /* Pattern */
  ParseMasksArgs(mask, nmask, pattern);
  InitStatTables();                             /* distr. uniforme du "fond" */

  SetNtCodes();        /* tableau de codage des nt pour le calcul des scores */

  if (pcflag) {                                            /* pseudo-comptes */
      ReadSUM(sumfname, hpcw, spcw, trset);
      GetPatternProfilesSM(trset, pattern);                      /* Profiles */
  }
  else GetPatternProfiles(trset, pattern);                       /* Profiles */

  GetMasks(mask, nmask, pattern);                                   /* Masks */

                                                    /* Ctrl du nb de configs */
  CtrlCfgs(mask, nmask, maskproc, 'S', stdout);        /* avant le demarrage */
                                                   /* du calcul des tableaux */

  if (maskproc == DYNAMIC) GetMaskCfgs(mask);       /* Configs du 1er masque */
  else GetMasksCfgs(mask, nmask);                         /* Configs de tous */

  /*--------------- Seuils supposant une database uniforme ------------------*/

  GetThresholds(thresholds, trset, mask, nmask);               /* Thresholds */
  fPrintfThresholds(mask, nmask, stdout);            /* affichage des seuils */
  DelTrset(trset);

  /*------------------- Statistique sur la Database -------------------------*/

  fprintf(stdout, "Database:\t\"%s\"\n", argv[2]);

  if (stat == GLOBAL) {
      datavol = fGetLongStat(argv[2], first, nseq, bgn, range,
                             mask[nmask-1].min_len, 'V', stdout);
      ChPatternStat(pattern);
  }
  else
      datavol = fGetShortStat(argv[2], first, nseq, bgn, range,
                              mask[nmask-1].min_len, 'V', stdout);

  /*------------- Calcul de la E-value du masque principal ------------------*/

  if (eval == ON) {
      double E, cutoff = mask[nmask-1].threshold;
      GetEvalues(mask + nmask - 1, datavol, dir, stat);
      E = Evalue(cutoff);
      fPrintfEvalue(stdout, E, cutoff, dir, datavol);
  }

  /*--------------- Sequences et Contextes ----------------------------------*/

  seq = InitSeq(GetTail(pattern));                              /*  Sequence */
  ctxt = SetupContext(mask, nmask);                             /* Contextes */
  InitSeqContexts(ctxt, nmask, seq);
  StoreInterval(ctxt, bgn, range);
  InitOutputContexts(ctxt, nmask, style, hist, eval, warnings, chrono, stdout);
  SetMasksScoresTabs(ctxt, nmask);                             /* Score Tabs */

  /*--------------- Exploration des Sequences -------------------------------*/


  StartStatus();           /* debut de l'affichage du nbre de sites explores */
  StartTimer(ctxt);                                   /* et du chronometrage */

  txt = fopen(argv[2], "r");

  while ( ReadSeq(txt, seq) != 0 )                                 /* Search */
  {
      if (seq->nb >= first && seq->nb < first + nseq)
      {
          PreProcessSeq(seq);

	  if (stat == UNIFORM) {
              if (dir == FORWARD || dir == FWD_AND_REV) {
                  ResetDetect(ctxt, nmask);
                  SetDir(ctxt, nmask, FORWARD);
                  total += Search(ctxt, nmask, maskproc);
              }
              if (dir == REV_CMPL || dir == FWD_AND_REV) {
                  ResetDetect(ctxt, nmask);
                  SetDir(ctxt, nmask, REV_CMPL);
                  RevCmpl(seq);
                  total += Search(ctxt, nmask, maskproc);
              }
	  }
          else                      /* stat est alors egal a LOCAL ou GLOBAL */
	  {
              if (stat == LOCAL) {
                  sGetStat(seq->data, seq->datalen);
                  ChPatternStat(pattern);
              }
              if (dir == FORWARD || dir == FWD_AND_REV) {
                  ResetDetect(ctxt, nmask);
                  SetDir(ctxt, nmask, FORWARD);
                  total += Search(ctxt, nmask, maskproc);
              }
              if (dir == REV_CMPL || dir == FWD_AND_REV) {
                  ResetDetect(ctxt, nmask);
                  SetDir(ctxt, nmask, REV_CMPL);
                  RevCmpl(seq);
                  ReverseStat();
                  ChPatternStat(pattern);
                  total += Search(ctxt, nmask, maskproc);

                  if (stat == GLOBAL) {        /* restaure les freq. du fond */
                      ReverseStat();
                      ChPatternStat(pattern);
                  }
	      }
          }
      }
  }
  fclose(txt);
  fPrintfProcessLog(ctxt, nmask, total);

  /*-------------------- Sortie ---------------------------------------------*/

  if (eval == ON) FreeEvalsTab();
  DelMasks(mask, nmask);
  DelPattern(pattern);
  DelSeq(seq);
  free(ctxt);
  free(thresholds);
  FreeNtCodes();
  exit(0);
}
/*=============================================================================
 ErpinHelp(): help
=============================================================================*/

void ErpinHelp(void)
{
  fprintf(stderr,
"\n"
"=========== erpin: Easy Rna Profile IdentificatioN, Version 4.2.5 ==========\n\n"

"Usage:\n"
"erpin [-h]                                                              help\n"
"erpin <training-set>                                  training set file name\n"
"      <input-file>                                database file name (fasta)\n"
"      <region>                                            region of interest\n"
"      -nomask|((-mask|-umask|-add) <elt1> ...)      level1,\n"
"      [-nomask|((-mask|-umask|-add) <elt1> ...)]    level2,    default: void\n"
"      [...]                                         level..             idem\n"
"      [-cutoff <cutoff1> <cutoff2> ..]                         default: 100%%\n"
"      [-dmp|-smp]                                              default: -dmp\n"
"      [-fwd|-rev|-fwd+rev]                                 default: -fwd+rev\n"
"      [-long|-short|-mute]                                    default: -long\n"
"      [-warnings]                                               default: OFF\n"
"      [-globstat|-locstat|-unifstat]                      default: -globstat\n"
"      [-Eon|-Eoff]                             E-value ON|OFF, default: -Eon\n"
"      [-hist]                                                   default: OFF\n"
"      [-seq1 <seqnb1>][-nseq <nseq>]            default: 1, SEQ_MAX_NB (all)\n"
"      [-bgn <seqbgn>][-len <range>]                  default: 1, SEQ_MAX_LEN\n"
"      [-logzero <logzero>]                                      default: -20\n"
"      [-tablen <tablen>]                                       default: 1024\n"
"      [-chrono]                                                 default: OFF\n"
"      [-sumf <fname>]           substitution matrix file name, default: none\n"
"      [-pcw <pcw>]                      pseudo-counts weight, default:   0.1\n"
"      [-hpcw <hpcw>]                        pseudo-counts weight for helices\n"
"      [-spcw <spcw>]                        pseudo-counts weight for strands\n\n"

"default command line: \"erpin <trset> <data> <region> <mask>\" stands for:\n\n"

"erpin <trset> <data> <region> <mask> -cutoff 100%% \\ \n"
"      -dmp -fwd+rev -long -globstat -Eon \\ \n"
"      -seq1 1 -nseq 10M -bgn 1 -len 300M -logzero -20. -tablen 1024\n\n"

"Example:\n"
"erpin trsets/trna.epn  sequences/ecoli.fna -2,2 -umask 20 11 -nomask \\ \n"
"      -cutoff 100%% 90%% -sumf ../sum/SUM.dat -fwd -short -chrono\n\n"

"============ D.Gautheret A.Lambert, Luminy Marseille, Jan 2005 =============\n\n"
);

  return;
}
/*===========================================================================*/
