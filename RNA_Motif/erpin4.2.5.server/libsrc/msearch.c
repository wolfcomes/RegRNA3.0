
/*=============================================================================
 msearch.c                         A.Lambert le 16/12/01 revu le 22/05/02

 fonctions impliquees dans une operation de recherche de masques, dans une
 sequence.
 L'operation se fait a 1 ou plusieurs niveaux emboites, le 2eme, optionnel, su-
 bordonne aux detections du 1er, etc ...

 Afin de conserver presents les parametres fixant les conditions de l'operation
 et son contexte a chaque niveau, on cree 1 nouvelle structure:

 - Context: parametres attributs et contraintes lies a l'operation.

 Il ne sera pas alloue de memoire aux pointeurs de la structure 'Context',
 ceux-ci pointeront des structures deja creees, il n'est donc pas prevu de
 fonction particuliere de liberation de memoire ('free' suffira).

 cc -O2 -Wall -c msearch.c -I../include -DDEBUG;

 cc -O2 -Wall -c msearch.c -I../include ;
 ar -rs ../lib/librnaIV.a msearch.o ;

=============================================================================*/

#include "rnaIV.h"

extern int ScoreTabLen;                   /* longueur des tableaux de scores */

Context *NewContext(void);
Context *SetupContext(Mask *mask, int nctxt);
void InitSeqContext(Context *ctxt, Sequence *seq);
void InitSeqContexts(Context *ctxt, int nctxt, Sequence *seq);
void InitOutputContext(Context *ctxt, int style, int hist, int eval,
                       int warnings, int chrono, FILE *outfile);
void InitOutputContexts(Context *ctxt, int nmask, int style, int hist, int eval,
                        int warnings, int chrono,  FILE *outfile);
void SetDir(Context *ctxt, int nmask, int dir);
void ResetDetect(Context *ctxt, int nmask);
void StartTimer(Context *ctxt);
void GetTime(Context *ctxt);
void PrintTime(Context *ctxt);
void ResetTimer(Context *ctxt, int nmask);
void GetBgn(Context *ctxt, int detect, int *bgn);
void StoreBgn(Context *ctxt, int bgn);
void StoreInterval(Context *ctxt, int bgn, int range);
int  SetSeqContext(Context *ctxt);

int  MaskSearch(Context *ctxt, int levels);

void StartStatus(void);
void PrintStatus(Context *ctxt, unsigned int *nscans);
void ClearStatus(void);
void fPrintfProcessLog(Context *ctxt, int nctxt, int total);

/*=============================================================================
 NewContext(): Cree une structure 'Context'.
=============================================================================*/

Context *NewContext(void)
{
  Context *ctxt = (Context *) malloc(sizeof(Context));

  return ctxt;
}
/*=============================================================================
 SetupContext(): Cree un tableau de 'nctxt' structures 'Context' chacune etant
                 associee au masque du meme indice dans le tableau de 'Mask'.
                 Des champs qui ne sont pas susceptibles de varier sont initia-
                 lises, d'autres prennent leur valeur par defaut.

  Note:
  'ctxt[i].range' est l'intervalle sur lequel un masque peut etre DETECTE,
  il sera RECHERCHE sur l'intervalle:
  ctxt[i].range - ctxt->mask->min_bgn - mask[i].min_len + 1;

  Si les 2 masques impliques dans la determination de 'ctxt[i].range' ont le
  meme debut que le pattern associe ou si
  aucun gap ne precede aucun des 2 masques, l'expression precedente a pour
  valeur:
  1 + 2 * SEARCH_EXT
=============================================================================*/

Context *SetupContext(Mask *mask, int nctxt)
{
  Context *ctxt;
  int     i, len = mask->pattern->max_len;

  extern int ScoreTabLen;

  if (ScoreTabLen < 5 * len)  ScoreTabLen = 5 * len;     /* borne inferieure */

  ctxt = (Context *) malloc(nctxt * sizeof(Context));

  for (i = 0; i < nctxt; i++)
  {
      ctxt[i].level = i;   
      ctxt[i].mask = mask + i;
      ctxt[i].bgn = 0;

      if (i == 0) {
          ctxt[i].range = SEQ_MAX_LEN;
          ctxt[i].scoretablen = ScoreTabLen;
      }
      else {                                         /* etendue plus reduite */

          ctxt[i].range = mask[0].max_bgn - mask[0].min_bgn +
                          mask[i].max_bgn +
                          mask[i].min_len + 2 * SEARCH_EXT ;

          ctxt[i].scoretablen = 3 * len + 2 * SEARCH_EXT ;
      }
      ctxt[i].dir = FORWARD;
      ctxt[i].detects = 0;
      ctxt[i].cumuls = 0;               /* mis a 0 une seule fois en general */
      ctxt[i].timer = OFF;
      ctxt[i].bgntime = 0;
      ctxt[i].duration = 0.0;
      ctxt[i].list = NULL;

                                          /* voir le code de 'SetSeqContext' */
      ctxt[i].tabscanlen = ctxt[i].scoretablen - mask[i].max_len + 1;
                             /* sera, en fait, seulement utilise au niveau 0 */

      if (ctxt[i].tabscanlen <= 0)
      {
          fprintf(stderr, "SetupContext: 'scoretablen' too short at level %d,",
                  ctxt[i].level + 1);          /* output: 'level' vaut 1,2.. */
          fprintf(stderr, " exit ..\n");
          exit(1);
      }
  }

  return ctxt;
}
/*=============================================================================
 InitSeqContext(): Initialise le champ 'ctxt->seq' d'un 'Context' pointe par
                   'ctxt' en lui associant une structure 'Sequence'.
=============================================================================*/

void InitSeqContext(Context *ctxt, Sequence *seq)
{
  ctxt->seq = seq;

  return;
}
/*=============================================================================
 InitSeqContexts(): Initialise le champ 'Sequence' d'un tableau de 'nctxt'
                    structures 'Context'. La sequence associee est commune a
                    l'ensemble des contextes.
=============================================================================*/

void InitSeqContexts(Context *ctxt, int nctxt, Sequence *seq)
{
  int i;

  for (i = 0; i < nctxt; i++)  ctxt[i].seq = seq;

  return;
}
/*=============================================================================
 InitOutputContext(): Charge depuis ses arguments le contexte de sortie des re-
                      sultats d'une recherche.
                      'style': LONG SHORT ou MUTE
                      'hist': ON ou OFF
                      'eval': ON ou OFF
                      'warnings': ON ou OFF
                      'chrono': ON ou OFF
                      'outfile': fichier de sortie: stdout, /dev/null ...
=============================================================================*/

void InitOutputContext(Context *ctxt, int style, int hist, int eval,
                       int warnings, int chrono, FILE *outfile)
{
  if (eval == OFF && hist == ON) {
      fprintf(stderr, "InitOutputContext: incompatible args 3 and 4, exit..\n");
      exit(1);
  }
  ctxt->style = style;
  ctxt->hist = hist;
  ctxt->eval = eval;
  ctxt->warnings = warnings;
  ctxt->timer = chrono;
  ctxt->outfile = outfile;

  if (hist == ON) InitDetectsHisto(ctxt);                 /* voir 'dhisto.c' */

  return;
}
/*=============================================================================
 InitOutputContexts(): Variante de 'InitOutputContext' pour un tableau de
                       structures 'Context'.
                       Tous sauf le dernier 'nmask - 1' ont la sortie MUTE et
                       les champs 'hist' et 'eval' mis a 'OFF',
                       tous sauf le 1er masque ont le chronometrage OFF.
                       Les autres arguments sont communs.
=============================================================================*/

void InitOutputContexts(Context *ctxt, int nmask, int style, int hist, int eval,
                        int warnings, int chrono, FILE *outfile)
{
  int i;

  if (eval == OFF && hist == ON) {
      fprintf(stderr, "InitOutputContexts: incompatible args 4 and 5, exit..\n");
      exit(1);
  }

  if (nmask == 1) {
      InitOutputContext(ctxt, style, hist, eval, warnings, chrono, outfile);
  }
  else {
      InitOutputContext(ctxt, MUTE, OFF, OFF, warnings, chrono, outfile);
      for (i = 1; i < nmask - 1; i++)
      {
          InitOutputContext(ctxt + i, MUTE, OFF, OFF, warnings, OFF, outfile);
      }
      InitOutputContext(ctxt + nmask - 1, style, hist, eval, warnings,
                        OFF, outfile);
  }

  return;
}
/*=============================================================================
 SetDir(): initialise les champs 'ctxt->dir' d'un tableau de contextes pointe
           par "ctxt'.
=============================================================================*/

void SetDir(Context *ctxt, int nmask, int dir)
{
  int i;

  for (i = 0; i < nmask; i++)  ctxt[i].dir = dir;

  return;
}
/*=============================================================================
 ResetDetect(): Remet a 0 le compteur des detections 'ctxt->detects' d'un ta-
                bleau de 'nmask' structures 'Context'.
=============================================================================*/

void ResetDetect(Context *ctxt, int nmask)
{
  int i;

  for (i = 0; i < nmask; i++)  ctxt[i].detects = 0;

  return;
}
/*=============================================================================
 StartTimer(): demarre le 'timer' d'une structure 'Context'.
=============================================================================*/

void StartTimer(Context *ctxt)
{
  if (ctxt->timer == ON) ctxt->bgntime = clock();
  return;
}
/*=============================================================================
 GetTime(): Enregistre le temps d'une structure 'Context'.
=============================================================================*/

void GetTime(Context *ctxt)
{
  if (ctxt->timer == ON) 
      ctxt->duration += CpuTime(&(ctxt->bgntime));
  return;
}
/*=============================================================================
 PrintTime(): Affiche le temps d'une structure 'Context'.
=============================================================================*/

void PrintTime(Context *ctxt)
{
  if (ctxt->timer == ON) fPrintfCpuTime(ctxt->duration, ctxt->outfile);
  return;
}
/*=============================================================================
 ResetTimer(): Remet a 0 le 'timer' d'un tableau de structures 'Context'.
=============================================================================*/

void ResetTimer(Context *ctxt, int nmask)
{
  int i;

  for (i = 0; i < nmask; i++) {
      ctxt[i].bgntime = 0;
      ctxt[i].duration = 0.0;
  }
  return;
}
/*=============================================================================
 GetBgn(): Determine l'abscisse voisine de la detection du 1er masque a la
           position 'detect' dans une sequence, cette abscisse pointee par
           'bgn' sera prise pour origine pour les tableaux de scores de la 
           recherche de 'ctxt[1].mask'.
           Tous les champs 'ctxt[i].bgn' sauf le premier auront la meme valeur,
           determinee au 1er niveau.
=============================================================================*/

void GetBgn(Context *ctxt, int detect, int *bgn)
{
  int level = ctxt->level;            /* niveau ou la detection est observee */

  if (level == 0)
      *bgn = detect - ctxt->mask->max_bgn - SEARCH_EXT ;
  else
      *bgn = ctxt->bgn;                                      /* simple copie */

  if (*bgn < 0) *bgn = 0;                    /* eviter de sortir des donnees */

  return;
}
/*=============================================================================
 StoreBgn(): Enregistre dans une structure 'Context' l'origine de l'intervalle
             de recherche dans une sequence.
=============================================================================*/

void StoreBgn(Context *ctxt, int bgn)
{
  ctxt->bgn = bgn;

  return;
}
/*=============================================================================
 StoreInterval(): Enregistre dans une structure 'Context' les parametres d'un
                  intervale de recherche ('bgn', 'range') dans une sequence.
=============================================================================*/

void StoreInterval(Context *ctxt, int bgn, int range)
{
  ctxt->bgn = bgn;
  ctxt->range = range;

  return;
}
/*=============================================================================
 SetSeqContext(): Initialise depuis la structure pointee par 'ctxt->seq' les
                  champs de la structure pointee par 'ctxt' utiles a son examen,
                  a l'aide de 'ctxt->bgn' et 'ctxt->range', respectivement le
                  debut et l'etendue d'une recherche.
                  On traite de maniere differente le premier niveau (recherche
                  sur toute une sequence) et les niveaux suivants.
                  Retourne 0 si le nombre de bases est insuffisant, sinon 1. 
=============================================================================*/

int SetSeqContext(Context *ctxt)
{
                                            /* origine des donnees a scanner */
                                       /* pour le 1er niveau 'ctxt->bgn' est */
                                           /* initialise dans 'SetupContext' */
  ctxt->data = ctxt->seq->data + ctxt->bgn;

                                    /* longueur a scanner, depend du  niveau */

  ctxt->datalen = MIN(ctxt->seq->datalen - ctxt->bgn, ctxt->range);

                                           /* longueur effectivement scannee */
  if (ctxt->level == 0)
  {
      ctxt->seqscanlen = ctxt->datalen - ctxt->mask->min_len + 1;

      ctxt->nscan = ctxt->seqscanlen / ctxt->tabscanlen;
      ctxt->lastscanlen = ctxt->seqscanlen % ctxt->tabscanlen;

      ctxt->lastscoretabLen = ctxt->lastscanlen + ctxt->mask->max_len - 1;
  }
  else
  {
      ctxt->seqscanlen = ctxt->datalen - ctxt->mask->min_bgn -
                         ctxt->mask->min_len + 1;

      ctxt->nscan = 0;
      ctxt->lastscanlen = ctxt->seqscanlen;

      ctxt->lastscoretabLen = ctxt->lastscanlen + ctxt->mask->min_bgn +
                              ctxt->mask->max_len - 1;
  }

  if (ctxt->seqscanlen <= 0)                       /* sequence trop courte ! */
  {
      if (ctxt->warnings == ON)
          fprintf(stderr, "SetSeqContext: too few bases in sequence '%d'\n",
                  ctxt->seq->nb);
      return 0;
  }

  return 1;
}
/*=============================================================================
 MaskSearch(): Fonction de recherche de masques operant en plusieurs etapes:

               - la 1ere, de preference sur un masque presentant peu de gaps,
               - la 2eme, qui demarre a chaque detection dans la 1ere etape,
               opere sur un masque plus structure et un intervalle reduit.
               - etc ..
               Cette organisation repond au besoin d'optimisation de la duree
               de l'operation de recherche.

               Les arguments: sequence, masques, parametres .. sont pointes par
               la structure Context 'ctxt'.
               'levels' est le nombre (generalement 1, 2 ou 3) de niveaux de la 
               recherche, c'est aussi le nombre de masques recherches.

               note: La variable 'detect' repere une detection depuis l'origine
               de la sequence visitee (et non pas depuis 'ctxt->bgn').
               Retourne le nombre d'occurences du masque recherche au niveau le
               plus eleve (d'indice 'levels - 1'). Les occurences sont comptees
               dans les champs 'ctxt->detects' qui doivent etre regulierement
               remis a 0 et cumulees dans 'ctxt->cumuls'.
=============================================================================*/

int MaskSearch(Context *ctxt, int levels)
{
  int     i, j, k, l, lo, detect, bgn, TabLen, ScanLen;

  extern unsigned long long int TotalNucScans;               /* voir 'env.c' */
  static unsigned int jobscans = 0;

  if ( SetSeqContext(ctxt) == 0 )  return 0;       /* sequence trop courte ! */

  TabLen  = ctxt->scoretablen;
  ScanLen = ctxt->tabscanlen;

  for (i = k = 0; i <= ctxt->nscan; i++, k += ctxt->tabscanlen)
  {
      if (i == ctxt->nscan)                         /* phase finale atteinte */
      {
          TabLen  = ctxt->lastscoretabLen;
          ScanLen = ctxt->lastscanlen;

          if (TabLen == 0) break;
      }
                                                       /* tableaux de scores */
      GetMaskScoresTab(ctxt->data + k, TabLen, ctxt->mask, ctxt->level);

      if (ctxt->level == 0) {
          lo = 0;
          TotalNucScans += ScanLen;           /* decompte des bases traitees */
      }
      else lo = ctxt->mask->min_bgn;

      for (j = 0, l = lo; j < ScanLen; j++, l++)
      {
          if (GetMaskScore(ctxt->mask, l, ctxt->level) > ctxt->mask->threshold)
          {
              ctxt->detects++;
              ctxt->cumuls++;
              detect = ctxt->bgn + k + l;

              if (levels == 1) {
                  AddLink(ctxt, detect, ctxt->mask->cfgindex, ctxt->mask->score);
              }
              else                   /* relance MaskSearch dans un voisinage */
              {                                           /* de la detection */
                  GetBgn(ctxt, detect, &bgn);
                  StoreBgn(ctxt + 1, bgn);

                  (int) MaskSearch(ctxt + 1, levels - 1);
              }
          }
      }
      PrintStatus(ctxt, &jobscans);
  }

  return ctxt[levels - 1].detects;
}
/*=============================================================================
 StartStatus(): demarre l'affichage de l'etat d'avancement d'une recherche.
                Il faut verifier la coherence entre les 3 fonctions suivantes,
                relatives au meme affichage.
=============================================================================*/

void StartStatus(void)
{
  fprintf(stderr, "Kb: %6d\r", 0);                      /* debut du comptage */
  return;
}
/*=============================================================================
 PrintStatus(): Affiche l'etat d'avancement du premier niveau d'une operation
                de recherche: le nombre de bases traitees a la periodicite de
                10 fois la longueur traitee par un tableau de scores precalcule.
                Le curseur de texte est replace en debut de ligne.
                A l'occasion le temps CPU est enregistre.
=============================================================================*/

void PrintStatus(Context *ctxt, unsigned int *nscans)
{
  if (ctxt->level == 0 && ++(*nscans) % 10 == 0)
  {
      fprintf(stderr, "Kb: %6d\r", *nscans * ctxt->tabscanlen / 1000);
      GetTime(ctxt);
  }
  return;
}
/*=============================================================================
 ClearStatus(): Efface l'etat d'avancement d'une recherche de 'erpin'.
                Fonction executee en entree de 'fPrintfPostProcess' de 'list.c'
=============================================================================*/

void ClearStatus(void)
{
  fprintf(stderr, "          \r");
  return;
}
/*=============================================================================
 fPrintfProcessLog(): Affichage recapitulatif d'une operation de recherche.
=============================================================================*/

void fPrintfProcessLog(Context *ctxt, int nctxt, int total)
{
  int i, n;
  extern unsigned long long int TotalNucScans;               /* voir 'env.c' */

  if (ctxt[nctxt-1].style != MUTE) fprintf(ctxt->outfile, "\n");

  for (i = 0; i < nctxt; i++)
  {
      fprintf(ctxt->outfile, "-------- at level %d --------\n", i+1);
      if (i == 0)
          fprintf(ctxt->outfile, "%lld bases processed\n", TotalNucScans);
      fprintf(ctxt->outfile, "cutoff: %.2f\n", ctxt[i].mask->threshold);

#ifdef DEBUG

      if (i > 0) 
          fprintf(ctxt->outfile, "range: %d\n", ctxt[i].seqscanlen);

fprintf(ctxt->outfile, "range: %d\n", ctxt[i].seqscanlen);
fprintf(ctxt->outfile, "DataLen: %d\n", ctxt[i].datalen);
fprintf(ctxt->outfile, "MaskMinLen: %d\n", ctxt[i].mask->min_len);
fprintf(ctxt->outfile, "MaskMinBgn: %d\n", ctxt[i].mask->min_bgn);
fprintf(ctxt->outfile, "CtxtBgn: %d\n", ctxt[i].bgn);
fprintf(ctxt->outfile, "CtxtRange: %d\n", ctxt[i].range);
fprintf(ctxt->outfile, "CtxtLevel: %d\n", ctxt[i].level);

#endif

      n = ctxt[i].mask->ncfg;
      fprintf(ctxt->outfile, "%d config. per site\n", n);

      n = ctxt[i].cumuls;
      fprintf(ctxt->outfile, "%d hit%s\n", n, n > 1 ? "s" : "");
  }
  if (total > 0)
      fprintf(ctxt->outfile,
              "%d independent hit%s\n", total, total > 1 ? "s" : "");
  
  PrintTime(ctxt);
  PrintScoresHisto(ctxt + nctxt - 1);       /* histogramme du dernier masque */
  printf("\n");

  return;
}
/*===========================================================================*/
