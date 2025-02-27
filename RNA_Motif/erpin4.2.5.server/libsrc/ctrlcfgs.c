

/*=============================================================================
 ctrlcfgs.c                                        A.Lambert le 04/12/02

 Fonctions concernant le denombrement et le controle du nombre de configura-
 tions des 'masks' et 'pattern'.

 cc -O2 -Wall -c ctrlcfgs.c -I../include ;

 ar -rs ../lib/librnaIV.a ctrlcfgs.o ;

=============================================================================*/

#include "rnaIV.h"

void  GetPatternCfgsNb(Pattern *pattern);
void  GetMaskCfgsNb(Mask *mask);
int   MaskCfgsNb(Mask *mask);
void  PrintMaskCfgsVol(Mask *mask, FILE *txt);
void  CtrlCfgs(Mask *mask, int nmask, int maskproc, char mode, FILE *txt);
void  CfgMsg(int maskproc, FILE *txt);
int   Rnd(int min, int max);
short *GetRandCfg(Mask *mask);

/*=============================================================================
 GetPatternCfgsNb: Calcule le nombre de configurations de la region pointee par
                   'pattern'.
                   Ce nombre, exprime en unites logarithmiques de base 2 (bits),
                   est enregistre dans le champ 'pattern->log2ncfg'.
                   Si ce nombre est inferieur a 'CFG_MAX_BITS', sa valeur deci-
                   male est enregistree dans le champ 'pattern->ncfg', sinon 0.
=============================================================================*/

void GetPatternCfgsNb(Pattern *pattern)
{
  double log2ncfg;
  int    i, n;

  for (i = 0, log2ncfg = 0.0; i < pattern->natom; i++) 
      log2ncfg += LOG2(pattern->atom[i].max_gaps + 1.);

  pattern->log2ncfg = log2ncfg;
  pattern->ncfg = 0;

  if (log2ncfg <= CFG_MAX_BITS)
  {
      for (i = 0, n = 1; i < pattern->natom; i++) 
          n *= (pattern->atom[i].max_gaps + 1);
      pattern->ncfg = n;
  }
  return ;
}
/*=============================================================================
 GetMaskCfgsNb(): Calcule le nombre de configurations du masque pointe par
                  'mask'. Retourne ce nombre exprime en unites logarithmiques
                  de base 2.
                  Si ce nombre est inferieur a 'CFG_MAX_BITS', en retour 'ncfg'
                  pointe sa valeur decimale sinon 0.
                  Les tableaux utilises sont crees par 'GetMaskGapList' qui
                  devra etre appelee avant, aucun controle n'est effectue.
                  ('GetMaskGapList' est executee par 'GetMask' lors de la crea-
                  tion d'un masque).
=============================================================================*/

void GetMaskCfgsNb(Mask *mask)
{
  double log2ncfg;
  int    j, n;

  for (j = 0, log2ncfg = 0.0; j < mask->ngaps; j++)
     log2ncfg += LOG2(mask->gapslist[2][j] - mask->gapslist[1][j] + 1.);

  mask->log2ncfg = log2ncfg;
  mask->ncfg = 0;

  if (log2ncfg <= CFG_MAX_BITS)
  {
      for (j = 0, n = 1; j < mask->ngaps; j++)
          n *= (int) (mask->gapslist[2][j] - mask->gapslist[1][j] + 1);
      mask->ncfg = n;
  }
  return ;
}
/*=============================================================================
 MaskCfgsNb(): Variante de 'GetMaskCfgsNb' qui suppose que les situations com-
               portant un nombre de configurations trop important ont ete ecar-
               tees. Le nombre en unites log2 n'est pas calcule.
               Le nombre de configurations est retourne.
               Cette fonction sera utilise dans 'ImportPatternCfgInfos' dans
               'dmp.c'.
=============================================================================*/

int MaskCfgsNb(Mask *mask)
{
  int  j, ncfg = 1;

  for (j = 0; j < mask->ngaps; j++)
      ncfg *= (int) (mask->gapslist[2][j] - mask->gapslist[1][j] + 1);

  return ncfg;
}
/*=============================================================================
 PrintMaskCfgsVol(): Affiche le nombre de configurations d'un masque et le vo-
                     lume occupe par le tableau des configurations.
                     La constuction effective du tableau n'est pas necessaire a
                     ce calcul.
                     Le volume est exprime en KB, si le nombre de configurations
                     depasse 'CFG_MAX_BITS' les unites logarithmiques de base 2
                     sont utilisees.
                     La sortie est effectuee dans le fichier pointe par 'txt'.
=============================================================================*/

void PrintMaskCfgsVol(Mask *mask, FILE *txt)
{
  double v;                       /* v: volume en Koctet d'une configuration */

  v = 1.*(sizeof(Config) + 3 * (mask->nhx + mask->nst + 2) * sizeof(short));
  v /= 1024.;

  if (mask->log2ncfg <= CFG_MAX_BITS)
  {
      fprintf(txt, "config.number: %d\n", mask->ncfg);
      fprintf(txt, "config.table volume: %.1f KB\n", mask->ncfg * v);
  }

  else          /* affichage en 'float' pour prevenir un overflow sur 'ncfg' */
  {
      fprintf(txt, "config.number: %.1e\n", pow(2., mask->log2ncfg));
      fprintf(txt, "config.table volume: %.1e MB\n",
              pow(2., mask->log2ncfg - 10) * v);
  }
  return;
}
/*=============================================================================
 CtrlCfgs(): Procede au controle du nombre de configurations de 'nmask' masques
             (le premier pointe par 'mask') dans l'ordre ou ils seront traites
             dans une operation de recherche.
             Les types de traitement consideres sont 'DYNAMIC' et 'STATIC', 
             l'argument 'maskproc' prenant une de ces 2 valeurs.
             Dans le cas 'DYNAMIC' on opere une simulation de recherche qui ne
             necessite pas la creation effective des tableaux de configurations.
             Le comportement va, en fonction du nombre de configurations trouve
             et du volume de memoire correspondant requis, du calcul 'silencieux'
             a l'arret inconditionnel, en passant par l'affichage d'informations
             et/ou de 'Warning'.
             Si l'argument 'mode' est 'V' l'information sera obligatoirement 
             fournie.
             Les sorties sont dirigees sur le fichier pointe par 'txt'.
=============================================================================*/

void CtrlCfgs(Mask *mask, int nmask, int maskproc, char mode, FILE *txt)
{
  Context  *ctxt;
  short    *gapcfg;
  int      i, stops, warns, infos;
  char     *warnstr, *stopstr;

  warnstr = (char *) calloc(64, 1);
  stopstr = (char *) calloc(64, 1);

  if (maskproc != STATIC && maskproc != DYNAMIC) {
      fprintf(stderr, "CtrlCfgs: Unknown argument #3, exit..\n");
      exit(1);
  }

  stops = warns = infos = 0;

  ctxt = SetupContext(mask, nmask);

  if (maskproc == DYNAMIC && nmask > 1)        /* simulation de la reduction */
  {                                        /* progressive des configurations */
      gapcfg = GetRandCfg(mask);
      ExportMaskCfgInfos(ctxt, gapcfg);
      free(gapcfg);
      for (i = 1; i < nmask; i++)
      {
          ImportPatternCfgInfos2(ctxt + i);
          GetMaskCfgsNb(ctxt[i].mask);
          if (i != nmask - 1)
          {
              gapcfg = GetRandCfg(ctxt[i].mask);
              ExportMaskCfgInfos(ctxt + i, gapcfg);
              free(gapcfg);
          }
      }
  }

  for (i = 0; i < nmask; i++)                          /* examen des masques */
  {
      if (mask[i].log2ncfg > CFG_BITS_INFO) {
      infos++;
      }
      if (mask[i].log2ncfg > CFG_BITS_WARN) {
          sprintf(warnstr, "%s#%d, ", warnstr, i+1);
          warns++;
      }
      if (mask[i].log2ncfg > CFG_MAX_BITS) {
          sprintf(stopstr, "%s#%d, ", stopstr, i+1);
          stops++;
      }
  }
                                               /* affichage des informations */
  if (infos > 0 || toupper(mode) == 'V' )
  {
      fprintf(txt, "\nMEMORY USAGE CONTROL: %d step%s to be processed\n",
              nmask, nmask > 1 ? "s" : "");
      for (i = 0; i < nmask; i++) {
          fprintf(txt, "step #%d:\n", i+1);
          PrintMaskCfgsVol(mask + i, txt);
      }
  }                            /* affichage precedant la sortie du programme */
  if (stops > 0) {
      fprintf(stderr, "\nWARNING: Too many configurations for step%s %s\n",
              stops > 1 ? "s" : "", stopstr);
      CfgMsg(maskproc, stderr);
      fprintf(stderr, "exit..\n");
      exit(1);
  }                /* affichage precedant la sortie optionnelle du programme */
  if (warns > 0) {
      fprintf(stderr, "\nWARNING: many configurations for step%s %s\n",
              warns > 1 ? "s" : "", warnstr);
      CfgMsg(maskproc, stderr);
      fprintf(stderr, "may be you prefer to kill the program..\n");
      sleep(5);
  }
  if (maskproc == DYNAMIC)                /* restauration du contenu initial */
      for (i = 0; i < nmask; i++) {                           /* des masques */
          ResetMaskGapList(mask + i);                     /* avant la sortie */
          GetMaskCfgsNb(mask + i);
      }

  free(warnstr);
  free(stopstr);
  free(ctxt);
  return;
}
/*=============================================================================
 CfgMsg(): Affiche vers 'txt' un commentaire, suivant la valeur de 'maskproc'.
=============================================================================*/

void CfgMsg(int maskproc, FILE *txt)
{
  if (maskproc == DYNAMIC)
      fprintf(txt, "It is advisable to use extra search steps.\n");
  else
      fprintf(txt,"It is advisable to use multi-step search (masks).\n");

  fprintf(txt, "Please consult Erpin documentation about search strategies.\n");

  return;
}
/*============================================================================
 Rnd(): genere un nombre pseudo aleatoire dans [min, max] (bornes incluses).
============================================================================*/

int Rnd(int min, int max)
{
 return  (min + (int) ((float) (max - min + 1) * rand() / (RAND_MAX + 1.0)));
}
/*=============================================================================
 GetRandCfg(): Extrait une configuration aleatoire de gaps parmi celles d'un
               masque pointe par 'mask'.
               Retourne un pointeur sur la configuration creee.
=============================================================================*/

short *GetRandCfg(Mask *mask)
{
  int   i;
  short *gapcfg;

  gapcfg = (short *) malloc(mask->ngaps * sizeof(short));

  for (i = 0; i < mask->ngaps; i++)
      gapcfg[i] = Rnd(mask->gapslist[1][i], mask->gapslist[2][i]);

  return gapcfg;
}
/*===========================================================================*/
