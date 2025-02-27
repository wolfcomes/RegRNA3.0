
/*=============================================================================
 cfgs.c                                     A.Lambert  le 04/12/02

 Pour un ensemble de masques d'une meme region, ce programme calcule le nombre
 de configurations et la memoire requise dans les 2 types de gestion, statique
 et dynamique, lors d'une exploration de sequences par 'erpin'.

 Il permet ainsi de selectionner les groupes d'arguments les mieux adaptes,
 c'est a dire les moins exigeants en temps CPU et/ou ressources en memoire.

 make -f apps.mk cfgs ;
 mv cfgs ../bin ;

 cc -O3 -Wall -o cfgs cfgs.c -I../include -L../lib -lrnaIV -lm ;
 chmod 755 cfgs ;
 strip cfgs ;

 cfgs <trset> <region>
      -nomask|((-mask|-umask|-add) <elt1> <elt2>..)   level1,
      [-nomask|((-mask|-umask|-add) <elt1> <elt2>..)] level2,  default: void
      [...]                                           level..           idem

 cfgs ~/devc/projets/bio/data/trsetsII/trna.db -2,2 -umask 7 20 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/trna.db -2,2 -umask 7 8 9 20 11 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/trna.db -2,2 -mask 9 20 11 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/trna98.db 1,12 -mask 2 9 20 11 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/trna98.db 1,12 -umask 2 9 20 11 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/trna98.db 1,12 \
      -umask 7 20 -umask 3 4 5 6 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/trna98.db -2,2 \
      -umask 7 20 -mask 2 3 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/secis.db -2,2 -umask 3 4 -nomask ;

 cfgs ~/devc/projets/bio/data/trsetsII/mirna1.epn -2,2 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/mirna1.epn 1,27 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/16s.epn 1,289 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/mirna1.epn -2,2 \
      -umask 10 13 15 17 -add 8 11 19 -add 6 9 21 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/mirna1.epn -2,2 \
      -umask 10 -add 13 15 17 -add 8 11 19 -add 6 9 21 -nomask ;
 cfgs ~/devc/projets/bio/data/trsetsII/mirna2.epn -2,2 \
      -umask 10 -add 13 -add 8 11 19 -add 6 9 21 -nomask ;

=============================================================================*/

#include "rnaIV.h"

void  DemoCtrlCfgs(Mask *mask, int nmask, FILE *txt);
void  cfgsHelp(void);

/*===========================================================================*/

int main(int argc, char *argv[])
{
  Trset    *trset;
  Pattern  *pattern;
  Mask     *mask;
  int      nmask;

  ReadHelpArgs(argc, argv, 3, cfgsHelp);
  fTest(argv[1]);

  mask = ReadMasksArgs(argc, argv, &nmask);

  trset = ReadTrset(argv[1], 'S', stdout);
  pattern = ReadPattern(argv[2], trset);
  ParseMasksArgs(mask, nmask, pattern);

  GetMasks(mask, nmask, pattern);

  DemoCtrlCfgs(mask, nmask, stdout);

  DelTrset(trset);
  DelMasks(mask, nmask);
  DelPattern(pattern);
  exit(0);
}
/*=============================================================================
 DemoCtrlCfgs(): Variante de la fonction 'CtrlCfgs' dans 'ctrlcfgs.c'.
                 Procede au controle du nombre de configurations de 'nmask'
                 masques (le premier pointe par 'mask') dans l'ordre ou ils
                 seront traites dans une operation de recherche.

                 Les 2 types de traitement, 'DYNAMIC' et 'STATIC', sont
                 consideres.
                 Dans le cas 'DYNAMIC' on opere une simulation de recherche qui
                 ne necessite pas la creation effective des tableaux de confi-
                 gurations.
                 Les sorties sont dirigees sur le fichier pointe par 'txt'.
=============================================================================*/

void DemoCtrlCfgs(Mask *mask, int nmask, FILE *txt)
{
  Context  *ctxt;
  short    *gapcfg;
  int      i, stops, warns;
  char     *warnstr, *stopstr;

  warnstr = (char *) calloc(64, 1);
  stopstr = (char *) calloc(64, 1);

  stops = warns = 0;

  ctxt = SetupContext(mask, nmask);
  fprintf(txt, "\nMEMORY USAGE CONTROL IN STATIC AND DYNAMIC MASK PROCESSING\n");
  fprintf(txt, "%d step%s to be processed\n\n", nmask, nmask > 1 ? "s" : "");

  fprintf(txt, "===========================================================\n");
  fprintf(txt, "I/ STATIC MASK PROCESSING (Erpin '-smp' option):\n");
  fprintf(txt, "===========================================================\n");

  for (i = 0; i < nmask; i++) {
      fprintf(txt, "step #%d:\n", i+1);
      PrintMaskCfgsVol(mask + i, txt);
  }

  for (i = 0; i < nmask; i++)                          /* examen des masques */
  {
      if (mask[i].log2ncfg > CFG_BITS_WARN) {
          sprintf(warnstr, "%s#%d, ", warnstr, i+1);
          warns++;
      }
      if (mask[i].log2ncfg > CFG_MAX_BITS) {
          sprintf(stopstr, "%s#%d, ", stopstr, i+1);
          stops++;
      }
  }
  if (stops > 0) {
      fprintf(txt, "\nWARNING: Too many configurations for step%s %s\n",
              stops > 1 ? "s" : "", stopstr);
      CfgMsg(STATIC, txt);
  }
  else if (warns > 0) {
      fprintf(txt, "\nWARNING: many configurations for step%s %s\n",
              warns > 1 ? "s" : "", warnstr);
      CfgMsg(STATIC, txt);
  }
  fprintf(txt, "\n");
  fprintf(txt, "===========================================================\n");
  fprintf(txt, "II/ DYNAMIC MASK PROCESSING (Erpin '-dmp' default option):\n");
  fprintf(txt, "===========================================================\n");

  if (nmask > 1)                               /* simulation de la reduction */
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
  for (i = 0; i < nmask; i++) {
      fprintf(txt, "step #%d:\n", i+1);
      PrintMaskCfgsVol(mask + i, txt);
  }

  stops = warns = 0;
  warnstr[0] = stopstr[0] = '\0';

  for (i = 0; i < nmask; i++)                          /* examen des masques */
  {
      if (mask[i].log2ncfg > CFG_BITS_WARN) {
          sprintf(warnstr, "%s#%d, ", warnstr, i+1);
          warns++;
      }
      if (mask[i].log2ncfg > CFG_MAX_BITS) {
          sprintf(stopstr, "%s#%d, ", stopstr, i+1);
          stops++;
      }
  }
  if (stops > 0) {
      fprintf(txt, "\nWARNING: Too many configurations for step%s %s\n",
              stops > 1 ? "s" : "", stopstr);
      CfgMsg(DYNAMIC, txt);
  }
  else if (warns > 0) {
      fprintf(txt, "\nWARNING: many configurations for step%s %s\n",
              warns > 1 ? "s" : "", warnstr);
      CfgMsg(DYNAMIC, txt);
  } 
  fprintf(txt, "\n");
  free(warnstr);
  free(stopstr);
  free(ctxt);
  return;
}
/*=============================================================================
 cfgsHelp(): help
=============================================================================*/

void cfgsHelp(void)
{
  fprintf(stderr,

"cfgs: masks configurations and memory usage statistics\n"
"Usage:\n"
"cfgs [-h]\n"
"cfgs <training-set> <region> <mask>\n"
"     [(-mask|-umask|-add) <arg1> <arg2> ..][..]\n"
);

  return;
}
/*===========================================================================*/
