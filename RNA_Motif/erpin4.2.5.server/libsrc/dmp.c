
/*=============================================================================
 dmp.c                                       A.Lambert le 24/05/02

 Fonctions destinees a assurer la gestion dynamique des masques lors d'une ope-
 ration de recherche a plusieurs niveaux.

 cc -O2 -Wall -c dmp.c -I../include ;
 ar -rs ../lib/librnaIV.a dmp.o ;

=============================================================================*/

#include "rnaIV.h"

void  GetPatternGapList(Pattern *pattern, int nmask);
void  FreePatternGapList(Pattern *pattern, int nmask);
void  ResetPatternGapList(Pattern *pattern);
void  ResetGapLists(Context *ctxt, int bgn, int end);
void  ExportMaskCfgInfos(Context *ctxt, short *mcfg);
void  ImportPatternCfgInfos(Context *ctxt);
void  ImportPatternCfgInfos2(Context *ctxt);

int   RecMaskSearch(Context *ctxt, int levels);
int   Search(Context *ctxt, int nmask, int maskproc);

/*=============================================================================
 GetPatternGapList(): Cree les tableaux 'gapslist' et 'gapslist_ori' des gaps
                      de 'pattern', relatifs a un groupe de 'nmask' masques.
                      
                      Les elements de 'gapslist' seront ulterieurement modifies.
                      Cette fonction sera executee lors d'un appel a 'GetMasks'.
=============================================================================*/

void GetPatternGapList(Pattern *pattern, int nmask)
{
  int   i, j;

  if (pattern->nst == 0) return;

  pattern->gapslist = (short ***) malloc(nmask * sizeof(short **));

  for (i = 0; i < nmask; i++)
  {
      pattern->gapslist[i] = sMat(2, pattern->nst);
  }
  pattern->gapslist_ori = sMat(2, pattern->nst);


  for (i = j = 0; i < pattern->natom; i++)
  {
      if (pattern->atom[i].type == STD)
      {
          pattern->gapslist_ori[0][j] = 0;
          pattern->gapslist_ori[1][j] = (short) pattern->atom[i].max_gaps;
          j++;
      }
  }
                        /* initialise le 1er tableau 'gapslist' de 'pattern' */

  CopysMat(pattern->gapslist_ori, 2, pattern->nst, pattern->gapslist[0]);

  return;
}
/*=============================================================================
 FreePatternGapList(): Libere les tableaux de gaps associes a 'nmask' masques.
                       Cette fonction sera executee avec 'FreeMasks'.
=============================================================================*/

void FreePatternGapList(Pattern *pattern, int nmask)
{
  int i;

  if (pattern->gapslist_ori != NULL)
  {
      FreesMat(pattern->gapslist_ori);
      pattern->gapslist_ori = NULL;
  }
  if (pattern->gapslist != NULL)
  {
      for (i = 0; i < nmask; i++)
      {
          FreesMat(pattern->gapslist[i]);
      }
      free(pattern->gapslist);
      pattern->gapslist = NULL;
  }
  return;
}
/*=============================================================================
 ResetPatternGapList(): Remet les elements du 1er tableau 'gaplist' de 'pattern'
                        a leur valeur originale.
=============================================================================*/

void ResetPatternGapList(Pattern *pattern)
{

  CopysMat(pattern->gapslist_ori, 2, pattern->nst, pattern->gapslist[0]);

  return;
}
/*=============================================================================
 ResetGapLists(): Remet a leur valeur initiale l'ensemble des elts des champs
                  'gapslist' des masques et du pattern associes au tableau de
                  structures 'Context' pointe par 'ctxt'.
                  Plus precisement aux contextes d'indices allant de 'bgn' a
                  'end' inclus.
=============================================================================*/

void ResetGapLists(Context *ctxt, int bgn, int end)
{
  int  i;

  for (i = bgn; i <= end; i++) ResetMaskGapList((ctxt + i)->mask);

  ResetPatternGapList(ctxt->mask->pattern);

  return;
}
/*=============================================================================
 ExportMaskCfgInfos(): Modifie (reduit !) le tableau des gaps d'un pattern en
                       utilisant une configuration de gaps particuliere d'un
                       masque associe.
                       ATTENTION: 'ctxt' ne doit pas pointer le dernier element
                       d'un tableau, aucun controle n'est effectue.
=============================================================================*/

void ExportMaskCfgInfos(Context *ctxt, short *mcfg)
{
  int i, j, col;

  i = ctxt->level;

  for (j = 0; j < ctxt->mask->ngaps; j++)
  {
      col = (int) ctxt->mask->gapslist[0][j];
      if (col != -1 )
      {
          ctxt->mask->pattern->gapslist[i][0][col] =
          ctxt->mask->pattern->gapslist[i][1][col] = mcfg[j];
      }
  }
              /* copie: attention que 'i' n'atteigne pas le dernier indice ! */

  CopysMat(ctxt->mask->pattern->gapslist[i], 2, ctxt->mask->pattern->nst,
           ctxt->mask->pattern->gapslist[i+1]);

  return;
}
/*=============================================================================
 ImportPatternCfgInfos(): Modifie le tableau des gaps d'un masque en utilisant
                          les informations (cumulees) du tableau des gaps du
                          pattern associe.
                          le champ 'ctxt->mask->ncfg' est aussi actualise.
=============================================================================*/

void ImportPatternCfgInfos(Context *ctxt)
{
  int   i, j, col;

  i = ctxt->level;

  for (j = 0; j < ctxt->mask->ngaps; j++)
  {
      col = (int) ctxt->mask->gapslist[0][j];
      if (col != -1)
      {
          ctxt->mask->gapslist[1][j] = ctxt->mask->pattern->gapslist[i][0][col];
          ctxt->mask->gapslist[2][j] = ctxt->mask->pattern->gapslist[i][1][col];
      }
  }
  ctxt->mask->ncfg = MaskCfgsNb(ctxt->mask);

  return;
}
/*=============================================================================
 ImportPatternCfgInfos2(): Variante de 'ImportPatternCfgInfos' qui utilise la
                           fonction 'GetMaskCfgsNb' laquelle calcule d'abord le
                           nombre de configurations en unites log2.
=============================================================================*/

void ImportPatternCfgInfos2(Context *ctxt)
{
  int   i, j, col;

  i = ctxt->level;

  for (j = 0; j < ctxt->mask->ngaps; j++)
  {
      col = (int) ctxt->mask->gapslist[0][j];
      if (col != -1)
      {
          ctxt->mask->gapslist[1][j] = ctxt->mask->pattern->gapslist[i][0][col];
          ctxt->mask->gapslist[2][j] = ctxt->mask->pattern->gapslist[i][1][col];
      }
  }
  GetMaskCfgsNb(ctxt->mask);

  return;
}
/*=============================================================================
 RecMaskSearch(): Fonction de recherche de masques operant en plusieurs etapes:

               - la 1ere, de preference sur un masque presentant peu de gaps,
               - la 2eme, qui demarre a chaque detection dans la 1ere etape,
               opere sur un masque plus structure et un intervalle reduit,
               utilisant la configuration trouvee a l'etape precedente.
               - etc ..
               Cette organisation repond au besoin d'optimisation de la duree
               de l'operation de recherche.

               Les arguments: sequence, masques, parametres .. sont pointes par
               la structure Context 'ctxt'.
               'levels' est le nombre (generalement 1, 2 ou 3 mais eventuelle-
               ment plus) de niveaux de la recherche, c'est aussi le nombre de
               masques recherches.

               note: La variable 'detect' repere la position d'une detection
               depuis l'origine de la sequence visitee (et non pas depuis 
               'ctxt->bgn').
               Retourne le nombre d'occurences du masque recherche au niveau le
               plus eleve (d'indice 'levels - 1'). Les occurences sont comptees
               dans les champs 'ctxt->detects' qui doivent etre regulierement
               remis a 0 et cumulees dans 'ctxt->cumuls'.
=============================================================================*/

int RecMaskSearch(Context *ctxt, int levels)
{
  int     i, j, k, l, lo, detect, bgn, index, TabLen, ScanLen;

  extern unsigned long long int TotalNucScans;
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
              index = ctxt->mask->cfgindex;

              if (levels == 1) {
                  AddLink2(ctxt, detect, index, ctxt->mask->score);
              }
              else                /* relance RecMaskSearch dans un voisinage */
              {                                           /* de la detection */
                  Context *ctxt1 = ctxt + 1;

                  GetBgn(ctxt, detect, &bgn);
                  StoreBgn(ctxt1, bgn);

                  ExportMaskCfgInfos(ctxt, ctxt->mask->gapscfg[index]);

                  DelGapCfgTable(ctxt1->mask);
                  DelCfgTable(ctxt1->mask);

                  ImportPatternCfgInfos(ctxt1);
                  GetMaskCfgs(ctxt1->mask);

                  (int) RecMaskSearch(ctxt1, levels - 1);

               /* si la detection a ete observee au 1er niveau, les tableaux */
                      /* de gaps sont reinitialises, des que la recherche au */
                                 /* voisinage de cette detection est achevee */

                  if (ctxt->level == 0)
                      ResetGapLists(ctxt, 0, levels - 1);
              }
          }
      }
      PrintStatus(ctxt, &jobscans);
  }

  return ctxt[levels - 1].detects;
}
/*=============================================================================
 Search(): Interface des fonctions de recherche 'MaskSearch' et 'RecMaskSearch'
           pilote par le choix de traitement 'maskproc'(option '-dmp' ou '-smp'
           de 'erpin'
=============================================================================*/

int Search(Context *ctxt, int nmask, int maskproc)
{
  int ndetects;

  if (maskproc == DYNAMIC)                /* recherche dynamique des masques */
  {                                              /* la plus souvent utilisee */
      (int) RecMaskSearch(ctxt, nmask);

      ndetects = fPrintfPostProcess2(ctxt + nmask - 1);
  }
  else if (maskproc == STATIC)             /* recherche statique des masques */
  {                                          /* convient a des motifs courts */
      (int) MaskSearch(ctxt, nmask);

      ndetects = fPrintfPostProcess(ctxt + nmask - 1);
  }
  else                                                             /* erreur */
  {
      fprintf(stderr, "Search: Invalid argument #3, exit..\n");
      exit(1);
  }

  return ndetects;
}
/*===========================================================================*/
