
/*=============================================================================
 cfg.c                         A.Lambert  le 26/09/01  revu le 24/04/02

 Soit un masque constitue de 'nhx' helices et 'nst' brins dont certains ont
 des longueurs variables, lues dans un 'trset';
 Une configuration de ce masque est donnee par un positionnement particulier
 de l'ensemble de ses elements dans les limites imposees par les nombres de gaps
 possibles dans les brins variables.
 Une configuration est identifiee de facon univoque par la donnee pour tous les
 constituants des elements suivant:
 Helice: absc. du debut et dist. entre les brins, soit 'bgn' et 'dist';
 Brin:   absc. du debut et longueur du brin, soit 'bgn' et 'len'.

 La donnee d'une configuration doit permettre le calcul immediat d'un score de
 sequence, et, les tableaux de scores des elements etant calcules, de permettre
 de trouver par simple enumeration celle qui a le score maximal.

 On etudie ici l'ensemble des configurations possibles d'un masque dont le
 nombre total est le produit des (max_gaps + 1) effectue sur l'ensemble des
 brins variables, soit 'pattern->ncfg'.
 On utilisera pour les calculs le tableau 'atom' d'un 'pattern' auquel le masque
 est associe.
 
 Une configuration sera stockee dans une structure 'Config'.

 cc -O2 -Wall -c cfg.c -I../include ;

 ar -rs ../lib/librnaIV.a cfg.o ;

=============================================================================*/

#include "rnaIV.h"

void SetCfgTable(Mask *mask);
void FreeCfgFields(Config *cfg);
void DelCfgTable(Mask *mask);
void SetGapCfgTable(Mask *mask);
void DelGapCfgTable(Mask *mask);

Config *CloneCfg(Mask *mask, Config *cfg);

/*=============================================================================
 SetCfgTable(): Cree un tableau de configurations associees au masque pointe 
                par 'mask' qui doit etre correctement initialise.
=============================================================================*/

void SetCfgTable(Mask *mask)
{
  int i;

  mask->cfg = (Config *) malloc(mask->ncfg * sizeof(Config));

  if (mask->cfg == NULL) {
      fprintf(stderr, "SetCfgTable: allocation failure, exit..\n");
      exit(1);
  }

  for (i = 0; i < mask->ncfg; i++)
  {
      mask->cfg[i].hxbgn  = (short *) calloc(mask->nhx + 1, sizeof(short));
      mask->cfg[i].hxgaps = (short *) calloc(mask->nhx + 1, sizeof(short));
      mask->cfg[i].hxdist = (short *) calloc(mask->nhx + 1, sizeof(short));
      mask->cfg[i].stbgn  = (short *) calloc(mask->nst + 1, sizeof(short));
      mask->cfg[i].stgaps = (short *) calloc(mask->nst + 1, sizeof(short));
      mask->cfg[i].stlen  = (short *) calloc(mask->nst + 1, sizeof(short));
  }
  return;
}
/*=============================================================================
 FreeCfgFields(): Libere la memoire allouee aux champs de la configuration
                  pointee par 'cfg';
=============================================================================*/

void FreeCfgFields(Config *cfg)
{
  free(cfg->hxbgn);
  free(cfg->hxgaps);
  free(cfg->hxdist);
  free(cfg->stbgn);
  free(cfg->stgaps);
  free(cfg->stlen); 

  return;
}
/*=============================================================================
 DelCfgTable(): Detruit un tableau de structures de type 'Config' associe a la
                structure pointee par 'mask'.
=============================================================================*/

void DelCfgTable(Mask *mask)
{
  int i;

  if (mask->cfg == NULL) return;

  for (i = 0; i < mask->ncfg; i++)   FreeCfgFields(mask->cfg + i);

  free(mask->cfg);
  mask->cfg = NULL;

  return;
}
/*=============================================================================
 SetGapCfgTable(): Cree le tableau des configurations de gaps d'un masque.
=============================================================================*/

void SetGapCfgTable(Mask *mask)
{
  mask->gapscfg = sMat(mask->ncfg, mask->ngaps);

  return;
}
/*=============================================================================
 DelGapCfgTable(): Detruit le tableau des configurations de gaps d'un masque.
=============================================================================*/

void DelGapCfgTable(Mask *mask)
{
  if (mask->gapscfg == NULL) return;

  FreesMat(mask->gapscfg);
  mask->gapscfg = NULL;

  return;
}
/*=============================================================================
 CloneCfg(): Cree une copie de la configuration du masque pointe par 'mask'
             pointee par 'cfg'.
             Retourne un pointeur sur cette copie.
=============================================================================*/

Config *CloneCfg(Mask *mask, Config *cfg)
{
  int    i;
  Config *copy;

  copy = (Config *) malloc(sizeof(Config));

  copy->hxbgn  = (short *) calloc(mask->nhx + 1, sizeof(short));
  copy->hxgaps = (short *) calloc(mask->nhx + 1, sizeof(short));
  copy->hxdist = (short *) calloc(mask->nhx + 1, sizeof(short));
  copy->stbgn  = (short *) calloc(mask->nst + 1, sizeof(short));
  copy->stgaps = (short *) calloc(mask->nst + 1, sizeof(short));
  copy->stlen  = (short *) calloc(mask->nst + 1, sizeof(short));

  for (i = 0; i < mask->nhx; i++)
  {
      copy->hxbgn[i]  = cfg->hxbgn[i];
      copy->hxgaps[i] = cfg->hxgaps[i];
      copy->hxdist[i] = cfg->hxdist[i];
  }

  for (i = 0; i < mask->nst; i++)
  {
      copy->stbgn[i]  = cfg->stbgn[i];
      copy->stgaps[i] = cfg->stgaps[i];
      copy->stlen[i]  = cfg->stlen[i];
  }
  copy->len = cfg->len;

  return copy;
}
/*===========================================================================*/
