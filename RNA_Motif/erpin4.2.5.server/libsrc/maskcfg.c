
 
/*=============================================================================
 maskcfg.c                       A.Lambert le 28/01/02    revu le 10/04/02

 Ce code est destine au traitement des configurations de motifs (masques) d'ARN.
 On gere de facon differente les atomes qui appartiennent au masque et ceux qui
 ne lui appartiennent pas, afin de calculer de facon exacte les differentes con-
 figurations du masque etudie;
 2 atomes consecutifs n'appartenant pas au masque, contenant respectivement g1
 et g2 gaps, donnent lieu a (g1 + g2 + 1) combinaisons alors qu'ils intervien-
 nent pour (g1 + 1)x(g2 + 1) s'ils appartiennent au masque.

 La chaine 'mask->atomstr' code les atomes de l'enveloppe du masque 'mask->atom'
 par '1' ou '0' suivant qu'ils appartiennent ou non au masque.
 Le tableau 'mask->gapslist' enregistre les gaps dans le meme ordre.

 cc -Wall -O2 -c maskcfg.c -I../include  -DDEBUG ;

 cc -Wall -O2 -c maskcfg.c -I../include ;
 ar -rs ../lib/librnaIV.a maskcfg.o ;

=============================================================================*/

#include "rnaIV.h"

void  GetMaskAtoms(Mask *mask);
void  GetMaskAtomStr(Mask *mask);
void  GetMaskGapList(Mask *mask);
void  ResetMaskGapList(Mask *mask);

void  EnumGapCfgs(short **mgaplist, short nvarst, short col, int *nb,
                  short *tmpgaps, short **gapscfglist);
int   GetGapsCfgs(Mask *mask);
void  PrtGapsCfg(Mask *mask, int cfgindex, FILE *txt);
void  PrtGapsCfgs(Mask *mask, FILE *txt);

void  GetCfg(Config *cfg, Mask *mask, short *gaps_cfg);
void  GetMaskCfgs(Mask *mask);

int   GetMaskCfgLen(Config *cfg, Mask *mask);
void  GetMaskCfgsLens(Mask *mask);

/*=============================================================================
 GetMaskAtoms(): Determine dans le 'pattern' associe a un masque le 1er atome
                 et le nombre d'atomes correspondant a la longueur du masque.
                 Initialise les champs 'mask->atom' et 'mask->natom' de 'mask'.

                 C'est la determination des atomes de l'ENVELOPPE du masque.

                 Cette fonction permettra de calculer les configurations d'un
                 masque, en excluant le debut et la fin du 'pattern' qui ne
                 pas partie du masque.
=============================================================================*/

void GetMaskAtoms(Mask *mask)
{
  int  *m, i, min_pos1, max_pos2, atom1 = 0,
       id1, id2;                      /* identificateurs des atomes extremes */

  min_pos1 = 100000;
  max_pos2 = 0;
  id1 = id2 = 0;

              /* recherche des identificateurs des atomes extremes du masque */

  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
  {
      if (mask->pattern->hlx[*m].db_bgn1 < min_pos1) {
          min_pos1 = mask->pattern->hlx[*m].db_bgn1;
          id1 = mask->pattern->hlx[*m].id;
      }
      if (mask->pattern->hlx[*m].db_bgn2 > max_pos2) {
          max_pos2 = mask->pattern->hlx[*m].db_bgn2;
          id2 = mask->pattern->hlx[*m].id;
      }

  }
  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
  { 
      if (mask->pattern->std[*m].db_bgn < min_pos1) {
          min_pos1 = mask->pattern->std[*m].db_bgn;
          id1 = mask->pattern->std[*m].id;
      }
      if (mask->pattern->std[*m].db_bgn > max_pos2) {
          max_pos2 = mask->pattern->std[*m].db_bgn;
          id2 = mask->pattern->std[*m].id;
      }
  }

     /* passage des identificateurs aux indices des atomes dans le 'pattern' */

  for (i = 0; i < mask->pattern->natom; i++)
  {
      if (id1 == mask->pattern->oatlist[i]) {     /* 1ere occurence de 'id1' */
          atom1 = i;
          mask->atom = mask->pattern->atom + atom1;
          break;
      }
  }
  for (i = mask->pattern->natom - 1; i >= 0; i--)
  {
      if (id2 == mask->pattern->oatlist[i]) { /* derniere occurence de 'id2' */
          mask->natom = i - atom1 + 1;
          break;
      }
  }

  return;
}
/*=============================================================================
 GetMaskAtomStr(): Cree une chaine de longueur 'mask->natom' constituee des
                   caracteres '1' et '0' suivant que l'atome correspondant ap-
                   partient ou non au masque pointe par 'mask'.
                   Cette chaine sera utilisee pour le calcul des configurations
                   d'un masque.
=============================================================================*/

void GetMaskAtomStr(Mask *mask)
{
  int  i, j, offset, *m;

  mask->atomstr = (char *) malloc(mask->natom + 1);

  offset = mask->atom - mask->pattern->atom;

  for (i = 0; i < mask->natom; i++)
  {
      switch (mask->atom[i].type) {
          case HLX1: 
          case HLX2:
              mask->atomstr[i] = '0'; 
              for (j = 0, m = mask->hxindex; j < mask->nhx; j++, m++)
                  if (mask->pattern->hlx[*m].id == mask->pattern->oatlist[offset+i])
                  {
                      mask->atomstr[i] = '1';   /* l'atome appart. au masque */
                      break;
                  }
              break;
          case STD:
             mask->atomstr[i] = '0'; 
             for (j = 0, m = mask->stindex; j < mask->nst; j++, m++)
                  if (mask->pattern->std[*m].id == mask->pattern->oatlist[offset+i])
                  {
                      mask->atomstr[i] = '1';   /* l'atome appart. au masque */
                      break;
                  }
              break;
           default: break;
      }
  }
  mask->atomstr[mask->natom] = '\0';

  return;
}
/*=============================================================================
 GetMaskGapList(): Retourne un tableau contenant les nombres extremes de gaps
                   des atomes successifs de 'mask' qui en contiennent et, dans
                   l'ordre des atomes, pour les atomes n'appartenant pas au
                   masque:
                   - si le groupe d'atomes ne contient qu'un brin avec gaps,
                   parmi d'autres sans gaps, celui-ci est traite comme les brins
                   du masque (l'indice du brin est enregistre en 1ere ligne).
                   - si le groupe d'atomes contient plusieurs brins avec gaps,
                   ceux-ci sont cumules en un seul element afin de ne pas creer
                   des repetitions dans les configurations du masque.
                   La 1ere ligne du tableau retourne contient l'indice du brin
                   associe dans les premiers cas et -1 dans le dernier.

                   Le nombre de configurations est enregistre dans 'mask->ncfg'.
                   'mask->ngaps' donne le nombre d'elements du tableau des gaps.
                   Si aucun brin ne contient de gaps 'mask->ngaps' est mis a 0
                   et 'mask->gapslist' a 'NULL'.
=============================================================================*/

void GetMaskGapList(Mask *mask)
{
  int    i, j, g, sum, nvarst;
  short  tmpindex = 0, **tmpgapslist;

  for (i = j = 0; i < mask->natom; i++)
      if (mask->atom[i].max_gaps > 0) j++;    /* compte les atomes avec gaps */

  if (j == 0) {
      mask->ngaps = 0;                /* initialisation a 0 de 'mask->ngaps' */
      mask->gapslist = NULL;
      mask->ncfg = 1;

      return;
  }

  tmpgapslist = sMat(2, j);
  FillsMat(tmpgapslist, 2, j, 0);                             /* mise a zero */
  mask->ngaps = 0;                    /* initialisation a 0 de 'mask->ngaps' */

  for (i = j = 0; i < mask->natom; i++)
  {
      if (mask->atom[i].type == STD)
      {
          switch (mask->atomstr[i])
          {
              case '1':                       /* atome appartenant au masque */
                  if (mask->atom[i].max_gaps > 0)
                  {
                      tmpgapslist[0][j] = (short) mask->atom[i].index;
                      tmpgapslist[1][j] = (short) mask->atom[i].max_gaps;
                      j++;
                  }
                  break;

              case '0':                 /* atome n'appartenant pas au masque */
                  nvarst = 0;
                  sum = 0;
                  while (mask->atomstr[i] == '0')
                  {
                      if (mask->atom[i].type == STD)
                      {
                          g = mask->atom[i].max_gaps;
                          if (g > 0) {
                              tmpindex = (short) mask->atom[i].index;
                              sum += g;                             /* cumul */
                              nvarst++;             /* nb de brins avec gaps */
                          }
                      }
                      i++;
                  }
                  if (sum > 0)
                  {
                      if (nvarst == 1)              /* 1 seul brin avec gaps */
                          tmpgapslist[0][j] = tmpindex;
                      else 
                          tmpgapslist[0][j] = -1;    /* pas d'indice de brin */
                                                               /* enregistre */
                      tmpgapslist[1][j] = (short) sum;
                      j++;
                  }
                  i--;
                  break;
          }
      }
  }
  mask->ngaps = j;    /* reinitialisation de 'mask->ngaps' (<= au precedent) */

  mask->gapslist_ori = sMat(3, mask->ngaps);
  mask->gapslist     = sMat(3, mask->ngaps);

  for (j = 0; j < mask->ngaps; j++) {
      mask->gapslist[0][j] = mask->gapslist_ori[0][j] = tmpgapslist[0][j];
      mask->gapslist[1][j] = mask->gapslist_ori[1][j] = 0;
      mask->gapslist[2][j] = mask->gapslist_ori[2][j] = tmpgapslist[1][j];
  }

  FreesMat(tmpgapslist);

  GetMaskCfgsNb(mask);       /* compte le nombre de configurations possibles */

  return;
}
/*=============================================================================
 ResetMaskGapList(): Remet les elements du champ 'gaplist' de 'mask' a leur
                     valeur originale.
=============================================================================*/

void ResetMaskGapList(Mask *mask)
{
  CopysMat(mask->gapslist_ori, 3, mask->ngaps, mask->gapslist);
  return;
}
/*=============================================================================
         ICI COMMENCE LA CONSTRUCTION DES CONFIGURATIONS D'UN MASQUE.
=============================================================================*/

/*=============================================================================
 EnumGapCfgs(): Enumere recursivement les configurations de gaps construites
                depuis le tableau 'mgaplist' de longueur 'nvarst', 'col' est
                l'indice de la colonne No 'col' dont les elements prennent les
                valeurs de mgaplist[1][col] a mgaplist[2][col].
                Le tableau 'gaps_cfg_list' est recursivement rempli des confi-
                gurations visitees.
                'nb' pointe le nombre de configurations trouvees.
                'tmpgaps' est un tableau temporaire ou s'inscrit une configu-
                ration.
=============================================================================*/

void EnumGapCfgs(short **mgaplist, short nvarst, short col, int *nb,
                 short *tmpgaps, short **gapscfglist)
{
  short  i, j;

  if (nvarst == 0) {
      gapscfglist[0][0] = 0;
      *nb = 1;
      return;
  }
  for (i = mgaplist[1][col]; i <= mgaplist[2][col]; i++)
  {
      tmpgaps[col] = i;
      if (col == nvarst - 1)                    /* derniere colonne atteinte */
      {
          for (j = 0; j < nvarst; j++) 
              gapscfglist[*nb][j] = tmpgaps[j];
          (*nb)++;
      }
      else
      EnumGapCfgs(mgaplist, nvarst, col + 1, nb, tmpgaps, gapscfglist);
  }
  return;
}
/*=============================================================================
 GetGapsCfgs(): Interface de la fonction recursive 'EnumGapCfgs'.
                En retour le tableau des configurations possibles concernant
                les brins variables de la structure secondaire pointee par
                'mask' est charge a l'adresse de 'mask->gapscfg'.
                Retourne, pour controle, le nombre de configurations.
=============================================================================*/

int GetGapsCfgs(Mask *mask)
{
  short  *tmpgaps;
  int    nb;

  tmpgaps = (short *) calloc(mask->ngaps, sizeof(short));
  SetGapCfgTable(mask);

  nb = 0;
  EnumGapCfgs(mask->gapslist, mask->ngaps, 0, &nb, tmpgaps, mask->gapscfg);

  free(tmpgaps);

#ifdef DEBUG
printf("\nconfigs: %d expected, %d found\n\n", mask->ncfg, nb);
#endif

  return nb;
}
/*=============================================================================
 PrtGapsCfg(): Imprime le contenu de la configuration de gaps No 'cfgindex'
               du masque pointe par 'mask' dans le fichier pointe par 'txt'.
=============================================================================*/

void PrtGapsCfg(Mask *mask, int cfgindex, FILE *txt)
{
  int j;

  for (j = 0; j < mask->ngaps; j++)
      fprintf(txt, "%-2d ", mask->gapscfg[cfgindex][j]);
  fprintf(txt, "\n");

  return;
}
/*=============================================================================
 PrtGapsCfgs(): Imprime le contenu du tableau des configurations de gaps du
                masque pointe par 'mask' dans le fichier pointe par 'txt'.
=============================================================================*/

void PrtGapsCfgs(Mask *mask, FILE *txt)
{
  int i;

  for (i = 0; i < mask->ncfg; i++)  PrtGapsCfg(mask, i, txt);
  fprintf(txt, "\n");

  return;
}
/*=============================================================================
 GetCfg(): Retourne, pointee par 'cfg', la configuration du masque pointe par
           'mask' calculee a l'aide des parametres lus dans la configuration
           de gaps pointee par 'gaps_cfg'.
           Le debut de la configuration a pour abscisse zero.

           La configuration, ainsi que 'gaps_cfg', prennent en compte de facon
           appropriee tous les atomes de l'enveloppe du masque, meme ceux qui
           ne lui appartiennent pas.

   rappel: atom[i].h2index vaut 0 ou, dans le cas ou atom[i].type = HLX1 (1er
           brin d'helice), a pour valeur l'indice de l'atome du 2eme brin de
           l'helice correspondante (a ete initialise dans 'ReadHelix').

=============================================================================*/

void GetCfg(Config *cfg, Mask *mask, short *gaps_cfg)
{
  int  i, j, h2indx, dist, gapsum, min_dist, offset,

       stindx, hxindx;        /* indices dans les tab. de brins et d'helices */

  offset = mask->atom - mask->pattern->atom;


             /* --- REINITIALISATION DES PARAMETRES VARIABLES DES ATOMES --- */


  mask->atom[0].len = mask->atom[0].min_len;          /* 1er atome du masque */
  mask->atom[0].gaps = 0;
  mask->atom[0].bgn  = 0;
 
  for (i = 1; i < mask->natom; i++)                   /* les atomes suivants */
  {
      mask->atom[i].gaps = 0;
      mask->atom[i].len = mask->atom[i].min_len;
      mask->atom[i].bgn = mask->atom[i-1].bgn + mask->atom[i-1].len;
  }


 /* -- ACTUALISATION DES ATOMES DU MASQUE A L'AIDE DE LA CONFIG. DES GAPS -- */


  for (i = j = gapsum = 0; i < mask->natom; i++)
  {
      mask->atom[i].bgn += gapsum;
                                          /* restriction aux brins avec gaps */
      if (mask->atom[i].type == STD && mask->atom[i].max_gaps > 0)
      {
          switch (mask->atomstr[i]) {
              case '1':                           /* atome appart. au masque */
                      mask->atom[i].gaps = gaps_cfg[j];
                      mask->atom[i].len += gaps_cfg[j];
                      gapsum += gaps_cfg[j++];
                  break;
              case '0':                                 /* atome hors masque */
                      while (mask->atomstr[i] == '0') i++;          /* passe */
                      gapsum += gaps_cfg[j++];           /* les brins exclus */
                      i--;                     /* et s'arrete sur le dernier */
                  break;
          }
      }
  }

  cfg->len = (short) (mask->min_len + gapsum);     /* longueur de la config. */


     /* --- ACTUALISATION DE LA GEOMETRIE DES HELICES ET BRINS DU MASQUE --- */


  for (i = stindx = hxindx = 0; i < mask->natom; i++)
  {
      if (mask->atomstr[i] == '1')       /* restriction aux atomes du masque */
      {
          switch (mask->atom[i].type)
          {
              case STD:                                              /* brin */
                  cfg->stbgn[stindx]  = (short) mask->atom[i].bgn;
                  cfg->stgaps[stindx] = (short) mask->atom[i].gaps;
                  cfg->stlen[stindx]  = (short) mask->atom[i].len;
                  stindx++;
                  break;

              case HLX1:                                /* 1er brin d'helice */
                  cfg->hxbgn[hxindx] = (short) mask->atom[i].bgn;

                  h2indx = mask->atom[i].h2index - offset;      /* 2eme brin */

                  dist = mask->atom[h2indx].bgn - mask->atom[i + 1].bgn;
                  min_dist = mask->atom[h2indx].min_bgn -
                             mask->atom[i + 1].min_bgn;

                  cfg->hxdist[hxindx] = (short) dist;
                  cfg->hxgaps[hxindx] = (short) (dist - min_dist);
                  hxindx++;
                  break;
              case HLX2:
                  break;
              default:
                  fprintf(stderr,
                          "GetCfg: invalid 'atom.type' for arg 2, exit..\n");
                  exit(1);
          }
      }

  }
  return;
}
/*=============================================================================
 GetMaskCfgs(): Dresse le tableau des configurations possibles du masque pointe
                par 'mask' a partir des configurations de gaps trouvees par la
                fonction 'GetGapsCfgs' qui est appelee.
                Le masque doit avoir ete completement initialise au prealable.
                Le tableau retourne doit posseder 'mask->ncfg' elements, un
                controle est effectue.
=============================================================================*/

void GetMaskCfgs(Mask *mask)
{
  int     i, cfgs;

  if (mask->nhx + mask->nst == 0) {
      mask->ncfg = 0;
      return;
  }
  if (mask->cfg != NULL)  DelCfgTable(mask);

  SetCfgTable(mask);               /* creation du tableau des configurations */

  cfgs = GetGapsCfgs(mask);                        /* configurations de gaps */

  if (cfgs != mask->ncfg) {                                      /* controle */
      fprintf(stderr, "GetMaskCfgs: %d config. expected, %d found! exit..\n",
              mask->ncfg, cfgs);
      exit(1);
  }

  for (i = 0; i < mask->ncfg; i++)               /* configurations du masque */
      GetCfg(mask->cfg + i, mask, mask->gapscfg[i]);

  return;
}
/*=============================================================================
 GetMaskCfgLen(): Depuis les donnees d'une configuration pointee par 'cfg' pi-
                  lotee par un masque pointe par 'mask' retourne la longueur
                  de cette configuration.
            note: Les longueurs des configurations ont ete determinees lors de
                  leur creation.
=============================================================================*/

int GetMaskCfgLen(Config *cfg, Mask *mask)
{
  int   *m, i, pos2, min_pos1, max_pos2;

  min_pos1 = 100000;
  max_pos2 = 0;

  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
  {
      if (cfg->hxbgn[i] < min_pos1)  min_pos1 = cfg->hxbgn[i];

      pos2 = cfg->hxbgn[i] + 2 * mask->pattern->hlx[*m].helix_len +
             cfg->hxdist[i];

      if (pos2 > max_pos2)  max_pos2 = pos2;
  }
  for (i = 0; i < mask->nst; i++)
  {
      if (cfg->stbgn[i] < min_pos1)  min_pos1 = cfg->stbgn[i];
 
      pos2 = cfg->stbgn[i] + cfg->stlen[i];

      if (pos2 > max_pos2)  max_pos2 = pos2;
  }

  return (max_pos2 - min_pos1);
}
/*=============================================================================
 GetMaskCfgsLens(): Applique 'GetMaskCfgLen' a l'ensemble des configurations
                    associees a un masque, les longueurs calculees sont chargees
                    dans le champ 'len' des configurations.
=============================================================================*/

void GetMaskCfgsLens(Mask *mask)
{
  int  i;

  for (i = 0; i < mask->ncfg; i++)
      mask->cfg[i].len = (short) GetMaskCfgLen(mask->cfg + i, mask);

  return;
}
/*===========================================================================*/
