
 
/*=============================================================================
 list.c                          A.Lambert le 03/01/02   revu le 17/04/02

 Creation d'une liste chainee, organisee en 'pile', destinee a enregistrer les
 donnees, en nombre aleatoire, lors d'une operation de recherche dans une sequ-
 ence, en vue d'un traitement global des resultats.

 On opere de 2 facons differentes suivant que les configurations des motifs
 sont gerees dynamiquement ou non.
 On utilise le champ 'cfg' ou l'indice 'cfgindex' suivant la gestion adoptee.

 Cette liste est associee a un 'Context': champ 'ctxt->list'.

 cc -O2 -Wall -c list.c -I../include ;

 ar -rs ../lib/librnaIV.a list.o ;

=============================================================================*/

#include "rnaIV.h"


void AddLink(Context *ctxt, int offset, int cfgindex, double score);
void AddLink2(Context *ctxt, int offset, int cfgindex, double score);
void FreeList(Context *ctxt);
Data *CopyList(Context *ctxt, int *links);
void OverlapFilter(Context *ctxt);
void OverlapFilter2(Context *ctxt);
int  fPrintfPostProcess(Context *ctxt);
int  fPrintfPostProcess2(Context *ctxt);

/*=============================================================================
 AddLink(): Ajoute un element a la liste pointee par 'toplist' et initialise
            ses champs de donnees avec ses arguments.
            Retourne un pointeur sur le nouvel element de la pile.
=============================================================================*/

void AddLink(Context *ctxt, int offset, int cfgindex, double score)
{
  List *newlink;

  newlink = (List *) malloc(sizeof(List));

  newlink->data.offset = offset;
  newlink->data.cfgindex = cfgindex;
  newlink->data.score = score;
  newlink->data.cfg = NULL;

  newlink->next = ctxt->list;
  ctxt->list = newlink;

  return;
}
/*=============================================================================
 AddLink2(): Ajoute un element a la liste pointee par 'toplist' et initialise
             ses champs de donnees avec ses arguments.
             Le champ 'data.cfg' est cree et initialise mais 'data.cfgindex'
             est mis a 0;
             Retourne un pointeur sur le nouvel element de la pile.
=============================================================================*/

void AddLink2(Context *ctxt, int offset, int cfgindex, double score)
{
  List *newlink;

  newlink = (List *) malloc(sizeof(List));

  newlink->data.offset = offset;
  newlink->data.cfgindex = 0;
  newlink->data.score = score;
  newlink->data.cfg = CloneCfg(ctxt->mask, ctxt->mask->cfg + cfgindex);

  newlink->next = ctxt->list;
  ctxt->list = newlink;

  return;
}
/*=============================================================================
 FreeList(): Libere la memoire allouee a l'ensemble des elements de la liste
             pointee par 'ctxt->list', et la reinitialise a 'NULL'.
=============================================================================*/

void FreeList(Context *ctxt)
{
  List  *tmp;

  while (ctxt->list != NULL)
  {
      tmp = ctxt->list->next;

      if (ctxt->list->data.cfg != NULL) {
          FreeCfgFields(ctxt->list->data.cfg);
          free(ctxt->list->data.cfg);
      }
      free((char *) ctxt->list);
      ctxt->list = tmp;
  }
  ctxt->list = NULL;

  return;
}
/*=============================================================================
 CopyList(): Copie les donnees contenues dans les elements de la liste pointee
             par 'ctxt->list', dans un tableau.
             si 'ctxt->dir == FORWARD' l'ordre de la pile source est inverse.
             Retourne un pointeur sur le tableau cree.
=============================================================================*/

Data *CopyList(Context *ctxt, int *links)
{
  int   i;
  Data  *data;
  List  *list = ctxt->list;

  *links = 0;
  while (list != NULL) {
      (*links)++;
      list = list->next;
  }
  data = (Data *) malloc(*links * sizeof(Data));

  if (ctxt->dir == FORWARD)
  {
      i = *links - 1;
      while (ctxt->list != NULL)
      {
          data[i--] = ctxt->list->data;
          ctxt->list = ctxt->list->next;
      }
  }
  else
  {
      i = 0;
      while (ctxt->list != NULL)
      {
          data[i++] = ctxt->list->data;
          ctxt->list = ctxt->list->next;
      }
  }

  return data;
}
/*=============================================================================
 OverlapFilter(): Filtre, pour le 'Context' pointe par 'ctxt', les detections
                  dont les sequences se chevauchent pour ne retenir que celle
                  de plus fort score.
=============================================================================*/

void OverlapFilter(Context *ctxt)
{
  int    i, j, index, previous, links;
  double scoremax;
  Data   *data;

  ctxt->detects = 0;                                     /* reinitialisation */
  data = CopyList(ctxt, &links);

  if (ctxt->dir == FORWARD)                                   /* brin direct */
  {
      for (i = 0; i < links; i++)
      {
          scoremax = data[i].score;
          index = i;
          previous = data[i].offset + ctxt->mask->cfg[data[i].cfgindex].len;

          for (j = i+1; j < links && data[j].offset < previous; j++)
              if (data[j].score > scoremax)
              {
                  scoremax = data[j].score;
                  index = j;
              }
          ctxt->detects++;
          fPrintfData(ctxt, data + index);
          i = --j;
      }
  }
  else                           /* dans ce cas on a 'ctxt->dir == REV_CMPL' */
  {
      for (i = 0; i < links; i++)
      {
          scoremax = data[i].score;
          index = i;

          for (j = i+1; j < links &&
              data[j].offset +
              ctxt->mask->cfg[data[j].cfgindex].len > data[i].offset;
              j++)
              if (data[j].score > scoremax)
              {
                  scoremax = data[j].score;
                  index = j;
              }
          ctxt->detects++;
          fPrintfData(ctxt, data + index);
          i = --j;
      }
  }
  free(data);

  return;
}
/*=============================================================================
 OverlapFilter2(): Filtre, pour le 'Context' pointe par 'ctxt', les detections
                   dont les sequences se chevauchent pour ne retenir que celle
                   de plus fort score.
                   Version adaptee a la gestion dynamique des configurations.
=============================================================================*/

void OverlapFilter2(Context *ctxt)
{
  int    i, j, index, previous, links;
  double scoremax;
  Data   *data;

  ctxt->detects = 0;                                     /* reinitialisation */
  data = CopyList(ctxt, &links);

  if (ctxt->dir == FORWARD)                                   /* brin direct */
  {
      for (i = 0; i < links; i++)
      {
          scoremax = data[i].score;
          index = i;
          previous = data[i].offset + data[i].cfg->len;

          for (j = i+1; j < links && data[j].offset < previous; j++)
              if (data[j].score > scoremax)
              {
                  scoremax = data[j].score;
                  index = j;
              }
          ctxt->detects++;
          fPrintfData(ctxt, data + index);
          i = --j;
      }
  }
  else                           /* dans ce cas on a 'ctxt->dir == REV_CMPL' */
  {
      for (i = 0; i < links; i++)
      {
          scoremax = data[i].score;
          index = i;

          for (j = i+1; j < links &&
              data[j].offset + data[j].cfg->len > data[i].offset; j++)
              if (data[j].score > scoremax)
              {
                  scoremax = data[j].score;
                  index = j;
              }
          ctxt->detects++;
          fPrintfData(ctxt, data + index);
          i = --j;
      }
  }
  free(data);

  return;
}
/*=============================================================================
 fPrintfPostProcess(): Affiche les resultats traites.
                       Retourne le nombre de motifs trouves dont sont exclues
                       les repetitions et 'overlap'.
=============================================================================*/

int fPrintfPostProcess(Context *ctxt)
{
  ClearStatus();               /* Efface l'etat d'avancement de la recherche */
  OverlapFilter(ctxt);
  FreeList(ctxt);

  return ctxt->detects;
}
/*=============================================================================
 fPrintfPostProcess2(): Affiche les resultats traites.
                        Retourne le nombre de motifs trouves dont sont exclues
                        les repetitions et 'overlap'.
=============================================================================*/

int fPrintfPostProcess2(Context *ctxt)
{
  ClearStatus();               /* Efface l'etat d'avancement de la recherche */
  OverlapFilter2(ctxt);
  FreeList(ctxt);

  return ctxt->detects;
}
/*===========================================================================*/
