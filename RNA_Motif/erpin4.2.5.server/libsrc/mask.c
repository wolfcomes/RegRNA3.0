
/*=============================================================================
 mask.c                         A.Lambert le 26/11/01   revu le 05/05/02

 Code des fonctions destinees a gerer les masques de pattern pour le calcul des
 scores.

 Par "masque" on entend un ensemble d'helices et brins d'un meme 'pattern'.

 Un masque est entre par l'enumeration de symboles de ses elements.
 Ces arguments pourront etre entres de 3 facons (voir 'ReadMasksArgs'):

   -mask: les arguments entres sont exclus du traitement,
   -umask: les arguments entres seront les seuls a etre traites.
   -nomask: pas d'arguments, le masque coincide avec le pattern.

 Un symbole d'helice conduit a la selection d'une helice seulement si le pat-
 tern associe contient les 2 brins de celle-ci, sinon l'element est interprete
 comme un simple brin.


 cc -O2 -Wall -c mask.c -I ../include ;

 ar -rs ../lib/librnaIV.a mask.o ;

=============================================================================*/

#include "rnaIV.h"

Mask  *NewMask(void);
void  FreeMask(Mask *mask);
void  DelMask(Mask *mask);
void  SetUMask(Mask *mask, Pattern *pattern);
void  CtrlMaskErrors(Mask *mask, Pattern *pattern);
void  GetMaskGeom(Mask *mask);

void  GetNoMask(Mask *mask, Pattern *pattern);
void  GetMask(Mask *mask, Pattern *pattern);

/*=============================================================================
 NewMask(): Cree une structure 'Mask' et initialise ses champs pointeurs sus-
            ceptibles de se voir allouer de la memoire a NULL.
=============================================================================*/

Mask *NewMask(void)
{
  Mask *mask;

  mask = (Mask *) malloc(sizeof(Mask));
  mask->args = NULL;
  mask->hxindex = NULL;
  mask->stindex = NULL;
  mask->atomstr = NULL;
  mask->gapslist = NULL;
  mask->gapslist_ori = NULL;
  mask->gapscfg = NULL;
  mask->cfg = NULL;
  mask->str = NULL;

  return mask;  
}
/*=============================================================================
 FreeMask(): Libere la memoire allouee aux champs d'une structure Mask et re-
             initialise ses elements pointeurs a NULL;
             Si le masque coincide avec le pattern 'mask->cfg' et 'mask->str'
             ne sont pas alloues, ils pointent sur les elements du 'pattern',
             ceci pour economiser la memoire: ils sont exclus de l'operation.
             Cette eventualite est reperee par: mask->mode == NOMASK
=============================================================================*/

void FreeMask(Mask *mask)
{
  if (mask->args != NULL)    { free(mask->args); mask->args = NULL; }
  if (mask->hxindex != NULL) { free(mask->hxindex); mask->hxindex = NULL; }
  if (mask->stindex != NULL) { free(mask->stindex); mask->stindex = NULL; }
  if (mask->atomstr != NULL) { free(mask->atomstr); mask->atomstr = NULL; }
  if (mask->str != NULL)     { free(mask->str); mask->str = NULL; }

  if (mask->gapslist != NULL) {
      FreesMat(mask->gapslist);
      mask->gapslist = NULL;
  }
  if (mask->gapslist_ori != NULL) {
      FreesMat(mask->gapslist_ori);
      mask->gapslist_ori = NULL;
  }
  if (mask->gapscfg != NULL) DelGapCfgTable(mask);

  if (mask->cfg != NULL) DelCfgTable(mask);

  return;
}
/*=============================================================================
 DelMask(): detruit une structure Mask.
=============================================================================*/

void DelMask(Mask *mask)
{
  FreeMask(mask);
  free(mask);
  return;
}
/*=============================================================================
 SetUMask(): Initialise une structure 'Mask' a partir des arguments lus par 
             'ReadMasksArgs' et une structure pointee par 'pattern'.
             Le tableau pointe par 'mask->cfg' sera cree plus bas, ou sera aussi
             determine 'mask->ncg', dans 'GetMaskCfgs'.
             note1: 
             La taille des tableaux est augmentee d'une unite afin d'eviter le
             debordement lors de l'incrementation de pointeurs.
=============================================================================*/

void SetUMask(Mask *mask, Pattern *pattern)
{
  int  i, j, k;

  mask->pattern = pattern;            /* pointe le pattern associe au masque */

  for (i = k = 0; i < pattern->nhx; i++)     /* compte les indices d'helices */
  {
      for (j = 0; j < mask->nargs; j++)
          if (pattern->hlx[i].id == mask->args[j]) k++;
  }
  mask->nhx = k;

  for (i = k = 0; i < pattern->nst; i++)      /* compte des indices de brins */
  {
      for (j = 0; j < mask->nargs; j++)
          if (pattern->std[i].id == mask->args[j]) k++;
  }
  mask->nst = k;

  mask->hxindex = (int *) malloc((mask->nhx + 1) * sizeof(int));    /* note1 */
  mask->stindex = (int *) malloc((mask->nst + 1) * sizeof(int));

  for (i = k = 0; i < pattern->nhx; i++)     /* saisie des indices d'helices */
  {
      for (j = 0; j < mask->nargs; j++)
          if (pattern->hlx[i].id == mask->args[j])
              mask->hxindex[k++] = i;
  }

  for (i = k = 0; i < pattern->nst; i++)      /* saisie des indices de brins */
  {
      for (j = 0; j < mask->nargs; j++)
          if (pattern->std[i].id == mask->args[j])
              mask->stindex[k++] = i;
  }

  CtrlMaskErrors(mask, pattern);

  return;
}
/*=============================================================================
 CtrlMaskErrors(): Controle si tous les arguments saisis du masque pointe par
                   'mask' figurent dans ceux dans la structure pointee par
                   'pattern'.
                   En cas d'erreur, affiche un message et provoque la sortie.
=============================================================================*/

void CtrlMaskErrors(Mask *mask, Pattern *pattern)
{
  int  i, j, unknown, errors;

  errors = 0;
  for (i = 0; i < mask->nargs; i++)
  {
      unknown = 1;

      for (j = 0; j < pattern->nhx; j++)                 /* scan des helices */
          if (mask->args[i] == pattern->hlx[j].id) {  /* detection -> sortie */
              unknown = 0;
              break;
          }
      for (j = 0; j < pattern->nst; j++)                   /* scan des brins */
          if (mask->args[i] == pattern->std[j].id) {  /* detection -> sortie */
              unknown = 0;
              break;
          }
      if (unknown != 0) {
          errors++;
          fprintf(stderr,
          "CtrlMaskErrors: mask argument '%d' not found in pattern '%s'\n", 
          mask->args[i], pattern->id);
      }  
  }
  if (errors != 0) {
      fprintf(stderr,
          "CtrlMaskErrors: %d error%s detected, exit..\n",
          errors, (errors > 1 ? "s" : ""));
          exit(1);
  }
  return;
}
/*=============================================================================
 GetMaskGeom(): Initialise les champs 'db_bgn' et 'max_len' d'un masque, en
                utilisant les donnees des elements helices et brins du pattern
                associe: 'mask->pattern'
=============================================================================*/

void GetMaskGeom(Mask *mask)
{
  int  *m, i, db_bgn, db_end, end, first, last;

  db_bgn = last = mask->pattern->db_bgn + mask->pattern->max_len - 1;
  db_end = first = 0;

  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
  {
      if (mask->pattern->hlx[*m].db_bgn1 < db_bgn)
      { 
          db_bgn = mask->pattern->hlx[*m].db_bgn1;
          first  = mask->pattern->hlx[*m].min_bgn;
      }

      end = mask->pattern->hlx[*m].db_bgn1 + mask->pattern->hlx[*m].max_len - 1;

      if (end > db_end)
      {
          db_end = end;
          last   = mask->pattern->hlx[*m].min_bgn +
                   mask->pattern->hlx[*m].min_len;
      }
  }

  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
  {
      if (mask->pattern->std[*m].db_bgn < db_bgn)
      {
          db_bgn = mask->pattern->std[*m].db_bgn;
          first  = mask->pattern->std[*m].min_bgn;
      }

      end = mask->pattern->std[*m].db_bgn + mask->pattern->std[*m].max_len - 1;

      if (end > db_end)
      {
          db_end = end;
          last   = mask->pattern->std[*m].min_bgn +
                   mask->pattern->std[*m].min_len;
      }
  }
  mask->db_bgn = db_bgn;
  mask->max_len = db_end - db_bgn + 1;
  mask->min_len = last - first;
  mask->min_bgn = first;
  mask->max_bgn = mask->db_bgn - mask->pattern->db_bgn;

  return;
}
/*=============================================================================
 GetNoMask(): En retour 'mask' pointe sur un masque constitue de la totalite
              des elements du pattern pointe par l'argument, par simple copie
              des champs du masque depuis ceux de correspondant du 'pattern'
              associe.

              L'appel de cette fonction suit en general un appel a 'ReadMasks-
              Args' qui alloue la memoire et detecte le mode "NOMASK'.
=============================================================================*/

void GetNoMask(Mask *mask, Pattern *pattern)
{
  int  i;

  mask->mode = NOMASK;
  mask->nargs = pattern->nhx + pattern->nst;

  mask->nhx  = pattern->nhx;
  mask->nst  = pattern->nst;

  mask->args = (int *) malloc(mask->nargs * sizeof(int));
  mask->hxindex = (int *) malloc((mask->nhx + 1) * sizeof(int));
  mask->stindex = (int *) malloc((mask->nst + 1) * sizeof(int));

  for (i = 0; i < pattern->nhx; i++)
  {
      mask->args[i] = pattern->hlx[i].id;
      mask->hxindex[i] = i;
  }

  for (i = 0; i < pattern->nst; i++)
  {
      mask->args[pattern->nhx + i] = pattern->std[i].id;
      mask->stindex[i] = i;
  }

  mask->pattern = pattern;

  mask->db_bgn  = pattern->db_bgn;
  mask->min_bgn = 0;
  mask->max_bgn = 0;
  mask->max_len = pattern->max_len;
  mask->min_len = pattern->min_len;

  return;
}
/*=============================================================================
 GetMask(): Interface assurant la creation d'un masque et son initialisation
            complete, a l'exception du tableau de ses configurations, depuis
            l'examen du contenu d'un 'pattern'.
=============================================================================*/

void GetMask(Mask *mask, Pattern *pattern)
{

  switch (mask->mode)
  {
      case MASK:
      case UMASK:
      case ADDMASK: SetUMask(mask, pattern);
                    GetMaskGeom(mask);
                    break;
      case NOMASK:  GetNoMask(mask, pattern);
                    break;
      default:
          fprintf(stderr, "GetMask: incorrect field 'mode': %d, exit..\n",
                  mask->mode);
          exit(1);
  }

  mask->str = (char *) malloc(mask->pattern->max_len +
                              mask->pattern->natom + 1);
                              /* extention necessaire pour des "separateurs" */
  GetMaskAtoms(mask);                             /* fonction de 'maskcfg.c' */
  GetMaskAtomStr(mask);                                              /* idem */
  GetMaskGapList(mask);                                              /* idem */

  return;
}
/*===========================================================================*/
