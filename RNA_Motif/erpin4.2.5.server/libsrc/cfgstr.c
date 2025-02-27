
/*=============================================================================
 cfgstr.c                                   A.Lambert  revu le 11/06/02

 fonctions destinees a reconstituer une sequence "alignee" sur la base d'entrai-
 nement associee, a partir de la configuration de score optimal enregistree dans
 un masque.

 cc -O2 -Wall -c cfgstr.c -I../include ;

 ar -rs ../lib/librnaIV.a cfgstr.o ;

=============================================================================*/
#include "rnaIV.h"

void RecordMaskSeq(Mask *mask, char *seq, Config *cfg);
void CharInsert(char *str, int pos, char insert);
void ToLower(char *str, int len);
void GetCfgAtoms(Config *cfg, Mask *mask);

/*=============================================================================
 RecordMaskSeq(): Depuis les donnees de la configuration pointee par 'cfg'
                  d'une structure pointee par 'mask', cette fonction enregis-
                  tre a partir de l'abscisse d'une sequence pointee par 'seq'
                  le contenu du brin correspondant.
                  Celui-ci est en retour pointe par 'pattern->str'.
                  Variante de 'RecordMaskSeq'.
                  ATTENTION: Les arguments 'mask' et 'cfg' doivent etre com-
                  patibles, aucun controle n'est effectue.
                  Le contenu des brins non traites car n'appartenant pas au
                  masque est egalement affiche, en lettres minuscules, sans
                  donner lieu a une operation d'alignement.
=============================================================================*/

void RecordMaskSeq(Mask *mask, char *seq, Config *cfg)
{
  int     *m, i, size;
  char    *dest, *src;

  FillStr(mask->str, mask->max_len + mask->natom, '-');

  for (i = 0, m = mask->hxindex; i < mask->nhx; i++, m++)
  {
      dest = mask->str + mask->pattern->hlx[*m].db_bgn1 - mask->db_bgn;
      src  = seq + (int) cfg->hxbgn[i];
      size = mask->pattern->hlx[*m].helix_len;

      strncpy(dest, src, size);                      /* 1er brin de l'helice */

      dest = mask->str + mask->pattern->hlx[*m].db_bgn2 - mask->db_bgn;
      src += size + (int) cfg->hxdist[i];

      strncpy(dest, src, size);                     /* 2eme brin de l'helice */
  }

  for (i = 0, m = mask->stindex; i < mask->nst; i++, m++)
  {
      dest = mask->str + mask->pattern->std[*m].db_bgn - mask->db_bgn;
      src  = seq + (int) cfg->stbgn[i];
      size = (int) cfg->stlen[i];

      if (mask->pattern->std[*m].max_gaps == 0)      /* brin sans gap: copie */
      {
          strncpy(dest, src, size);
      }
      else                                     /* brin avec gaps: alignement */
      {
          mask->pattern->std[*m].len = size;     /* init. nec. a 'AlignBack' */

          AlignSProfile(src, mask->pattern->std + *m);
          AlignBack(src, mask->pattern->std + *m);

          src  = mask->pattern->std[*m].str;
          size = mask->pattern->std[*m].max_len;
 
          strncpy(dest, src, size);
      }
  }

  mask->str[mask->max_len] = '\0';

  if (mask->mode != NOMASK)               /* affichage des brins non traites */
  {
      GetCfgAtoms(cfg, mask);

      for (i = 0; i < mask->natom; i++)
      {
          if (mask->atomstr[i] == '0') /* restriction aux atomes hors masque */
          {
              int bgn = mask->atom[i-1].bgn + mask->atom[i-1].len;

              src = seq + bgn;
              dest = mask->str + mask->atom[i].db_bgn - mask->db_bgn;

              while (mask->atomstr[i] == '0') i++;
              size = mask->atom[i].bgn - bgn;
              i--;

              strncpy(dest, src, size);
              ToLower(dest, size);
          }
      }
  }
                                 /* place des "separateurs" entre les atomes */
  for (i = 0; i < mask->natom - 1; i++)
  {
      if (mask->atomstr[i] == '1')
          CharInsert(mask->str + mask->atom[i].db_bgn -
                     mask->db_bgn + i,
                     mask->atom[i].max_len, '.');
      else
      {
          int j, separs = 0;
          while (mask->atomstr[i] == '0') {  /* compte les brins hors masque */
              i++;
              separs++;
          }
          i--;            /* les "separateurs" sont regroupes en fin de brin */
          for (j = 0; j < separs; j++)
              CharInsert(mask->str + mask->atom[i].db_bgn -
                         mask->db_bgn + i - separs + 1,
                         mask->atom[i].max_len, '.');
      }
  }

  return;
}
/*=============================================================================
 CharInsert(): Inserre le caractere 'insert' a la position 'pos' de la chaine
               pointee par 'str', qui doit se terminer par '\0' et disposer
               d'un volume suffisant.
               'pos' doit etre inferieur ou egal a la longueur de 'str', aucun
               controle n'est effectue.
=============================================================================*/

void CharInsert(char *str, int pos, char insert)
{
  int  i, len = strlen(str);
  char *s;

  for (i = len, s = str + len; i >= pos; i--, s--)  *(s + 1) = *s;

  str[pos] = insert;

  return;
}
/*=============================================================================
 ToLower(): Convertit aux minuscules la chaine pointee par 'seq' sur une lon-
            gueur egale a 'len'.
=============================================================================*/

void ToLower(char *str, int len)
{
  int i;

  for (i = 0; i < len; i++) str[i] = tolower(str[i]);

  return;
}
/*=============================================================================
 GetCfgAtoms(): A partir d'une configuration de masque pointee par 'cfg' cette
                fonction intialise les champs 'bgn' et 'len' des atomes du
                masque pointe par 'mask' (le code est inspire de 'GetCfg' dans
                'maskcfg.c').
=============================================================================*/

void GetCfgAtoms(Config *cfg, Mask *mask)
{
  int  i, offset, h2indx,
       stindx, hxindx;        /* indices dans les tab. de brins et d'helices */ 

  offset = mask->atom - mask->pattern->atom;

  for (i = stindx = hxindx = 0; i < mask->natom; i++)
  {
      if (mask->atomstr[i] == '1')       /* restriction aux atomes du masque */
      {
          switch (mask->atom[i].type)
          {
              case STD:                                              /* brin */   
                  mask->atom[i].bgn = (int) cfg->stbgn[stindx];
                  mask->atom[i].len = (int) cfg->stlen[stindx];
                  stindx++;
                  break;
              case HLX1:                                /* 1er brin d'helice */
                  mask->atom[i].bgn = (int) cfg->hxbgn[hxindx];
                  mask->atom[i].len = (int) mask->atom[i].min_len;

                  h2indx = mask->atom[i].h2index - offset;      /* 2eme brin */
                  mask->atom[h2indx].bgn = mask->atom[i].bgn + mask->atom[i].len
                                           + cfg->hxdist[hxindx];
                  mask->atom[h2indx].len = mask->atom[i].len;
                  hxindx++;
                  break;
              default:
                  break;
          }
      }
  }
  return;
}
/*===========================================================================*/
