
/*=============================================================================
 helix.c                       A.Lambert le 13/10/01  revu le 13/11/01

 Ce code concerne les fonctions relatives aux manipulations des structures
 'Helix'.

 L'identificateur d'une helice est le nombre entier qui lui est affecte dans le
 codage original de la structure, enregitre dans le champ 'ohxlist' de 'pattern'
 (ceci facilitera les operations d'entree/sortie, comme pour 'pattern').

 cc -O2 -Wall -c helix.c -I../include ;

 ar -r ../lib/librnaIV.a helix.o ;

=============================================================================*/

#include "rnaIV.h"

Helix *SetHlxTable(int nhx);
void  FreeHlxTable(Helix *Hlx, int nhx);
void  DelHlxTable(Helix *Hlx, int nhx);
void  ReadHelix(Helix *Hlx, Pattern *pattern, int hlx_id);
void  ReadHelices(Pattern *pattern);

/*=============================================================================
 SetHlxTable(): Cree un tableau de 'nhx' helices en initialisant a 'NULL' les 
                pointeurs, si 'nhx' est nul retourne 'NULL'.
=============================================================================*/

Helix *SetHlxTable(int nhx)
{
  Helix *Hlx;
  int     i;

  if (nhx == 0) return NULL;

  Hlx = (Helix *) malloc(nhx * sizeof(Helix));

  for (i = 0; i < nhx; i++)
  {
      Hlx[i].Profile = NULL;
      Hlx[i].Scores  = NULL;
      Hlx[i].ScoresBis = NULL;
  }
  return Hlx;
}
/*=============================================================================
 FreeHlxTable(): Libere la memoire allouee au tableau d'helices pointe par 'Hlx'.
=============================================================================*/

void FreeHlxTable(Helix *Hlx, int nhx)
{
  int     i;

  if (nhx == 0) return;

  for (i = 0; i < nhx; i++)
  {
      if (Hlx[i].Profile != NULL){ 
          free(Hlx[i].Profile); Hlx[i].Profile = NULL;
      }
      if (Hlx[i].Scores != NULL) {
          free(Hlx[i].Scores);  Hlx[i].Scores = NULL;
      }
      if (Hlx[i].ScoresBis != NULL) {
          free(Hlx[i].ScoresBis);  Hlx[i].ScoresBis = NULL;
      }
  }
  return;
}
/*=============================================================================
 DelHlxTable(): Detruit le tableau d'helices pointe par 'Hlx'.
=============================================================================*/

void DelHlxTable(Helix *Hlx, int nhx)
{
  FreeHlxTable(Hlx, nhx);
  free(Hlx);
  return;
}
/*=============================================================================
 ReadHelix(): Initialise les champs de la structure pointee par 'Hlx' a partir
              des donnees pointees par 'pattern'.
              'hlx_id' est l'identificateur, dans la liste 'ohxlist' des 
              helices de 'pattern', de l'helice a enregistrer.
=============================================================================*/

void ReadHelix(Helix *Hlx, Pattern *pattern, int hlx_id)
{
  int   i, hlx_found, hlx_index, H1_index, H2_index;

  hlx_index = H1_index = H2_index = 0;

  for (i = hlx_found = 0; i < pattern->nhx; i++)
     if (hlx_id == pattern->ohxlist[i]) {
         hlx_found = 1;
         hlx_index = i;
         break;
     }
  if (!(hlx_found)) {
      fprintf(stderr, "ReadHelix: helix '%d' not found, exit..", hlx_id);
      exit(1);
  }
                                          /* cherche le 1er brin de l'helice */

  for (i = hlx_found = 0; i < pattern->natom; i++) 
      if (pattern->atom[i].type == HLX1 && 
          pattern->atom[i].index == hlx_index)
      {
          hlx_found = 1;
          H1_index = i;
          break;
      }
  if (!(hlx_found)) {
      fprintf(stderr, "ReadHelix: strand 1 of helix '%d' not found, exit..",
                       hlx_id);
      exit(1);
  }
                                     /* cherche le deuxieme brin de l'helice */

  for (i = H1_index + 1, hlx_found = 0; i < pattern->natom; i++) 
  {
      if (pattern->atom[i].type == HLX2 && 
          pattern->atom[i].index == hlx_index)
      {
          hlx_found = 1;
          H2_index = i;
          pattern->atom[H1_index].h2index = H2_index;
          break;
      }
  }
  if (!(hlx_found))  {
      fprintf(stderr, "ReadHelix: strand 2 of helix '%d' not found, exit..",
                       hlx_id);
      exit(1);
  }
  Hlx->id = hlx_id;
  Hlx->helix_len = pattern->atom[H1_index].max_len;

  if (Hlx->helix_len > HLX_MAX_LEN) {
      fprintf(stderr, "ReadHelix: helix too long (max: %d), exit..\n", HLX_MAX_LEN);
      exit(1);
  }

  Hlx->db_bgn1 = pattern->atom[H1_index].db_bgn;
  Hlx->db_bgn2 = pattern->atom[H2_index].db_bgn;
  Hlx->min_bgn = pattern->atom[H1_index].min_bgn;

  Hlx->max_dist = Hlx->db_bgn2 - Hlx->db_bgn1 - Hlx->helix_len;
  Hlx->max_len  = 2 * Hlx->helix_len + Hlx->max_dist;

                   /* decompte dans  les 'atoms' du nombre maximal de gaps */

  for (i = H1_index + 1, Hlx->max_gaps = 0; i < H2_index; i++)
  {
      Hlx->max_gaps += pattern->atom[i].max_gaps;
  }

  Hlx->min_dist = Hlx->max_dist - Hlx->max_gaps;
  Hlx->min_len  = Hlx->max_len - Hlx->max_gaps;

  if (Hlx->max_gaps == 0) Hlx->dist = Hlx->max_dist;

  return;
}
/*=============================================================================
 ReadHelices(): Enregistre les donnees du tableau d'helices contenues dans la
                structure pointee par 'pattern'.
                Ce tableau est enregistre dans le champ 'hlx' de 'pattern'.
=============================================================================*/

void ReadHelices(Pattern *pattern)
{
  int i;

  if (pattern->nhx == 0) return;

  pattern->hlx = SetHlxTable(pattern->nhx);

  for (i = 0; i < pattern->nhx; i++)
  {
      ReadHelix(pattern->hlx + i, pattern, pattern->ohxlist[i]);
  }
  return;
}
/*===========================================================================*/
