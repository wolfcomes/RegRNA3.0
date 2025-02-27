
/*=============================================================================
 strand.c                          A.Lambert le 15/10/01

 Ce code concerne les fonctions relatives aux manipulations des structures
 'Strand'.

 L'identificateur d'un brin est le nombre entier qui lui est affecte dans le
 codage original de la structure, enregitre dans le champ 'ostlist' de 'pattern'
 (ceci facilitera les operations d'entree/sortie, comme pour 'pattern').

 cc -O2 -Wall -c strand.c -I../include ;

 ar -r ../lib/librnaIV.a strand.o ;

=============================================================================*/

#include "rnaIV.h"

Strand *SetStdTable(int nst);
void   FreeStdTable(Strand *Std, int nst);
void   DelStdTable(Strand *Std, int nst);
void   ReadStrand(Strand *Std, Pattern *pattern, int std_id);
void   ReadStrands(Pattern *pattern);

/*=============================================================================
 SetStdTable(): Cree un tableau de 'nst' brins en initialisant a 'NULL' les 
                pointeurs, si 'nst' est nul retourne 'NULL'.
=============================================================================*/

Strand *SetStdTable(int nst)
{
  Strand *Std;
  int    i;

  if (nst == 0) return NULL;

  Std = (Strand *) malloc(nst * sizeof(Strand));

  for (i = 0; i < nst; i++)
  {
      Std[i].str     = NULL;
      Std[i].Profile = NULL;
      Std[i].Align   = NULL;
      Std[i].Scores  = NULL;
      Std[i].ScoresBis  = NULL;
  }
  return Std;
}
/*=============================================================================
 FreeStdTable(): Libere la memoire allouee au tableau de brins pointe par 'Std'.
=============================================================================*/

void FreeStdTable(Strand *Std, int nst)
{
  int   i;

  if (nst == 0) return;

  for (i = 0; i < nst; i++)
  {
      if (Std[i].str != NULL)     {
          free(Std[i].str);     Std[i].str = NULL;
      }
      if (Std[i].Profile != NULL) {
          free(Std[i].Profile); Std[i].Profile = NULL;
      }
      if (Std[i].Align != NULL)   {
          free(Std[i].Align);   Std[i].Align = NULL;
      }
      if (Std[i].Scores != NULL)  {
          free(Std[i].Scores);  Std[i].Scores = NULL;
      }
      if (Std[i].ScoresBis != NULL)  {
          free(Std[i].ScoresBis);  Std[i].ScoresBis = NULL;
      }
  }
  return;
}
/*=============================================================================
 DelStdTable(): Detruit le tableau de brins pointe par 'Std'.
=============================================================================*/

void DelStdTable(Strand *Std, int nst)
{
  FreeStdTable(Std, nst);
  free(Std);
  return;
}
/*=============================================================================
 ReadStrand(): Initialise les champs de la structure pointee par 'Std' a par-
               tir des donnees pointees par 'pattern' et 'trset'.
               'std_id' est l'identificateur, dans la liste 'ostlist' des 
               brins de 'pattern', du brin a enregistrer.
               Si le brin contient des gaps le tableau 'St->Align' sera cree et
               mis a 0 par la fonction 'GetStProfile' dans 'profs.c'.
=============================================================================*/

void ReadStrand(Strand *Std, Pattern *pattern, int std_id)
{
  int i, std_found, std_index, atom_index;

  std_index = atom_index = 0;

  for (i = std_found = 0; i < pattern->nst; i++)
     if (std_id == pattern->ostlist[i]) {
         std_found = 1;
         std_index = i;
         break;
     }
  if (!(std_found)) {
      fprintf(stderr, "ReadStrand: strand '%d' not found, exit..", std_id);
      exit(1);
  }
                                         /* cherche le brin parmi les atomes */

  for (i = std_found = 0; i < pattern->natom; i++)
      if (pattern->atom[i].type == STD && 
          pattern->atom[i].index == std_index)
      {
          std_found = 1;
          atom_index = i;
          break;
      }
  if (!(std_found)) {
      fprintf(stderr, "ReadStrand: strand '%d' not found, exit..", std_id);
      exit(1);
  }
  Std->id = std_id;
  Std->db_bgn = pattern->atom[atom_index].db_bgn;
  Std->min_bgn = pattern->atom[atom_index].min_bgn;
  Std->min_len = pattern->atom[atom_index].min_len;
  Std->max_len = pattern->atom[atom_index].max_len;
  Std->max_gaps = pattern->atom[atom_index].max_gaps;

  Std->len = Std->max_len;             /* pourra etre modifie ulterieurement */

  return;
}
/*=============================================================================
 ReadStrands(): Enregistre les donnees du tableau ds brins contenus dans la
                structure pointee par 'pattern'.
                Ce tableau est enregistre dans le champ 'std' de 'pattern'.
=============================================================================*/

void ReadStrands(Pattern *pattern)
{
  int i;

  if (pattern->nst == 0) return;

  pattern->std = SetStdTable(pattern->nst);

  for (i = 0; i < pattern->nst; i++)
  {
      ReadStrand(pattern->std + i, pattern, pattern->ostlist[i]);
  }
  return;
}
/*===========================================================================*/
