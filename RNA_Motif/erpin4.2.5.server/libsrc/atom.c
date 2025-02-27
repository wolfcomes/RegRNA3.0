/*=============================================================================
 atom.c                         A.Lambert le 17/09/01    revu le 13/11/01

 Fonctions concernant la gestion des 'atomes'.
 La structure 'Atom' decrit l'element elementaire dans un alignement multiple
 constitue d'un brin isole (code par un seul nombre entier) ou de l'un des 2
 brins d'une helice.
 Cette structure sera utilisee dans de nombreuses operations utilisant les pa-
 rametres geometriques des divers brins et leur codage, notamment dans l'etude
 des configurations d'ARN.

 cc -O2 -Wall -c atom.c -I../include ;

 ar -rs ../lib/librnaIV.a atom.o ;

=============================================================================*/

#include "rnaIV.h"

Atom  *ReadAtoms(int *atomlist, Trset *trset);

/*=============================================================================
 ReadAtoms(): Initialise depuis les donnees pointees par 'trset' le tableau
              d'atomes identifiees par le tableau pointe par 'atomlist' dont la
              fonction suppose qu'il est un sous-ensemble du modele de 'trset'
              d'elements contigus et ranges dans le meme odre.
              Un controle est effectue, si il est negatif la fonction retourne
              le pointeur NULL.
              Les champs "variables" 'bgn, len, gaps' sont initialises a des
              valeurs utiles pour les manipulations ulterieures.
              Rappel: les tableaux d'entiers doivent se terminer par 0.
=============================================================================*/

Atom *ReadAtoms(int *atomlist, Trset *trset)
{  
  Atom  *atoms;
  int   *pti, i, j, l, list_len, max, gaps;

  if (trset->model == NULL || trset->data == NULL) {
      fprintf(stderr, "ReadAtoms: trset data not allocated, exit..\n");
      exit(1);
  }
  if ((list_len = TabLen(atomlist)) == 0) {
      fprintf(stderr, "ReadAtoms: void argument 1, exit..\n");
      exit(1);
  }
                                    /* verifie que atomlist coincide avec un */
                                   /* sous-ensemble connexe de trset->atlist */

  if ((pti = TabSearch(trset->atlist, atomlist[0])) == NULL)  return NULL;

  for (i = 1; i < list_len; i++)
      if (atomlist[i] != pti[i])  return NULL;

  atoms = (Atom *) malloc(list_len * sizeof(Atom));

  for (i = 0; i < list_len; i++)                        /* lecture du modele */
  {
      pti = TabSearch(trset->model, atomlist[i]);

      atoms[i].db_bgn  = pti - trset->model;
      atoms[i].max_len = TabRepeats(pti, atomlist[i]);
      atoms[i].id = atomlist[i];
  }

  for (i = 0; i < list_len; i++)                  /* lecture de l'alignement */
  {
      max = 0;
      for (j = 0; j < trset->nseq; j++)
      {                                                   /* compte les gaps */
          for (l = 0, gaps = 0; l < atoms[i].max_len; l++)
          {
              if (trset->data[j][atoms[i].db_bgn + l] == '-') gaps++;
          }
          if (gaps > max) max = gaps;
      }
      atoms[i].max_gaps = max;
      atoms[i].min_len = atoms[i].max_len - atoms[i].max_gaps;
  }

  atoms[0].min_bgn = 0;
  for (i = 1; i < list_len; i++)
  {
      atoms[i].min_bgn = atoms[i-1].min_bgn + atoms[i-1].min_len;
  }

  for (i = 0; i < list_len; i++)      /* initialisation des 'bgn, len, gaps' */
  {
      atoms[i].len = atoms[i].min_len;
      atoms[i].gaps = 0;
  }
  atoms[0].bgn = 0;

  return atoms;
}
/*===========================================================================*/
