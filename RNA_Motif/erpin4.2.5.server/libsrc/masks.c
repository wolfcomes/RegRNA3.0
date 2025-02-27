
/*=============================================================================
 masks.c                        A.Lambert le 19/12/01  revu le 05/05/02

 On cree ici des fonctions permettant la creation et la gestion de tableaux de
 masques. il s'agit surtout d'interfaces ou de variantes de fonctions ecrites
 dans 'mask.c'.

 cc -O2 -Wall -c masks.c -I ../include ;

 ar -rs ../lib/librnaIV.a masks.o ;

=============================================================================*/

#include "rnaIV.h"


Mask *NewMasks(int nmask);
void FreeMasks(Mask *mask, int nmask);
void DelMasks(Mask *mask, int nmask);
int  Read1MaskArgs(int argc, char *argv[], int *index, Mask *mask);
void ParseMasksArgs(Mask *mask, int nmask, Pattern *pattern);

Mask *ReadMasksArgs(int argc, char *argv[], int *nmask);

void GetMasks(Mask *mask, int nmask, Pattern *pattern);
void GetMasksCfgs(Mask *mask, int nmask);
void SetMasksScoresTabs(Context *ctxt, int nctxt);

/*=============================================================================
 NewMasks(): Cree un tableau de 'nmask' structures 'Mask', initialise les
             champs pointeurs susceptibles de se voir allouer de la memoire a
             NULL.
=============================================================================*/

Mask *NewMasks(int nmask)
{
  Mask *mask;
  int  i;

  mask = (Mask *) malloc(nmask * sizeof(Mask));

  if (mask == NULL) {
      fprintf(stderr, "NewMasks: allocation failure, exit..\n");
      exit(1);
  }
  
  for (i = 0; i < nmask; i++)
  {
      mask[i].args = NULL;
      mask[i].hxindex = NULL;
      mask[i].stindex = NULL;
      mask[i].atomstr = NULL;
      mask[i].gapslist = NULL;
      mask[i].gapslist_ori = NULL;
      mask[i].gapscfg = NULL;
      mask[i].cfg = NULL;
      mask[i].str = NULL;
  }

  return mask;  
}
/*=============================================================================
 FreeMasks(): Libere la memoire allouee aux champs d'un tableau de 'nmask' 
              structures Mask et reinitialise les elements pointeurs a NULL;
=============================================================================*/

void FreeMasks(Mask *mask, int nmask)
{
  int i;

  FreePatternGapList(mask->pattern, nmask);

  for (i = 0; i < nmask; i++)  FreeMask(mask + i);

  return;
}
/*=============================================================================
 DelMasks(): detruit un tableau de 'nmask' structures Mask.
=============================================================================*/

void DelMasks(Mask *mask, int nmask)
{
  FreeMasks(mask, nmask);
  free(mask);
  return;
}
/*=============================================================================
 Read1MaskArgs(): Saisit dans les arguments de '*argv[]', ex: parmi les argu-
                  ments de 'main', les elements d'un masque de pattern, sous
                  forme d'une sequences de chaines representant des symboles
                  d'helices et brins.
                  En retour la structure pointee par 'mask' est partiellement
                  initialisee: le tableau 'mask->args' des arguments convertis
                  en entiers, si il y en a et la fonction retourne 1 sinon 0
                  si le mode n'est pas NOMASK (cad: detection d'une erreur).
                  Le nombre 'mask->nargs' d'arguments saisis, ainsi que le mode
                  'mask->mode' de saisie (MASK, UMASK, ADDMASK ou NOMASK) sont
                  enregistres.
                  L'indice dans 'argv' de la fin des arguments est pointe en
                  sortie pour une relance ulterieure de la fonction pour la
                  saisie des arguments d'un autre masque avec les arguments:
                  'argc - index' et 'argv + index'.
=============================================================================*/

int Read1MaskArgs(int argc, char *argv[], int *index, Mask *mask)
{
  int   i, j, k, nargs, mode = 0, bgnmask = 0, diff;


  for (i = nargs = 0; i < argc; i++)         /* compte les elements du masque */
  {
      if (strcmp(argv[i], "-mask") == 0   ||
          strcmp(argv[i], "-umask") == 0  ||
          strcmp(argv[i], "-add") == 0    ||
          strcmp(argv[i], "-nomask") == 0
         )
      {
          if (strcmp(argv[i], "-mask") == 0)
              mode = MASK;
          else
          if (strcmp(argv[i], "-umask") == 0)
              mode = UMASK;
          else
          if (strcmp(argv[i], "-add") == 0)
              mode = ADDMASK;
          else
          {
              mode = NOMASK;
              *index = ++i;
              break;                  /* il n'y a pas d'autre argument a lire */
          }

          *index = bgnmask = ++i;

          for (j = i; j < argc && argv[j][0] != '-'; j++)     /* jusqu'a une */
          {                                               /* nouvelle option */
              *index = j;
              diff = 1;
              for (k = bgnmask; k < j; k++)        /* exclue les repetitions */
                  if ((diff = strcmp(argv[k], argv[j])) == 0)
                      break;
              if (diff != 0) nargs++;
          }          
          break;                                /* enregistre un seul masque */
      }
  }

  if (nargs == 0 && mode != NOMASK) return 0;             /* erreur detectee */

  mask->mode = mode;

  if (mask->mode != NOMASK)  /* les arg. de NOMASK seront lus ulterieurement */
  {
      mask->args = (int *) malloc(nargs * sizeof(int));
      for (i = 0; i < nargs; i++)
          mask->args[i] = atoi(argv[bgnmask + i]);

      mask->nargs = nargs;
  }
  return 1;
}
/*=============================================================================
 ReadMasksArgs(): Interface de 'Read1MaskArgs' destine a operer sur les argu-
                  ments de 'main' pour y saisir les arguments d'un nombre de
                  masques pointe, en retour par 'nmask'.
                  Retourne un pointeur sur le tableau de masques cree, dont les
                  champs seront initialises ulterieurement.

                  Les options precedant les arguments des masques sont:
                  '-mask', '-umask', '-add' ou '-nomask' dans le cas ou le 
                  masque represente la totalite du pattern associe dans les
                  arguments.

                  Les elements du tableau de masques sont ranges suivant leur
                  ordre dans la lecture de 'argv'.
                  Si aucun argument '-mask', '-umask', '-add' ou '-nomask' un
                  message est emis precedant la sortie du programme.
=============================================================================*/

Mask *ReadMasksArgs(int argc, char *argv[], int *nmask)
{
  int   i, index;
  Mask  *mask;

  for (i = *nmask = 0; i < argc; i++)         /* compte le nombre de masques */
  {
      if (argv[i][0] == '-') {
          if (strcmp(argv[i], "-mask") == 0   ||
              strcmp(argv[i], "-umask") == 0  ||
              strcmp(argv[i], "-add") == 0    ||
              strcmp(argv[i], "-nomask") == 0
             )
              (*nmask)++;
      }      
  }

  if (*nmask == 0)                                  /* pas d'argument 'Mask' */ 
  {
      fprintf(stderr, "ReadMasksArgs: needs at least one argument\n");
      fprintf(stderr, "               -umask |-mask |-add |-nomask, exit..\n");
      exit(1);
  }  

  mask = NewMasks(*nmask);         /* creation et init. des pointeurs a NULL */

  for (i = index = 0; i < *nmask; i++)
  {
      argc -= index;
      argv += index; 
      if (Read1MaskArgs(argc, argv, &index, mask + i) == 0)
      {
          fprintf(stderr,
                  "ReadMasksArgs: Error reading 'mask' arguments, exit..\n");
          exit(1);
      }
  }

  return mask;
}
/*=============================================================================
 ParseMasksArgs(): Analyse les arguments d'un tableau de 'namsk' masques poin-
                   te par 'mask' provenant des options '-mask', '-umask' et
                   '-add' lues sur une ligne de commande et subordonnes a la
                   region pointee par 'pattern'.
                   Pour '-mask' et '-add' le tableau des arguments est revu et
                   enumere tous les elements consideres, en vue d'un traitement
                   unifie notamment par 'SetUMask'.
                   Le masque cree apres un argument '-mask' aura pour elements
                   ceux du complementaire dans 'pattern' des elements cites.
                   Le masque cree apres un argument '-add' cumulera ses propres
                   elements et ceux du masque qui le precede s'il y en a un,
                   les repetitions eventuelles sont eliminees.
                   Cette fonction utilise les informations reunies par la fonc-
                   tion 'ReadMasksArgs'.
=============================================================================*/

void ParseMasksArgs(Mask *mask, int nmask, Pattern *pattern)
{
  int i, j, k, nargs, diffs;

  for (i = 0; i < nmask; i++)
  {
      if (mask[i].mode == MASK)  /* select. du complementaire dans 'pattern' */
      {
          int *tmp1;

          nargs = pattern->nhx + pattern->nst - mask[i].nargs;
          tmp1 = (int *) malloc(nargs * sizeof(int));
          nargs = 0;

          for (j = 0; j < pattern->nhx; j++)
          {
               for (k = diffs = 0; k < mask[i].nargs; k++)
                   if (pattern->hlx[j].id != mask[i].args[k]) diffs++;

              if (diffs == mask[i].nargs) {
                  tmp1[nargs] = pattern->hlx[j].id;
                  nargs++;
              }
          }
          for (j = 0; j < pattern->nst; j++)
          {
               for (k = diffs = 0; k < mask[i].nargs; k++)
                   if (pattern->std[j].id != mask[i].args[k]) diffs++;

              if (diffs == mask[i].nargs) {
                  tmp1[nargs] = pattern->std[j].id;
                  nargs++;
              }
          }
          free(mask[i].args);
          mask[i].args = tmp1;
          mask[i].nargs = nargs;
      }

      else if (mask[i].mode == ADDMASK  &&  i != 0)

      {
          int *tmp2;

          nargs = mask[i-1].nargs + mask[i].nargs;
          tmp2 = (int *) malloc(nargs * sizeof(int));
          for (j = 0; j < mask[i-1].nargs; j++) tmp2[j] = mask[i-1].args[j];
          nargs = mask[i-1].nargs;

          for (k = 0; k < mask[i].nargs; k++)     /* elimine les repetitions */
          {
              for (j = diffs = 0; j < mask[i-1].nargs; j++) {
                  if (mask[i].args[k] != mask[i-1].args[j]) diffs++;
              }
              if (diffs == mask[i-1].nargs) {
                  tmp2[nargs] = mask[i].args[k];
                  nargs++;
              }
          }
          free(mask[i].args);
          mask[i].args = tmp2;
          mask[i].nargs = nargs;
      }
  }

  return;
}
/*=============================================================================
 GetMasks(): Interface assurant l'initialisation complete d'un tableau de mas-
             ques, a l'exception de l'ensemble de ses configurations, depuis
             l'examen du contenu d'un 'pattern'.
             Variante de 'GetMask' de 'mask.c'.
=============================================================================*/

void GetMasks(Mask *mask, int nmask, Pattern *pattern)
{
  int i;

  for (i = 0; i < nmask; i++)  GetMask(mask + i, pattern);

  GetPatternGapList(pattern, nmask);                 /* voir le code 'dmp.c' */

  return;
}
/*=============================================================================
 GetMasksCfgs(): Interface assurant le calcul des tableaux de configurations 
                 d'un tableau de 'nmask' pointe par 'mask'.
                 Les configurations de gaps ne sont pas retenues.
=============================================================================*/

void GetMasksCfgs(Mask *mask, int nmask)
{
  int i;

  for (i = 0; i < nmask; i++) {
      GetMaskCfgs(mask + i);
      DelGapCfgTable(mask + i);
  }

  return;
}
/*=============================================================================
 SetMasksScoresTabs(): Cree l'ensemble des tableaux de scores associes au ta-
                       bleau de 'nctxt' contextes pointe par 'ctxt'.
                       Chaque tableau de scores a une dimension a priori inde-
                       pendante.
                       Suivant la valeur de 'i' les tableaux 'Scores' (i = 0)
                       seront crees ou 'ScoresBis' (i != 0) dans les structures
                       'Helix' et 'Strand'.
=============================================================================*/

void SetMasksScoresTabs(Context *ctxt, int nctxt)
{
  int i;

  for (i = 0; i < nctxt; i++)
      SetMaskScoresTab(ctxt[i].mask, ctxt[i].scoretablen, i);

  return;
}
/*===========================================================================*/
