
/*=============================================================================
 pattern.c                    A.Lambert le 23/09/01   revu le 05/10/01

 Fonctions concernant la lecture et l'analyse d'un 'pattern', structure consti-
 tuee d'une succession de brins, helices ou parties d'helices.

 cc -O2 -Wall -c pattern.c -I../include ;

 ar -rs ../lib/librnaIV.a pattern.o ;

=============================================================================*/

#include "rnaIV.h"

Pattern *NewPattern(void);
void    FreePattern(Pattern *pattern);
void    DelPattern(Pattern *pattern);

void    ReadPatternId(char *pattern_id, int *bgn, int *end);

Pattern *ReadPatternAtoms(char *pattern_id, Trset *trset);
Pattern *ReadPattern(char *pattern_id, Trset *trset);

void    sPrintfPattern(char *str, Pattern *pattern);
void    sPrintfPatternStruct(char *str, Pattern *pattern);

/*=============================================================================
 NewPattern(): Cree une nouvelle structure 'pattern', initialise les pointeurs
               a NULL;
=============================================================================*/

Pattern *NewPattern(void)
{
  Pattern *pattern;

  pattern = (Pattern *) malloc(sizeof(Pattern));

  pattern->id = NULL;
  pattern->model = NULL;
  pattern->atlist = pattern->hxlist = pattern->stlist = NULL;
  pattern->oatlist = pattern->ohxlist = pattern->ostlist = NULL;

  pattern->atom = NULL;
  pattern->hlx  = NULL;
  pattern->std  = NULL;

  pattern->gapslist_ori = NULL;
  pattern->gapslist = NULL;

  return pattern;
}
/*=============================================================================
 FreePattern(): libere la memoire d'une structure 'pattern'.        
=============================================================================*/

void FreePattern(Pattern *pattern)
{

  if (pattern->id != NULL)      { free(pattern->id);      pattern->id    = NULL;}
  if (pattern->model != NULL)   { free(pattern->model);   pattern->model = NULL;}
 
  if (pattern->atlist != NULL)  { free(pattern->atlist);  pattern->atlist = NULL;}
  if (pattern->hxlist != NULL)  { free(pattern->hxlist);  pattern->hxlist = NULL;}
  if (pattern->stlist != NULL)  { free(pattern->stlist);  pattern->stlist = NULL;}

  if (pattern->oatlist != NULL) { free(pattern->oatlist); pattern->oatlist = NULL;}
  if (pattern->ohxlist != NULL) { free(pattern->ohxlist); pattern->ohxlist = NULL;}
  if (pattern->ostlist != NULL) { free(pattern->ostlist); pattern->ostlist = NULL;}

  if (pattern->atom != NULL) { free(pattern->atom);  pattern->atom = NULL;}
  if (pattern->hlx != NULL)   FreeHlxTable(pattern->hlx, pattern->nhx);
  if (pattern->std != NULL)   FreeStdTable(pattern->std, pattern->nst);

  return;
}
/*=============================================================================
 DelPattern(): detruit la structure pointee par 'pattern'.        
=============================================================================*/

void DelPattern(Pattern *pattern)
{
  FreePattern(pattern);
  free(pattern);
  return;
}
/*=============================================================================
 ReadPatternId(): Controle l'argument 'pattern_id' qui doit etre une chaine,
                  lue parmi les arguments de 'main', formee de 2 entiers sepa-
                  res par une virgule seulement.
                  En retour 'bgn' et 'end' pointent les entiers codant les
                  les atomes extremes d'un pattern.
                  Exemples: "-12,17" ou "+3,-5" ou "3,-5" etc..
=============================================================================*/

void ReadPatternId(char *pattern_id, int *bgn, int *end)
{
  char *ptc1, *ptc2;

  if ((ptc1 = strchr(pattern_id, ',')) == NULL || *(++ptc1) == '\0') {
      fprintf(stderr, "ReadPatternId: invalid argument '%s', exit..\n",
              pattern_id);
      exit(1);
  }

  *bgn = (int) strtol(pattern_id, &ptc1, 10);

  if (ptc1 == pattern_id || *ptc1 != ',') {
      fprintf(stderr, "ReadPatternId: invalid argument '%s', exit..\n",
              pattern_id);
      exit(1);
  }
  ptc2 = ++ptc1;
  *end = (int) strtol(ptc2, &ptc1, 10);

  if (ptc1 == ptc2) {
      fprintf(stderr, "ReadPatternId: invalid argument '%s', exit..\n",
              pattern_id);
      exit(1);
  }
  return;
}
/*=============================================================================
 ReadPatternAtoms(): Enregistre les elements structurels d'un 'pattern' depuis
                     une structure 'trset'.
                     Cette fonction ne cree pas les tableaux d'helices, strands
                     et configurations.
=============================================================================*/

Pattern *ReadPatternAtoms(char *pattern_id, Trset *trset)
{
  Pattern *pattern;
  int     i, j, h1, h2, *pti, *pti1, *pti2, bgn, end, bgn_index;

  ReadPatternId(pattern_id, &bgn, &end);     /* lecture et controle de l'id. */

                                                 /* detection dans 'oatlist' */
  if (bgn < 0)        
     pti1 = TabSearch(trset->oatlist, - bgn);                    /* 1er brin */
  else
  if ((pti = TabSearchn(trset->oatlist, bgn, 2)) != NULL)
     pti1 = pti;                                                /* 2eme brin */
  else
     pti1 = TabSearch(trset->oatlist, bgn);                   /* simple brin */

  if (end < 0)        
     pti2 = TabSearch(trset->oatlist, - end);                    /* 1er brin */
  else
  if ((pti = TabSearchn(trset->oatlist, end, 2)) != NULL)
     pti2 = pti;                                                /* 2eme brin */
  else
     pti2 = TabSearch(trset->oatlist, end);                   /* simple brin */

  if (pti1 == NULL || pti2 == NULL) {
      fprintf(stderr, 
             "ReadPattern: invalid argument '%s', exit..\n", pattern_id);
      exit(1);
  }

  pattern = NewPattern();                        /* creation de la structure */

  pattern->id = (char *) calloc(strlen(pattern_id) + 1, 1);
  strcpy(pattern->id, pattern_id);

  bgn_index = pti1 - trset->oatlist;           /* recodage de 'bgn' et 'end' */
  pattern->bgnid = *(trset->atlist + bgn_index);
  pattern->endid = *(trset->atlist + (pti2 - trset->oatlist));/* et de 'end' */

  pattern->trsetbgnindex = bgn_index;  /* indice dans la liste des atomes de */
                                        /* 'trset' du 1er atome de 'pattern' */
  pattern->natom = pti2 - pti1 + 1;

  pattern->atlist = (int *) calloc(pattern->natom + 1, sizeof(int));
  pattern->oatlist = (int *) calloc(pattern->natom + 1, sizeof(int));

                           /* dresse la liste des identificateurs des atomes */
  
  for (i = 0; i < pattern->natom; i++) 
  {
      pattern->atlist[i] = trset->atlist[i + bgn_index];
      pattern->oatlist[i] = trset->oatlist[i + bgn_index];
  }
                                       /* et initialise les stuctures 'atom' */

  pattern->atom = ReadAtoms(pattern->atlist, trset);

  if (pattern->atom == NULL) {
      fprintf(stderr, "ReadPatternAtoms: error reading atom list, exit..\n");
      exit(1);
  }

  for (i = 0; i < pattern->natom; i++) pattern->atom[i].h2index = 0;

                                                     /* et les autres champs */

  pattern->db_bgn = pattern->atom[0].db_bgn;
  pattern->max_len = 0;
  pattern->max_gaps = 0;

  for (i = 0; i < pattern->natom; i++)
  {
      pattern->max_len += pattern->atom[i].max_len;
      pattern->max_gaps += pattern->atom[i].max_gaps;
  }
  pattern->min_len = pattern->max_len - pattern->max_gaps;

                            /* compte les brins de longueur variables (gaps) */

  pattern->nvarst = 0;
  for (i = 0; i < pattern->natom; i++)
      if (pattern->atom[i].max_gaps != 0) (pattern->nvarst)++;
 
                   /* compte les helices de 'trset' presentes dans 'pattern' */
  pattern->nhx = 0;
  for (i = 0; i < trset->nhx; i++)
  {
      if (TabSearch(pattern->atlist, - trset->hxlist[i]) != NULL &&
          TabSearch(pattern->atlist,   trset->hxlist[i]) != NULL)
          (pattern->nhx)++;
  }
  pattern->hxlist  = (int *) calloc(pattern->nhx + 1, sizeof(int));
  pattern->ohxlist = (int *) calloc(pattern->nhx + 1, sizeof(int));
                     
                                                        /* et les enregistre */

  for (i = 0, j = 0; i < trset->nhx; i++)
  {                                                          /* si helice .. */
      if ((pti1 = TabSearch(pattern->atlist, - trset->hxlist[i])) != NULL &&
          (pti2 = TabSearch(pattern->atlist,   trset->hxlist[i])) != NULL)
      {
          pattern->hxlist[j]  = trset->hxlist[i];
          pattern->ohxlist[j] = trset->ohxlist[i];

          h1 = pti1 - pattern->atlist;
          h2 = pti2 - pattern->atlist;
          pattern->atom[h1].type = HLX1;
          pattern->atom[h2].type = HLX2;
          pattern->atom[h1].index = pattern->atom[h2].index = j;
          j++;
      }
  }
                                   /* compte les brins isoles dans 'pattern' */
                       /* qui peuvent etre des parties d'helices incompletes */

  pattern->nst = pattern->natom - 2 * pattern->nhx;

  pattern->stlist = (int *) calloc(pattern->nst + 1, sizeof(int));
  pattern->ostlist = (int *) calloc(pattern->nst + 1, sizeof(int));

                                                        /* et les enregistre */

  for (i = 0, j = 0; i < pattern->natom; i++)
  {
      pti1 = TabSearch(pattern->atlist, - pattern->atlist[i]);
      pti2 = TabSearch(pattern->atlist,   pattern->atlist[i]);
      if (pti1 == NULL || pti2 == NULL)                        /* si brin .. */
      {
          pattern->stlist[j]  = pattern->atlist[i];
          pattern->ostlist[j] = pattern->oatlist[i];
          h1 = pti1 - pattern->atlist;
          h2 = pti2 - pattern->atlist;

          h1 = (pti1 == NULL ? h2 : h1);
          pattern->atom[h1].type = STD;
          pattern->atom[h1].index = j;
          j++;
      }
  }
                         /* compte le nbre total de configurations possibles */

  GetPatternCfgsNb(pattern);
                                /* extrait le modele de structure secondaire */

  pattern->model = (int *) calloc(pattern->max_len + 1, sizeof(int));

  pti1 = trset->model + pattern->db_bgn;
  for (i = 0; i < pattern->max_len; i++)
      pattern->model[i] = *(pti1++);

  return pattern;
}
/*=============================================================================
 ReadPattern(): Enregistre les elements structurels d'un 'pattern' incluant les
                atomes, helices et brins (mais pas les configs.) depuis une
                structure 'trset'.
=============================================================================*/

Pattern *ReadPattern(char *pattern_id, Trset *trset)
{
  Pattern *pattern;

  pattern = ReadPatternAtoms(pattern_id, trset);
  ReadHelices(pattern);
  ReadStrands(pattern);

  return pattern;
}
/*=============================================================================
 sPrintfPattern(): Charge dans la chaine 'str' la sequence des atomes dans 
                   l'ordre ou ils apparaissent dans la structure pointee par
                   'pattern'.
                   Les indices de brins et d'helices dans leur tableau respec-
                   tif sont indiques pour chaque atome.
                   Destinee a visualiser la composition en helices et brins de
                   la structure 'pattern'.
=============================================================================*/

void sPrintfPattern(char *str, Pattern *pattern)
{
  int  i, t;

  for (i = 0, str[0] = '\0'; i < pattern->natom; i++)
  {
      t = pattern->atom[i].type;
      sprintf(str, "%s%s%s[%d]",
              str, i != 0 ? "." : "",                          /* separateur */
              t == STD ? "St" : t == HLX1 ? "H1" : "H2",             /* type */
              pattern->atom[i].index);                             /* indice */
  }
  return;
}
/*=============================================================================
 sPrintfPatternStruct(): Charge dans la chaine 'str' la sequence des symboles
                         'S' ou 'H' suivant que le caractere correspondant
                         appartient a un brin ou une helice dans le 'pattern'
                         pointe en argument.
                         Destinee a visualiser la composition en helices et
                         brins de la structure 'pattern'.
=============================================================================*/

void sPrintfPatternStruct(char *str, Pattern *pattern)
{
  int   i, j, t, l;
  char  symbol;

  for (i = 0, str[0] = '\0'; i < pattern->natom; i++)
  {
      t = pattern->atom[i].type;
      l = pattern->atom[i].max_len;
      symbol = (t == STD ? 'S' : 'H');

      sprintf(str, "%s%s", str, i != 0 ? " " : "");            /* separateur */

      for (j = 0; j < l; j++)
          sprintf(str, "%s%c", str, symbol);   
  }

  return;
}
/*===========================================================================*/
