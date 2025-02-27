
/*=============================================================================
 trset.c                       A.Lambert le 17/09/01     revu le 20/06/02

 Fonctions concernant la lecture et l'analyse d'une base d'entrainement consti-
 tuee d'un alignement multiple de sequences d'ARN.

 cc -O2 -Wall -c trset.c -I../include ;

 ar -rs ../lib/librnaIV.a trset.o ;

=============================================================================*/

#include "rnaIV.h"

Trset *NewTrset(void);
void  FreeTrset(Trset *trset);
void  DelTrset(Trset *trset);
void  FreeTrsetSeqs(Trset *trset);

int   GetFileType(char *filename);
void  GetTrsetGeom(char *filename, Trset *trset, FILE *fout);
void  GetTrsetData(char *filename, Trset *trset);
void  GetTrsetStruct(Trset *trset);
void  CtrlTrsetStruct(Trset *trset, int mode, FILE *fout);
void  ProcessTrset(Trset *trset);

Trset *ReadTrset(char *filename, int mode, FILE *fout);

int   GetInputCode(Trset *trset, int col);
void  GetInputModel(int *model, Trset *trset);
char  *GapsStr(Trset *trset);
char  *StructStr(Trset *trset);
void  fPrintfAtoms(char *str, Trset *trset, int bgn, int len, FILE *txt);
void  fPrintfSubTrset(Trset *trset, int bgn, int len, FILE *txt);
void  fPrintfTrset(Trset *trset, FILE *txt);
void  fPrintfTrsetStat(char *filename, int nseq, int len, FILE *fout);
void  GetTrsetPatternId(Trset *trset, char *id);

int   StripTrset(char *filename, Trset *trset, char *outputname);

/*=============================================================================
 NewTrset(): Creation d'une structure de type 'Trset' dont les elements
             pointeurs sont initialises a NULL;
=============================================================================*/

Trset *NewTrset(void)
{
  Trset *trset = (Trset *) malloc(sizeof(Trset));

  trset->model   = NULL;
  trset->atlist  = trset->hxlist = trset->stlist = NULL;
  trset->oatlist = trset->ohxlist = trset->ostlist = NULL;
  trset->omodel  = NULL;
  trset->data    = NULL;

  trset->hlxsum = NULL;        /* pointeurs sur les matrices de substitution */
  trset->stsum = NULL;

  return trset;
}
/*=============================================================================
 FreeTrset(): libere la memoire d'une structure 'Trset' et celle allouee a ses
              champs.
=============================================================================*/

void  FreeTrset(Trset *trset)
{
  if (trset->model != NULL)   { free(trset->model);    trset->model   = NULL; }
  if (trset->atlist != NULL)  { free(trset->atlist);   trset->atlist  = NULL; }
  if (trset->oatlist != NULL) { free(trset->oatlist);  trset->oatlist = NULL; }
  if (trset->hxlist != NULL)  { free(trset->hxlist);   trset->hxlist  = NULL; }
  if (trset->ohxlist != NULL) { free(trset->ohxlist);  trset->ohxlist = NULL; }
  if (trset->stlist != NULL)  { free(trset->stlist);   trset->stlist  = NULL; }
  if (trset->ostlist != NULL) { free(trset->ostlist);  trset->ostlist = NULL; }
  if (trset->omodel != NULL)  FreecMat(trset->omodel); 
  if (trset->data != NULL)    FreecMat(trset->data);

  if (trset->hlxsum != NULL)  FreedMat(trset->hlxsum);
  if (trset->stsum != NULL)  FreedMat(trset->stsum);

  return;
}
/*=============================================================================
 DelTrset(): Detruit la structure pointee par'Trset'.
=============================================================================*/

void  DelTrset(Trset *trset)
{
  FreeTrset(trset);
  free(trset);
  return;
}
/*=============================================================================
 FreeTrsetSeqs(): Libere les champs 'data' et 'omodel' d'une structure 'Trset',
                  c'est a dire le tableau des sequences alignees et le codage 
                  original (entre au programme),
                  sera utilise pour librer la memoire de donnees devenues
                  inutiles.
=============================================================================*/

void FreeTrsetSeqs(Trset *trset)
{
  if (trset->data != NULL) {
      FreecMat(trset->data);   trset->data = NULL;
  }
  if (trset->omodel != NULL) {
      FreecMat(trset->omodel);  trset->omodel = NULL;
  }
  return;
}
/*=============================================================================
 GetFileType(): Indique si un fichier est un 'TRSET' ou 'DATASET' ou 'UNKNOWN'
                suivant qu'il commence ou non par une chaine de caracteres ex-
                clusivement numeriques (codage d'une struture secondaire d'ARN)
=============================================================================*/

int GetFileType(char *filename)
{
  FILE *txt;
  char *str;
  int  i, len, alphas, digits, buffer_sz = TRSET_MAX_LEN;

  if ((txt = fopen(filename, "r")) == NULL) {
      fprintf(stderr, "%s: file not found, exit.. \n", filename);
      exit(1);
  }
  str = (char *) calloc(buffer_sz + 1, sizeof(char));

                                   /* passe les 1eres lignes de commentaires */
  do {                                                /* et les lignes vides */
      if ( fgets(str, buffer_sz, txt) == NULL )
      {
          fprintf(stderr, "GetFileType: no data in '%s' file, exit..\n",
          filename);
          exit(1);
      }
      len = DataLen(str, strlen(str));
  }
  while (str[0] == '>' || len == 0);

  for (i = alphas = digits = 0; i < len; i++)
  {
      if (isalpha(str[i])) alphas++;
      else
      if (isdigit(str[i])) digits++;
  }

  free(str);
  fclose(txt);

  if (digits == len)           /* caracteres numeriques exclusivement: TRSET */
      return TRSET;
  if (alphas == len)      /* caracteres alphabetiques exclusivement: DATASET */
      return DATASET;            

      return UNKNOWN;                                      /* autre: UNKNOWN */
}
/*=============================================================================
 GetTrsetGeom(): Dans le fichier 'filename' d'une base d'entrainement cette 
                 fonction saisit les elements geometriques du modele de struc-
                 ture secondaire et des sequences de l'alignement multiple.
                 Un controle de la largeur de l'alignement est effectue.
                 - Les resultats sont retournes dans une structure 'Trset' 
                 pointee par 'trset',
                 Des indications sont sortis sur le fichier pointe par 'fout'.
=============================================================================*/

void GetTrsetGeom(char *filename, Trset *trset, FILE *fout)
{
  FILE   *txt;
  char   *str;
  int    model_len = 0, model_rows, nseq, len, nerrors, lines,
         buffer_sz = TRSET_MAX_LEN;

  if ((txt = fopen(filename, "r")) == NULL) {
      fprintf(stderr, "%s: file not found, exit.. \n", filename);
      exit(1);
  }
  str = (char *) malloc(buffer_sz);
  nerrors = 0;                             /* compteur des erreurs detectees */
  lines = 0;                               /* compteur des lignes du fichier */
  model_rows = 0;                           /* compteur des lignes du modele */
  nseq = 0;                                        /* compteur des sequences */

  while (fgets(str, buffer_sz, txt) != NULL)
  {
      lines++;
      len = DataLen(str, strlen(str));
      if (str[0] != '>' && len != 0)
      {
          if (isdigit(str[0]))                  /* ligne du modele de codage */
          {
              model_rows++;
              if (model_rows == 1) model_len = len;
              if (len != model_len) {
                  fprintf(stderr,
                  "GetTrsetGeom: Invalid model (unequal lengths: line %d)\n",
                   lines);
                  nerrors++;
              }
          }
          else                                   /* sequence de l'alignement */
          {
              nseq++;
              if (len != model_len) {
                  fprintf(stderr, 
                  "GetTrsetGeom: alignment error in seq. %d, line %d\n",
                  nseq, lines);
                  nerrors++;
              }
          }
      }

  }
  if (nerrors > 0)                   /* au moins une erreur detectee: sortie */
  {
      fprintf(stderr, "GetTrsetGeom: %d error%s found, exit.. \n",
              nerrors, (nerrors > 1 ? "s" : ""));
      free(str);
      fclose(txt); 
      exit(1);
  }
              /* on dispose desormais de 'model_rows', 'model_len' et 'nseq' */


  fPrintfTrsetStat(filename, nseq, model_len, fout);

  trset->nseq   = nseq;
  trset->len    = model_len;
  trset->digits = model_rows;

  free(str);
  fclose(txt);

  return;
}
/*=============================================================================
 GetTrsetData(): Dans le fichier 'filename' d'une base d'entrainement cette 
                 fonction enregistre les donnees du modele de structure secon-
                 daire et des sequences de l'alignement multiple, a partir des
                 indications geometriques contenues dans 'trset'.
                 - Le tout est retourne dans les champs d'une structure 'Trset'
                 pointee par 'trset', dont les tableaux sont crees.
                 Cette fonction utilise les donnees lues par 'GetTrsetGeom' qui
                 devra donc etre appelee avant.
=============================================================================*/

void GetTrsetData(char *filename, Trset *trset)
{
  FILE   *txt;
  char   *str;
  int    i, lines, buffer_sz = TRSET_MAX_LEN;

  if ((txt = fopen(filename, "r")) == NULL) {
      fprintf(stderr, "%s: file not found, exit.. \n", filename);
      exit(1);
  }
  str = (char *) malloc(buffer_sz);
  trset->omodel = cMat(trset->digits, trset->len + 1);
  trset->model  = (int *) calloc((trset->len + 1), sizeof(int));
  trset->data   = cMat(trset->nseq, trset->len + 1);

  lines = 0;                               /* compteur des lignes du fichier */

  for (i = 0; i < trset->digits;  )  /* copie des lignes de codage du modele */
  {
      lines++;
      fgets(str, buffer_sz, txt);
      if (isdigit(str[0]))
      {
          strncpy(trset->omodel[i], str, trset->len);
          trset->omodel[i][trset->len] = '\0';
          i++;                              /* compteur des lignes du modele */
      }      
  }
                     /* transfert des valeurs dans une seule ligne d'entiers */

  GetInputModel(trset->model, trset);

               /* controle qu'aucun element de 'trset->model' n'est egal a 0 */

  for (i = 0; i < trset->len; i++)
      if (trset->model[i] == 0) {
          fprintf(stderr, "GetTrsetData: Secondary structure code error\n");
          fprintf(stderr, "invalid column #%d\n", i);
          fprintf(stderr, "columns filled with 0 are not allowed, exit..\n");
          exit(1);
      }

                             /* copie des sequences de l'alignement multiple */
  i = 0;
  while (fgets(str, buffer_sz, txt) != NULL)
  {
      lines++;
      if (str[0] != '>' && DataLen(str, strlen(str)) != 0)  
      {
          strncpy(trset->data[i], str, trset->len);
          trset->data[i][trset->len] = '\0';
          i++;                                       /* compte les sequences */
      }
  }

  free(str);
  fclose(txt); 

  return;
}
/*=============================================================================
 GetTrsetStruct(): Etudie et recode la structure secondaire de 'trset'.

                   Le codage lu dans la base d'entrainement obeit aux regles
                   suivantes:
                   - Chaque colonne de l'alignement est code par un entier PO-
                   SITIF ou NUL, le meme sur toute la longueur d'un meme brin
                   ou d'une partie d'helice.
                   - Dans le cas d'une helice c'est le meme nombre qui code ses
                   deux brins: la repetition permet la detection des helices.

                   Le nouveau codage, issu de l'original, est destine a fournir
                   plus d'informations aux fonctions:
                   - Chaque colonne de l'alignement est code par un entier
                   STRICTEMENT POSITIF ce qui permet de reserver le nombre 0 
                   pour la detection de la fin d'un tableau.

                   CODAGE DES BRINS ISOLES:
                   - Les 'nst' brins isoles, indexes par 'i' dans l'ordre d'ap-
                   parition dans la liste originale des atomes 'oatlist' sont
                   codes par le nombre IMPAIR STRICTEMENT POSITIF: 2*i + 1.
                   L'index 'i' dans 'oatlist' du brin identifie par 'stdid' est
                   donc donne par: 'atlist[(stdid - 1)/2]'.
                   Trest *trset;
                   #define STDINDEX(trset, stdid)  (trset->atlist[(stdid-1)/2])

                   CODAGE DES HELICES:
                   - Les 'nhx' helices, indexees par 'i' dans l'ordre d'appa-
                   rition du 1er brin dans la liste originale 'oatlist' sont
                   codees par un nombre PAIR:
                   1er brin:  - (2 * (i + 1))       entier PAIR NEGATIF
                   2eme brin:   (2 * (i + 1))       entier PAIR POSITIF
                   helice:      (2 * (i + 1))   
                   L'index 'i' dans 'oatlist' du 1er brin de l'helice identifiee
                   par 'hlxid' est donc donne par: 'atlist[hlxid/2 - 1]'.
                   Trest *trset;
                   #define HLX1INDEX(trset, hlxid) (trset->atlist[(hlxid)/2-1])

                   Cette fonction utilise les donnees lues par 'GetTrsetData'
                   qui devra donc etre appelee avant.
=============================================================================*/

void GetTrsetStruct(Trset *trset)
{
  int i, j, count, *pti, errors;

  if (trset->model == NULL) {
      fprintf(stderr, "GetTrsetStruct: argument not allocated, exit..\n");
      exit(1);
  }
                                       /* compte le nombre total de symboles */

  for (i = count = 1; i < trset->len; i++)
  {
     if (trset->model[i] != trset->model[i - 1])
        count++;
  }
  trset->natom = count;
                                   /* cree et initialise la liste des atomes */

  trset->oatlist = (int *) calloc(trset->natom + 1, sizeof(int));

  trset->oatlist[0] = trset->model[0];
  for (i = 1, count = 0; i < trset->len; i++)
  {
     if (trset->model[i] != trset->model[i - 1])
         trset->oatlist[++count] = trset->model[i];
  }

                                         /* compte le nombre total d'helices */

  for (pti = trset->oatlist, count = 0; *pti != 0; pti++)
  {
      if (TabSearchn(pti, *pti, 2) != NULL)
          count++;
  }
  trset->nhx = count;
  trset->ohxlist = (int *) calloc(trset->nhx + 1, sizeof(int));
  trset->hxlist = (int *) calloc(trset->nhx + 1, sizeof(int));

          /* controle que les symboles d'helices ne sont detectes que 2 fois */

  for (pti = trset->oatlist, errors = 0; *pti != 0; pti++)
  {
      if (TabSearchn(pti, *pti, 3) != NULL) {           /* si 3 detections.. */
          fprintf(stderr, "GetTrsetStruct: Secondary structure code error\n");
          fprintf(stderr, "3 symbols '%d' (at least) detected\n", *pti);
          errors++;
      }
  }
  if (errors != 0) {
          fprintf(stderr, "exit..\n");
          exit(1);
  }
                                        /* initialise le tableau des helices */

  for (i = 0, pti = trset->oatlist; *pti != 0; pti++)
  {
      if (TabSearch(trset->ohxlist, *pti) == NULL && 
          TabSearchn(pti, *pti, 2) != NULL
         )
          trset->ohxlist[i++] = *pti;
  }

                          /* cree et initialise le tableau des simples brins */

  trset->nst = trset->natom - 2 * trset->nhx;             /* nombre de brins */

  trset->ostlist = (int *) calloc(trset->nst + 1, sizeof(int));
  trset->stlist = (int *) calloc(trset->nst + 1, sizeof(int));

  for (i = 0, pti = trset->oatlist; *pti != 0; pti++)
  {
      if (TabSearch(trset->ohxlist, *pti) == NULL)
          trset->ostlist[i++] = *pti;
  }

                         /* -- RECODAGE DU MODELE de structure secondaire -- */


  trset->atlist = (int *) calloc(trset->natom + 1, sizeof(int));

  for (i = 0; i < trset->nhx; i++)                                /* helices */
  {
      pti = TabSearch(trset->oatlist, trset->ohxlist[i]);        /* 1er brin */
      j = pti - trset->oatlist; 
      trset->atlist[j] = - 2 * (j + 1);

      pti = TabSearchn(trset->oatlist, trset->ohxlist[i], 2);   /* 2eme brin */
      trset->atlist[pti - trset->oatlist] = 2 * (j + 1);      /* de l'helice */

      trset->hxlist[i] = 2 * (j + 1);                              /* helice */
  }

  for (i = 0; i < trset->nst; i++)                                /* strands */
  {
      pti = TabSearch(trset->oatlist, trset->ostlist[i]);
      j = pti - trset->oatlist;       
      trset->atlist[j] = 2 * j + 1;

      trset->stlist[i] = 2 * j + 1;
  }

                      /* actualise 'trset->model' a l'aide du nouveau codage */

  for (i = count = 0; i < trset->len; )
  {
      if (trset->model[i] == trset->oatlist[count])
          trset->model[i++] = trset->atlist[count];
      else count++;
  }

  return;
}
/*=============================================================================
 CtrlTrsetStruct(): Controle la structure d'un alignement multiple pointe par
                    'trset': Verifie l'egalite des longueurs des 2 brins des
                    helices, et qu'il n'y ait pas de gaps dans les helices.
                    En cas d'echec, si mode == 'V', les donnees de 'trset' sont
                    dirigees sur le fichier 'fout'.
=============================================================================*/

void CtrlTrsetStruct(Trset *trset, int mode, FILE *fout)
{
  int i, j, l1, l2, nerrors = 0, symbol;

  for (i = 0; i < trset->nhx; i++)      /* controle des longueurs d'helices */
  {
      l1 = TabRepeats(trset->model, - (trset->hxlist[i]));      /* 1er brin */
      l2 = TabRepeats(trset->model, trset->hxlist[i]);         /* 2eme brin */
      if (l1 != l2) {
          nerrors++;
          fprintf(stderr,
          "CtrlTrsetStruct: non equal lengths (%d & %d) in helix '%d' code\n",
          l1, l2, trset->ohxlist[i]);
      }
  }

  for (j = 0; j < trset->len; j++)               /* colonnes de l'alignement */
  {
      symbol = trset->model[j];           /* symbol % 2 == 0 en cas d'helice */

      for (i = 0; i < trset->nseq; i++)          /* parcourt de la colonne j */
      {
          if (symbol % 2 == 0  &&  trset->data[i][j] == '-')
          {
              nerrors++;
              fprintf(stderr,
              "CtrlTrsetStruct: a gap in helix '%d', col. %d, seq. %d\n",
              GetInputCode(trset, j), j + 1, i + 1);
          }
      }
  }
  if (nerrors > 0)
  {
      if (toupper(mode) == 'V')  fPrintfTrset(trset, fout);

      fprintf(stderr, "CtrlTrsetStruct: %d alignment error%s found, exit..\n",
      nerrors, nerrors > 1 ? "s" : "");
      exit(1);
  }
  return;
}
/*=============================================================================
 ProcessTrset(): Traite les donnees saisies dans 'trset->data':
                 Passage aux majuscules.
                 Convertion de l'ARN en ADN.
                 Convertion des lettres differentes de A,T,G,C,- en N.
                 Preservation des gaps.
=============================================================================*/

void ProcessTrset(Trset *trset)
{
  int   i, j;
  char  c;

  for (i = 0; i < trset->nseq; i++)
      for (j = 0; j < trset->len; j++)
      {
          if (trset->data[i][j] != '-')
          {
              c = toupper(trset->data[i][j]);
              if (c == 'U')  c = 'T';
              if (c != 'A' && c != 'G' && c != 'C' && c != 'T')  c = 'N';

              trset->data[i][j] = c;
          }
      }
  return;
}
/*=============================================================================
 ReadTrset(): Interface des fonctions precedentes concernant la lecture des
              donnees d'une base d'entrainement et leur controle.
              En cas d'echec, si mode == 'V', les donnees de 'trset' sont diri-
              gees sur le fichier 'fout'.
=============================================================================*/

Trset *ReadTrset(char *filename, int mode, FILE *fout)
{
  Trset *trset;

  if (GetFileType(filename) != TRSET) {
      fprintf(stderr, "ReadTrset: inappropriate file '%s', exit..\n", filename);
      exit(1);
  }
  trset = NewTrset();
  GetTrsetGeom(filename, trset, fout);
  GetTrsetData(filename, trset);
  GetTrsetStruct(trset);
  CtrlTrsetStruct(trset, mode, fout);
  ProcessTrset(trset);

  return trset;
}
/*=============================================================================
 GetInputCode(): Retourne l'entier codant une colonne de 'trset', a l'origine
                 des operations (dans le fichier lui-meme).
                 S'applique a un 'trset' dont au moins la geometrie et le code 
                 ont ete enregistres.
=============================================================================*/

int GetInputCode(Trset *trset, int col)
{
  int i, code;
  char *str;

  str = (char *) malloc(trset->digits + 1);

  for (i = 0; i < trset->digits; i++)  str[i] = trset->omodel[i][col];
  str[trset->digits] = '\0';

  code = atoi(str);
  free(str);

  return code;
}
/*=============================================================================
 GetInputModel(): Charge dans le tableau pointe par 'model', le tableau 
                  d'entiers codant 'trset', a l'origine des operations (dans le
                  fichier lui-meme).
                  S'applique a un 'trset' dont au moins la geometrie et le code 
                  ont ete enregistres.
                  'model' doit pointer un volume suffisant.
=============================================================================*/

void GetInputModel(int *model, Trset *trset)
{
  int  i, j;
  char *str;

  str = (char *) malloc(trset->digits + 1);

  for (j = 0; j < trset->len; j++)
  {
      for (i = 0; i < trset->digits; i++)  str[i] = trset->omodel[i][j];
      str[trset->digits] = '\0';
      model[j] = atoi(str);
  }
  model[trset->len] = '\0';
  free(str);

  return;
}
/*=============================================================================
 GapsStr(): Retourne une chaine ou sont montres les gaps de 'trset'.
            Les colonnes ou des gaps sont presents sont indiques par '-'.
            Les colonnes constituees seulement de gaps sont indiques par '#',
            leur nombre est enregistre dans le champ 'trset->voidcols' en vue
            d'un traitement ulterieur du fichier des donnees.
            Les colonnes sans gaps sont indiques par '*'.
=============================================================================*/

char *GapsStr(Trset *trset)
{
  int  i, j, gaps;
  char *str = (char *) calloc(trset->len + 1, 1);
  
  trset->voidcols = 0;
  for (j = 0; j < trset->len; j++)               /* colonnes de l'alignement */
  { 
      str[j] = '*';
      for (i = gaps = 0; i < trset->nseq; i++)   /* parcourt de la colonne j */
          if (trset->data[i][j] == '-') gaps++;           /* compte les gaps */

      if (gaps != 0) str[j] = '-';
      if (gaps == trset->nseq) {             /* que des gaps dans la colonne */
          str[j] = '#';
          trset->voidcols++;
      }
  }

  return str;
}
/*=============================================================================
 StructStr(): Retourne une chaine ou, par 'H' ou 'S' sont montres les parties
              d'helices et les brins de 'trset'.
=============================================================================*/

char *StructStr(Trset *trset)
{
  int  j;
  char *str = (char *) calloc(trset->len + 1, 1);

  for (j = 0; j < trset->len; j++)               /* colonnes de l'alignement */
  { 
      if (trset->model[j] % 2 == 0) str[j] = 'H';
      else
      str[j] = 'S';
  }

  return str;
}
/*=============================================================================
 fPrintfAtoms(): Imprime les donnees d'une partie, limitee par 'bgn' et 'len',
                 d'une sequence 'str' de 'trset' en regroupant les 'atoms' de
                 l'alignement multiple.
=============================================================================*/

void fPrintfAtoms(char *str, Trset *trset, int bgn, int len, FILE *txt)
{
  int j, k;

  fprintf(txt, "%c", str[bgn]);
 
  for (j = 1, k = bgn + 1; j < len; j++, k++)    /* colonnes de l'alignement */
  {
      if (trset->model[k] == trset->model[k - 1])
          fprintf(txt, "%c", str[k]);
      else
          fprintf(txt, " %c", str[k]);          /* separation entre 2 atomes */
  }
  fprintf(txt, "\n");

  return;
}
/*=============================================================================
 fPrintfSubTrset(): Imprime dans le fichier 'txt' les donnees d'un fragment de
                    la structure pointee par 'trset' commencant a 'bgn' et de
                    longueur 'len'.
=============================================================================*/

void fPrintfSubTrset(Trset *trset, int bgn, int len, FILE *txt)
{
  int  i; 
  char *strg;

  strg = GapsStr(trset);

  for (i = 0; i < trset->digits; i++) {
      fprintf(txt, "     ");
      fPrintfAtoms(trset->omodel[i], trset, bgn, len, txt);
  }

  fprintf(txt, "     ");
  fPrintfAtoms(strg, trset, bgn, len, txt);

  if (trset->voidcols != 0)
      fprintf(txt, "     warning: %d column%s found void (#) in trset data\n",
              trset->voidcols, trset->voidcols > 1 ? "s" : "");

  fprintf(txt, "\n");

  for (i = 0; i < trset->nseq; i++) {
      fprintf(txt, "%-4d ", i+1);
      fPrintfAtoms(trset->data[i], trset, bgn, len, txt);
  }
  fprintf(txt, "\n");

  free(strg);

  return;
}
/*=============================================================================
 fPrintfTrset(): Imprime dans le fichier 'txt' l'ensemble des donnees de la
                 structure pointee par 'trset'.
=============================================================================*/

void fPrintfTrset(Trset *trset, FILE *txt)
{
  fPrintfSubTrset(trset, 0, trset->len, txt);
  return;
}
/*=============================================================================
 fPrintfTrsetStat(): Affiche les resultats de la statistique faite sur une base
                   d'entrainement.
=============================================================================*/

void fPrintfTrsetStat(char *filename, int nseq, int len, FILE *fout)
{

  fprintf(fout, "\nTraining set:\t\"%s\":\n\t\t%d sequence%s of length %d\n",
          filename, nseq, (nseq > 1 ? "s" : ""), len);

  return;
}
/*=============================================================================
 GetTrsetPatternId(): En retour 'id' pointe la chaine identifiant le 'pattern'
                      constitue de la structure complete pointee par 'trset'.
                      'id' doit disposer d'un volume suffisant.
=============================================================================*/

void GetTrsetPatternId(Trset *trset, char *id)
{
  sprintf(id, "%s%d%c",
         (trset->atlist[0] < 0 ? "-" : ""), trset->oatlist[0], ',');
  sprintf(id, "%s%d", id, trset->oatlist[trset->natom - 1]);
  
  return;
}
/*=============================================================================
 StripTrset(): Cree une copie du fichier pointe par 'filename' contenant un
               alignement multiple, en supprimant dans la copie les colonnes de
               l'alignement ne contenant que des gaps, les sequences codant la
               structure secondaire sont aussi modifiees.
               L'argument 'trset' pointe la structure supposee prealablement
               lue dans le fichier source.
               Au fichier cree est attribue le nom pointe par 'outputname'.
               La fonction retourne 1 si le fichier est cree, sinon 0.
=============================================================================*/

int StripTrset(char *filename, Trset *trset, char *outputname)
{
  FILE  *input, *output;
  char  *gapstr, *str;
  int   i, buffer_sz = TRSET_MAX_LEN;

  gapstr = GapsStr(trset);                    /* analyse des gaps de 'trset' */

  if (trset->voidcols == 0)  return 0;                       /* rien a faire */

  if ((input = fopen(filename, "r")) == NULL) {
      fprintf(stderr, "%s: file not found, exit.. \n", filename);
      exit(1);
  }
  str = (char *) malloc(buffer_sz);

  output = fopen(outputname, "w");

  while (fgets(str, buffer_sz, input) != NULL)
  {
      if (str[0] == '>')
          fprintf(output, "%s", str);    /* copie d'une ligne de commentaire */
      else
      if (DataLen(str, strlen(str)) != 0)
      {
          for (i = 0; i < trset->len; i++)
              if (gapstr[i] != '#')  fputc(str[i], output);
          fputc('\n', output);
      }
  }

  trset->len -= trset->voidcols;    /* actualisation de la structure 'trset' */
  trset->voidcols = 0;

  free(gapstr);
  free(str);
  fclose(input); 
  fclose(output); 

  return 1;
}
/*===========================================================================*/
