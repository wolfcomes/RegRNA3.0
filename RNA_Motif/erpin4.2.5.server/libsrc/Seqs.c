
/*=============================================================================
 Seqs.c                                   A.Lambert le 12/12/01

 fonctions impliquees dans des operations sur des sequence, generalement issues
 d'un fichier au format 'fasta', commode dans le sens ou ce format permet de
 conserver dans un meme fichier un nombre quelconque de sequences.

 Afin de conserver toujours presents les parametres fixant les conditions de
 l'operation et son contexte, on cree 1 nouvelle structure:
 - Seqence: donnees, parametres et attributs de la sequence cible.

 cc -O2 -Wall -c Seqs.c -I../include ;

 ar -rs ../librnaIV.a Seqs.o ;

=============================================================================*/

#include "rnaIV.h"


Sequence *NewSeq(void);
void     FreeSeqData(Sequence *seq);
void     FreeSeq(Sequence *seq);
void     DelSeq(Sequence *seq);
Sequence *InitSeq(int tailsize);
int      GetTail(Pattern *pattern);
void     AddTail(Sequence *seq, char car);
void     CutTail(Sequence *seq);
void     PreProcessSeq(Sequence *seq);
void     RevCmpl(Sequence *seq);

int      ReadSeq(FILE *src, Sequence *seq);
int      GetSeqPos(FILE *src, int *seqbgn, int *seqlen, char *comment);
void     CopySeq(FILE *src, int seqbgn, int bgn, int len, int cpl, FILE *target);

/*=============================================================================
 NewSeq(): Alloue une structure 'Sequence' et met a "NULL' les pointeurs,
           ulterieurement, seuls les champs "comment' et 'data' seront alloues.
=============================================================================*/

Sequence *NewSeq(void)
{
  Sequence *seq = (Sequence *) malloc(sizeof(Sequence));

  if (seq == NULL) {
      fprintf(stderr, "NewSeq: allocation failure, exit..\n");
      exit(1);
  }

  seq->comment = NULL;
  seq->data = NULL;

  return seq;
}
/*=============================================================================
 FreeSeqData(): Libere la memoire du champ 'data', le champ 'comment' pourra le
                plus souvent subsister. le champ 'data' est remis a 'NULL'.
=============================================================================*/

void FreeSeqData(Sequence *seq)
{
  if (seq->data != NULL) {
      free(seq->data);
      seq->data = NULL;
  }
  return;
}
/*=============================================================================
 FreeSeq(): Libere la memoire des champs 'data' et 'comment' et remet les poin-
            teurs a 'NULL'.
=============================================================================*/

void FreeSeq(Sequence *seq)
{
  if (seq->data != NULL) {
      free(seq->data);
      seq->data = NULL;
  }
  if (seq->comment != NULL) {
      free(seq->comment);
      seq->comment = NULL;
  }
  return;
}
/*=============================================================================
 DelSeq(): Detruit unestructure 'Sequence'.
=============================================================================*/

void DelSeq(Sequence *seq)
{
  FreeSeq(seq);
  free(seq);
  return;
}
/*=============================================================================
 InitSeq(): Avant le demarrage de la lecture de sequences le champ 'comment'
            est alloue, 'tail' initialise et le compteur mis a 0.
=============================================================================*/

Sequence *InitSeq(int tailsize)
{
  Sequence *seq;

  seq = NewSeq();
  seq->nb = 0;
  seq->tail = tailsize;
  seq->comment = (char *) malloc(COMMENT_MAX_LEN + 1);

  return seq;
}
/*=============================================================================
 GetTail(): determine depuis les donnees de 'pattern' le prolongement 'tail'
            qui sera charge dans le champ 'seq->tail' d'une sequence.
=============================================================================*/

int GetTail(Pattern *pattern)
{
  return (pattern->max_len - pattern->min_len + SEQ_SECURE); 
}
/*=============================================================================
 AddTail(): ajoute un groupe de caracteres 'car' a la fin d'une sequence (de
            puis le caractere '\0' en position 'seq->datalen') afin de permet-
            tre la terminaison correcte d'une operation de recherche.
=============================================================================*/

void AddTail(Sequence *seq, char car)
{
  int j;

  for (j = 0; j < seq->tail; j++)  seq->data[seq->datalen + j] = car; 
 
  seq->data[seq->datalen + seq->tail] = '\0';

  return;
}
/*=============================================================================
 CutTail(): supprime des caracteres (en general des 'E') a la fin d'une se-
	    quence en placant en tete de ce groupe le caractere '\0'.
=============================================================================*/

void CutTail(Sequence *seq)
{
  seq->data[seq->datalen] = '\0';
  return;
}
/*=============================================================================
 PreProcessSeq(): Traite une structure 'Sequence' prealablement aux operations
                  de recherche:
                  Convertion aux majuscules, changement des 'U' en 'T', toutes
                  les lettres differentes de ATGC sont changees en 'N'.
                  Ajoute a la fin des caracteres 'E' qui auront un traitement
                  approprie.
=============================================================================*/

void PreProcessSeq(Sequence *seq)
{
  int   i;
  char  c, *ptc;

  for (i = 0, ptc = seq->data; i < seq->datalen; i++, ptc++)
  {
      c = toupper(*ptc);
      if (c == 'U')  c = 'T';
      if (c != 'A' && c != 'G' && c != 'C' && c != 'T')  c = 'N';
      *ptc = c;
  }
  AddTail(seq, 'E');

  return;
}
/*=============================================================================
 RevCmpl(): Inverse, sur une longueur egale a 'seq->datalen' l'ordre des carac-
            teres d'une sequence et procede aux substitutions:
            A -> T, G -> C, C -> G et T -> A.
=============================================================================*/

void RevCmpl(Sequence *seq)
{
  int   i;
  char  c, *ptc, *ptc1, *ptc2;

  ptc = seq->data + seq->datalen - 1;

  for (ptc1 = seq->data, ptc2 = ptc;  ptc1 < ptc2;  ptc1++, ptc2--)
  {
      c = *ptc1;
      *ptc1 = *ptc2;
      *ptc2 = c;
  }

  for (i = 0, ptc1 = seq->data; i < seq->datalen; i++, ptc1++)
      switch(*ptc1) {
          case 'A': *ptc1 = 'T'; break;
          case 'T': *ptc1 = 'A'; break;
          case 'G': *ptc1 = 'C'; break;
          case 'C': *ptc1 = 'G'; break;
          default:  break;
      } 

  return;
}
/*=============================================================================
 ReadSeq(): Lit une structure 'Sequence' depuis un fichier 'ascii' au format 
            'fasta', variante de 'fReadSeq' dans 'seqs.c'.

            Dans le code la variable 'cur_pos' est declaree 'static' afin que
            les appels successifs utilisent la valeur fixee a la fin de l'appel
            precedent, c'est a dire la position finale du curseur de fichier.

            La derniere ligne de commentaire precedant la sequence est chargee
            dans le champ 'comment' de le structure pointee par 'seq'.
            Le champ 'tail' est suppose prealablement initialise (InitSeq).  
            La fonction retourne le nombre d'appels soit le No d'ordre de la
            sequence dans le fichier, ou 0 si le tableau n'a pas pu etre cree
            (EOF atteint, auquel cas la var.statique 'cur_pos' est remise a 0). 
=============================================================================*/

int ReadSeq(FILE *src, Sequence *seq)
{
  short  c, bgn; 
  char   *ptc;
  int    l, seq_bgn, seq_len, line_len;

  static int cur_pos = 0;                  /* position du curseur de fichier */
  static int call_nb = 0;                             /* compteur des appels */

  if (cur_pos < 0) {                   /* prend en compte un appel precedent */
      cur_pos = 0;                             /* arrivant en fin de fichier */
      return 0;
  }

  call_nb++;                            /* compteur des appels a la fonction */
  seq->nb++;                                      /* compteur  des sequences */

  fseek(src, cur_pos, SEEK_SET);         /* positionne le curseur de fichier */

  while ((c = getc(src)) != '>' && c != EOF) /* avance jusqu'au prochain '>' */
      /* rien */ ;
  if (c == EOF) {
      cur_pos = 0;
      return 0;
  }                                                 /* ou EOF, alors: sortie */


       /* --- analyse le bloc des commentaires precedant chaque sequence --- */


  bgn = c;
  while (bgn == '>')              /* conservera la derniere ligne rencontree */
  {                                                    /* commencant par '>' */
      seq->comment[0] = '>';
      ptc = seq->comment;
      line_len = 0;

      while ((c = getc(src)) != '\n')
          if (++line_len < COMMENT_MAX_LEN)  *(++ptc) = c;

      *(++ptc) = '\0';
  
      if ((c = getc(src)) == EOF) {
          cur_pos = 0;
          return 0;
      }
      bgn = c;
  } 
  fseek(src, -1, SEEK_CUR);                         /* recule d'un caractere */


                     /* ------ demarrage de la lecture de la sequence ------ */


  seq_bgn = ftell(src);    /* enregistre la position du debut de la sequence */

  seq_len = 0;                /* mise a 0 du compteur des car. alphabetiques */

  while ((c = getc(src)) != '>' && c != EOF)       /* compte les car. alpha. */
      if (isalpha(c)) seq_len++;                           /* de la sequence */ 

  cur_pos = (c != EOF ? ftell(src) - 1 : EOF);        /* position du curseur */

  fseek(src, seq_bgn, SEEK_SET);              /* retour en debut de sequence */
 

                       /* ----- enregistrement de la sequence analysee ----- */
                                                /* limite superieure imposee */

  if (seq_len > SEQ_MAX_LEN)
  {
      fprintf(stderr, "\nWarning:\n");
      fprintf(stderr, "in '%s'\n", seq->comment);
      fprintf(stderr, "sequence length exceeds limit of %d bases\n", SEQ_MAX_LEN);
      fprintf(stderr, "truncating sequence at this value..\n\n");

      seq_len = SEQ_MAX_LEN;
  }

  seq->datalen = seq_len;
  seq->len = seq_len + seq->tail;

  if ((seq->data = (char *) malloc(seq->len + 1)) == NULL) {
      fprintf(stderr, "ReadSeq: allocation failure, exit..\n");
      exit(1);
  } 

  ptc = seq->data;
  l = 0; 
  while ((c = getc(src)) != '>' && c != EOF)        /* jusqu'au prochain '>' */
      if (isalpha(c) && l++ < seq->datalen)
          *(ptc++) = c;

  seq->data[seq_len] = '\0';    /* terminaison par '\0' a la fin des donnees */

  return call_nb;
}
/*=============================================================================
 GetSeqPos(): Saisit les parametres geometriques d'une sequence dans un fichier.
              - Variante de 'ReadSeq' qui ne requiert pas de stockage a l'ex-
              ception du commentaire 'comment' saisi.

              'seqbgn' pointe en sortie sur la position du debut de la sequence
              concernee, de longueur '*seqlen'.

              La fonction retourne le nombre d'appels soit le No de la sequence
              ou 0 si le tableau n'a pas pu etre cree (EOF atteint), les va-
              leurs pointees par 'seqbgn' et 'seqlen' sont alors indeterminees.
=============================================================================*/

int GetSeqPos(FILE *src, int *seqbgn, int *seqlen, char *comment)
{
  short  c, bgn; 
  char   *ptc;
  int    line_len;

  static int cur_pos = 0;                  /* position du curseur de fichier */
  static int call_nb = 0;                             /* compteur des appels */

  if (cur_pos < 0) {                   /* prend en compte un appel precedent */
      cur_pos = 0;                             /* arrivant en fin de fichier */
      return 0;
  }

  call_nb++;                            /* compteur des appels a la fonction */
 
  fseek(src, cur_pos, SEEK_SET);         /* positionne le curseur de fichier */

  while ((c = getc(src)) != '>' && c != EOF) /* avance jusqu'au prochain '>' */
      /* rien */ ;
  if (c == EOF) {
      cur_pos = 0;
      return 0;
  }                                                 /* ou EOF, alors: sortie */


       /* --- analyse le bloc des commentaires precedant chaque sequence --- */


  bgn = c;
  while (bgn == '>')              /* conservera la derniere ligne rencontree */
  {                                                    /* commencant par '>' */     
      comment[0] = '>';
      ptc = comment;
      line_len = 0;

      while ((c = getc(src)) != '\n')
          if (++line_len < COMMENT_MAX_LEN)  *(++ptc) = c;

      *(++ptc) = '\0'; 
   
      if ((c = getc(src)) == EOF) {
          cur_pos = 0;
          return 0;
      }
      bgn = c;
  } 
  fseek(src, -1, SEEK_CUR);                         /* recule d'un caractere */


                     /* ------ demarrage de la lecture de la sequence ------ */


  *seqbgn = ftell(src);    /* enregistre la position du debut de la sequence */

  *seqlen = 0;                /* mise a 0 du compteur des car. alphabetiques */

  while ((c = getc(src)) != '>' && c != EOF)       /* compte les car. alpha. */
      if (isalpha(c)) (*seqlen)++;                         /* de la sequence */ 

  cur_pos = (c != EOF ? ftell(src) - 1 : EOF);        /* position du curseur */
                                                       /* en fin de sequence */
  return call_nb;
}
/*=============================================================================
 CopySeq(): Transfere tout ou partie d'une sequence reperee par 'GetSeqPos'
            d'un fichier dans un autre.
            'seqbgn' est la position du curseur de fichier dans 'src', pris
            comme point de depart dans la sequence a copier.
            'bgn' marque le debut de la lecture dans la sequence (0 pour le 1er
            caractere).
            'len' est la longueur lue dans 'src' et copiee dans 'target'.

            L'operation est stoppee si la fin du fichier est atteinte.
            L'ecriture dans 'target' se fait sous un format de 'cpl' caracteres
            par ligne.
=============================================================================*/

void CopySeq(FILE *src, int seqbgn, int bgn, int len, int cpl, FILE *target)
{
  short  c; 
  int    i = 0;

  fseek(src, seqbgn, SEEK_SET);          /* positionne le curseur de fichier */

  i = 0;
  while (i < bgn && (c = getc(src)) != EOF) {
      if (isalpha(c)) i++;           /* passe 'bgn' caracteres alphabetiques */
  }

  i = 0;
  while (i < len && (c = getc(src)) != EOF)        /* copie 'len' caracteres */
  {                                     /* alphabetiques de 'src' a 'target' */
      if (isalpha(c)) {      
          putc(c, target);
          if (++i % cpl == 0) putc('\n', target);          /* nouvelle ligne */
      }
  }
  putc('\n', target);

  return;
}
/*===========================================================================*/
