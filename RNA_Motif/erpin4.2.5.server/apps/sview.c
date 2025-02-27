
/*=============================================================================
 sview.c                                      A.Lambert  le 06/02/02


 Affiche tout ou partie d'une sequence extraite d'un fichier au format 'fasta'.

 La sequence est recherchee, soit par entree de son commentaire (prive ou non
 de son 1er caractere: '>' ainsi que des espaces qui peuvent suivre), soit par
 son numero, l'indexation partant de 1.
 Si le commentaire entre convient a plusieurs sequences (trop court !) seule la
 premiere occurence est retenue.

 Note: Les sequences sont lues par 'GetSeqPos' qui ne retient que les caracteres
       alphabetiques, supprimant les '-', les caracteres numeriques, les espaces
       et lignes vides.

 cc -O3 -Wall -o sview sview.c -I../include -L../lib -lrnaIV -lm ;
 strip sview ;
 chmod 755 sview ;
 mv sview ../bin ;

 sview   [-h]                                                            help
 sview   <inputfile>                     nom du fichier ou saisir la sequence
         [-nb <numero>]             indice de la sequence a saisir, debut a 1
         [-comment <comment>]          ou commentaire de la sequence a saisir
         [-bgn <indice du debut>]                 debut de la saisie, sinon 1
         [-len <longueur de la saisie>] longueur de la saisie, sinon complete

 Notes:  -nb, -comment    :   un seul de ces deux arguments peut etre present
         -n, -c, -b, -l   :    arguments raccourcis acceptes par le programme
         sview <inputfile>:     affiche le contenu entier de la 1ere sequence

 Exemples:
 sview  seqs.fna  -comment seq3
 sview  ~/devc/projets/bio/data/trsets/trna.db \
        -comment  "DA0340 TGC ARCHAEGLOBUS FULG. ARCHAE"
 sview  ~/devc/projets/bio/data/trsets/trna.db -nb 3 -b 10 -l 20
 sview  ~/devc/projets/bio/data/sequences/ecoli.fna -b 779000 -l 2000

=============================================================================*/

#include "rnaIV.h"

#define  SEARCH_NUMBER  0                           /* recherche par l'index */
#define  SEARCH_COMMENT 1                    /* recherche par le commentaire */
#define  CHAR_PER_LINE  70                 /* nombre de caracteres par ligne */

static void read_args(int argc, char *argv[]);
static char *StripFastaCmt(char *str, int *len);
static void usage(void);

static char *cmt;
static int  cmt_len, nb, bgn, length, search;

/*=============================================================================
 main():
=============================================================================*/

int main(int argc, char *argv[])
{
  char   *comment, *pc;
  int    seq_len, seq_bgn, cpl;
  FILE   *txt;
  int    nul, count, seq_found;

  read_args(argc, argv);

  txt = fopen(argv[1], "r");

  comment = (char *) malloc(COMMENT_MAX_LEN + 1);

  count = seq_found = 0;

  while (GetSeqPos(txt, &seq_bgn, &seq_len, comment))
  {

      count++;                                   /* commence a 1, non pas 0 */
      pc = StripFastaCmt(comment, &nul);

      if (
         (search == SEARCH_NUMBER && nb == count) ||
         (search == SEARCH_COMMENT && strncmp(pc, cmt, cmt_len) == 0)
         )
      {
          if (bgn > seq_len) {
              fprintf(stderr, "\nseq %d: %d bases, bgn = %d is out of range !\n",
                      count, seq_len, bgn);
          }
          else
          {
              length = bgn + length > seq_len + 1 ? seq_len - bgn + 1: length;

              if (bgn == 1 && length == seq_len)             
                  fprintf(stdout, "\n%s (seq %d: %d bases)\n",
                          comment, count, seq_len);
              else
                  fprintf(stdout, "\n%s (seq %d: Extract %d+%d)\n",
                          comment, count, bgn, length);

              cpl = (length <= 80 ? 80 : CHAR_PER_LINE);
              CopySeq(txt, seq_bgn, bgn - 1, length, cpl, stdout);

              seq_found = 1;
          }
          break;
      }
  }
                               /* informations sur les echecs de l'operation */

  if (search == SEARCH_NUMBER &&  nb > count)
  {
      char *tail;
      GetPathTail(argv[1], &tail);
      fprintf(stderr, "\n%d sequence%s read in '%s': Nb %d is out of range\n",
      count, (count > 1 ? "s" : ""), tail, nb);
      free(tail);
  }
  if (search == SEARCH_COMMENT && !seq_found)
  {
      char *tail;
      GetPathTail(argv[1], &tail);
      fprintf(stderr, "\n%d sequence%s read in '%s'\n",
              count, (count > 1 ? "s" : ""), tail);
      fprintf(stderr, "arg \"%s\" not found\n", cmt);    
      free(tail);
  }

  fclose(txt);
  free(comment);
  exit(0);
}
/*=============================================================================
 read_args(): lit et analyse les arguments entres a 'main'
=============================================================================*/

void read_args(int argc, char *argv[])
{
  int  i, test;
  char arg;
                             /* valeur par defaut des arguments du programme */

  nb = 1;                            /* le 1er indice est mis a 1, non pas 0 */
  bgn = 1;                           /* la saisie demarre du debut, fixe a 1 */
  length = SEQ_MAX_LEN;                      /* toute la sequence est saisie */
  search = SEARCH_NUMBER;                 /* pour etre coherent avec: nb = 1 */

  if (argc < 2 || strcmp(argv[1], "-h") == 0)  usage();

  fTest(argv[1]);                            /* teste la presence du fichier */


  for (i = 2, test = 0; i < argc; i++)    /* evite les args. contradictoires */
      if (argv[i][0] == '-') {
          arg = argv[i][1];
          if (arg == 'n' || arg == 'c') test++;
      }
  if (test > 1)  usage();


  for (i = 2; i < argc; i++)               /* fin de l'analyse des arguments */
  {
      if (argv[i][0] == '-')
      {
          arg = argv[i][1];
          switch (arg)
          {
              case 'n':                          /* -nb: recherche par index */
                  search = SEARCH_NUMBER;
                  nb = atoi(argv[++i]);
                  break;
              case 'c':             /* -comment: recherche par 'commentaire' */
                  search = SEARCH_COMMENT;
                  cmt = StripFastaCmt(argv[++i], &cmt_len);
                  break;
              case 'b':                          /* -bgn: debut de la saisie */
                  bgn = atoi(argv[++i]);
                  break;
              case 'l':                       /* -len: longueur de la saisie */
                  length = atoi(argv[++i]);
                  break;
          }
      }
  }
  bgn = (bgn <= 0 ? 1 : bgn);                     /* numerotation depuis 1 ! */
  nb  = (nb <= 0 ? 1 : nb);

  return;
}
/*=============================================================================
 StripFastaCmt(): Prepare la chaine de caracteres pointee par 'str' a la compa-
                  raison avec un commentaire de fichier 'Fasta', en selection-
                  nant sa partie "informative" (ignorant les 'espaces' en debut
                  et fin).
                  Retourne un pointeur sur la partie selectionnee de 'str' dont
                  la longueur est, en retour, pointee par 'len'. 
=============================================================================*/

char *StripFastaCmt(char *str, int *len)
{
  char *ptc = str;

  while (*ptc == '>' || isspace(*ptc)) ptc++;     /* elimine '\n', space ... */
                                                       /* en debut de chaine */
  *len = strlen(ptc);

  while (isspace(ptc[*len - 1])) (*len)--;          /* idem en fin de chaine */

  return ptc;
}
/*=============================================================================
 usage(): affiche l'usage du programme et provoque la sortie.
=============================================================================*/

void usage(void)
{
  fprintf(stderr,

"sview: Sequence VIEWer: prints to 'stdout' a sequence read in a 'fasta' file\n\n"

"Usage:\n"
"sview   [-h]                                                            help\n"
"sview   <inputfile>              name of file containing the target sequence\n"
"        [-nb <index>]              index of the target sequence, starts at 1\n"
"        [-comment <comment>]               or comment of the target sequence\n"
"        [-bgn <bgn of capture>]        begining of the capture, otherwise: 1\n"
"        [-len <length of capture>]     length of the capture, otherwise: end\n\n"

"Notes:  -nb, -comment    :         only one of these 2 arguments may be used\n"
"        -n, -c, -b, -l   :                short cuts accepted by the program\n"
"        sview <inputfile>:     prints the whole contents of the 1st sequence\n\n"

"Examples:\n\n"

"sview data/trsets/trna.db -comment \"DA0340 TGC ARCHAEGLOBUS FULG. ARCHAE\"\n"
"sview data/trsets/trna.db -nb 3 -bgn 10 -len 20\n"
"sview data/sequences/ecoli.fna -b 779000 -l 2000\n"
);

  exit(0);
}
/*===========================================================================*/




