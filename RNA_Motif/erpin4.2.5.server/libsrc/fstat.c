
/*=============================================================================
 fstat.c                                       A.Lambert le 10/05/03

 Ce code concerne un ensemble de fonctions destinees a faire la statistique des
 bases dans des sequences, afin de completer les profils statistiques:

 - denombrer les frequences des bases A,T,G,C,

 Il est issu de 'sstat.c' dedouble pour en diminuer le volume.
 Il contient des fonctions concernant des mesures statistiques sur les sequen-
 ces de fichiers 'fasta'.

 On ecrit ici des variantes de 'fGetStat' afin d'adapter le calcul des E-value
 aux differentes options de 'erpin':
 '-globstat', '-locstat', '-unifstat' et aux options concernant le nombre et
 les portions de sequences a traiter.
 Le but est d'obtenir avant le debut des recherches le nombre le plus precis
 possible de nucleotides a visiter (ormis le double passage FWD&REV qui se re-
 duit a une multiplication par 2.

 cc -Wall -O2 -c fstat.c -I../include ;

 ar -rs ../lib/librnaIV.a fstat.o ;

=============================================================================*/

#include "rnaIV.h"

int  fGetStat1(FILE *fd, int *cur_pos, double *freqs, int pr);
int  fGetStat(char *data_name, int seqnb1, int nseq, int mode, FILE *fout);
void PrintfStat(int nseq, int sum, int filetype, FILE *fout);

int    fGetShortStat1(FILE *fd, int *cur_pos, int *length, int proc);
double fGetShortStat(char *data_name, int seqnb1, int nseq, int bgn, int range,
                     int masklen, int mode, FILE *fout);
int    fGetLongStat1(FILE *fd, int *cur_pos, double *freqs, int *length, int proc);
double fGetLongStat(char *data_name, int seqnb1, int nseq, int bgn, int range,
                    int masklen, int mode, FILE *fout);
void   PrintfLongStat(int nseq, double nucs_to_scan, FILE *fout);
unsigned int fGetFreqs(char *data, double *freqs);

extern double *NewLogDataFreqs, *LogDataFreqs;
extern double *LogDataFreqs1, *LogDataFreqs2;

/*=============================================================================
 fGetStat1(): Lit une sequence dans un fichier 'ascii' au format 'fasta' et
        opere le decompte des bases ATGC (ou AUGC, mais ignore les 'N') et les
        transitions entre bases, les sequences ne sont pas copiees.

        Cette fonction etant destinee a etre iteree: le tableau 'freqs' est
	incremente et initialise a l'exterieur.

        L'argument 'cur_pos' pointe la position du curseur de texte a partir de
        laquelle la lecture commence (mais pas forcement la lecture de la
        sequence qui, en general, est precedee de commentaires).

        Si l'argument 'pr' (process) vaut 0 la sequence est seulement lue.

        En sortie 'cur_pos' pointe la position finale du curseur de fichier
        (utilisee comme position de depart de la lecture suivante).
        La fonction retourne EOF (-1) si la fin du fichier a ete atteinte, ou
        1 ou 0 suivant que la sequence a ete traitee ou seulement lue.
=============================================================================*/

int fGetStat1(FILE *fd, int *cur_pos, double *freqs, int pr)
{
  short  c, bgn;

  fseek(fd, *cur_pos, SEEK_SET);         /* positionne le curseur de fichier */

  while ((c = getc(fd)) != '>' && c != EOF)  /* avance jusqu'au prochain '>' */
      /* rien */ ;

  if (c == EOF)  return EOF;                        /* ou EOF, alors: sortie */


         /* --- passe le bloc des commentaires precedant chaque sequence --- */


  bgn = c;
  while (bgn == '>')                  /* passe les lignes commencant par '>' */
  {
      while (getc(fd) != '\n')
          /* rien */ ;

      if ((c = getc(fd)) == EOF) return EOF;
      bgn = c;
  } 
  fseek(fd, -1, SEEK_CUR);                          /* recule d'un caractere */


                               /* ------ lecture seule de la sequence ------ */

  if (pr == 0) {
      while ((c = getc(fd)) != '>' && c != EOF)
          /* rien */ ;

      *cur_pos = (c != EOF ? ftell(fd) - 1 : EOF);    /* position du curseur */
      if (c == EOF)  return EOF;
      return 0;
  }
                       /* ------ lecture et traitement de la sequence ------ */

  while ((c = getc(fd)) != '>' && c != EOF)
  {
      switch(toupper(c))
      {
          case 'A': freqs[_A_]++; break;
          case 'G': freqs[_G_]++; break;
          case 'C': freqs[_C_]++; break;
          case 'U':
          case 'T': freqs[_T_]++; break;
          default :               break;
      }
  }
  *cur_pos = (c != EOF ? ftell(fd) - 1 : EOF);        /* position du curseur */

  if (c == EOF)  return EOF;
  return 1;
}
/*=============================================================================
 fGetStat(): Opere sur les sequences d'un fichier 'ascii' au format 'fasta' le
             decompte des bases ATGC (ou AUGC) et des transitions.
             Les sequences d'indice < 'seqnb1' (seqnb1 >= 1) sont ignorees:
             en cas de 'training-set' il faudra que: seqnb1 >= 2, une securite
             est prevue pour exclure la sequence d'annotation.
             'nseq' est le nombre de sequences a analyser.
             Les proportions de ATGC sont chargees dans le tab. 'NewLogDataFreqs'
             (initialises a 1 pour eviter les divisions par 0).
             Si mode == 'V' la statistique est affichee.

             Retourne le nombre total de caracteres alphabetiques decomptes dans
             les sequences.
             note: 'fGetStat' est destinee a interfacer la fonction 'fGetStat1'.
=============================================================================*/

int fGetStat(char *data_name, int seqnb1, int nseq, int mode, FILE *fout)
{
  FILE   *txt;
  int    i, j, filetype, cur_pos, count_seq;
  double sum;

  extern double *NewLogDataFreqs;

                         /* exclusion de la sequence de codage d'un 'treset' */
  if ((filetype = GetFileType(data_name)) == TRSET) 
      seqnb1 = MAX(2, seqnb1);

  if ((txt = fopen(data_name, "r")) == NULL) {
      fprintf(stderr, "%s: file not found, exit.. \n", data_name);
      exit(1);
  }
                          /* initialisation a 0.01 pour eviter les elts nuls */

  FilldTab(NewLogDataFreqs, ALPHA_LEN, 0.01);

                                      /* lecture et traitement des sequences */
  cur_pos = 0;

  for (i = j = 1; i < seqnb1 && j != EOF; i++)          /* lecture seulement */
      j = fGetStat1(txt, &cur_pos, NewLogDataFreqs, 0);

  for (i = 0; i < nseq && j != EOF; i++)            /* lecture et traitement */
      j = fGetStat1(txt, &cur_pos, NewLogDataFreqs, 1);

  count_seq = i;
  fclose(txt);
                                               /* fin du calcul des tableaux */

  for (i = 0, sum = 0.; i < ALPHA_LEN; i++)  sum += NewLogDataFreqs[i];

  for (i = 0; i < ALPHA_LEN; i++)  NewLogDataFreqs[i] /= sum;

  if (mode == 'V')
      PrintfStat(count_seq, sum, filetype, fout);
                                                  /* passage aux logarithmes */
  for (i = 0; i < ALPHA_LEN; i++)  NewLogDataFreqs[i] = log(NewLogDataFreqs[i]);

  return (int) sum;
}
/*=============================================================================
 PrintfStat(): Affiche les resultats de la statistique faite sur un fichier de
               donnees.
=============================================================================*/

void PrintfStat(int nseq, int sum, int filetype, FILE *fout)
{
  extern double *NewLogDataFreqs;

  if (filetype == DATASET)
      fprintf(fout, "\t\t%d sequence%s, %d nucleotides\n",
              nseq, nseq > 1 ? "s" : "", sum);

  fprintf(fout, "\t\tATGC ratios: %.3f  %.3f  %.3f  %.3f\n",
          NewLogDataFreqs[_A_], NewLogDataFreqs[_T_],
          NewLogDataFreqs[_G_], NewLogDataFreqs[_C_]
         );

  return;
}
/*=============================================================================
 fGetShortStat1(): variante allegee de 'fGetStat1()'
                   On ne calcule que le nombre de nucleotides d'une sequence,
                   pointe en retour par 'length' si 'proc' != 0.
=============================================================================*/

int fGetShortStat1(FILE *fd, int *cur_pos, int *length, int proc)
{
  short  c, bgn;

  *length = 0;
  fseek(fd, *cur_pos, SEEK_SET);         /* positionne le curseur de fichier */

  while ((c = getc(fd)) != '>' && c != EOF)  /* avance jusqu'au prochain '>' */
      /* rien */ ;

  if (c == EOF)  return EOF;                        /* ou EOF, alors: sortie */

         /* --- passe le bloc des commentaires precedant chaque sequence --- */

  bgn = c;
  while (bgn == '>')                  /* passe les lignes commencant par '>' */
  {
      while (getc(fd) != '\n')
          /* rien */ ;

      if ((c = getc(fd)) == EOF) return EOF;
      bgn = c;
  }
  fseek(fd, -1, SEEK_CUR);                          /* recule d'un caractere */

          /* ------ lecture seule d'une sequence exclue du traitement ------ */

  if (proc == 0)
  {
      while ((c = getc(fd)) != '>' && c != EOF)
          /* rien */ ;

      *cur_pos = (c != EOF ? ftell(fd) - 1 : EOF);    /* position du curseur */
      if (c == EOF)  return EOF;

      return 0;
  }

                /* ------ decompte des bases d'une sequence a traiter ------ */

  while ((c = getc(fd)) != '>' && c != EOF)
      if (isalpha(c))
          (*length)++;             /* selection des caracteres alphabetiques */

  *cur_pos = (c != EOF ? ftell(fd) - 1 : EOF);   /* positionement du curseur */
  if (c == EOF)  return EOF;

  return 1;
}
/*=============================================================================
 fGetShortStat(): Variante de fGetStat() adaptee a '-locstat' et '-unifstat'.
                  rappel: en reperage "interne" bgn = 0 si 1 a ete entre.
                  'range' est la longueur du fragment a visiter.
                  'masklen' est la longueur minimale du masque principal re-
                  cherche (le dernier de la liste), a deduire du decompte des
                  bases a explorer sur chaque sequence.
                  Cette version ne fait pas intervenir les var. 'LogDataFreqs'.
                  Retourne le nombre de nucleotides a visiter.
=============================================================================*/

double fGetShortStat(char *data_name, int seqnb1, int nseq, int bgn, int range,
                     int masklen, int mode, FILE *fout)
{
  FILE   *txt;
  int    i, j, cur_pos, count_seq, length, nucs;
  double nucs_to_scan;

                         /* exclusion de la sequence de codage d'un 'treset' */
  if (GetFileType(data_name) == TRSET)
      seqnb1 = MAX(2, seqnb1);

  if ((txt = fopen(data_name, "r")) == NULL) {
      fprintf(stderr, "%s: file not found, exit.. \n", data_name);
      exit(1);
  }

                                      /* lecture et traitement des sequences */
  cur_pos = 0;
  nucs = 0;
  nucs_to_scan = 0.;
                                     /* lecture seule des sequences ignorees */
  for (i = j = 1; i < seqnb1 && j != EOF; i++)          /* lecture seulement */
      j = fGetShortStat1(txt, &cur_pos, &length, 0);

  for (i = 0; i < nseq && j != EOF; i++)            /* lecture et traitement */
  {
      j = fGetShortStat1(txt, &cur_pos, &length, 1);
      nucs = MIN(range, (length - bgn));
      nucs -= (masklen-1);                       /* neutralisation de la fin */
      nucs = MAX(nucs, 0);                                         /* butoir */
      nucs_to_scan += (double)nucs;       /* total des nucleotides a visiter */
  }

  count_seq = i;
  fclose(txt);
  
  if (mode == 'V')
      fprintf(fout, "\t\t%lld nucleotides to be processed in %d sequence%s\n",
              (long long int)nucs_to_scan, count_seq, count_seq > 1 ? "s" : "");

  return nucs_to_scan;
}
/*=============================================================================
 fGetLongStat1(): Variante de 'fGetStat1' ou 'length' pointe en retour le nbre
                  de nucleotides de la sequence.
=============================================================================*/

int fGetLongStat1(FILE *fd, int *cur_pos, double *freqs, int *length, int proc)
{
  short  c, bgn;

  *length = 0;
  fseek(fd, *cur_pos, SEEK_SET);         /* positionne le curseur de fichier */

  while ((c = getc(fd)) != '>' && c != EOF)  /* avance jusqu'au prochain '>' */
      /* rien */ ;

  if (c == EOF)  return EOF;                        /* ou EOF, alors: sortie */


         /* --- passe le bloc des commentaires precedant chaque sequence --- */


  bgn = c;
  while (bgn == '>')                  /* passe les lignes commencant par '>' */
  {
      while (getc(fd) != '\n')
          /* rien */ ;

      if ((c = getc(fd)) == EOF) return EOF;
      bgn = c;
  }
  fseek(fd, -1, SEEK_CUR);                          /* recule d'un caractere */


                               /* ------ lecture seule de la sequence ------ */

  if (proc == 0)
  {
      while ((c = getc(fd)) != '>' && c != EOF)
          /* rien */ ;

      *cur_pos = (c != EOF ? ftell(fd) - 1 : EOF);    /* position du curseur */
      if (c == EOF)  return EOF;
      return 0;
  }
                       /* ------ lecture et traitement de la sequence ------ */
		       
  while ((c = getc(fd)) != '>' && c != EOF)
  {
      if (isalpha(c)) {            /* selection des caracteres alphabetiques */
          switch(toupper(c))
          {
              case 'A': freqs[_A_]++; break;
              case 'G': freqs[_G_]++; break;
              case 'C': freqs[_C_]++; break;
              case 'U':
              case 'T': freqs[_T_]++; break;
              default :               break;
          }
          (*length)++;                        /* decompte du total des bases */
      }
  }
  *cur_pos = (c != EOF ? ftell(fd) - 1 : EOF);        /* position du curseur */

  if (c == EOF)  return EOF;
  return 1;
}
/*=============================================================================
 fGetLongStat(): Variante de 'fGetStat' adaptee a l'option '-globstat' dont les
                 arguments sont ceux de 'fGetStat' et 'fGetShortStat'
                 Retourne le nombre total de nucleotides a visiter.
=============================================================================*/

double fGetLongStat(char *data_name, int seqnb1, int nseq, int bgn, int range,
                    int masklen, int mode, FILE *fout)
{
  FILE   *txt;
  int    i, j, cur_pos, count_seq, length, nucs;
  double sum, nucs_to_scan;

  extern double *NewLogDataFreqs;
                                  
                         /* exclusion de la sequence de codage d'un 'treset' */
  if (GetFileType(data_name) == TRSET)
      seqnb1 = MAX(2, seqnb1);

  if ((txt = fopen(data_name, "r")) == NULL) {
      fprintf(stderr, "%s: file not found, exit.. \n", data_name);
      exit(1);
  }

                          /* initialisation a 0.01 pour eviter les elts nuls */

  FilldTab(NewLogDataFreqs, ALPHA_LEN, 0.01);

                                      /* lecture et traitement des sequences */
  cur_pos = 0;
  nucs = 0;
  nucs_to_scan = 0.;

  for (i = j = 1; i < seqnb1 && j != EOF; i++)          /* lecture seulement */
      j = fGetLongStat1(txt, &cur_pos, NewLogDataFreqs, &length, 0);

  for (i = 0; i < nseq && j != EOF; i++)            /* lecture et traitement */
  {
      j = fGetLongStat1(txt, &cur_pos, NewLogDataFreqs, &length, 1);
      nucs = MIN(range, (length - bgn));
      nucs -= (masklen-1);                       /* neutralisation de la fin */
      nucs = MAX(nucs, 0);                                         /* butoir */
      nucs_to_scan += (double)nucs;       /* total des nucleotides a visiter */
  }
  count_seq = i;
  fclose(txt);
                                               /* fin du calcul des tableaux */

  for (i = 0, sum = 0.; i < ALPHA_LEN; i++)  sum += NewLogDataFreqs[i];

  for (i = 0; i < ALPHA_LEN; i++)  NewLogDataFreqs[i] /= sum;

  if (mode == 'V')
      PrintfLongStat(count_seq, nucs_to_scan, fout);
                                                  /* passage aux logarithmes */
  for (i = 0; i < ALPHA_LEN; i++)   NewLogDataFreqs[i] = log(NewLogDataFreqs[i]);

  return nucs_to_scan;
}
/*=============================================================================
 PrintfLongStat(): Affiche les resultats de la statistique faite sur un fichier
                   de donnees.
=============================================================================*/

void PrintfLongStat(int nseq, double nucs_to_scan, FILE *fout)
{
  extern double *NewLogDataFreqs;
  unsigned long long int N;

  N = (long long int) nucs_to_scan;

  fprintf(fout, "\t\t%lld nucleotides to be processed in %d sequence%s\n",
                N, nseq, nseq > 1 ? "s" : "");

  fprintf(fout, "\t\tATGC ratios: %.3f  %.3f  %.3f  %.3f\n",
          NewLogDataFreqs[_A_], NewLogDataFreqs[_T_],
          NewLogDataFreqs[_G_], NewLogDataFreqs[_C_]
         );

  return;
}
/*=============================================================================
 fGetFreqs(): Mesure dans le fichier 'fasta' dont le nom est pointe par 'data'
              les frequences des ATGC (ou AUGC) pointees en retour par 'freqs'
	      dont la memoire necessaire est supposee allouee.
	      Les autres caracteres (comme 'N') participent au decompte mais
	      pas au calcul des frequences.
	      La fonction retourne le nombre total de nucleotides lues.
=============================================================================*/

unsigned int fGetFreqs(char *data, double *freqs)
{
  FILE   *txt;
  int    i, c;
  double sum, nn = 0.;

  if ((txt = fopen(data, "r")) == NULL) {
      fprintf(stderr, "fGetFreqs: file '%s' not found, exit.. \n", data);
      exit(1);
  }
  FilldTab(freqs, ALPHA_LEN, 0.0);

  while ((c = getc(txt)) != EOF)
  {
      if (c == '>')
          while (getc(txt) != '\n')  ;             /* passe la ligne entiere */
      else
      switch(toupper(c)) {
          case 'A': freqs[0]++;  break;
          case 'G': freqs[1]++;  break;
          case 'C': freqs[2]++;  break;
          case 'U':
          case 'T': freqs[3]++;  break;
          default : if (isalpha(c)) nn++;  break;
      }
  }

  fclose(txt);

  for (i = 0, sum = 0.; i < ALPHA_LEN; i++) sum += freqs[i];
  for (i = 0; i < ALPHA_LEN; i++)  freqs[i] /= sum;

  return (unsigned int) (sum + nn);
}
/*===========================================================================*/
