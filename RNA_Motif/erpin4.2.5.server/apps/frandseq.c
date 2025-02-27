
/*=============================================================================
 frandseq.c                           A.Lambert le 21/02/03 revu le 20/04/04

 Generation de sequences aleatoires de nucleotides dont les bases sont distri-
 buees suivant les pourcentages entres en argument, par defaut: 25.

 La "semence" du generateur de nombres pseudo-aleatoires, la longueur et le
 nombre de sequences a creer peuvent etre entres en argument au programme.
 Les valeurs par defaut sont respectivement 1, 600 et 1.
 Le programme utilise la sortie standard, sous un format de type 'fasta'.
 Les sequences sont directement creees dans un fichier, sans passer par l'in-
 termediaire d'un tableau.
 Le programme n'utilise que la bibliotheque standard du C.

 cc -O3 -Wall -o frandseq frandseq.c ;
 strip frandseq ;
 chmod 755 frandseq ;
 mv frandseq ../bin;

 usage:   tous les arguments sont optionnels
 frandseq  [-len <length>]           longueur des sequences creees, defaut: 600
           [-nseq <nseq>]                 nombre de sequences creees, defaut: 1
	   [-ATG <%A> <%T> <%G> [%C]]        pourcentages des bases, defaut: 25
	   [-seed <seed>]   initialisation du gen. de nb. aleatoires, defaut: 1

 defaut:   frandseq -len 600 -nseq 1 -ATG 25 25 25 -seed 1
 exemples: frandseq -len 1000000 > ~/devc/projets/bio/data/sequences/rnd1M.fna
           frandseq -len 1200 -nseq 4 > randseqs.fna
           frandseq -nseq 2 -ATG 20 40 15 -seed 12

=============================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void fRandSeq(FILE *txt, int len, double *freqs, int cpl);
void fPrintCmt(FILE *txt, int no, int length, double *freqs);

/*===========================================================================*/

int main(int argc, char *argv[])
{
  double sum, freqs[4] = {0.25, 0.25, 0.25, 0.25};
  int    i, seed = 1, len = 600, nseq = 1;

  for (i = 1; i < argc; i++)
  {
      if (strcmp(argv[i], "-len") == 0)   len  = atoi(argv[++i]);
      if (strcmp(argv[i], "-nseq") == 0)  nseq = atoi(argv[++i]);
      if (strcmp(argv[i], "-seed") == 0)  seed = atoi(argv[++i]);
      if (strcmp(argv[i], "-ATG") == 0) {
          if (argc - i - 1 < 3) {
              fprintf(stderr, "'ATG' option needs 3 arguments, exit..\n");
              exit(1);
          }
          freqs[0] = atof(argv[++i])/100.;                  /* freq. des 'A' */
          freqs[1] = atof(argv[++i])/100.;                  /* freq. des 'T' */
          freqs[2] = atof(argv[++i])/100.;                  /* freq. des 'G' */
      }
  }
  for (i = 0, sum = 0.0; i < 3; i++) {
      sum += freqs[i];
      if (freqs[i] < 0.0) {
          fprintf(stderr, "arguments of 'ATG' option must be >= 0, exit..\n");
          exit(1);
      }
  }
  if (sum > 1.) {
      fprintf(stderr, "invalid arguments of 'ATG' option, exit..\n");
      exit(1);
  }
  srand(seed);

  for (i = 0; i < nseq; i++)  {
      fPrintCmt(stdout, i+1, len, freqs);
      fRandSeq(stdout, len, freqs, 60);
  }
  return 0;
}
/*=============================================================================
 fRandSeq(): Genere une sequence aleatoire de longueur 'len', dont la composi-
             tion en bases ATGC est issue des elements du tableau 'freqs'.
             La sequence est creee dans le fichier pointe par 'txt', sous un
             format de 'cpl' caracteres par ligne.
=============================================================================*/

void fRandSeq(FILE *txt, int len, double *freqs, int cpl)
{
  int    i;
  double rnd, freqA, freqAT, freqATG;

  freqA = freqs[0];                                         /* freq. des 'A' */
  freqAT = freqA + freqs[1];                              /* freq. des 'A+T' */
  freqATG = freqAT + freqs[2];                          /* freq. des 'A+T+G' */

  for (i = 0; i < len;  )
  {
      rnd = rand() / (RAND_MAX + 1.0);

      if (rnd < freqA)   fputc('A', txt);
      else
      if (rnd < freqAT)  fputc('T', txt);
      else
      if (rnd < freqATG) fputc('G', txt);
      else
          fputc('C', txt);

      if ((++i) % cpl == 0) fputc('\n', txt);           /* retour a la ligne */
  }
  fputc('\n', txt);                               /* retour a la ligne final */
}
/*=============================================================================
 fPrintCmt(): Enregistre dans le fichier 'txt' un commentaire de type 'fasta'
              pour une sequence de longueur 'length'.
	      L'indice 'no' est joint ainsi que les frequences des bases.
=============================================================================*/

void fPrintCmt(FILE *txt, int no, int length, double *freqs)
{
  freqs[3] = 1.0 - (freqs[0] + freqs[1] + freqs[2]);

  fprintf(txt, ">random seq %d: %d bases, ATGC-freqs: %.2f %.2f %.2f %.2f\n",
          no, length, freqs[0], freqs[1], freqs[2], freqs[3]);
}
/*===========================================================================*/
