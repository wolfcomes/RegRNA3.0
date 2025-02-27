
/*=============================================================================
 mksum.c                                        A.Lambert le 24/03/04

 Ce programme calcule les matrices de substitutions depuis une base d'entraine-
 ment dont les helices peuvent contenir des gaps.
 Il n'y a pas ici de ponderation des sequences.

 cc -O2 -Wall -c stsum.c -I../include ;
 cc -O2 -Wall -c hlxsum.c -I../include ;

 cc -O3 -Wall -o mksum mksum.c stsum.o hlxsum.o \
    -I../include -L../lib -lrnaIV -lm ;
 strip mksum ;
 chmod 755 mksum ;

 mksum <trset-name> <region> <mask>
       [-V]                            affichage des profils d'helices et brins

 mksum ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -2,2 -nomask

 tstrip2 ~/devc/projets/bio/data/16S.3.epn
 mksum ~/devc/projets/bio/data/16S.3.epn.strip -2,700 -nomask -V > mksum.dat

=============================================================================*/

#include "rnaIV.h"

void   PrdMat(FILE *txt, double **M, int nl, int nc);
void   PrhSUM(FILE *txt, double **M);
void   PrsSUM(FILE *txt, double **M);
double *SetSeqWeights(Trset *trset, double val);
Trset  *ReadTrset2(char *filename, int mode, FILE *fout);

double **HlxSuMat(Trset *trset, Mask *mask, double *sweights, int Vmode);
double **StSuMat(Trset *trset, Mask *mask, double *sweights, int Vmode);


/*===========================================================================*/

int main(int argc, char *argv[])
{
  Trset   *trset;
  Pattern *pattern;
  Mask    *mask;
  double  *sWeights;
  int     i, Vmode = 0, nmask;
  double  **SbstMat;
  FILE    *txt1, *txt2;

  if (argc < 4) {
      fprintf(stderr, "'%s' needs at least 3 arguments, exit..\n", argv[0]);
      exit(1);
  }
  fTest(argv[1]);

  for (i = 1; i < argc; i++)
      if (strcmp(argv[i], "-V") == 0)  Vmode = ON;

  mask = ReadMasksArgs(argc, argv, &nmask);
  if (nmask != 1) {
      fprintf(stderr, "'%s' needs 1 (and only 1) mask arg, exit..\n", argv[0]);
      exit(1);
  }

  trset = ReadTrset2(argv[1], 'S', stderr);
  pattern = ReadPattern(argv[2], trset);
  ParseMasksArgs(mask, nmask, pattern);
  sWeights = SetSeqWeights(trset, 1.0);
  GetMask(mask, pattern);

  txt1 = fopen("xSUM.dat", "w");
  txt2 = fopen("xSUM.txt", "w");

  if (mask->nhx != 0) {
      SbstMat = HlxSuMat(trset, mask, sWeights, Vmode);
      PrdMat(txt1, SbstMat, SQR_ALPHA_LEN, SQR_ALPHA_LEN);
      PrhSUM(txt2, SbstMat);
      FreedMat(SbstMat);
  }
  if (mask->nst != 0) {
      SbstMat = StSuMat(trset, mask, sWeights, Vmode);
      PrdMat(txt1, SbstMat, ALPHA_LEN, ALPHA_LEN);
      PrsSUM(txt2, SbstMat);
      FreedMat(SbstMat);
  }
  fclose(txt1);
  fclose(txt2);

  free(sWeights);
  DelTrset(trset);
  DelMasks(mask, nmask);
  DelPattern(pattern);

  return 0;
}
/*=============================================================================
 PrdMat(): Enregistre dans le fichier pointe par 'txt' le tableau 'M' de 'nl'
           lignes et 'nc' colonnes.
=============================================================================*/

void PrdMat(FILE *txt, double **M, int nl, int nc)
{
  int  i, j;

  fprintf(txt, "\n");

  for (i = 0; i < nl; i++)  {
      for (j = 0; j < nc; j++)   fprintf(txt, " %.3e", M[i][j]);
      fprintf(txt, "\n");
  }
  fprintf(txt, "\n");
}
/*=============================================================================
 PrhSUM(): Affiche dans le fichier pointe par 'txt' la matrice de substitution
           des helices, sous un format adapte a la lecture directe.
=============================================================================*/

void PrhSUM(FILE *txt, double **M)
{
  int  i, j, nl, nc;
  char AUGC[4] = "AUGC";

  nl = nc = SQR_ALPHA_LEN;

  fprintf(txt, "\n");

  for (j = 0; j < nc; j++)
      fprintf(txt, "       %c%c", AUGC[j/ALPHA_LEN], AUGC[j%ALPHA_LEN]);
  fprintf(txt, "\n");

  for (i = 0; i < nl; )
  {
      fprintf(txt, " %c%c ", AUGC[i/ALPHA_LEN], AUGC[i%ALPHA_LEN]);
      for (j = 0; j < nc; j++)   fprintf(txt, " %.1e ", M[i][j]);
      fprintf(txt, "\n");
      if ((++i) % ALPHA_LEN == 0) fprintf(txt, "\n");
  }
  fprintf(txt, "\n");
}
/*=============================================================================
 PrsSUM(): Affiche dans le fichier pointe par 'txt' la matrice de substitution
           des brins, sous un format adapte a la lecture directe.
=============================================================================*/

void PrsSUM(FILE *txt, double **M)
{
  int  i, j, nl, nc;
  char AUGC[4] = "AUGC";

  nl = nc = ALPHA_LEN;

  fprintf(txt, "\n");

  for (j = 0; j < nc; j++)
      fprintf(txt, "       %c ", AUGC[j%ALPHA_LEN]);
  fprintf(txt, "\n");

  for (i = 0; i < nl; i++)
  {
      fprintf(txt, " %c ", AUGC[i%ALPHA_LEN]);
      for (j = 0; j < nc; j++)   fprintf(txt, " %.1e ", M[i][j]);
      fprintf(txt, "\n");
  }
  fprintf(txt, "\n");
}
/*=============================================================================
 SetSeqWeights(): Cree le tableau des poids associes aux sequences de la base
                  d'entrainement pointe par 'trset', et le retourne initialise
		  uniformement a la valeurs 'val'.
=============================================================================*/

double *SetSeqWeights(Trset *trset, double val)
{
  double *sweights;
  int    i;

  if ((sweights = (double *) malloc(trset->nseq * sizeof(double))) == NULL) {
      fprintf(stderr, "SetSeqWeights: allocation failure, exit..\n");
      exit(1);
  }
  for (i = 0; i < trset->nseq; i++)  sweights[i] = val;

  return sweights;
}
/*=============================================================================
 ReadTrset2(): Interface des fonctions precedentes concernant la lecture des
              donnees d'une base d'entrainement et leur controle.
              En cas d'echec, si mode == 'V', les donnees de 'trset' sont diri-
              gees sur le fichier 'fout'.
=============================================================================*/

Trset *ReadTrset2(char *filename, int mode, FILE *fout)
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
/*
  CtrlTrsetStruct(trset, mode, fout);
*/
  ProcessTrset(trset);

  return trset;
}
