
/*=============================================================================
 tstrip2.c                                    A.Lambert  le  19/06/02

 Procede au nettoyage d'une base d'entrainement des colonnes vides dans l'ali-
 gnement multiple.
 Un fichier contenant la nouvelle base d'entrainement est cree qui, par defaut
 porte le meme nom augmente du suffixe '.strip', ou le nom entre en 2eme argu-
 ment du programme.
 Dans cette version les gaps sont admis dans les helices.

 cc -O3 -Wall -o tstrip2 tstrip2.c -I../include -L../lib -lrnaIV -lm ;
 strip tstrip2 ;
 chmod 755 tstrip2 ;


 tstrip2 <trset> [<ouputname>]

 tstrip2 ~/devc/projets/bio/data/trsetsII/rnaseP.db rnaseP.db2

 tstrip2 ~/devc/projets/bio/data/16S.3.epn  16S.3.epn.strip

=============================================================================*/

#include "rnaIV.h"

void  TstripHelp(void);

Trset *ReadTrset2(char *filename, int mode, FILE *fout);

/*===========================================================================*/

int main(int argc, char *argv[])
{
  Trset   *trset;
  char    *ouputname;
  int     status;

  ReadHelpArgs(argc, argv, 1, TstripHelp);
  fTest(argv[1]);
  trset = ReadTrset2(argv[1], 'S', stdout);

  if (argc == 3) ouputname = argv[2];
  else
  {
      ouputname = (char *) malloc(strlen(argv[1]) + 10);
      sprintf(ouputname, "%s%s", argv[1], ".strip");  
  }

  status = StripTrset(argv[1], trset, ouputname);

  if (status == 1)
      fprintf(stdout, "file '%s' has been created..\n", ouputname);
  else
      fprintf(stdout, "file '%s' needs not changes..\n", argv[1]);

  DelTrset(trset);
  exit(0);
} 
/*=============================================================================
 TstripHelp(): help
=============================================================================*/

void TstripHelp(void)
{
  fprintf(stderr,

"tstrip: Delete void columns from a training set\n"
"Usage:\n"
"tstrip [-h]\n"
"tstrip <training-set>    training set file name\n"
"       [<ouputname>]     output file name, default: <training-set>.strip\n"
);

  return;
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
