
/*=============================================================================
 tstrip.c                                       A.Lambert  le  19/06/02

 Procede au nettoyage d'une base d'entrainement des colonnes vides dans l'ali-
 gnement multiple.
 Un fichier contenant la nouvelle base d'entrainement est cree qui, par defaut,
 porte le meme nom augmente du suffixe '.strip' ou le nom optionnellement entre
 en 2eme argument du programme.

 make -f apps.mk tstrip ;
 mv tstrip ../bin ;

 tstrip <trset> [<ouputname>]

 tstrip ~/devc/projets/bio/data/trsetsII/rnaseP.db rnaseP.db2

=============================================================================*/

#include "rnaIV.h"

void  TstripHelp(void);

/*===========================================================================*/

int main(int argc, char *argv[])
{
  Trset   *trset;
  char    *ouputname;
  int     status;

  ReadHelpArgs(argc, argv, 1, TstripHelp);
  fTest(argv[1]);
  trset = ReadTrset(argv[1], 'S', stdout);

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
/*===========================================================================*/

