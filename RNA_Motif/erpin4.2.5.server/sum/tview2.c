
/*=============================================================================
 tview2.c                    A.Lambert le 03/10/01  revu le  20/06/02

 Visualisation d'un training-set ou d'une region de celui-ci.
 Dans cette version les gaps sont admis dans les helices.

 cc -O3 -Wall -o tview2 tview2.c -I../include -L../lib -lrnaIV -lm ;
 chmod 755 tview2 ;
 strip tview2 ;

 tview2 <trset> [<pattern>]

 tview2  ~/devc/projets/bio/data/trsetsII/trna.db  -2,4 | more
 tview2  ~/devc/projets/bio/data/trsetsII/sarcin.db -4,4
 tview2  ~/devc/projets/bio/data/trsetsII/23S301proPK.db
 tview2  ~/devc/projets/bio/data/trsetsII/uniq+forts-100-100.epn > \
         ~/devc/projets/bio/erpin5/jobs/uniq+forts-100-100.dat
 tview2  ~/devc/projets/bio/data/trsetsII/rnaseP.db 67,14 | more
 tview2  ~/devc/projets/bio/data/trsetsII/rnaseP.db.strip 67,14 | more
 tview2  ~/devc/projets/bio/data/trsetsII/mirna1.epn -2,2
 tview2  ~/devc/projets/bio/data/trsetsII/mirna2.epn -2,2
 tview2  ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -2,2

 tview2 ~/devc/projets/bio/data/16S.3.epn > 16S.3.aln
 tview2 16S.3.epn.strip > 16S.3.aln


=============================================================================*/

#include "rnaIV.h"

void  get_geom(char *pattern_id, Trset *trset, int *bgn, int *len);
void  TviewHelp(void);

Trset *ReadTrset2(char *filename, int mode, FILE *fout);


/*===========================================================================*/

int main(int argc, char *argv[])
{
  Trset   *trset;
  int     bgn, len;
  char    id[32];

  ReadHelpArgs(argc, argv, 1, TviewHelp);
  fTest(argv[1]);
  trset = ReadTrset2(argv[1], 'V', stdout);

  if (argc < 3) GetTrsetPatternId(trset, id);
  else strcpy(id, argv[2]);

  get_geom(id, trset, &bgn, &len);

  fPrintfSubTrset(trset, bgn, len, stdout);

  DelTrset(trset);
  exit(0);
} 
/*=============================================================================
 get_geom(): Saisit dans le 'modele' de 'trset' les indices 'bgn' et 'len' d'un
             'pattern' identifie par 'pattern_id'.
=============================================================================*/

void get_geom(char *pattern_id, Trset *trset, int *bgn, int *len)
{
  Pattern  *P;
  char     *str = (char *) malloc(TRSET_MAX_LEN);
  
  P = ReadPatternAtoms(pattern_id, trset);
  *bgn = P->db_bgn;
  *len = P->max_len;
  sPrintfPattern(str, P);
  fprintf(stdout, "\n%s: %s\n\n", P->id, str);
  sPrintfPatternStruct(str, P);
  fprintf(stdout, "     %s\n", str);

  DelPattern(P);
  free(str);
  
  return;
}
/*=============================================================================
 TviewHelp(): help
=============================================================================*/

void TviewHelp(void)
{
  fprintf(stderr,

"tview: Training set VIEWer\n"
"Usage:\n"
"tview [-h]\n"
"tview <training-set> [<region>]\n"
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
