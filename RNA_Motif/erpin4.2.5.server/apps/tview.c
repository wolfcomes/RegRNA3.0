
/*=============================================================================
 tview.c                        A.Lambert le 03/10/01  revu le  16/04/04

 Visualisation d'un training-set ou d'une region de celui-ci.

 make -f apps.mk tview ;
 mv tview ../bin ;

 tview <trset> [<pattern>]

 tview  ~/devc/projets/bio/data/trsetsII/trna.db  -2,4 | more
 tview  ~/devc/projets/bio/data/trsetsII/sarcin.db -4,4
 tview  ~/devc/projets/bio/data/trsetsII/23S301proPK.db
 tview  ~/devc/projets/bio/data/trsetsII/uniq+forts-100-100.epn > \
        ~/devc/projets/bio/erpin5/jobs/uniq+forts-100-100.dat
 tview  ~/devc/projets/bio/data/trsetsII/rnaseP.db 67,14 | more
 tview  ~/devc/projets/bio/data/trsetsII/rnaseP.db.strip 67,14 | more
 tview  ~/devc/projets/bio/data/trsetsII/mirna1.epn -2,2
 tview  ~/devc/projets/bio/data/trsetsII/mirna2.epn -2,2
 tview  ~/devc/projets/bio/data/trsetsII/trna-typeI.epn -2,2

=============================================================================*/

#include "rnaIV.h"

void  get_geom(char *pattern_id, Trset *trset, int *bgn, int *len);
void  TviewHelp(void);

/*===========================================================================*/

int main(int argc, char *argv[])
{
  Trset   *trset;
  int     bgn, len;
  char    id[32];

  ReadHelpArgs(argc, argv, 1, TviewHelp);
  fTest(argv[1]);
  trset = ReadTrset(argv[1], 'V', stdout);

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
/*===========================================================================*/
