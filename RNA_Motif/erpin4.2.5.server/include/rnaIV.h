
/*=============================================================================
 rnaIV.h            A.Lambert le 17/09/01      derniere revision le 15/04/04
=============================================================================*/

#ifndef __RNAIV__
#define __RNAIV__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

/********************* macros a modifier si necessaire ***********************/


#define TRSET_MAX_LEN    12000    /* longueur maximale de l'alignement d'ARN */
#define COMMENT_MAX_LEN  300     /* longueur max. d'une ligne de commentaire */
#define HLX_MAX_LEN      64                /* longueur maximale d'une helice */

#define SEQ_MAX_LEN      300000000   /* borne sup. a la long. d'une sequence */
#define SEQ_MAX_NB       10000000       /* borne sup. au nombre de sequences */
#define SEQ_SECURE       10       /* securite: octets ajoutes en fin de seq. */
                                  /* est mis a 0 lors des tests de fonctions */
#define SEARCH_EXT       0     /* proplongement d'un intervalle de recherche */
                                           /* d'un masque a ses 2 extremites */

#define CFG_MAX_BITS     24.0     /* nombre maximal de bits admis pour coder */
                                 /* le nombre des configurations d'un masque */
                                          /* au dela de ce nombre -> Exit    */
#define CFG_BITS_WARN    16.0             /* au dela de ce nombre -> Warning */
#define CFG_BITS_INFO    10.0             /* au dela de ce nombre -> Info    */

#define HISTO_FNAME  "epnhist.dat"    /* fichier de l'histogramme des scores */
#define MIN_DETECTS  10    /* nbre minimum de detections pour creer l'histo. */

#define PSEUDO_COUNTS_FACTOR 2.e-3            /* le poids des pseudo-comptes */
#define PSEUDO_COUNTS_USER   0.1          /* par defaut est le produit des 2 */


                   /* nombre d'echantillons pour la 1ere colonne des profils */
#define SAMPLES1_SINGLE 300                  /* dans le cas d'un simple brin */
#define SAMPLES1_DOUBLE 150               /* et dans le cas d'un double brin */

#define DELTA_H   0.05 /* mesure de la division elementaire des histogrammes */


/*****************************************************************************/

#define  ON   1                                             /* usage general */
#define  OFF  0
                               /* trois types de fichiers 'FASTA' consideres */
#define TRSET   0                /* base d'entrainement: modele + alignement */
#define DATASET 1                        /* ensemble de sequences a analyser */
#define UNKNOWN 2                /* fichier au contenu non gerable: erreur ? */

#define DYNAMIC 0                /* gestion dynamique des masques successifs */
#define STATIC  1                 /* gestion statique des masques successifs */

#define FORWARD       0                                       /* brin direct */
#define REV_CMPL      1                               /* brin complementaire */
#define FWD_AND_REV   2            /* signale le traitement dans les 2 brins */

#define LONG    0                      /* version longue des infos affichees */
#define SHORT   1                                     /* version plus courte */
#define MUTE    2                                         /* aucun affichage */

#define GLOBAL  0            /* statistique globale des poids et transitions */
#define LOCAL   1                      /* statistique locale (seq. par seq.) */
#define UNIFORM 2                                     /*statistique uniforme */

#define PERCENT 0        /* 2 types d'entree des seuils: entree de la valeur */
#define VALUE   1              /* ou du pourcentage de captures dans 'trset' */

#define HLX     0                                 /* type des atomes: helice */
#define HLX1    1                                       /* 1er brin d'helice */
#define HLX2    2                                      /* 2eme brin d'helice */
#define STD     3                                                  /* strand */

#define MASK    0                     /* 2 facons complementaires de masquer */
#define UMASK   1                  /* des elements d'une structure 'Pattern' */
#define NOMASK  2                  /* ici le masque ne cache rien du pattern */
#define ADDMASK 3              /* ici on ajoute des elts au masque precedent */

#define ALPHA_LEN       4                  /* longueur de l'alphabet utilise */
#define SQR_ALPHA_LEN   16
#define _A_             0                   /* indices des lettres utilisees */
#define _T_             1
#define _G_             2
#define _C_             3
#define _X_             4                               /* concerne les gaps */
#define _N_             5                       /* nucleotide non identifiee */

#define ALPHA_LENxA     0
#define ALPHA_LENxT     4
#define ALPHA_LENxG     8
#define ALPHA_LENxC     12

#define MAX(a, b)   ((a) > (b) ? (a) : (b))
#define MIN(a, b)   ((a) < (b) ? (a) : (b))
#define LOG2(x)     (log((double) x) / log(2.0))

extern double       LOG_ZERO;                  /* borne inferieure de log(x) */
extern int          ScoreTabLen;          /* longueur des tableaux de scores */
extern unsigned long long int TotalNucScans;     /* total des bases visitees */
extern short  *NtStCode;                           /* codage des nucleotides */
extern short  **NtHlxCode;                                /* voir 'ntcode.c' */

/*-----------------------------------------------------------------------------
 Point: Description d'un point du plan
-----------------------------------------------------------------------------*/

typedef struct {        /* structure servant a stocker les coord. d'un point */
   double x, y;
} Point;

/*-----------------------------------------------------------------------------
 Map: structure destinee a construire l'etalonnage d'un axe
-----------------------------------------------------------------------------*/

typedef struct {
  double  min, max, ratio;
  int     Pixmin, Pixmax;
  } Map;

/*-----------------------------------------------------------------------------
 Histo: description simplifiee d'un histogramme
-----------------------------------------------------------------------------*/

typedef struct {
  double  hmin, hmax;
  int     bins, samples;
  double  *vals;
  } Histo;

Histo  MainMaskEvals;             /* histogramme des E-values des detections */
                                                      /* du masque principal */
Histo  MainMaskDetects;    /* histogramme des detections du masque principal */
Map    MainMaskMap;       /* pour l'etalonnage de l'axe de 'MainMaskDetects' */

/*-----------------------------------------------------------------------------
  Config: description des configurations d'un 'Pattern' ou d'un 'Mask',
          en general, sera construit comme un champ de ces structures et
          utilisera les donnees qui y sont enregistrees
-----------------------------------------------------------------------------*/

typedef struct {
    short len,                               /* longueur de la configuration */
          *hxbgn, *hxgaps, *hxdist,        /* helices: tabl. de taille 'nhx' */
                                   /* debut, gaps et dist. entre les 2 brins */
          *stbgn, *stgaps, *stlen;           /* brins: tabl. de taille 'nst' */
                                        /* debut, gaps et longueur des brins */
  }
  Config;

/*-----------------------------------------------------------------------------
  List: description d'une liste pour la collecte de resultats d'une recherche.
-----------------------------------------------------------------------------*/

typedef struct {
      int    offset;  /* position d'une detection depuis le debut de la seq. */
      int    cfgindex;                /* indice de la configuration detectee */
      double score;                    /* score de la configuration detectee */
      Config *cfg;                   /* contenu de la configuration detectee */
  }
  Data;

struct lklist {
      Data          data;
      struct lklist *next;
  };

typedef struct lklist List;

/*-----------------------------------------------------------------------------
  Trset: description des donnees d'une base d'entrainement, constituee d'un
         alignement multiple de sequences.
-----------------------------------------------------------------------------*/

typedef struct {
    int   nseq,                       /* nombre de sequences de l'alignement */
          len,                                     /* longueur des sequences */
          voidcols,               /* nombre de colonnes vides (que des gaps) */
          digits,                /* nbre de caracteres codant chaque colonne */
          natom, nhx, nst;       /* nbre d'atomes, d'helices et brins isoles */

    int   *model;                 /* tableau (len) codant la struct. second. */
    int   *atlist, *hxlist, *stlist;  /* listes des atomes, helices et brins */
    int   *oatlist, *ohxlist, *ostlist;  /* les meme dans le codage original */
    char  **omodel,                 /* codage original de la struct. second. */
          **data;                       /* sequences: tableau nseq * (len+1) */

    double **hlxsum, **stsum;      /* pointent les matrices de substitution  */
    double hpcw, spcw;                           /* poids des pseudo-comptes */
    int    pcflag;            /* signal de l'introduction des pseudo-comptes */
  }
  Trset;

/*-----------------------------------------------------------------------------
  Atom: description des brins elementaires d'une base d'entrainement.
-----------------------------------------------------------------------------*/

typedef struct {
    int   id,                                              /* identificateur */
          type,                              /* valeurs: HLX1 ou HLX2 ou STD */
          index,                 /* index dans les listes d'helices ou brins */
          h2index,        /* en cas d'helice, indice de l'atome du 2eme brin */
          min_len, max_len, max_gaps,                 /* param. geometriques */
          db_bgn,               /* abscisse du debut dans la base de donnees */
          min_bgn,                         /* abscisse du debut, gaps exclus */
                         /* en general repere depuis le debut d'un 'pattern' */
          bgn, len, gaps;          /* variables a l'interieur d'un 'pattern' */
  }
  Atom;

/*-----------------------------------------------------------------------------
  Helix: description des helices, contituees de 2 brins sans gaps
-----------------------------------------------------------------------------*/

typedef struct {
    int     id,                    /* identificateur dans le modele original */
            helix_len, min_len, max_len,
            min_dist, max_dist, max_gaps,
            db_bgn1, db_bgn2,      /* debuts reperes dans la base de donnees */
            min_bgn,                       /* abscisse du debut, gaps exclus */
                                    /* repere depuis le debut d'un 'pattern' */
            bgn, dist;           /* variables: debut de l'helice et distance */
                                                      /* entre les 2 parties */
    double  score,
            **Profile;
    float   **Scores,        /* pointeur sur un tab. des scores d'1 sequence */
            **ScoresBis;                    /* 2eme tab. completant 'Scores' */
  }
  Helix;

/*-----------------------------------------------------------------------------
  Strand: description de simples brins contenant ou non des gaps.
-----------------------------------------------------------------------------*/

typedef struct {
    int     id,                    /* identificateur dans le modele original */
            min_len, max_len, max_gaps,
            db_bgn,             /* abscisse du debut dans la base de donnees */
            min_bgn,                       /* abscisse du debut, gaps exclus */
                                    /* repere depuis le debut d'un 'pattern' */
            bgn, len;                /* variables: debut et longueur du brin */
    char    *str;                                   /* chaine du brin aligne */
    double  score,
            **Profile, **Align;
    float   **Scores,        /* pointeur sur un tab. des scores d'1 sequence */
            **ScoresBis;                    /* 2eme tab. completant 'Scores' */
  }
  Strand;

/*-----------------------------------------------------------------------------
 Pattern: description d'un motif complexe, constitue d'une sequence CONNEXE de
          structures 'Atom'.
-----------------------------------------------------------------------------*/

typedef struct {
    char  *id;                                /* identificateur du 'pattern' */
    int    bgnid, endid,               /* identif. des 1er et dernier atomes */
           trsetbgnindex,        /* index DANS TRSET du 1er atome de pattern */
           natom,                         /* nombre d'elts 'Atom' du pattern */
           nhx, nst,              /* nombres d'helices et strands du pattern */
           nvarst,                /* nombre de strands variables (avec gaps) */
           ncfg;                 /* nombre total de configurations possibles */
    double log2ncfg;           /* log2 de 'ncfg', pour prevenir les overflow */

    int    *model,                      /* codage de la structure secondaire */
           *atlist, *hxlist, *stlist, /* listes des atomes, helices et brins */
           *oatlist, *ohxlist, *ostlist; /* les meme dans le codage original */

    Atom   *atom;      /* pointeur sur le tableau d'elts ordonnes du pattern */
    Helix  *hlx;           /* pointeur sur le tableau des helices du pattern */
    Strand *std;           /* pointeur sur le tableau des strands du pattern */

    short  ***gapslist,        /* pointeurs sur les tableaux des nombres min */
           **gapslist_ori;           /* et max des gaps dans les 'nst' brins */

    int    db_bgn,      /* absc. du debut de pattern dans la base de donnees */
           min_len,                /* somme des longueurs minimum des atomes */
           max_gaps,              /* somme des nbres de gaps max. des atomes */
           max_len;
  }
  Pattern;

/*-----------------------------------------------------------------------------
  Mask: Description d'un 'masque', sous-ensemble d'un meme 'pattern',
        constitue de certains des ses brins et helices completes
-----------------------------------------------------------------------------*/

typedef struct {
    int      nargs, mode,         /* nombre et mode de saisie des arguments: */
                                           /* MASK, UMASK, ADDMASK ou NOMASK */
             *args,                  /* pointeur sur les arguments du masque */
             nhx, nst, ncfg,    /* nombre d'helices, brins et configurations */
             natom,                                       /* nombre d'atomes */
             *hxindex, *stindex;      /* pointeurs sur les indices d'helices */
    double   log2ncfg;         /* log2 de 'ncfg', pour prevenir les overflow */
                                            /* et brins du 'pattern' associe */
    Atom     *atom;                   /* pointeur sur le 1er atome du masque */
    char     *atomstr;        /* chaine de '1' et '0' suivant que les atomes */
                                           /* appartiennent ou non au masque */
    short    **gapslist,      /* tableau des min. et max. de gaps des atomes */
             **gapslist_ori,      /* le meme, destine a restaurer 'gapslist' */
             **gapscfg;                /* tableau des configurations de gaps */
    int      ngaps;                     /* nb. de colonnes de ces 2 tableaux */
                                  /* du masque et du complementaire, et dim. */
    Pattern  *pattern;             /* structure 'pattern' associee au masque */
    Config   *cfg;              /* pointeur sur les configurations du masque */
                                        /* parametres lies a une exploration */
    double   score, threshold;          /* enregistrement des score et seuil */
    int      cfgindex;          /* indice d'1 config., par ex. de score max. */
    char     *str;        /* pointe une seq. reconstituee depuis une config. */

    int      db_bgn,     /* absc. du debut du masque dans la base de donnees */
             min_bgn,    /* somme des long. minimales des elts qui precedent */
             max_bgn,        /* abscisse maximale du debut, dans le pattern  */
             max_len,                     /* longueur maximale avec les gaps */
             min_len;                 /* ecart minimal des elements extremes */
  }
  Mask;

/*-----------------------------------------------------------------------------
  Sequence: description des donnees et attributs d'une sequence
-----------------------------------------------------------------------------*/

typedef struct {
    char    *comment;                  /* commentaire attache a une sequence */
    int     nb,                     /* nomero d'ordre dans le fichier source */
            datalen,                          /* volume original des donnees */
            tail, len;              /* rajout et volume total de la sequence */
    char    *data;                          /* tableau des donnees A,T,G,C.. */
  }
  Sequence;

/*-----------------------------------------------------------------------------
  Context: description des parametres, attributs d'un process de recherche
           d'un masque en 1 ou plusieurs etapes emboitees, dans une sequence.
-----------------------------------------------------------------------------*/

typedef struct {
    int     level,                        /* niveau de la recherche: 0, 1,.. */
                          /* pour l'utilisateur 'level' est decalle: 1, 2,.. */
            bgn,                             /* debut dans la seq. originale */
            range,             /* longueur consideree dans la seq. originale */
            datalen,                          /* long. de sequence a visiter */
            scoretablen,                   /* longueur des tableaux de score */
            tabscanlen,            /* longueur "utile" des tableaux de score */
            seqscanlen,             /* longueur "utile" de  sequence visitee */
            nscan,        /* nombre de cycles complets de tableaux de scores */
            lastscanlen,                     /* longueur de la partie finale */
            lastscoretabLen;                    /* longueur du tableau final */
    char    *data;                    /* pointeur sur la 1ere base a visiter */

    Mask    *mask;             /* pointeur sur le masque associe a 'Context' */
    Sequence *seq;          /* pointeur sur la sequence associee a 'Context' */

    int     dir,              /* FORWARD ou REV_CMPL: brin direct ou inverse */
            style,                   /* sortie des infos: LONG SHORT ou MUTE */
            hist,       /* ON ou OFF: sortie de l'histogramme des detections */
            eval,          /* ON ou OFF: sortie de la E-value des detections */
            warnings,               /* ON ou OFF: sortie ou non des messages */
            detects,            /* compteur des detections du masque associe */
            cumuls,                      /* compteur des detections cumulees */
            timer;                /* ON ou OFF: chronometrage des recherches */
    unsigned long  bgntime;                        /* debut du chronometrage */
    double  duration;                                               /* duree */
    List    *list;                    /* liste pour la collecte de resultats */
    FILE    *outfile;             /* fichier de sortie: stdout, /dev/null .. */
  }
  Context;

/*-----------------------------------------------------------------------------
  Threshold: description de l'entree d'un seuil pour une operation de recherche.
-----------------------------------------------------------------------------*/

typedef struct {
      int    input;                       /* type d'entree: VALUE ou PERCENT */
      int    percent;             /* pourcentage de captures dans un 'trset' */
      double val;               /* valeur du seuil entre ou apres convertion */
  }
  Threshold;

/*=============================================================================
 fonctions de 'env.c', environnement
=============================================================================*/

void  ChScoreTabLen(int val);
void  ChLogZero(double val);

/*=============================================================================
 fonctions de 'io.c', entrees/sorties.
=============================================================================*/

void   fPrintfSeq(FILE *txt, char *seq, int length, int cpl);
void   fPrintfdMat(double **m, int nrow, int ncol, FILE *txt);
void   fPrintfiMat(int **m, int nrow, int ncol, FILE *txt);
double CpuTime(unsigned long *bgntime);
void   fPrintfCpuTime(double duration, FILE *txt);
void   GetPathTail(char *path, char **tail);
void   fTest(char *fname);
double StrToD(char *str, char *status);
void   ioError(char *func_name, char c);

/*=============================================================================
 fonctions de 'tab1.c', gestion de tableaux unidimensionnels.
=============================================================================*/

int    DataLen(char *str, int len);
int    TabLen(const int *t);
int    *TabSearch(const int *t, int val);
int    *TabSearchn(int *t, int val, int p);
int    TabRepeats(int *t, int n);

void   FilldTab(double *tab, int len, double val);
void   FilliTab(int *tab, int len, int val);
void   FillStr(char *str, int len, char val);

double *padd(double val, int n);

/*=============================================================================
 fonctions de 'tab2.c', gestion de tableaux bidimensionnels.
=============================================================================*/

char   **cMat(int nrow, int ncol);
short  **sMat(int nrow, int ncol);
int    **iMat(int nrow, int ncol);
double **dMat(int nrow, int ncol);
float  **fMat(int nrow, int ncol);
void   FreecMat(char **m);
void   FreesMat(short **m);
void   FreeiMat(int **m);
void   FreedMat(double **m);
void   FreefMat(float **m);
void   FillcMat(char **m, int nrow, int ncol, char val);
void   FillsMat(short **m, int nrow, int ncol, short val);
void   FilliMat(int **m, int nrow, int ncol, int val);
void   FilldMat(double **m, int nrow, int ncol, double val);
void   FillfMat(float **m, int nrow, int ncol, float val);
void   CopysMat(short **src, int nrow, int ncol, short **target);

/*=============================================================================
 fonctions de 'Seqs.c': manipulations de sequences.
=============================================================================*/

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
 fonctions de 'trset.c' concernant une base d'entrainement.
=============================================================================*/

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
 fonctions de 'sstat.c' concernant la statistique des bases dans des sequences.
=============================================================================*/

void sGetStat(char *seq, int len);
void InitStatTables(void);
void ReInitStatTables(double *freqs);
void SwapStatTables(void);
void ReverseStat(void);
char *MakeRandSeq(int len);
double *Cumuls(double *freqs);
void SetRandSeq(char *seq, int len, double *cumuls);
char *RandSeq(int len, double *freqs);

void ChStStat(Strand *St);
void ChHlxStat(Helix *Hlx);
void ChPatternStat(Pattern *pattern);

/*=============================================================================
 fonctions de 'fstat.c' concernant la statistique des bases dans un fichier.
=============================================================================*/

int  fGetStat1(FILE *fd, int *cur_pos, double *freqs, int pr);
int  fGetStat(char *data_name, int seqnb1, int nseq, int mode, FILE *fout);
void PrintfStat(int nseq, int sum, int filetype, FILE *fout);

int    fGetShortStat1(FILE *fd, int *cur_pos, int *length, int proc);
double fGetShortStat(char *data_name, int seqnb1, int nseq, int bgn, int range,
                     int masklen, int mode, FILE *fout);
int    fGetLongStat1(FILE *fd, int *cur_pos, double *freqs, int *length,
                     int proc);
double fGetLongStat(char *data_name, int seqnb1, int nseq, int bgn, int range,
                    int masklen, int mode, FILE *fout);
void   PrintfLongStat(int nseq, double nucs_to_scan, FILE *fout);
unsigned int fGetFreqs(char *data, double *freqs);

/*=============================================================================
 fonctions de 'atom.c'
=============================================================================*/

Atom  *ReadAtoms(int *atomlist, Trset *trset);

/*=============================================================================
 fonctions de 'pattern.c'
=============================================================================*/

Pattern *NewPattern(void);
void    FreePattern(Pattern *pattern);
void    DelPattern(Pattern *pattern);
void    ReadPatternId(char *pattern_id, int *bgn, int *end);

Pattern *ReadPatternAtoms(char *pattern_id, Trset *trset);
Pattern *ReadPattern(char *pattern_id, Trset *trset);

void    sPrintfPattern(char *str, Pattern *pattern);
void    sPrintfPatternStruct(char *str, Pattern *pattern);

/*=============================================================================
 fonctions de 'helix.c'.
=============================================================================*/

Helix  *SetHlxTable(int nhx);
void   FreeHlxTable(Helix *Hlx, int nhx);
void   DelHlxTable(Helix *Hlx, int nhx);
void   ReadHelix(Helix *Hlx, Pattern *pattern, int hlx_id);
void   ReadHelices(Pattern *pattern);

/*=============================================================================
 fonctions de 'strand.c'.
=============================================================================*/

Strand *SetStdTable(int nst);
void   FreeStdTable(Strand *Std, int nst);
void   DelStdTable(Strand *Std, int nst);
void   ReadStrand(Strand *Std, Pattern *pattern, int std_id);
void   ReadStrands(Pattern *pattern);

/*=============================================================================
 fonctions de 'profs.c': construction des profils de brins et helices.
=============================================================================*/

double **GetWeights(char **data, int nbstr, int bgn, int len);
void   GetStWeights(Trset *trset, Strand *St);

double **GetCorrels(char **data, int nbstr, int bgn1, int bgn2, int len);

void   GetStProfile(Trset *trset, Strand *St);
void   GetHlxProfile(Trset *trset, Helix *Hlx);
void   GetPatternProfiles(Trset *trset, Pattern *pattern);

/*=============================================================================
 fonctions de 'sum.c': introduction des matrices de substitution.
=============================================================================*/

void ReadSUM(char *fname, double hpcw, double spcw, Trset *trset);
void SumImg(double **Sbstm, int dim, double **prof, int len, double weight);

/*=============================================================================
 fonctions de 'profsSM.c': profils utilisant une matrice de substitution.
=============================================================================*/

double **GetWeightsSM(char **data, int nbstr, int bgn, int len,     /* brins */
                     double **Sbstm, double pcw);
void   GetStWeightsSM(Trset *trset, Strand *St);
                                                                  /* helices */
double **GetCorrelsSM(char **data, int nbstr, int bgn1, int bgn2, int len,
                      double **Sbstm, double pcw);
                                                               /* interfaces */
void   GetStProfileSM(Trset *trset, Strand *St);
void   GetHlxProfileSM(Trset *trset, Helix *Hlx);
void   GetPatternProfilesSM(Trset *trset, Pattern *pattern);
void   GetMaskProfilesSM(Trset *trset, Mask *mask);

void   FreePatternProfiles(Pattern *pattern);
void   FreeMaskProfiles(Mask *mask);

/*=============================================================================
 fonctions de 'scores.c': calcul de scores de brins et helices.
=============================================================================*/

double GetStNoGScore(char *seq, Strand *St);
void   GetHlxScores(char *seq, Helix *Hlx, double *scores);
double GetHlxBestScore(char *seq, Helix *Hlx);

void   GetStScores(char *seq, Strand *St, double *scores);
double GetStBestScore(char *seq, Strand *St);

void   GetHlxScoresTab(char *seq, int length, Helix *Hlx);
void   GetStScoresTab(char *seq, int length, Strand *St);
void   GetHlxScoresTabBis(char *seq, int length, Helix *Hlx);
void   GetStScoresTabBis(char *seq, int length, Strand *St);

/*=============================================================================
 fonctions de 'align.c': alignement d'un brin sur un profil.
=============================================================================*/

void   AlignSProfile(char *seq, Strand *St);
void   AlignBack(char *seq, Strand *St);

/*=============================================================================
 fonctions de 'cfg.c': configurations d'un masque.
=============================================================================*/

void SetCfgTable(Mask *mask);
void FreeCfgFields(Config *cfg);
void DelCfgTable(Mask *mask);
void SetGapCfgTable(Mask *mask);
void DelGapCfgTable(Mask *mask);

Config *CloneCfg(Mask *mask, Config *cfg);

/*=============================================================================
 fonctions de 'mask.c': gestion d'un masque de pattern.
=============================================================================*/

Mask  *NewMask(void);
void  FreeMask(Mask *mask);
void  DelMask(Mask *mask);
void  SetUMask(Mask *mask, Pattern *pattern);
void  CtrlMaskErrors(Mask *mask, Pattern *pattern);
void  GetMaskGeom(Mask *mask);

void  GetNoMask(Mask *mask, Pattern *pattern);
void  GetMask(Mask *mask, Pattern *pattern);

/*=============================================================================
 fonctions de 'maskcfg.c': configurations d'un masque.
=============================================================================*/

void  GetMaskAtoms(Mask *mask);
void  GetMaskAtomStr(Mask *mask);
void  GetMaskGapList(Mask *mask);
void  ResetMaskGapList(Mask *mask);

void  EnumGapCfgs(short **mgaplist, short nvarst, short col, int *nb,
                  short *tmpgaps, short **gapscfglist);
int   GetGapsCfgs(Mask *mask);
void  PrtGapsCfg(Mask *mask, int cfgindex, FILE *txt);
void  PrtGapsCfgs(Mask *mask, FILE *txt);

void  GetCfg(Config *cfg, Mask *mask, short *gaps_cfg);
void  GetMaskCfgs(Mask *mask);

int   GetMaskCfgLen(Config *cfg, Mask *mask);
void  GetMaskCfgsLens(Mask *mask);

/*=============================================================================
 fonctions de 'masks.c': gestion d'un tableau de masques.
=============================================================================*/

Mask *NewMasks(int nmask);
void FreeMasks(Mask *mask, int nmask);
void DelMasks(Mask *mask, int nmask);
int  Read1MaskArgs(int argc, char *argv[], int *index, Mask *mask);

Mask *ReadMasksArgs(int argc, char *argv[], int *nmask);
void ParseMasksArgs(Mask *mask, int nmask, Pattern *pattern);

void GetMasks(Mask *mask, int nmask, Pattern *pattern);
void GetMasksCfgs(Mask *mask, int nmask);
void SetMasksScoresTabs(Context *ctxt, int nctxt);

/*=============================================================================
 fonctions de 'mscores.c': calcul de scores d'un masque de pattern.
=============================================================================*/

void  GetMaskProfiles(Trset *trset, Mask *mask);
void  SetMaskScoresTab(Mask *mask, int len, int level);
void  FreeMaskScoresTab(Mask *mask);
void  GetMaskScoresTab(char *seq, int len, Mask *mask, int level);
float GetMaskScore(Mask *mask, int tab_index, int level);

/*=============================================================================
 fonctions de 'tscores.c': calcul et statistique des scores d'un 'trset'.
=============================================================================*/

double GetHlxTScore(char *seq, Helix *Hlx);
double GetStTScore(char *seq, Strand *St);
double GetPatternTScore(char *seq, Pattern *pattern);
double GetMaskTScore(char *seq, Mask *mask);

double *GetPatternTScores(Trset *trset, Pattern *pattern);
double *GetMaskTScores(Trset *trset, Mask *mask);

int    cmpdbl(const void *x1, const void *x2);
void   fPrintfCaptures(double *ts_scores, int nscores, FILE *txt);
double ConvertRatio(int percent, double *ts_scores, int nscores);
void   fPrintfScoresStat(double *scores, int nscores, FILE *txt);

void   fPrintfPatternTStat(Trset *trset, Pattern *pattern, FILE *txt);
void   fPrintfMaskTStat(Trset *trset, Mask *mask, FILE *txt);
double GetPatternThreshold(int percent, Trset *trset, Pattern *pattern);
double GetMaskThreshold(int percent, Trset *trset, Mask *mask);
void   GetMasksThresholds(int *percent, Trset *trset, Mask *mask, int nmask);

/*=============================================================================
 fonctions de 'cfgstr.c': sortie de sequence depuis la configuration.
=============================================================================*/

void RecordMaskSeq(Mask *mask, char *seq, Config *cfg);
void CharInsert(char *str, int pos, char insert);
void ToLower(char *str, int len);
void GetCfgAtoms(Config *cfg, Mask *mask);

/*=============================================================================
 fonctions de 'outputs.c': sortie des resultats d'une recherche de masque.
=============================================================================*/

void PrintOutput(Context *ctxt, int detect, Config *cfg, double score);
void fPrintfOutput(Context *ctxt, int detect);
void fPrintfData(Context *ctxt, Data *data);

/*=============================================================================
 fonctions de 'thresholds.c': entree des seuils.
=============================================================================*/

Threshold *GetThresholdsArgs(int argc, char *argv[], int nmask);
void  GetThresholds(Threshold *threshold, Trset *trset, Mask *mask, int nmask);
void  fPrintfThresholds(Mask *mask, int nmask, FILE *txt);

/*=============================================================================
 fonctions de 'msearch.c': gestion des 'Context' et recherche de masques.
=============================================================================*/

Context *NewContext(void);
Context *SetupContext(Mask *mask, int nctxt);
void  InitSeqContext(Context *ctxt, Sequence *seq);
void  InitSeqContexts(Context *ctxt, int nctxt, Sequence *seq);

void  InitOutputContext(Context *ctxt, int style, int hist, int eval,
                        int warnings, int chrono, FILE *outfile);
void  InitOutputContexts(Context *ctxt, int nmask, int style, int hist, int eval,
                         int warnings, int chrono,  FILE *outfile);
void  SetDir(Context *ctxt, int nmask, int dir);
void  ResetDetect(Context *ctxt, int nmask);
void  StartTimer(Context *ctxt);
void  GetTime(Context *ctxt);
void  PrintTime(Context *ctxt);
void  ResetTimer(Context *ctxt, int nmask);

void  GetBgn(Context *ctxt, int detect, int *bgn);
void  StoreBgn(Context *ctxt, int bgn);
void  StoreInterval(Context *ctxt, int bgn, int range);
int   SetSeqContext(Context *ctxt);

int   MaskSearch(Context *ctxt, int levels);

void  StartStatus(void);
void  PrintStatus(Context *ctxt, unsigned int *nscans);
void  ClearStatus(void);
void  fPrintfProcessLog(Context *ctxt, int nctxt, int total);

/*=============================================================================
 fonctions de 'args.c': saisie d'arguments entres aux programmes.
=============================================================================*/

void ReadHelpArgs(int argc, char *argv[], int argcmin, void (*helpfunc)(void));
void GetOutputArgs(int argc, char *argv[], int *style, int *warnings,
                   int *evflag, int *histflag);
void GetDirArg(int argc, char *argv[], int *dir);
void GetStatArg(int argc, char *argv[], int *stat);
void GetSeqArgs(int argc, char *argv[],
                       int *seqnb1, int *nseq, int *bgn, int *range);
void GetEnvArgs(int argc, char *argv[], int *chrono);
void GetProcArg(int argc, char *argv[], int *maskproc);
char *GetSumArgs(int argc, char *argv[], int *pcflag, double *hpcw, double *spcw);

/*=============================================================================
 fonctions de 'list.c': liste pour la collecte des resultats
=============================================================================*/

void AddLink(Context *ctxt, int offset, int cfgindex, double score);
void AddLink2(Context *ctxt, int offset, int cfgindex, double score);
void FreeList(Context *ctxt);
Data *CopyList(Context *ctxt, int *links);
void OverlapFilter(Context *ctxt);
void OverlapFilter2(Context *ctxt);
int  fPrintfPostProcess(Context *ctxt);
int  fPrintfPostProcess2(Context *ctxt);

/*=============================================================================
 fonctions de 'dmp.c': gestion dynamique des configurations de masques
=============================================================================*/

void  GetPatternGapList(Pattern *pattern, int nmask);
void  FreePatternGapList(Pattern *pattern, int nmask);
void  ResetPatternGapList(Pattern *pattern);
void  ResetGapLists(Context *ctxt, int bgn, int end);
void  ExportMaskCfgInfos(Context *ctxt, short *mcfg);
void  ImportPatternCfgInfos(Context *ctxt);
void  ImportPatternCfgInfos2(Context *ctxt);

int   RecMaskSearch(Context *ctxt, int levels);
int   Search(Context *ctxt, int nmask, int maskproc);

/*=============================================================================
 fonctions de 'ctrlcfgs.c': gestion et controle du nbre de configs. des masques
=============================================================================*/

void  GetPatternCfgsNb(Pattern *pattern);
void  GetMaskCfgsNb(Mask *mask);
int   MaskCfgsNb(Mask *mask);
void  PrintMaskCfgsVol(Mask *mask, FILE *txt);
void  CtrlCfgs(Mask *mask, int nmask, int maskproc, char mode, FILE *txt);
void  CfgMsg(int maskproc, FILE *txt);
int   Rnd(int min, int max);
short *GetRandCfg(Mask *mask);

/*=============================================================================
 fonctions de 'maps.c': correspondance entre nombres reels et entiers.
=============================================================================*/

Map    MapSetup(double x1, double x2, int X1, int X2);
double Getx(Map map, int X);
int    GetX(Map map, double x);
int    GetXFloor(Map map, double x);
int    GetXCeil(Map map, double x);

/*=============================================================================
 fonctions de 'dhisto.c': histogramme des scores des detections.
=============================================================================*/

void   InitDetectsHisto(Context *ctxt);
void   AddToDetectsHisto(double x);
void   PrintHisto(Histo hist, char *filename);
Histo  ReadHisto(char *filename);
void   PrintScoresHisto(Context *ctxt);

/*=============================================================================
 fonctions de 'histools.c': histogrammes de scores.
=============================================================================*/

Histo  SetupHist(double min, double max, double dh, int samples);
Histo  GetNHist(float *scores, int samples, double dx);
Histo  GetNHistW(double *scores, double *weights, int samples, double dx);
void   NormalizeHist(Histo hist);
double *SetHistAxis(Histo hist);
void   PrHistInfo(FILE *txt, char *cmt, Histo hist);
void   PrHistDat(FILE *txt, Histo hist);
void   GetHistStat(Histo hist, double *mean, double *std);
double Gauss(double x, double mean, double std);
int    RandNucl(double p0, double p1, double p2);
char   RndNucl(double p0, double p1, double p2);
void   SortAsc(int *T, int n);
int    ResetHistBins(Histo *hist, int bins);

/*=============================================================================
 fonctions de 'hshisto.c': histogrammes de scores de brins et helices.
=============================================================================*/

double **CatStProfiles(Mask *mask, int *width);
double **CatStNoGProfs(Mask *mask, int *width);
double **CatHlxProfiles(Mask *mask, int *width);

Histo  StsNoGHist(double **Prof, int width, double dh, double *Prfsc);
Histo  GetStNoGHist(Strand *St, double dh, double *Prfsc);

double **GetWFScoresProb(Strand *St, double *bkgfreqs);
float  *RndAlnScores(Strand *St, int samples, double *Prfsc);
Histo  GetAlnHist(Strand *St, double dx, double *Pfs, int samples);

void   GetSSProfileStat(double **SSprof, int width, double *Prfsc,
                        double *min, double *max, double *mean, double *std);
Histo  GetSSHist(double **SSProf, int width, double dh, double *Prfsc);

/*=============================================================================
 fonctions de 'mhisto.c': histogramme d'helices, brins et masques.
=============================================================================*/

void   GetSSStat(Mask *mask, double *pfsc, double *min, double *max,
                 double *mean, double *std);
void   GetTScoresStat(double *scores, int nscores,
                      double *min, double *max , double *mean , double *std);
void   GetMaskTStat(Trset *trset, Mask *mask,
                    double *min, double *max , double *mean , double *std);

Histo  GetHlxHist(Mask *mask, double dx, double *Prfsc);
Histo  GetStHist(Mask *mask, double dx, double *Prfsc, int spercol);
Histo  GetStsHist(Mask *mask, double dx, double *Prfsc, int spercol);
Histo  GetMaskHist(Mask *mask, double dx, double *Prfsc, int spercol);

/*=============================================================================
 fonctions de 'conv.c': produit de convolution.
=============================================================================*/

double *convlv(double *a, double *b, int na, int nb, int *ncnv);
void   conv(double **a, int *na, double *b, int nb);
void   ConvHist(Histo *h1, Histo h2);

/*=============================================================================
 fonctions de 'cdf.c': "cumulative distribution function"
=============================================================================*/

Histo  GetHistCdf(Histo hist);
double Interpol(Histo cdf, double x);

/*=============================================================================
 fonctions de 'lfit.c': fit lineaire par les moindres carres.
=============================================================================*/

void fit(Point *Pts, int ndata, double *a, double *b,
         double *siga, double *sigb, double *chi2);
void lfit(Point *Pts, int ndata, double *a, double *b);

/*=============================================================================
 fonctions de 'bkgstat.c': statistique des bases de sequences aleatoires.
=============================================================================*/

void ReadBkgFreqs(double *freqs, int offset, char *argv[]);
void ResetBkgFreqs(double *freqs);
void ChMaskStat(Mask *mask);

/*=============================================================================
 fonctions de 'Eval.c': calcul de la E-value des scores.
=============================================================================*/

Histo  GetEvals1(Mask *mask, double datavol, int smplc1);
Histo  GetEvals(Mask *mask, double datavol, int dir, int stat);
void   GetEvalues(Mask *mask, double datavol, int dir, int stat);
double Evalue(double x);
void   FreeEvalsTab(void);
double evprob(double p, unsigned int N);
void   fPrintfEvalue(FILE *fout, double E, double cutoff, int dir, double datavol);

/*=============================================================================
 fonctions de 'ntcode.c': tableaux de codage des nucleotides.
=============================================================================*/

void  SetNtStCode(void);
void  SetNtHlxCode(void);
void  SetNtCodes(void);
void  FreeNtCodes(void);

/*=============================================================================
=============================================================================*/

#endif
