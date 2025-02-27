
/*=============================================================================
 io.c                             A.Lambert  le 12/11/01

 Fonctions concernant les operations d'entree/sortie.

 cc -O2 -Wall -c io.c -I../include ;

 ar -r ../lib/librnaIV.a io.o ;

=============================================================================*/

#include "rnaIV.h"

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
 fPrintfSeq(): ecrit dans le fichier 'txt' une sequence entiere ou une partie
               'length' sous un format de 'cpl' caracteres par ligne.
=============================================================================*/

void fPrintfSeq(FILE *txt, char *seq, int length, int cpl)
{
  int  i, len = strlen(seq);

  len = length < len ? length : len; 

  fputc(seq[0], txt);          /* sorti de la boucle pour eviter le 1er '\n' */

  for (i = 1; i < len; i++)
  {
      if (i % cpl == 0) fputc('\n', txt);
      fputc(seq[i], txt);
  }
  fputc('\n', txt);

  return;
}
/*=============================================================================
 fPrintfdMat(): Enregistre le contenu de la matrice 'm' de doubles dans le 
                fichier 'txt'. ('stdout' pour un affichge au shell).
=============================================================================*/

void fPrintfdMat(double **m, int nrow, int ncol, FILE *txt)
{
  int    i, j;

  for (i = 0; i < nrow; i++)
  {
      for (j = 0; j < ncol; j++)   fprintf(txt, "%6.2f ", m[i][j]);
      fprintf(txt, "\n");
  }
  fprintf(txt, "\n");

  return;
}
/*=============================================================================
 fPrintfiMat(): Enregistre le contenu de la matrice 'm' d'entiers dans le 
                fichier 'txt'. ('stdout' pour un affichge au shell).
=============================================================================*/

void fPrintfiMat(int **m, int nrow, int ncol, FILE *txt)
{
  int    i, j;

  for (i = 0; i < nrow; i++)
  {
      for (j = 0; j < ncol; j++)   fprintf(txt, "%4d  ", m[i][j]);
      fprintf(txt, "\n");
  }
  fprintf(txt, "\n");

  return;
}
/*=============================================================================
 CpuTime(): retourne, converti en secondes le temps Cpu ecoule depuis la valeur
            pointee par 'bgntime', laquelle est reinitialisee.
=============================================================================*/

double CpuTime(unsigned long *bgntime)
{
  double          duration;
  unsigned long   endtime = clock();

  duration = (double) (endtime - *bgntime) / CLOCKS_PER_SEC;
  *bgntime = endtime;

  return duration;
}
/*=============================================================================
 fPrintfCpuTime(): affiche, converti en heures, minutes et secondes le temps 
                   Cpu entre, exprime en secondes.
=============================================================================*/

void fPrintfCpuTime(double duration, FILE *txt)
{
  int    hours, mins;
  double secs;

  fprintf(txt, "cpu-time: ");  

  hours = (int) duration / 3600;
  if (hours != 0)
      fprintf(txt, "%dh ", hours);
  
  mins = ((int) duration % 3600) / 60;
  if (mins != 0)
      fprintf(txt, "%dm ", mins);

  secs = duration - hours * 3600 - mins * 60;
  fprintf(txt, "%.2fsec\n", secs);

  return;
}
/*=============================================================================
 GetPathTail(): charge a l'adresse de 'tail' la chaine qui suit la derniere 
	occurence de '/' dans la chaine 'path'. 
        Exemple: char  *myname; ...  GetPathTail(argv[0], &myname);
=============================================================================*/

void GetPathTail(char *path, char **tail)
{
  char *p;
    
  p = strrchr(path, '/');
  if (p) {
      p++;
  }
  else {
      p = path;
  }
    
  *tail = (char *) malloc(strlen(p) + 1);
  strcpy(*tail, p);

  return;
}
/*=============================================================================
 fTest(): teste l'existence du fichier dont le nom est entre en argument.
          Provoque la sortie en cas d'echec.
=============================================================================*/

void fTest(char *fname)
{
  FILE *txt;

  if ((txt = fopen(fname, "r")) == NULL)  {
      fprintf(stderr, "%s: file not found, exit.. \n", fname);
      exit(1);
  }
  fclose(txt);

  return;
}
/*=============================================================================
 StrToD(): convertion controlee de la chaine pointee par 'str' en double, in-
           terface de 'strtod' de la lib. standard. 
           Cette fonction est construite pour controler les arguments d'un 
           programme qui en accepte un grand nombre dont certains ordonnes,
           ce qui accroit le risque d'erreurs.
           Le controle impose que la convertion soit possible (end != str) et
           qu'elle porte sur tous les caracteres de la chaine (*end == '\0'),
           en cas de succes 'status' vaut 1 en retour, sinon 0, laissant a la
           fonction appelante le soin de traiter le probleme.
=============================================================================*/

double StrToD(char *str, char *status)
{
  char    *end;
  double  x;

  *status = 1;
  x = strtod(str, &end);

  if (end == str || *end != '\0') *status = 0;

  return x;
}
/*=============================================================================
 ioError(): Signale une erreur de lecture d'un caractere non gere et provoque
            la sortie du programme.
            Arguments: le nom de la fonction appelante et le caractere concerne.
=============================================================================*/

void ioError(char *func_name, char c)
{
  fprintf(stderr, "%s: character '%c' not managed, exit..\n", func_name, c);
  exit(1);
}
/*===========================================================================*/
