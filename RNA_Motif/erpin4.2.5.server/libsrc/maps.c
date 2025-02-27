
/*=============================================================================
 maps.c                                              A.Lambert le 30/04/03
 
 Fonctions concernant la correspondance entre un axe reel et les indices d'un
 tableau qui en representent une version discretisee.
 
 cc -O2 -Wall -c maps.c -I../include ;
 ar -rs ../lib/librnaIV.a maps.o ;

=============================================================================*/

#include "rnaIV.h"

/*
voir "rnaIV.h"

typedef struct {
  double  min, max, ratio;
  int     Pixmin, Pixmax;
  } Map;

*/

Map    MapSetup(double x1, double x2, int X1, int X2);
double Getx(Map map, int X);
int    GetX(Map map, double x);
int    GetXFloor(Map map, double x);
int    GetXCeil(Map map, double x);

/*=============================================================================
 MapSetup(): Initialise une structure 'Map'.
=============================================================================*/

Map MapSetup(double x1, double x2, int X1, int X2)
{
  Map map;
  map.min = x1;
  map.max = x2;
  map.Pixmin = X1;
  map.Pixmax = X2;
  map.ratio = (double) (X2 - X1)/(x2 - x1);
  return map;
}
/*=============================================================================
 Getx(): Retourne le nombre reel correspondant a l'entier 'X'.
=============================================================================*/

double Getx(Map map, int X)
{
  return (map.min + (double) (X - map.Pixmin) / map.ratio);
}
/*=============================================================================
 GetX(): Retourne le nombre entier correspondant au nombre reel 'x'.
=============================================================================*/

int GetX(Map map, double x)
{
  return (int) rint((double)map.Pixmin + (x - map.min) * map.ratio);
}
/*=============================================================================
 GetXFloor(): Retourne le nombre entier le plus grand inferieur a 'x', dans la
              correspondance controlee par 'map'.
=============================================================================*/

int GetXFloor(Map map, double x)
{
  return (int) floor((double) map.Pixmin + (x - map.min) * map.ratio);
}
/*=============================================================================
 GetXCeil(): Retourne le nombre entier le plus petit superieur a 'x', dans la
             correspondance controlee par 'map'.
=============================================================================*/

int GetXCeil(Map map, double x)
{
  return (int) ceil((double) map.Pixmin + (x - map.min) * map.ratio);
}
/*===========================================================================*/
