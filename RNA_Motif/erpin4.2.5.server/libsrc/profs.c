
/*=============================================================================
 profs.c                           A.Lambert le 17/11/01 revu le 12/02/04

 Ce code regroupe les fonctions procedant a la construction des profils statis-
 tiques relatifs aux structures helices et brins depuis le contenu d'un aligne-
 ment multiple constituant la "base d'apprentissage".

 Les tableaux construits sont stockes dans les champs des structures 'Helix' et
 'Strand', la memoire est restituee par appel de 'FreeStrand', 'FreeHelix'.

 les brins sont traites de la meme maniere, qu'ils contiennent ou non des gaps:

 Brin avec gaps: --> Matrice des poids respectifs des A,T,G,C et Gaps.
 Brin sans gaps: --> Matrice des poids respectifs des A,T,G,C et 0 Gaps.
 Helice:         --> Statistique des correlations entre 1/2 helices.

 Fonctions concernees:

 GetWeights       1 brin non structure avec gaps.
 GetStWeights     2 brin structure avec gaps, interface de 1.
 GetCorrels       3 helice non structuree.
 GetStProfile     4 brin structure avec ou sans gaps, interface de 2.
 GetHlxProfile    5 helice structuree.

 cc -O2 -Wall -c profs.c -I../include ;

 ar -r ../lib/librnaIV.a profs.o ;

=============================================================================*/

#include "rnaIV.h"

extern double LOG_ZERO;
extern double *NewLogDataFreqs;

double **GetWeights(char **data, int nbstr, int bgn, int len);      /* brins */
void   GetStWeights(Trset *trset, Strand *St);
                                                                  /* helices */
double **GetCorrels(char **data, int nbstr, int bgn1, int bgn2, int len);
                                                               /* interfaces */
void   GetStProfile(Trset *trset, Strand *St);
void   GetHlxProfile(Trset *trset, Helix *Hlx);
void   GetPatternProfiles(Trset *trset, Pattern *pattern);

/*=============================================================================
 GetWeights(): Mesure sur les sequences contenues dans un tableau 'data' la
               frequence des caracteres rencontres (A, T, G, C, - et N) a
               partir de 'bgn' et sur une longueur 'len'.
               Stocke les resultats dans un profil dont le nombre de lignes est
               egal a la longueur de l'alphabet utilise plus 2:
                 - une ligne pour la frequence des gaps,
                 - une ligne pour la frequence moyenne des lettres ATGC, ou
               l'on note dans chaque colonne (1 - frequence des gaps)/ALPHALEN,
               afin de gerer, lors du calcul de scores, les caracteres inconnus
               'N' comme une distribution equiprobable de A, T, G et C:
               (1/ALPHA_LEN = 0.25).
               Les resultats sont retournes dans un tableau de doubles.
               Remarque:
               Cette fonction n'a pas de structure 'Strand' ou 'Helix' comme
               argument car elle est aussi destinee a etre utilisee sur des
               brins de sequences non structures.
               Pour les brins definis par une structure 'Strand' on construit
               plus bas 'GetStProfile()', interface de celle-ci.

               Note sur la determination de Po et Pg:

               "Une facon d'exprimer qu'on ne sait rien sur une structure secon-
               daire est de supposer que toutes ses realisations sont equipro-
               bables".

               Lors de l'alignement d'une sequence aleatoire (recherche de la
               structure secondaire), pour un brin de longueur L, le nombre de
               gaps est une variable aleatoire X qui peut prendre de facon equi-
               probable les (L+1) valeurs possibles: n = 0,1,2,..,L soit une
               probabilite P(X = n) = 1/(L+1);
               Soit l'evenement Ei: "le caractere i est un gap".

               P(Ei) = Sum { P(Ei, X = n) } = Sum { P(X = n) x P(Ei | X = n) }
                        n = 0,1,..,L           n = 0,1,..,L
               P(Ei) = Sum { 1/(L+1) x n/L } = 1/(L+1) x 1/L x L.(L+1)/2 = 1/2
                        n = 0,1,..,L

               Par hypothese Ei ne depend pas de i,
               donc: Pg = 1/2.
               Par ailleurs si les |A| lettres A,T,G,C sont equiprobables:
               Po = (1 - Pg) / |A|,
               soit: Po = 1/8.      est utilise pour les 'N'

               Les elements du profil concernant ATGC sont divises par le log
               de la frequence des bases observees dans les sequences a explo-
               rer.
               Le tableau 'NewLogDataFreqs' est une variable globale declaree
               dans 'env.c'.
               Pour simplifier, les 'N' continuent d'etre interpretes comme
               une distribution uniforme de ATGC.
	       Le cas ou il n'y a pas de gaps dans le profil est gere.
=============================================================================*/

double **GetWeights(char **data, int nbstr, int bgn, int len)
{
  double **profile, Po, Pg, logPo, logPg, log1mPg, u, eps, sum;
  int    i, j, totalgaps, gaps, na, nt, ng, nc, nn, ncar, height;

  extern double *NewLogDataFreqs;
  extern double LOG_ZERO;

  eps = exp(LOG_ZERO);
  Pg = 0.5;
  Po = 1.0 / ALPHA_LEN;
  logPg   = log(Pg);
  logPo   = log(Po);
  log1mPg = log(1.0 - Pg);


  height = ALPHA_LEN + 2;                 /* hauteur du profil: nb de lignes */
  totalgaps = 0;                        /* compteur global du nombre de gaps */

  profile = dMat(height + 2, len);

  for (j = 0; j < len; j++) {     /* 2 lignes de + pour le calcul des scores */
      profile[6][j] = LOG_ZERO;
      profile[7][j] = 0.0;
  }

  for (j = 0; j < len; j++)           /* boucle sur les colonnes successives */
  {
      na = nt = ng = nc = nn = gaps = 0;

      for (i = 0; i < nbstr; i++)                   /* boucle sur les lignes */
          switch (data[i][j + bgn]) {
              case 'A': na++;   break;
              case 'T': nt++;   break;
              case 'G': ng++;   break;
              case 'C': nc++;   break;
              case '-': gaps++; break;
              case 'N': nn++;   break;
              default : ioError("GetWeights", data[i][j + bgn]);
          }
      u = nn / (double) ALPHA_LEN;
      ncar = nbstr - gaps;
      profile[_A_][j] = (na == 0 && nn == 0 ? eps : na + u);
      profile[_T_][j] = (nt == 0 && nn == 0 ? eps : nt + u);
      profile[_G_][j] = (ng == 0 && nn == 0 ? eps : ng + u);
      profile[_C_][j] = (nc == 0 && nn == 0 ? eps : nc + u);
      profile[_X_][j] = (gaps == 0 ? eps : (double) gaps);
      profile[_N_][j] = (ncar == 0 ? eps : ncar / (double) ALPHA_LEN);

      totalgaps += gaps;
  }

                         /* du decompte aux frequences et elements du profil */

  if (totalgaps == 0) log1mPg = logPg = 0.0;   /* pas de gaps dans le profil */

  for (j = 0; j < len; j++)
  {
      for (i = 0, sum = 0.0; i <= ALPHA_LEN; i++)           /* normalisation */
          sum += profile[i][j];                           /* somme sur ATGC- */

      for (i = 0; i < ALPHA_LEN; i++)
          profile[i][j] = log(profile[i][j] / sum) - log1mPg - NewLogDataFreqs[i];

      profile[_X_][j] = log(profile[_X_][j] / sum) - logPg;
      profile[_N_][j] = log(profile[_N_][j] / sum) - log1mPg - logPo;
  }

  return profile;
}
/*=============================================================================
 GetStWeights(): Interface de 'GetWeights' prenant pour argument les adresses
                 de 'trset' et d'une structure 'Strand'.
                 Le profil n'est cree que si St->Profile == NULL.
=============================================================================*/

void GetStWeights(Trset *trset, Strand *St)
{
  if (St->Profile == NULL)
      St->Profile = GetWeights(trset->data, trset->nseq,
                               St->db_bgn, St->max_len);
  return;
}
/*=============================================================================
 GetCorrels(): calcule le tableau des correlations entre 2 brins d'egale longu-
               eur 'len' debutant a 'bgn1' et 'bgn2' structures comme 2 parties
               d'helices, sur un alignement de sequences pointe par 'data'.
               Si la longueur de l'helice est 'len' le tableau est un rectangle
               de 'len' colonnes et ALPHA_LEN * ALPHA_LEN + 2 * ALPHA_LEN, soit
               24 lignes (dont les 8 dernieres destinees la gestion des 'N' lors
               du calcul des scores).
               - On gere les caracteres inconnus 'N' en leur associant une dis-
               tribution equiprobable de ATGC (1/ALPHA_LEN).
=============================================================================*/

double **GetCorrels(char **data, int nbstr, int bgn1, int bgn2, int len)
{
  int     i, j, k, n, i1 = 0, i2 = 0, findN, height, end2;
  char    *h1, *h2;
  double  **hprofile, sum, eps, logPo, u, alphalen_inv, sqr_alphalen_inv;

  extern double *NewLogDataFreqs;
  extern double LOG_ZERO;

  alphalen_inv = 1.0 / ALPHA_LEN;
  sqr_alphalen_inv = 1.0 / SQR_ALPHA_LEN;

  eps = exp(LOG_ZERO);
  logPo = log(sqr_alphalen_inv);            /* si distr. uniforme: Po = 1/16 */
  end2 = bgn2 + len - 1;
  height = SQR_ALPHA_LEN + 2 * ALPHA_LEN;         /* nb de lignes du tableau */

  hprofile = dMat(height + 2, len);
  FilldMat(hprofile, height + 2, len, 0.0);
  for (j = 0; j < len; j++)  hprofile[24][j] = LOG_ZERO;    /* 2 lignes de + */
                                                /* pour le calcul des scores */

  for (n = 0; n < nbstr; n++)              /* boucle sur la base des donnees */
  {
      h1 = data[n] + bgn1;           /* adresse du 1er caractere du 1er brin */
      h2 = data[n] + end2;      /* adresse du dernier caractere du 2eme brin */

      for (j = 0; j < len; j++)
      {
          findN = 0;
          switch(h1[j]) {
              case 'A': i1 = ALPHA_LENxA; break;
              case 'T': i1 = ALPHA_LENxT; break;
              case 'G': i1 = ALPHA_LENxG; break;
              case 'C': i1 = ALPHA_LENxC; break;
              case 'N': findN += 1;       break;
              default : ioError("GetCorrels", h1[j]);
          }
          switch(h2[-j]) {
              case 'A': i2 = _A_;   break;
              case 'T': i2 = _T_;   break;
              case 'G': i2 = _G_;   break;
              case 'C': i2 = _C_;   break;
              case 'N': findN += 2; break;
              default : ioError("GetCorrels", h2[-j]);
          }
          switch (findN) {
              case 0:                       /* h1[j] != 'N' && h2[-j] != 'N' */
                  hprofile[i1 + i2][j]++;
                  break;
              case 1:                       /* h1[j] == 'N' && h2[-j] != 'N' */
                  for (i = 0; i < ALPHA_LEN; i++)
                      hprofile[i * ALPHA_LEN + i2][j] += alphalen_inv;
                  break;
              case 2:                       /* h1[j] != 'N' && h2[-j] == 'N' */
                  for (i = 0; i < ALPHA_LEN; i++)
                      hprofile[i1 + i][j] += alphalen_inv;
                  break;
              case 3:                       /* h1[j] == 'N' && h2[-j] == 'N' */
                  for (i = 0; i < SQR_ALPHA_LEN; i++)
                      hprofile[i][j] += sqr_alphalen_inv;
                  break;
          }
      }
  }
          /* ------ traite les elements vides des 16 premieres lignes ------ */

  for (i = 0; i < SQR_ALPHA_LEN; i++)
      for (j = 0; j < len; j++)
          if (hprofile[i][j] < eps) hprofile[i][j] = eps;

  for (j = 0; j < len; j++)                    /* normalisation des colonnes */
  {
      for (i = 0, sum = 0.0; i < SQR_ALPHA_LEN; i++)
          sum += hprofile[i][j];
      for (i = 0; i < SQR_ALPHA_LEN; i++)
          hprofile[i][j] /= sum;
  }

      /* ----- lignes 16 a 19 du profil gerant les couples (N, A|T|G|C) ---- */

  for (j = 0; j < len; j++)                /* boucle sur toutes les colonnes */
  {
      i1 = SQR_ALPHA_LEN;
      for ( i = 0; i < ALPHA_LEN; i++)       /* boucle sur les lignes 16..19 */
      {
          for (k = 0; k < ALPHA_LEN; k++)
              hprofile[i1 + i][j] += hprofile[ALPHA_LEN * k + i][j];

          hprofile[i1 + i][j] *= alphalen_inv;                    /* moyenne */
      }
  }

      /* ----- lignes 20 a 23 du profil gerant les couples (A|T|G|C, N) ---- */

  for (j = 0; j < len; j++)                /* boucle sur toutes les colonnes */
  {
      i1 = SQR_ALPHA_LEN + ALPHA_LEN;
      for ( i = 0; i < ALPHA_LEN; i++)       /* boucle sur les lignes 20..23 */
      {
          for (k = 0; k < ALPHA_LEN; k++)
              hprofile[i1 + i][j] += hprofile[ALPHA_LEN * i + k][j];

          hprofile[i1 + i][j] *= alphalen_inv;                    /* moyenne */
      }
  }

                          /* ------ passage a l'echelle logarithmique ------ */

  for (i = 0; i < SQR_ALPHA_LEN; i++)             /* les 16 premieres lignes */
  {
      u = NewLogDataFreqs[i / ALPHA_LEN] + NewLogDataFreqs[i % ALPHA_LEN];
      for (j = 0; j < len; j++)
      {
          hprofile[i][j] = log(hprofile[i][j]) - u;
      }
  }
  for (i = SQR_ALPHA_LEN; i < height; i++)         /* les 8 lignes suivantes */
      for (j = 0; j < len; j++)
      {
          hprofile[i][j] = log(hprofile[i][j]) - logPo;
      }

  return hprofile;
}
/*=============================================================================
 GetStProfile(): Calcule le profil du brin pointe par 'St'.
=============================================================================*/

void GetStProfile(Trset *trset, Strand *St)
{
      GetStWeights(trset, St);

  if (St->max_gaps != 0)                      /* necessaire seulement si gaps */
  {
      St->str = (char *) malloc(St->max_len + 1);
      St->Align = dMat(St->max_len + 1, St->max_len + 1);
      FilldMat(St->Align, St->max_len + 1, St->max_len + 1, 0.0);
  }
  return;
}
/*=============================================================================
 GetHlxProfile(): Interface de 'GetCorrels' prenant pour argument les adresses
                  de 'trset' et d'une structure 'Helix'.
                  Le profil n'est cree que si Hlx->Profile == NULL.
=============================================================================*/

void GetHlxProfile(Trset *trset, Helix *Hlx)
{
  if (Hlx->Profile == NULL)
      Hlx->Profile = GetCorrels(trset->data, trset->nseq,
                            Hlx->db_bgn1, Hlx->db_bgn2, Hlx->helix_len);

  return;
}
/*=============================================================================
 GetPatternProfiles(): Cree l'ensemble des profils d'un pattern.
=============================================================================*/

void GetPatternProfiles(Trset *trset, Pattern *pattern)
{
  int i;

  for (i = 0; i < pattern->nst; i++)
  {
      GetStProfile(trset, pattern->std + i);
  }
  for (i = 0; i < pattern->nhx; i++)
  {
      GetHlxProfile(trset, pattern->hlx + i);
  }
  return;
}
/*===========================================================================*/
