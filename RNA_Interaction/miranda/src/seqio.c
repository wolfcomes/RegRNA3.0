/* seqio.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/*                                                                                */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details                              */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*                                                                                */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:                                 */
/* miranda@cbio.mskcc.org (reaches both).                                         */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either 
 * version 2 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "miranda.h"

long
readinseq (long pos, FILE * fp, char *sequence, char *description, char *name)
{

  char c;
  int i;
  int j;
  int flag = 0;

  while ((c = fgetc (fp)) != EOF)
    {

      if (c == '>')
	{
	  /* ungetc(c, fp); */
	  i = 0;
	  j = 0;
	  while ((c = fgetc (fp)) != '\n')
	    {
	      if (c == '\t')
		{
		  c = ' ';
		}

	      if (c == ' ')
		{
		  flag = 1;
		}

	      if (!flag)
		{
		  name[i] = c;
		  i++;
	      } else
		{
		  description[j] = c;
		  j++;
		}

	    }
	  name[i] = '\0';
	  description[j] = '\0';

	}
      i = 0;
      while (((c = fgetc (fp)) != EOF) && (c != '>'))
	{
	  if ((c != '\n') && (c != ' ') && (c != '\r'))
	    {

	      if (((c >= 'a') && (c <= 'z')))
		{
		  c += 'A' - 'a';
		}
	      if (((c >= 'A') && (c <= 'Z')) || (c == '*'))
		{
		  sequence[i] = c;
		  i++;
	      } else
		{
		  sequence[i] = 'X';
		  i++;
		  fprintf (stderr,
			   "Error: Sequence file contains non-standard or lowercase characters %c\n",
			   c);
		}
	    }
	}
      sequence[i] = '\0';
      ungetc (c, fp);
      return (1);
    }

  return (0);
}

int
readinpairs (FILE * fp , pair_struct **pairs)
{
  char BUFFER[400];
  char ident1[400];
  char ident2[400];
  int  n_pairs = 0;
  int  N_pairs = 16; 

  if (!(*pairs = malloc(N_pairs * sizeof (pair_struct)))) {
      fprintf(stderr, "Cannot obtain %ld pair structs\n", (long) N_pairs);
      exit;
  }

  while ( (fgets(BUFFER, 400, fp)) ){
       if (n_pairs >= N_pairs) {
           N_pairs *= 1.44;
           if (!(*pairs = realloc(*pairs, N_pairs * sizeof (pair_struct)))) {
               fprintf(stderr, "Cannot realloc to %ld pair structs\n", (long) N_pairs);
               exit;
           }
       }
       sscanf(BUFFER,"%s\t%s",(*pairs)[n_pairs].identifier1,(*pairs)[n_pairs].identifier2);
        n_pairs++;
  }

   return(n_pairs);
}

int read_profiles(FILE *fpin, mirna_profile *profile_list){

int n_profiles=0;
char mirna_id[400];
char c;
double location=0;
double scale=0;

while(fscanf(fpin,"%s\t%lf\t%lf\n",&mirna_id,&location,&scale) != EOF) {
    strcpy(profile_list[n_profiles].mirna_id,mirna_id);
    profile_list[n_profiles].location=location;
    profile_list[n_profiles].scale=scale;
    n_profiles++;
    if (n_profiles > max_input_profs){
         max_input_profs *= 1.4;
         profile_list = realloc(profile_list, max_input_profs * sizeof(mirna_profile));
    }
}

return(n_profiles);
}
