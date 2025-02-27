/* miranda.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/* 										  */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details				  */
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
/*									          */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:				  */
/* miranda@cbio.mskcc.org (reaches both).					  */
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

int
main (int argc, char *argv[])
{
	
	char filename1[200];
	char filename2[200];
	char pairs_file[200];
	char fileout[200];
	char profileout[200];
	char profilein[200];
	FILE *fp1 = 0;
	FILE *fp2 = 0;
	FILE *fp_pairs = 0;
	FILE *fp_profile = 0;
	FILE *fp_read_profile =0;
	FILE *fpout = stdout;
	
	/* hit_list *hit_matrix=NULL; */
	
	mirna_profile *profile_list=NULL;
	
	/* Set Default Parameter Values */
	
	length_5p_for_weighting=8.0;  /* The 5' sequence length to be weighed  except for the last residue*/ 
	scale = 4.0;			/* The 5' miRNA scaling parameter             */
	strict = 0;			/* Strict seed model on/off                   */
	debug = 0;
	norm_scores = 0;	        /* Normalise scores according to a maximum    */
	key_value_pairs = 0;
	gap_open = -8;		/* Gap-open Penalty                           */
	gap_extend = -2;		/* Gap-extend Penalty                         */
	score_threshold = 50;		/* SW Score Threshold for reporting hits      */
	norm_score_threshold = 0;
	energy_threshold = -20;	/* Energy Threshold (DG) for reporting hits   */
	verbosity = 1;		/* Verbose mode on/off                        */
	outfile = 0;			/* Dump to file on/off                        */
	truncated = 0;		/* Truncate sequences on/off                  */
	do_shuffle = 0;		/* Generate statistics using seq shuffling    */
	no_energy = 0;		/* Turn off Vienna Energy Calcs - FASTER      */
	average = 0;			/* Some statistics for shuffled searches      */
	stdev = 0;
	z_threshold = 5.0;		/* Z-Score threshold >=           */
	p_threshold = 0.1;
	shuffle_window = 10;		/* Size of shuffling window       */
	total_shuffles = 100;		/* Total number of shuffles       */
	uniform = 0;			/* Uniform Shuffling mode on/off  */
	total_hits = 0;		/* Generic counter for alignments */
	restricted = 0;		/* Perform restricted search space */
	do_profile = 0;               /* Statistical Profiling Mode */
	do_pvalues = 0;
	evd_location = 0.0;
	evd_scale = 0.0;
	max_input_profs=1000;
	
	/* Command-line parsing begins here */
	strcpy(profileout,"profile.out");
	parse_command_line (argc, argv, &filename1, &filename2, &fileout, &profileout, &profilein, &pairs_file);
	
	printf ("%s\n",argv[3]);
	/*Doron - when running under profile mode or -evd, -norm the alignment score is normalized*/
	if (norm_scores){
        score_threshold=norm_score_threshold;
	}
	
	/* Now check our input and output files can be accessed / created */
	
	if ((fp1 = fopen (filename1, "r")) == NULL)
    {
		fprintf (stderr, "Error: Cannot open file %s\n", filename1);
		exit (1);
    }
	
	if ((fp2 = fopen (filename2, "r")) == NULL)
    {
		fprintf (stderr, "Error: Cannot open file %s\n", filename2);
		exit (1);
    }
	
	if ((do_profile) && ((fp_profile = fopen (profileout, "w")) == NULL))
	{
		fprintf (stderr, "Error: Cannot open profile file %s for output\n", profileout);
		exit (1);
	}
	
	if (((do_pvalues) && (!evd_location) && (!evd_scale)) && ((fp_read_profile = fopen (profilein, "r")) == NULL))
	{
		fprintf (stderr, "Error: Cannot open profile file %s for input\n", profilein);
		exit (1);
	}
	
	/* Load a pairs file if neccessary */
	if ((restricted) && ((fp_pairs = fopen (pairs_file, "r")) == NULL))
    {
		fprintf (stderr, "Error: Cannot open restrict pairs file %s\n", pairs_file);
		exit (1);
    }
	
	/* Output file, if needed */
	if ((outfile) && ((fpout = fopen (fileout, "w")) == NULL))
    {
		fprintf (stderr, "Error: Cannot create output file %s\n", fileout);
		exit (1);
    }
	

	print_parameters (filename1, filename2, fpout);
	
	/* Everything looks good.... Start the Scan! */
	
	
    /*
	 hit_matrix=calloc(10,sizeof(hit_list));
	 */
	
	int n_profs=0;
	
	/*Doron- Load profiles if used -evd option. */
	if ((do_pvalues) && (!evd_location) && (!evd_scale)){
        profile_list=calloc(max_input_profs,sizeof(mirna_profile)); 
        n_profs=read_profiles(fp_read_profile,profile_list+0);
        qsort (profile_list, n_profs, sizeof (mirna_profile), cmpprofiles);
    }    
	
	/*Doron - Find miRNA targets*/
	char *percent_file=argv[3];
	find_targets (fp1, fp2, fpout, fp_pairs, fp_profile, filename2,profile_list,n_profs, do_profile, percent_file);
	
    /*
	 {  hit_list* list = hit_matrix+0;
		 int i;
		 for (i=0;i<list->n_hits;i++) {
			 hit_struct* hit = list->hits+i;
			 fprintf(stdout, "score %f %s\n", hit->score,hit->alignment[2]);
		 }
	 }
	 */
	
	fclose(fp1);
	if (outfile){
        fclose(fpout);
	}
	if (fp_profile){
		fclose(fp_profile);
	}
	if (fp_read_profile){
		fclose(fp_read_profile);
	}
	if (fp_pairs){
		fclose(fp_pairs);
	}
	exit (0);
}

