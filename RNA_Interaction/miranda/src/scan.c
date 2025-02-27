/* scan.c */
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
#include <math.h>
#include <ctype.h>
#include "utils.h"
#include "miranda.h"


/* FIXME:
 * this is still not sufficient if there are many small sequence stretches inbetween
 * N stretches.
*/

int sequence_length_effective(const char* sequence)
{  const char* x = sequence;
	int n_base = 0;
	while (*x)
		if (*(x++) != 'N')
			n_base++;
	return n_base;
}

/* Load Sequences and Set-up the Alignment Run */

int
find_targets (FILE * fp1, FILE * fp2, FILE * fpout, FILE * fp_pairs, FILE * fp_profile, char *filename, mirna_profile *profile_list, int n_profs, int do_profile, char *percent_file)
{
	
	/* The three key alignment matrices */
	int **matrix1;		/* Best score of all three states (nt-nt,nt-gap,gap-nt  */
	int ***matrix2;		/* Traceback Matrix           */
	int **matrix3;		/* best score for state nt-nt */
	int **matrix4;		/* best score for state gap-nt */
	int **matrix5;		/* best score for state nt-gap */
	int **nt_nt_score;            /* nt nt match matrix for easy lookup */
	
	int utr_processed = 0;
	int mirna_processed = 0;
	int i = 0;
	int max_possible_score=0;
	int seqlen1;
	int seqlen2;
	int total_pairs;
	hit_struct hit;		/*Struct to store hit information*/
	hit_struct *hit_ptr;
	score_struct *scores;	/*Score struct for non-optimal path detection*/
	final_score finalscore;
	final_score *fscore_ptr;
	pair_struct *pairs;	/* Score structure for non-optimal path detection (?)*/
	double profile_location=0;
	double profile_scale=0;
	
	/* Sequence Information, IDs and Descriptions */
	char *query;			/* The Query Sequence (miRNA)             */
	char *query_des;
	char *query_id;
	char *rev_query;
	char *reference;		/* The Reference Sequence (UTR)           */
	char *reference2;		/* A Second Copy of the Ref for Shuffling */
	char *reference_des;
	char *reference_id;
	
	/* Scoring Information */
	
	double end_score = 0.0;
	double max_score = 0.0;
	double maximum = 0;
	double *dist;
	
	/* File IO            */
	long stream_pos_q = 0;
	long current_pos_q = 0;
	long stream_pos_r = 0;
	long current_pos_r = 0;
	
	species_profile sp_prof;
	
	if (do_profile) {
		sp_prof.n_profiles_max = 1000;
		sp_prof.n_profiles = 0;
		sp_prof.profiles = calloc(sp_prof.n_profiles_max,sizeof(mirna_profile));  
	}
	
	hit_ptr = &hit;
	fscore_ptr = &finalscore;
	
	/* Memory Allocation for Sequences */
	query = (char *) calloc (SEQUENCE_LENGTH, sizeof (char));
	query_des = (char *) calloc (IDENTIFIER_LENGTH, sizeof (char));
	query_id = (char *) calloc (IDENTIFIER_LENGTH, sizeof (char));
	
	reference = (char *) calloc (SEQUENCE_LENGTH, sizeof (char));
	reference_des = (char *) calloc (IDENTIFIER_LENGTH, sizeof (char));
	reference2 = (char *) calloc (SEQUENCE_LENGTH, sizeof (char));
	reference_id = (char *) calloc (IDENTIFIER_LENGTH, sizeof (char));
	
	rev_query = (char *) calloc (SEQUENCE_LENGTH, sizeof (char));
	
	/* Array to store distribution of shuffled alignments */
	dist = (double *) calloc (total_shuffles, sizeof (double));
	
	/* Initialize the pairs list for restrived searches */
	if (restricted) {
        pairs = NULL;
        total_pairs=load_pairs(fp_pairs,&pairs);
        if (verbosity){
			printf("Performing Restricted Scan on:%d pairs\n",total_pairs);
        }
	}
	
	/* Prepare the generic base lookup array */
	initialize_bases ();
	
	/* Read the query sequence(s) (microRNA(s)) from a FASTA file */
	while ((current_pos_q =
			readinseq (stream_pos_q, fp1, query, query_des, query_id)))
    {
		mirna_profile * mpr = NULL;
		
		if (do_profile) {
			if (sp_prof.n_profiles >= sp_prof.n_profiles_max) {
                sp_prof.n_profiles_max *= 1.4;
                sp_prof.profiles = realloc(sp_prof.profiles, sp_prof.n_profiles_max * sizeof(mirna_profile));
			}
			mpr = sp_prof.profiles+sp_prof.n_profiles;
			mpr->scores = calloc(40000,sizeof(double));
			mpr->cdf =    calloc(40000,sizeof(double));
			mpr->location = 0.0;
			mpr->scale    = 0.0;
			mpr->n_scores = 0;
			mpr->n_scores_max = 40000;
			strncpy(mpr->mirna_id, query_id, 99);
		}
		
		if ( (verbosity) || (debug))
		{
			fprintf (fpout, "Read Sequence:%s %s(%d nt)\n", query_id, query_des,
					 (int) strlen (query));
			//printf ("%s\n",percent_file);
  			FILE * pFile;
			pFile = fopen (percent_file,"w");
			fputs (query_des,pFile);
			fclose (pFile);



		}
		
		
		/* We are doing alignments like this:
			* 
			*           microRNA
			*   3'-<<<<<<<<<<<<<<<<<<<-5'
			*        |||o|||||  ||||| 
			*   5'->>>>>>>>>>>>>>>>>>>-3'
			*      Reference Sequence 
			*
			*
			* Hence we should reverse one of the two sequences 
			*/
		
		/* Reverse the wquery (microRNA) sequence */
		strcpy (rev_query, query);
		revstring (query);
		
		/* Loop over all reference sequences in FASTA file */
		/* Do full scan for each                           */
		
		fclose (fp2);
		if ((fp2 = fopen (filename, "r")) == NULL)
		{
			fprintf (stderr, "Error: Cannot open file %s\n", filename);
			exit (1);
		}
		
		stream_pos_r = 0;
		
		while
			(  (  current_pos_r
				  =  readinseq (stream_pos_r, fp2, reference, reference_des, reference_id)
				  )
			   )
		{ int reference_len = sequence_length_effective(reference);
				
				if
					(  (!restricted || find_pair(query_id,reference_id,total_pairs,pairs))
					   && (!do_profile || (reference_len >= 300))
					   )
				{
						/* Keep track of the number of sequences scanned so far */
						utr_processed++;
						
						if ((verbosity) || (debug))
						{
							fprintf (fpout, "Read Sequence:%s %s(%d nt)\n", reference_id,
									 reference_des, (int) (reference_len));
						}
						
						if (truncated)
						{
							reference[truncated] = '\0';
						}
						
						/* Get sequence lengths for query and reference */
						/* we recompute strlen(reference); it might have been truncated */
						
						seqlen1 = strlen (query);
						seqlen1_global=seqlen1;
						seqlen2 = strlen (reference);
						
						/*Doron - the length of the miRNA 3'-end not weighted*/
						length_3p_for_weighting=seqlen1-length_5p_for_weighting;
						overlap_cutoff=(seqlen1);
						
						strcpy (reference2, reference);
						
						/*
						 if (norm_scores){
							 max_possible_score=get_max_score(query,length_3p_for_weighting);
						 }
						 */ 
						
						
						/* Initialize the hit / alignment constructs for this sequence */
						
						
						hit.alignment[0] =
							calloc (seqlen1 + seqlen2, sizeof (char));
						hit.alignment[1] =
							calloc (seqlen1 + seqlen2, sizeof (char));
						hit.alignment[2] =
							calloc (seqlen1 + seqlen2, sizeof (char));
						hit.rest[0] =  calloc (30, sizeof (char));
						hit.rest[1] =  calloc (30, sizeof (char));
						hit.rest[2] =  calloc (30, sizeof (char));
						hit.rest[3] =  calloc (30, sizeof (char));
						hit.rest[4] =  calloc (30, sizeof (char));
						hit.rest[5] =  calloc (30, sizeof (char));
						
						
						
						
						/* Structure for sub-optimal score list */
						scores =
							(score_struct *) calloc (seqlen1 * seqlen2,
													 sizeof (score_struct));
						
						/* Initialize the three alignment matrices */
						matrix1 = calloc ((seqlen1 + 1), sizeof (int *));
						matrix2 = (int ***) calloc (4, sizeof (int **));
						matrix3 = calloc ((seqlen1 + 1), sizeof (int *));
						matrix4 = calloc ((seqlen1 + 1), sizeof (int *));
						matrix5 = calloc ((seqlen1 + 1), sizeof (int *));
						nt_nt_score=calloc ((seqlen1 + 1), sizeof (int *));
						
						/*Initialize 4-D call-back matrix*/
						for (i=0; i <4;i++){
							matrix2[i]= (int **) calloc ((seqlen1 + 1), sizeof (int *));
						}
						
						
						for (i = 0; i < seqlen1 + 1; i++)
						{
							matrix1[i] = calloc ((seqlen2 + 1), sizeof (int));
							
							matrix2[0][i] = (int *) calloc ((seqlen2 + 1), sizeof (int));
							matrix2[1][i] = (int *) calloc ((seqlen2 + 1), sizeof (int));
							matrix2[2][i] = (int *) calloc ((seqlen2 + 1), sizeof (int));
							matrix2[3][i] = (int *) calloc ((seqlen2 + 1), sizeof (int));
							
							matrix3[i] = calloc ((seqlen2 + 1), sizeof (int));
							matrix4[i] = calloc ((seqlen2 + 1), sizeof (int));
							matrix5[i] = calloc ((seqlen2 + 1), sizeof (int));
							nt_nt_score[i] = calloc ((seqlen2 + 1), sizeof (int));
							matrix1[i][0] =  matrix3[i][0] = nt_nt_score[i][0]=0;
							
							matrix2[0][i][0] =matrix2[1][i][0]=matrix2[2][i][0]=matrix2[3][i][0]=0;
							matrix4[i][0]=0;matrix5[i][0]=0;
							
							/*	       matrix4[i][0]=-10000;matrix5[i][0]=gap_open+(i*gap_extend); */
						}
						
						for (i = 0; i < seqlen2 + 1; i++)
						{
							matrix1[0][i] = matrix3[0][i] = nt_nt_score[0][i]=0;
							matrix4[0][i] =0;matrix5[0][i] =0;
							/*	      matrix4[0][i] =gap_open+(i*gap_extend);matrix5[0][i] =-10000; */
							matrix2[0][0][i] =matrix2[1][0][i]=matrix2[2][0][i]=matrix2[3][0][i]=0;
						}
						
						if (verbosity && do_shuffle)
						{
							fprintf
							(fpout,
							 "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
							fprintf (fpout,
									 "Generating Alignment Distribution of Shuffled Sequences\n");
							fprintf (fpout,
									 "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
						}
						
						if (uniform)
						{
							shuffle_window = seqlen2;
						}
						
						
						if (do_shuffle)
						{
							irand (0);
							for (i = 0; i < total_shuffles; i++)
							{
								end_score = 0;
								shuffle (reference2, seqlen2, shuffle_window);
								/*	  dist[i] =
									build_matrix_quick (matrix1, matrix2, query, reference2,
														seqlen1, seqlen2);
								
								this functionality is disabled as a quick fix  */
							}
							
							for (i = 0; i < total_shuffles; i++)
							{
								average += (dist[i]);
								if (dist[i] > maximum)
								{
									maximum = dist[i];
								}
							}
							average = average / (double) total_shuffles;
							
							for (i = 0; i < total_shuffles; i++)
							{
								stdev += ((dist[i] - average) * (dist[i] - average));
							}
							stdev = stdev / (double) (total_shuffles - 1);
							stdev = sqrt (stdev);
						}
						
						
						
						if (do_shuffle)
						{
							fprintf (fpout, "done\t");
							fprintf (fpout,
									 "Average: %3.2f\tSt. Dev: %3.2f\tMax: %3.2f\n",
									 average, stdev, maximum);
						}
						
						if(verbosity || debug){
							fprintf
							(fpout,
							 "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
							fprintf (fpout, "Performing Scan: %s vs %s\n", query_id,
									 reference_id);
							fprintf (fpout,
									 "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
							fflush(fpout);
						}
						
						
						
						if (do_pvalues) {
							
							if ((!evd_location) && (!evd_scale)){
								mirna_profile *result;
								mirna_profile search;
								strcpy(search.mirna_id,query_id);
								result=(mirna_profile *) bsearch(&search,profile_list,n_profs,sizeof(mirna_profile), cmpprofiles);
								if (result == NULL){
									fprintf(stderr,"Error: Profile NOT FOUND for %s\n",query_id);
									exit(1);
								} else {
									if (debug){
										fprintf(fpout,"FOUND %s %f %f\n",result->mirna_id,result->location,result->scale);
									}
									profile_location=result->location;
									profile_scale=result->scale;
								}
							} else {
								profile_location=evd_location;
								profile_scale=evd_location;
							}
						} /*End of do_pvalues*/
          
          if (do_profile) {
			  double score =
			  do_mirna_utr_profile (matrix1, matrix2, matrix3, matrix4, matrix5, nt_nt_score,query, reference, scores,
									max_possible_score, seqlen1, seqlen2, 1,
									fpout); 
			  if (mpr->n_scores >= mpr->n_scores_max) {
				  mpr->n_scores_max *= 1.4;
				  mpr->scores = realloc(mpr->scores, mpr->n_scores_max * sizeof (double));
				  mpr->cdf =    realloc(mpr->cdf, mpr->n_scores_max * sizeof (double));
			  }
			  mpr->scores[mpr->n_scores] = score+0.5;         /* extreme int/double/int conversion protection */
			  mpr->n_scores++;
			  /*Doron - run alignment*/
		  } else {
			  end_score =
			  do_alignment (matrix1, matrix2, matrix3, matrix4, matrix5, nt_nt_score,query, reference, scores,
							max_possible_score, seqlen1, seqlen2, 1, fscore_ptr, FORWARD,
							query_id, reference_id, &hit, fpout, do_profile, profile_location, profile_scale);
          }

	  if ( (verbosity) || (debug))
	  {
	      fprintf (fpout, "Score for this Scan:\n");
	  }

	  if (end_score > 0.0)
	  {
	      fprintf
		  (fpout,
		   "Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions\n");
		  
	      if (!no_energy)
		  {
			  fprintf
			  (fpout,
			   ">>%s\t%s\t%2.2f\t-%2.2f\t%2.2f\t%2.2f\t%d\t%d\t%d\t%s\n",
			   query_id, reference_id, finalscore.total_score,
			   end_score, finalscore.max_score, finalscore.max_hit,
			   utr_processed, seqlen1, seqlen2, finalscore.positional);
		  } else
		  {
			  fprintf
			  (fpout,
			   ">>%s\t%s\t%2.2f\t0.0\t%2.2f\t0.0\t%d\t%d\t%d\t%s\n",
			   query_id, reference_id, finalscore.total_score, 
			   finalscore.max_score, utr_processed, seqlen1, seqlen2,
			   finalscore.positional);
			  
		  }
	      fflush (fpout);
	  } else
	  {
		  if ((verbosity) || (debug) ){
			  fprintf (fpout, "No Hits Found above Threshold\n");
		  }
	  }
	  
	 if (verbosity || debug || end_score > 0.0) {
		 
		 fprintf (fpout, "Complete\n\n");
	 }
		/*Doron - clean up and free memory*/
	  fflush (fpout);
	  stream_pos_r = current_pos_r;

	  for (i = 0; i < seqlen1 + 1; i++)
	  {
	      free (matrix1[i]);
	      free (matrix2[0][i]);
	      free (matrix2[1][i]);
	      free (matrix2[2][i]);
	      free (matrix2[3][i]);
	      free (matrix3[i]);
	      free (matrix4[i]);
	      free (matrix5[i]);
	      free (nt_nt_score[i]);
	  }
            free (matrix2[0]);
            free (matrix2[1]);
            free (matrix2[2]);
            free (matrix2[3]);

	  free (matrix1);
	  free (matrix2);
	  free (matrix3);
	  free (matrix4);
	  free (matrix5);
	  free (nt_nt_score);


    
	  free (hit.alignment[0]);
	  free (hit.alignment[1]);
	  free (hit.alignment[2]);
	  free (hit.rest[0]);
	  free (hit.rest[1]);
	  free (hit.rest[2]);
	  free (hit.rest[3]);
	  free (hit.rest[4]);
	  free (hit.rest[5]);
      
		

	  free (scores);

				} else {
					/* printf("Skipped: %s vs %s (length %d)\n",query_id,reference_id, reference_len); */
				} 

		}
      stream_pos_q = current_pos_q;
      stream_pos_r = 0;
      current_pos_r = 0;
      mirna_processed++;

      if (do_profile) {
		  sp_prof.n_profiles++;
      }
    }
  
  fprintf (fpout, "Scan Complete\n\n");


  if (do_profile){
	  fprintf (fpout, "Generating Alignment Statistics\n");
	  fprintf (fpout, "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
	  
	  for (i=0;i<sp_prof.n_profiles;i++){
		  mirna_profile* mpr=sp_prof.profiles+i;
		  int j=0;
		  double slope=0.0;
		  double intercept=0.0;
		  double correlation=0.0;
		  /*
		   fprintf (fpout, "Calibrating: %s\n",mpr->mirna_id);
		   fprintf (fpout, "\tSorting Score Distribution\n");
		   */
		  
		  qsort(mpr->scores, mpr->n_scores, sizeof(double), cmp_dbl_asc); 
		  
		  /*
		   fprintf (fpout, "\tGenerating Empirical Cumulative Density Function\n");
		   */
		  
		  for (j=0;j<mpr->n_scores;j++){
			  mpr->cdf[j] =   j + 1 < mpr->n_scores
			  ? log (-log((j+1) / (double) mpr->n_scores))
			  : mpr->cdf[j-1];
		  }
		  
		  
		  /*fprintf (fpout, "\tPerforming Linear Regression\n");*/
		  correlation=linear_regression(mpr->scores,mpr->cdf,mpr->n_scores,&slope,&intercept);
		  
		  fprintf(fpout,"\tSlope: %f\tIntercept: %f\tCorrelation [R^2]: %f\n\n",slope,intercept,correlation*correlation);
		  mpr->scale = -(1/slope);
		  mpr->location    = mpr->scale * intercept;
		  fprintf(fpout,"%s\t%f\t%f\n",mpr->mirna_id,mpr->location,mpr->scale);
		  fprintf(fp_profile,"%s\t%f\t%f\n",mpr->mirna_id,mpr->location,mpr->scale);
		  
	      if (1)
			  graph_cdf(mpr->scores,mpr->cdf,mpr->n_scores,mpr->location,mpr->scale,fpout);  
		  
		  
		  free(mpr->scores);
		  free(mpr->cdf); 
		  
	  }
  } /*End of do_profile*/




    fflush (fpout);
    if (do_profile){
        free(sp_prof.profiles); 
    }

    free(query);
    free(query_des);
    free(query_id);
    free(reference);
    free(reference_des);
    free(reference_id);
    free(reference2);
    free(rev_query);
    free(dist);

  return (1);
}


double normalise_score (int score, int mirna_len, int utr_len) {
     return (score/( 0 + log(mirna_len * (utr_len < 300 ? 300 : utr_len))));
}

/*Doron - Main alignment function*/
double
do_alignment (int **best, int ***track, int **a_nt_nt, int **b_gap_nt, int **c_nt_gap, int **nt_nt_score, 
			  char *query, char *reference, score_struct * scores, int max_possible_score, 
			  int seqlen1,int seqlen2, int verbose, 
			  final_score * finalscore,int direction, char *query_id, char *reference_id,
			  hit_struct *hit, FILE * fpout, int do_profile, double location, double scale)
{
	
	int i = 0;
	int j = 0;
	int z = 0;
	int utr_offset3p;
	int utr_offset5p;
	double energy = 0;
	double scan_score = 0;
	double z_score = 0;
	double identity = 0;
	int hit_cluster[100];
	int valid_hits = 0;
	int good_call = 0;
	int cmin = 0;
	int cmax = 0;
	int diff = 0;
	char strop1[200];
	char strop2[200];
	int non_gap_count = 0; /*Counters for gaps and matches in the alignment*/
	int perfect_match_count = 0;
	int gap_count = 0;
	
	int fail = 0;
	int mypos = 0;
	int seed_gaps = 0;
	double hit_score;
	int tmp_integer;
	double p_value=0;
	
	
	int *good_ones_starts_j, *good_ones_ends_j,good_ones_count;
	
	good_ones_count=-1;
	good_ones_starts_j=(int *) calloc (seqlen2, sizeof (int));
	good_ones_ends_j=(int *) calloc (seqlen2, sizeof (int));
	
	total_hits = 0;
	finalscore->no_hits = 0;
	finalscore->max_hit = 0;
	finalscore->max_score = 0;
	finalscore->scan_score = 0;
	finalscore->total_score = 0;
	finalscore->positional[0] = '\0';
	
	get_nt_nt_seq_scores (nt_nt_score,query,reference,seqlen1,seqlen2);
	
	/*Doron - build_matrix core alignment function */
	build_matrix (best, track, a_nt_nt, b_gap_nt,c_nt_gap,nt_nt_score, query, 
				  reference, seqlen1, seqlen2, scores);
	
	qsort (scores, total_hits, sizeof (score_struct), cmpscores);
	
	for (i = 0; i <=total_hits; i++)
	{
		utr_offset3p=0;
		utr_offset5p=0;
		good_call = 1;
		clear_hit(hit,seqlen1,seqlen2);
		hit_score=scores[i].score;
		
		
		if (norm_scores){
			/* hit_score=(int) ((double) ((double) hit_score / (double) max_possible_score) * 100); */
			int seqlen2_eff = sequence_length_effective(reference);
			hit_score = (double) normalise_score((double) hit_score,seqlen1,seqlen2_eff);
		}
		
		/*Doron - does miRanda reports the optimal alignment only or all alignments>threshold ? */
		if (hit_score >= score_threshold)
		{
			
			if (debug){
				fprintf (fpout,"[%d] score i j path ->> %d %d %d %d\n",i,scores[i].score,scores[i].i,scores[i].j,scores[i].path); 
				fflush(fpout);
			}
			
			/*Doron - traceback (?)*/
			traceback (best, track, query, reference, scores[i].i, scores[i].j, 1,
					   hit, hit_score);
			
			/*Doron - testfor_overlap (?) , good_ones *(?)*/
			good_call=testfor_overlap (good_ones_starts_j,good_ones_ends_j,&good_ones_count,
									   hit->ref_start,hit->ref_end);
			if (good_call==1){
				good_ones_starts_j[good_ones_count]=hit->ref_start;
				good_ones_ends_j[good_ones_count]=hit->ref_end;
			}
			 
							                   
			/*Doron - miRNA alignment, convert un-aligned nt to lowercase
			* hit->rest[0,1,2] = The 5' unaligned regions in the query alignments and ref
			* hit->rest[3,4,5] = The 3' unaligned regions in the query, alignment and ref*/
			if (hit->query_start >= 1)
			{
				for (j = 0; j <= hit->query_start - 1; j++)
				{
					diff = hit->query_start - j;
					hit->rest[0][j] = tolower(query[j]);
					tmp_integer=hit->ref_start - diff;
					if (tmp_integer >=0){
						hit->rest[1][j] = tolower(reference[tmp_integer]);
						utr_offset3p++;
					}else{
						
						hit->rest[1][j] = '-';
					}
					
					hit->rest[2][j] = ' ';
				}
			} /* if hit->query_start >= 1) */

			if ((hit->query_end) < seqlen1)
			{
				for (j = hit->query_end; j < seqlen1; j++)
				{
					diff = j - hit->query_end;
					hit->rest[3][j - hit->query_end] = tolower(query[j]);
					tmp_integer=hit->ref_end + diff;
					if (tmp_integer < seqlen2){
						hit->rest[4][j - hit->query_end] =
						tolower(reference[tmp_integer]);
						utr_offset5p++;
					} else{
						hit->rest[4][j - hit->query_end] = '-';
					} /* if (tmp_integer < */
					hit->rest[5][j - hit->query_end] = ' ';
				}
			}/* if hit->query_end <seqlen1) */

			/* Adjusting for offset due to local alignment in next two lines */

			hit->ref_end+=(utr_offset5p-1);
			hit->ref_start-=utr_offset3p;

			/*strop1[0] = '\0';
			strop2[0] = '\0'; */
			sprintf (strop1, "%s%s%s", hit->rest[3], hit->alignment[0],
					 hit->rest[0]);
			sprintf (strop2, "%s%s%s", hit->rest[5], hit->alignment[1],
					 hit->rest[2]);

			/*
			 printf ("   QueTE:    3' %s5'\n                %s\n   RTT:      5' %s 3'\n\n",
					 hit->alignment[0], 
					 hit->alignment[1], 
					 hit->alignment[2]);
			 */
			/*
			 printf("strop1 >%s:\n",strop1);
			 printf("strop2 >%s:\n",strop2);
			 */

			mypos = 0;
			fail = 0;
			non_gap_count = perfect_match_count = gap_count = 0; /*Doron - non-gap positions, matches,gaps*/

			/* This code now looks for strict seed matches */
			if (strict)
			{
				/*Doron - traverse the alignment*/
				for (j = 0; j < strlen (strop1); j++)
				{
					if (strop1[j] != '-') /*Doron - if no gaps in the miRNA alignment*/
					{
						mypos++;
					}
					
					if ((mypos >= 2) && (mypos <= 8)) /*Doron - extended changed 7->8*/
					{
						if (strop2[j] != ' ') /*Doron - if no gaps in the alignment i.e. if `|` or `:`*/
						{
							non_gap_count++;
						}
						
						if (strop2[j] == '|')
						{
							perfect_match_count++;
						}
						
						if (strop1[j] == '-'){
							gap_count++;
						}
					}/* if  ((mypos >= 1) && (mypos <= 3)) */

				/*printf("On Pos: %d [%d]-->%c %c\t[%d:%d]\n",mypos,j,strop1[j],strop2[j],count1,count2);*/
				}
				
				/*Doron -fail if less than 7 `|:` char in the alignment
				* if less than 7 perfect matches and any gaps*/
				
				if ( non_gap_count <7 || perfect_match_count< 7 || gap_count>0)
				{
					fail=1;
					good_call=0;
				}  /* if !strict) 1 */

			} /* if !strict) 2 */

			/*printf("Model Analysis: %d %d %d %d\n",count1,count2,count3,count4); */
			/*Doron - compute percent identity*/
			identity=0;
			for (j = 0; j < strlen (hit->alignment[0]); j++)
			{
				if (hit->alignment[1][j] == '|')
				{
					identity++;
				}
			}/* for (j = 0; j < strlen (hit->alignment[0]); j++) */

			/*Doron - this %identity is not used in anyway. The printed %identity and similarity is 
			recomputed in printhit function*/
			identity = (identity / strlen (hit->alignment[0])) * 100;

			if (do_pvalues){
				p_value=1-((exp(-exp(-(((hit->score)-location)/scale)))));
				if (debug){
					printf("Score: %f\tP-Value: %f [%e] [location=%f scale=%f]\n",hit->score,p_value,p_value,location,scale);
				}
			}

			if (do_shuffle){
				z_score = (hit->score - average) / stdev;
			} else {
				z_score = 1000000;
			}

			if (!no_energy)
			{
				energy = get_energy (hit);
			} else
			{
				energy = -1000000;
			}

			if ((energy < energy_threshold) && (z_score >= z_threshold) && (p_value <= p_threshold) )
			{			
				if (valid_hits == 0)
				{
					hit_cluster[valid_hits] = hit->ref_start;
					valid_hits++;
					cmax = hit->ref_start;
					cmin = hit->ref_start;
				
				} else {			
					for (z = 0; z < valid_hits; z++)
					{
						if (hit_cluster[z] > cmax)
						{
							cmax = hit_cluster[z];
						}
						
						if (hit_cluster[z] < cmin)
						{
							cmin = hit_cluster[z];
						}
						
					} /* for (z = 0; z < valid_hits; z++) */


				} /* if (valid_hits == 0) */

				/*Doron - good_call, a good alignment that passes score, energy, etc. requirements*/
				if ( (good_call) )
				{
					
					hit_cluster[valid_hits] = hit->ref_start;
					valid_hits++;
					
					scan_score += (energy * -1);
					finalscore->no_hits++;
					sprintf (finalscore->positional, "%s %d",
							 finalscore->positional, hit->ref_start);
					
					if (energy < finalscore->max_hit)
					{
						finalscore->max_hit = energy;
					}
					
					
					finalscore->total_score += hit->score;
					if (hit->score > finalscore->max_score)
					{
						finalscore->max_score = hit->score;
					}
					
					if (!key_value_pairs){
						printhit (query_id, reference_id, hit, query, reference,
								  direction, z_score, p_value, energy, fpout);
					} else {
						printhit_tskv (query_id, reference_id, hit, query, reference,
									   direction, z_score, p_value, energy, fpout);
					}
					
				} /* if (good_call) */
			} /* if ((energy < energy_threshold) && (z_score >= z_threshold)) */
		} /* if (scores[i].score > score_threshold) */
		scores[i].score = scores[i].path = scores[i].i = scores[i].j = 0;
	} /* for i =0 to total hits */

	free (good_ones_starts_j);
	free (good_ones_ends_j);

	return (scan_score);
}

double
do_mirna_utr_profile (int **best, int ***track, int **a_nt_nt, int **b_gap_nt, int **c_nt_gap, int **nt_nt_score, 
					  char *query, char *reference,score_struct * scores, int max_possible_score, 
					  int seqlen1,int seqlen2, int verbose, 
					  FILE * fpout )
{
	
	double score=0;
	
	/*
	 hit_struct hit1;
	 hit_struct *hit=&hit1;
	 
	 hit1.alignment[0] =
	 calloc (seqlen1 + seqlen2, sizeof (char));
	 hit1.alignment[1] =
	 calloc (seqlen1 + seqlen2, sizeof (char));
	 hit1.alignment[2] =
	 calloc (seqlen1 + seqlen2, sizeof (char));
	 hit1.rest[0] =  calloc (30, sizeof (char));
	 hit1.rest[1] =  calloc (30, sizeof (char));
	 hit1.rest[2] =  calloc (30, sizeof (char));
	 hit1.rest[3] =  calloc (30, sizeof (char));
	 hit1.rest[4] =  calloc (30, sizeof (char));
	 hit1.rest[5] =  calloc (30, sizeof (char));
	 */
	
	
	int seqlen2_eff = sequence_length_effective(reference);
	
	get_nt_nt_seq_scores (nt_nt_score,query,reference,seqlen1,seqlen2);
	score=build_matrix_prof (best, track, a_nt_nt, b_gap_nt,c_nt_gap,nt_nt_score, query, 
							 reference, seqlen1, seqlen2, scores); 
	
	if (debug){
		fprintf(fpout,"Score Raw: %f",score);
	}
	if (norm_scores){
        score=normalise_score(score,seqlen1,seqlen2_eff);
	}
	if (debug){
		fprintf(fpout,"\tScore Norm: %f [Sequence1: %d  Sequence2: %d] [%d]\n",score,seqlen1,seqlen2_eff);
	}
	return score;
}
