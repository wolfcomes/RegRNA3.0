/* miranda.h									  */
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
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   
*/
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


/* Setup Traceback Flags */
#define DIAG 0
#define UP 1
#define LEFT -1
int CURR;
/* Setup Strands */
#define FORWARD 0
#define REVERSE 1

/* Setup lengths */
#define IDENTIFIER_LENGTH 1000
#define SEQUENCE_LENGTH 100000

/* Adjustable Algorithm Parameters */
double scale;
int strict;
int debug;
int norm_scores;
double gap_open;
double gap_extend;
double average;
double stdev;
double z_threshold;
double p_threshold;
double score_threshold;
double norm_score_threshold;
double energy_threshold;
double evd_location;
double evd_scale;
int length_5p_for_weighting;
int length_3p_for_weighting;
int key_value_pairs;
int do_shuffle;
int do_profile;
int do_pvalues;
int no_energy;
int shuffle_window;
int total_shuffles;
int verbosity;
unsigned int uniform;
int outfile;
int truncated;
int total_hits;
int restricted;
int seqlen1_global;
int overlap_cutoff;
int max_input_profs;
/* Output Filehandle */
FILE *fpout;


/* Generic structure to store individual hit information */
typedef struct hit_struct
{

  double score;
  double cdf;
  int query_start;
  int query_end;
  int ref_start;
  int ref_end;
  char *alignment[3];
  char *rest[6];

}
hit_struct;


typedef struct mirna_profile
{
        double  *scores;
        double  *cdf;
        int     n_scores;
        int     n_scores_max;
        double  location;
        double  scale;
        char    mirna_id[IDENTIFIER_LENGTH];
}
mirna_profile;

typedef struct species_profile
{
        mirna_profile *profiles;
        int            n_profiles;
        int            n_profiles_max;
        char          *species_id; 
}
species_profile;

typedef struct hit_list
{
	hit_struct *hits;
	long	    n_hits;
	long 	    n_hits_max;
}
hit_list;

/* Generic structure to store consecutive UTR hit information */
typedef struct final_score
{
  double scan_score;
  int no_hits;
  double max_hit;
  double max_score;
  double total_score;
  char positional[1000];
}
final_score;

/* Score structure for non-optimal path detection */
typedef struct score_struct
{

  int score;
  int path;
  int i;
  int j;

}
score_struct;

/* Structure for pair-wise restriction */
typedef struct pairs_struct
{

  char identifier1[200];
  char identifier2[200];

}
pair_struct;

/* Declare Functions */
double max (double, double, double);
int max_finder_fourstates (int, int, int);
int max_finder_threestates (int, int);
int max_finder_and_track_threestates (int, int, int);

int score5p (char, char);
int score (char, char);

void print2dmatrix (int **, int ,int);
void print1dmatrix (int *, int);

void traceback (int **, int ***, char *, char *, int, int, int,
		hit_struct *, double);
int testfor_overlap (int *,int *,int *,int,int);

long readinseq (long, FILE *, char *, char *, char *);
int  readinpairs(FILE *, pair_struct **);
void shuffle (char *, int, int);
int build_matrix (int **, int ***, int **, int **,int **,int **, char *, char *, int, int,
		   score_struct *);
int build_matrix_prof (int **, int ***, int **, int **,int **,int **, char *, char *, int, int,
           score_struct *);

double linear_regression (double *, double *, int, double *, double *);
double graph_cdf(double *, double *, int, double , double , FILE *);


void get_nt_nt_seq_scores (int **,char *, char *, int, int);
double build_matrix_quick (int **, int **, char *, char *, int, int);



void revstring (char s[]);
int get_max_score (char *,int );
unsigned char complementary(unsigned char);
void irand (int);
int nrand (int);
void printhit (char *, char *, hit_struct *, char *, char *, int, double, double, double, FILE *);
void printhit_tskv (char *, char *, hit_struct *, char *, char *, int, double, double, double, FILE *);
void print_license();
void nullify_j_overlaps (int,score_struct *);
void print_small_license();
double get_energy (hit_struct *);
double vfold (char *);
void clear_hit (hit_struct *, int, int);
void clear_matrix (double **, int, int, int, int);
void initialize_bases ();
int cmpscores (const void *, const void *);
int cmp_dbl_asc (const void *p1, const void *p2);
int cmpprofiles(const void *, const void *);
int find_targets (FILE *, FILE *, FILE *, FILE *, FILE *, char *, mirna_profile *, int, int, char *);
int sort_pairs(const void *, const void *);
double do_alignment (int **, int ***, int **, int **, int **, int **, char *, char *,
		     score_struct *, int, int, int, int,
		     final_score *, int, char *, char *, hit_struct *, FILE *, int, double, double);
double do_mirna_utr_profile (int **, int ***, int **, int **, int **, int **, char *, char *,
		     score_struct *, int, int, int, int, FILE *);
void print_banner();
void print_usage();
void print_options();
