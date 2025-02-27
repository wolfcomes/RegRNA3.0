/* statistics.c */
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
#include <ctype.h>
#include <time.h>
#include <math.h>
#include "utils.h"
#include "miranda.h"

void
shuffle (char *sequence, int seqlen, int wsize)
{
  int k = 0;
  char tmp;
  char *top;
  int i, j = 0;
  int mm = seqlen % wsize;

  for (k = 0; k < (seqlen - wsize) + 1; k += wsize)
    {
      top = &sequence[k];
      for (i = wsize; i > 0; i--)
	{
	  j = nrand (i);
	  tmp = top[j];
	  top[j] = top[i - 1];
	  top[i - 1] = tmp;
	}
    }
  top = &sequence[seqlen - mm];
  for (i = mm; i > 0; i--)
    {
      j = nrand (i);
      tmp = top[j];
      top[j] = top[i - 1];
      top[i - 1] = tmp;
    }


}

double
linear_regression(double *x_vals, double *y_vals, int n_vals, double *slopep, double *interceptp){

        double sumx = 0;
        double sumy = 0;
        double sumxx = 0;
        double sumyy = 0;
        double sumxy = 0;
        double censor = 2;
        int n_crap =  (n_vals * censor) / 100; /* number considered */

        int nc = 0;
        double correlation = 0.0;
        double slope = 0.0;
        double intercept = 0.0;
        int i;
        double x=0;
        double y=0;

        for (i=n_crap;i<n_vals-n_crap;i++) { 

                x = x_vals[i];
                y = y_vals[i];

                if (x < 4.0) {
                  continue;
                }

                y=y;
                sumx += x;
                sumy += y;
                sumxx += x * x;
                sumyy += y * y;
                sumxy += x * y;
                nc++;
        }


        if (nc) {
                slope =     (nc * sumxy - sumx * sumy)
                         /  (nc*sumxx - sumx * sumx ); 
                intercept   =    (sumy - slope * sumx) / nc; 

                correlation =        ( nc * sumxy - sumx * sumy )
                                /  sqrt( (nc*sumxx - sumx* sumx) * ( nc * sumyy - sumy * sumy ) );
        }
        else {
                fprintf(stderr,"Panic: Data Bad, Me not like, Me terminate now\n");
        }

        *slopep = slope;
        *interceptp = intercept;
        return (correlation);
}

double graph_cdf(double *x_vals, double *y_vals, int n_vals, double location, double scale, FILE *fpout){

int i=0;
int j=0;
int height=16;
int width=100;
int hist_height=height;
double censor = 2;
int n_crap = (n_vals * censor) / 100; /* number considered */
char plot_space[height+1][width];
char hist_space[hist_height+1][width];

double histogram_emp[width];
double histogram_evd[width];

double max_x=0;
double max_y=0;
double min_x=10000;
double min_y=10000;


int    n_bins=20;
double bin=0;
double bin_size=0;
double bin_start=0;
double bin_end=0;

double max_p=0;
int ticks=0;
int ticks_done=0;
double tick_size=0;
for (i=0;i< n_vals;i++){
       if (x_vals[i] >= max_x){
                max_x=x_vals[i];
       } 
       if (x_vals[i] <= min_x){
                min_x=x_vals[i];
       }
       if (y_vals[i] >= max_y){
                max_y=y_vals[i];
       } 
       if (y_vals[i] <= min_y){
                min_y=y_vals[i];
       }
}

bin_start=min_x;
bin_end=max_x;
bin_size=(bin_end-bin_start)/(n_bins);


fprintf(fpout,"\n");
fprintf(fpout,"\tEmpirical Cumulative Density Plot [(1-log(-log(P(s>=x)))) vs x]\n");
fprintf(fpout,"\t=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
for (i=0;i<=height;i++){
        for (j=0;j<width;j++){
                plot_space[i][j]=' ';
        }
}

for (i=0;i<n_vals;i++){
        int new_x=((x_vals[i]-min_x)/(max_x-min_x))*(height-1);
        int new_y=((y_vals[i]-min_y)/(max_y-min_y))*(width-1);
        if ((i > n_crap) && (i < (n_vals - n_crap))){
                plot_space[new_x][new_y]='+';
        } else {
                plot_space[new_x][new_y]='.';
        }
}

ticks=(height/4);
ticks_done=0;
tick_size=(max_y-min_y)/ticks;

for (i=0;i<height;i++){

        if (!(i%ticks) || (i==(height-1))){
                fprintf(fpout,"\t %5.2f |",min_y+((double) ticks_done*tick_size));
                ticks_done++;
        } else {
                fprintf(fpout,"\t       |");
        }

        for (j=(width-1);j>=0;j--){
                fprintf(fpout,"%c",plot_space[i][j]);
        }
        fprintf(fpout,"\n");
}

fprintf(fpout,"\t       |");
for (j=0;j<width;j++){
        fprintf(fpout,"=");
}
fprintf(fpout,"\n");

ticks=(width/10);
ticks_done=0;
tick_size=(max_x-min_x)/ticks;

fprintf(fpout,"\t       ");
for (j=0;j<width;j++){
        if (!(j%ticks) || (j==(width-1))){
                fprintf(fpout,"|%6.3f   ",min_x+((double) ticks_done*tick_size));
                ticks_done++;
        } else {
                /* fprintf(fpout,""); */
        }

}

fprintf(fpout,"\n\n\n\n");
fprintf(fpout,"\tExtreme Value Distribution Fit [Location: %6.3f Scale: %5.3f]\n",location,scale);
fprintf(fpout,"\t=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
for (i=0;i<n_bins;i++){
        double t1=0;
        double t2=0;
        bin=bin_start+(i*bin_size);
        histogram_emp[i]=0;
        histogram_evd[i]=0;
        for (j=0;j<n_vals;j++){
                if ( (x_vals[j] >= bin) && (x_vals[j] <= (bin+bin_size)) ){ 
                       histogram_emp[i]++; 
                }
        }
        histogram_emp[i]=(histogram_emp[i]/n_vals);
	    t1=1-((exp(-exp(-(((bin)-location)/scale)))));
	    t2=((exp(-exp(-(((bin+bin_size)-location)/scale)))));
        histogram_evd[i]=(t1+t2)-1;

	if (histogram_evd[i] >= max_p){
		max_p=histogram_evd[i];
	}
	if (histogram_emp[i] >= max_p){
		max_p=histogram_emp[i];
	}
}

for (j=0;j<n_bins;j++){

	int this_height=(histogram_emp[j]/max_p) * hist_height;
        int evd_height= (histogram_evd[j]/max_p) * hist_height;
	int next_height=evd_height;
    char symbol;
	if (j<(n_bins-1)){
		next_height=(histogram_evd[j+1]/max_p) * hist_height;
	}
	for (i=0;i<=hist_height;i++){
		if (i<=this_height){
			hist_space[i][j]='#';
		} else {
                	hist_space[i][j]=' ';
		}
        }

	symbol='#';
	if (next_height > evd_height){
		symbol='/';
	} 
	if (next_height == evd_height){
		symbol='-';
	} 
	if (next_height < evd_height){
		symbol='\\';
	}
        hist_space[evd_height][j]=symbol;
}

for (i=(hist_height);i>=0;i--){
	int ticks=(hist_height/4);
	if (!(i%ticks)){
		fprintf(fpout,"\t %5.3f |",((double)i/(double) hist_height)*max_p);
	} else {
		fprintf(fpout,"\t       |");
	}
	for (j=0;j<n_bins;j++){
                fprintf(fpout,"%-5c", hist_space[i][j]);
        }
        fprintf(fpout,"\n");
}
fprintf(fpout,"\t       |");
for (j=0;j<width;j++){
        fprintf(fpout,"=");
}
fprintf(fpout,"\n");
fprintf(fpout,"\t       ");
for (j=0;j<n_bins;j++){
	fprintf(fpout,"%-5.1f",bin_start+(j*bin_size));
}
fprintf(fpout,"\n\n\n\n");

return(1);
}

void
irand (int n)
{				/* initialize random number generator */
  if (n == 0)
    {
      n = time (NULL);
      n = n % 16381;
      if ((n % 2) == 0)
	n++;
    }
  srand48 (n);
}


int
nrand (int n)
{				/* returns a random number between 0 and n-1
				 * where n < 64K) */
  int rn;
  rn = lrand48 ();
  rn = rn >> 16;
  rn = (rn % n);
  return rn;
}

int
getfreq (char *sequence, int seqlen, double *frequency)
{

  int i = 0;
  for (i = 0; i < 256; i++)
    {
      frequency[i] = 0;
    }

  for (i = 0; i < seqlen; i++)
    {
      frequency[toupper (sequence[i])]++;
    }

  for (i = 0; i < 256; i++)
    {
      frequency[i] = frequency[i] / seqlen;
    }

  return (1);
}
