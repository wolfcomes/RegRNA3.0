/* scmatrix.h									  */
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


/* Scoring Matrix - Only Takes Normal Nucleotides
* Doron - changed G:U wobble weights 2->0 */
//int match[6][6] = {
//		/* C   G   A   T   U   X  */
//  /* C */ {-3, +5, -3, -3, -3, -1},
//  /* G */ {+5, -3, -3, 2, 2, -1},
//  /* A */ {-3, -3, -3, +4, +4, -1},
//  /* T */ {-3, 2, +4, -3, -3, -1},
//  /* U */ {-3, 2, +4, -3, -3, -1},
//  /* X */ {-1, -1, -1, -1, -1, -1}
//};

/* testing matrix */
int match[6][6] = {
	  /* C   G   A   T   U   X  */
  /* C */ {-3, +5, -3, -3, -3, -1},
  /* G */ {+5, -3, -3, 1, 1, -1},
  /* A */ {-3, -3, -3, +5, +5, -1},
  /* T */ {-3, 1, +5, -3, -3, -1},
  /* U */ {-3, 1, +5, -3, -3, -1},
  /* X */ {-1, -1, -1, -1, -1, -1}
};


/*Doron - as far as I can see this matrix is never used*/
int match5p_strict[6][6] = {
         /* C   G   A   T   U   X  */
  /* C */ {-1000, +5, -1000, -1000, -1000, -1},
  /* G */ {+5, -1000, -1000,  2,  2, -1},
  /* A */ {-1000, -1000, -1000, +5, +5, -1},
  /* T */ {-1000,  2, +5, -1000, -1000, -1},
  /* U */ {-1000,  2, +5, -1000, -1000, -1},
  /* X */ {-1, -1, -1, -1, -1, -1}
};

/* Doron - Changed match5p to disfavour G:U wobbles at 5' end
//
//int match5p[6][6] = {
//		/* C   G   A   T   U   X  */
//  /* C */ {-3, +5, -3, -3, -3, -1},
//  /* G */ {+5, -3, -3, +2, +2, -1},
//  /* A */ {-3, -3, -3, +5, +5, -1},
//  /* T */ {-3, +2, +5, -3, -3, -1},
//  /* U */ {-3, +2, +5, -3, -3, -1},
//  /* X */ {-1, -1, -1, -1, -1, -1}
//};


int match5p[6][6] = {
	/* C   G   A   T   U   X  */
  /* C */ {-3, +5, -3, -3, -3, -1},
  /* G */ {+5, -3, -3, 1, 1, -1},
  /* A */ {-3, -3, -3, +5, +5, -1},
  /* T */ {-3, 1, +5, -3, -3, -1},
  /* U */ {-3, 1, +5, -3, -3, -1},
  /* X */ {-1, -1, -1, -1, -1, -1}
};


char ali_rep[6][6] = {
         /* C   G   A   T   U   X  */
  /* C */ {' ','|',' ',' ',' ',' '},
  /* G */ {'|',' ',' ',':',':',' '},
  /* A */ {' ',' ',' ','|','|',' '},
  /* T */ {' ',':','|',' ',' ',' '},
  /* U */ {' ',':','|',' ',' ',' '},
  /* X */ {' ',' ',' ',' ',' ',' '}
};



/* Master list of Nucleotide Bases */
unsigned char baselist[6] = { 'C', 'G', 'A', 'T', 'U', 'X' };
unsigned char bases[1000];
