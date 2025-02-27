/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     SYM_PARMS = 258,
     SYM_DESCR = 259,
     SYM_SITES = 260,
     SYM_SCORE = 261,
     SYM_SE = 262,
     SYM_CTX = 263,
     SYM_SS = 264,
     SYM_H5 = 265,
     SYM_H3 = 266,
     SYM_P5 = 267,
     SYM_P3 = 268,
     SYM_T1 = 269,
     SYM_T2 = 270,
     SYM_T3 = 271,
     SYM_Q1 = 272,
     SYM_Q2 = 273,
     SYM_Q3 = 274,
     SYM_Q4 = 275,
     SYM_ACCEPT = 276,
     SYM_BEGIN = 277,
     SYM_BREAK = 278,
     SYM_CONTINUE = 279,
     SYM_ELSE = 280,
     SYM_END = 281,
     SYM_FOR = 282,
     SYM_HOLD = 283,
     SYM_IF = 284,
     SYM_IN = 285,
     SYM_REJECT = 286,
     SYM_RELEASE = 287,
     SYM_WHILE = 288,
     SYM_IDENT = 289,
     SYM_INT = 290,
     SYM_FLOAT = 291,
     SYM_STRING = 292,
     SYM_PAIRSET = 293,
     SYM_AND = 294,
     SYM_ASSIGN = 295,
     SYM_DOLLAR = 296,
     SYM_DONT_MATCH = 297,
     SYM_EQUAL = 298,
     SYM_GREATER = 299,
     SYM_GREATER_EQUAL = 300,
     SYM_LESS = 301,
     SYM_LESS_EQUAL = 302,
     SYM_MATCH = 303,
     SYM_MINUS = 304,
     SYM_MINUS_ASSIGN = 305,
     SYM_MINUS_MINUS = 306,
     SYM_NEGATE = 307,
     SYM_NOT = 308,
     SYM_NOT_EQUAL = 309,
     SYM_OR = 310,
     SYM_PERCENT = 311,
     SYM_PERCENT_ASSIGN = 312,
     SYM_PLUS = 313,
     SYM_PLUS_ASSIGN = 314,
     SYM_PLUS_PLUS = 315,
     SYM_STAR = 316,
     SYM_STAR_ASSIGN = 317,
     SYM_SLASH = 318,
     SYM_SLASH_ASSIGN = 319,
     SYM_LPAREN = 320,
     SYM_RPAREN = 321,
     SYM_LBRACK = 322,
     SYM_RBRACK = 323,
     SYM_LCURLY = 324,
     SYM_RCURLY = 325,
     SYM_COLON = 326,
     SYM_COMMA = 327,
     SYM_SEMICOLON = 328,
     SYM_CALL = 329,
     SYM_LIST = 330,
     SYM_KW_STREF = 331,
     SYM_IX_STREF = 332,
     SYM_ERROR = 333
   };
#endif
/* Tokens.  */
#define SYM_PARMS 258
#define SYM_DESCR 259
#define SYM_SITES 260
#define SYM_SCORE 261
#define SYM_SE 262
#define SYM_CTX 263
#define SYM_SS 264
#define SYM_H5 265
#define SYM_H3 266
#define SYM_P5 267
#define SYM_P3 268
#define SYM_T1 269
#define SYM_T2 270
#define SYM_T3 271
#define SYM_Q1 272
#define SYM_Q2 273
#define SYM_Q3 274
#define SYM_Q4 275
#define SYM_ACCEPT 276
#define SYM_BEGIN 277
#define SYM_BREAK 278
#define SYM_CONTINUE 279
#define SYM_ELSE 280
#define SYM_END 281
#define SYM_FOR 282
#define SYM_HOLD 283
#define SYM_IF 284
#define SYM_IN 285
#define SYM_REJECT 286
#define SYM_RELEASE 287
#define SYM_WHILE 288
#define SYM_IDENT 289
#define SYM_INT 290
#define SYM_FLOAT 291
#define SYM_STRING 292
#define SYM_PAIRSET 293
#define SYM_AND 294
#define SYM_ASSIGN 295
#define SYM_DOLLAR 296
#define SYM_DONT_MATCH 297
#define SYM_EQUAL 298
#define SYM_GREATER 299
#define SYM_GREATER_EQUAL 300
#define SYM_LESS 301
#define SYM_LESS_EQUAL 302
#define SYM_MATCH 303
#define SYM_MINUS 304
#define SYM_MINUS_ASSIGN 305
#define SYM_MINUS_MINUS 306
#define SYM_NEGATE 307
#define SYM_NOT 308
#define SYM_NOT_EQUAL 309
#define SYM_OR 310
#define SYM_PERCENT 311
#define SYM_PERCENT_ASSIGN 312
#define SYM_PLUS 313
#define SYM_PLUS_ASSIGN 314
#define SYM_PLUS_PLUS 315
#define SYM_STAR 316
#define SYM_STAR_ASSIGN 317
#define SYM_SLASH 318
#define SYM_SLASH_ASSIGN 319
#define SYM_LPAREN 320
#define SYM_RPAREN 321
#define SYM_LBRACK 322
#define SYM_RBRACK 323
#define SYM_LCURLY 324
#define SYM_RCURLY 325
#define SYM_COLON 326
#define SYM_COMMA 327
#define SYM_SEMICOLON 328
#define SYM_CALL 329
#define SYM_LIST 330
#define SYM_KW_STREF 331
#define SYM_IX_STREF 332
#define SYM_ERROR 333




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 20 "rmgrm.y"
{
	int	ival;
	NODE_T	*npval;
}
/* Line 1529 of yacc.c.  */
#line 210 "y.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

