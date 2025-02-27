/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



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




/* Copy the first part of user declarations.  */
#line 1 "rmgrm.y"


#include <stdio.h>
#include "rnamot.h"

extern	VALUE_T	rm_tokval;
extern	int	rm_context;

static	NODE_T	*np;

/*
typedef	union	{
	int	ival;
	NODE_T	*npval;
} YYSTYPE;
*/



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 20 "rmgrm.y"
{
	int	ival;
	NODE_T	*npval;
}
/* Line 193 of yacc.c.  */
#line 276 "y.tab.c"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 289 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  6
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   382

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  79
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  83
/* YYNRULES -- Number of rules.  */
#define YYNRULES  171
/* YYNRULES -- Number of states.  */
#define YYNSTATES  247

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   333

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     8,     9,    13,    15,    17,    19,    20,
      24,    25,    29,    30,    31,    35,    36,    38,    41,    44,
      46,    49,    51,    53,    55,    57,    59,    61,    63,    65,
      67,    69,    71,    73,    75,    77,    79,    81,    83,    85,
      88,    92,    96,    98,   101,   102,   106,   108,   110,   112,
     114,   118,   120,   123,   125,   127,   129,   131,   133,   135,
     137,   139,   141,   143,   145,   147,   149,   151,   154,   157,
     160,   164,   167,   171,   175,   178,   181,   185,   188,   189,
     195,   198,   202,   203,   210,   212,   213,   214,   220,   225,
     226,   227,   235,   237,   239,   241,   243,   245,   247,   249,
     251,   253,   257,   261,   263,   265,   267,   269,   271,   273,
     275,   279,   281,   285,   287,   289,   293,   295,   297,   299,
     301,   303,   305,   307,   309,   311,   315,   317,   319,   321,
     325,   327,   329,   331,   333,   336,   339,   341,   343,   347,
     349,   353,   355,   357,   359,   363,   368,   370,   372,   377,
     382,   384,   386,   389,   392,   394,   396,   398,   400,   402,
     404,   406,   408,   410,   414,   416,   420,   421,   426,   428,
     432,   434
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
      80,     0,    -1,    81,    84,    86,    88,    -1,    -1,    83,
      82,    90,    -1,   161,    -1,     3,    -1,   161,    -1,    -1,
       4,    85,    92,    -1,    -1,     5,    87,    96,    -1,    -1,
      -1,     6,    89,    99,    -1,    -1,    91,    -1,    91,    90,
      -1,   132,    73,    -1,    93,    -1,    93,    92,    -1,    94,
      -1,   148,    -1,    95,    -1,     7,    -1,     8,    -1,     9,
      -1,    10,    -1,    11,    -1,    12,    -1,    13,    -1,    14,
      -1,    15,    -1,    16,    -1,    17,    -1,    18,    -1,    19,
      -1,    20,    -1,    97,    -1,    96,    97,    -1,   144,    30,
     157,    -1,   143,    30,   157,    -1,   100,    -1,   100,    99,
      -1,    -1,   102,   101,   103,    -1,   103,    -1,    22,    -1,
      26,    -1,   134,    -1,    69,   104,    70,    -1,   105,    -1,
     105,   104,    -1,   106,    -1,   107,    -1,   108,    -1,   109,
      -1,   110,    -1,   111,    -1,   112,    -1,   113,    -1,   114,
      -1,   115,    -1,   116,    -1,   118,    -1,   119,    -1,   120,
      -1,    21,    73,    -1,   132,    73,    -1,   151,    73,    -1,
      23,   122,    73,    -1,   146,    73,    -1,    69,   104,    70,
      -1,    24,   122,    73,    -1,   161,    73,    -1,   125,   105,
      -1,    28,   153,    73,    -1,   123,   105,    -1,    -1,   123,
     105,    25,   117,   105,    -1,    31,    73,    -1,    32,   153,
      73,    -1,    -1,    33,    65,   134,   121,    66,   105,    -1,
      35,    -1,    -1,    -1,    29,    65,   134,   124,    66,    -1,
      27,    65,   126,    66,    -1,    -1,    -1,   129,   127,    73,
     130,   128,    73,   131,    -1,   132,    -1,   151,    -1,   161,
      -1,   132,    -1,   134,    -1,   161,    -1,   132,    -1,   151,
      -1,   161,    -1,   150,   133,   132,    -1,   150,   133,   134,
      -1,    40,    -1,    50,    -1,    59,    -1,    57,    -1,    64,
      -1,    62,    -1,   135,    -1,   134,    55,   135,    -1,   136,
      -1,   136,    39,   135,    -1,    98,    -1,   138,    -1,   138,
     137,   138,    -1,    42,    -1,    43,    -1,    44,    -1,    45,
      -1,    46,    -1,    47,    -1,    48,    -1,    54,    -1,   140,
      -1,   138,   139,   140,    -1,    58,    -1,    49,    -1,   142,
      -1,   140,   141,   142,    -1,    56,    -1,    63,    -1,    61,
      -1,   145,    -1,    49,   145,    -1,    53,   145,    -1,   147,
      -1,   147,    -1,   147,    71,   143,    -1,   148,    -1,   148,
      71,   144,    -1,   150,    -1,   152,    -1,   146,    -1,    65,
     134,    66,    -1,   153,    65,   155,    66,    -1,   148,    -1,
     149,    -1,    94,    65,   156,    66,    -1,    94,    67,   155,
      68,    -1,   153,    -1,   151,    -1,   154,   153,    -1,   153,
     154,    -1,    35,    -1,    36,    -1,    41,    -1,   160,    -1,
     157,    -1,    34,    -1,    51,    -1,    60,    -1,   134,    -1,
     134,    72,   155,    -1,   132,    -1,   132,    72,   156,    -1,
      -1,    69,   158,   159,    70,    -1,   160,    -1,   160,    72,
     159,    -1,    37,    -1,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   182,   182,   185,   185,   187,   189,   190,   192,   192,
     196,   196,   199,   201,   201,   204,   207,   208,   210,   213,
     214,   215,   220,   222,   230,   231,   232,   233,   234,   235,
     236,   237,   238,   239,   240,   241,   242,   243,   246,   247,
     249,   256,   264,   265,   267,   267,   269,   271,   272,   273,
     275,   278,   279,   281,   282,   283,   284,   285,   286,   287,
     288,   289,   290,   291,   292,   293,   294,   296,   299,   305,
     311,   314,   319,   322,   325,   327,   329,   332,   334,   333,
     338,   341,   345,   344,   349,   350,   353,   352,   357,   360,
     362,   360,   366,   367,   368,   370,   371,   372,   374,   375,
     376,   379,   387,   396,   397,   399,   401,   403,   405,   408,
     409,   412,   413,   416,   417,   418,   421,   423,   424,   425,
     427,   428,   430,   431,   433,   434,   437,   438,   440,   441,
     444,   445,   446,   448,   449,   451,   453,   457,   460,   465,
     468,   473,   474,   475,   476,   479,   482,   483,   485,   494,
     497,   498,   500,   501,   503,   504,   505,   507,   508,   510,
     513,   515,   517,   518,   521,   524,   529,   529,   533,   534,
     537,   540
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "SYM_PARMS", "SYM_DESCR", "SYM_SITES",
  "SYM_SCORE", "SYM_SE", "SYM_CTX", "SYM_SS", "SYM_H5", "SYM_H3", "SYM_P5",
  "SYM_P3", "SYM_T1", "SYM_T2", "SYM_T3", "SYM_Q1", "SYM_Q2", "SYM_Q3",
  "SYM_Q4", "SYM_ACCEPT", "SYM_BEGIN", "SYM_BREAK", "SYM_CONTINUE",
  "SYM_ELSE", "SYM_END", "SYM_FOR", "SYM_HOLD", "SYM_IF", "SYM_IN",
  "SYM_REJECT", "SYM_RELEASE", "SYM_WHILE", "SYM_IDENT", "SYM_INT",
  "SYM_FLOAT", "SYM_STRING", "SYM_PAIRSET", "SYM_AND", "SYM_ASSIGN",
  "SYM_DOLLAR", "SYM_DONT_MATCH", "SYM_EQUAL", "SYM_GREATER",
  "SYM_GREATER_EQUAL", "SYM_LESS", "SYM_LESS_EQUAL", "SYM_MATCH",
  "SYM_MINUS", "SYM_MINUS_ASSIGN", "SYM_MINUS_MINUS", "SYM_NEGATE",
  "SYM_NOT", "SYM_NOT_EQUAL", "SYM_OR", "SYM_PERCENT",
  "SYM_PERCENT_ASSIGN", "SYM_PLUS", "SYM_PLUS_ASSIGN", "SYM_PLUS_PLUS",
  "SYM_STAR", "SYM_STAR_ASSIGN", "SYM_SLASH", "SYM_SLASH_ASSIGN",
  "SYM_LPAREN", "SYM_RPAREN", "SYM_LBRACK", "SYM_RBRACK", "SYM_LCURLY",
  "SYM_RCURLY", "SYM_COLON", "SYM_COMMA", "SYM_SEMICOLON", "SYM_CALL",
  "SYM_LIST", "SYM_KW_STREF", "SYM_IX_STREF", "SYM_ERROR", "$accept",
  "program", "parm_part", "@1", "parm_hdr", "descr_part", "@2",
  "site_part", "@3", "score_part", "@4", "pd_list", "pdef", "se_list",
  "strel", "strhdr", "strtype", "kw_site_list", "kw_site", "site",
  "rule_list", "rule", "@5", "pattern", "action", "stmt_list", "stmt",
  "accept_stmt", "asgn_stmt", "auto_stmt", "break_stmt", "call_stmt",
  "cmpd_stmt", "continue_stmt", "empty_stmt", "for_stmt", "hold_stmt",
  "if_stmt", "@6", "reject_stmt", "release_stmt", "while_stmt", "@7",
  "loop_level", "if_hdr", "@8", "for_hdr", "for_ctrl", "@9", "@10",
  "for_init", "for_test", "for_incr", "asgn", "asgn_op", "expr", "conj",
  "compare", "comp_op", "a_expr", "add_op", "term", "mul_op", "factor",
  "pairing", "kw_pairing", "primary", "fcall", "stref", "kw_stref",
  "ix_stref", "lval", "auto_lval", "literal", "ident", "incr_op", "e_list",
  "a_list", "pairset", "@11", "s_list", "string", "empty", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    79,    80,    82,    81,    81,    83,    83,    85,    84,
      87,    86,    86,    89,    88,    88,    90,    90,    91,    92,
      92,    93,    93,    94,    95,    95,    95,    95,    95,    95,
      95,    95,    95,    95,    95,    95,    95,    95,    96,    96,
      97,    98,    99,    99,   101,   100,   100,   102,   102,   102,
     103,   104,   104,   105,   105,   105,   105,   105,   105,   105,
     105,   105,   105,   105,   105,   105,   105,   106,   107,   108,
     109,   110,   111,   112,   113,   114,   115,   116,   117,   116,
     118,   119,   121,   120,   122,   122,   124,   123,   125,   127,
     128,   126,   129,   129,   129,   130,   130,   130,   131,   131,
     131,   132,   132,   133,   133,   133,   133,   133,   133,   134,
     134,   135,   135,   136,   136,   136,   137,   137,   137,   137,
     137,   137,   137,   137,   138,   138,   139,   139,   140,   140,
     141,   141,   141,   142,   142,   142,   142,   143,   143,   144,
     144,   145,   145,   145,   145,   146,   147,   147,   148,   149,
     150,   150,   151,   151,   152,   152,   152,   152,   152,   153,
     154,   154,   155,   155,   156,   156,   158,   157,   159,   159,
     160,   161
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     4,     0,     3,     1,     1,     1,     0,     3,
       0,     3,     0,     0,     3,     0,     1,     2,     2,     1,
       2,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     2,
       3,     3,     1,     2,     0,     3,     1,     1,     1,     1,
       3,     1,     2,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     2,     2,     2,
       3,     2,     3,     3,     2,     2,     3,     2,     0,     5,
       2,     3,     0,     6,     1,     0,     0,     5,     4,     0,
       0,     7,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     3,     3,     1,     1,     1,     1,     1,     1,     1,
       3,     1,     3,     1,     1,     3,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     3,     1,     1,     1,     3,
       1,     1,     1,     1,     2,     2,     1,     1,     3,     1,
       3,     1,     1,     1,     3,     4,     1,     1,     4,     4,
       1,     1,     2,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     3,     1,     3,     0,     4,     1,     3,
       1,     0
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
     171,     6,     0,     0,     3,     7,     1,     8,    12,     0,
       0,    10,    15,   159,   160,   161,     4,    16,     0,     0,
     151,   150,     0,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,     9,    19,    21,
      23,    22,     0,    13,     2,    17,    18,   103,   104,   106,
     105,   108,   107,     0,   153,   152,    20,     0,     0,    11,
      38,     0,   139,     0,   154,   155,   170,   156,     0,     0,
       0,   166,     0,   113,   101,   102,   109,   111,   114,   124,
     128,     0,   133,   143,   136,   146,   147,   141,   142,   150,
     158,   157,   164,     0,    39,     0,     0,    47,    48,   166,
      14,    42,    44,    46,    49,   141,   134,   135,     0,     0,
       0,     0,     0,   116,   117,   118,   119,   120,   121,   122,
     127,   123,   126,     0,     0,   130,   132,   131,     0,     0,
       0,     0,     0,   148,    40,   140,     0,    85,    85,     0,
       0,     0,     0,     0,     0,   171,     0,    51,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,   171,   171,     0,     0,   151,     0,    43,     0,
     144,     0,   168,   162,     0,   110,   112,   115,   136,   125,
     129,    41,   138,   137,     0,   165,    67,    84,     0,     0,
     171,     0,     0,    80,     0,     0,     0,    50,    52,    77,
      75,    68,    71,    69,    74,   171,    45,   167,     0,     0,
     149,   145,    70,    73,     0,    89,    92,   151,    94,    76,
      86,    81,    82,    72,    78,   169,   163,    88,     0,     0,
       0,   171,   171,    87,   171,    79,    90,    95,    96,    97,
      83,     0,   171,    91,    98,   151,   100
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     2,     3,     9,     4,     8,    10,    12,    42,    44,
      63,    16,    17,    37,    38,    72,    40,    59,    60,    73,
     100,   101,   169,   102,   103,   146,   147,   148,   149,   150,
     151,   152,   153,   154,   155,   156,   157,   158,   231,   159,
     160,   161,   230,   188,   162,   229,   163,   214,   228,   241,
     215,   236,   243,   164,    53,   173,    76,    77,   123,    78,
     124,    79,   128,    80,    81,    61,    82,    83,    84,    85,
      86,   105,    20,    88,    89,    22,   174,    93,    90,   109,
     171,    91,   167
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -145
static const yytype_int16 yypact[] =
{
      30,  -145,    40,    54,  -145,    62,  -145,  -145,    65,   -13,
     362,  -145,    38,  -145,  -145,  -145,  -145,   -13,    -6,    39,
    -145,   -25,    44,  -145,  -145,  -145,  -145,  -145,  -145,  -145,
    -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,   362,    15,
    -145,  -145,   362,  -145,  -145,  -145,  -145,  -145,  -145,  -145,
    -145,  -145,  -145,   299,  -145,  -145,  -145,   -13,    15,   362,
    -145,    52,    13,   236,  -145,  -145,  -145,  -145,   140,   140,
     299,  -145,   -40,  -145,  -145,    37,  -145,    46,   171,    25,
    -145,    63,  -145,  -145,   -14,  -145,  -145,    39,  -145,   -29,
    -145,  -145,    22,    29,  -145,    35,   362,  -145,  -145,    84,
    -145,   236,  -145,  -145,    37,  -145,  -145,  -145,   -32,    69,
     299,   299,   299,  -145,  -145,  -145,  -145,  -145,  -145,  -145,
    -145,  -145,  -145,   299,   299,  -145,  -145,  -145,   299,    35,
     362,   299,   -13,  -145,  -145,  -145,    36,    75,    75,    55,
      44,    56,    49,    44,    58,   138,    59,    84,  -145,  -145,
    -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,
    -145,  -145,   138,   138,    53,    57,    61,    68,  -145,    73,
    -145,    77,    60,   -35,    70,  -145,  -145,   -17,  -145,    25,
    -145,  -145,  -145,    79,    86,  -145,  -145,  -145,    83,    85,
     -13,    87,   299,  -145,    90,   299,    94,  -145,  -145,   103,
    -145,  -145,  -145,  -145,  -145,   138,  -145,  -145,    69,   299,
    -145,  -145,  -145,  -145,   102,  -145,  -145,   100,  -145,  -145,
      37,  -145,    37,  -145,  -145,  -145,  -145,  -145,   105,   113,
     114,   138,   299,  -145,   138,  -145,  -145,  -145,    37,  -145,
    -145,   111,   -13,  -145,  -145,   121,  -145
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,
    -145,   175,  -145,    93,  -145,     1,  -145,  -145,   129,  -145,
      92,  -145,  -145,  -145,    26,   -94,  -144,  -145,  -145,  -145,
    -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,
    -145,  -145,  -145,    64,  -145,  -145,  -145,  -145,  -145,  -145,
    -145,  -145,  -145,    -5,  -145,   -46,   -82,  -145,  -145,    71,
    -145,    72,  -145,    76,    78,   107,    -4,   -86,   -55,     4,
    -145,    -8,   -91,  -145,    -7,   -15,  -126,    67,   -67,  -145,
      -2,  -106,     0
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -172
static const yytype_int16 yytable[] =
{
       5,    19,    21,   172,    18,   184,    54,    75,   166,    19,
      21,    39,    18,   165,    41,    55,  -137,   104,   199,   200,
     111,    13,    14,   111,   108,    57,    14,   110,   134,   175,
     176,    15,   120,     1,   170,    15,   131,   209,    14,    39,
       6,   122,    41,    58,    43,    87,    62,    15,    74,    19,
      21,   196,    92,   198,   166,   104,   166,   130,     7,   165,
      58,   165,   181,    62,   106,   107,    -5,    46,   178,   178,
      11,   166,   166,   178,    54,   183,   165,   165,    13,    47,
      57,   125,    95,   226,    96,   112,   126,   235,   127,    48,
     240,    19,   111,   129,   132,   133,    49,    58,    50,   217,
      62,    51,   172,    52,    71,   136,    66,   137,   138,   186,
     187,   139,   140,   141,   166,   142,   143,   144,    13,   165,
     190,   192,   193,   195,    19,    21,   201,    92,   224,   197,
     202,    56,   208,   191,   203,    14,   194,    19,   210,    19,
     166,   204,   205,   166,    15,   165,   220,   207,   165,   222,
     130,   245,   211,   145,    19,    19,   212,  -171,   213,   136,
     219,   137,   138,   221,   223,   139,   140,   141,   227,   142,
     143,   144,    13,   -93,    13,    64,    65,    66,   232,   233,
     234,    67,    19,    21,   242,   216,   238,   -99,    94,    14,
     218,    14,    45,   168,   177,   206,   179,    19,    15,   185,
      15,     0,   189,   135,   180,    70,   225,   145,   182,    71,
       0,     0,     0,   113,   114,   115,   116,   117,   118,   119,
     120,     0,     0,    19,    87,   121,    19,   237,     0,   122,
       0,     0,   239,     0,    19,    21,     0,   244,     0,     0,
       0,     0,   246,    23,    24,    25,    26,    27,    28,    29,
      30,    31,    32,    33,    34,    35,    36,     0,    97,     0,
       0,     0,    98,     0,     0,     0,     0,     0,     0,     0,
      13,    64,    65,    66,     0,     0,     0,    67,     0,     0,
       0,     0,     0,     0,     0,    68,     0,    14,     0,    69,
       0,     0,     0,     0,     0,     0,    15,     0,     0,     0,
       0,    70,     0,     0,     0,    99,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    33,    34,    35,    36,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    13,    64,    65,    66,     0,     0,     0,
      67,     0,     0,     0,     0,     0,     0,     0,    68,     0,
      14,     0,    69,     0,     0,     0,     0,     0,     0,    15,
       0,     0,     0,     0,    70,     0,     0,     0,    71,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36
};

static const yytype_int16 yycheck[] =
{
       0,     9,     9,   109,     9,   131,    21,    53,    99,    17,
      17,    10,    17,    99,    10,    22,    30,    63,   162,   163,
      55,    34,    51,    55,    70,    65,    51,    67,    95,   111,
     112,    60,    49,     3,    66,    60,    65,    72,    51,    38,
       0,    58,    38,    42,     6,    53,    42,    60,    53,    57,
      57,   145,    57,   147,   145,   101,   147,    71,     4,   145,
      59,   147,   129,    59,    68,    69,     4,    73,   123,   124,
       5,   162,   163,   128,    89,   130,   162,   163,    34,    40,
      65,    56,    30,   209,    71,    39,    61,   231,    63,    50,
     234,    99,    55,    30,    72,    66,    57,    96,    59,   190,
      96,    62,   208,    64,    69,    21,    37,    23,    24,    73,
      35,    27,    28,    29,   205,    31,    32,    33,    34,   205,
      65,    65,    73,    65,   132,   132,    73,   132,    25,    70,
      73,    38,    72,   140,    73,    51,   143,   145,    68,   147,
     231,    73,    69,   234,    60,   231,   192,    70,   234,   195,
      71,   242,    66,    69,   162,   163,    73,    73,    73,    21,
      73,    23,    24,    73,    70,    27,    28,    29,    66,    31,
      32,    33,    34,    73,    34,    35,    36,    37,    73,    66,
      66,    41,   190,   190,    73,   190,   232,    66,    59,    51,
     190,    51,    17,   101,   123,   169,   124,   205,    60,   132,
      60,    -1,   138,    96,   128,    65,   208,    69,   130,    69,
      -1,    -1,    -1,    42,    43,    44,    45,    46,    47,    48,
      49,    -1,    -1,   231,   232,    54,   234,   232,    -1,    58,
      -1,    -1,   232,    -1,   242,   242,    -1,   242,    -1,    -1,
      -1,    -1,   242,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    -1,    22,    -1,
      -1,    -1,    26,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      34,    35,    36,    37,    -1,    -1,    -1,    41,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    49,    -1,    51,    -1,    53,
      -1,    -1,    -1,    -1,    -1,    -1,    60,    -1,    -1,    -1,
      -1,    65,    -1,    -1,    -1,    69,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    34,    35,    36,    37,    -1,    -1,    -1,
      41,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    49,    -1,
      51,    -1,    53,    -1,    -1,    -1,    -1,    -1,    -1,    60,
      -1,    -1,    -1,    -1,    65,    -1,    -1,    -1,    69,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,    80,    81,    83,   161,     0,     4,    84,    82,
      85,     5,    86,    34,    51,    60,    90,    91,   132,   150,
     151,   153,   154,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    92,    93,    94,
      95,   148,    87,     6,    88,    90,    73,    40,    50,    57,
      59,    62,    64,   133,   154,   153,    92,    65,    94,    96,
      97,   144,   148,    89,    35,    36,    37,    41,    49,    53,
      65,    69,    94,    98,   132,   134,   135,   136,   138,   140,
     142,   143,   145,   146,   147,   148,   149,   150,   152,   153,
     157,   160,   132,   156,    97,    30,    71,    22,    26,    69,
      99,   100,   102,   103,   134,   150,   145,   145,   134,   158,
      67,    55,    39,    42,    43,    44,    45,    46,    47,    48,
      49,    54,    58,   137,   139,    56,    61,    63,   141,    30,
      71,    65,    72,    66,   157,   144,    21,    23,    24,    27,
      28,    29,    31,    32,    33,    69,   104,   105,   106,   107,
     108,   109,   110,   111,   112,   113,   114,   115,   116,   118,
     119,   120,   123,   125,   132,   146,   151,   161,    99,   101,
      66,   159,   160,   134,   155,   135,   135,   138,   147,   140,
     142,   157,   143,   147,   155,   156,    73,    35,   122,   122,
      65,   153,    65,    73,   153,    65,   104,    70,   104,   105,
     105,    73,    73,    73,    73,    69,   103,    70,    72,    72,
      68,    66,    73,    73,   126,   129,   132,   151,   161,    73,
     134,    73,   134,    70,    25,   159,   155,    66,   127,   124,
     121,   117,    73,    66,    66,   105,   130,   132,   134,   161,
     105,   128,    73,   131,   132,   151,   161
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 3:
#line 185 "rmgrm.y"
    { rm_context = CTX_PARMS; }
    break;

  case 6:
#line 189 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;

  case 7:
#line 190 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;

  case 8:
#line 192 "rmgrm.y"
    { rm_context = CTX_DESCR; }
    break;

  case 9:
#line 194 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;

  case 10:
#line 196 "rmgrm.y"
    { rm_context = CTX_SITES; }
    break;

  case 11:
#line 198 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;

  case 12:
#line 199 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;

  case 13:
#line 201 "rmgrm.y"
    { rm_context = CTX_SCORE; }
    break;

  case 14:
#line 203 "rmgrm.y"
    { RM_accept(); (yyval.npval) = NULL; }
    break;

  case 15:
#line 204 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;

  case 21:
#line 215 "rmgrm.y"
    { if( rm_context == CTX_DESCR )
					SE_close();
				  else if( rm_context == CTX_SITES )
					POS_close();
				}
    break;

  case 23:
#line 222 "rmgrm.y"
    { if( rm_context == CTX_DESCR )
					SE_open( (yyvsp[(1) - (1)].ival) );
				  else if( rm_context == CTX_SITES )
					POS_open( (yyvsp[(1) - (1)].ival) );
				  else
					(yyval.npval) = RM_node( (yyvsp[(1) - (1)].ival), 0, 0, 0 );
				}
    break;

  case 24:
#line 230 "rmgrm.y"
    { (yyval.ival) = SYM_SE; }
    break;

  case 25:
#line 231 "rmgrm.y"
    { (yyval.ival) = SYM_CTX; }
    break;

  case 26:
#line 232 "rmgrm.y"
    { (yyval.ival) = SYM_SS; }
    break;

  case 27:
#line 233 "rmgrm.y"
    { (yyval.ival) = SYM_H5; }
    break;

  case 28:
#line 234 "rmgrm.y"
    { (yyval.ival) = SYM_H3; }
    break;

  case 29:
#line 235 "rmgrm.y"
    { (yyval.ival) = SYM_P5; }
    break;

  case 30:
#line 236 "rmgrm.y"
    { (yyval.ival) = SYM_P3; }
    break;

  case 31:
#line 237 "rmgrm.y"
    { (yyval.ival) = SYM_T1; }
    break;

  case 32:
#line 238 "rmgrm.y"
    { (yyval.ival) = SYM_T2; }
    break;

  case 33:
#line 239 "rmgrm.y"
    { (yyval.ival) = SYM_T3; }
    break;

  case 34:
#line 240 "rmgrm.y"
    { (yyval.ival) = SYM_Q1; }
    break;

  case 35:
#line 241 "rmgrm.y"
    { (yyval.ival) = SYM_Q2; }
    break;

  case 36:
#line 242 "rmgrm.y"
    { (yyval.ival) = SYM_Q3; }
    break;

  case 37:
#line 243 "rmgrm.y"
    { (yyval.ival) = SYM_Q4; }
    break;

  case 40:
#line 250 "rmgrm.y"
    { if( rm_context == CTX_SITES )
					SI_close( (yyvsp[(3) - (3)].npval) );
				  else if( rm_context == CTX_SCORE )
					(yyval.npval) = RM_node( SYM_IN, 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) );
				}
    break;

  case 41:
#line 257 "rmgrm.y"
    { if( rm_context == CTX_SITES )
					SI_close( (yyvsp[(3) - (3)].npval) );
				  else if( rm_context == CTX_SCORE )
					(yyval.npval) = RM_node( SYM_IN, 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) );
				}
    break;

  case 44:
#line 267 "rmgrm.y"
    { RM_action( (yyvsp[(1) - (1)].npval) ); }
    break;

  case 45:
#line 268 "rmgrm.y"
    { RM_endaction(); }
    break;

  case 47:
#line 271 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_BEGIN, 0, 0, 0 ); }
    break;

  case 48:
#line 272 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_END, 0, 0, 0 ); }
    break;

  case 49:
#line 273 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 50:
#line 276 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;

  case 67:
#line 297 "rmgrm.y"
    { RM_accept(); (yyval.npval) = NULL; }
    break;

  case 68:
#line 300 "rmgrm.y"
    { RM_mark();
				  RM_expr( 0, (yyvsp[(1) - (2)].npval) );
				  RM_clear();
				}
    break;

  case 69:
#line 306 "rmgrm.y"
    { RM_mark();
				  RM_expr( 0, (yyvsp[(1) - (2)].npval) );
				  RM_clear();
				}
    break;

  case 70:
#line 312 "rmgrm.y"
    { RM_break( (yyvsp[(2) - (3)].npval) ); (yyval.npval) = NULL; }
    break;

  case 71:
#line 315 "rmgrm.y"
    { RM_expr( 0, (yyvsp[(1) - (2)].npval) );
				  RM_clear();
				}
    break;

  case 72:
#line 320 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;

  case 73:
#line 323 "rmgrm.y"
    { RM_continue( (yyvsp[(2) - (3)].npval) ); (yyval.npval) = NULL; }
    break;

  case 75:
#line 327 "rmgrm.y"
    { RM_endfor(); }
    break;

  case 76:
#line 330 "rmgrm.y"
    { RM_hold( (yyvsp[(2) - (3)].npval) ); }
    break;

  case 77:
#line 332 "rmgrm.y"
    { RM_endif(); }
    break;

  case 78:
#line 334 "rmgrm.y"
    { RM_else(); }
    break;

  case 79:
#line 336 "rmgrm.y"
    { RM_endelse(); }
    break;

  case 80:
#line 339 "rmgrm.y"
    { RM_reject(); (yyval.npval) = NULL; }
    break;

  case 81:
#line 342 "rmgrm.y"
    { RM_release( (yyvsp[(2) - (3)].npval) ); }
    break;

  case 82:
#line 345 "rmgrm.y"
    { RM_while( (yyvsp[(3) - (3)].npval) ); (yyval.npval) = NULL; }
    break;

  case 83:
#line 347 "rmgrm.y"
    { RM_endwhile(); (yyval.npval) = NULL; }
    break;

  case 84:
#line 349 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_INT, &rm_tokval, 0, 0 ); }
    break;

  case 85:
#line 350 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;

  case 86:
#line 353 "rmgrm.y"
    { RM_if( (yyvsp[(3) - (3)].npval) ); }
    break;

  case 87:
#line 355 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;

  case 88:
#line 358 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;

  case 89:
#line 360 "rmgrm.y"
    {  RM_forinit( (yyvsp[(1) - (1)].npval) ); }
    break;

  case 90:
#line 362 "rmgrm.y"
    { RM_fortest( (yyvsp[(4) - (4)].npval) ); }
    break;

  case 91:
#line 364 "rmgrm.y"
    { RM_forincr( (yyvsp[(7) - (7)].npval) ); }
    break;

  case 92:
#line 366 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 93:
#line 367 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 94:
#line 368 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 95:
#line 370 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 96:
#line 371 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 97:
#line 372 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 98:
#line 374 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 99:
#line 375 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 100:
#line 376 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 101:
#line 380 "rmgrm.y"
    { (yyval.npval) = RM_node( (yyvsp[(2) - (3)].ival), 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) );
				  if( rm_context == CTX_PARMS )
					PARM_add( (yyval.npval) );
				  else if( rm_context == CTX_DESCR ||
					rm_context == CTX_SITES )
					SE_addval( (yyval.npval) );
				}
    break;

  case 102:
#line 388 "rmgrm.y"
    { (yyval.npval) = RM_node( (yyvsp[(2) - (3)].ival), 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) );
				  if( rm_context == CTX_PARMS )
					PARM_add( (yyval.npval) );
				  else if( rm_context == CTX_DESCR ||
					rm_context == CTX_SITES )
					SE_addval( (yyval.npval) );
				}
    break;

  case 103:
#line 396 "rmgrm.y"
    { (yyval.ival) = SYM_ASSIGN; }
    break;

  case 104:
#line 398 "rmgrm.y"
    { (yyval.ival) = SYM_MINUS_ASSIGN; }
    break;

  case 105:
#line 400 "rmgrm.y"
    { (yyval.ival) = SYM_PLUS_ASSIGN; }
    break;

  case 106:
#line 402 "rmgrm.y"
    { (yyval.ival) = SYM_PERCENT_ASSIGN; }
    break;

  case 107:
#line 404 "rmgrm.y"
    { (yyval.ival) = SYM_SLASH_ASSIGN; }
    break;

  case 108:
#line 406 "rmgrm.y"
    { (yyval.ival) = SYM_STAR_ASSIGN; }
    break;

  case 109:
#line 408 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 110:
#line 410 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_OR, 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) ); }
    break;

  case 111:
#line 412 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 112:
#line 414 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_AND, 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) ); }
    break;

  case 113:
#line 416 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 114:
#line 417 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 115:
#line 419 "rmgrm.y"
    { (yyval.npval) = RM_node( (yyvsp[(2) - (3)].ival), 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) ); }
    break;

  case 116:
#line 422 "rmgrm.y"
    { (yyval.ival) = SYM_DONT_MATCH; }
    break;

  case 117:
#line 423 "rmgrm.y"
    { (yyval.ival) = SYM_EQUAL; }
    break;

  case 118:
#line 424 "rmgrm.y"
    { (yyval.ival) = SYM_GREATER; }
    break;

  case 119:
#line 426 "rmgrm.y"
    { (yyval.ival) = SYM_GREATER_EQUAL; }
    break;

  case 120:
#line 427 "rmgrm.y"
    { (yyval.ival) = SYM_LESS; }
    break;

  case 121:
#line 429 "rmgrm.y"
    { (yyval.ival) = SYM_LESS_EQUAL; }
    break;

  case 122:
#line 430 "rmgrm.y"
    { (yyval.ival) = SYM_MATCH; }
    break;

  case 123:
#line 431 "rmgrm.y"
    { (yyval.ival) = SYM_NOT_EQUAL; }
    break;

  case 124:
#line 433 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 125:
#line 435 "rmgrm.y"
    { (yyval.npval) = RM_node( (yyvsp[(2) - (3)].ival), 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) ); }
    break;

  case 126:
#line 437 "rmgrm.y"
    { (yyval.ival) = SYM_PLUS; }
    break;

  case 127:
#line 438 "rmgrm.y"
    { (yyval.ival) = SYM_MINUS; }
    break;

  case 128:
#line 440 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 129:
#line 442 "rmgrm.y"
    { (yyval.npval) = RM_node( (yyvsp[(2) - (3)].ival), 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) ); }
    break;

  case 130:
#line 444 "rmgrm.y"
    { (yyval.ival) = SYM_PERCENT; }
    break;

  case 131:
#line 445 "rmgrm.y"
    { (yyval.ival) = SYM_SLASH; }
    break;

  case 132:
#line 446 "rmgrm.y"
    { (yyval.ival) = SYM_STAR; }
    break;

  case 133:
#line 448 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 134:
#line 450 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_NEGATE, 0, 0, (yyvsp[(2) - (2)].npval) ); }
    break;

  case 135:
#line 452 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_NOT, 0, 0, (yyvsp[(2) - (2)].npval) ); }
    break;

  case 136:
#line 453 "rmgrm.y"
    { if( rm_context == CTX_SCORE )
					(yyval.npval) = (yyvsp[(1) - (1)].npval);
				}
    break;

  case 137:
#line 457 "rmgrm.y"
    { if( rm_context == CTX_SCORE )
					(yyval.npval) = (yyvsp[(1) - (1)].npval);
				}
    break;

  case 138:
#line 461 "rmgrm.y"
    { if( rm_context == CTX_SCORE )
					(yyval.npval) = RM_node( SYM_COLON, 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) );
				}
    break;

  case 139:
#line 465 "rmgrm.y"
    { if( rm_context == CTX_SCORE )
					(yyval.npval) = (yyvsp[(1) - (1)].npval);
				}
    break;

  case 140:
#line 469 "rmgrm.y"
    { if( rm_context == CTX_SCORE )
					(yyval.npval) = RM_node( SYM_COLON, 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) );
				}
    break;

  case 141:
#line 473 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 142:
#line 474 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 143:
#line 475 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 144:
#line 477 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(2) - (3)].npval); }
    break;

  case 145:
#line 480 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_CALL, 0, (yyvsp[(1) - (4)].npval), (yyvsp[(3) - (4)].npval) ); }
    break;

  case 146:
#line 482 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 147:
#line 483 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 148:
#line 486 "rmgrm.y"
    { if( rm_context == CTX_DESCR )
					SE_close();
				  else if( rm_context == CTX_SITES )
					POS_close();
				  else if( rm_context == CTX_SCORE )
					(yyval.npval) = RM_node( SYM_KW_STREF, 0, (yyvsp[(1) - (4)].npval), (yyvsp[(3) - (4)].npval) );
				}
    break;

  case 149:
#line 495 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_IX_STREF, 0, (yyvsp[(1) - (4)].npval), (yyvsp[(3) - (4)].npval) ); }
    break;

  case 150:
#line 497 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 151:
#line 498 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 152:
#line 500 "rmgrm.y"
    { (yyval.npval) = RM_node( (yyvsp[(1) - (2)].ival), 0, 0, (yyvsp[(2) - (2)].npval) ); }
    break;

  case 153:
#line 501 "rmgrm.y"
    { (yyval.npval) = RM_node( (yyvsp[(2) - (2)].ival), 0, (yyvsp[(1) - (2)].npval), 0 ); }
    break;

  case 154:
#line 503 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_INT, &rm_tokval, 0, 0 ); }
    break;

  case 155:
#line 504 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_FLOAT, &rm_tokval, 0, 0 ); }
    break;

  case 156:
#line 505 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_DOLLAR,
					&rm_tokval, 0, 0 ); }
    break;

  case 157:
#line 507 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 158:
#line 508 "rmgrm.y"
    { (yyval.npval) = (yyvsp[(1) - (1)].npval); }
    break;

  case 159:
#line 510 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_IDENT,
					&rm_tokval, 0, 0 ); }
    break;

  case 160:
#line 514 "rmgrm.y"
    { (yyval.ival) = SYM_MINUS_MINUS; }
    break;

  case 161:
#line 515 "rmgrm.y"
    { (yyval.ival) = SYM_PLUS_PLUS; }
    break;

  case 162:
#line 517 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_LIST, 0, (yyvsp[(1) - (1)].npval), 0 ); }
    break;

  case 163:
#line 519 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_LIST, 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) ); }
    break;

  case 164:
#line 521 "rmgrm.y"
    { if( rm_context == CTX_SCORE )
					(yyval.npval) = RM_node( SYM_LIST, 0, (yyvsp[(1) - (1)].npval), 0 );
				}
    break;

  case 165:
#line 525 "rmgrm.y"
    { if( rm_context == CTX_SCORE )
					(yyval.npval) = RM_node( SYM_LIST, 0, (yyvsp[(1) - (3)].npval), (yyvsp[(3) - (3)].npval) );
				}
    break;

  case 166:
#line 529 "rmgrm.y"
    { PR_open(); }
    break;

  case 167:
#line 531 "rmgrm.y"
    { (yyval.npval) = PR_close(); }
    break;

  case 168:
#line 533 "rmgrm.y"
    { PR_add( (yyval.npval) ); }
    break;

  case 169:
#line 535 "rmgrm.y"
    { PR_add( (yyvsp[(1) - (3)].npval) ) ; }
    break;

  case 170:
#line 537 "rmgrm.y"
    { (yyval.npval) = RM_node( SYM_STRING,
					&rm_tokval, 0, 0 ); }
    break;

  case 171:
#line 540 "rmgrm.y"
    { (yyval.npval) = NULL; }
    break;


/* Line 1267 of yacc.c.  */
#line 2541 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 542 "rmgrm.y"


#include "lex.yy.c"

int	yyerror( msg )
char	msg[];
{

	fprintf( stderr, "yyerror: %s\n", msg );
	return( 0 );
}

