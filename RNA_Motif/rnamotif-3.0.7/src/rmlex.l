%{

#include <string.h>
#include "rnamot.h"
#include "y.tab.h"

double	atof();
char	*RM_str2seq();

extern	int	rm_context;
extern	char	*rm_wdfname;
extern	int	rm_lineno;
extern	int	rm_emsg_lineno;
extern	VALUE_T	rm_tokval;

static	char	*sp;
static	POS_T	*posp;
static	char	emsg[ 256 ];

typedef	struct	fname_t	{
	struct	fname_t	*f_left;
	struct	fname_t	*f_right;
	char	*f_name;
} FNAME_T;

static	FNAME_T	*ftree;

static	void	setfileinfo( char [] );
static	FNAME_T	*enter_fname( FNAME_T **, char [] );

%}
%%

^#\ line.*	setfileinfo( yytext );
#.*		;

\n		;
\r		;

[\ \t\f]+	;

parms		return( SYM_PARMS );
descr		return( SYM_DESCR );
sites		return( SYM_SITES );
score		return( SYM_SCORE );

se		return( SYM_SE );
ctx		return( SYM_CTX );
ss		return( SYM_SS );
h5		return( SYM_H5 );
h3		return( SYM_H3 );
p5		return( SYM_P5 );
p3		return( SYM_P3 );
t1		return( SYM_T1 );
t2		return( SYM_T2 );
t3		return( SYM_T3 );
q1		return( SYM_Q1 );
q2		return( SYM_Q2 );
q3		return( SYM_Q3 );
q4		return( SYM_Q4 );

ACCEPT		return( SYM_ACCEPT );
BEGIN		return( SYM_BEGIN );
END		return( SYM_END );
HOLD		return( SYM_HOLD );
REJECT		return( SYM_REJECT );
RELEASE		return( SYM_RELEASE );
break		return( SYM_BREAK );
continue	return( SYM_CONTINUE );
else		return( SYM_ELSE );
for		return( SYM_FOR );
if		return( SYM_IF );
in		return( SYM_IN );
while		return( SYM_WHILE );

[a-zA-Z][_a-zA-Z0-9]*	{
		  sp = ( char * )malloc(( yyleng + 1 ) *sizeof( char ));
		  if( sp == NULL ){
			RM_errormsg( 1, 
				"lex.yy.c: can't alloc sp for SYM_IDENT." );
		  }
		  strcpy( sp, yytext );
		  rm_tokval.v_type = T_STRING;
		  rm_tokval.v_value.v_pval = sp;
		  return( SYM_IDENT );
		}
[0-9]+		{ rm_tokval.v_type = T_INT;
		  rm_tokval.v_value.v_ival = atoi( yytext );
		  return( SYM_INT );
		}
([0-9]+[eE][\+\-]?[0-9]+)|((([0-9]+\.[0-9]*)|(\.[0-9]+))([eE][\+\-]?[0-9]+)?) {
		  rm_tokval.v_type = T_FLOAT;
		  rm_tokval.v_value.v_dval = atof( yytext );
		  return( SYM_FLOAT );
		}
\"[^"\n]*	{ if( yytext[ yyleng - 1 ] == '\\' )
			yymore();
		  else{
			yytext[ yyleng++ ] = input();
			yytext[ yyleng ] = '\0';
			sp = ( char * )malloc( ( yyleng - 1 )*sizeof( char ) );
			if( sp == NULL ){
				RM_errormsg( 1,
				"lex.yy.c: can't alloc sp for SYM_STRING." );
			}
			strncpy( sp, &yytext[1], yyleng - 2 );
			sp[ yyleng - 2 ] = '\0';
			rm_tokval.v_type = T_STRING;
/*
			rm_tokval.v_value.v_pval = sp;
*/
			rm_tokval.v_value.v_pval = RM_str2seq( sp );
			return( SYM_STRING );
		  }
		}
\'[^'\n]*	{ if( yytext[ yyleng - 1 ] == '\\' )
			yymore();
		  else{
			yytext[ yyleng++ ] = input();
			yytext[ yyleng ] = '\0';
			sp = ( char * )malloc( ( yyleng - 1 )*sizeof( char ) );
			if( sp == NULL ){
				RM_errormsg( 1,
				"lex.yy.c: can't alloc sp for SYM_STRING." );
			}
			strncpy( sp, &yytext[1], yyleng - 2 );
			sp[ yyleng - 2 ] = '\0';
			rm_tokval.v_type = T_STRING;
			rm_tokval.v_value.v_pval = sp;
			return( SYM_STRING );
		  }
		}

&&		return( SYM_AND );
=		return( SYM_ASSIGN );
\$		{ posp = ( POS_T * )malloc( sizeof( POS_T ) );
		  if( posp == NULL ){
			RM_errormsg( 1,
				"lex.yy.c: can't alloc posp for SYM_DOLLAR." );
		  }
		  posp->p_type = SYM_DOLLAR;
		  posp->p_lineno = rm_lineno;
		  posp->p_tag = NULL;
		  posp->p_addr.a_l2r = 0;
		  posp->p_addr.a_offset = 0;
		  rm_tokval.v_type = T_POS;
		  rm_tokval.v_value.v_pval = posp;
		  return( SYM_DOLLAR );
		}
!~		return( SYM_DONT_MATCH );
==		return( SYM_EQUAL );
\>		return( SYM_GREATER );
\>=		return( SYM_GREATER_EQUAL );
\<		return( SYM_LESS );
\<=		return( SYM_LESS_EQUAL );
=~		return( SYM_MATCH );
\-		return( SYM_MINUS );
\-=		return( SYM_MINUS_ASSIGN );
\-\-		return( SYM_MINUS_MINUS );
!		return( SYM_NOT );
!=		return( SYM_NOT_EQUAL );
\|\|		return( SYM_OR );
%		return( SYM_PERCENT );
%=		return( SYM_PERCENT_ASSIGN );
\+		return( SYM_PLUS );
\+=		return( SYM_PLUS_ASSIGN );
\+\+		return( SYM_PLUS_PLUS );
\*		return( SYM_STAR );
\*=		return( SYM_STAR_ASSIGN );
\/		return( SYM_SLASH );
\/=		return( SYM_SLASH_ASSIGN );

\(		return( SYM_LPAREN );
\)		return( SYM_RPAREN );
\[		return( SYM_LBRACK );
\]		return( SYM_RBRACK );
\{		return( SYM_LCURLY );
\}		return( SYM_RCURLY );
\,		return( SYM_COMMA );
:		return( SYM_COLON );
;		return( SYM_SEMICOLON );

.		return( SYM_ERROR );

%%

int	yywrap()
{

	return( 1 );
}

static	void	setfileinfo( char data[] )
{
	char	*dp, *qp;
	int	lnum;
	char	work[ 256 ];
	FNAME_T	*fnp;

	dp = &data[ 6 ];	/* skip # line	*/
	if( !isspace( *dp ) )
		return;
	dp += strspn( dp, " \t" );
	if( !isdigit( *dp ) )
		return;
	for( lnum = *dp++ - '0'; isdigit( *dp ); dp++ )
		lnum = 10 * lnum + *dp - '0';
	if( !isspace( *dp ) )
		return;
	for( ; isspace( *dp ); dp++ )
		;
	if( *dp != '\'' )
		return;
	dp++;
	if( qp = strchr( dp, '\'' ) ){
		rm_lineno = lnum;
		strncpy( work, dp, qp - dp );
		work[ qp - dp ] = '\0';
		fnp = enter_fname( &ftree, work );
		rm_wdfname = fnp->f_name;
	}
}

static	FNAME_T	*enter_fname( FNAME_T **root, char fname[] )
{
	char	*sp;
	FNAME_T	*fnp;
	int	cv;

	if( *root == NULL ){
		fnp = ( FNAME_T * )malloc( sizeof( FNAME_T ) );
		if( fnp == NULL ){
			sprintf( emsg,
			"enter_fname: can't allocate space for fnp for '%s'.",
				fname );
			rm_emsg_lineno = rm_lineno;
			RM_errormsg( 1, emsg );
		}
		fnp->f_left = NULL;
		fnp->f_right = NULL;
		sp = ( char * )malloc( strlen( fname ) + 1 );
		if( sp == NULL ){
			sprintf( emsg,
		"enter_fname: can't allocate space for f_name for '%s'.",
				fname );
			rm_emsg_lineno = rm_lineno;
			RM_errormsg( 1, emsg );
		}
		strcpy( sp, fname );
		fnp->f_name = sp;
		*root = fnp;
	}else if( ( cv = strcmp( ( *root )->f_name, fname ) ) < 0 )
		fnp = enter_fname( &( *root )->f_right, fname );
	else if( cv > 0 )
		fnp = enter_fname( &( *root )->f_left, fname );
	else
		fnp = *root;
	return( fnp );
}
