	case SYM_PARMS :
		fprintf( fp, "SYM_PARMS\n" );
		break;
	case SYM_DESCR :
		fprintf( fp, "SYM_DESCR\n" );
		break;
	case SYM_SITES :
		fprintf( fp, "SYM_SITES\n" );
		break;
	case SYM_SCORE :
		fprintf( fp, "SYM_SCORE\n" );
		break;
	case SYM_SE :
		fprintf( fp, "SYM_SE\n" );
		break;
	case SYM_CTX :
		fprintf( fp, "SYM_CTX\n" );
		break;
	case SYM_SS :
		fprintf( fp, "SYM_SS\n" );
		break;
	case SYM_H5 :
		fprintf( fp, "SYM_H5\n" );
		break;
	case SYM_H3 :
		fprintf( fp, "SYM_H3\n" );
		break;
	case SYM_P5 :
		fprintf( fp, "SYM_P5\n" );
		break;
	case SYM_P3 :
		fprintf( fp, "SYM_P3\n" );
		break;
	case SYM_T1 :
		fprintf( fp, "SYM_T1\n" );
		break;
	case SYM_T2 :
		fprintf( fp, "SYM_T2\n" );
		break;
	case SYM_T3 :
		fprintf( fp, "SYM_T3\n" );
		break;
	case SYM_Q1 :
		fprintf( fp, "SYM_Q1\n" );
		break;
	case SYM_Q2 :
		fprintf( fp, "SYM_Q2\n" );
		break;
	case SYM_Q3 :
		fprintf( fp, "SYM_Q3\n" );
		break;
	case SYM_Q4 :
		fprintf( fp, "SYM_Q4\n" );
		break;
	case SYM_ACCEPT :
		fprintf( fp, "SYM_ACCEPT\n" );
		break;
	case SYM_BEGIN :
		fprintf( fp, "SYM_BEGIN\n" );
		break;
	case SYM_BREAK :
		fprintf( fp, "SYM_BREAK\n" );
		break;
	case SYM_CONTINUE :
		fprintf( fp, "SYM_CONTINUE\n" );
		break;
	case SYM_ELSE :
		fprintf( fp, "SYM_ELSE\n" );
		break;
	case SYM_END :
		fprintf( fp, "SYM_END\n" );
		break;
	case SYM_FOR :
		fprintf( fp, "SYM_FOR\n" );
		break;
	case SYM_HOLD :
		fprintf( fp, "SYM_HOLD\n" );
		break;
	case SYM_IF :
		fprintf( fp, "SYM_IF\n" );
		break;
	case SYM_IN :
		fprintf( fp, "SYM_IN\n" );
		break;
	case SYM_REJECT :
		fprintf( fp, "SYM_REJECT\n" );
		break;
	case SYM_RELEASE :
		fprintf( fp, "SYM_RELEASE\n" );
		break;
	case SYM_WHILE :
		fprintf( fp, "SYM_WHILE\n" );
		break;
	case SYM_IDENT :
		fprintf( fp, "SYM_IDENT = \"%s\"\n", (char *)np->n_val.v_value.v_pval );
		break;
	case SYM_INT :
		fprintf( fp, "SYM_INT = %d\n", np->n_val.v_value.v_ival );
		break;
	case SYM_FLOAT :
		fprintf( fp, "SYM_FLOAT = %lg\n", np->n_val.v_value.v_dval );
		break;
	case SYM_STRING :
		fprintf( fp, "SYM_STRING = \"%s\"\n", (char *)np->n_val.v_value.v_pval );
		break;
	case SYM_PAIRSET :
		fprintf( fp, "SYM_PAIRSET = RM_dump_pairset( fp, np->n_val.v_value.v_pval );\n" );
		break;
	case SYM_AND :
		fprintf( fp, "SYM_AND\n" );
		break;
	case SYM_ASSIGN :
		fprintf( fp, "SYM_ASSIGN\n" );
		break;
	case SYM_DOLLAR :
		fprintf( fp, "SYM_DOLLAR\n" );
		break;
	case SYM_DONT_MATCH :
		fprintf( fp, "SYM_DONT_MATCH\n" );
		break;
	case SYM_EQUAL :
		fprintf( fp, "SYM_EQUAL\n" );
		break;
	case SYM_GREATER :
		fprintf( fp, "SYM_GREATER\n" );
		break;
	case SYM_GREATER_EQUAL :
		fprintf( fp, "SYM_GREATER_EQUAL\n" );
		break;
	case SYM_LESS :
		fprintf( fp, "SYM_LESS\n" );
		break;
	case SYM_LESS_EQUAL :
		fprintf( fp, "SYM_LESS_EQUAL\n" );
		break;
	case SYM_MATCH :
		fprintf( fp, "SYM_MATCH\n" );
		break;
	case SYM_MINUS :
		fprintf( fp, "SYM_MINUS\n" );
		break;
	case SYM_MINUS_ASSIGN :
		fprintf( fp, "SYM_MINUS_ASSIGN\n" );
		break;
	case SYM_MINUS_MINUS :
		fprintf( fp, "SYM_MINUS_MINUS\n" );
		break;
	case SYM_NEGATE :
		fprintf( fp, "SYM_NEGATE\n" );
		break;
	case SYM_NOT :
		fprintf( fp, "SYM_NOT\n" );
		break;
	case SYM_NOT_EQUAL :
		fprintf( fp, "SYM_NOT_EQUAL\n" );
		break;
	case SYM_OR :
		fprintf( fp, "SYM_OR\n" );
		break;
	case SYM_PERCENT :
		fprintf( fp, "SYM_PERCENT\n" );
		break;
	case SYM_PERCENT_ASSIGN :
		fprintf( fp, "SYM_PERCENT_ASSIGN\n" );
		break;
	case SYM_PLUS :
		fprintf( fp, "SYM_PLUS\n" );
		break;
	case SYM_PLUS_ASSIGN :
		fprintf( fp, "SYM_PLUS_ASSIGN\n" );
		break;
	case SYM_PLUS_PLUS :
		fprintf( fp, "SYM_PLUS_PLUS\n" );
		break;
	case SYM_STAR :
		fprintf( fp, "SYM_STAR\n" );
		break;
	case SYM_STAR_ASSIGN :
		fprintf( fp, "SYM_STAR_ASSIGN\n" );
		break;
	case SYM_SLASH :
		fprintf( fp, "SYM_SLASH\n" );
		break;
	case SYM_SLASH_ASSIGN :
		fprintf( fp, "SYM_SLASH_ASSIGN\n" );
		break;
	case SYM_LPAREN :
		fprintf( fp, "SYM_LPAREN\n" );
		break;
	case SYM_RPAREN :
		fprintf( fp, "SYM_RPAREN\n" );
		break;
	case SYM_LBRACK :
		fprintf( fp, "SYM_LBRACK\n" );
		break;
	case SYM_RBRACK :
		fprintf( fp, "SYM_RBRACK\n" );
		break;
	case SYM_LCURLY :
		fprintf( fp, "SYM_LCURLY\n" );
		break;
	case SYM_RCURLY :
		fprintf( fp, "SYM_RCURLY\n" );
		break;
	case SYM_COLON :
		fprintf( fp, "SYM_COLON\n" );
		break;
	case SYM_COMMA :
		fprintf( fp, "SYM_COMMA\n" );
		break;
	case SYM_SEMICOLON :
		fprintf( fp, "SYM_SEMICOLON\n" );
		break;
	case SYM_CALL :
		fprintf( fp, "SYM_CALL\n" );
		break;
	case SYM_LIST :
		fprintf( fp, "SYM_LIST\n" );
		break;
	case SYM_KW_STREF :
		fprintf( fp, "SYM_KW_STREF\n" );
		break;
	case SYM_IX_STREF :
		fprintf( fp, "SYM_IX_STREF\n" );
		break;
	case SYM_ERROR :
		fprintf( fp, "SYM_ERROR\n" );
		break;
	default :
		fprintf( fp, "RM_dumpnode: Unknown symbol %d\n", np->n_sym );
		break;
