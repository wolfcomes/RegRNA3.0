
/*=============================================================================
 outputs.c                          A.Lambert le 21/12/01 revu le 12/05/03

 fonctions concernant la sortie de resultats lors d'une operation de recherche
 d'un masque (recherche a 1 niveau ou plus), ou en fin d'operation.
 
 Position affichee des motif detectes:
 Brin direct: elle indique le debut du motif depuis l'origine de la sequence.
 Brin complementaire: elle indique la trace du motif sur le brin direct (anno-
 tation "genbank").
 Optionnellement la E-value des detections est aussi affichee.

 cc -O2 -Wall -c outputs.c -I../include ;
 ar -rs ../lib/librnaIV.a outputs.o ;

=============================================================================*/

#include "rnaIV.h"

void PrintOutput(Context *ctxt, int detect, Config *cfg, double score);
void fPrintfOutput(Context *ctxt, int detect);
void fPrintfData(Context *ctxt, Data *data);

/*=============================================================================
 PrintOutput(): Affiche, suivant le contexte pointe par 'ctxt', le resultat
                d'une detection du masque associe en position 'detect' reperee
                depuis l'origine de la sequence, et dont la configuration est
                pointee par 'cfg' et de score 'score'.
                En output le debut des sequences est repere par 1 (et non pas 0
                comme en "interne"): 'detect' est decalle de 1.
                ATTENTION: 'cfg' doit etre une des configurations possibles de
                'ctxt->mask', aucun controle n'est effectue.
=============================================================================*/

void PrintOutput(Context *ctxt, int detect, Config *cfg, double score)
{
  int pos, end, len;
  
  if (ctxt->hist == ON)  AddToDetectsHisto(score);      /* actualise l'histo */

  if (ctxt->style == MUTE) return;

  len = cfg->len;

  if (ctxt->style == LONG && ctxt->detects == 1)  /* une seule fois par seq. */
  {
      fprintf(ctxt->outfile, "%s\n", ctxt->seq->comment);
  }

  pos = detect + 1;                              /* l'origine passe de 0 a 1 */
  end = pos + len - 1;

  if (ctxt->dir == FORWARD)
  {                                                           /* brin direct */
      fprintf(ctxt->outfile, "FW %3d %7d..%-7d  %.2f",
              ctxt->detects, pos, end, score);
      if (ctxt->eval == ON) {
          fprintf(ctxt->outfile, "  %.2e", Evalue(score));      /* + E-value */
      }
      fprintf(ctxt->outfile, "\n");                           /* terminaison */
  }
  else
  {                                                    /* reverse complement */
      pos = ctxt->seq->datalen - pos + 1;     /* inversion de 'pos' et 'end' */
      end = ctxt->seq->datalen - end + 1;

      fprintf(ctxt->outfile, "RC %3d %7d..%-7d  %.2f",
              ctxt->detects, end, pos, score);
      if (ctxt->eval == ON) {
          fprintf(ctxt->outfile, "  %.2e", Evalue(score));      /* + E-value */
      }
      fprintf(ctxt->outfile, "\n");                           /* terminaison */
  }

  if (ctxt->style == LONG)
  {
      RecordMaskSeq(ctxt->mask, ctxt->seq->data + detect, cfg);
      fprintf(ctxt->outfile, "%s\n", ctxt->mask->str);
  }

  return;
}
/*=============================================================================
 fPrintfOutput(): Interface de 'PrintOutput' adaptee a une sortie de resultat
                  en cours d'exploration: 'cfgindex' est enregistre dans la
                  structure 'ctxt->mask' et 'score dans 'ctxt->mask->score'.
=============================================================================*/

void fPrintfOutput(Context *ctxt, int detect)
{
  Config  *cfg = ctxt->mask->cfg + ctxt->mask->cfgindex;

  PrintOutput(ctxt, detect, cfg, ctxt->mask->score);

  return;
}
/*=============================================================================
 fPrintfData(): Interface de 'PrintOutput' adaptee a une sortie de resultats
                en fin de recherche: les donnees ont ete enregistrees dans une 
                structure 'Data' pointee par 'data'.
=============================================================================*/

void fPrintfData(Context *ctxt, Data *data)
{

  if (data->cfg == NULL)
  {
      Config  *cfg = ctxt->mask->cfg + data->cfgindex;
      PrintOutput(ctxt, data->offset, cfg, data->score);
  }
  else
      PrintOutput(ctxt, data->offset, data->cfg, data->score);

  return;
}
/*===========================================================================*/
