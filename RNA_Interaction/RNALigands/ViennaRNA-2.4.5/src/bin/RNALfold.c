/*
 *                Ineractive Access to folding Routines
 *
 *                c Ivo L Hofacker
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/mfe.h"
#include "ViennaRNA/Lfold.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/commands.h"
#include "ViennaRNA/constraints_SHAPE.h"
#include "RNALfold_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"

#include "ViennaRNA/color_output.inc"

typedef struct {
  FILE  *output;
  int   dangle_model;
} hit_data;


#ifdef VRNA_WITH_SVM
PRIVATE void
default_callback_z(int        start,
                   int        end,
                   const char *structure,
                   float      en,
                   float      zscore,
                   void       *data);


#endif

PRIVATE void
default_callback(int        start,
                 int        end,
                 const char *structure,
                 float      en,
                 void       *data);


int
main(int  argc,
     char *argv[])
{
  FILE                        *input, *output;
  struct  RNALfold_args_info  args_info;
  char                        *ParamFile, *ns_bases, *rec_sequence, *rec_id, **rec_rest,
                              *command_file, *orig_sequence, *infile, *outfile, *filename_delim,
                              *shape_file, *shape_method, *shape_conversion;
  unsigned int                rec_type, read_opt;
  int                         length, istty, noconv, maxdist, zsc, tofile, filename_full,
                              with_shapes, verbose;
  double                      min_en, min_z;
  vrna_md_t                   md;
  vrna_cmd_t                  *commands;
  dataset_id                  id_control;

  ParamFile     = ns_bases = NULL;
  do_backtrack  = 1;
  noconv        = 0;
  dangles       = 2;
  maxdist       = 150;
  zsc           = 0;
  min_z         = -2.0;
  gquad         = 0;
  rec_type      = read_opt = 0;
  rec_id        = rec_sequence = orig_sequence = NULL;
  rec_rest      = NULL;
  outfile       = NULL;
  infile        = NULL;
  input         = NULL;
  output        = NULL;
  tofile        = 0;
  filename_full = 0;
  command_file  = NULL;
  commands      = NULL;

  /* apply default model details */
  vrna_md_set_default(&md);

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNALfold_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* parse options for ID manipulation */
  ggo_get_id_control(args_info, id_control, "Sequence", "sequence", "_", 4, 1);

  /* temperature */
  if (args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;

  /* do not take special tetra loop energies into account */
  if (args_info.noTetra_given)
    md.special_hp = tetra_loop = 0;

  /* set dangle model */
  if (args_info.dangles_given) {
    if ((args_info.dangles_arg < 0) || (args_info.dangles_arg > 3))
      vrna_message_warning(
        "required dangle model not implemented, falling back to default dangles=2");
    else
      md.dangles = dangles = args_info.dangles_arg;
  }

  /* do not allow weak pairs */
  if (args_info.noLP_given)
    md.noLP = noLonelyPairs = 1;

  /* do not allow wobble pairs (GU) */
  if (args_info.noGU_given)
    md.noGU = noGU = 1;

  /* do not allow weak closing pairs (AU,GU) */
  if (args_info.noClosingGU_given)
    md.noGUclosure = no_closingGU = 1;

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (args_info.noconv_given)
    noconv = 1;

  /* set energy model */
  if (args_info.energyModel_given)
    md.energy_set = energy_set = args_info.energyModel_arg;

  /* take another energy parameter set */
  if (args_info.paramFile_given)
    ParamFile = strdup(args_info.paramFile_arg);

  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if (args_info.nsp_given)
    ns_bases = strdup(args_info.nsp_arg);

  /* set the maximum base pair span */
  if (args_info.span_given)
    maxdist = args_info.span_arg;

  if (args_info.zscore_given) {
#ifdef VRNA_WITH_SVM
    zsc = 1;
    if (args_info.zscore_arg != -2)
      min_z = args_info.zscore_arg;

#else
    vrna_message_error("\'z\' option is available only if compiled with SVM support!");
#endif
  }

  /* gquadruplex support */
  if (args_info.gquad_given)
    md.gquad = gquad = 1;

  if (args_info.verbose_given)
    verbose = 1;

  /* SHAPE reactivity data */
  ggo_get_SHAPE(args_info, with_shapes, shape_file, shape_method, shape_conversion);

  if (args_info.outfile_given) {
    tofile = 1;
    if (args_info.outfile_arg)
      outfile = strdup(args_info.outfile_arg);
  }

  if (args_info.infile_given)
    infile = strdup(args_info.infile_arg);

  /* filename sanitize delimiter */
  if (args_info.filename_delim_given)
    filename_delim = strdup(args_info.filename_delim_arg);
  else if (get_id_delim(id_control))
    filename_delim = strdup(get_id_delim(id_control));

  if ((filename_delim) && isspace(*filename_delim)) {
    free(filename_delim);
    filename_delim = NULL;
  }

  /* full filename from FASTA header support */
  if (args_info.filename_full_given)
    filename_full = 1;

  if (args_info.commands_given)
    command_file = strdup(args_info.commands_arg);

  /* check for errorneous parameter options */
  if (maxdist <= 0) {
    RNALfold_cmdline_parser_print_help();
    exit(EXIT_FAILURE);
  }

  /* free allocated memory of command line data structure */
  RNALfold_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */

  md.max_bp_span = md.window_size = maxdist;

  if (infile) {
    input = fopen((const char *)infile, "r");
    if (!input)
      vrna_message_error("Could not read input file");
  } else {
    input = stdin;
  }

  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (command_file != NULL)
    commands = vrna_file_commands_read(command_file, VRNA_CMD_PARSE_HC | VRNA_CMD_PARSE_SC);

  if (ns_bases != NULL)
    vrna_md_set_nonstandards(&md, ns_bases);

  istty     = (!infile) && isatty(fileno(stdout)) && isatty(fileno(stdin));
  read_opt  |= VRNA_INPUT_NO_REST;
  if (istty) {
    vrna_message_input_seq_simple();
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  }

  /*
   #############################################
   # main loop: continue until end of file
   #############################################
   */
  while (
    !((rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, input, read_opt))
      & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))) {
    /*
     ########################################################
     # init everything according to the data we've read
     ########################################################
     */
    char  *SEQ_ID       = NULL;
    char  *v_file_name  = NULL;
    char  *tmp_string   = NULL;
    /*
     ########################################################
     # init everything according to the data we've read
     ########################################################
     */
    if (rec_id) /* remove '>' from FASTA header */
      rec_id = memmove(rec_id, rec_id + 1, strlen(rec_id));

    /* construct the sequence ID */
    set_next_id(&rec_id, id_control);
    SEQ_ID = fileprefix_from_id(rec_id, id_control, filename_full);

    if (tofile) {
      /* prepare the file name */
      if (outfile)
        v_file_name = vrna_strdup_printf("%s", outfile);
      else
        v_file_name = (SEQ_ID) ?
                      vrna_strdup_printf("%s.lfold", SEQ_ID) :
                      vrna_strdup_printf("RNALfold_output.lfold");

      tmp_string = vrna_filename_sanitize(v_file_name, filename_delim);
      free(v_file_name);
      v_file_name = tmp_string;

      if (infile && !strcmp(infile, v_file_name))
        vrna_message_error("Input and output file names are identical");

      output = fopen((const char *)v_file_name, "a");
      if (!output)
        vrna_message_error("Failed to open file for writing");
    } else {
      output = stdout;
    }

    if (!istty)
      print_fasta_header(output, rec_id);

    length = (int)strlen(rec_sequence);

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if (!noconv)
      vrna_seq_toRNA(rec_sequence);

    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(rec_sequence);

    if (!tofile && istty)
      vrna_message_info(output, "length = %d", length);

    /*
     ########################################################
     # done with 'stdin' handling
     # begin actual computations
     ########################################################
     */

    vrna_fold_compound_t *vc = vrna_fold_compound((const char *)rec_sequence,
                                                  &md,
                                                  VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

    if (commands)
      vrna_commands_apply(vc, commands, VRNA_CMD_PARSE_HC | VRNA_CMD_PARSE_SC);

    if (with_shapes) {
      vrna_constraints_add_SHAPE(vc,
                                 shape_file,
                                 shape_method,
                                 shape_conversion,
                                 verbose,
                                 VRNA_OPTION_WINDOW);
    }

    hit_data data;
    data.output       = output;
    data.dangle_model = md.dangles;

#ifdef VRNA_WITH_SVM
    min_en =
      (zsc) ? vrna_mfe_window_zscore_cb(vc, min_z, &default_callback_z,
                                        (void *)&data) : vrna_mfe_window_cb(vc, &default_callback,
                                                                            (void *)&data);
#else
    min_en = vrna_mfe_window_cb(vc, &default_callback, (void *)&data);
#endif
    fprintf(output, "%s\n", orig_sequence);

    char *msg = NULL;
    if (!tofile && istty)
      msg = vrna_strdup_printf(" minimum free energy = %6.2f kcal/mol", min_en);
    else
      msg = vrna_strdup_printf(" (%6.2f)", min_en);

    print_structure(output, NULL, msg);
    free(msg);

    if (output)
      (void)fflush(output);

    if (tofile && output) {
      fclose(output);
      output = NULL;
    }

    /* clean up */
    vrna_fold_compound_free(vc);
    free(rec_id);
    free(SEQ_ID);
    free(rec_sequence);
    free(orig_sequence);
    rec_id    = rec_sequence = orig_sequence = NULL;
    rec_rest  = NULL;

    free(v_file_name);

    if (with_shapes)
      break;

    /* print user help for the next round if we get input from tty */

    if (istty)
      vrna_message_input_seq_simple();
  }

  if (infile && input)
    fclose(input);

  free(filename_delim);
  free(command_file);
  vrna_commands_free(commands);

  free_id_data(id_control);

  return EXIT_SUCCESS;
}


PRIVATE void
default_callback(int        start,
                 int        end,
                 const char *structure,
                 float      en,
                 void       *data)
{
  FILE  *output       = ((hit_data *)data)->output;
  int   dangle_model  = ((hit_data *)data)->dangle_model;
  char  *struct_d2    = NULL;
  char  *msg          = NULL;

  if ((dangle_model == 2) && (start > 1)) {
    msg       = vrna_strdup_printf(" (%6.2f) %4d", en, start - 1);
    struct_d2 = vrna_strdup_printf(".%s", structure);
    print_structure(output, struct_d2, msg);
    free(struct_d2);
  } else {
    msg = vrna_strdup_printf(" (%6.2f) %4d", en, start);
    print_structure(output, structure, msg);
  }

  free(msg);
}


#ifdef VRNA_WITH_SVM
PRIVATE void
default_callback_z(int        start,
                   int        end,
                   const char *structure,
                   float      en,
                   float      zscore,
                   void       *data)
{
  FILE  *output       = ((hit_data *)data)->output;
  int   dangle_model  = ((hit_data *)data)->dangle_model;
  char  *struct_d2    = NULL;
  char  *msg          = NULL;

  if ((dangle_model == 2) && (start > 1)) {
    msg       = vrna_strdup_printf(" (%6.2f) %4d z= %.3f", en, start - 1, zscore);
    struct_d2 = vrna_strdup_printf(".%s", structure);
    print_structure(output, struct_d2, msg);
    free(struct_d2);
  } else {
    msg = vrna_strdup_printf(" (%6.2f) %4d z= %.3f", en, start, zscore);
    print_structure(output, structure, msg);
  }

  free(msg);
}


#endif
