/**********************************************/
/* BEGIN interface for sliding-window MFE     */
/* prediction                                 */
/**********************************************/

%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") mfe_window;
%feature("kwargs") mfe_window;
#endif

  float mfe_window(FILE *nullfile = NULL){

    return vrna_mfe_window($self, nullfile);
  }

  /* ONLY possible if VRNA_WITH_SVM is set
  float mfe_window_zscore(double min_z,FILE *file=NULL){

    return vrna_mfe_window_zscore($self,min_z,file);
  }*/
  
}

/* tell swig that these functions return objects that require memory management */

%include <ViennaRNA/mfe.h>


%include <ViennaRNA/Lfold.h>
