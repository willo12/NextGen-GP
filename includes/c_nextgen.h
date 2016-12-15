#ifndef c_nextgen_h__
#define c_nextgen_h__
 

extern void c_map_f(char *treestring, double *A, double *B, int shape_i, int shape_j, double S0_start, double S1_start, double S0_end, double S1_end, double *forcing, int forcing_dim);
extern void c_single_tree(double* ffs,double* result,double* obs, long*  I_indices, char treestr[], int ffs0, int ffs1, int result0,int obs0,int obs1, double *model, int model0, int ts_factor, int startscore_i, double S_init);
extern void c_nextgen(double* ffs,double* result, double* old_score_vals, double* score_vals,double* obs, long*  I_indices, char treefile[], char treeoutfile[], int ffs0, int ffs1, int result0,int my_number,int qsubs,int obs0,int obs1,int runlen,int popsize, double *model, int model0, int ts_factor, int startscore_i, double S_init);

extern int test_llvm(void);

#endif
