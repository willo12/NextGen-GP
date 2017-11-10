#ifndef BASICOPS_H
#define BASICOPS_H

#ifdef INTSTATES
extern int reg[MAXPAR];
#else
extern double reg[MAXPAR];
#endif

extern int arg_table[OPTABLE];

extern char* instrsetnodes; /* for printing */
extern char* instrset; /* used at runtime */
extern char* instrset_rare; /* used at runtime */
extern char* instrsetswap; /* instr that can swap */
extern char* instrsetzero; /* function subleaf nodes */
/*extern char* instrsetzero;  */

extern char* instrsetstem;
extern char* instrsetsuperleaf;
extern char* instrsetleaf;

extern char const_op_char;
extern char par_op_char;


double nopFun(Node *tree);
double constFun(Node *tree);
double parFun(Node *tree);

double parFun0(Node *tree);
double parFun1(Node *tree);
double parFun2(Node *tree);
double parFun3(Node *tree);
double parFun4(Node *tree);


double addFun(Node *tree);
double subFun(Node *tree);
double mulFun(Node *tree);
double divFun(Node *tree);
double isgreaterFun(Node *tree);
double iseqFun(Node *tree);
double ifFun(Node *tree);
double sqrtFun(Node *tree);
double tanhFun(Node *tree);
double markerFun(Node *tree);

double stemFun(Node *tree);
double funFun(Node *tree);

int nopFunInt(Node *tree);
int constFunInt(Node *tree);
int parFunInt(Node *tree);

int addFunInt(Node *tree);
int subFunInt(Node *tree);
int mulFunInt(Node *tree);
int divFunInt(Node *tree);
int isgreaterFunInt(Node *tree);
int iseqFunInt(Node *tree);
int ifFunInt(Node *tree);
int markerFunInt(Node *tree);
int stemFunInt(Node *tree);


Node *init_node(char human_op);
Node *make_par_node(char human_op, int i);
Node *make_par_fun_node(char human_op, int i, double a);
Node *make_const_node_int(char human_op, int value);
Node *make_const_node(char human_op, double value);
Node *make_zero_fun_node(void);

double randomdouble(void);
double randomdouble_slope(double avg, double amp);

Node *makerandomparnode(int max_par);
Node *makerandomparfunnode(int max_par);
Node *makerandomconstnode(void);
Node *makerandomconstnodeint(void);
Node *makerandomterminal_basic(int maxdepth,double fpr,double ppr, int max_par, Node* (*f_const)(void));

void free_param_node(Node *t);
void free_const_node(Node *t);
void free_param_fun_node(Node *t);

// COPY FUNCTIONS
Node *copy_const_node(Node *tree);
Node *copy_par_node(Node *tree);
Node *copy_par_fun_node(Node *tree);
Node *copy_non_leaf(Node *tree);
Node *copy_node(Node *tree);
Node *copy_node_nojump(Node *tree);
Node *make_run_node(Node *tree);

// compare functions
int compare_const_nodes(char op, Node *t1, Node *t2);
int compare_par_nodes(char op, Node *t1, Node *t2);
int compare_par_fun_nodes(char op, Node *t1, Node *t2);


int init_tables(void);

#endif
