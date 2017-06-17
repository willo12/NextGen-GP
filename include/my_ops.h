
#ifndef OPS_H
#define OPS_H

#ifdef INTSTATES
extern int (*op_table[OPTABLE])(Node *tree);
#else
extern double (*op_table[OPTABLE])(Node *tree);
#endif


extern char* (*code2c_str_table[OPTABLE])(Node *tree, char trstr[]);
extern char* (*code2str_table[OPTABLE])(Node *tree, char trstr[]);

extern Node* (*copy_table[OPTABLE])(Node *tree);
extern Node* (*str2leaf_table[OPTABLE])(char text[], char left_bracket, char delim);

Node *makerandomterminal(int maxdepth,double fpr,double ppr, int max_par);

int tree_height(Node *tree, int height);

void free_node(Node *t);
int check_node(Node *t);
int tree_consts(Node *t, double *consts[], int *current, int maxconsts);

int leaves_equal(Node *t1,Node *t2);
void user_functions(NodeScore ns);

int init_tables_ops(void);
int init_tables_copy(void);
int init_io_tables(int arg_table[]);

void evaluate_ns(NodeScore ns, State S_return);
void evaluate_tree(Node* tree, State S_return);
double run_node(Node* node);

#endif
