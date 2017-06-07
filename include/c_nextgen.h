#ifndef NEXTGENH
#define NEXTGENH


// global preprocessor directives for c_nextgen.h


#define MAXCONSTS 100

#define TREEBUFFER 120000000 // tree array length

// shared vars for c_nextgen.h

extern char treebuffer[TREEBUFFER];
extern int node_count;


//const char* delims = "\n";
//const char delim = ':';  // is this really used?


// global function definitions
void split(char text[],int *I, char left_bracket, char delim);

Point i2point(int i, int max_x);
int *make_int(int value);
double *make_double(double value);

int NodeScore2c_string(NodeScore ns, char treebuffer[] );
int read_text(char *filepath,char *dest_str );
int free_int_pars(IntPars ip);
int free_node_score(NodeScore ns);
int free_pop(Population pop);
IntPars make_int_pars(int size);
NodeScore make_node_score(int size);
int check_numpar(int numpar);
int is_null_node(Node *tree);
Node *talloc(void);
char* node2str(Node *tree, char trstr[]);
Node *child_choice(Node *t);
char string_choice(char str[]);
int selectindex(int ntrees, double pexp);
Node *makerandomtree(int maxdepth,double fpr,double ppr, int max_par);
NodeScore makerandomns(int maxdepth,double fpr,double ppr, Experiment Exp);
int is_terminal(Node *t);
IntPars crossover_int_pars(IntPars ip1, IntPars ip2);
NodeScore crossover_ns(NodeScore ns1,NodeScore ns2, double probswap, int top);
Node *unifcross_node(Node *t1,Node *t2, double probswap, int top);
Node *crossover(Node *t1,Node *t2,double probswap,int top);
int mutconsts(double *consts[],int size, double mutationrate);
Node *mutate(Node *t, double probchange);
Node *swapmut(Node *t, double probchange);
Node *hoistmut(Node *t, double probreturn,int top);
NodeScore allmut(NodeScore ns, double probchange, Experiment Exp);
int compare_nodes(Node *t1, Node *t2);
int compare_ns(NodeScore ns1, NodeScore ns2);
int does_tree_exist(NodeScore newns, Population pop, Point my_loc, int compgridsize);
int conparcount(Node *t, int *current);
int init_buffer(void);
double my_strtod(const char *token);
NodeScore bufferline2NodeScore(char tmp_str[], char * infile);
int check_brackets(char *tmp_str, char left_bracket);
int check_ascii(char *tmp_str);
int read_data(const char *filepath,Node *trees[],int maxtrees);
void c_nextgen(int my_number,int qsubs,int runlen,int popsize, int compgridsize, double mutationrate, int tour);

#endif
