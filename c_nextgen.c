/* ISSUES: MAXTREESTR may be exceeded during runtime: need to account for this.
One way is to calculate expected string length during tree construction?

*/

//#define SPECIAL_MUL

#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* MAXCHILD must be >= SPACEDIM */
#define TANHONLY
#define INTERP_FORCING
#define HI_MUT_SPOT

#define SPACEDIM 2
#define RESPWINDOW 500
#define MAXLEAF 25
#define MAXSIZE 12
#define HILLTHRESH 0.02
#define INITNODES 1
#define MAXHILL 8
#define CONSTSCUTOFF 8
#define MIGRATE 0.02
#define TOUR 7
#define MIGTOUR 5
#define RNDGRAIN 0.005
#define MAXSEARCH 0
#define MAXCONSTS 100
#define MAXCHILD 3
#define MAXPAR 10
#define NUMPAR 7
#define OPTABLE 26
#define MAXTREESTR 10000
#define TREEBUFFER 120000000
#define ERRORMARG 0.1
#define MIN_SCORE 40

#define BLOWUP 10
#define MAXCHANGE 0.1
#define PARSIMONY 0.5
#define MAXDEPTH 4
#define FPR 0.5
#define PPR 0.6
#define PNEW 0.001
#define BREEDINGRATE 0.3
#define MUTATIONRATE 0.08

#define PEXP 0.9 
#define PROBSWAP 0.9

#define PULSETIME 10000
#define MAXTREES 360000

#define FNAMESIZE 100
#define STACKSIZE 2000

int numpar;
double buffers[STACKSIZE][SPACEDIM];

struct node *talloc(void);
int *make_int(int);
double *make_double(double);

void interp_node(struct node *tree, int i_tofill);
char* node2str(struct node *tree, char trstr[]);

char* instrsetnodes = "QEOTLASMIV"; /* for printing */
char* instrset = "ASM"; /* used at runtime */
char* instrsetswap = "ASM"; /* instr that can swap */
char* instrsetzero = "T"; /* function subleaf nodes */
/* char* instrsetzero = "QEOTL";  */

double reg[MAXPAR];

double accum, xreg;

double * fixedforc;
double * fixedforc2;
double * emph;

char treebuffer[TREEBUFFER];

char treestr1[MAXTREESTR];
char treestr2[MAXTREESTR];

char* (*code2str_table[OPTABLE])(struct node *tree, char trstr[]);
void (*op_table[OPTABLE])(struct node *, int);

struct node *copy_node(struct node *tree);
void free_node(struct node *t);


double get_score(double* ffs,double* result, double *obs, long* I, struct node *newtree, int ffs0, int ffs1, int obs0, int obs1, double aspect, double* S_init, int ts_factor, int startscore_i);
int mutconsts(double *consts[],int size);
double hillclimb(double* ffs, double* result,double*  obs,long*  I,struct node *newtree, int ffs0, int ffs1, int obs0, int obs1, double *consts[] , int size, double error, double aspect, double* S_init, int ts_factor, int startscore_i);
int sectortour(double *score_vals,int popsize,int tour, int sector, int adjacent[]);
int tournament_size(struct node *trees[], int is,int popsize,int tour);
int tournament_scsize(struct node *trees[],double *score_vals, int is,int popsize,int tour);
int reverse_tour_size(struct node *trees[],double *score_vals, int is,int popsize,int tour);
int tournament(double *score_vals, int is,int popsize,int tour);
int reversetour(double *score_vals,int is,int popsize,int tour);
int * mapsector(int sector,int popsize);
void  init_fixedforc(double *result, int result0, double *model,int model0);
int check_brackets(char *tmp_str);
int check_node(struct node *t);
int init_tables(void);
void handle_S_init(double S_init, double* obs, int obs1, double* S_init_array);

static int arg_table[OPTABLE];

const char* delims = "\n";

const char delim = ':';

int node_count = 0;

int once_flag=0;

int nhillsearches=0;

FILE *fp;

/*

To create a const node:
node op is 'C'
first child points to a value

*/

struct fun {
  void * function;
  unsigned char argtypes[MAXCHILD];  
};

struct node {
  int op;
  void * children[MAXCHILD];
};

int mod (int a, int b)
{
  if (b<0)
    return mod(-a,-b);

  int ret = a%b;
  if (ret < 0)
    ret+=b;
  return ret;

}

void split(char text[],int *I)
{
/* 
split string on ',' except for comma's within bracketed parts.
split is achieved by replacing those ',' with \0 and maintaining
indices into string in I, indicating where the pieces start.

(p1,2.56)S,p2

becomes 

'(p1,2.56)S', 'p2'

*/

  int counter=0;
  int cursor=0;
  int icomma=0;

  
  while (*(text+cursor) !='\0')
  {

    if (*(text+cursor)=='(')
      counter += 1;
    else if (*(text+cursor)==')')
      counter -= 1;

    if (*(text+cursor)==',' && counter==0)
    {
      
      text[cursor]='\0';
      *(I+icomma)=cursor+1;      
      icomma++;
    } 

    cursor+=1;
  }

  *(I+icomma)=0;  
}

struct node *copy_node(struct node *tree)
{
  struct node *newnode=talloc();
  newnode->op = tree->op;

  /* C and P are special as their children are non-standard nodes */
  if (tree->op == 'C'-'A')
  {
    newnode->children[0] = make_double(*((double *) tree->children[0]));

  }
  else if (tree->op == 'P'-'A')
  {
    newnode->children[0] = make_int(*((int *) tree->children[0]));
  }
  else
  { /* copy children. copies no children in case of zero-child nodes */
    int i;
    for (i=0;i<arg_table[tree->op];i++)
    {
      newnode->children[i] = copy_node(tree->children[i]);
    };
  };

  return newnode;
};


struct node *str2node(char text[])
{

  /* recursively fill node from text containing tree in bracket notation */

  int l,i;
  int I[4];

  if (text[0] == '*')
  {
    return NULL;
  }

  struct node *newnode=talloc();
  
  if (text[0] == '(') /* this is not a leaf node */
  {
    l=strlen(text);
    newnode->op=text[l-1]-'A';  /* found op. op code is relative to 'A' */

/* trim */
    text[l-2]='\0'; /* remove trailing op and ) bracket */
    text = text+1; /* remove beginning ( bracket */

/* split string on ',' except for comma's within bracketed parts.
split is achieved by replacing those ',' with \0 and maintaining
indices into string in I, indicating where the pieces start. */
    split(text,I);

    /* fill children of newnode, using str2node recursively */
    newnode->children[0] = str2node(text);    /* first child at start of string (\0 have now been inserted) */

    /* collect remaining children */
    i=0;
    while (I[i] !=0)
    {
      newnode->children[i+1] = str2node(text+I[i]);    
      i++;
    }
  }    
  else /* this is a sub-leaf node (identified by not having a bracket) */
  {
    if (text[0]=='p') /* it's a parameter */
    {
      newnode->op='P'-'A'; 
      newnode->children[0]=make_int(atoi(text+1));
    }
    else if ( isdigit(text[0]) || (text[0]=='-') )
    {
      newnode->op='C'-'A'; /* it's a  constant */
      newnode->children[0]=make_double(atof(text));
    }
    else
    {
      newnode->op=text[0]-'A'; /* it's a sub-leaf function, no children */
    }
  }
  return newnode;
}


struct node *talloc(void)
{
  struct node *newtree =  (struct node *) calloc(1,sizeof(struct node));
  if (newtree == NULL)
  {
    fprintf(stderr,"Error: cannot calloc newtree in talloc.\n");
    exit(0);
  }

  int i;
  for (i=0;i<MAXCHILD;i++)
  {
    newtree->children[i]=NULL;
  }

  node_count++;
  return newtree;
};

int *make_int(int value)
{
/*
  Create new integer constant node with malloc. Returns pointer to new int.
*/
  int *ptr = (int *)malloc(sizeof(int));

  if (ptr == NULL)
  {
    fprintf(stderr,"Error: cannot malloc new int in make_int.\n");
    exit(0);
  }

  *ptr=value;
  return ptr;
};

double *make_double(double value)
{
/*
  Create new double constant node with malloc. Returns pointer to new double
*/

  double *ptr = (double *)malloc(sizeof(double));

  if (ptr == NULL)
  {
    fprintf(stderr,"Error: cannot malloc new double in make_double.\n");
    exit(0);
  }

  *ptr=value;
  return ptr;
};

/* Functions ending in Chr print node to string. Used in code2str_table. 
   Names correspond to node names.
*/

char* nopChr(struct node *tree, char trstr[])
{
  return trstr;
};

char* constChr(struct node *tree, char trstr[])
{
/* convert tree to str where tree is a constant leaf */
  sprintf(trstr,"%g",*((double *) tree->children[0]));

  return trstr; 
};


char* parChr(struct node *tree, char trstr[])
{
/* convert tree to str where tree is a par leaf */

  sprintf(trstr,"p%d",*((int *) tree->children[0]));

  return trstr; 
};





/*
  The following four functions print functions taking 0,1,2 and 3 argument nodes.
*/

char* zeroFunChr(struct node *tree, char trstr[])
{
/* convert tree to str recursively, where root node is an op that takes 0 arguments */

  sprintf(trstr,"%c", tree->op+'A'); /* print Q etc to trstr */

  return trstr;

};

char* oneFunChr(struct node *tree, char trstr[])
{
/* convert tree to str recursively, where root node is an op that takes 1 argument */

  char tmp_str[MAXTREESTR];

  /* insert tmp_str between brackets, e.g. ()M */
  strcpy(trstr,"(");  

/* code2str_table table linking op codes to function pointers, functions convert node to str in brack notation.
   code2str_table points to oneFunChr (this function), etc and leaf chr function parChr, constChr 
   for example, ((p1)Q)Q would call oneFunChr again on (p1)Q and yield (result)Q, where result is obtained by oneFunChr calling parChr on p1 and enclosing it in ()Q.

*/

/*   apply code2str_table[op] to children and put result in tmp_str */
  (*code2str_table[((struct node *) tree->children[0])->op])( ((struct node *) tree->children[0]),tmp_str );
  
  strcat(trstr,  tmp_str   );  /* append tmp_str to trstr */
  sprintf(tmp_str,")%c", tree->op+'A'); /* print )M or )A etc to tmp_str */
  strcat(trstr, tmp_str); 

  return trstr;

};

/* new generalized 2,3,..FunChr to be used later */
char* nFunChr(struct node *tree, char trstr[], int n)
{
/* convert tree to str recursively, where root node is an op that takes n arguments */
 
  int i;
  char tmp_str[MAXTREESTR];
  strcpy(trstr,"(");

  for (i=0;i<n;i++){  
    (*code2str_table[((struct node *) tree->children[0])->op])( ((struct node *) tree->children[0]),tmp_str );
  strcat(trstr,tmp_str);
    if (i<n-1) /* join on , delimiter */
    {
      strcat(trstr,",");
    }
  } 


  sprintf(tmp_str,")%c", tree->op+'A');
  strcat(trstr, tmp_str);

  return trstr;

};

char* twoFunChr(struct node *tree, char trstr[])
{
/* convert tree to str recursively, where root node is an op that takes 2 arguments */
  char tmp_str[MAXTREESTR];
  strcpy(trstr,"(");
  
  (*code2str_table[((struct node *) tree->children[0])->op])( ((struct node *) tree->children[0]),tmp_str );
  strcat(trstr,tmp_str);
  strcat(trstr,",");
  (*code2str_table[((struct node *) tree->children[1])->op])( ((struct node *) tree->children[1]),tmp_str );
  strcat(trstr,tmp_str);  
  sprintf(tmp_str,")%c", tree->op+'A');
  strcat(trstr, tmp_str);

  return trstr;

};


char* threeFunChr(struct node *tree, char trstr[])
{
/* convert tree to str recursively, where root node is an op that takes 3 arguments */
  char tmp_str[MAXTREESTR];
  strcpy(trstr,"(");  
  (*code2str_table[((struct node *) tree->children[0])->op])( ((struct node *) tree->children[0]),tmp_str );
  strcat(trstr,tmp_str);
  strcat(trstr,",");
  (*code2str_table[((struct node *) tree->children[1])->op])( ((struct node *) tree->children[1]),tmp_str );
  strcat(trstr,tmp_str);    
  strcat(trstr,",");
  (*code2str_table[((struct node *) tree->children[2])->op])( ((struct node *) tree->children[2]),tmp_str );
  strcat(trstr,tmp_str);    
  sprintf(tmp_str,")%c", tree->op+'A');
  strcat(trstr, tmp_str);

  return trstr;
};

/* The following functions are called from the table of operations to execute nodes recursively */


double constFun(struct node *tree)
{
/* Yield a constant leaf. Triggered by op code 'C' */
  return *((double *) tree->children[0]); 

};

void nopFun(struct node *tree, int i_tofill)
{
/* no operation function to initialize op_table */

};


double parFun(struct node *tree)
{
/* Yield a parameter leaf. Triggered by op code 'D' */
  return reg[*((int *) tree->children[0])];
};

double zeroFun(struct node *tree, int i, double input)
{
/* subleaf functions. Apply scalar function. Triggered by subleaf op codes e.g. 'Q' 
   int i: determines what input to use for the scalar function: >-1 signals reg[i], -1 signals input argument.

   Unlike branch (above-leaf) functions, which deposit in the buffer, this function returns its value.

*/

  if (i>-1)
  {

    input = reg[i];
  
  }

  if (tree->op=='T'-'A')
  {
    return tanh(input);
  }
  else if (tree->op=='E'-'A')
  {
    return exp(input);
  }
  else if (tree->op=='O'-'A')
  {
    return 1.0/input;
  }
  else if (tree->op=='Q'-'A')
  {
 
    return sqrt(input);

  }
  else if (tree->op=='L'-'A')
  {
    return log(input);
  }

  else
  {
    return 0;
  }
};

double zeroFunInput(struct node *tree, double input)
{
/* subleaf functions. Apply scalar function to input argument. Triggered by subleaf op codes e.g. 'Q' 

*/


  if (tree->op=='T'-'A')
  {
    return tanh(input);
  }
  else if (tree->op=='E'-'A')
  {
    return exp(input);
  }
  else if (tree->op=='O'-'A')
  {
    return 1.0/input;
  }
  else if (tree->op=='Q'-'A')
  {
 
    return sqrt(input);

  }
  else if (tree->op=='L'-'A')
  {
    return log(input);
  }

  else
  {
    return 0;
  }
};


double zeroFunReg(struct node *tree, int i)
{
/* subleaf functions. Apply scalar function. Triggered by subleaf op codes e.g. 'Q' 
   use reg[i] as input

*/


  if (tree->op=='T'-'A')
  {
    return tanh(reg[i]);
  }
  else if (tree->op=='E'-'A')
  {
    return exp(reg[i]);
  }
  else if (tree->op=='O'-'A')
  {
    return 1.0/reg[i];
  }
  else if (tree->op=='Q'-'A')
  {
 
    return sqrt(reg[i]);

  }
  else if (tree->op=='L'-'A')
  {
    return log(reg[i]);
  }

  else
  {
    return 0;
  }
};





// Jump functions

void leafFun(struct node *tree, int i_tofill)
{
/* Represents vector of parameters and constants.

   Calling it fills buffers[i_tofill][:] with constant values and current values of parameters.

   leaf nodes have SPACEDIM amount of (sub-leaf) children.
*/

  int i;  

  for (i=0;i<SPACEDIM;i++)
  { /* interpret child nodes (vector components) here locally instead of via op_table */
    if ( ((struct node *) tree->children[i])->op == 'C'-'A')
    {
      buffers[i_tofill][i] = *((double *) ((struct node *) tree->children[i])->children[0]);
    }
    else if ( ((struct node *) tree->children[i])->op == 'P'-'A')
    {
      buffers[i_tofill][i] = reg[*((int *) ((struct node *) tree->children[i])->children[0])];

                            
    }
    else
    {
#ifdef TANHONLY
      buffers[i_tofill][i] = tanh(reg[i]);
#else
      buffers[i_tofill][i] = zeroFunReg(tree->children[i], i);
#endif

    }
  }
}

void addFun(struct node *tree, int i_tofill)
{
  /* execute addition node: add the values of the two child nodes. */

  /* double buffers[4][SPACEDIM] */

  int i;  
  (*op_table[(  (struct node *) tree->children[0])->op])(((struct node *) tree->children[0]),i_tofill+1); /* fills buffer i_tofill+1 */
  (*op_table[(  (struct node *) tree->children[1])->op])(((struct node *) tree->children[1]),i_tofill+2); /* fills buffer i_tofill+2 */

  for (i=0;i<SPACEDIM;i++)
  {
    buffers[i_tofill][i] = buffers[i_tofill+1][i] + buffers[i_tofill+2][i];
  }

};

void subFun(struct node *tree, int i_tofill)
{
  /* execute subtraction node: subtract the values of the two child nodes. */

  int i;  
  (*op_table[(  (struct node *) tree->children[0])->op])(((struct node *) tree->children[0]),i_tofill+1); /* fills buffer i_tofill+1 */
  (*op_table[(  (struct node *) tree->children[1])->op])(((struct node *) tree->children[1]),i_tofill+2); /* fills buffer i_tofill+2 */

  for (i=0;i<SPACEDIM;i++)
  {
    buffers[i_tofill][i] = buffers[i_tofill+1][i] - buffers[i_tofill+2][i];
  }

};

void mulFun(struct node *tree, int i_tofill)
{
  /* execute multiplication node: multiply the values of the two child nodes. If left node is a function, apply that function to right multiplicant */

#ifdef SPECIAL_MUL
  struct node *child;
#endif
  int i; 

#ifdef SPECIAL_MUL
  if ( ((struct node *) tree->children[0])->op == 'V'-'A' ) /* left multiplicant is a vector V: look at each of its components */
  { 

    /* interpret right multiplicant, and then decide how left multiplicant acts on resulting buffer values for each component */
    (*op_table[(  (struct node *) tree->children[1])->op])(((struct node *) tree->children[1]),i_tofill+1); /* fills buffer i_tofill+1 */

    for (i=0;i<SPACEDIM;i++)
    { /* get children of left multiplicant */
      child = ((struct node *) ((struct node *) tree->children[0])->children[i]);
      if ( child->op == 'C'-'A')
      {
        buffers[i_tofill][i] = (*((double *) child->children[0]))*buffers[i_tofill+1][i];  /* left multiplicant is a constant */
      }
      else if ( child->op == 'P'-'A')
      {
        buffers[i_tofill][i] = reg[*((int *) child->children[0])]*buffers[i_tofill+1][i]; /* it's a parameter node */
      }
      else
      {

#ifdef TANHONLY
        buffers[i_tofill][i] = tanh(buffers[i_tofill+1][i]);
#else

        buffers[i_tofill][i] =  zeroFunInput(child, buffers[i_tofill+1][i]); /* child of left node is function: apply function to right multiplicant */
#endif

      }
    }

  }
  else
  {
#endif

    (*op_table[(  (struct node *) tree->children[0])->op])(((struct node *) tree->children[0]),i_tofill+1); /* fills buffer i_tofill+1 */
    (*op_table[(  (struct node *) tree->children[1])->op])(((struct node *) tree->children[1]),i_tofill+2); /* fills buffer i_tofill+2 */

    for (i=0;i<SPACEDIM;i++)
    {
      buffers[i_tofill][i] = buffers[i_tofill+1][i] * buffers[i_tofill+2][i];
    }

#ifdef SPECIAL_MUL
  }
#endif

};

void divFun(struct node *tree, int i_tofill)
{
  /* execute division node: divide the values of the two child nodes. */

  int i;  
  (*op_table[(  (struct node *) tree->children[0])->op])(((struct node *) tree->children[0]),i_tofill+1); /* fills buffer i_tofill+1 */
  (*op_table[(  (struct node *) tree->children[1])->op])(((struct node *) tree->children[1]),i_tofill+2); /* fills buffer i_tofill+2 */

  for (i=0;i<SPACEDIM;i++)
  {
    buffers[i_tofill][i] = buffers[i_tofill+1][i] / buffers[i_tofill+2][i];
  }

};

void ifFun(struct node *tree, int i_tofill)
{
  /* execute if-then branching node: execute child 1 if child 0 > 0, child 2 otherwise */

  int i;  

  (*op_table[(  (struct node *) tree->children[0])->op])(((struct node *) tree->children[0]),i_tofill+1); /* fills buffer i_tofill+1 */
  (*op_table[(  (struct node *) tree->children[1])->op])(((struct node *) tree->children[1]),i_tofill+2); /* fills buffer i_tofill+2 */
  (*op_table[(  (struct node *) tree->children[0])->op])(((struct node *) tree->children[2]),i_tofill+3); /* fills buffer i_tofill+1 */


  for (i=0;i<SPACEDIM;i++)
  {

    if (buffers[i_tofill+1][i]>0)
      buffers[i_tofill][i] = buffers[i_tofill+2][i];      
    else 
      buffers[i_tofill][i] = buffers[i_tofill+3][i];      
  }

};


void sqrtFun(struct node *tree, int i_tofill)
{
  /* execute square root node: take the square root of the child node. */

  int i;  
  (*op_table[(  (struct node *) tree->children[0])->op])(((struct node *) tree->children[0]),i_tofill+1); /* fills buffer i_tofill+1 */

  for (i=0;i<SPACEDIM;i++)
  {
    buffers[i_tofill][i] = sqrt( ( double ) buffers[i_tofill+1][i] );
  }

};



/* The following functions deal with intpereting trees and converting them to strings */

void interp_node(struct node *tree, int i_tofill)
{
/* Interpret the argument node recursively using the table of operations op_table */

  (*op_table[tree->op])(tree, i_tofill);
};

char* node2str(struct node *tree, char trstr[])
{
/* convert the argument node to bracket string trstr recursively using the table of codes code2str_table 
   convenience function for code2str_table
*/
  return (*code2str_table[tree->op])(tree,trstr);
};

/* Random selection utilities */

struct node *child_choice(struct node *t)
{
  /* Choose a child node randomly */
  return ((struct node *) t->children[rand()%arg_table[t->op]]);
};

char string_choice(char str[])
{
  /* Choose a character inside a string randomly */
  int n = strlen(str);
  int i = rand()%n;

  return str[i];
};


int selectindex(int ntrees, double pexp)
{
  /* Select an index between 0 and ntrees, but with bias towards 0. Used in elitism */
  int result = (int)floor(log( ( (double) rand() )/RAND_MAX ) /log(pexp));
  if (result > ntrees-2)
    result = ntrees-2;
  return result;
};

double randomdouble(void)
{
  /* Produce a random floating point number, with some biases */

  int i = rand()%8;

  if (i<1)
  {
    if (rand()%2==0)
    {
      return 1.0;
    }
    else
    {
      return 0.0;
    }
  }
  else if (i<2)
  {
    return (( (double) (2*rand() -RAND_MAX) )/RAND_MAX)*pow(10,rand()%11-5);
  }
  else if (i<5)
  {
    return (( (double) (2*rand() -RAND_MAX) )/RAND_MAX)*pow(10,rand()%5-2);
  }

  else
  {
    return (( (double) (2*rand() -RAND_MAX) )/RAND_MAX)*pow(10,rand()%7-3);
  }
}

struct node *makerandomtree(int maxdepth,double fpr,double ppr)
{
  /* Create a random tree 
     fpr probability of creating normal node (recursively)
     ppr probability of creating parameter node

  */

  char op;
  int i=0;
  int rnd;
 
  struct node *child;
  struct node *result=talloc();

  rnd = rand();

  if ( (rnd<fpr*RAND_MAX) && (maxdepth>0) )
  { /* only recursion: create random node and call makerandomtree to make children */
    op = string_choice(instrset)-'A';
    result->op = op;

    while ( i<arg_table[(int) op] )  
    {  
      result->children[i] = ((struct node *)  makerandomtree(maxdepth-1,fpr,ppr) );
      i++;
    };
  }
  else
  {  /* Terminal. Create random V node. */
    result->op = 'V'-'A';
    
    for (i=0; i<SPACEDIM;i++ )  
    {
      child = talloc();
      rnd = rand();
      if ( rnd<ppr*RAND_MAX  )
      {
        if (rand()%2 == 0)
        {        
  /* make parameter node */
          child->op = 'P'-'A';
          if (i==0)
          {
            child->children[0] = make_int((int)  rand()%numpar );
          }
          else
          { /* exclude forcing in other state variables */
            child->children[0] = make_int((int)  rand()%SPACEDIM );
          }


        }
        else
        {        
  /* make function node */
          child->op = string_choice(instrsetzero)-'A';
        }


      }
      else
      {
  /* make constant node */
        child->op = 'C'-'A';
        child->children[0] = make_double(  randomdouble()  );
      };
      result->children[i] = ((struct node *) child);
     
    };    
  };
  return result;  
};


struct node *mrt_steered(int maxdepth,double fpr,double ppr)
{
  char op;
  int i=0;
  int rnd;

  char tmp_str[MAXTREESTR];


  struct node *result;

  rnd = rand();

  if ( (rnd<fpr*RAND_MAX) && (maxdepth>0) )
  { /* make normal node */

    if (rand()%3 != 0)
    {
      sprintf(tmp_str,"(%g,(p0,%g)S)M",randomdouble(),randomdouble());
      result = str2node(tmp_str);

 /*     ( (struct node *)  result->children[1])->children[1] = ((struct node *)  mrt_steered(maxdepth-1,fpr,ppr) );
*/
    }
    else
    {
      result=talloc();
      op = string_choice(instrset)-'A';
      result->op = op;

      while ( i<arg_table[(int) op] )  
      {  
        result->children[i] = ((struct node *)  mrt_steered(maxdepth-1,fpr,ppr) );
        i++;
      };
    }

  }
  else if ( rnd<ppr*RAND_MAX  )
  {
  /* make parameter node */

    result=talloc();

    result->op = 'P'-'A';
    result->children[0] = make_int((int)  rand()%numpar );
  }
  else
  {
  /* make constant node */
    result=talloc();

    result->op = 'C'-'A';
 
    result->children[0] = make_double(  randomdouble()  );
  };

  return result;  
};



int leaves_equal(struct node *t1,struct node *t2)
{
// return 0 if t1 and t2 are leaf nodes that are equal, 1 otherwise


//  struct node *child;

//  int return_val = 1;

  int i;

  if ( (  ((struct node *) t1)->op =='V'-'A')  &&  (    ( (struct node *) t2)->op =='V'-'A')   )
  {


    for (i=0;i<SPACEDIM;i++)
    { /* get children of left multiplicant */

      if (  ((struct node *) t1->children[i])->op !=  ((struct node *) t2->children[i])->op )
      {
        return 1;
      }
      else if ( ( ( (struct node *) t1->children[i])->op == 'C'-'A') && ( (*((double *) ( (struct node *) t1->children[i])->children[0])) != (*((double *) ( (struct node *) t2->children[i])->children[0])) ) )
      {

        return 1;
      }
      else if ( ( ( (struct node *) t1->children[i])->op == 'P'-'A') && ( *((int *) ( (struct node *) t1->children[i])->children[0]) != *((int *) ( (struct node *) t2->children[i])->children[0])  )  )
      {
        return 1;
      }
  

    }

  return 0;
  }

  return 1;


}

struct node *unifcross(struct node *t1,struct node *t2, double probswap, int top)
{
  /* Uniform crossover

  */
  int i=0;
  struct node *result;

  if ( (top==0) && (rand()>probswap*RAND_MAX) )
  { /* terminal. swap t1 for t2 */
 
    return ((struct node *) copy_node(t2));
  }
  else
  {
    if ( (t1->op =='V'-'A')  ||  (t2->op =='V'-'A')   )
    { /* terminal: t1 or t2 is a leaf node, then randomly return a copy of t1 or t2, this to achieve the mixing */


      if ( (rand()>probswap*RAND_MAX) && (leaves_equal(t1,t2) != 0) )
      { 

        return copy_node(t2);
      }
      else
      {
        return copy_node(t1);
      }
    }
    else if (arg_table[t1->op] == arg_table[t2->op])
    { /* only recursion: t1 and t2 are not leaf nodes, an take same number or args 
         now commit space to result node pointer (goes out of scope in other if-then branches) 
         heads or tails whether t1 or t2 op is used for result, and apply unifcross pairwise to children
      */

      result = talloc();
      if (  (rand()>probswap*RAND_MAX) && (t1->op != t2->op  )  )
      {

        result->op = t2->op; 
      }
      else
      {
        result->op = t1->op; 
      }

      while (i<arg_table[t1->op])
      {
        result->children[i] = ((struct node *)  unifcross(((struct node *) t1->children[i]) , ((struct node *) t2->children[i]),probswap, 0) );
        i++;
      }
      return result;
    }
    else
    {
      /* t1 and t2 are not leaf nodes and do not take same number or args, return copy of t1 */
      return copy_node(t1);
    }
  }
}




struct node *crossover(struct node *t1,struct node *t2,double probswap,int top)
{
  /*
    
  */

  int i=0;
  struct node *result;

  if ((rand()<probswap*RAND_MAX) && (top==0))
  {
    return ((struct node *) copy_node(t2));
  }
  else
  {
  
    if ((t1->op !='C'-'A') && (t1->op !='P'-'A') && (t2->op !='C'-'A') && (t2->op !='P'-'A'))
    {
      result = talloc();
      result->op = t1->op;

      while ( (i<arg_table[t1->op])  )
      {
  
        result->children[i] = ((struct node *)  crossover(((struct node *) t1->children[i]) , child_choice(t2),probswap, 0) ) ;
        i++;
      };
    }
    else
    { /* leaf nodes return copies of node t1 depending on dice and top!=0 */ 
      result=((struct node *) copy_node(t1));
    }
    return ((struct node *) result);
  };
};




struct node *mutsteered(struct node *t, double probchange)
{
  int i=0;
  struct node *result;

  if ( rand()<probchange*RAND_MAX )
  {
    return mrt_steered(MAXDEPTH,FPR,PPR);  
/*    return makerandomtree(MAXDEPTH,FPR,PPR); */
  }
  else
  {
    if ((t->op !='C'-'A') && (t->op !='P'-'A') )
    {
      result = talloc();
      result->op = t->op;

      while (i<arg_table[t->op])  
      {
  
        result->children[i] = ((struct node *)  mutsteered(((struct node *) t->children[i]) ,probchange) );

        i++;
      };  
    }     
    else
    {
      result = copy_node(t);
    }
    return result;
  }
};


struct node *mutate(struct node *t, double probchange)
{
  /* 
    Go through tree recursively and roll dice to replace a node and graft new random tree. 
    C/ P leaf nodes are not touched.

    e.g. in (p1,2.0)M p1 and 2.0 can be replaced with a subtree
    
  */

  int i=0;
  struct node *result;

  if ( rand()<probchange*RAND_MAX )
  {
    /* dice successful: grow random tree here */
    return makerandomtree(MAXDEPTH,FPR,PPR);
  }
  else
  {
    if (t->op !='V'-'A') 
    { /* only recursion: it is not a leaf node, preserve op but apply mutate to children */
      result = talloc();
      result->op = t->op;

      while (i<arg_table[t->op])  
      { /* recursively throw dice to replace children with subtrees */
        result->children[i] = ((struct node *)  mutate(((struct node *) t->children[i]) ,probchange) );

        i++;
      };  
    }     
    else
    { /* terminal: return copied node t for leaf */
      result = copy_node(t);
    }

    return result;
  }
};

struct node *swapmut(struct node *t, double probchange)
{
/* Swap mutation, rolls dice to swap out a suitable operation somewhere in the tree with one from instrsetswap. 
   If node t takes 2 args, swapmut replaces op in node with a randomly chosen 2-arg op 
   at most 1 node is swapped.
*/

  int i=0;
  struct node *result;

  if (t->op !='V'-'A') 
  {
    if ( (rand()<probchange*RAND_MAX) && (arg_table[t->op]==2) )
    {
      /* terminal: dice successful and arg condition met: copy t with changed op */
      result = copy_node(t);
      result->op = string_choice(instrsetswap)-'A';
      return result;
    }
    else
    {
        /* only recursion: create new node preserving op and applying swapmut to children */
        result = talloc();
        result->op = t->op;

        while (i<arg_table[t->op])  
        {
          result->children[i] = ((struct node *)  swapmut(((struct node *) t->children[i]) ,probchange) );
          i++;
        }
    }
  }
  else
  {
    /* terminal: leaf node */
    result = copy_node(t);
  }

  return result;
};



struct node *hoistmut(struct node *t, double probreturn,int top)
{
  /* hoist mutation: throw dice to replace nodes with one of their children.
     acts to shorten/ simplify trees again.
 */
 
  struct node *result;

  if (top==1)
  { /* only applied on root: root node always copied */
    return copy_node(  hoistmut(t,probreturn,0) ); 
  }

  if (rand()<probreturn*RAND_MAX )
  { /* terminal: unsuccessful dice return t */

    if (t->op !='V'-'A')  // not leaf node
    { /*  only recursion. Hoist up one of the child nodes */
      result = ((struct node *)  hoistmut(((struct node *) child_choice(t) ) ,probreturn,0) );
    }     
    else // leaf node
    { /* terminal: leaf nodes returned as is  */
      result = t;
    }
  }
  else
  {

    result = t;

  }
  return result;
};

struct node *allmut(struct node *t, double probchange)
{
  int rnd;

  rnd = rand()%4;
  if (rnd < 1)
  {
    return hoistmut(t, probchange,1);
  }
  else if (rnd < 2)
  {
    return swapmut(t, probchange);
  }
  else
  {
    return mutate(t, probchange);
  }
};



/*
To create int:

ptr = (int *)malloc(sizeof(int));

create the doubles and ints this way for constant and param leaf nodes

The context of the children must come from the function of the node,
yielding the type

*/


void c_doloop(double* ffs,double* result, struct node *tree, int m, int n, double* S_init, int ts_factor)
{
/* integration loop using simple time stepping.

 result is on the forcing grid */


  int t0=ffs[0];
  int tend=ffs[(m-1)*n];
  double dt;
  int steps=0;
  int forcing_t=0;
  int t_elapsed=0;
  int i;
  double dt_forc;
  int steps_forc=0;

/*  double S_init=0.3; 
  double S_init=0.345;
*/

  double S[SPACEDIM];

/*  double S_mem;
  double Sdotdt=0;
*/

/*  double lincoef = ((double) 1.0)/((double) (0-t0)); */
/*  S[0]=S_init; */
/*  S_mem=S;
*/

/*
  for (i=1;i<SPACEDIM;i++)
  {
    S[i] = 0;
  }
*/

  for (i=0;i<SPACEDIM;i++)
  {
    S[i] = S_init[i];
  }

  dt_forc = ffs[n]-ffs[0];   /* time index to forcing array */
  dt = dt_forc/ts_factor;

/*  fprintf(stderr,"%d %d ", dt_forc, dt);
*/
  
  while ((t_elapsed<tend-t0) && (steps_forc<m))
  {
    forcing_t= (int) round(t_elapsed/dt_forc);    /* time index to forcing array */
    reg[0]=S[0];
    reg[1]=S[1];
/*    reg[2]=ffs[forcing_t*n+1]; */

  /* remaining registers reserved for forcing */
    for (i=0;i<n-1;i++)
    {
      reg[2+i]= ffs[forcing_t*n+1+i];
    }


/*    reg[2]=fixedforc2[steps];
    reg[3]=d;
    reg[2]=fixedforc[steps];

    reg[3]=Sdotdt;
    reg[3]=t_elapsed;
    reg[3]=0.0;
*/
    for (i=0;i<SPACEDIM;i++)
    {
      *(result+(steps_forc*SPACEDIM)+i ) = S[i];
    }


    t_elapsed += dt;


    (*op_table[tree->op])(tree,0);

/*    if (fabs(buffers[0][0])>0.0004 )
    {
      *result=1e19;
      break;
    }  
*/

    for (i=0;i<SPACEDIM;i++)
    {
      S[i] += buffers[0][i]*dt;
    }

/* degeneracy checks */
    for (i=0;i<SPACEDIM;i++)
    {
      if (fabs(S[i])>BLOWUP || isnan(S[i]) || (fabs(buffers[0][i])>MAXCHANGE/dt ) )
      {
        t_elapsed = tend-t0;
        i = SPACEDIM;
       
        *result=1e19;
        break;
      }
    }

    if (steps%ts_factor == 0)
    {
      steps_forc++;
    }

    steps++;
  };
  fprintf(stderr,"steps_forc: %d\n",steps_forc);

}


void f_t(double *S,int n, double t_elapsed, double dt_forc, double *ffs, struct node *tree)
{ /* deposit S values and forcing in registers reg[] and run tree */
  /* n: number of forcing terms plus time (so n=2 for one forcing term) */

  int i, forcing_t;
  double remainder;

  /* first SPACEDIM registers reserved for state vector S */
  for (i=0;i<SPACEDIM;i++)
  {
    reg[i]=S[i];
  }


  /* explicit t dependence */
/*  reg[2]=ffs[forcing_t*n+1]; */

  forcing_t = ((int) floor(t_elapsed/dt_forc));  /* time index to forcing array */

# ifdef INTERP_FORCING
  remainder = t_elapsed - forcing_t*dt_forc;
# endif

  /* remaining registers reserved for forcing. Fill them with forcing */
  for (i=0;i<n-1;i++)
  {
# ifdef INTERP_FORCING
    reg[SPACEDIM+i]= ffs[forcing_t*n+1+i] + remainder*(ffs[(forcing_t+1)*n+1+i] - ffs[forcing_t*n+1+i])/dt_forc;
# else
    reg[SPACEDIM+i]= ffs[forcing_t*n+1+i];
#endif
  }

  (*op_table[tree->op])(tree,0); /* run tree and deposit in buffers[:] */
}


void f_t_2D(double *S,int n, double t_elapsed, double dt_forc, double *ffs, struct node *tree)
{ /* deposit S values and forcing in registers reg[] and run tree */

  /* n: number of forcing terms plus time (so n=2 for one forcing term) */

  int i, forcing_t;
  double remainder;

  /* first 2 registers reserved for state vector S */
  reg[0]=S[0];
  reg[1]=S[1];
  /* explicit t dependence */
/*  reg[2]=ffs[forcing_t*n+1]; */

  forcing_t = ((int) floor(t_elapsed/dt_forc));  /* time index to forcing array */

//  remainder = fmod(t_elapsed, dt_forc);

#ifdef INTERP_FORCING
  remainder = t_elapsed - forcing_t*dt_forc;
#endif

  /* remaining registers reserved for forcing */
  for (i=0;i<n-1;i++)
  {

#ifdef INTERP_FORCING
    reg[2+i]= ffs[forcing_t*n+1+i] + remainder*(ffs[(forcing_t+1)*n+1+i] - ffs[forcing_t*n+1+i])/dt_forc;
#else
    reg[2+i]= ffs[forcing_t*n+1+i];
#endif

  }


  (*op_table[tree->op])(tree,0); /* run tree and deposit in buffers[:] */
}

void RK4(double* ffs,double* result, struct node *tree, int m, int n, double *S_init, int ts_factor)
{
/*
  Runge Kutta 4 solver of tree equation. 

  Trees are executed recursively by calling (*op_table[tree->op])(tree,0), where the 0 arguments tells functions to deposit results in vector buffers[0][:]. Calling the tree reads values from reg[:].

  reg[2] is forcing, i.e. the explicit time dependence of the function f.

  result is on the forcing grid, indexed with steps_forc

  m: shape[0] of forcing ffs, time dimension array length of forcing
  n: shape[1] of forcing ffs, the amount of forcing terms plus time (so n=2 for one forcing term)


*/

  int t0=ffs[0];
  int tend=ffs[(m-1)*n];
  double dt;
  int steps=0;
  double t_elapsed=0;
  int i;
  double dt_forc;
  int steps_forc=0;
  double S_dot_dt[SPACEDIM];

/*  double S_init=0.3; 
  double S_init=0.345;
*/

  double S[SPACEDIM];
  double S_arg[SPACEDIM];

/* Runge Kutta k's */
  double k1[SPACEDIM]; 
  double k2[SPACEDIM]; 
  double k3[SPACEDIM]; 
  double k4[SPACEDIM];

/*  S[0]=S_init; */

  for (i=0;i<SPACEDIM;i++)
  {
    S[i] = S_init[i];
  }

  dt_forc = ffs[n]-ffs[0]; /* time step for forcing */
  dt=dt_forc/ts_factor; /* time step for integration loop */

/*  fprintf(stderr,"%d %d ", dt_forc, dt);
*/
  
  while (t_elapsed<tend-t0) /* integration loop */
  {
    for (i=0;i<SPACEDIM;i++)
    {
      *(result+(steps_forc*SPACEDIM)+i ) = S[i];
    }

    t_elapsed += dt;

    /* RK4 term k1 */
    f_t_2D(S, n, t_elapsed, dt_forc, ffs, tree);
    for (i=0;i<SPACEDIM;i++)
    { 
      k1[i] = buffers[0][i];
    }

    /* RK4 term k2 */
    for (i=0;i<SPACEDIM;i++)
    {
      S_arg[i] = S[i] + k1[i]*dt/2;
    }
    f_t_2D(S_arg, n, t_elapsed + dt/2, dt_forc, ffs, tree);  /* deposit S values and forcing in registers reg[] and run tree */
    for (i=0;i<SPACEDIM;i++)
    {
      k2[i] = buffers[0][i];
    }

    /* RK4 term k3 */
    for (i=0;i<SPACEDIM;i++)
    {
      S_arg[i] = S[i] + k2[i]*dt/2;
    }
    f_t_2D(S_arg, n, t_elapsed + dt/2, dt_forc, ffs, tree); /* run tree */
    for (i=0;i<SPACEDIM;i++)
    {
      k3[i] = buffers[0][i];
    }

    /* RK4 term k4 */
    for (i=0;i<SPACEDIM;i++)
    {
      S_arg[i] = S[i] + k3[i]*dt;
    }
    f_t_2D(S_arg, n, t_elapsed + dt, dt_forc, ffs, tree); /* run tree */
    for (i=0;i<SPACEDIM;i++)
    {
      k4[i] = buffers[0][i];
    }

    /* increment S */
    for (i=0;i<SPACEDIM;i++)
    {
      S_dot_dt[i] = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])*dt/6.0;
      S[i] += S_dot_dt[i];
    }


/* degeneracy checks */
    for (i=0;i<SPACEDIM;i++)
    {
      if (fabs(S[i])>BLOWUP || isnan(S[i]) || (fabs( S_dot_dt[i] ) > MAXCHANGE  ) )
      {

  /*      fprintf(stderr,"degeneracy at i=%d",i);  */
        t_elapsed = tend-t0;
        i = SPACEDIM;
    
        *result=1e19;
        break;
      }
    }

    if (steps%ts_factor == 0)
    {
      steps_forc++;
    }

    steps++;
  };
}



void fill_reg2(double* ffs, int steps, int ts_factor, int steps_forc, int n)
{
  reg[2]= ffs[steps_forc*n+1]; 
}

void fill_reg2_mod(double* ffs, int steps, int ts_factor, int steps_forc, int n)
{ /* this may need to be replaced with dt_forc and fmod code, a la RK4 */
  reg[2]= ffs[steps_forc*n+1] + (steps%ts_factor)*(ffs[(steps_forc+1)*n+1] - ffs[steps_forc*n+1])/ts_factor;
}

void c_doloop_interp(double* ffs,double* result, struct node *tree, int m, int n, double* S_init, int ts_factor)
{
/* Numerical integration loop. result[0.0] = 1e19 signals failure */


  int t0=ffs[0];
  int tend=ffs[(m-1)*n];
  double dt;
  int steps=0;  /* incremented by 1 during run loop to index result array */
  int steps_forc=0;
  double t_elapsed=0;
  int i;
  double dt_forc;
  int forcing_t;
  double remainder;

/*  double S_init=0.3; 
  double S_init=0.345;
*/

  double S[SPACEDIM];

/*  double S_mem;
  double Sdotdt=0;
*/


/*  double lincoef = ((double) 1.0)/((double) (0-t0)); */
/*  S[0]=S_init; */
/*  S_mem=S;
*/

  for (i=0;i<SPACEDIM;i++)
  {
    S[i] = S_init[i];
  }

  for (i=1;i<SPACEDIM;i++)
  {
    S[i] = 0;
  }

  dt_forc = ffs[n]-ffs[0]; /* time step forcing */
  dt = dt_forc/ts_factor;  /* time step of integration */

  fprintf(stderr,"%g %g ", dt_forc, dt);

  while (t_elapsed<tend-t0-dt) /* to keep step_forc in bounds */
  { /* steps index to result array initialized to 0, incremented by 1 */


   
    reg[0]=S[0];
    reg[1]=S[1];
    forcing_t = ((int) floor(t_elapsed/dt_forc));  /* time index to forcing array */

    remainder = fmod(t_elapsed, dt_forc);

    reg[2]= ffs[forcing_t*n+1] + remainder*(ffs[(forcing_t+1)*n+1] - ffs[forcing_t*n+1])/dt_forc;



/*    fprintf(stderr,"yo: %g",(ffs[(steps_forc+1)*n+1] - ffs[steps_forc*n+1])/ts_factor);
*/


    for (i=0;i<SPACEDIM;i++)
    {
      *(result+(steps_forc*SPACEDIM)+i ) = S[i];
    }

    t_elapsed += dt;
    

    (*op_table[tree->op])(tree,0);

/* check on S_dot. */
/*    omit for now
    if (fabs(buffers[0][0])>0.02/dt ) 
    {
     
      *result=1e19;
      break;
    }  
*/

    for (i=0;i<SPACEDIM;i++)
    {
      S[i] += buffers[0][i]*dt;
    }

/* degeneracy checks */
    for (i=0;i<SPACEDIM;i++)
    {
      if (fabs(S[i])>BLOWUP || isnan(S[i])  )
      {
        t_elapsed = tend-t0-dt;
        i = SPACEDIM;

        *result=1e19;
        break;
      }
    }


    steps++;
    if (steps%ts_factor == 0)
    {
      steps_forc++;
    }
 
  };
/*  fprintf(stderr,"yo: %d, m: %d \n",steps_forc,m);   */
}


void c_map_f(char *treestring, double *A, double *B, int shape_i, int shape_j, double S0_start, double S1_start, double S0_end, double S1_end, double *forcing, int forcing_dim)
{
  /* For specific forcing value, create two 2D arrays of tree values for S vector values between -2 and 2. */

  struct node *tree;

  int i,j;
  double dS0, dS1;
 
  char tmp_str[MAXTREESTR];

  init_tables();

  reg[0] = -2.0;
  reg[1] = -2.0;
  for (i=0;i<forcing_dim;i++)
  {
    reg[2+i] = *(forcing+i);
  }


  dS0 = ((double)  (S0_end-S0_start)/shape_i);
  dS1 = ((double)  (S1_end-S1_start)/shape_j);

/*  fprintf(stderr,"%s\n",treestring);
*/
  strcpy(tmp_str,treestring);  

  if (check_brackets(treestring) == -1)
  {
    fprintf(stderr,"Error in tree string: bracket tree integrity check failed. \n");
    /* abort and trigger failure */
    return;    
  }

  tree = str2node(tmp_str);

  if (check_node(tree) == 1)
  {
    fprintf(stderr,"Error: tree integrity check failed.");
        /* abort and trigger failure */
    return;  
  }

/*  printf("dS0: %g, dS1: %g \n",dS0,dS1); */

  for (j=0;j<shape_j;j++)
  {
    reg[0] = -2.0;
    for (i=0;i<shape_i;i++)
    {
      
/*      printf("%g %g %g \n",reg[0],reg[1], reg[2]); */
      (*op_table[tree->op])(tree,0); /* run tree and deposit in buffers[:] */

      A[j*shape_i+i] =  buffers[0][0];
      B[j*shape_i+i] =  buffers[0][1];

      reg[0] += dS0;
    }
    reg[1] += dS1;
  }
  free_node(tree);
}

void f(double *S, struct node *tree)
{

  reg[0]=S[0];
  reg[1]=S[1];
  reg[2]=S[2];

  (*op_table[tree->op])(tree,0); /* run tree and deposit in buffers[:] */
}



int compare_trees(struct node *t1, struct node *t2)
{
/* determine whether 2 trees are identical. returns 0 when identical, 1 otherwise*/

  int j;

  if (t1->op != t2->op)  
  {
    /* unequal ops -> finished */
    return 1;
  }
  else
  {

    if ( (t1->op =='V'-'A') && (leaves_equal(t1,t2) != 0) )
    {
      return 1;

    }
    else
    { /* compare children pair-wise if any. if one of the children returns 1, return 1 */

      for (j=0;j<arg_table[t1->op];j++)  
      {
        if (compare_trees( ((struct node *) t1->children[j]) , ((struct node *) t2->children[j]) ) == 1)
        {
          return 1;
        }
      }
    }
    return 0;
  }
}


int tree_in_old(struct node *newtree,  struct node *trees[], int popsize)
{
  // checks if tree is in array of trees. returns index if yes, -1 otherwise

  int i;

  for (i=0;i<popsize;i++)
  {
    if (compare_trees(newtree, trees[i]) == 0)
    {
      return i;
    }

  }

  return -1;
}

int conparcount(struct node *t, int *current)
{
  /* counts constants and parameters in a tree */
  int j=0;

  if (t->op =='C'-'A')  
  {
    (*current)++;
    return 1;
  }
  else if (t->op =='P'-'A')
  {
    (*current)++;
    return 1;
  }
  else
  {
    while (j<arg_table[t->op])  
    {  
      if ( conparcount( ((struct node *) t->children[j]),current) > -1 )     
      {
        j++;
      }
      else
      {
        return -1;
      }
    }; 
    return *current;  
  }     
  return -1;
};



int tree_consts(struct node *t, double *consts[], int *current, int maxconsts)
{
 
  int j=0;
  
  if (t->op =='C'-'A')  
  {
    if (*current<maxconsts)
    {
      consts[(*current)++] =  ((double *) t->children[0]); 
      return 1;
    }
    else
    {
      return -1;
    }
  }
  else if (t->op =='P'-'A')
  {
    return 0;
  }
  else
  {
    while (j<arg_table[t->op])  
    {  
      if ( tree_consts( ((struct node *) t->children[j]), consts,current,maxconsts) > -1 )     
      {
        j++;
      }
      else
      {
        return -1;
      }
    };
  
    if (*current>=maxconsts)
    {
      return -1;
    }
    return *current;  
  }     
  return -1;
};






void free_node(struct node *t)
{

  if (t == NULL)
  {
    return;
  }
 
  int i;

  if  ((t->op == 'C'-'A') || (t->op == 'P'-'A'))
  {

    if (t->children[0] != NULL)
    {
      free(t->children[0]);
    }
  }
  else
  {
    for (i=0;i<arg_table[(int) t->op];i++)
    {
      free_node(t->children[i]);
    }
  }

  if (t != NULL)
  {
    free(t);
  }

  node_count--;
}

int check_node(struct node *t)
{
  struct node *child;

  if (t == NULL)
  {
    return 1;
  }
 
  int i;

  if  (t->op == 'C'-'A')
  { /* terminal */

    if (t->children[0] != NULL)
    {
      return 0;
    }
    else
    {
      fprintf(stderr,"Warning: illegal const encountered: NULL child!\n");
      return 1;
    }
  }
  else if (t->op == 'P'-'A')
  { /* terminal */
    if ( ( t->children[0] != NULL ) && ( *((int *) t->children[0]) < MAXPAR) && ( *((int *) t->children[0]) >=0 ) )
    {
      return 0;
    }
    else
    {
      fprintf(stderr,"Warning: illegal par encountered!\n");
      return 1;
    }  

  }
  else
  {
    if (strchr(instrsetnodes, t->op + 'A') == NULL)
    {
      fprintf(stderr,"Warning: illegal op encountered!\n");
      return 1;
    }

    if (t->op == 'V'-'A') 
    { /* terminal: check for non sub-leaf children to V */

      for (i=0;i<arg_table[(int) t->op];i++)
      {
        child = ((struct node *) t->children[i]);
        if (child == NULL)
        {
          fprintf(stderr,"Warning: NULL child node encountered!\n");
          return 1;
        }

        if (  (strchr(instrsetzero, child->op + 'A') == NULL) && (child->op != 'C'-'A') && (child->op != 'P'-'A') )
        {
          fprintf(stderr,"Warning: illegal V node encountered!\n");
          return 1;
        }

      }

    }
    else
    { /* check for illegal children in super-leaf nodes */
      for (i=0;i<arg_table[(int) t->op];i++) /* zero-nodes will always return 0  */
      {
        child = ((struct node *) t->children[i]);
        if (child == NULL)
        {
          fprintf(stderr,"Warning: NULL child node encountered!\n");
          return 1;
        }
    
        if (   (strchr(instrsetzero, child->op + 'A') != NULL) || (child->op == 'C'-'A') || (child->op == 'P'-'A') )
        {
          fprintf(stderr,"Warning: illegal super-leaf node encountered!\n");
          return 1;
        }        
        else if (check_node(child) == 1) /* the only recursion: recursively check children */
        { 
          return 1;
        }

      }
    }

    return 0;
  }

  if (t != NULL)
  {
    free(t);
  }

  node_count--;
}


int init_reg(int ffs1)
{
  int i;
  for (i=0;i<MAXPAR;i++)
  {
    reg[i] = 0.0;
  }

  numpar = ffs1 - 1 + SPACEDIM;

  if (numpar > MAXPAR)
  {
    fprintf(stderr,"Warning! Number of variables numpar=%d, with SPACEDIM=%d, forcing terms=%d exceeds MAXPAR=%d \n",numpar,SPACEDIM, ffs1-1,MAXPAR);
    return 1;
  }

  return 0;
}

int init_tables(void)
{
/* Initialize op_table, arg_table and code2str_table
   op_table		table linking op codes to function pointers executed at runtime
   arg_table		table recording number of arguments for each op
   code2str_table 	table linking op codes to function pointers converting nodes to strings

    code2str_table requires instrsetnodes, including the leaf node V.
*/

  int i;

  for (i=0;i<OPTABLE;i++)
  {
    op_table[i] = &nopFun;
  };

  /* linking ops to recursive functions, leafFun is terminal as it produces stored values */
  op_table['A'-'A'] = &addFun;  
  op_table['S'-'A'] = &subFun;
  op_table['M'-'A'] = &mulFun;
  op_table['D'-'A'] = &divFun;
  op_table['I'-'A'] = &ifFun;
  op_table['Q'-'A'] = &sqrtFun;
  op_table['V'-'A'] = &leafFun;


  for (i=0;i<OPTABLE;i++)
  {
    arg_table[i] = 0;
  };

  /* record number of fun args 
     Highest value must be less than MAXCHILD
  */

  i = 0;
  while (instrsetzero[i] != '\0')
  {
    arg_table[instrsetzero[i]-'A'] = 0;
    i++;
  }

  arg_table['A'-'A'] = 2;  
  arg_table['S'-'A'] = 2;
  arg_table['M'-'A'] = 2;
  arg_table['D'-'A'] = 2;
  arg_table['I'-'A'] = 3;

  arg_table['V'-'A'] = SPACEDIM;


/* code2str_table table linking op codes to function pointers, functions convert node to str in brack notation 
   core of the node2str convenience function  */

  for (i=0;i<OPTABLE;i++)
  {
    code2str_table[i] = &nopChr;
  };

  /* these are the true leaf nodes for printing purposes, not V (V is leaf during execution) */
  code2str_table['C'-'A'] = &constChr;  
  code2str_table['P'-'A'] = &parChr;  

  i=0;
  while (instrsetnodes[i] != '\0')
  {
    switch(arg_table[instrsetnodes[i]-'A'])
    {
    case 0 :
      code2str_table[instrsetnodes[i]-'A'] = &zeroFunChr;
      break;
    case 1 :
      code2str_table[instrsetnodes[i]-'A'] = &oneFunChr;
      break;
    case 2 :
      code2str_table[instrsetnodes[i]-'A'] = &twoFunChr;
      break;
    case 3 :
      code2str_table[instrsetnodes[i]-'A'] = &threeFunChr;
      break;
    default :
      code2str_table[instrsetnodes[i]-'A'] = &nopChr;
      break;
    }

    i++;
  }; 
  return 0;
}


int init_buffer(int popsize)
{
/* initialize the buffer with text representations of random trees*/

  int i;
  char tmp_str[MAXTREESTR];
  struct node *yo;

  treebuffer[0] = '\0';
        for (i=0;i<popsize;i++)
  {
    yo = makerandomtree(MAXDEPTH,FPR,PPR); 
    
    (*code2str_table[(int) yo->op])(yo,tmp_str);

    strcat(treebuffer,  tmp_str   );
    strcat(treebuffer,  delims  );
    free_node(yo);
  }  
  return 0;

}


int migrants2buffer(struct node *trees[],double *score_vals, double proportion, int low, int hi ,int popsize, double* ffs,double* result, double *obs, long* I, int ffs0, int ffs1, int obs0, int obs1, double aspect, double* S_init_array,int ts_factor, int startscore_i)
{

  int * bound_low;
  int *  bound_hi;

  int i,k;
  char tmp_str[MAXTREESTR];
  int migrants;

#ifdef CHECK_MIG_IC_STABILITY
  double S_init_array2[SPACEDIM];
#endif

#if defined (CHECK_MIG_STABILITY) || defined (CHECK_MIG_IC_STABILITY)
  double error;
#endif

  bound_low = mapsector(low,popsize);
  bound_hi = mapsector(hi,popsize);

  migrants = ((int) (proportion*((double) (bound_low[1]-bound_low[0] + bound_hi[1]-bound_hi[0]  )     )  ));

/*  fprintf(stderr,"low: %d, hi: %d, migrants: %d, popsize: %d \n",low,hi,migrants,popsize);
  fprintf(stderr,"is_low: %d, is_hi: %d, ie_low: %d, ie_hi: %d \n",is_low,is_hi, ie_low, ie_hi  );
*/

  treebuffer[0] = '\0';
  for (i=0;i<migrants;i++)
  {
    if (i%2==0)
    {
      k = tournament_scsize(trees,score_vals,bound_low[0],bound_low[1]-bound_low[0],MIGTOUR);
    }
    else
    {
      k = tournament_scsize(trees,score_vals,bound_hi[0],bound_hi[1]-bound_hi[0],MIGTOUR);
    }

#ifdef CHECK_MIG_STABILITY
    /* check time step stability to export only stable trees */
    error = get_score(ffs, result, obs, I, trees[k], ffs0, ffs1, obs0,obs1,aspect, S_init_array, ts_factor*4, startscore_i);
    if (fabs(score_vals[k]-error)>5.0)
    {
      score_vals[k] = 1e19;
    }
    else
    {
#endif
#ifdef CHECK_MIG_IC_STABILITY
      handle_S_init(4000.0, obs, obs1, S_init_array2); /* stability test w.r.t. i.c.*/
      error = get_score(ffs, result, obs, I, trees[k], ffs0, ffs1, obs0,obs1,aspect, S_init_array2, ts_factor, startscore_i);
 /*     fprintf(stderr,"TESTING: %g vs %g. S_init0: %g, S_init1: %g ",error,score_vals[k],S_init_array2[0],S_init_array2[1]); */

      if (fabs(score_vals[k]-error)>10.0)
      {
        score_vals[k] = (fabs(score_vals[k]-error)*score_vals[k])/5.0; 
      }
#endif

      /* add strings containing score and tree in bracket notation to treebuffer string */
      /* add score output_score */
      sprintf(tmp_str,"%g",score_vals[k]);
      strcat(treebuffer,  tmp_str );
      strcat(treebuffer,  delims );

      /* add tree */
      (*code2str_table[(int) trees[k]->op])(trees[k],tmp_str); /* convert tree to string */
      strcat(treebuffer,  tmp_str );
      strcat(treebuffer,  delims  );  
#ifdef CHECK_MIG_STABILITY
    }
#endif

  }  

  free(bound_low);
  free(bound_hi);
  return 0;
}


int readmigrants(const char *filepath, struct node *trees[], double *score_vals, int low, int hi , int popsize)
{
/* fill an array of tree structure pointers by translating text tree representations to structures  

   Produces warning to stderr when file reading error occurs.

   This function sometimes fails due to file corruption. At the moment, this can lead to corrupted trees entering the population, possibly with P/ C nodes with parents other than V nodes. floating point exceptions have been observed in conjunction with file reading errors. Possible fixes include a tree vetting function.

   opposite of trees2buffer
*/



  int * bound_low;
  int * bound_hi;

  char tmp_str[MAXTREESTR];
  int i, i_tour;
  double tmp_score;
  
  bound_low = mapsector(low,popsize);
  bound_hi = mapsector(hi,popsize);

  i=0;
  fp = fopen(filepath,"r");
  if (fp != NULL)
  {
    while((fgets (tmp_str,MAXTREESTR,fp)!=NULL ) && (i<MAXTREES)  )
    {
      strtok(tmp_str,"\n");   
      tmp_score = atof(tmp_str);

      if (fgets (tmp_str,MAXTREESTR,fp) != NULL) /* there should be a tree on the line after a score*/
      {

        if (( tmp_str[0] == '\n' ) || ( tmp_str[0] == '\0' ))
        {
          fprintf(stderr,"Error reading mig file: empty line \n");
          i=-1; /* abort and trigger failure */
          break;
        }
     
        if (check_brackets(tmp_str) == -1)
        {
          fprintf(stderr,"Error reading mig file: bracket tree integrity check failed for %s \n",tmp_str);
          i=-1; /* abort and trigger failure */
          break;
        }

        strtok(tmp_str,"\n");   

        if (i%2==0)
        {
          i_tour = reverse_tour_size(trees,score_vals,bound_low[0], bound_low[1]-bound_low[0], TOUR);
        }
        else
        {
          i_tour = reverse_tour_size(trees,score_vals,bound_hi[0], bound_hi[1]-bound_hi[0], TOUR);
        }

        free_node(trees[i_tour]); /* Note: free_node checks for NULL */

        score_vals[i_tour] = tmp_score;
        trees[i_tour] = str2node(tmp_str);
        if (check_node(trees[i_tour]) == 1) /* tree integrity check */
        {
          i = -1; /* abort and trigger failure */
          break;
        }
        i++;
      }
      else
      { /* if there's no tree after the line with a score, trigger failure */
        fprintf(stderr,"Error in reading migrant file %s (wrong format?). \n",filepath);   
        i = -1; /* abort and trigger failure */
        break;
      }
    }

    fclose(fp);
    free(bound_low);
    free(bound_hi);
    return i;
  }
  free(bound_low);
  free(bound_hi);
  return(-2); /* can't open file, not a big problem */
}


int readmigrants_buffered(const char *filepath, struct node *trees[], double *score_vals, int low, int hi , int popsize)
{
/* fill an array of tree structure pointers by translating text tree representations to structures  

   Produces warning to stderr when file reading error occurs.

   This function sometimes fails due to file corruption. At the moment, this can lead to corrupted trees entering the population, possibly with P/ C nodes with parents other than V nodes. floating point exceptions have been observed in conjunction with file reading errors. Possible fixes include a tree vetting function.

   opposite of trees2buffer
*/

  struct node *nodes_buffer[MAXTREES+5]; /* can we allocate a smaller amount?? */
  double scores[MAXTREES+5];

  int * bound_low;
  int * bound_hi;

  char tmp_str[MAXTREESTR];
  int i, num_trees;
  int i_tour=0;
  double tmp_score;

  
  bound_low = mapsector(low,popsize);
  bound_hi = mapsector(hi,popsize);

  i=0;
  fp = fopen(filepath,"r");
  if (fp != NULL)
  {
    while((fgets (tmp_str,MAXTREESTR,fp)!=NULL ) && (i<MAXTREES)  )
    {
      strtok(tmp_str,"\n");   
      tmp_score = atof(tmp_str);

      if (fgets (tmp_str,MAXTREESTR,fp) != NULL)
      {
        strtok(tmp_str,"\n");   

        scores[i] = tmp_score;
        nodes_buffer[i] = str2node(tmp_str);
        if (check_node(nodes_buffer[i]) == 1) /* tree integrity check */
        {
          i = -1; /* abort and trigger failure */
          break;
        }

        i++;
      }
      else
      {
        fprintf(stderr,"Error in reading migrant file %s (wrong format?). \n",filepath);   
        i=-1; /* abort and trigger failure */
        break;
      }
    }

    if (i>-1)
    { /* tree reading succesful, now introduce buffer trees into population via tournament */
      num_trees = i; /* previous while loop read this many trees */

      fprintf(stderr,"yo: %d \n",num_trees);
     
      for (i=0;i<num_trees;i++)
      {     

        if (i%2==0)
        {
          i_tour = reverse_tour_size(trees,score_vals,bound_low[0], bound_low[1]-bound_low[0], TOUR);
        }
        else
        {
          i_tour = reverse_tour_size(trees,score_vals,bound_hi[0], bound_hi[1]-bound_hi[0], TOUR);
        }

        fprintf(stderr,"i_tour %d \n",i_tour);
        free_node(trees[i_tour]); /* remove old tree (note: free_node checks for NULL) */

        score_vals[i_tour] = scores[i];
        trees[i_tour] = nodes_buffer[i];
        fprintf(stderr,"after\n");

      }
    }

    fclose(fp);
    free(bound_low);
    free(bound_hi);
    return i;
  }

  free(bound_low);
  free(bound_hi);
  return(-2);
}



int trees2buffer(struct node *trees[],int popsize)
{
/* fill the buffer with text representations of tree structures from an array 

   opposite of buffer2trees
*/
  int i;
  char tmp_str[MAXTREESTR];

  treebuffer[0] = '\0';
  for (i=0;i<popsize;i++)
  {
    (*code2str_table[(int) trees[i]->op])(trees[i],tmp_str);

    strcat(treebuffer,  tmp_str );
    strcat(treebuffer,  delims  );  
  }  
  return 0;
}

int buffer2trees(struct node *trees[])
{
/* fill an array of tree structure pointers by translating text tree representations to structures  

   opposite of trees2buffer
*/

  int i=0;
  char * pch;
  char tmp_str[MAXTREESTR];
  
  pch = strtok(treebuffer,delims);
  while ( pch != NULL)
  {
    strcpy(tmp_str,pch);
    trees[i] = str2node(tmp_str);  
    pch = strtok(NULL, delims);
    i++;
  };
  return 0;
}

void store_data(const char *filepath,const char *data, char *mode)
{
/* add return value in future */

  fp = fopen(filepath,mode);

  if (fp != NULL)
  {
    fprintf(fp,"%s",data);
    fclose(fp);
  }
}

int check_brackets(char *tmp_str)
{
  int i=0;
  int b=0;
  int t=0;

  while (tmp_str[i] != '\0')
  {
    if (tmp_str[i] == '(')
    {
      b++;
      t++;
    }
    else if (tmp_str[i] == ')')
    {
      b--;
      t++;
    }
    i++;
  }
  if ((t==0) || (b!=0)) 
  {
    return -1;
  }
  return 0;   
}

int check_ascii(char *tmp_str)
{
  int i=0;


  while ((tmp_str[i] != '\0') && (tmp_str[i] != '\n'))
  {
    if ((tmp_str[i] < '(') || (tmp_str[i] > 'z'))
    {
      return -1;
    }
    i++;
  }
  return 0;   
}

int read_data(const char *filepath,struct node *trees[],int maxtrees)
{
  char tmp_str[MAXTREESTR];
  int i=0;

  fp = fopen(filepath,"r");
  if (fp != NULL)
  {
    while((fgets (tmp_str,MAXTREESTR,fp)!=NULL ) && (i<maxtrees)  )
    {
      if (( tmp_str[0] == '\n' ) || ( tmp_str[0] == '\0' ))
      {
        fprintf(stderr,"Error reading pop: empty line \n");
        i=-1; /* abort and trigger failure */
        break;
      }
     
      if (check_brackets(tmp_str) == -1)
      {
        fprintf(stderr,"Error reading pop: bracket tree integrity check failed for %s \n",tmp_str);
        i=-1; /* abort and trigger failure */
        break;
      }

      strtok(tmp_str,"\n");
      trees[i] = str2node(tmp_str);
      if (check_node(trees[i]) == 1)
      {
        fprintf(stderr,"Error reading pop: tree integrity check failed.");
        i=-1; /* abort and trigger failure */
        break;
      }

      i++;
    }

    fclose(fp);
    return i;
  }
  return(-2); /* signal file failure */
}

double score_fun_mag(int startscore_i,int obs0,int obs1,double *obs, long* I, double* result, double magnification)
{
  int j;
  double error, tmp_error0;

  double upper_bound, lower_bound;

  lower_bound = -1.0;

  error = 0;
  for (j=startscore_i;j<obs0;j++)
  {      

    if (j-startscore_i < 220)
    {
      upper_bound = 0.1;
    }
    else
    {
      upper_bound = 0.5;
    }

    tmp_error0 = obs[j*obs1] - result[I[j]*SPACEDIM];

    if ( ( (obs[j*obs1] > upper_bound) && (result[I[j]]<upper_bound) ) ||  ( (obs[j*obs1] < lower_bound) && (result[I[j]]>lower_bound) )  ){
      tmp_error0 *= magnification; 
    }

    error += tmp_error0*tmp_error0;
  }
  
  return error;
}

double score_fun_basic(int startscore_i,int obs0,int obs1,double *obs, long* I, double* result)
{
  int j,k;
  double error, tmp_error0;
 
  double weights[10];

  weights[0] = 0.6;
  weights[1] = 0.4; 
  weights[2] = 0;
  weights[3] = 0;
  weights[4] = 0;

  error = 0;

  for (k=0;k<obs1;k++)  /* number of obs columns must not exceed SPACEDIM! */
  {
    for (j=startscore_i;j<obs0;j++)
    {      
 /*   fprintf(stderr,"%g # ",obs[j*obs1+1]);  */

      tmp_error0 = obs[j*obs1+k] - result[I[j]*SPACEDIM+k];
      error += tmp_error0*tmp_error0*weights[k];
    }
  }

  return 1000.0*error/(obs1*obs0);
}

double score_fun_pow4(int startscore_i,int obs0,int obs1,double *obs, long* I, double* result)
{
  int j,k;
  double error, tmp_error0;

  error = 0;

  for (k=0;k<obs1;k++)  /* number of obs columns must not exceed SPACEDIM! */
  {
    for (j=startscore_i;j<obs0;j++)
    {      
 /*   fprintf(stderr,"%g # ",obs[j*obs1+1]);  */

      tmp_error0 = obs[j*obs1+k] - result[I[j]*SPACEDIM+k];
      error += pow(1.25*tmp_error0,4);
    }
  }

  return 1000.0*error/(obs1*obs0);
}


double score_fun_aspect(int startscore_i,int obs0,int obs1,double *obs, long* I, double* result, double aspect)
{
  int j;
  double error, tmp_error0, tmp_error1;

  error = 0;
  for (j=startscore_i;j<obs0;j++)
  {      

    tmp_error0 = obs[j*obs1] - result[I[j]*SPACEDIM];
    tmp_error1 = (obs[j*obs1+1]-(result[I[j]*SPACEDIM]-result[I[j-1]*SPACEDIM]))*aspect;

    error += tmp_error0*tmp_error0 + tmp_error1*tmp_error1;
  }

  return error;
}


double get_score(double* ffs,double* result, double *obs, long* I, struct node *newtree, int ffs0, int ffs1, int obs0, int obs1, double aspect, double* S_init,int ts_factor, int startscore_i)
{
  double error;

/*  response_loop(ffs,result, newtree, ffs0, ffs1,S_init, ts_factor);  */

  if (check_node(newtree) == 1)
  {
    return 1e19;
  }

  RK4(ffs,result, newtree, ffs0, ffs1,S_init, ts_factor);

/*  fprintf(stderr,"yoyo %g\n",result[0]);
*/

  if (result[0] <1e18) 
  {

    if (aspect > 0)
    {
      error = score_fun_aspect( startscore_i, obs0, obs1, obs, I, result, aspect);
    }
    else
    {
      error = score_fun_basic( startscore_i, obs0, obs1, obs, I, result);
    }
/*    fprintf(stderr,"%f \n", sqrt(error));
*/
    return error;

  }
  else
  {
    return 1e19;
  }
  return 1e19;
}

int mutconsts(double *consts[],int size)
{
  int k, k2;
  double tmpval=0;
 
  if (size>0)
  {
    tmpval=rand();
    if ( tmpval<0.1*MUTATIONRATE*RAND_MAX )
    {
      k2=rand()%size;

      if (rand()%2 == 0)
      {
        *consts[k2] =   randomdouble();
      }
      else
      {
        *consts[k2] =   randomdouble()*(*consts[k2]);
      }


      for (k=0;(k<size) && (k != k2);k++)
      {
        
        *consts[k] =   (1.0+((double) (2*rand() -RAND_MAX) /RAND_MAX)*
            RNDGRAIN)*(*consts[k]);
      }  
    }
  }
return 0;
}


double hillclimb(double* ffs, double* result,double*  obs,long*  I,struct node *newtree, int ffs0, int ffs1, int obs0, int obs1, double *consts[] , int size, double error, double aspect, double* S_init, int ts_factor, int startscore_i)
{
  int k, km, oldkm, searches;
  double minerror;

  double delta[MAXCONSTS];
  double oldpos[MAXCONSTS];

  searches=0;

  if (size>0)
  {

    if (size>MAXHILL)
    {
      size = MAXHILL;
    }

    oldkm=-1;       
    for (k=0;k<size;k++)
    {
/*      delta[k] = (( (double) rand()  )/RAND_MAX)*pow(10,rand()%5-2); */

      delta[k] = (*consts[k])*( 1.0 + (( (double) rand() )/RAND_MAX)*0.5)*RNDGRAIN;
    }

    while(searches<MAXSEARCH)
    {
      km=0;
      minerror=error;
      for (k=0;k<size;k++)
      {
        oldpos[k] = *consts[k];  
        if (k != -oldkm)
        {
          (*consts[k]) = oldpos[k]+delta[k];
          error = get_score(ffs, result, obs, I, newtree, ffs0, ffs1, obs0, obs1,aspect, S_init, ts_factor, startscore_i);
          if (error < minerror-HILLTHRESH)
          {
            km=k;
            minerror=error;
          }
          (*consts[k]) = oldpos[k];
        }  
        else if (k != oldkm)
        {      
          (*consts[k]) = oldpos[k]-delta[k];
          error = get_score(ffs, result, obs, I, newtree, ffs0, ffs1, obs0, obs1,aspect, S_init, ts_factor, startscore_i);
          if (error < minerror-HILLTHRESH)
          {
            km=-k;
            minerror=error;
          }
                  (*consts[k]) = oldpos[k];
        }
      }
      if (km == 0)
      {
        break;
      }
      else
      {
        if (km>0)
        {
          (*consts[km]) = oldpos[km]+delta[km];
        }
        else
        {
          (*consts[-km]) = oldpos[-km]-delta[-km];
        }
        error = minerror;
        oldkm = km;
      }
      searches++;
    }  
  }

  nhillsearches += searches;
  return error;
}


int sectortour(double *score_vals,int popsize,int tour, int sector, int adjacent[])
{
/* tournament selection by sector, with process-internal migration. returns index of chosen tree
*/

  int * bound;
  int kmin=0;

  if (rand()<2*MIGRATE*RAND_MAX)
  {
    if (rand()%2==0)
    {
      bound = mapsector(adjacent[sector],popsize);
      kmin = tournament(score_vals, bound[0],bound[1]-bound[0],tour);
    }
    else
    {
      bound = mapsector(adjacent[sector+4],popsize);
      kmin = tournament(score_vals, bound[0],bound[1]-bound[0],tour);
    }
  }
  else
  {   
    bound = mapsector(sector,popsize);
    kmin = tournament(score_vals, bound[0],bound[1]-bound[0],tour);      

  }
 
  free(bound);

  return kmin;
}


int tournament_size(struct node *trees[], int is,int popsize,int tour)
{
/* tournament selection on size alone. returns index of chosen tree
*/
  int i,k,leafcount;
  int kmin=0;
  double M=1e20; /* minimum leafcount */

  int *current = make_int(0);

  for (i=0;i<tour;i++)
  {
    k = is + rand()%popsize; /* choose population wide within segment is to is+popsize */

    leafcount = conparcount(trees[k], current);
    *current = 0; /* always reset after use */
    if ( leafcount < M )
    {
      M = leafcount;
      kmin = k;
    }
  }
  free(current);

  return kmin;
}



int tournament_scsize(struct node *trees[],double *score_vals, int is,int popsize,int tour)
{
/* tournament selection based on score and tree size
*/
  int i,k,l,leafcount;
  int kmin=0;
  double M=1e20;  /* minimum value of leafcount */
  double sizescore; /* amalgemation of size and score */

  int *current = make_int(0);

  for (i=0;i<tour;i++)
  {
    l=0;
    do
    { /* look for a tree with leaves between 2 and MAXLEAF and try no more than 10 times */
      k = is + rand()%popsize; /* choose population wide in segment is to is+popsize*/
      leafcount = conparcount(trees[k], current);
      *current = 0;  /* always reset leafcount current afterwards */
      l++;
    }while( ( (leafcount < 2) || (leafcount > MAXLEAF) ) && (l<10) );

    if ( ( sizescore = ( score_vals[k] + 2*PARSIMONY*( (double) leafcount) )  ) < M)  
    { /* find the minimum value over the tournament */
      M = sizescore;
      kmin = k;
    }
/*    fprintf(stderr,"score: %f, parsimony: %f \n",score_vals[k],PARSIMONY*( (double) leafcount) ); 
*/

  }
  free(current);

  return kmin;
}

int reverse_tour_size(struct node *trees[],double *score_vals, int is,int popsize,int tour)
{
/* tournament selection on size and score, seeking maximum values. Used for tree elimination in migration.
*/
  int i,k,leafcount;
  int kmax=0;
  double M=0; /* maximum sizescore value */
  double sizescore; /* combination of size and score */

  int *current = make_int(0);

  for (i=0;i<tour;i++)
  {
    k = is + rand()%popsize; /* choose population wide inside segment is to is+ popsize */

    leafcount = conparcount(trees[k], current); /* calculate number of consts and pars (leaves) */
    *current = 0;      /* always reset leafcount current afterwards */

/*    if ( ( sizescore = score_vals[k]*( (double) leafcount) ) > M)   */
    if ( ( sizescore = ( score_vals[k] + 2*PARSIMONY*( (double) leafcount) )  ) > M) 
    { /* calculate maximum sizescore */
      M = sizescore;
      kmax = k;
    }
  }
  free(current);

  return kmax;
}





int tournament(double *score_vals, int is,int popsize,int tour)
{
/* tournament selection
*/
  int i,k;
  int kmin=0;
  double M=1e20;

  for (i=0;i<tour;i++)
  {
    k = is + rand()%popsize; /* choose population wide */

    if (score_vals[k] < M)
    {
      M = score_vals[k];
      kmin = k;
    }
  }
  return kmin;
}

int reversetour(double *score_vals,int is,int popsize,int tour)
{
/* inverse tournament selection
*/
  int i,k;
  int kmax=0;
  double M=0;

  for (i=0;i<tour;i++)
  {
    k = is + rand()%popsize;

    if (score_vals[k] > M)
    {
      M = score_vals[k];
      kmax = k;
    }
  }
  return kmax;
}

int * mapsector(int sector,int popsize)
{

  
  int * bound = malloc(sizeof(int)*2);
  bound[0] = sector*popsize/4;
  bound[1] = (sector+1)*popsize/4;

  return bound;
}


void  init_fixedforc(double *result, int result0, double *model,int model0)
{
  int j,k;
  int t0=-2500000;
  int dt=50;
 

  fixedforc = (double*)calloc(result0,sizeof(double));
  fixedforc2 = (double*)calloc(result0,sizeof(double));
  emph = (double*)calloc(result0,sizeof(double));
 
  for (j=0;j<result0;j++)
  {
    for (k=0;k<model0;k++)
    {
      fixedforc[j] += model[k]*pow(t0 + j*dt, model0-k-1);
    }

    for (k=0;k<model0-1;k++)
    {
      fixedforc2[j] += (model0-1-k)*model[k]*pow(t0 + j*dt, model0-k-2);
    }

    emph[j] = result[j*SPACEDIM];

/*    fixedforc[j] =result[j];
    fixedforc2[j] = (1.70188-13)*(t0 + j*dt) - 2.932128e-07;
*/
  }
}


void random_init(int popsize, struct node *trees1[], double* old_score_vals, double* ffs,double* result, double *obs, long* I, int ffs0, int ffs1, int obs0, int obs1, double aspect, double* S_init, int ts_factor, int startscore_i)
{
  int tries=0;
  int i=0;
  double error;
  char tmp_str[MAXTREESTR];

  struct node *newtree;

  while (i<popsize)
  {        
       
    node_count = 0;  /* node_count is a global var */

    if (tries<-2) /* avenue to introduce specific trees for testing */
    {
      sprintf(tmp_str,"((T,(p1,p2)M)V,(-5.69596e-05,-5.69596e-05)V)M");
      newtree = str2node(tmp_str);
    }
    else
    {
      newtree = makerandomtree(MAXDEPTH,FPR,PPR);
    }
    tries++;

    if (node_count > INITNODES)
    {
      /* score calculation for random init */
      error = get_score(ffs, result, obs, I, newtree, ffs0, ffs1, obs0,obs1,aspect, S_init, ts_factor, startscore_i);   
      if (error <1e18) 
      {
        trees1[i] = newtree;
        old_score_vals[i] = error;
        i++;
      }
      else
      {
        free_node(newtree);
      }
    }
    else
    {
      free_node(newtree);
    }
  }
}

void handle_S_init(double S_init, double* obs, int obs1, double* S_init_array)
{ /* set initial value array S_init_array for S array, based on signal S_init 

     S_init_array must be initialized to values before calling this function (e.g. [0,0] )

     Set S_init in config.ini to following values:
     S_init < -10 indicates that initial value of obs must be used to initialize all state variables
     S_init > 100 that all state variables init must have a random component around initial obs value
     S_init > 10 that all state variables init init must be entirely random

*/

  int i;

  if (obs1 > SPACEDIM)
  {
    obs1 = SPACEDIM;  
  }
  else if (obs1 < SPACEDIM)
  {
    for (i=obs1;i<SPACEDIM;i++)
    {
      S_init_array[i] = 0.0;
    }    
  }


  if (S_init < -10.0)
  { 
    /* remember time is removed from obs in python code! hence obs[0]*/
/*    S_init = obs[0];  */

    for (i=0;i<obs1;i++)
    {
      S_init_array[i] = obs[i];
    }

  }
  else if (S_init > 100.0)
  {
    for (i=0;i<obs1;i++)
    {
      S_init_array[i] =  obs[i] + (S_init*1e-4)*(( (double) (rand() -0.5*RAND_MAX) )/RAND_MAX);

    }
  }
  else if (S_init > 10.0)
  {
    for (i=0;i<SPACEDIM;i++)
    {
      S_init_array[i] = 4*(( (double) (rand() -0.5*RAND_MAX) )/RAND_MAX);
    }
  }
  else
  {
    for (i=0;i<SPACEDIM;i++)
    {
      S_init_array[i] = S_init;
    }
  }

/*
*/ 

}

void c_single_tree(double* ffs,double* result,double* obs, long* I, char treestr[], int ffs0, int ffs1, int result0,int obs0,int obs1, double *model, int model0, int ts_factor, int startscore_i, double S_init)
{
  int i;
  char tmp_str[MAXTREESTR];

  struct node *newtree;
  double aspect;
  double error;
  double S_init_array[SPACEDIM];

  aspect = model[0];
  srand( getpid()+time(NULL)  );

  init_tables();
  init_reg(ffs1);

/*  sprintf(tmp_str,"(-7.75295e-07,E)V"); */

  strcpy(tmp_str,treestr); /* leave original in tact: string will be altered on interpretation */
  newtree = str2node(tmp_str);

  init_fixedforc(result,result0,model,model0);

  handle_S_init(S_init, obs, obs1, S_init_array);

  for (i=0;i<SPACEDIM;i++)
  {
   

    printf("Init value S%d: %g\n",i,S_init_array[i]);
  }

  error = get_score(ffs, result, obs, I, newtree, ffs0, ffs1, obs0,obs1,aspect, S_init_array, ts_factor, startscore_i);
  model[1] = error;
/*  printf("tree score: %g with S_init: %g, numpar: %d \n", error,S_init, numpar);
*/
}



void c_nextgen(double* ffs,double* result, double* old_score_vals, double* score_vals,double* obs, long* I, char treefile[], char treeoutfile[], int ffs0, int ffs1, int result0,int my_number,int qsubs,int obs0,int obs1,int runlen,int popsize, double *model, int model0, int ts_factor, int startscore_i, double S_init)
{

/*  fprintf(stderr,"Aspect: %f\n",model[0]);
*/

  struct node *trees1[MAXTREES+5]; /* the old trees */
  struct node *trees2[MAXTREES+5]; /* the new trees */

/*
  struct node *(*previousgen)[MAXTREES+5] = &trees1;
  struct node *(*currentgen)[MAXTREES+5] = &trees2;
  struct node *(*tmp_pointer)[]; 
*/

  int *current = make_int(0);

  int *nconsts;
  int size=0;
  int sector;

  int i,j,run, totsize,qrows,qcols,my_number_x,my_number_y, i_min, leafcount, trees_reused;
  int error_flag=0;

  double S_init_array[SPACEDIM];
  double aspect;
  double avghill, avgsize, avgtreesreused;
  
  double mutationrate = MUTATIONRATE; 
  double error, minerror, toterror,varerror;
  double *consts[MAXCONSTS];

  char tmp_str[MAXTREESTR]; /* used only once */
  char line_str[MAXTREESTR];

  char * infiles[4];
  char * outfiles[4];

  char outfile1[FNAMESIZE];
  char outfile2[FNAMESIZE];
  char outfile3[FNAMESIZE];
  char outfile4[FNAMESIZE];

  char infile1[FNAMESIZE];
  char infile2[FNAMESIZE];
  char infile3[FNAMESIZE];
  char infile4[FNAMESIZE];

  int lowmig[4] = {0,1,0,2};
  int himig[4] = {2,3,1,3};

  int adjacent[8] = {1,0,0,1,2,3,3,2};

  char repfile[FNAMESIZE];
  char elitefile[FNAMESIZE];

  struct node *newtree;
  struct node *tmptree1;
  struct node *tmptree2;
  struct node *tmptree3;

  int i_tree1, i_tree2, i_tree_old;
  double existing_score = -1; 

/*
  gettimeofday(&time,NULL);
  srand((unsigned) getpid()+(time.tv_sec*100) + (time.tv_usec/100)  );
*/
  once_flag = 0;

  aspect = model[0];

  srand( getpid()+time(NULL)  );

/*  printf("random seed: %d \n",getpid()+(time.tv_sec*100) + (time.tv_usec/100));
*/
  init_tables();
  init_reg(ffs1);

  printf("runlen: %d, my_number: %d, qsubs: %d, S_init: %g, numpar: %d \n",runlen,my_number, qsubs, S_init, numpar);


  init_fixedforc(result,result0,model,model0);

  handle_S_init(S_init, obs, obs1, S_init_array);

  qrows = (int) sqrt(qsubs);
  qcols = (int) qsubs/qrows;

  my_number_y = ((int) mod(my_number/qcols,qrows));
  my_number_x = ((int) mod(my_number, qcols) );

#ifdef HI_MUT_SPOT
  // Create a small spot with increased mutation rate
  if (  (qsubs > 16) &&   (my_number_y <= 1) && ( my_number_x <= 1)    ) 
  {
    mutationrate = 0.5;
  }
#endif


  for (j=0; j < popsize; j++)
  {
    trees1[j]=NULL;
  }

  nconsts = make_int(0);

  if (popsize > MAXTREES)
  {
    fprintf(stderr,"Error: popsize %d exceeds MAXTREES %d.\n",popsize,MAXTREES);
  }

  int ntrees = read_data(treefile, trees1,MAXTREES);
  if (ntrees <= -1)
  {
    if (ntrees == -1)
    { 
      fprintf(stderr,"Bad pop file detected, populating from random init.\n");
    }

    printf("No data file or bad file. Random init \n");
/*    mutationrate = 4.0*mutationrate;   */
    
    random_init(popsize, trees1, old_score_vals, ffs, result, obs, I, ffs0, ffs1, obs0, obs1, aspect, S_init_array, ts_factor, startscore_i);
    ntrees = popsize;
  }



/* the trees in the ordered pointer array is used to create a new sequence of trees (now not ordered) */
/* create file name strings for migrant files */
/* in from-to format */

  sprintf(infile1,"migr%d-%dxp",mod(my_number_x-1,qcols),mod(my_number_y,qrows)  );
  sprintf(infile2,"migr%d-%dxm",mod(my_number_x+1,qcols),mod(my_number_y,qrows) );

  sprintf(infile3,"migr%d-%dyp",mod(my_number_x,qcols),mod(my_number_y-1,qrows) );
  sprintf(infile4,"migr%d-%dym",mod(my_number_x,qcols),mod(my_number_y+1,qrows)  );

  sprintf(outfile1,"migr%d-%dxm",my_number_x,my_number_y );
  sprintf(outfile2,"migr%d-%dxp",my_number_x,my_number_y );

  sprintf(outfile3,"migr%d-%dym",my_number_x,my_number_y );
  sprintf(outfile4,"migr%d-%dyp",my_number_x,my_number_y );

  sprintf(repfile,"report%d",my_number);
  sprintf(elitefile,"elite%d",my_number);

  infiles[0] = infile1;
  infiles[1] = infile2;
  infiles[2] = infile3;
  infiles[3] = infile4;

  outfiles[0] = outfile1;
  outfiles[1] = outfile2;
  outfiles[2] = outfile3;
  outfiles[3] = outfile4;


  for (run=0;run<runlen;run++)
  {

    for (j=0;j<4;j++) /* read migrant files */
    { 
      error_flag = readmigrants(infiles[j],trees1,old_score_vals, lowmig[j], himig[j] , popsize); /* j arg is sector */
      if (error_flag == -1)
      {
        break;
      }
    }
    if (error_flag == -1) 
    { 
      /* One bad file will cause pop to be dropped. Using function readmigrants_buffered instead (above) allows keeping the pop. */
      fprintf(stderr, "Error in migrant io detected, replacing entire population.\n");
      random_init(popsize, trees1, old_score_vals, ffs, result, obs, I, ffs0, ffs1, obs0, obs1, aspect, S_init_array, ts_factor, startscore_i);
    }

    totsize=0;
    minerror=1e20;
    i_min=0;
    toterror=0;
    trees_reused=0;

    i=0;
    while (i<popsize) /* fill trees2 and score_vals  */
    {
      node_count = 0; /* count how many nodes are produced for newtree for cost function*/

      existing_score = -1.0; 
      if ( rand() > PNEW*RAND_MAX)   
      {
        sector = 4*i/popsize;

        i_tree1 = sectortour(old_score_vals,ntrees, TOUR , sector,adjacent);
        i_tree2 = sectortour(old_score_vals,ntrees, TOUR , sector,adjacent);

        tmptree1 = trees1[i_tree1];
        tmptree2 = trees1[i_tree2];

    
        tmptree3 = unifcross(tmptree1,tmptree2, PROBSWAP,1);

        newtree = allmut(tmptree3, mutationrate); 
 //       fprintf(stderr,"%d",compare_trees(tmptree1,newtree));

        if (compare_trees(tmptree1,newtree) == 0)
        {
          existing_score = old_score_vals[i_tree1];
          trees_reused++;
        }
        else
        {
          i_tree_old = tree_in_old(newtree,trees1, popsize);
          if (i_tree_old > -1)
          {
            existing_score = old_score_vals[i_tree_old];
            trees_reused++;            
          } 
        }

        free_node(tmptree3);
      }
      else
      {
        newtree = makerandomtree(MAXDEPTH,FPR,PPR);
      }


      size = tree_consts(newtree, consts,nconsts,CONSTSCUTOFF); /* fill consts array */
      *nconsts=0;
/*      fprintf(stderr,"size=%d\n",size);  */

      if ( (size > -1) && (size < MAXSIZE) )
      {  
    /* perform optimization here */    

   /*     mutconsts(consts,size); */

        if (existing_score > MIN_SCORE)
        {
          error = existing_score;
        }
        else
        {
          error = get_score(ffs, result, obs, I, newtree, ffs0, ffs1, obs0,obs1,aspect, S_init_array, ts_factor, startscore_i);
        }
 
        if (error <1e18) 
        {


            /* note that hillclimb checks size>0 */
/*          error = hillclimb(ffs, result, obs, I, newtree, ffs0, ffs1, obs0, obs1, consts , 
                             size, error,aspect, S_init, ts_factor, startscore_i);
*/ 

/*        parsimony not included in error reporting to report and elite output files */

          leafcount = conparcount(newtree, current);
          *current = 0;


     /*     fprintf(stderr,"%f  %f %d \n",error,((double) PARSIMONY*((double) leafcount)) , leafcount );  */

          // fill score vals in new generation (the old generation has old_score_vals).
          score_vals[i] = error + ((double) PARSIMONY*((double) leafcount)); 


          // fill tree array for new generation
          trees2[i] = newtree; /* found an acceptable tree */
          totsize +=size;
          toterror += error;
          if (error < minerror)
          {
            minerror = error;
            i_min = i;
          }
          i++;
        }
        else
        {
          free_node(newtree);
        }
      }
      else
      {
        free_node(newtree);
      }
    }; /* end while loop on popsize */


  /*
    copy score_vals to old_score_vals and copy trees2 to trees1 (to be improved later). 
  */

    varerror=0;
    toterror = ((double) toterror/popsize);
    for (j=0; j < popsize; j++)
    {
      varerror += (score_vals[j] - toterror)*(score_vals[j] - toterror);
    }

    avghill = ((double) nhillsearches)/popsize;    
    nhillsearches = 0;

    avgsize = ((double) totsize)/popsize;
    
    avgtreesreused = ((double) trees_reused)/popsize;  

    for (j=0;j<4;j++) /* output migrants to adjacent processes to files */
    {
      migrants2buffer(trees2, score_vals, MIGRATE, lowmig[j], himig[j],popsize, ffs, result, obs, I, ffs0, ffs1, obs0,obs1,aspect, S_init_array, ts_factor, startscore_i);
      store_data(outfiles[j],treebuffer,"w"); /* write migrant tree files */
    }
  
    /* report stats */
    sprintf(line_str,"%g %g %g %g %g %g\n",minerror,toterror,varerror,avgsize,avghill, avgtreesreused);
    store_data(repfile,line_str,"a"); /* append stats to report file (filename repfile) */

    (*code2str_table[(int) trees2[i_min]->op])(trees2[i_min],tmp_str);

    sprintf(line_str,"%g %s\n",minerror, tmp_str);
    store_data(elitefile,line_str,"a"); /* append stats to report file */

    for (j=0; j < popsize; j++)
    {
      old_score_vals[j] = score_vals[j];
      free_node(trees1[j]);
      trees1[j] = trees2[j];
    }


  }  /* end of runlen for loop */

  trees2buffer(trees2,popsize);

  for (j=0; j < popsize; j++)
  {
    free_node(trees1[j]);
  }
  free(fixedforc);
  free(nconsts);
  store_data(treeoutfile,treebuffer,"w");
}

