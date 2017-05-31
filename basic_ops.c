#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <stdint.h>
#include <sys/mman.h>

#include <compile_options.h>
#include <states.h>
#include <node_structures.h>
#include <fields.h>
#include <basic_ops.h>
#include <my_ops.h>
#include <tree_io.h>
#include <model.h>
#include <c_nextgen.h>

#ifdef INTSTATES
int reg[MAXPAR];
int reg_bck[MAXPAR];
#else
double reg[MAXPAR];
double reg_bck[MAXPAR];
#endif

void replace_p(Node *node);



// integer leaf-related functions





// make_par_node can be used for par nodes in the int setting

Node *make_const_node_int(char op, int value)
{
  Node *node = init_node(op);

  node->children[0] = make_int(  value  );

  return node;
}

Node *makerandomconstnodeint(void)
{
#ifdef MAXINT
  return make_const_node_int(const_op_char,  (int)  rand()%MAXINT  );
#else
  return make_const_node_int(const_op_char,  (int)  rand()%10  );
#endif
}


int p2index(Node *tree)
{ // extract integer index from par node
  return *((int *) tree->children[0]);
}


//typedef double (*FP)();

//FP make_const_node_run(double value)
//{

// on 64 bit NX protected systems:
//  uint8_t *buf = mmap(NULL, 1000, PROT_EXEC | PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);


//  buf[0] = 0xb8;
//  double u32 = value;
//  memcpy(buf + 1, &u32, 8);
//  buf[5] = 0xc3;

//  double (*ptr)(void) = (double (*)(void)) buf;

//  return ptr;

//}




// op_table functions

double nopFun(Node *tree)
{
/* no operation function to initialize op_table */

  fprintf(stderr,"NOPFUN\n ");
  exit(-1);

  return 0;

};


double constFun(Node *tree)
{
/* Yield a constant leaf.  */
  return *((double *) tree->children[0]); 

};


double parFun(Node *tree)
{
/* Yield a parameter leaf.  */

  return reg[*((int *) tree->children[0]) ];
};


double parFun0(Node *tree)
{
/* Yield a parameter leaf.  */

  return reg[0];
};


double parFun1(Node *tree)
{
/* Yield a parameter leaf.  */

  return reg[1];
};

double parFun2(Node *tree)
{
/* Yield a parameter leaf.  */

  return reg[2];
};

double parFun3(Node *tree)
{
/* Yield a parameter leaf.  */

  return reg[3];
};

double parFun4(Node *tree)
{
/* Yield a parameter leaf.  */

  return reg[4];
};




double addFun(Node *tree)
{
  /* execute addition node: add the values of the two child nodes. */

  return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) + (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))) ;

};

double subFun(Node *tree)
{
  /* execute subtraction node: subtract the values of the two child nodes. */

  return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) - (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))) ;

};

double mulFun(Node *tree)
{
  /* execute multiplication node */

  return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) * (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))) ;

};

double divFun(Node *tree)
{
  /* execute division node: divide the values of the two child nodes. */

  return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) / (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))) ;

};

double isgreaterFun(Node *tree)
{

  if ((*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) > (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))))
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

double iseqFun(Node *tree)
{
  if ((*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) == (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))))
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

double ifFun(Node *tree)
{
  /* execute if-then branching node: execute child 1 if child 0 > 0, child 2 otherwise */

  if ((*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0])))>0)
    return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))) ;      
  else 
    return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[2]))->op)])(((Node *) ((Node *) tree->children[2]))) ;  
};


double markerFun(Node *tree)
{
  /* marker node, returns evaluation of single child unmodified */
  return (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0])));
};




double sqrtFun(Node *tree)
{
  /* execute square root node: take the square root of the child node. */
  return  sqrt( (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0])))   );
};


double tanhFun(Node *tree)
{
  /* execute square root node: take the square root of the child node. */
  return  rational_tanh( (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0])))  );
};

double stemFun(Node *tree)
{
  // for now it's a nop
  return 0;
}


double funFun(Node *tree)
{
#if NUMFUN>0
  int i;
  double result;
  Node *fun_node;
  Node *fun_arg;

  for (i=0;i<NUMFUN;i++)
  {
    fun_arg = tree->children[i+1];
    reg_bck[i] = reg[i];
    reg[i] = run_node(fun_arg);
  }

  for (i=NUMFUN;i<numpar;i++)
  {
    reg_bck[i] = reg[i];
    reg[i] = 0;
  }

  fun_node = tree->children[0];
  result = run_node(fun_node);

  for (i=0;i<numpar;i++)
  {
    reg[i] = reg_bck[i];
  }

  return result;
#else
  return 0;
#endif

};


// Integer run functions


int nopFunInt(Node *tree)
{
/* no operation function to initialize op_table */

  fprintf(stderr,"NOPFUNINT\n ");
  exit(-1);

  return 0;

};

int constFunInt(Node *tree)
{
/* Yield a constant leaf.  */
  return *((int *) tree->children[0]); 

};

int parFunInt(Node *tree)
{
/* Yield a parameter leaf.  */

  return reg[*((int *) tree->children[0]) ];
};

int addFunInt(Node *tree)
{
  /* execute addition node: add the values of the two child nodes. */

  return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) + (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))) ;

};

int subFunInt(Node *tree)
{
  /* execute subtraction node: subtract the values of the two child nodes. */

  return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) - (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))) ;

};

int mulFunInt(Node *tree)
{
  /* execute multiplication node */

  return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) * (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))) ;

};

int divFunInt(Node *tree)
{
  /* execute division node: divide the values of the two child nodes. */

  return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) / (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))) ;

};

int isgreaterFunInt(Node *tree)
{

  if ((*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) > (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))))
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int iseqFunInt(Node *tree)
{
  if ((*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0]))) == (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))))
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int ifFunInt(Node *tree)
{
  /* execute if-then branching node: execute child 1 if child 0 > 0, child 2 otherwise */

  if ((*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0])))>0)
    return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[1]))->op)])(((Node *) ((Node *) tree->children[1]))) ;      
  else 
    return  (*op_table[( (int)  (  (Node *) ((Node *) tree->children[2]))->op)])(((Node *) ((Node *) tree->children[2]))) ;  
};


int markerFunInt(Node *tree)
{
  /* marker node, returns evaluation of single child unmodified */
  return (*op_table[( (int)  (  (Node *) ((Node *) tree->children[0]))->op)])(((Node *) ((Node *) tree->children[0])));
};

int stemFunInt(Node *tree)
{
  // for now it's a nop
  return 0;
}












// node production functions

Node *init_node(char op)
{
  Node *node;
  node = talloc();
  node->op = op ;

  return node;
}


Node *make_par_node(char op, int i)
{
  Node *node = init_node(op);
  
  node->children[0] = make_int((int)  i );  

  return node;
}

Node *make_par_fun_node(char op, int i, double a)
{
  Node *node = init_node(op);
  
  node->children[0] = make_int((int)  i );  
  node->children[1] = make_double(  a  );

  return node;
}


Node *make_const_node(char op, double value)
{
  Node *node = init_node(op);

  node->children[0] = make_double(  value  );

  return node;
}



// random node functions

double randomdouble(void)
{
  /* Produce a random floating point number, with some biases */

  int i = rand()%3;
  int p;
  double t = DOUBLETRUNC;

  double value = (( (double) (2*rand() -RAND_MAX) )/RAND_MAX);

  value = trunc(t*value)/t;


  if (i<1)
  {
    p = rand()%4-2;
  }  
  else if (i<2)
  {
    p = rand()%6-3;
  }  
  else
  {
    p = rand()%14-7;
  }

  return value*pow(10,p);
}

Node *makerandomparnode(int max_par)
{
#if RESTRICTED_PARS > 0
  return make_par_node(par_op_char, (int)  numpar - RESTRICTED_PARS + rand()%RESTRICTED_PARS );
#else
  return make_par_node(par_op_char, (int)  rand()%max_par );
#endif

}


Node *makerandomparfunnode(int max_par)
{
  int step = 8;
  double a;
  int i_a = (rand()%step+1);

  a = ((double) i_a)/( (double) step );


  return make_par_fun_node('H', (int)  rand()%max_par,  a );
}


Node *makerandomconstnode(void)
{
  return make_const_node(const_op_char,  randomdouble()  );
}


Node *makerandomterminal_basic(int maxdepth,double fpr,double ppr, int max_par, Node* (*f_const)(void) )
{

  int rnd;

  rnd = rand();
  if ( rnd<ppr*RAND_MAX  )
  {
  /* make parameter node */
    return makerandomparnode(max_par);
  }
  else
  {
  /* make constant node */
    return (*f_const)();
  };
}



// free the leaf nodes

void free_param_node(Node *t)
{
    if (t->children[0] != NULL)
    {
      free(((Node *) t->children[0]));
      t->children[0] = NULL;
    }
}

void free_const_node(Node *t)
{
    if (t->children[0] != NULL)
    {
      free(((Node *) t->children[0]));
      t->children[0] = NULL;
    }
}


void free_param_fun_node(Node *t)
{

  if (t->children[0] != NULL)
  {
    free(((Node *) t->children[0]));
    t->children[0] = NULL;
  }

  if (t->children[1] != NULL)
  {
    free(t->children[1]);
    t->children[1] = NULL;
  }
  

}




// node copy functions

Node *copy_const_node(Node *tree)
{
  return make_const_node(tree->op, *((double *) tree->children[0]) );
}

Node *copy_par_node(Node *tree)
{
  return make_par_node(tree->op, *((int *) tree->children[0]) );
}

Node *copy_par_fun_node(Node *tree)
{
  return make_par_fun_node(tree->op,*((int *) tree->children[0]) , *((double *) tree->children[1]) );
}


Node *copy_non_leaf(Node *tree)
{
  Node *newnode;
  int i;
  Node *child;

  newnode=talloc();
  newnode->op = tree->op;

  for (i=0;i<arg_table[tree->op];i++) // no looping for childless nodes
  {
    child = ((Node *) tree->children[i]);
    newnode->children[i] = (*copy_table[child->op])(child);
  }; 

  return newnode;
}


void replace_p(Node *node)
{ // recursively replace p nodes with p0-5 run nodes

  int i;

  if (node->op == par_op_char)
  {
    node->op = 'a'+*((int *) node->children[0]); 

//    free(node->children[0]);
//    node->children[0] = NULL;

  }
  else
  {
    for (i=0;i<arg_table[node->op];i++)
    {
      replace_p((Node *) node->children[i]);
    }
  }
}

Node *make_run_node(Node *tree)
{
  Node * runtree = copy_node(tree);
  replace_p(runtree);

  return runtree;

}

Node *copy_node(Node *tree)
{
  /* Copy a node including its subtrees, and allocate memory using talloc for all nodes involved.

  */

  // jump on copy_table defined in my_ops.c, providing appropriate recursion.
  return (*copy_table[tree->op])(tree);

};


Node *copy_node_nojump(Node *tree)
{

  if (tree->op == const_op_char)
  {
    return make_const_node(tree->op, *((double *) tree->children[0]) );
  }
  else if (tree->op == par_op_char)
  {
    return make_par_node(tree->op, *((int *) tree->children[0]) );
  }
  else
  {
    Node *newnode;
    int i;
    Node *child;

    newnode=talloc();
    newnode->op = tree->op;

    for (i=0;i<arg_table[tree->op];i++) // no looping for childless nodes
    {
      child = ((Node *) tree->children[i]);
      newnode->children[i] = copy_node_nojump(child);
    }; 

    return newnode;
  }

}





int compare_const_nodes(char op, Node *t1, Node *t2)
{

  if (  ((Node *) t1)->op != op )
  {
    return 1;
  }

  if (t1->op != t2->op)
  {
    return 1;
  }

  if (  ((double *) t1->children[0]) != ((double *) t2->children[0])   )
  {
    return 1;
  }

  return 0;
}

int compare_par_nodes(char op, Node *t1, Node *t2)
{
  if (  ((Node *) t1)->op != op )
  {
    return 1;
  }

  if (t1->op != t2->op)
  {
    return 1;
  }

  if (  ((int *) t1->children[0]) != ((int *) t2->children[0])   )
  {
    return 1;
  }

  return 0;
}

int compare_par_fun_nodes(char op, Node *t1, Node *t2)
{
  if (  ((Node *) t1)->op != op )
  {
    return 1;
  }

  if (t1->op != t2->op)
  {
    return 1;
  }

  if (  ((int *) t1->children[0]) != ((int *) t2->children[0])   )
  {
    return 1;
  }

  if (  ((double *) t1->children[1]) != ((double *) t2->children[1])   )
  {
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

    define init_tables_* functions in my_ops module

*/

  init_tables_ops();
  init_tables_copy();
  init_io_tables(arg_table);

  
  return 0;
}



