#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <compile_options.h>
#include <states.h>
#include <node_structures.h>
#include <fields.h>
#include <basic_ops.h>
#include <my_ops.h>
#include <tree_io.h>
#include <model.h>
#include <c_nextgen.h>

// globals

#ifdef INTSTATES
int (*op_table[OPTABLE])(Node *tree);
#else
double (*op_table[OPTABLE])(Node *tree);
#endif

char* (*code2str_table[OPTABLE])(Node *tree, char trstr[]);
char* (*code2c_str_table[OPTABLE])(Node *tree, char trstr[]);

Node* (*copy_table[OPTABLE])(Node *tree);
Node* (*str2leaf_table[OPTABLE])(char text[], char left_bracket, char delim);


#ifdef DEBUG
char*instrsetnodes = "QEOTLASMIGEVFY"; /* for printing */
#endif

int arg_table[OPTABLE];
char treebuffer[TREEBUFFER];

int node_count;
void remove_op(Node *tree, char op, int first_encounter);
int force_params_node(Node *tree, int i_param);


// make random functions

Node *makerandomterminal(int maxdepth,double fpr,double ppr, int max_par)
{
#ifdef INTCONSTS
  return makerandomterminal_basic(maxdepth, fpr, ppr, max_par, makerandomconstnodeint);
#else
  return makerandomterminal_basic(maxdepth, fpr, ppr, max_par, makerandomconstnode);
#endif
}


void free_node(Node *t)
{

  if (t == NULL)
  {
    return;
  }
 
  int i, i_at;

  if (t->op == const_op_char)
  {
    free_const_node(t);
  }
  else if ( (t->op == par_op_char) || ( (t->op>='a' ) && (t->op<='a'+5 )  )  )
  {
    free_param_node(t);
  }
//  else if  ((t->op == 'H'))
//  {
//    free_param_fun_node(t);
//  }
  else
  {

    i_at = ( (int) t->op );

    if (i_at > OPTABLE)
    {
      fprintf(stderr,"Error: node op %d greater than OPTABLE %d.\n",i_at, OPTABLE);
      exit(-1);
    }

    for (i=0;i<arg_table[i_at];i++)
    {
      free_node(( (Node *) t->children[i]));
    }
  }

  if (t != NULL)
  {
    free(t);
  }
  t = NULL;


  node_count--;
}



int check_node(Node *t)
{
  Node *child;

  if (t == NULL)
  {
    return 1;
  }
 
  int i;

  if  (t->op == const_op_char)
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
  else if (t->op == par_op_char)
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

#ifdef DEBUG
    if (strchr(instrsetnodes, t->op ) == NULL)
    {
      fprintf(stderr,"Warning: illegal op \'%c\' encountered!  \n",t->op );
      return 1;
    }
#endif


    for (i=0;i<arg_table[(int) t->op];i++) /* zero-nodes will always return 0  */
    {
      child = ((Node *) ( (Node *) t->children[i]));
      if (child == NULL)
      {
        fprintf(stderr,"Warning: NULL child node encountered under %c at i=%d! Is spacedim=%d correct? Or check meta attribute options. \n",t->op ,i,SPACEDIM);
        return 1;
      }
      else
      {
        if ( check_node(child) == 1 )
        {
          return 1;
        }
      }
    }
    
    return 0;
  }

}


int tree_consts(Node *t, double *consts[], int *current, int maxconsts)
{
 
  int j=0;

  if (t->op ==const_op_char)  
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
//  else if ( (t->op ==const_op_char) || (t->op =='H') )
  else if  (t->op ==par_op_char) 
  {
    return 0;
  }
  else
  {
    while (j<arg_table[t->op])  
    {  
      if ( tree_consts( ((Node *) t->children[j]), consts,current,maxconsts) > -1 )     
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

void remove_op(Node *tree, char op, int first_encounter)
{
/* remove op by linking parent to grandchild.
   returns number of replacements made.

   so far only works on 1-child ops

  */

  int i;
  Node *child;
  Node *grandchild;

  for (i=0;i<arg_table[tree->op];i++) // no looping for childless nodes
  {
    child = tree->children[i];

    if (child->op == op)
    {
      if (first_encounter == 1)
      {
        remove_op(child, op, 0);
      }
      else
      {
        do
        { // there might be multiple nested nodes with op
          grandchild = child->children[0];
          tree->children[i] = grandchild; 
          free(child);
          child = grandchild;
        } while (child->op == op);

        remove_op(grandchild, op, 0);   
      }
    } 
    else
    {
      remove_op(child, op, first_encounter );
    }
  }
  return;
};



int tree_height(Node *tree, int depth)
{
/* Determine tree height (maximum path from leaf to root)

  */

  int i;
  Node *child;

  if (arg_table[tree->op] > 0)
  {
    depth++; // every time a node has children, depth is increased by 1
    int max_child_depth = 0;
    for (i=0;i<arg_table[tree->op];i++) 
    {
      child = tree->children[i];
      int depth_plus_child = tree_height(child, depth); // current depth

      if (depth_plus_child > max_child_depth)
      {
        max_child_depth = depth_plus_child;
      }
    }  
    return max_child_depth;
  }

  return depth;

};



int force_params_node(Node *tree, int i_param)
{
/*  alter existing tree to force all parameters with index greater than SPACEDIM to be set to i_param */

  int i_value;
  int i;

  if (tree->op == par_op_char)
  {
    i_value = *((int *) tree->children[0]);
    if (i_value >= SPACEDIM)
    {

      *((int *) tree->children[0]) = i_param;
    }
  }
  else
  { 
    for (i=0;i<arg_table[tree->op];i++) // no looping for childless nodes
    {
      force_params_node(( (Node *) tree->children[i]) , i_param);
    };
  };

  return 0;
};

int force_at_marker_node(Node *tree, char marker, int i_param)
{
/* Walk tree and look for marker op and use force_params_node for subtrees of those markers */

  int i;

  if (tree->op == marker)
  {
    force_params_node(tree , i_param);
  }
  else
  { 
    for (i=0;i<arg_table[tree->op];i++) // no looping for childless nodes
    {
      force_at_marker_node(( (Node *) tree->children[i]), marker , i_param);
    };
  };

  return 0;
}




int leaves_equal(Node *t1,Node *t2)
{
// return 0 if t1 and t2 are leaf nodes that are equal, 1 otherwise

  if  ( (is_terminal(t1) == 1) && (is_terminal(t2) == 1)  ) // is_terminal returns 1 if t is a terminal
  {
    if (compare_par_nodes(par_op_char,t1, t2) == 0)
    {
      return 0;
    }
//    else if (compare_par_fun_nodes('H', t1, t2) == 0)
//    {
//      return 0;
//    }
    else if (compare_const_nodes(const_op_char, t1, t2) == 0)
    {
      return 0;
    }
  }
  return 1;
}

void curtail_par_nodes(Node *tree, int max_num)
{
 
  int i;
  int *int_child;

  if (tree->op == par_op_char)
  {
    int_child = ((int *) tree->children[0]);
    if (*int_child >= max_num)
    {
      *int_child = (*int_child)%max_num;
    }
  }
  else
  {
    for (i=0;i<arg_table[tree->op];i++)
    {
      curtail_par_nodes(( (Node *) tree->children[i]), max_num);  
    }
  }
}


void user_functions(NodeScore ns)
{


#if INTPARS > 0
  force_params_node(ns.node, *ns.ip.data);
#endif

#if MARKER > 0
  force_at_marker_node(ns.node, 'Y' , MARKER);
#endif

#ifdef TANH_ONCE
  remove_op(ns.node, 'T', 1);
#endif

#if defined(STEM_NODES) && defined(CURTAILFORCING)
// avoid having forcing terms driving state variables other than S1
  
  int i;
  Node *child;

  for (i=1;i<arg_table[((int) ((Node *) ns.node)->op)];i++)
  {
    child = ns.node->children[i];
    curtail_par_nodes(child,SPACEDIM); // this limits forcing terms to first dimension only (reg[0])
  }
#endif
}


int init_tables_ops(void)
{
/* Initialize op_table
*/

  int i;

  // operations to be executed at runtime


#ifdef INTSTATES
  for (i=0;i<OPTABLE;i++)
  {
    op_table[i] = &nopFunInt;
  };

  op_table['A'] = &addFunInt;  
  op_table['S'] = &subFunInt;
  op_table['M'] = &mulFunInt;
  op_table['D'] = &divFunInt;
  op_table['I'] = &ifFunInt;
  op_table['G'] = &isgreaterFunInt;
  op_table['E'] = &iseqFunInt;
  op_table['V'] = &stemFunInt;
  op_table['Y'] = &markerFunInt;

  op_table[((int) const_op_char)] = &constFunInt;
  op_table[((int) par_op_char)] = &parFunInt;  

#else
  /* linking ops to recursive functions, superleafFun is terminal as it produces stored values */

  for (i=0;i<OPTABLE;i++)
  {
    op_table[i] = &nopFun;
  };

  op_table['A'] = &addFun;  
  op_table['S'] = &subFun;
  op_table['M'] = &mulFun;
  op_table['D'] = &divFun;
  op_table['I'] = &ifFun;
  op_table['G'] = &isgreaterFun;
  op_table['E'] = &iseqFun;
  op_table['Q'] = &sqrtFun;
  op_table['T'] = &tanhFun;
  op_table['V'] = &stemFun;

  op_table['Y'] = &markerFun;

  op_table[((int) const_op_char)] = &constFun;
  op_table[((int) par_op_char)] = &parFun;  

#if FORCING_TERMS > 0
// fast parFuns
  op_table['a'] = &parFun0;  
  op_table['b'] = &parFun1;  
  op_table['c'] = &parFun2;  
  op_table['d'] = &parFun3;  
  op_table['e'] = &parFun4;  
#endif

#endif

  // argument table
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
    arg_table[((int) instrsetzero[i])] = 0;
    i++;
  }

  
  arg_table['Q'] = 1;
  arg_table['T'] = 1;
  arg_table['Y'] = 1;

  arg_table['A'] = 2;  
  arg_table['S'] = 2;
  arg_table['M'] = 2;
  arg_table['D'] = 2;
  arg_table['G'] = 2;
  arg_table['E'] = 2;
  arg_table['N'] = 2;
  arg_table['R'] = 2;
 
  arg_table['I'] = 3;

  arg_table['F'] = NUMFUN+1;

  arg_table['V'] = SPACEDIM;

//  arg_table['H'] = 2;

  return 0;
}


int init_tables_copy(void)
{

  int i;

  // node copy functions

  for (i=0;i<OPTABLE;i++)
  {
    copy_table[i] = &copy_non_leaf;
  };

  // leaf nodes
  copy_table[((int) const_op_char)] = &copy_const_node;  
  copy_table[((int) par_op_char)] = &copy_par_node;
//  copy_table['H'] = &copy_par_fun_node;
  
  return 0;
}

int init_io_tables(int arg_table[])
{
  int i;
//  int op;

  for (i=0;i<OPTABLE;i++)
  {
    str2leaf_table[i] = NULL;
  };

  // leaf nodes
#ifdef INTCONSTS
  str2leaf_table[0] = &str2const_node_int;  
#else
  str2leaf_table[0] = &str2const_node;  
#endif
  str2leaf_table[1] = &str2par_node;
  str2leaf_table[2] = &str2par_fun_node;


/* fill table for rendering trees into bracket string notation

   code2str_table table linking op codes to function pointers, functions convert node to str in brack notation 
   core of the node2str convenience function  */

  for (i=0;i<OPTABLE;i++)
  {
    code2str_table[i] = &nFunChr;
  };

  /* these are the true leaf nodes for printing purposes, not V (V is leaf during execution) */
#ifdef INTCONSTS
  code2str_table[((int) const_op_char)] = &constChrInt;  
#else
  code2str_table[((int) const_op_char)] = &constChr;  
#endif

  code2str_table[((int) par_op_char)] = &parChr;  
  code2str_table['H'] = &parfunChr;  

/* initialize c string table */

  for (i=0;i<OPTABLE;i++)
  {
    code2c_str_table[i] = &nFunChr_c;
  };

  /* these are the true leaf nodes for printing purposes, not V (V is leaf during execution) */
#ifdef INTCONSTS
  code2c_str_table[((int) const_op_char)] = &constChrInt;  
#else
  code2c_str_table[((int) const_op_char)] = &constChr;  
#endif

  code2c_str_table[((int) par_op_char)] = &parChr_c;  
  code2c_str_table['H'] = &parfunChr;  


  return 0;
}


void evaluate_tree(Node *tree, State S_return)
{ /* evaluates three and deposits result in S_return.data */

#if defined(STEM_NODES)
  Node *child;
  int i;

  // only evaluate the children of the stem node
  for (i=0;i<SPACEDIM;i++)
  { 
    child =  ((Node *) ( (Node *) tree->children[i]));
    S_return.data[i] = (*op_table[( (int)  (  (Node *) child)->op)])(((Node *) child));
  }
 
#else

  S_return.data[0] = (*op_table[( (int)  (  (Node *) tree)->op)])(((Node *) tree)); 

#endif

  return;
}


double run_node(Node* node)
{
/* typically run on tree->children[0]. E.g. run_node(((Node *) tree->children[0]))
*/

  return (*op_table[( (int)  (  (Node *) node)->op)])(((Node *) node));
}



