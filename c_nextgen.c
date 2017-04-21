/* NextGen Computational Core.

  By Will Sijp.


*/

/* ISSUES: MAXTREESTR may be exceeded during runtime: need to account for this.
One way is to calculate expected string length during tree construction?

Trees inside one thread are defined on 2D grid.
When EXPLICITMIG is selected, trees are also grouped inside internal grid cells of width and height $compgridsize.

Tournament selection picks random trees from moving 2D rectangular window around current tree being created by crossing if EXPLICITMIG is not defined. If it is defined, tournament selection picks inside 2D cell enclosing current tree, and tree migration between internal cells takes place explicitly via exchange similar to migration between threads. 

*/

//#define SPECIAL_MUL // multiplying with left multiplicant vector member functions acting as operators on the right member

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


/* MAXCHILD must be greater than or equal to SPACEDIM 
*/

// local preprocessor directives

#define MAXLEAF 25
#define MAXSIZE 12
#define INITNODES 1
#define CONSTSCUTOFF 8
#define MIGRATE 0.08 // proportion of population migrating 
#define MIGRATE_INTERNAL 0.08 // proportion of population migrating internally 
#define MIGTOUR 5
//#define EXPLICITMIG // exchange migrant trees explicitly between cells internal to CPU process

#define MIN_SCORE -1
#define MIGPARSFACT 3 // amplification of parsimony in assessing score migrants
#define MAXDEPTH 4
#define FPR 0.5
#define PPR 0.6
#define PNEW 0.02
#define PROBSWAP 0.9
#define MAXTREES 400000
#define FNAMESIZE 100

// local vars

typedef struct migr {

  char* infiles;
  char* outfiles;

} Migr;




// local function definitions

int mod (int a, int b);
char find_right_bracket(char left_bracket);

Point make_point(int x, int y);
int point2i(Point p, int max_x);
Point calc_max_p(int max_i);

Node* node_from_str_lit(char tree_str[]);
int test_tree_io(void);

Node *pick_F(Node *t, int max_i);
int export_internal(Population pop, Point my_loc, Point dest_loc, int tour, int compgridsize, int migrants);
int internal_migration(Population pop, int tour, int compgridsize, int migrants);
Point get_boundary(int i_iofile, int i, Point max_loc, int migrants);
IntPars str2int_pars(char * text);
int int_pars2bufferline(NodeScore ns, char treebuffer[], IntPars ip );
int NodeScore2bufferline(NodeScore ns, char treebuffer[] );
int migrants2buffer(Population pop, int i_iofile , int migrants ,int compgridsize, Experiment Exp, Point my_loc_external );
int numbers_dot_only(const char *s);

int readmigrants(Migr mig,int i_iofile, Population pop, int migrants, int compgridsize, int tour );
int pop2buffer(Population pop);
int buffer2trees(Node *trees[]);

int read_pop(char *filepath,Population pop);
int tournament_size(Population pop, int is,int tour);
int tournament_scsize(Population pop, int is,int tour);
int reverse_tour_size(Population pop, int is,int tour);
int tournament(Population pop, int is,int tour);
int pick_random_neighbor(int i_tree, Point max_loc , int compgridsize);
int pick_random_neighbor_point(Point my_loc, Point max_loc , int compgridsize);
int rnd_in_area(int i_tree, Point max_loc , int compgridsize);
Rect point2cell(Point my_loc, Point max_loc, int compgridsize);
int rnd_in_area_point(Point my_loc, Point max_loc , int compgridsize);
int tournament_scsize2D(Population pop, Point my_loc,int tour, int compgridsize);
int tournament_size2D(Population pop, Point my_loc,int tour, int compgridsize);
int reverse_tour_scsize2D(Population pop, Point my_loc,int tour, int compgridsize);
int reverse_tour_size2D(Population pop, Point my_loc,int tour, int compgridsize);
int tournament2D(Population pop, Point my_loc,int tour, int compgridsize);
int reversetour2D(Population pop, int i_tree,int tour, int compgridsize);
int reversetour(Population pop,int is,int tour);
void random_init(Population pop,  Experiment Exp);
Population copy_pop(Population pop, Population pop_next);
Population make_pop(int popsize);
Migr make_migr(void);
Migr init_mig_files(Point p, Point max_p);


// global vars

//int nhillsearches=0;
int node_count = 0;



// functions

int mod (int a, int b)
{
  if (b<0)
    return mod(-a,-b);

  int ret = a%b;
  if (ret < 0)
    ret+=b;
  return ret;

}

char find_right_bracket(char left_bracket)
{

  if (left_bracket == '(')
  {
    return left_bracket + 1;
  }
  else
  {
    return left_bracket + 2;
  }  
}




void split(char text[],int *I, char left_bracket, char delim)
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

  char right_bracket = find_right_bracket(left_bracket);
  
  while (*(text+cursor) !='\0')
  {

    if (*(text+cursor)==left_bracket)
      counter += 1;
    else if (*(text+cursor)==right_bracket)
      counter -= 1;

    if (*(text+cursor)==delim && counter==0)
    {
      
      text[cursor]='\0';
      *(I+icomma)=cursor+1;      
      icomma++;
    } 

    cursor+=1;
  }

  *(I+icomma)=0;  
}


/*Function to find minimum of x and y*/
int min(int x, int y)
{
  return y ^ ((x ^ y) & -(x < y));
}
 
/*Function to find maximum of x and y*/
int max(int x, int y)
{
  return x ^ ((x ^ y) & -(x < y)); 
}

Point make_point(int x, int y)
{
  Point p;
  p.x = x;
  p.y = y;

  return p;
}


Point i2point(int i, int max_x)
{
  // translates index i into x and y values, modulo the max values
  // deposits values into xy 

  Point p;

  p.y = i/max_x;
  p.x = i - max_x*p.y;

  return p;
}

int point2i(Point p, int max_x)
{
  // translates index i into x and y values, modulo the max values
  // deposits values into xy 
 
  return p.y*max_x + p.x;
}

Point calc_max_p(int max_i) 
{
  // Calculate grid boundaries based on total linear size.
  // assumes number of form floor(s(max_i))*floor(max_i/floor(s(max_i))) are given

  Point max_p;
  int max_x = (int) sqrt(max_i);

  max_p.x = max_x;
  max_p.y = max_i/max_x;

  return max_p;
}


int *make_int(int value)
{
/*
  Create new integer constant node with malloc. Returns pointer to new int.
*/
  int *ptr = (int *)malloc(sizeof(int));

  if (ptr == NULL)
  {
    fprintf(stderr,"Error: cannot malloc new int in make_int.\n");
    exit(-1);
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
    exit(-1);
  }

  *ptr=value;
  return ptr;
};



int free_int_pars(IntPars ip)
{
  free(ip.data);
  ip.data = NULL;

  return 0;
}


int free_node_score(NodeScore ns)
{
  free_node(ns.node);
#ifdef EVOLVEIC
  free_state(ns.S_i);
#endif

#if INTPARS > 0
  free_int_pars(ns.ip);
#endif

  return 0;
}

int free_pop(Population pop)
{

  int j;

  for (j=0; j < pop.popsize; j++)
  {
    free_node_score(pop.pop[j]);
  }

  free(pop.pop);

  return 0;
}





IntPars make_int_pars(int size)
{
  IntPars ip;

  ip.data = (int *)malloc(sizeof(int)*size);

  return ip;
}







NodeScore make_node_score(int size)
{

  NodeScore ns;

  ns.score = 1e19;
  ns.node = talloc();
#ifdef EVOLVEIC
  ns.S_i = make_state(size);
#endif

#if INTPARS > 0
  ns.ip = make_int_pars(INTPARS);
#endif 

  return ns;
}


NodeScore copy_node_score(NodeScore ns)
{

  NodeScore new_ns;

  new_ns.score = ns.score;
  new_ns.node = copy_node(ns.node);
#ifdef EVOLVEIC
  new_ns.S_i = make_state(ns.S_i.size);
#endif

#if INTPARS > 0
  new_ns.ip = make_int_pars(INTPARS);
#endif 

  return new_ns;

}


int check_numpar(int numpar)
{
  if (numpar > MAXPAR)
  {
    fprintf(stderr,"Warning! Number of variables numpar=%d, with SPACEDIM=%d, exceeds MAXPAR=%d \n",numpar,SPACEDIM,MAXPAR);
    return 1;
  }

  return 0;
}







int is_null_node(Node *tree)
{ /* check if this is a non-initialized Node */

  int i;
  for (i=0;i<MAXCHILD;i++)
  {
    if ( (Node *)  (tree->children[i]) != NULL)
    {
      return 1;
    }

  }
  return 0;

}


Node *talloc(void)
{
  Node *newtree =  (Node *) calloc(1,sizeof(Node));
  if (newtree == NULL)
  {
    fprintf(stderr,"Error: cannot calloc newtree in talloc.\n");
    exit(-1);
  }

  int i;
  for (i=0;i<MAXCHILD;i++)
  {
    newtree->children[i]=NULL;
  }

  node_count++;
  return newtree;
};


int make_test_tree(char test_tree_str[])
{
  char tmp_str[MAXTREESTR];
  Node *player;
 
  State S_result = make_state(SPACEDIM);

  sprintf(tmp_str, "%s",test_tree_str  );  
  player = str2node(tmp_str,'(',',');

  (*code2str_table[(int) player->op])(player,tmp_str);

  evaluate_tree(player, S_result);

#ifdef INTSTATES
  fprintf(stderr,"%s %d\n",tmp_str, S_result.data[0]);
#else
  fprintf(stderr,"%s %g\n",tmp_str, S_result.data[0]); // double
#endif

  free_node(player);
  return 0;
}


Node* node_from_str_lit(char tree_str[])
{
  char tmp_str[MAXTREESTR];
 
  sprintf(tmp_str, "%s",tree_str  );  
  return str2node(tmp_str,'(',',');

}

int test_tree_io(void)
{

#ifdef STEM_NODES

#ifdef INTSTATES

#else

#endif

//  make_test_tree("((5,4)A,(5,4)A)V"  );  
#else

#ifdef INTSTATES


#endif

//  make_test_tree("((5,4)A,(5,4)A)V"  );  
  make_test_tree("(5,4)A"  );  

  make_test_tree("(4,5)A"  );  
  make_test_tree("(5,4)S"  );  
  make_test_tree("(5,4)M"  );  

#endif

  return 0;
}


char* node2str(Node *tree, char trstr[])
{
/* convert the argument node to bracket string trstr recursively using the table of codes code2str_table 
   convenience function for code2str_table
*/
  return (*code2str_table[tree->op])(tree,trstr);
};

/* Random selection utilities */

Node *child_choice(Node *t)
{
  /* Choose a child node randomly */
  return ((Node *) t->children[rand()%arg_table[t->op]]);
};

char string_choice(char str[])
{
  int i;

  /* Choose a character inside a string randomly */
  int n = strlen(str);
  if (n>1)
  {
    i = rand()%n;
  }
  else
  { 
    i = 0;
  }

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







Node *makerandomtree(int maxdepth,double fpr,double ppr, int max_par)
{
  /* Create a random tree 
     fpr probability of creating normal node (recursively)
     ppr probability of creating parameter node

  */

  char op;
  int i=0;
  int rnd;
 
  Node *result;

  rnd = rand();

  if ( (rnd<fpr*RAND_MAX) && (maxdepth>0) )
  { /* only recursion: create random node and call makerandomtree to make children */


#ifdef INSTRSET_RARE
    if (rnd<0.1*fpr*RAND_MAX)
    {
      op = string_choice(instrset_rare);
    }
    else
    {
      op = string_choice(instrset);
    }
#else
    op = string_choice(instrset);
#endif

    result = init_node(op);

    while ( i<arg_table[(int) op] )  
    {  
      result->children[i] = ((Node *)  makerandomtree(maxdepth-1,fpr,ppr, max_par) );
      i++;
    };
  }
  else
  {  /* Terminal. Create random leaf node. */

    result=makerandomterminal(maxdepth, fpr, ppr, max_par);
  
  };
  return result;  
};



NodeScore makerandomns(int maxdepth,double fpr,double ppr, Experiment Exp)
{
  NodeScore ns;

#ifdef STEM_NODES
  int i;

  Node * newtree = init_node(instrsetstem[0]);
 
  for (i=0;i<arg_table[newtree->op];i++)
  {
    newtree->children[i] = ( (Node *)  makerandomtree(maxdepth,fpr,ppr,numpar) );
  }

#else
  Node * newtree = makerandomtree(maxdepth,fpr,ppr, numpar);
#endif

  ns.score = 1e19;
  ns.node = newtree;

#ifdef EVOLVEIC
  ns.S_i = makerandomstate(Exp.S_i);
#endif

#if INTPARS > 0
  ns.ip = makerandomint_pars(Exp);
#endif

  return ns;

}









int is_terminal(Node *t)
{
/* Determines whether node is a terminal node. Returns 1 if so, 0 if not. */


//  if (t == NULL)
//  {
//    fprintf(stderr,"YO! \n");
//    exit(-1);
//  }

//  fprintf(stderr,"(%s, %c ) ", instrsetleaf, t->op);


# ifdef SUPERLEAF
  if (strchr(instrsetsuperleaf, t->op)  != NULL)
#else
  if (strchr(instrsetleaf, t->op)  != NULL)
//  if ((t->op == 'P') || (t->op == 'C'))  // faster?

#endif

  {
    return 1;
  }
  else
  {
    return 0;
  }
}



IntPars crossover_int_pars(IntPars ip1, IntPars ip2)
{

  int i;
  int divisor;

  IntPars ip = make_int_pars(INTPARS);

  divisor = rand()%(INTPARS+1);

  for (i=0;i<divisor;i++)
  {
    *(ip.data+i) = *(ip1.data+i);
  }

  for (i=divisor;i<INTPARS;i++)
  {
    *(ip.data+i) = *(ip2.data+i);
  }

  return ip;
}



NodeScore crossover_ns(NodeScore ns1,NodeScore ns2, double probswap, int top)
{
  NodeScore newns;

  Node * newtree = unifcross_node(ns1.node,ns2.node, probswap, top);
//  State S_i = make_copy_state(ns1.S_i);  // should become evolutionary crossover

//  Node * newtree = crossover(ns1.node,ns2.node, probswap, top);

  newns.score = 1e19;
  newns.node = newtree;

#ifdef EVOLVEIC
  newns.S_i = crossover_states(ns1.S_i, ns2.S_i);
# endif

#if INTPARS > 0
  newns.ip = crossover_int_pars(ns1.ip, ns2.ip);

#endif

  return newns;

}

Node *unifcross_node(Node *t1,Node *t2, double probswap, int top)
{
  /* Uniform crossover

  */
  int i=0;
  Node *result;

  if ( (top==0) && (rand()>probswap*RAND_MAX) )
  { /* recursion terminal. swap t1 for t2 */
    return ((Node *) copy_node(t2));
  }
  else
  {
    if (  (is_terminal(t1) == 1) || (is_terminal(t2) == 1) )
    { /* recursion terminal. t1 or t2 is a leaf node: randomly return a copy of t1 or t2, this to achieve the mixing */

      if (rand()>probswap*RAND_MAX) 
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

      if (rand()>probswap*RAND_MAX)
      {
        result = init_node(t2->op );
      }
      else
      {
        result = init_node(t1->op );
      }

      while (i<arg_table[t1->op])
      {
        result->children[i] = ((Node *)  unifcross_node(((Node *) t1->children[i]) , ((Node *) t2->children[i]),probswap, 0) );
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

Node *crossover(Node *t1,Node *t2,double probswap,int top)
{
  /* to be updated later, with NodeScore return value instead of Node
  */

  int i=0;
  Node *result;

  if ((rand()<probswap*RAND_MAX) && (top==0))
  {
    return ((Node *) copy_node(t2));
  }
  else
  {
  
    if (  (is_terminal(t1) == 0) && (is_terminal(t2) == 0) ) // t1 and t2 both not leaf nodes
    {
      result = init_node(t1->op );

      while ( (i<arg_table[t1->op])  )
      {
  
        result->children[i] = ((Node *)  crossover(((Node *) t1->children[i]) , child_choice(t2),probswap, 0) ) ;
        i++;
      };
    }
    else
    { /* leaf nodes return copies of node t1 depending on dice and top!=0 */ 
      result=((Node *) copy_node(t1));
    }
    return ((Node *) result);
  };
};



int mutconsts(double *consts[],int size, double mutationrate)
{
  int k2;
  double tmpval=0;
 
  if (size>0)
  {
    tmpval=rand();
    if ( tmpval<0.1*mutationrate*RAND_MAX )
    {
      k2=rand()%size;

      if (rand()%2 == 0)
      {
        *consts[k2] =   randomdouble();
      }

    }
  }
return 0;
}


Node *pick_F(Node *t, int max_i)
{
  // picks the first function it finds

  Node *child;
  Node *result;
  int i;

  if (t->op == 'F')
  {
    return t;
  }


  else if (is_terminal(t) == 0) // not a leaf node  
  {
    for (i=0;i<arg_table[t->op];i++)
    {
      child = ((Node *) t->children[i]);

      if ( (result = pick_F(child, max_i)) != NULL)
      {
        return result;
      } 
    }
    return NULL;

  }
  else
  {
    return NULL;
  }


}

Node *mutate(Node *t, double probchange)
{
  /* 
    Go through tree recursively and roll dice to replace a node and graft new random tree. 
    C/ P leaf nodes are not touched.

    e.g. in (p1,2.0)M p1 and 2.0 can be replaced with a subtree
 
    Creates new tree
   
  */

  int i=0;
  Node *result;

  if ( rand()<probchange*RAND_MAX )
  {
    /* dice successful: grow random tree here */
    return makerandomtree(MAXDEPTH,FPR,PPR, numpar);
  }
  else
  {
    if (is_terminal(t) == 0) // not a leaf node
    { /* only recursion: it is not a leaf node, preserve op but apply mutate to children */

      result = init_node(t->op );

      while (i<arg_table[t->op])  
      { /* recursively throw dice to replace children with subtrees */
        result->children[i] = ((Node *)  mutate(((Node *) t->children[i]) ,probchange) );

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

Node *swapmut(Node *t, double probchange)
{
/* Swap mutation, rolls dice to swap out a suitable operation somewhere in the tree with one from instrsetswap. 
   If node t takes 2 args, swapmut replaces op in node with a randomly chosen 2-arg op 
   at most 1 node is swapped.
*/

  int i=0;
  Node *result;

  if (is_terminal(t) == 0) // not a leaf node
  {
    if ( (rand()<probchange*RAND_MAX) && (arg_table[t->op]==2) )
    {
      /* terminal: dice successful and arg condition met: copy t with changed op */
      result = copy_node(t);
      result->op = string_choice(instrsetswap);
      return result;
    }
    else
    {
        /* only recursion: create new node preserving op and applying swapmut to children */

        result = init_node(t->op );

        while (i<arg_table[t->op])  
        {
          result->children[i] = ((Node *)  swapmut(((Node *) t->children[i]) ,probchange) );
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



Node *hoistmut(Node *t, double probreturn,int top)
{
  /* hoist mutation: throw dice to replace nodes with one of their children.
     acts to shorten/ simplify trees again.
 */
 
  Node *result;

  if (top==1)
  { /* only applied on root: root node always copied */
    return copy_node(  hoistmut(t,probreturn,0) ); 
  }

  if (rand()<probreturn*RAND_MAX )
  { /* dice result induces hoisting on non-terminal node */

    if (is_terminal(t) == 0)  // not leaf node
    { /*  only recursion. Hoist up one of the child nodes */
      result = ((Node *)  hoistmut(((Node *) child_choice(t) ) ,probreturn,0) );
    }     
    else // leaf node
    { /* terminal: leaf nodes returned as is  */
      result = t;
    }
  }
  else
  { /* dice result induces no hoisting */

    result = t;

  }
  return result;
};


Node* allmut_node(Node *tree, double probchange)
{
  int rnd;
  Node * newtree;

  rnd = rand()%4;
  if (rnd < 1)
  {
    newtree = hoistmut(tree, probchange,1);
  }
  else if (rnd < 2)
  {
    newtree = swapmut(tree, probchange);
  }
  else
  {
    newtree = mutate(tree, probchange);
  }

  return newtree;
}


NodeScore allmut(NodeScore ns, double probchange, Experiment Exp)
{
  Node * newtree;
  NodeScore newns;

#ifdef STEM_NODES
  int i;

  newtree = init_node(instrsetstem[0]);

  for (i=0;i<arg_table[newtree->op];i++)
  {
    newtree->children[i] =  allmut_node( ((Node *) ns.node->children[i]) , probchange );
  }

#else
  newtree = allmut_node(ns.node, probchange);


#endif

  newns.score = 1e19;
  newns.node = newtree;

#ifdef EVOLVEIC
  newns.S_i = mut_state(ns.S_i);  // will become mutation on ns.S_i
#endif

#if INTPARS > 0
  newns.ip = mut_int_pars(ns.ip, Exp);
#endif 

  return newns;

};



/*
To create int:

ptr = (int *)malloc(sizeof(int));

create the doubles and ints this way for constant and param leaf nodes

The context of the children must come from the function of the node,
yielding the type

*/




int compare_nodes(Node *t1, Node *t2)
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
    if  (leaves_equal(t1,t2) == 0)
    {
      return 0;
    }
    else
    { /* compare children pair-wise if any. if one of the children returns 1, return 1 */

      for (j=0;j<arg_table[t1->op];j++)  
      {
        if (compare_nodes( ((Node *) t1->children[j]) , ((Node *) t2->children[j]) ) == 1)
        {
          return 1;
        }
      }
    }
    return 0;
  }
}



int compare_ns(NodeScore ns1, NodeScore ns2)
{
/* determine whether 2 NodeScores are identical. returns 0 when identical*/

  int S_i_flag=0;
  int ip_flag = 0;

  int ns_flag = compare_nodes(ns1.node,ns2.node);

#ifdef EVOLVEIC

  S_i_flag = compare_state(ns1.S_i, ns2.S_i);

#endif

#if INTPARS > 0
 
  int i;
  for (i=0;i<INTPARS;i++)
  {
    if ( ns1.ip.data[i] != ns2.ip.data[i]   )  
    {
      ip_flag = 1; 
    }
  }
#endif

  return ns_flag + S_i_flag + ip_flag;

}




int does_tree_exist(NodeScore newns, Population pop, Point my_loc, int compgridsize)
{
  // checks if tree is in array of trees. returns index if yes, -1 otherwise
  // spiral around my_loc, scanning most likely locations first


  int i,j, i_tree;
 
  Rect cell = point2cell(my_loc, pop.max_loc, compgridsize);

  for (i=cell.min.x;i<cell.max.x;i++)
  {  
    for (j=cell.min.y;j<cell.max.y;j++)
    {  
      my_loc.x = i; // reusing my_loc for different purpose
      my_loc.y = j;

      i_tree = point2i(my_loc, pop.max_loc.x);
      if ((compare_ns(newns, pop.pop[i_tree]) == 0) && (i_tree < pop.popsize))
      {
        return i_tree;
      }
    } 
  }
  return -1;
}

int conparcount(Node *t, int *current)
{
  /* counts constants and parameters in a tree */
  int j=0;

  if   ( strchr(instrsetleaf, t->op ) != NULL) 
  {
    (*current)++;
    return 1;
  }
  else
  {
    while (j<arg_table[t->op])  
    {  
      if ( conparcount( ((Node *) t->children[j]),current) > -1 )     
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





int init_buffer(void)
{
/* initialize the buffer */


  treebuffer[0] = '\0';
 
  return 0;
}


int export_internal(Population pop, Point my_loc, Point dest_loc, int tour, int compgridsize, int migrants)
{ /* export $migrants internal migrants from cell around Point my_loc to cell around Point dest_loc
     identify pairs of indices into trees array for tree pointers to be swapped     

     origin index determined via tournament_scsize2D
     destination index via reverse_tour_size2D


*/

  Node *tmptree_ptr; // temporary storage for tree swap

  int i, i_orig, i_dest;

  for (i=0;i<migrants;i++)
  {

    i_orig = tournament_scsize2D(pop,my_loc,MIGTOUR, compgridsize);
    i_dest = reverse_tour_scsize2D(pop,dest_loc, tour, compgridsize);


    tmptree_ptr = pop.pop[i_orig].node;
    pop.pop[i_orig].node = pop.pop[i_dest].node;
    pop.pop[i_dest].node = tmptree_ptr;
  }

  return 0;
}

int internal_migration(Population pop, int tour, int compgridsize, int migrants)
{
/* loop through all cells
	for each cell
		exports $migrants trees to neighboring cells

*/

  int i_cell, j_cell, max_i_cell, max_j_cell;
  Point my_loc, dest_loc;

// determine cells
// determine what is a neighbor

  max_i_cell = pop.max_loc.x/compgridsize;
  max_j_cell = pop.max_loc.y/compgridsize;

  for (j_cell=0;j_cell<max_j_cell;j_cell++)
  {
    for (i_cell=0;i_cell<max_i_cell;i_cell++)
    {
      // identify cell by upper left corner. Note: this creates some repeated calcs in tournament.
      my_loc.x = i_cell*compgridsize; 
      my_loc.y = j_cell*compgridsize;

      if (i_cell>0) // do not export into the wall
      { // export towards left neighbor
        dest_loc.x = my_loc.x - 1;
        dest_loc.y = my_loc.y;  
        export_internal(pop, my_loc, dest_loc, tour, compgridsize,migrants); // export to one neighbor
      }

      if (i_cell<max_i_cell-1) // do not export into the wall
      { // export towards right neighbor
        dest_loc.x = my_loc.x + 1;
        dest_loc.y = my_loc.y;  
        export_internal(pop, my_loc, dest_loc, tour, compgridsize,migrants); // export to one neighbor
      }

      if (j_cell>0) // do not export into the wall
      { // export towards upward neighbor
        dest_loc.x = my_loc.x;
        dest_loc.y = my_loc.y-1;  
        export_internal(pop, my_loc, dest_loc, tour, compgridsize,migrants); // export to one neighbor
      }

      if (j_cell<max_j_cell-1) // do not export into the wall
      { // export towards downward neighbor
        dest_loc.x = my_loc.x;
        dest_loc.y = my_loc.y+1;  
        export_internal(pop, my_loc, dest_loc, tour, compgridsize,migrants); // export to one neighbor
      }
    }
  }

  return 0;
}

Point get_boundary(int i_iofile, int i, Point max_loc, int migrants)
{   
/*
 Return point inside appropriate internal compute grid boundary at height proportional to i. 
 This point can be used elsewhere to select prospective migrants in the external grid, or local trees to be replaced by imported migrants.  Locatation depends on direction of origin or export (l,r,u,d)

 i_iofile: index to iofile corresponding to left-right, up-down directions, values 0,1,2,3 
   i: tree index in migrant file: should run from 0 to migrants
   max_loc: corner with highest coordinate values in rectangle of internal compute grid
   migrants: number of migrants exchanged
*/

  Point boundary_loc;

  if (i_iofile==0){ // xm case: export towards/ import from left, pick from left vertical boundary
    boundary_loc.x = 0;
    boundary_loc.y = (max_loc.y * min(i,migrants)/migrants); // locate on axis
  }
  else if (i_iofile==1){ // xp case: export towards right/ import from, pick from right vertical boundary
    boundary_loc.x = max_loc.x - 1;
    boundary_loc.y = (max_loc.y * min(i,migrants)/migrants); // locate on axis
  }
  else if (i_iofile==2){ // ym case: export upward/ import from above (towards lower y), pick from top horizontal boundary
    boundary_loc.x = (max_loc.x * min(i,migrants)/migrants); // locate on axis
    boundary_loc.y = 0;
  }
  else { //  ym case: export downward/ import from below (towards higher y), pick from bottom horizontal boundary
    boundary_loc.x = (max_loc.x * min(i,migrants)/migrants); // locate on axis
    boundary_loc.y = max_loc.y - 1;
  }

  return boundary_loc;

}


IntPars str2int_pars(char * text)
{

  int i;
  char * p_end;
  char * p_next;
  char * tmp_str = (char *)malloc(sizeof(char)*(MAXTREESTR));
  IntPars ip = make_int_pars(INTPARS);
  int end_flag = 0;

  const char start_ch = '(';
  const char end_ch = ')';

  const char s[2] = ",";
  char *token;
  char *saveptr;


  strcpy(tmp_str, text);

//  fprintf(stderr, " %s\n", tmp_str );

  if (*tmp_str != start_ch)
  {
    fprintf(stderr,"Problem reading int pars. No opening bracket. \n"); 
    free_int_pars(ip);
    ip.data = NULL;
    return ip;
  }

  p_next = tmp_str + 1;

  /* get the first token */
  token = strtok_r(p_next, s, &saveptr);
   
  /* walk through other tokens */
  i=0;
  while( (token != NULL) && (end_flag == 0) && (i<INTPARS) ) 
  {

    if ( (p_end=strchr(token, end_ch)) != NULL  )
    {
      end_flag = 1;
      *p_end = '\0';
    }

    *(ip.data+i) = atoi(token);
    
    token = strtok_r(NULL, s, &saveptr);

    i++;
  }

  if ( (i < INTPARS) || (end_flag == 0) )
  {
    fprintf(stderr,"Problem reading int pars. Saved data length %d.\n",i);
    free_int_pars(ip);
    ip.data = NULL;
    return ip;    
  }

  free(tmp_str);

  return ip;
}






int int_pars2bufferline(NodeScore ns, char treebuffer[], IntPars ip )
{
  char tmp_str[MAXTREESTR];
  int i;

  sprintf(tmp_str,"(");
  strcat(treebuffer,  tmp_str ); 

  for (i=0;i<INTPARS;i++)
  {
    sprintf(tmp_str,"%d", *(ip.data+i) );
    strcat(treebuffer,  tmp_str );

    if (i<INTPARS-1)
    {
      sprintf(tmp_str,",");
      strcat(treebuffer,  tmp_str ); 
    }
  }

  sprintf(tmp_str,")");
  strcat(treebuffer,  tmp_str ); 

  return 0;
}


int NodeScore2bufferline(NodeScore ns, char treebuffer[] )
{
  /* add string containing score and tree in bracket notation to treebuffer string */

  char tmp_str[MAXTREESTR];
  const char hor_delim[2] = " ";

  /* add score output_score */
  sprintf(tmp_str,"%.10g", ns.score); // print score value into string
  strcat(treebuffer,  tmp_str ); // add string to treebuffer string
  strcat(treebuffer,  hor_delim ); // add space


#if INTPARS > 0
  int_pars2bufferline(ns, treebuffer, ns.ip);
  strcat(treebuffer,  hor_delim ); 
#endif

#ifdef EVOLVEIC
  state2bufferline(treebuffer, ns.S_i);
  strcat(treebuffer,  hor_delim ); 
#endif

  /* add tree */
  (*code2str_table[(int) ns.node->op])(ns.node,tmp_str); /* convert tree to string */
  strcat(treebuffer,  tmp_str );  // add string to treebuffer string
  strcat(treebuffer,  "\n"  );  // add \n

  return 0;
}


int migrants2buffer(Population pop, int i_iofile , int migrants ,int compgridsize, Experiment Exp, Point my_loc_external )
{
/* creates output file with spatially ordered list of trees, depending on the direction of the output file.

  j: output direction (1-4: xm,xp,ym,yp)

  File IO not included. Works on string buffer treebuffer.

*/

  Point boundary_loc; // used to locate trees within grid

  int i,k;

#if defined (CHECK_MIG_STABILITY)
  double error;
#endif

  init_buffer(); // buffer string to fill with tree representations


#if COMPARTMENTS > 0 

  if ( ( (mod(my_loc_external.x,COMPARTMENTS) == 0) && (i_iofile == 1) ) ||
       ( (mod(my_loc_external.x,COMPARTMENTS) == 1) && (i_iofile == 0) ) ||
       ( (mod(my_loc_external.y,COMPARTMENTS) == 0) && (i_iofile == 3) ) ||
       ( (mod(my_loc_external.y,COMPARTMENTS) == 1) && (i_iofile == 2) ) )

  {  // modified exchange towards x+1

    for (i=0;i<migrants/16;i++)
    { // loop over all migrants to output

      boundary_loc = get_boundary(i_iofile, i, pop.max_loc, migrants); // internal cell boundary

      k = tournament_size2D(pop,boundary_loc,MIGTOUR, compgridsize);
 
 
#ifdef CHECK_MIG_STABILITY
    /* check time step stability to export only stable trees */
      error = get_score(pop.pop[k].node, Exp.S_i, Exp);
      if (fabs(pop.pop[k].score-error)>5.0)
      {
        pop.pop[k].score = 1e19;
      }
      else
      {
#endif
      /* add strings containing score and tree in bracket notation to treebuffer string */

        NodeScore2bufferline(pop.pop[k], treebuffer );

#ifdef CHECK_MIG_STABILITY
      }
#endif

    } // end migrants loop  

  }
  else
  {

    for (i=0;i<migrants;i++)
    { // loop over all migrants to output

      boundary_loc = get_boundary(i_iofile, i, pop.max_loc, migrants); // internal cell boundary

      k = tournament_scsize2D(pop,boundary_loc,MIGTOUR, compgridsize);

#ifdef CHECK_MIG_STABILITY
    /* check time step stability to export only stable trees */
      error = get_score(pop.pop[k].node, Exp.S_i, Exp);
      if (fabs(pop.pop[k].score-error)>5.0)
      {
        pop.pop[k].score = 1e19;
      }
      else
      {
#endif
      /* add strings containing score and tree in bracket notation to treebuffer string */

        NodeScore2bufferline(pop.pop[k], treebuffer );

#ifdef CHECK_MIG_STABILITY
      }
#endif

    } // end migrants loop  

  }


#else


  for (i=0;i<migrants;i++)
  { // loop over all migrants to output

    boundary_loc = get_boundary(i_iofile, i, pop.max_loc, migrants); // internal cell boundary

    k = tournament_scsize2D(pop,boundary_loc,MIGTOUR, compgridsize);
 
 
#ifdef CHECK_MIG_STABILITY
    /* check time step stability to export only stable trees */
    error = get_score(pop.pop[k].node, Exp.S_i, Exp);
    if (fabs(pop.pop[k].score-error)>5.0)
    {
      pop.pop[k].score = 1e19;
    }
    else
    {
#endif
      /* add strings containing score and tree in bracket notation to treebuffer string */

      NodeScore2bufferline(pop.pop[k], treebuffer );

#ifdef CHECK_MIG_STABILITY
    }
#endif

  } // end migrants loop  
#endif




  return 0;
}

int numbers_dot_only(const char *s)
{
  char *valid_chars = ".e+-";

    while (*s) {

        if ((isdigit(*s) == 0) && ( strchr(valid_chars, *s) == NULL  ))
          return 0; // signals it's not a number

        s++;
    }

    return 1; // signals it's a number
}

double my_strtod(const char *token)
{

  double tmp_score;

  if (  numbers_dot_only(token) == 0 ) // separate check for garbled number (containing chars not number or .)
  {
    fprintf(stderr,"Error reading file: non-numerical string for score detected: %s. \n",token);
    return -1.0; /* abort and trigger failure. */          
  }


  tmp_score = strtod(token,NULL);

  if (tmp_score == 0.0) // indicates atof failed
  {
    fprintf(stderr,"Error reading file: could not convert float \n");
    return 0.0; /* abort and trigger failure. */          
  }

  return tmp_score;

}


NodeScore bufferline2NodeScore(char tmp_str[], char * infile)
{ /* read NodeScore from a line in the buffer 
     argument tmp_str is buffer line 
     infile argument used for error reporting only */


  double tmp_score;

  const char hor_delim[2] = " ";
  char *token;  
  char *token_wrong;  
  char *pos;

  NodeScore ns;

#if INTPARS > 0
  IntPars ip;
#endif

#ifdef EVOLVEIC
  State S;
#endif

  ns.score = 1e19;
  ns.node = NULL;
  
  token = strtok(tmp_str,hor_delim);   

  tmp_score = my_strtod(token);

  if (tmp_score <= 0.0) // indicates string to double conversion failed
  {
    //fprintf(stderr,"Error reading file: could not convert float \n");
    return ns; /* abort and trigger failure via ns values */          
  }

  if (tmp_score < MIN_SCORE) // separate check for implausible values
  {
    fprintf(stderr,"Error reading file: implausibly low score detected. \n");
    return ns; /* abort and trigger failure via ns values */          
  }


#if INTPARS > 0
  token = strtok(NULL,hor_delim);
  if (token != NULL) 
  {
    ip = str2int_pars(token);
       
    if (ip.data == NULL)
    {
      fprintf(stderr,"Error reading state. \n");
      return ns; /* abort and trigger failure via ns values */         
    }
  }
  else
  { 
    fprintf(stderr,"Error in reading file %s (wrong format?). \n",infile);   
    /* abort and trigger failure via ns values */       
  }  
#endif

#ifdef EVOLVEIC
  token = strtok(NULL,hor_delim);
  if (token != NULL) 
  {
    S = str2state(token);
       
    if (S.data == NULL)
    {
      fprintf(stderr,"Error reading state. \n");
      return ns; /* abort and trigger failure via ns values */         
    }
  }
  else
  { 
    fprintf(stderr,"Error in reading file %s (wrong format?). \n",infile);   
    /* abort and trigger failure via ns values */       
  }  
#endif

  
  token = strtok(NULL,hor_delim);
  if (token != NULL) 
  {

//    fprintf(stderr,"%s \n",token);

    token_wrong = strtok(NULL,hor_delim);  // is there anything else on the line?
    if (token_wrong != NULL)  // if so, error
    {
      fprintf(stderr,"Error reading file: to many items on line \n");
      return ns; /* abort and trigger failure via ns values */      
    }


    // file reading checks
    if (( token[0] == '\n' ) || ( token[0] == '\0' )) // empty line? this can't be an empty line
    {
      fprintf(stderr,"Error reading file: empty line \n");
      return ns; /* abort and trigger failure via ns values */
    }

    if ((check_brackets(token,'(') == -1) || (check_brackets(token,'[') == -1) )  // brackets wrong? bracket tree integrity check
    {
      fprintf(stderr,"Error reading file: bracket tree integrity check failed for %s \n",token);
      return ns; /* abort and trigger failure via ns values */
    }
    // end integrity checks

    if ( (pos=strchr(token,'\n')) != NULL)
      *pos = '\0';

    ns.score = tmp_score;
    ns.node = str2node(token,'(',',');

#ifdef EVOLVEIC
    ns.S_i = S;
#endif

#if INTPARS > 0
    ns.ip = ip;
#endif

    if (check_node(ns.node) == 1) /* tree integrity check after translation to node */
    {
      fprintf(stderr,"Triggering abort from bufferline2NodeScore \n");
      ns.score=1e19;
      ns.node=NULL; /* abort and trigger failure */
      return ns;   
    }
  }
  else
  { /* if there's no tree after the line with a score, trigger failure */
    fprintf(stderr,"Error in reading file %s (wrong format?). \n",infile);   
    /* abort and trigger failure via ns values */       
  }

  return ns;
}




int readmigrants(Migr mig,int i_iofile, Population pop, int migrants, int compgridsize, int tour )
{
/* fill an array of tree structure pointers *trees[] by translating text tree representations (from file) to structures  

   arg "migrants" is the number of migrants expected, used for location in tree placement on boundary. This function will read as many migrant trees as are available in the file. 

   Produces warning to stderr when file reading error occurs.

   This function sometimes fails due to file corruption. At the moment, this can lead to corrupted trees entering the population, possibly with P/ C nodes with parents other than V nodes. floating point exceptions have been observed in conjunction with file reading errors. Possible fixes include a tree vetting function.

   Includes file IO

   opposite of migrants2buffer
*/

  int i, i_tour;
  Point boundary_loc;
  NodeScore ns;

  char tmp_str[MAXTREESTR];
  
  i=0; // counts the migrant trees

  FILE *fp;

  fp = fopen(mig.infiles + i_iofile*FNAMESIZE,"r");  // file pointer fp
  if (fp != NULL)
  {
    while((fgets (tmp_str,MAXTREESTR,fp)!=NULL ) && (i<migrants)  )  // read line by line from file
    {
      ns = bufferline2NodeScore(tmp_str,mig.infiles + i_iofile*FNAMESIZE); // convert to NodeScore, file arg only for error reporting 

      if (ns.node == NULL)
      {
        i = -1;
        break;
      }
      i++;
      boundary_loc = get_boundary(i_iofile, i , pop.max_loc, migrants);
      i_tour = reverse_tour_scsize2D(pop,boundary_loc, tour, compgridsize);

      free_node_score(pop.pop[i_tour]); /* This tree will be replaced. Note: free_node checks for NULL, no need to do it here */
      pop.pop[i_tour] = ns;
 
    }

    fclose(fp);

#if COMPARTMENTS == 0
    if (i<migrants)
    {
      fprintf(stderr,"Not enough lines in file %s. \n", mig.infiles + i_iofile*FNAMESIZE);
      return -1;
    }
#endif
 
    return i;
  }
//  fprintf(stderr,"Warning: cannot open %s\n",mig.infiles + i_iofile*FNAMESIZE);
  return(-2); /* can't open file, not a big problem */
}







int pop2buffer(Population pop)
{
/* fill the buffer with scores and text representations of tree structures from an array 

*/
  int i;

  init_buffer();
  for (i=0;i<pop.popsize;i++)
  {
  //  fprintf(stderr,"%d ",i);
    NodeScore2bufferline(*(pop.pop+i), treebuffer ); // write line to buffer. Same call as from migrants2buffer
  }  
  return 0;
}




int check_brackets(char *tmp_str, char left_bracket)
{
  int i=0;
  int b=0;
  int t=0;

  char right_bracket = find_right_bracket(left_bracket);

//  int l = strlen(tmp_str);


  while (tmp_str[i] != '\0')
  {
 
    if (tmp_str[i] == left_bracket)
    {
      b++;
      t++;
    }
    else if (tmp_str[i] == right_bracket)
    {
      b--;
      t++;
    }
    else if ( tmp_str[i] == ' '  )
    {
      fprintf(stderr,"Error: tree str cannot contain spaces");
      return -1;    
    }

    i++;
  }
  if (b!=0) 
  {
    return -1;
  }

  return t;   
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

int read_data(const char *filepath,Node *trees[],int maxtrees)
{
  char tmp_str[MAXTREESTR];
  int i=0;

  FILE *fp;

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
     
      if ((check_brackets(tmp_str,'(') == -1) || (check_brackets(tmp_str,'[') == -1) )
      {
        fprintf(stderr,"Error reading pop: bracket tree integrity check failed for %s \n",tmp_str);
        i=-1; /* abort and trigger failure */
        break;
      }

      strtok(tmp_str,"\n");
      trees[i] = str2node(tmp_str,'(',',');
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




int read_pop(char *filepath,Population pop)
{

  int i;
  NodeScore ns;

  char tmp_str[MAXTREESTR];

  i=0; // counts the trees

  FILE *fp;

  fp = fopen(filepath,"r");  // file pointer fp
  if (fp != NULL)
  {
    while((fgets (tmp_str,MAXTREESTR,fp)!=NULL ) && (i<pop.popsize)  )  // read line by line from file
    {
      ns = bufferline2NodeScore(tmp_str,filepath); // convert to NodeScore, file arg only for error reporting 

      if (check_node(ns.node) == 1) /* tree integrity check after translation to node */
      {
        i = -1;
        break;
      }

      pop.pop[i] = ns;
      i++;
    }

    fclose(fp);
 
    if (i<pop.popsize)
    {
      fprintf(stderr,"Not enough lines in %s. \n",filepath);
      return -1;
    }

    return i;
  }
  
  return(-2); /* can't open file. */
}







int tournament_size(Population pop, int is,int tour)
{
/* tournament selection on size alone. returns index of chosen tree
*/
  int i,k,leafcount;
  int kmin=0;
  double M=1e20; /* minimum leafcount */

  int *current = make_int(0);

  for (i=0;i<tour;i++)
  {
    k = is + rand()%pop.popsize; /* choose population wide within segment is to is+popsize */

    leafcount = conparcount(pop.pop[k].node, current);

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



int tournament_scsize(Population pop, int is,int tour)
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
      k = is + rand()%pop.popsize; /* choose population wide in segment is to is+popsize*/
      leafcount = conparcount(pop.pop[k].node, current);
      *current = 0;  /* always reset leafcount current afterwards */
      l++;
    }while( ( (leafcount < 2) || (leafcount > MAXLEAF) ) && (l<10) );

    if ( ( sizescore = ( pop.pop[k].score + 2*PARSIMONY*( (double) leafcount) )  ) < M)  
    { /* find the minimum value over the tournament */
      M = sizescore;
      kmin = k;
    }
  }
  free(current);

  return kmin;
}




int reverse_tour_size(Population pop, int is,int tour)
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
    k = is + rand()%pop.popsize; /* choose population wide inside segment is to is+ popsize */

    leafcount = conparcount(pop.pop[k].node, current); /* calculate number of consts and pars (leaves) */
    *current = 0;      /* always reset leafcount current afterwards */

    if ( ( sizescore = ( pop.pop[k].score + 2*PARSIMONY*( (double) leafcount) )  ) > M) 
    { /* calculate maximum sizescore */
      M = sizescore;
      kmax = k;
    }
  }
  free(current);

  return kmax;
}





int tournament(Population pop, int is,int tour)
{
/* tournament selection
*/
  int i,k;
  int kmin=0;
  double M=1e20;

  for (i=0;i<tour;i++)
  {
    k = is + rand()%pop.popsize; /* choose population wide */

    if (pop.pop[k].score < M)
    {
      M = pop.pop[k].score;
      kmin = k;
    }
  }
  return kmin;
}


int pick_random_neighbor(int i_tree, Point max_loc , int compgridsize)
{
/* Picks random neighbor. Used for dispersed tree computation grid.

Same as pick_random_neighbor_point, but with index i_tree instead of Point my_loc as arg. */

  return pick_random_neighbor_point(i2point(i_tree, max_loc.x), max_loc , compgridsize);
}

int pick_random_neighbor_point(Point my_loc, Point max_loc , int compgridsize)
{
/* Picks random neighbor. Used for dispersed tree computation grid.

Same as pick_random_neighbor, but with Point my_loc as arg instead of index */

  Point rnd_loc; // point to return

  int x_m,y_m; // minimum x, y

  if (my_loc.x < 2*compgridsize)
  {
    x_m = 0;   
  }
  else if (my_loc.x > max_loc.x - 2*compgridsize)
  {
    x_m = max_loc.x - 2*compgridsize;
  }
  else  // an inner point
  {
    x_m = my_loc.x - compgridsize;
  }

  if (my_loc.y < 2*compgridsize)
  {
    y_m = 0;   
  }
  else if (my_loc.y > max_loc.y - 2*compgridsize)
  {
    y_m = max_loc.y - 2*compgridsize;
  }
  else  // an inner point
  {
    y_m = my_loc.y - compgridsize;
  }

  rnd_loc.x = x_m + rand()%(2*compgridsize);
  rnd_loc.y = y_m + rand()%(2*compgridsize);

  return point2i(rnd_loc, max_loc.x);
}


int rnd_in_area(int i_tree, Point max_loc , int compgridsize)
{
/* Picks random neighbor. Used for dispersed tree computation grid.

Same as pick_random_neighbor_point, but with index i_tree instead of Point my_loc as arg. */

  return rnd_in_area_point(i2point(i_tree, max_loc.x), max_loc , compgridsize);
}


Rect point2cell(Point my_loc, Point max_loc, int compgridsize)
{ // find enclosing cell around point

  Rect cell;

  // upper left corner
  cell.min.x = compgridsize*(my_loc.x/compgridsize);
  cell.min.y = compgridsize*(my_loc.y/compgridsize);

  // lower right corner of cell. Cutoff depending on max_loc in ambient grid
  cell.max.x = min(cell.min.x + compgridsize,max_loc.x);
  cell.max.y = min(cell.min.y + compgridsize,max_loc.y);

  return cell;
}

int rnd_in_area_point(Point my_loc, Point max_loc , int compgridsize)
{
/* Determine which compgrid cell the point my_loc is in, and pick rnd point in area. Used for dispersed tree computation grid.

   compgridsize: width of a compgrid cell

  */

  Point rnd_loc; // point to return
  Rect cell = point2cell(my_loc, max_loc, compgridsize);

  rnd_loc.x = cell.min.x + rand()%(cell.max.x - cell.min.x );
  rnd_loc.y = cell.min.y + rand()%(cell.max.y - cell.min.y );

  return point2i(rnd_loc, max_loc.x);
}


int tournament_scsize2D(Population pop, Point my_loc,int tour, int compgridsize)
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

#ifdef EXPLICITMIG
      k = rnd_in_area_point(my_loc, pop.max_loc , compgridsize);
#else
      k = pick_random_neighbor_point(my_loc, pop.max_loc , compgridsize);
#endif

      leafcount = conparcount(pop.pop[k].node, current);
      *current = 0;  /* always reset leafcount current afterwards */
      l++;
    }while( ( (leafcount < 2) || (leafcount > MAXLEAF) ) && (l<10) );

    if ( ( sizescore = ( pop.pop[k].score + MIGPARSFACT*PARSIMONY*( (double) leafcount) )  ) < M)  
    { /* find the minimum value over the tournament */
      M = sizescore;
      kmin = k;
    }
  }
  free(current);

  return kmin;
}


int tournament_size2D(Population pop, Point my_loc,int tour, int compgridsize)
{
/* tournament selection based on score and tree size
*/
  int i,k,l,leafcount;
  int kmin=0;
  int M=100000;  /* minimum value of leafcount */

  int *current = make_int(0);

  for (i=0;i<tour;i++)
  {
    l=0;
    do
    { /* look for a tree with leaves between 2 and MAXLEAF and try no more than 10 times */

#ifdef EXPLICITMIG
      k = rnd_in_area_point(my_loc, pop.max_loc , compgridsize);
#else
      k = pick_random_neighbor_point(my_loc, pop.max_loc , compgridsize);
#endif

      leafcount = conparcount(pop.pop[k].node, current);

      *current = 0;  /* always reset leafcount current afterwards */
      l++;
    }while( ( (leafcount < 2) || (leafcount > MAXLEAF) ) && (l<10) );

    if ( leafcount < M)  
    { /* find the minimum value over the tournament */
      M = leafcount;
      kmin = k;
    }
  }

  free(current);

  return kmin;
}

int reverse_tour_scsize2D(Population pop, Point my_loc,int tour, int compgridsize)
{
/* tournament selection based on score and tree size
*/
  int i,k,l,leafcount;
  int kmax=0;
  double M=0;  /* minimum value of leafcount */
  double sizescore; /* amalgemation of size and score */

  int *current = make_int(0);

  for (i=0;i<tour;i++)
  {
    l=0;
    do
    { /* look for a tree with leaves between 2 and MAXLEAF and try no more than 10 times */

#ifdef EXPLICITMIG
      k = rnd_in_area_point(my_loc, pop.max_loc , compgridsize);
#else
      k = pick_random_neighbor_point(my_loc, pop.max_loc , compgridsize);
#endif

      leafcount = conparcount(pop.pop[k].node, current);
      *current = 0;  /* always reset leafcount current afterwards */
      l++;
    }while( ( (leafcount < 2) || (leafcount > MAXLEAF) ) && (l<10) );

    if ( ( sizescore = ( pop.pop[k].score + MIGPARSFACT*PARSIMONY*( (double) leafcount) )  ) > M)  
    { /* find the minimum value over the tournament */
      M = sizescore;
      kmax = k;
    }
  }
  free(current);

  return kmax;
}



int reverse_tour_size2D(Population pop, Point my_loc,int tour, int compgridsize)
{
/* tournament selection on size and score, seeking maximum values. Used for tree elimination in migration.
   selection takes place around Point my_loc

*/
  int i,k,leafcount;
  int kmax=0;
  double M=0; /* maximum sizescore value */
  double sizescore; /* combination of size and score */

  int *current = make_int(0);

  for (i=0;i<tour;i++)
  {
#ifdef EXPLICITMIG
    k = rnd_in_area_point(my_loc, pop.max_loc , compgridsize);
#else
    k = pick_random_neighbor_point(my_loc, pop.max_loc , compgridsize);
#endif

    leafcount = conparcount(pop.pop[k].node, current); /* calculate number of consts and pars (leaves) */
    *current = 0;      /* always reset leafcount current afterwards */

/*    if ( ( sizescore = pop.pop[k].score*( (double) leafcount) ) > M)   */
    if ( ( sizescore = ( pop.pop[k].score + 2*PARSIMONY*( (double) leafcount) )  ) > M) 
    { /* calculate maximum sizescore */
      M = sizescore;
      kmax = k;
    }
  }
  free(current);

  return kmax;
}

int tournament2D(Population pop, Point my_loc,int tour, int compgridsize)
{
/* tournament selection inside a neighborhood on a 2D grid
 
  for points i close to the left edge: left-most box is 0 to 2*compgridsize
  for inner points i: i-compgridsize to i+compgridsize


  returns linear index i (to be used in score_vals and trees1)

*/

  int i,k;
  int kmin=0;
  double M=1e20;

  for (i=0;i<tour;i++)
  {

#ifdef EXPLICITMIG
    k = rnd_in_area_point(my_loc, pop.max_loc , compgridsize);
#else
    k = pick_random_neighbor_point(my_loc, pop.max_loc , compgridsize);
#endif

//    k = is + rand()%popsize; /* choose population wide */

    if (pop.pop[k].score < M)
    {
      M = pop.pop[k].score;
      kmin = k;
    }
  }
  return kmin;
}

int reversetour2D(Population pop, int i_tree,int tour, int compgridsize)
{
/* tournament selection inside a neighborhood on a 2D grid
  
  for points i close to the left edge: left-most box is 0 to 2*compgridsize
  for inner points i: i-compgridsize to i+compgridsize


  returns linear index i (to be used in score_vals and trees1)

*/

  int i,k;
  int kmax=0;
  double M=0;


  for (i=0;i<tour;i++)
  {

#ifdef EXPLICITMIG
    k = rnd_in_area(i_tree, pop.max_loc , compgridsize);
#else
    k = pick_random_neighbor(i_tree, pop.max_loc , compgridsize);
#endif

//    k = is + rand()%popsize; /* choose population wide */

    if (pop.pop[k].score > M)
    {
      M = pop.pop[k].score;
      kmax = k;
    }
  }
  return kmax;
}



int reversetour(Population pop,int is,int tour)
{
/* inverse tournament selection
*/
  int i,k;
  int kmax=0;
  double M=0;

  for (i=0;i<tour;i++)
  {
    k = is + rand()%pop.popsize;

    if (pop.pop[k].score > M)
    {
      M = pop.pop[k].score;
      kmax = k;
    }
  }
  return kmax;
}





void random_init(Population pop,  Experiment Exp)
{

  int i=0;
  double error;
//  char tmp_str[MAXTREESTR];

  NodeScore newns;

  while (i<pop.popsize)
  {        
    node_count = 0;  /* node_count is a global var */
    newns = makerandomns(MAXDEPTH,FPR,PPR, Exp);

    if (node_count > INITNODES)
    {
      /* score calculation for random init */
 //     fprintf(stderr,"(%g , %g) ", newns.S_i[0], newns.S_i[1]);


      user_functions(newns);

#ifdef EVOLVEIC
      error = get_score(newns.node, newns.S_i, Exp);   
#else 
      error = get_score(newns.node, Exp.S_i, Exp); 
#endif

      if (error <1e18) 
      {
        newns.score = error;

        pop.pop[i] = newns;
       
        i++;
      }
      else
      {
        free_node_score(newns);
      }
    }
    else
    {
      free_node_score(newns);
    }
  }
}



Population copy_pop(Population pop, Population pop_next)
{ /* copy pop_next into pop and free previous NodeScore elements of pop.pop */

  int j;

  if (pop.popsize != pop_next.popsize)
  {
    fprintf(stderr,"Error copying pop: pop sizes differ. \n");
    exit(-1);
  } 

  for (j=0; j < pop.popsize; j++)
  { 
 //   pop.pop[j].score = pop_next.pop[j].score;
    free_node_score(pop.pop[j]);  // free up memory previously allocated to NodeScore members of pop.pop array
    pop.pop[j] = pop_next.pop[j];

  }
  return pop;
}

Population reset_counters_pop(Population population)
{
  population.trees_reused = 0;
  population.totsize = 0; 
  population.toterror = 0.0; 
  population.minerror = 1e20; 
  population.i_min = 0; 

  return population;
}


Population make_pop(int popsize)
{ /* create population array and initialize with NULL NodeScore entries 
     usage:

     pop = make_pop(popsize);
 */

  int i;
  //NodeScore *pop; // population array

  Population population;

  population.popsize = popsize;
  population.max_loc = calc_max_p(popsize);

  population.pop = (NodeScore *)malloc(sizeof(NodeScore)*(popsize+5));  

  for (i=0;i<popsize+5;i++)
  { 
    population.pop[i].score = 1e19;
    population.pop[i].node = NULL;
  }  

  population.i = 0;
  population.my_loc.x = 0;
  population.my_loc.y = 0;

  population = reset_counters_pop(population);

  return population;
}


Migr make_migr(void)
{
  Migr mig;
  mig.infiles = (char *)malloc(sizeof(char)*(4*FNAMESIZE));
  mig.outfiles = (char *)malloc(sizeof(char)*(4*FNAMESIZE));

  return mig;
}

Migr init_mig_files(Point p, Point max_p)
{
  Migr mig = make_migr();

/* create file name strings for migrant files: in from-to format */
 
  // Migrant input files. x increases rightwards, y increases downwards
  sprintf(mig.infiles,"migr%d-%dxp",mod(p.x-1,max_p.x),mod(p.y,max_p.y) ); // from left process
  sprintf(mig.infiles+FNAMESIZE,"migr%d-%dxm",mod(p.x+1,max_p.x),mod(p.y, max_p.y ) ); // from right process
  sprintf(mig.infiles+2*FNAMESIZE,"migr%d-%dyp",mod(p.x,max_p.x),mod(p.y-1,max_p.y) ); // from above process
  sprintf(mig.infiles+3*FNAMESIZE,"migr%d-%dym",mod(p.x,max_p.x),mod(p.y+1,max_p.y) ); // from below process

  // Migrant output files
  sprintf(mig.outfiles,"migr%d-%dxm",p.x,p.y ); // leftward
  sprintf((mig.outfiles+FNAMESIZE),"migr%d-%dxp",p.x,p.y ); // rightward

  sprintf((mig.outfiles+2*FNAMESIZE),"migr%d-%dym",p.x,p.y ); // upward
  sprintf((mig.outfiles+3*FNAMESIZE),"migr%d-%dyp",p.x,p.y ); // downward

  return mig;

}

void c_nextgen(int my_number,int qsubs,int runlen,int popsize, int compgridsize, double mutationrate, int tour)
{

  Population pop; // struct containing population array
  Population pop_next; // next population

  char restartfile[20];

  int *current = make_int(0);

  int *nconsts;
  int size=0;

  int i,j,run, leafcount, new_flag;
  int error_flag=0;

  double error, avgsize, avgerror, avgtreesreused;
  
  double varerror;
  double *consts[MAXCONSTS];

  char tmp_str[MAXTREESTR]; /* used only once */
  char line_str[MAXTREESTR];

  char repfile[FNAMESIZE];
  char elitefile[FNAMESIZE];

  int migrants = MIGRATE*popsize/4;  // number of migrants to output in each of the 4 migrant files.

  NodeScore newns;
  NodeScore tmpns1;
  NodeScore tmpns2;
  NodeScore tmpns3;

  int i_tree1, i_tree2, i_tree_old;

  Point my_loc; // to be used in internal grid
  Point my_loc_external;
  Point max_loc_external; // based on total qsubs

  double existing_score;

  // report file and elite tree output files
  sprintf(repfile,"report%d",my_number);
  sprintf(elitefile,"elite%d",my_number);
  sprintf(restartfile,"pop%d",my_number);

  Experiment Exp = make_experiment(-1); // also initializes Exp.S_i

  pop = make_pop(popsize);
  pop_next = make_pop(popsize);

  if (pop.max_loc.x < 2*compgridsize)
  {
    fprintf(stderr,"Error: compgridsize %d must be at least twice max_loc.x %d. \n", compgridsize, pop.max_loc.x);
    exit(-1);

  }

  srand( getpid()+time(NULL)  );

/*  printf("random seed: %d \n",getpid()+(time.tv_sec*100) + (time.tv_usec/100));
*/

  model_init(Exp);

#ifdef DEBUG
  test_tree_io();
#endif

  //printf("runlen: %d, my_number: %d, qsubs: %d, S_init: %g, numpar: %d \n",runlen,my_number, qsubs, *Exp.S_i.data, numpar);

  max_loc_external = calc_max_p(qsubs);  // max x is qrows and max y is qcols
  my_loc_external = i2point(my_number, max_loc_external.x); // my_number_x is p.x, my_number_y is p.y

  Migr mig = init_mig_files(my_loc_external, max_loc_external);

  nconsts = make_int(0);

  if (popsize > MAXTREES)
  {
    fprintf(stderr,"Error: popsize %d exceeds MAXTREES %d.\n",popsize,MAXTREES);
  }

  int ntrees = read_pop(restartfile, pop);

  if (ntrees <= -1)
  {
    if (ntrees == -1)
    { 
      fprintf(stderr,"Bad pop file detected, populating from random init.\n");
    }

    printf("No data file. Random init \n");
/*    mutationrate = 4.0*mutationrate;   */
    
    random_init(pop,  Exp);
    ntrees = popsize;
  }

  for (run=0;run<runlen;run++)
  {

    for (j=0;j<4;j++) /* Import trees from external processes. Read migrant files. */
    { 
      error_flag = readmigrants(mig ,j ,pop, migrants, compgridsize, tour); /* j arg is sector */
      if (error_flag == -1)
      {
        break;
      }
    }
    if (error_flag == -1) 
    { 
      /* One bad file will cause pop to be dropped. Using function readmigrants_buffered instead (above) allows keeping the pop. */
      fprintf(stderr, "Error in migrant io detected, replacing entire population.\n");
      random_init(pop,  Exp);
    }

#ifdef DEBUG
    for (i=0;i<popsize;i++)
    {
      if (check_node(pop.pop[i].node) == 1)
      {
        fprintf(stderr,"bad pop at i=%d\n",i);
        exit(-1);
      }
    }
#endif


// --------- MAIN LOOP CREATING THE NEXT GENERATION ---------------
    i=0;
    while (i<popsize) /* fill pop_next. i is only incremented once new_ns.node is succesfully added, not discarded.  */
    {
      node_count = 0; /* count how many nodes are produced for new_ns.node for cost function*/

      // TREE MUTATIONS AND CROSSING using evolve

      new_flag = 1;
      existing_score = -1.0; 
      if ( rand() > PNEW*RAND_MAX)   
      {
        my_loc = i2point(i, pop.max_loc.x);

        i_tree1 = min(tournament2D(pop, my_loc, tour, compgridsize),pop.popsize-1);
        i_tree2 = min(tournament2D(pop, my_loc, tour, compgridsize),pop.popsize-1);

        tmpns1 = pop.pop[i_tree1];
        tmpns2 = pop.pop[i_tree2];

        tmpns3 = crossover_ns(tmpns1,tmpns2, PROBSWAP,1);  // crossover. this creates a new tree

        newns = allmut(tmpns3, mutationrate, Exp);  // mutate. this also creates a new tree

//    newns = copy_node_score(tmpns1);  // TESTING

#ifdef DEBUG
        if (check_node(newns.node) == 1)
        {
          fprintf(stderr,"newns not sound after unifcross and allmut, aborting. \n");
          exit(-1);
        }
#endif
        user_functions(newns);

#ifdef INLINESCORING
        i_tree_old = does_tree_exist(newns,pop, my_loc, compgridsize);
        if (i_tree_old > -1)
        {
          new_flag = 0;
          existing_score = pop.pop[i_tree_old].score;
          pop_next.trees_reused++;            
        } 
#endif

        free_node_score(tmpns3);
      }
      else
      {
        newns = makerandomns(MAXDEPTH,FPR,PPR, Exp);
      }

      size = tree_consts(newns.node, consts,nconsts,CONSTSCUTOFF); /* fill consts array */
      *nconsts=0;
 
 //     fprintf(stderr,"size=%d    \n",size);  

      if ( (size > -1) && (size < MAXSIZE) )
      {  
    /* perform optimization here */    

  //      mutconsts(consts,size, mutationrate/4.0); 

#ifdef INLINESCORING

        if (new_flag == 0)
        {
          error = existing_score;
        }
        else
        {

 #ifdef EVOLVEIC
          error = get_score(newns.node, newns.S_i, Exp);
 #else
          error = get_score(newns.node, Exp.S_i, Exp);
 #endif
        }
 
        if (error <1e18) 
        {

/*        parsimony not included in error reporting to report and elite output files */

          leafcount = conparcount(newns.node, current);
          *current = 0;

          // fill score vals in new generation (the old generation has old_score_vals).

          if (new_flag == 0)  // this indicates we are using cached score: no parsimony added.
          {
            //score_vals[i] 
            newns.score = error;          
          }
          else
          {
            //score_vals[i] = error + ((double) PARSIMONY*((double) leafcount)); 
            newns.score = error + ((double) PARSIMONY*((double) leafcount)); 
          }

          // fill tree array for new generation
          //trees2[i] = newns; /* found an acceptable tree */
          pop_next.pop[i] = newns;

          pop_next.totsize +=size;
          pop_next.toterror += error;
          if (error < pop_next.minerror)
          {
            pop_next.minerror = error;
            pop_next.i_min = i; // keep track of best score
          }
          i++; // increment loop counter
        }
        else
        {
          free_node_score(newns);  // tree not useable with error 1e19, discard
        }
#else

        pop_next.pop[i] = newns;

        pop_next.totsize +=size;

        i++; // increment loop counter

#endif  

      }
      else
      {
        free_node_score(newns); // tree not useable due to size limits, discard and try again
      }
    }; /* end while loop on popsize */


#ifndef INLINESCORING
// score all the new trees

    pop_next = score_pop(pop, pop_next, Exp, compgridsize);

#endif

  /* the next generation, trees2, is now ready
    copy score_vals to old_score_vals and copy trees2 to trees1 (to be improved later). 
  */

#ifdef EXPLICITMIG
    // exchange trees among neighboring cells internal to this CPU process.
    internal_migration(pop_next, tour, compgridsize, ((int) MIGRATE_INTERNAL*(compgridsize*compgridsize))/4  );
#endif

    varerror=0;
    avgerror = ((double) pop_next.toterror/pop_next.popsize);
    for (j=0; j < pop_next.popsize; j++)
    {
      varerror += (pop_next.pop[j].score - avgerror)*(pop_next.pop[j].score - avgerror);
    }

    varerror = varerror/pop_next.popsize;

 //   avghill = ((double) nhillsearches)/popsize;    

//    nhillsearches = 0;

    avgsize = ((double) pop_next.totsize)/pop_next.popsize;
    
    avgtreesreused = ((double) pop_next.trees_reused)/pop_next.popsize;  

    for (j=0;j<4;j++) /* output migrants to adjacent processes to files. j represents xm,xp,ym,yp */
    {
      migrants2buffer(pop_next, j, migrants, compgridsize, Exp, my_loc_external);
      
      sprintf(tmp_str,"%s_tmp",mig.outfiles+j*FNAMESIZE);
      store_data(tmp_str,treebuffer,"w"); /* write migrant tree files: xm,xp,ym,yp */
      rename(tmp_str,mig.outfiles+j*FNAMESIZE); // used tmp file as other process also read from these files
    }

    /* report stats */
    sprintf(line_str,"%g %g %g %g %g\n",pop_next.minerror, avgerror,sqrt(varerror),avgsize, avgtreesreused);
    store_data(repfile,line_str,"a"); /* append stats to report file (filename repfile) */

    (*code2str_table[(int) pop_next.pop[pop_next.i_min].node->op])(pop_next.pop[pop_next.i_min].node,tmp_str);

//    sprintf(line_str,"%g %s\n",minerror, tmp_str);

    line_str[0] = '\0';
    NodeScore2bufferline(pop_next.pop[pop_next.i_min],line_str);


    store_data(elitefile,line_str,"a"); /* append stats to report file */

    if (pop_next.pop[pop_next.i_min].score == 0.0)
    {
      printf("Found exactly matching solution: %sFinished.\n", line_str);
      exit(0);
    }

    pop = copy_pop(pop, pop_next);
    pop_next = reset_counters_pop(pop_next); // filling pop_next again

  }  /* end of runlen for loop */

  pop2buffer(pop_next);


  free_pop(pop);
  free(pop_next.pop);  // pop_next nodes have just been copied to pop so not using free_pop

  free_experiment(Exp);

  free(mig.infiles);
  free(mig.outfiles);  

  free(nconsts);
  free(current);

  store_data(restartfile,treebuffer,"w"); // buffer contains tree population, is written to file.
}



int main ( int argc, char *argv[] )
{

    int my_number;
    int qsubs;
    int runlen;
    int popsize; 
    int compgridsize; 
    double mutationrate; 
    int tour; 
 
    if ( argc != 9 ) /* argc should be 7 for correct execution */
    {
        printf( "usage: %s my_number qsubs runlen popsize compgridsize mutationrate tour\n", argv[0] );
    }
    else 
    {
        my_number = atoi(argv[1]);
        qsubs = atoi(argv[2]);
        runlen = atoi(argv[3]);
        popsize = atoi(argv[4]);
        compgridsize = atoi(argv[5]);
        mutationrate = strtod(argv[6],NULL);
        tour = atoi(argv[7]);
    
        printf("Starting my_number %d, qsubs %d, runlen %d, popsize %d, compgridsize %d, mutationrate %g, tour %d\n",my_number, qsubs , runlen, popsize, compgridsize, mutationrate, tour);

        c_nextgen(my_number, qsubs, runlen, popsize, compgridsize, mutationrate, tour);

    }

    return 0;
}


