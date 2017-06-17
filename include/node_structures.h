
#ifndef NODESTRUCTH
#define NODESTRUCTH

#define MAXCHILD 3


typedef struct point {
  int x;
  int y;
} Point;

typedef struct rect {
  Point min;
  Point max;
} Rect;


typedef struct node {
  int op;
  void * children[MAXCHILD];
} Node;


typedef struct run_node {
  int op;
  void * children[MAXCHILD];

  double (*func)(Node *tree);

} RunNode; // node executed during tree runtime


typedef struct node_score {
  double score;
  Node * node;
#ifdef EVOLVEIC
  State S_i;
#endif

#if SCALARPARS > 0
  ScalarPars ip;
#endif

} NodeScore;


typedef struct population {
  // general
  NodeScore * pop;
  int popsize;
  Point max_loc;

  // specific counters
 
  int i; // if you want to track index in pop
  Point my_loc; 

  int trees_reused;
  int totsize;
  double toterror; // total error
  double minerror; // keep track of best score
  int i_min; // keep track of index of best score

} Population;

#endif
