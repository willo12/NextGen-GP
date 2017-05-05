
#ifndef COMPOP_H
#define COMPOP_H

#define COMPARTMENTS 4

#define TANH_ONCE
#define INTPARS 0
#define EVOLVEIC

//#define CUSTOM_NODES
#define TS_FACTOR 2
#define FORCING_TERMS 2
#define PARSIMONY 0.0 // to calculate effect of tree size on score
#define DOUBLETRUNC 10000
#define MAXTREESTR 10000  // max size tree string
#define OPTABLE 200 // size optable
#define MAXPAR 40 // max number of registers
#define STACKSIZE 2000

#define RESTRICTED_PARS 0

#define MARKER 0
#define NUMFUN 0
//#define CURTAILFORCING
#define INLINESCORING
//#define INTCONSTS
//#define INTSTATES
//#define MACHINE_FUNCTIONS

//#define DEBUG

//#define SUPERLEAF


#endif
