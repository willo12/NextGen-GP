
#ifndef COMPOP_H
#define COMPOP_H

#define SPACEDIM 2   // spatial dimension of dynamical system
#define COMPARTMENTS 4

#define I_START_POLY 0
//#define DOUBLEPARS
#define SCALARPARS 0

#define USE_RK4
#define TANH_ONCE
#define OBSCOLS 1
#define EVOLVEIC
#define OBSERROR
//#define COMPILE

//#define CUSTOM_NODES
#define TS_FACTOR 2
#define FORCING_TERMS 1
#define PARSIMONY 0.04 // to calculate effect of tree size on score
#define DOUBLETRUNC 10000
#define MAXTREESTR 10000  // max size tree string
#define MAXTEMPLATESIZE 10000
#define OPTABLE 200 // size optable
#define MAXPAR 40 // max number of registers
#define STACKSIZE 2000

#define MAXHEIGHT 8
#define MAXHEIGHT2 3

#define RESTRICTED_PARS 0

#define MARKER 0
#define NUMFUN 0
#define CURTAILFORCING
#define INLINESCORING

//#define INTCONSTS
//#define INTSTATES
//#define MACHINE_FUNCTIONS

//#define DEBUG

//#define SUPERLEAF


#endif
