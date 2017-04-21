#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <compile_options.h>
#include <states.h>
#include <node_structures.h>
#include <fields.h>

#include <model.h>


int free_state(State S)
{
  if (S.data == NULL)
  {
    fprintf(stderr,"Warning, free_state trying to free NULL pointer. \n");
    return -1;
  }

  free(S.data);
  S.data = NULL;

  return 0;
}


State make_state(int size)
{
  State S;

  S.size = size;
#ifdef INTSTATES
  S.data = (int *)malloc(sizeof(int)*size);
#else
  S.data = (double *)malloc(sizeof(double)*size);
#endif

  return S;
}

#ifdef INTSTATES
State fill_state(State S, int *A)
#else
State fill_state(State S, double *A)
#endif
{
  int i;

  if (A == NULL)
  {
    for (i=0;i<SPACEDIM;i++) // SPACEDIM instead of S.size for speed
    {
#ifdef INTSTATES
      (*(S.data+i)) = 0;
#else
      (*(S.data+i)) = 0.0;
#endif
    }
  }
  else
  {
    for (i=0;i<SPACEDIM;i++)
    {
      (*(S.data+i)) = *(A+i);
    }
  }

  return S;
}

#ifdef INTSTATES
State init_state(int size, int *A)
#else
State init_state(int size, double *A)
#endif
{
/* allocate memory for state and fill data with zeros */

  State S = make_state(size);

  return fill_state(S, A);

}

//State make_copy_state(State S)
//{
//  return init_state(SPACEDIM, S.data); 
//}


State copy_state(State S_from, State S_to)
{
  int i;

  if (S_from.size != S_to.size)
  {
    fprintf(stderr,"Error copying states: sizes differ. %d vs %d \n", S_from.size , S_to.size);
    exit(-1);
    return S_from; 
  }
  else
  {

    for (i=0;i<S_from.size;i++)
    {
      (*(S_to.data+i)) = (*(S_from.data+i));
    }  

    return S_to;
  }
}

int compare_state(State S1, State S2)
{
  int i;

  for (i=0;i<SPACEDIM;i++)
  {
    if ( S1.data[i] != S2.data[i]   )  
    {
      return 1;
    }
  }

  return 0;

}


#ifdef INTSTATES

State makerandomstate(State S_i)
{
  int i;
 
  State S = make_state(SPACEDIM);

  for (i=0;i<SPACEDIM;i++)
  {
    (*(S.data+i)) = rand()%MAXINT;
  }
 
  return S;
}

#else

State makerandomstate(State S_i)
{
  int i;
  double pert;

  double ampl = 0.5;

  State S = make_state(SPACEDIM);

  for (i=0;i<SPACEDIM;i++)
  {
    pert = (( (double) (2*rand() -RAND_MAX) )*ampl/RAND_MAX);

    (*(S.data+i)) = trunc(   (pert + (*(S_i.data+i)) )*DOUBLETRUNC  )/DOUBLETRUNC;    
  }
 
  return S;
}

#endif

State crossover_states(State S1, State S2)
{
// for SPACEDIM 4: (a,b,c,d)
// division point at index 1,2,3

  int i;
  int divisor;

  State S = make_state(SPACEDIM);

  divisor = rand()%(SPACEDIM+1);

  for (i=0;i<divisor;i++)
  {
    (*(S.data+i)) = (*(S1.data+i));
  }

  for (i=divisor;i<SPACEDIM;i++)
  {
    (*(S.data+i)) = (*(S2.data+i));
  }

  return S;
}



State mut_state(State S)
{
  int i;
  int rnd;

  State new_S = make_state(SPACEDIM);

  for (i=0;i<SPACEDIM;i++)
  {
    rnd = rand()%5;

    if (rnd < 3)
    {
#ifdef INTSTATES
      new_S.data[i] = rand()%MAXINT;
#else
      new_S.data[i] = S.data[i] + ((double) (rnd-1))*WANDERSTEP;
#endif

    }

    else
    {
      new_S.data[i] = S.data[i];
    }
  }
  return new_S;
}


int state2bufferline(char treebuffer[], State S )
{
/* write out state to string for file io

*/
  char tmp_str[MAXTREESTR];
  int i;

  sprintf(tmp_str,"(");
  strcat(treebuffer,  tmp_str ); 

  for (i=0;i<SPACEDIM;i++)
  {

#ifdef INTSTATES
    sprintf(tmp_str,"%d", (*(S.data+i)) );
#else
    sprintf(tmp_str,"%.7g", (*(S.data+i)) );
#endif

    strcat(treebuffer,  tmp_str );

    if (i<SPACEDIM-1)
    {
      sprintf(tmp_str,",");
      strcat(treebuffer,  tmp_str ); 
    }
  }

  sprintf(tmp_str,")");
  strcat(treebuffer,  tmp_str ); 

  return 0;
}



State str2state(char * text)
{

  int i;
  char * p_end;
  char * p_next;
  char * tmp_str = (char *)malloc(sizeof(char)*(MAXTREESTR));
  State S = make_state(SPACEDIM);
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
    fprintf(stderr,"Problem reading IC. No opening bracket. \n"); 
    free_state(S);
    S.data = NULL;
    return S;
  }

  p_next = tmp_str + 1;

  /* get the first token */
  token = strtok_r(p_next, s, &saveptr);
   
  /* walk through other tokens */
  i=0;
  while( (token != NULL) && (end_flag == 0) && (i<SPACEDIM) ) 
  {

    if ( (p_end=strchr(token, end_ch)) != NULL  )
    {
      end_flag = 1;
      *p_end = '\0';
    }

    (*(S.data+i)) = strtod(token,NULL);
    
    token = strtok_r(NULL, s, &saveptr);

    i++;
  }

  if ( (i < SPACEDIM) || (end_flag == 0) )
  {
    fprintf(stderr,"Problem reading IC. Saved data read up to index %d. Also check SPACEDIM. \n",i);
    free_state(S);
    S.data = NULL;
    return S;    
  }

  free(tmp_str);

  return S;
}










