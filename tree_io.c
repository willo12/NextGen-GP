
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



Node *str2node(char text[], char left_bracket, char delim)
{

  /* recursively fill node from text containing tree in bracket notation */

  int j;
  int l,i;
  int I[4];

  char op;

  Node *newnode;

  j=0;
  // try all the leaves
  while ( str2leaf_table[j] != NULL  ) // scan through the leaf nodes to see whether text indicates leaf node
  {
    if ( (newnode=(*str2leaf_table[j])(text, left_bracket, delim)) != NULL)
    {
      return newnode;
    }
    j++;
  }

  if (text[0] == left_bracket) /* this is not a leaf node */
  {
    l=strlen(text);
    op=text[l-1];  /* found op. */

    text[l-2]='\0'; /* remove trailing op and ) bracket */
    text = text+1; /* remove beginning ( bracket */

/* split string on ',' except for comma's within bracketed parts.
split is achieved by replacing those ',' with \0 and maintaining
indices into string in I, indicating where the pieces start. */
    split(text,I, left_bracket, delim);

    newnode=talloc();
    newnode->op = op;
    /* fill children of newnode, using str2node recursively */
    newnode->children[0] = str2node(text, left_bracket, delim);    /* first child at start of string (\0 have now been inserted) */

    /* collect remaining children */
    i=0;
    while (I[i] !=0)
    {
      newnode->children[i+1] = str2node(text+I[i], left_bracket, delim);    
      i++;
    }
    return newnode;

  }

  return NULL; // there is a problem
}





Node *str2par_fun_node(char text[], char left_bracket, char delim)
{
  int l;
  int I[4];

  char op;

  l=strlen(text);
  op=text[l-1];  /* found op.  */

    if (op == 'H' )
    {
      text[l-2]='\0'; /* remove trailing op and ) bracket */
      text = text+1; /* remove beginning ( bracket */

      split(text,I, left_bracket, delim);

      return make_par_fun_node('H', atoi(text) , strtod(text+I[0],NULL) );

    }

  return NULL;
}


Node *str2par_node(char text[], char left_bracket, char delim)
{
  if (text[0]=='p') /* it's a parameter */
  {
    return make_par_node(par_op_char,atoi(text+1));
  }
  else
  {
    return NULL;
  }
}

Node *str2const_node_int(char text[], char left_bracket, char delim)
{
  if ( isdigit(text[0]) || (text[0]=='-') )
  {
    return make_const_node_int(const_op_char, atoi(text));
  }
  else
  {
    return NULL;
  }
}


Node *str2const_node(char text[], char left_bracket, char delim)
{
  if ( isdigit(text[0]) || (text[0]=='-') )
  {
    return make_const_node(const_op_char, strtod(text,NULL));
  }
  else
  {
    return NULL;
  }
}


char* constChrInt(Node *tree, char trstr[])
{
/* convert tree to str where tree is a constant leaf */
  sprintf(trstr,"%d",*((int *) tree->children[0]));

  return trstr; 
};


/* Functions ending in Chr print node to string. Used in code2str_table. 
   Names correspond to node names.
*/

char* nopChr(Node *tree, char trstr[])
{
  return trstr;
};

char* constChr(Node *tree, char trstr[])
{
/* convert tree to str where tree is a constant leaf */
  sprintf(trstr,"%0.8g",*((double *) tree->children[0]));

  return trstr; 
};


char* parChr(Node *tree, char trstr[])
{
/* convert tree to str where tree is a par leaf */

  sprintf(trstr,"p%d",*((int *) tree->children[0]));

  return trstr; 
};


char* parfunChr(Node *tree, char trstr[])
{
/* convert tree to str where tree is a par leaf */

  sprintf(trstr,"(%d,%0.8g)H",*((int *) tree->children[0])  ,  *((double *) tree->children[1])  );

  return trstr; 
};



char* nFunChr(Node *tree, char trstr[])
{
/* convert tree to str recursively, where root node is an op that takes n arguments */
 
  int i, n, op, child_op;
  char tmp_str[MAXTREESTR];
  Node * child;

  op = ((int) ((Node *) tree)->op);
  n = arg_table[op];

  strcpy(trstr,"(");

  for (i=0;i<n;i++){  
    child = ((Node *) tree->children[i]);
    child_op = ((int) child->op);

    (*code2str_table[child_op])(child, tmp_str );
    strcat(trstr,tmp_str);

    if (i<n-1) /* join on ',' delimiter */
    {
      strcat(trstr,",");
    }
  } 

  sprintf(tmp_str,")%c", tree->op);
  strcat(trstr, tmp_str);

  return trstr;

};

char* parChr_c(Node *tree, char trstr[])
{
/* convert tree to str where tree is a par leaf */

  sprintf(trstr,"p%d",*((int *) tree->children[0]));


  int i = *((int *) tree->children[0]);

  if (i<SPACEDIM)
  {
    sprintf(trstr,"S.data[%d]",i);    
  }
  else
  {
    sprintf(trstr,"Exp.Iffs.data[steps_forc*(FORCING_TERMS+1)+1+%d]",i);    
  }

  return trstr; 
};


char* nFunChr_c(Node *tree, char trstr[])
{
/* convert tree to str recursively, where root node is an op that takes n arguments */
 
  double tmp_double;
  int child_op;
  char tmp_str[MAXTREESTR];
  Node * child;

  if ( ((Node *) tree)->op == 'T' )
  {
    strcpy(trstr,"rth(");

    child = ((Node *) tree->children[0]);
    child_op = ((int) child->op);
    (*code2c_str_table[child_op])(child, tmp_str );
    strcat(trstr,tmp_str);
    sprintf(tmp_str,")" );
    strcat(trstr, tmp_str);
  }
  else
  {
    strcpy(trstr,"(");
 
    child = ((Node *) tree->children[0]);
    child_op = ((int) child->op);
    (*code2c_str_table[child_op])(child, tmp_str );
    strcat(trstr,tmp_str);

    child = ((Node *) tree->children[1]);
    child_op = ((int) child->op);

    if ( ((Node *) tree)->op == 'A' )
    {
      strcat(trstr,"+");
      (*code2c_str_table[child_op])(child, tmp_str );
    }
    else if ( ((Node *) tree)->op == 'S' )
    {
      if (child_op == const_op_char)
      {

        tmp_double = *((double *) child->children[0]);
 //       fprintf(stderr,"yo %g",tmp_double);

        if (tmp_double<0)        
        {
          strcat(trstr,"+");
          sprintf(tmp_str,"%0.8g",-tmp_double);
        }
        else if (tmp_double>0)
        {
          strcat(trstr,"-");
          sprintf(tmp_str,"%0.8g",tmp_double);
        }
        else
        {
          strcat(trstr,"+");
          sprintf(tmp_str,"%0.8g",0.0);
        }

      }
      else
      { 
        strcat(trstr,"-");
        (*code2c_str_table[child_op])(child, tmp_str );
      }
    }
    else if ( ((Node *) tree)->op == 'M' )
    {
      strcat(trstr,"*");
      (*code2c_str_table[child_op])(child, tmp_str );
    }

    strcat(trstr,tmp_str);

    sprintf(tmp_str,")" );
    strcat(trstr, tmp_str);
  }

  return trstr;

};






