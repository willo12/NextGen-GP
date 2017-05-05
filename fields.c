#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <compile_options.h>
#include <fields.h>

// should go into separate utilities file, or call fields.c fields_and_utilities.c

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


char* stradd(const char* a, const char* b){
/* from stack overflow  */

    size_t len = strlen(a) + strlen(b);
    char *ret = (char*)malloc(len * sizeof(char) + 1);
    *ret = '\0';
    return strcat(strcat(ret, a) ,b);
}



double rational_tanh(double x)
{

    if( x < -3 )
    {
        return -1;
    }
    else if( x > 3 )
    {
        return 1;
    }
    else
    {
        return x * ( 27 + x * x ) / ( 27 + 9 * x * x );
    }
} 


// ------- FIELD FUNCTIONS ------------------



void store_data(const char *filepath,const char *data, char *mode)
{
/* Write string buffer "data" to file. 


  add return value in future */

  FILE *fp;

  fp = fopen(filepath,mode);

  if (fp != NULL)
  {
    fprintf(fp,"%s",data);
    fclose(fp);
  }
  else
  {
    fprintf(stderr,"Could not open %s\n",filepath);
  }

}



Dims make_dims(int cols, int rows)
{
  Dims dims;

  dims.cols = cols;
  dims.rows = rows;

  return dims;
}

Dims bufferline2dims(char tmp_str[], char * infile)
{ // read cols, rows in that order from first line of file

  Dims dims;

  int cols, rows;

  const char hor_delim[2] = " ";
  char *token;  
  char *token_wrong;  
 
  cols = 0;
  rows = 0;

  dims = make_dims(0,0);

  token = strtok(tmp_str,hor_delim);   // start parsing

  if (*token != '#')
  { /* should start with # */
    fprintf(stderr,"Error: no header in file %s (wrong format?). \n",infile);   
    /* abort and trigger failure via ns values */
  }
  

  token = strtok(NULL,hor_delim);
  if (token != NULL) 
  {
    cols = atoi(token);
  }
  else
  { /* if there's no tree after the line with a score, trigger failure */
    fprintf(stderr,"Error in reading file %s (wrong format?). \n",infile);   
    /* abort and trigger failure via ns values */
       
  }


  if (cols == 0) // indicates string to int conversion failed
  {
    //fprintf(stderr,"Error reading file: could not convert cols int \n");
    return dims; /* abort and trigger failure via F values */          
  }

  token = strtok(NULL,hor_delim);
  if (token != NULL) 
  {
    token_wrong = strtok(NULL,hor_delim);
    if (token_wrong != NULL)
    {
      fprintf(stderr,"Error reading file: to many items on line \n");
      return dims; /* abort and trigger failure via ns values */      
    }

    // file reading checks
    if (( token[0] == '\n' ) || ( token[0] == '\0' )) // empty line? this can't be an empty line
    {
      fprintf(stderr,"Error reading file: empty line \n");
      return dims; /* abort and trigger failure via ns values */
    }

   // if ( (pos=strchr(token,'\n')) != NULL)
   //   *pos = '\0';

    rows = atoi(token);

    if (rows == 0) // indicates string to int conversion failed
    {
      //fprintf(stderr,"Error reading file: could not convert rows int \n");
      return dims; /* abort and trigger failure via F values */          
    }

    dims.cols = cols;
    dims.rows = rows;

  }
  else
  { /* if there's no tree after the line with a score, trigger failure */
    fprintf(stderr,"Error in reading file %s (wrong format?). \n",infile);   
    /* abort and trigger failure via ns values */
       
  }

  return dims;
}



Field bufferline2doubles(char tmp_str[], char * infile, Field F, int i)
{

  int j; // col
  double value;

  const char hor_delim[2] = " ";
  char *token;  

  token = strtok(tmp_str,hor_delim);   

  j = 0;
  while (( token != NULL ) && (j<F.dims.cols))
  {
    value = strtod(token, NULL);
    
    *(F.data + i*F.dims.cols + j  )  = value;
    token = strtok(NULL, hor_delim);
    j++;
  }
   
  return F;
}


FieldInt bufferline2ints(char tmp_str[], char * infile, FieldInt F, int i)
{
  int j; // col
  double value;

  const char hor_delim[2] = " ";
  char *token;  

  token = strtok(tmp_str,hor_delim);   

  j = 0;
  while (( token != NULL ) && (j<F.dims.cols))
  {
    value = ( (long int) strtod(token, NULL) );
    
    *(F.data + i*F.dims.cols + j  )  = value;
    token = strtok(NULL, hor_delim);
    j++;
  }
   
  return F;
}




FieldInt bufferline2config(char tmp_str[], char * infile, FieldInt F, int i)
{
  int j; // col
  double value;

  const char hor_delim[2] = " ";
  char *token;  

  token = strtok(tmp_str,hor_delim);   

  j = 0;
  while (( token != NULL ) && (j<F.dims.cols))
  {
    value = ( (long int) strtod(token, NULL) );
    
    *(F.data + i*F.dims.cols + j  )  = value;
    token = strtok(NULL, hor_delim);
    j++;
  }
   
  return F;
}




Field make_field(Dims dims)
{
  Field data;

  data.dims = dims;

  if (dims.cols*dims.rows > 0)
  {
    data.data = (double *)malloc(sizeof(double)*(dims.cols*(dims.rows+1) )); 
  }
  else
  {
    data.data =NULL;
  }
  return data;
}

Field make_field_empty(Dims dims)
{
  Field data;

  data.dims = dims;

  data.data =NULL;
  
  return data;
}

FieldInt make_field_int(Dims dims)
{
  FieldInt data;

  data.dims = dims;

  if (dims.cols*dims.rows > 0)
  {
    data.data = (long int *)malloc(sizeof(long int)*(dims.cols*(dims.rows+1) )); 
  }
  else
  {
    data.data =NULL;
  }
  return data;
}

FieldInt make_field_int_empty(Dims dims)
{
  FieldInt data;

  data.dims = dims;

  data.data =NULL;
  
  return data;
}

int free_field(Field F)
{
  if (F.data != NULL)
  {
    free(F.data);
    F.data = NULL;
  }

  return 0;
}

int free_int_field(FieldInt F)
{
  if (F.data != NULL)
  {
    free(F.data);
    F.data = NULL;
  }

  return 0;
}


Field read_array(char *filepath)
{
  int i,j;
  FILE *fp;

  Dims dims;
  Field F;

  //Field F = make_field(-1,-1);

  char tmp_str[MAXDATALINESTR];

  i=0; // counts the trees

  fp = fopen(filepath,"r");  // file pointer fp
  if (fp != NULL)
  {
    if (fgets (tmp_str,MAXDATALINESTR,fp)!=NULL )
    {
      dims = bufferline2dims(tmp_str, filepath); // add cols and rows to F
      if (dims.cols*dims.rows == 0)
      {
        return make_field(dims); // return field with NULL data 
      }
    }
    else
    {
      fprintf(stderr,"Not enough lines in %s. \n",filepath);
      return make_field(make_dims(0,0));
    }

    F = make_field(dims);
 
    while((fgets (tmp_str,MAXDATALINESTR,fp)!=NULL ) && (i<F.dims.rows)  )  // read line by line from file
    {
      F = bufferline2doubles(tmp_str, filepath, F, i); // insert line into F.data 

//      if (F.dims.cols == 0) // signal 
//      {
//        break;
//      }

      i++;
    }

    for (j=0;j<F.dims.cols;j++)
    {
      *(F.data+F.dims.cols*i+j) = *(F.data+F.dims.cols*(i-1)+j);
    }
   

    fclose(fp);
 
    if (i<F.dims.rows)
    {
      fprintf(stderr,"Not enough lines in %s. \n",filepath);
      F.dims.rows = 0;
      return F;
    }

    return F;
  }
  else
  {
    fprintf(stderr,"Error: could not open %s \n",filepath);
  }  


  return make_field(make_dims(0,0)); /* can't open file. */  

}

FieldInt read_array_int(char *filepath)
{
  int i;
  FILE *fp;

  Dims dims;
  FieldInt F;

  //Field F = make_field(-1,-1);

  char tmp_str[MAXDATALINESTR];

  i=0; // counts the trees

  fp = fopen(filepath,"r");  // file pointer fp
  if (fp != NULL)
  {
    if (fgets (tmp_str,MAXDATALINESTR,fp)!=NULL )
    {
      dims = bufferline2dims(tmp_str, filepath); // add cols and rows to F
      if (dims.cols*dims.rows == 0)
      {
        return make_field_int(dims); // return field with NULL data 
      }
    }
    else
    {
      fprintf(stderr,"Not enough lines in %s. \n",filepath);
      return make_field_int(make_dims(0,0));
    }

    F = make_field_int(dims);
 
    while((fgets (tmp_str,MAXDATALINESTR,fp)!=NULL ) && (i<F.dims.rows)  )  // read line by line from file
    {
      F = bufferline2ints(tmp_str, filepath, F, i); // insert line into F.data 

//      if (F.dims.cols == 0) // signal 
//      {
//        break;
//      }

      i++;
    }

    fclose(fp);
 
    if (i<F.dims.rows)
    {
      fprintf(stderr,"Not enough lines in %s. \n",filepath);
      F.dims.rows = 0;
      return F;
    }

    return F;
  }
  else
  {
    fprintf(stderr,"Error: could not open %s \n",filepath);
  }  
  
  return make_field_int(make_dims(0,0)); /* can't open file. */  

}



void copy_field(Field F_from, Field F_to)
{
  int i,j, rows, cols;

  rows = min(F_from.dims.rows,F_to.dims.rows);
  cols = min(F_from.dims.cols,F_to.dims.cols);

  for (i=0;i<rows;i++)
  {
    for (j=0;j<cols;j++)
    {
      *(F_to.data + i*F_to.dims.cols + j  ) = *(F_from.data + i*F_from.dims.cols + j  );
    }
  }

}

int write_array(Field F, char *filepath)
{

  int i, j;
  char tmp_str[MAXDATALINESTR];
  char *buffer; // string buffer for output

  int max_cell = 30; // max width of a printed cell (double). Make this big enough for most print cases.
  int buff_elts = (max_cell+1)*F.dims.cols*(F.dims.rows+1); // total amount of chars in the string buffer

//  fprintf(stderr, "%d bytes. \n", buff_elts );

  buffer = (char *)malloc( buff_elts*sizeof(char) ); 



  sprintf(tmp_str, "# %d %d\n",F.dims.cols,F.dims.rows); // write dimensions
  strcpy(buffer, tmp_str);

  for (i=0;i<F.dims.rows;i++)
  {
    for (j=0;j<F.dims.cols;j++)
    {
      sprintf(tmp_str, "%0.8g", *(F.data + i*F.dims.cols + j  ) );
      strcat(buffer, tmp_str);

      if (j < F.dims.cols-1)
      {
        sprintf(tmp_str," ");
        strcat(buffer, tmp_str);
      }

    }
    sprintf(tmp_str,"\n");
    strcat(buffer, tmp_str);

  }

  store_data(filepath,buffer,"w"); // buffer contains tree population, is written to file.

  free(buffer);
   
  return(0);
}


Field check_field(Field F)
{
  /* Check whether field is empty. If so, exit program.

  */

  if ( (F.data == NULL) || ( (F.dims.rows == 0) && (F.dims.cols == 0)   )   )
  {
    fprintf(stderr,"Exiting program ...\n");
    exit(-1);
  }


  return F;
}

FieldInt check_int_field(FieldInt F)
{
  /* Check whether field is empty. If so, exit program.

  */

  if ( (F.data == NULL) || ( (F.dims.rows == 0) && (F.dims.cols == 0)   )   )
  {
    fprintf(stderr,"Exiting program ...\n");
    exit(-1);
  }


  return F;
}
