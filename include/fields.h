#ifndef FIELDSH

#define MAXDATALINESTR 400

typedef struct dims {

  int cols;
  int rows;

} Dims;

typedef struct field {

  double* data;
  Dims dims;

} Field;

typedef struct field_int {

  long int* data;
  Dims dims;

} FieldInt;

int min(int x, int y);
int max(int x, int y);

char* stradd(const char* a, const char* b);
double rational_tanh(double x);

void store_data(const char *filepath,const char *data, char *mode);
Dims make_dims(int cols, int rows);
Dims bufferline2dims(char tmp_str[], char * infile);
Field bufferline2doubles(char tmp_str[], char * infile, Field F, int i);
FieldInt bufferline2ints(char tmp_str[], char * infile, FieldInt F, int i);
FieldInt bufferline2config(char tmp_str[], char * infile, FieldInt F, int i);
Field make_field(Dims dims);
Field make_field_empty(Dims dims);
FieldInt make_field_int(Dims dims);
FieldInt make_field_int_empty(Dims dims);
int free_field(Field F);
int free_int_field(FieldInt F);
Field read_array(char *filepath);
FieldInt read_array_int(char *filepath);
void copy_field(Field F_from, Field F_to);
int write_array(Field F, char *filepath);

Field check_field(Field F);
FieldInt check_int_field(FieldInt F);

#endif
