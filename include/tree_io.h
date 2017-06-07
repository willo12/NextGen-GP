#ifndef TREEIO_H
#define TREEIO_H


Node *str2node(char text[], char left_bracket, char delim);

Node *str2par_node(char text[], char left_bracket, char delim);
Node *str2const_node_int(char text[], char left_bracket, char delim);
Node *str2const_node(char text[], char left_bracket, char delim);
Node *str2par_fun_node(char text[], char left_bracket, char delim);

char* constChrInt(Node *tree, char trstr[]);
char* nopChr(Node *tree, char trstr[]);
char* constChr(Node *tree, char trstr[]);
char* parChr(Node *tree, char trstr[]);
char* parfunChr(Node *tree, char trstr[]);
char* nFunChr(Node *tree, char trstr[]);

char* parChr_c(Node *tree, char trstr[]);
char* nFunChr_c(Node *tree, char trstr[]);

char* zeroFunChrSL(Node *tree, char trstr[]);
char* oneFunChrSL(Node *tree, char trstr[]);
char* twoFunChrSL(Node *tree, char trstr[]);
char* threeFunChrSL(Node *tree, char trstr[]);

#endif
