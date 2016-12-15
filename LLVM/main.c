
#include <stdio.h>
#include "c_nextgen.h"
 
int main(void)
{
    puts("This is a shared library test...");
    test_llvm();
    return 0;
}
