
INC_DIR=./includes


CC=clang
CFLAGS=-g `/usr/bin/llvm-config-3.5 --cflags` -fPIC -Wall -Wextra -O2 -I$(INC_DIR)
DEPS=c_nextgen.h

LD=clang++
LDFLAGS=`/usr/bin/llvm-config-3.5 --cxxflags --ldflags --libs core executionengine jit interpreter analysis native bitwriter --system-libs`

#all: libc_nextgen.so
all: c_nextgen


c_nextgen: c_nextgen.o
	$(LD) $< $(LDFLAGS) -o $@

c_nextgen.o: c_nextgen.c
	$(CC) $(CFLAGS) -c $<

libc_nextgen.so: c_nextgen.o
	$(LD) $< $(LDFLAGS) -shared -o $@

c_nextgen.bc: c_nextgen
	./c_nextgen 0 0

c_nextgen.ll: c_nextgen.bc
	llvm-dis $<

clean:
	-rm -f c_nextgen.o c_nextgen c_nextgen.bc c_nextgen.ll
