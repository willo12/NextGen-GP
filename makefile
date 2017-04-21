IDIR =include
CC=gcc
CFLAGS=-I$(IDIR) -Wall -O3

ODIR=obj

_DEPS = c_nextgen.h model.h fields.h node_structures.h compile_options.h states.h basic_ops.h my_ops.h tree_io.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = model.o c_nextgen.o fields.o states.o basic_ops.o my_ops.o tree_io.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))



$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

c_nextgen: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) -lm -Wall -O3

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
