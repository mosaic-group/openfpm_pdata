include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o
OBJ_EXPR = main_expr.o

all: md_dyn md_dyn_expr

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

md_dyn: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

md_dyn_expr: $(OBJ_EXPR)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all
	source $$HOME/openfpm_vars; mpirun -np 3 ./md_dyn; mpirun -np 3 ./md_dyn_expr

.PHONY: clean all run

clean:
	rm -f *.o *~ core md_dyn md_dyn_expr

