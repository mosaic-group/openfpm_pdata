include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o
OBJ_EXPR = main_expr.o
OBJ_VL = main_vl.o
OBJ_VL_SYM = main_vl_sym.o

all: md_dyn md_dyn_expr md_dyn_vl md_dyn_vl_sym

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

md_dyn: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

md_dyn_expr: $(OBJ_EXPR)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

md_dyn_vl: $(OBJ_VL)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

md_dyn_vl_sym: $(OBJ_VL_SYM)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all
	source $$HOME/openfpm_vars; mpirun -np 3 ./md_dyn; mpirun -np 3 ./md_dyn_expr; mpirun -np 3 ./md_dyn_vl; mpirun -np 3 ./md_dyn_vl_sym;

.PHONY: clean all run

clean:
	rm -f *.o *~ core md_dyn md_dyn_expr md_dyn_vl md_dyn_vl_sym

