include ../../example.mk

CC=mpic++

LDIR =

OPT=
OBJ_VIC_PETSC = main_vic_petsc.o
OBJ_VIC_PETSC_OPT = main_vic_petsc_opt.o

all: vic_petsc vic_petsc_opt

vic_petsc_test: OPT += -DTEST_RUN
vic_petsc_test: vic_petsc

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

vic_petsc_opt: $(OBJ_VIC_PETSC_OPT)
	$(CC) -o $@ $^ $(LIBS_PATH) $(LIBS)

vic_petsc: $(OBJ_VIC_PETSC)
	$(CC) -o $@ $^ $(LIBS_PATH) $(LIBS)
	
run: vic_petsc_test
	mpirun -np 4 ./vic_petsc
	
.PHONY: clean all

clean:
	rm -f *.o *~ core vic

