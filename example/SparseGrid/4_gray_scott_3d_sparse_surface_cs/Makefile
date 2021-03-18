include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o
gray_scott_sparse_cs_surface: OPT += -DTEST_RUN
gray_scott_sparse_cs_surface: gray_scott_sparse_surface_cs

%.o: %.cpp
	$(CC) -O3 -g $(OPT) -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

gray_scott_sparse_surface_cs: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: gray_scott_sparse_surface_cs

run: gray_scott_sparse_cs_surface
	mpirun -np 4 ./gray_scott_sparse_surface_cs

.PHONY: clean all run

clean:
	rm -f *.o *~ core gray_scott_sparse_surface_cs

