include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o
gray_scott_sparse_cs_test: OPT += -DTEST_RUN
gray_scott_sparse_cs_test: gray_scott_sparse_cs

%.o: %.cpp
	$(CC) -O3 -g $(OPT) -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

gray_scott_sparse_cs: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: gray_scott_sparse_cs

run: gray_scott_sparse_cs_test
	mpirun --oversubscribe -np 4 ./gray_scott_sparse_cs

.PHONY: clean all run

clean:
	rm -f *.o *~ core gray_scott_sparse_cs

