include ../../common.mk
include ../../example.mk

OBJ = main.o
gray_scott_sparse_test: OPT += -DTEST_RUN
gray_scott_sparse_test: gray_scott_sparse

%.o: %.cpp
	$(CC) -O3 -g $(OPT) -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

gray_scott_sparse: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: gray_scott_sparse

run: gray_scott_sparse_test
	mpirun --oversubscribe -np 4 ./gray_scott_sparse

.PHONY: clean all run

clean:
	rm -f *.o *~ core gray_scott_sparse

