include ../../example.mk

CC=mpic++

LDIR =
OPT=

OBJ = main.o

cp_inter:
cp_inter_test: OPT += -DTEST_RUN
cp_inter_test: cp_inter

%.o: %.cpp
	$(CC) -O3 -g  $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

cp_inter: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: cp_inter

run: cp_inter_test
	mpirun --oversubscribe -np 2 ./cp_inter

.PHONY: clean all run

clean:
	rm -f *.o *~ core cp_inter

