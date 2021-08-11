include ../../example.mk

CC=mpic++

LDIR =
OPT=

OBJ = main.o

cp_circle:
cp_circle_test: OPT += -DTEST_RUN
cp_circle_test: cp_circle

%.o: %.cpp
	$(CC) -O3 -g  $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

cp_circle: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: cp_circle

run: cp_circle_test
	mpirun --oversubscribe -np 2 ./cp_circle

.PHONY: clean all run

clean:
	rm -f *.o *~ core cp_circle

