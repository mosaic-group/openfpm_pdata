include ../../example.mk

CC=mpic++

LDIR =
OPT=-I/Users/lschulze/gsl/include/

OBJ = main.o

cp_ellipse:
cp_ellipse_test: OPT += -DTEST_RUN
cp_ellipse_test: cp_ellipse

%.o: %.cpp
	$(CC) -O3 -g  $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

cp_ellipse: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: cp_ellipse

run: cp_ellipse_test
	mpirun --oversubscribe -np 2 ./cp_ellipse

.PHONY: clean all run

clean:
	rm -f *.o *~ core cp_ellipse
