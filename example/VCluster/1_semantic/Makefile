include ../../example.mk
include ../../common.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

vcluster: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: vcluster

run: all
	mpirun --oversubscribe -np 3 ./vcluster

.PHONY: clean all run

clean:
	rm -f *.o *~ core vcluster

