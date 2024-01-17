include ../../../common.mk
include ../../../example.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

diff_1d: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: diff_1d

run: all
	mpirun --oversubscribe -np 4 ./diff_1d

.PHONY: clean all

clean:
	rm -f *.o *~ core diff_1d

