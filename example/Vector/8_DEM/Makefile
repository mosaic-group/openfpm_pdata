include ../../example.mk
include ../../common.mk

OBJ = main.o

all: dem

dem_test: OPT += -DTEST_RUN
dem_test: all

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

dem: $(OBJ)
	$(CC) -o $@ $^ $(OPT) $(CFLAGS) $(LIBS_PATH) $(LIBS)


run: dem_test
	mpirun --oversubscribe -np 2 ./dem

.PHONY: clean all run

clean:
	rm -f *.o *~ core dem

