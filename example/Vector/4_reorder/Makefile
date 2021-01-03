include ../../example.mk

CC=mpic++

LDIR =

OBJ_DORD = main_data_ord.o
OBJ_CORD = main_comp_ord.o
OPT =

all: md_data_ord md_comp_ord
all_test: OPT += -DTEST_RUN
all_test: md_data_ord_test md_comp_ord_test

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++14 $(OPT) -o $@ $< $(INCLUDE_PATH)

md_data_ord: $(OBJ_DORD)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

md_comp_ord: $(OBJ_CORD)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

md_data_ord_test: $(OBJ_DORD)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

md_comp_ord_test: $(OBJ_CORD)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all_test
	mpirun -np 4 ./md_data_ord_test && mpirun -np 4 ./md_comp_ord_test

.PHONY: clean all run all_test on_test

clean:
	rm -f *.o *~ core md_data_ord md_comp_ord

