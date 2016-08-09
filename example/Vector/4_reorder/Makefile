include ../../example.mk

CC=mpic++

LDIR =

OBJ_DORD = main_data_ord.o
OBJ_CORD = main_comp_ord.o

all: md_data_ord md_comp_ord

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

md_data_ord: $(OBJ_DORD)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

md_comp_ord: $(OBJ_CORD)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

.PHONY: clean all

clean:
	rm -f *.o *~ core md_data_ord md_comp_ord

