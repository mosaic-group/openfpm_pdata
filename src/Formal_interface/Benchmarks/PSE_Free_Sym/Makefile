include paths.mk

CC=mpic++

all: PSE_Free_Sym

%.o: %.cpp %.hpp
	$(CC) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

PSE_Free_Sym: PSE_Free_Sym.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)


run: PSE_Free_Sym
	mpirun --oversubscribe -np 2 ./PSE_Free_Sym $(ARGS);

.PHONY: clean all run

clean:
	rm -f *.o *~ core *.pvtp *.vtp *.vtk

