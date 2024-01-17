include ../../example.mk
include ../../common.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

ps_cma_es: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: ps_cma_es

run: all
	mpirun --oversubscribe  -np 2 ./ps_cma_es

.PHONY: clean all run

clean:
	rm -f *.o *~ core ps_cma_es

