include ../../../common.mk
include ../../../example.mk

OBJ = mainDCPSE.o
OBJ2 = mainFD.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

DcpseLid: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
FDLid: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: DcpseLid FDLid

run: all
	mpirun -np 1 ./DcpseLid

.PHONY: clean all run

clean:
	rm -f *.o *~ core DcpseLid FDLid
