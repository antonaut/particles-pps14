#
# For computer with libs installed running Ubuntu Linux 12.04 or higher...
#
CC = g++
MPCC =  mpic++
OPENMP = -fopenmp
LIBS = -lm
CFLAGS = -g -Wall

TARGETS = serial pthreads openmp mpi linearserial linpthreads linopenmp linmpi

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o
pthreads: pthreads.o common.o
	$(CC) -o $@ $(LIBS) pthreads.o common.o -pthread
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o
linearserial: linearserial.o common.o
	$(CC) -o $@ $(LIBS) linearserial.o common.o
linopenmp: linopenmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) linopenmp.o common.o
linpthreads: linpthreads.o common.o
	$(CC) -o $@ $(LIBS) linpthreads.o common.o -pthread
linmpi: linmpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) linmpi.o common.o


openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
pthreads.o: pthreads.cpp common.h
	$(CC) -c $(CFLAGS) pthreads.cpp -pthread
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp
linearserial.o: linearserial.cpp common.h
	$(CC) -c $(CFLAGS) linearserial.cpp
linopenmp.o: linopenmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) linopenmp.cpp
linpthreads.o:
	$(CC) -c $(CLAGS) linpthreads.cpp -pthread
linmpi.o: linmpi.cpp common.h
	$(MPCC) -c $(CFLAGS) linmpi.cpp

clean:
	rm -f *.o $(TARGETS)
