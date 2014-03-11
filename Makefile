#
# Computers with Red Hat Enterprise Linux 5 in the computer room 648, KTH Forum, Kista
#
CC = g++
MPCC =  mpic++
OPENMP = -fopenmp
LIBS = -lm
CFLAGS = -g

TARGETS = serial pthreads openmp mpi linearserial

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
	

clean:
	rm -f *.o $(TARGETS)
