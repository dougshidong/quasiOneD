CXX	= g++
INCLGN  = ./gnuplot-iostream
INCLEI  = ./eigen/
DEBUG	= -g -Wall
#DEBUG	= -O3
CPPOBJ	= main.o quasiOneD.o grid.o optimizer.o flux.o adjoint.o timestep.o globals.o
CPPH	= quasiOneD.h
CFLAGS	= -Wall -c $(DEBUG) -I $(INCLEI)
LFLAGS  = -Wall $(DEBUG)

EXEC	= p1.exe

$(EXEC) : $(CPPOBJ)
	$(CXX) $(LFLAGS) $(CPPOBJ) -o $(EXEC)

main.o: main.cpp quasiOneD.h grid.h optimizer.h
	$(CC) $(CFLAGS) main.cpp

grid.o: grid.h grid.cpp
	$(CC) $(CFLAGS) grid.cpp

quasiOneD.o : quasiOneD.h quasiOneD.cpp flux.h timestep.h
	$(CC) $(CFLAGS) quasiOneD.cpp 

optimizer.o : optimizer.h optimizer.cpp quasiOneD.h grid.h adjoint.h
	$(CC) $(CFLAGS) optimizer.cpp

flux.o : flux.h flux.cpp
	$(CC) $(CFLAGS) flux.cpp

timestep.o : timestep.h timestep.cpp flux.h globals.h
	$(CC) $(CFLAGS) timestep.cpp

adjoint.o : adjoint.h adjoint.cpp
	$(CC) $(CFLAGS) adjoint.cpp

globals.o : globals.h globals.cpp
	$(CC) $(CFLAGS) globals.cpp

clean:
	\rm *.o p1.*
