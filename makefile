CXX	= g++
INCLGN  = ./gnuplot-iostream
INCLEI  = ./eigen/
DEBUG	= -g
CPPOBJ	= main.o quasiOneD.o grid.o optimizer.o flux.o
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

quasiOneD.o : quasiOneD.h quasiOneD.cpp flux.h
	$(CC) $(CFLAGS) quasiOneD.cpp 

optimizer.o : optimizer.h optimizer.cpp quasiOneD.h grid.h adjoint.h
	$(CC) $(CFLAGS) optimizer.cpp

flux.o : flux.h flux.cpp
	$(CC) $(CFLAGS) flux.cpp

adjoint.o : adjoint.h adjoint.cpp

clean:
	\rm *.o p1.*
