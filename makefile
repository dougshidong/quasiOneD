CXX	= g++
INCLEI	= ./eigen 
INCLGN  = ./gnuplot-iostream
DEBUG	= -g
CPPOBJ	= main.o quasiOneD.o grid.o
CPPH	= quasiOneD.h
CFLAGS	= -Wall -c $(DEBUG)
LFLAGS  = -Wall $(DEBUG)

EXEC	= p1.exe

$(EXEC) : $(CPPOBJ)
	$(CXX) $(LFLAGS) $(CPPOBJ) -o $(EXEC)

main.o: main.cpp quasiOneD.h grid.h
	$(CC) $(CFLAGS) main.cpp

grid.o: grid.h grid.cpp
	$(CC) $(CFLAGS) grid.cpp

quasiOneD.o : quasiOneD.h quasiOneD.cpp
	$(CC) $(CFLAGS) -I $(INCLEI) quasiOneD.cpp 

clean:
	\rm *.o p1.*
