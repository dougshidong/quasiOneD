CXX	= g++
INCL	= ./eigen
DEBUG	= -g
CPPOBJ	= main.o quasiOneD.o
CPPH	= quasiOneD.h
CFLAGS	= -Wall -c $(DEBUG)
LFLAGS  = -Wall $(DEBUG)

EXEC	= p1.exe

$(EXEC) : $(CPPOBJ)
	$(CXX) $(LFLAGS) $(CPPOBJ) -o $(EXEC)

main.o: main.cpp quasiOneD.h
	$(CC) $(CFLAGS) main.cpp

quasiOneD.o : quasiOneD.h quasiOneD.cpp
	$(CC) $(CFLAGS) -I $(INCL) quasiOneD.cpp

clean:
	\rm *.o p1.*
