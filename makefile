CXX	= g++
INCLEI  = ./eigen/
#DEBUG	= -g -Wall
DEBUG	= -O3
OBJDIR  = ./obj
CPP_FILES	=	$(wildcard *.cpp)
OBJ_FILES	=	$(addprefix	$(OBJDIR)/,	$(notdir	$(CPP_FILES:.cpp=.o)))
CFLAGS	= -Wall $(DEBUG) -I $(INCLEI)
LFLAGS  = -Wall $(DEBUG)

EXEC	= p1.exe

$(EXEC): $(OBJ_FILES)
	g++ $(LFLAGS) -o $@ $^

obj/%.o: %.cpp
	g++ $(CFLAGS) -c -o $@ $<

clean:
	@rm -f $(EXEC) $(wildcard *.o)
	@rm -rf $(OBJDIR)
