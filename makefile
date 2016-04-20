CXX	= g++
INCLEI  = ./eigen/
INCLPE  = -L ./petsc-3.6.3/include/ ${PETSC_CC_INCLUDES}
#DEBUG	= -g -Wall
DEBUG	= -O3
OBJDIR  = ./obj
CPP_FILES	=	$(wildcard *.cpp)
OBJ_FILES	=	$(addprefix	$(OBJDIR)/,	$(notdir	$(CPP_FILES:.cpp=.o)))
CFLAGS	= -Wall $(DEBUG) -I $(INCLEI)
LFLAGS	= -Wall $(DEBUG)
PETSC_DIR	=	./petsc-3.6.3
PETSC_ARCH	=	linux-cpp
PETSC_LIB2	=	${PETSC_DIR}/$(PETSC_ARCH)/lib

include ${PETSC_DIR}/lib/petsc/conf/variables

EXEC	= p1.exe

$(EXEC): $(OBJ_FILES)
	$(CXX) $(LFLAGS) -o $@ $^ -L $(PETSC_LIB2) $(INCLPE) $(PETSC_LIB)

obj/%.o: %.cpp
	mkdir -p $(OBJDIR)
	$(CXX) $(CFLAGS) -c -o $@ $< -L $(PETSC_LIB2) $(INCLPE)

clean:
	@rm -f $(EXEC) $(wildcard *.o)
	@rm -rf $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)
