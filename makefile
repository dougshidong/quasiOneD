CXX	= g++
#DEBUG	= -pg
DEBUG	= -g -Wall -Warray-bounds
#DEBUG	= -O3
OBJDIR  = ./obj
CPP_FILES	=	$(wildcard *.cpp)
OBJ_FILES	=	$(addprefix	$(OBJDIR)/,	$(notdir	$(CPP_FILES:.cpp=.o)))
CFLAGS	= -std=c++11 -Wall -fPIC $(DEBUG)
LFLAGS	= -Wall $(DEBUG)
#PETSC_DIR	=	~/Libraries/petsc/petsc-3.8.4
#PETSC_ARCH	=	linux-mkl-mpich-release
PETSC_LIB2	=	${PETSC_DIR}/$(PETSC_ARCH)/lib
INCLPE  = -L ${PETSC_DIR}/include/ ${PETSC_CC_INCLUDES}
#INCLAD  = -I/Users/ddong/adolc_base/include -I/Users/ddong/adolc_base/include/adolc -L/Users/ddong/adolc_base/lib64 -ladolc
INCLEI  = -I./eigen-git-mirror/ -I/eigen-git-mirror/unsupported/

include ${PETSC_DIR}/lib/petsc/conf/variables

EXEC	= p1.exe

$(EXEC): $(OBJ_FILES)
	$(CXX) $(LFLAGS) -o $@ $^ -L $(PETSC_LIB2) $(INCLPE) $(PETSC_LIB) ${INCLAD} ${INCLEI}

obj/%.o: %.cpp
	mkdir -p $(OBJDIR)
	$(CXX) $(CFLAGS) -c -o $@ $< -L $(PETSC_LIB2) $(INCLPE) ${INCLAD} ${INCLEI}

clean:
	@rm -f $(EXEC) $(wildcard *.o)
	@rm -rf $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)
