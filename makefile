CXX	= mpicxx
#CXX	= gcc-8
#DEBUG	= -pg
#DEBUG	= -g -Wall -Warray-bounds -march=native
DEBUG	= -O3
OBJDIR  = ./obj
CPP_FILES	=	$(wildcard *.cpp)
OBJ_FILES	=	$(addprefix	$(OBJDIR)/,	$(notdir	$(CPP_FILES:.cpp=.o)))
CFLAGS	= -std=c++17 -Wall -fPIC $(DEBUG)
LFLAGS	= -Wall $(DEBUG)
#PETSC_DIR	=	~/Libraries/petsc/petsc-3.8.4
#PETSC_ARCH	=	linux-mkl-mpich-release
LINKPE	=	-L ${PETSC_DIR}/$(PETSC_ARCH)/lib
LINKAD  = -L/${ADOLC_DIR}/lib64 -ladolc
INCLPE  = -I ${PETSC_DIR}/include/ ${PETSC_CC_INCLUDES}
INCLAD  = -I${ADOLC_DIR}/include -I${ADOLC_DIR}/include/adolc
INCLEI  = -I./eigen-git-mirror/
#INCLEI  = -I./eigen-git-mirror/ -I/eigen-git-mirror/unsupported/

include ${PETSC_DIR}/lib/petsc/conf/variables

EXEC	= p1.exe

$(EXEC): $(OBJ_FILES)
	$(CXX) $(CFLAGS) $(LFLAGS) -o $@ $^ $(LINKPE) ${PETSC_LIB} ${LINKAD}

obj/%.o: %.cpp
	mkdir -p $(OBJDIR)
	$(CXX) $(CFLAGS) -c -o $@ $< ${INCLAD} ${INCLEI} ${INCLPE}

clean:
	@rm -f $(EXEC) $(wildcard *.o)
	@rm -rf $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)
