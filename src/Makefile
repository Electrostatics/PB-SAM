## COMPILE FLAGS ###

CC = g++
         
CFLAGS = -O3 -ffast-math -funroll-loops -fopenmp -m64 -I$(MKL_HOME)/include

PARAMS = -D__MKL -D__OMP -D__ACCURATE__  #-D__FEXT # -D__TOOLS -D__PBC  -D__POT

LINK =  -L$(MKL_HOME)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm

BIN     = ../bin

### OBJECT FILES ###

## FindContacts ##

FC_C =  expansion.cpp readutil.cpp findContacts.cpp
FC_OBJS      =      $(FC_C:.cpp=.o)

## PB ##
PB_C =  getmemory.cpp expansion.cpp triexpan.cpp expcenter.cpp  \
						rotcoeff.cpp transcoeff.cpp xforms.cpp lautil.cpp  \
						readutil.cpp molecule.cpp moldynamics.cpp 	\
						system.cpp selfp_makemat.cpp
PB_OBJS = $(PB_C:.cpp=.o)

## BDM ##
BDM_C = getmemory.cpp expansion.cpp triexpan.cpp expcenter.cpp  \
                		rotcoeff.cpp transcoeff.cpp xforms.cpp	\
						lautil.cpp readutil.cpp  \
                		molecule.cpp moldynamics.cpp system.cpp \
						contact.cpp BDmulti.cpp BDmulti_main.cpp
BDM_OBJS = $(BDM_C:.cpp=.o)

## BDNAM ##
BDNAM_C =       getmemory.cpp expansion.cpp triexpan.cpp  \
						expcenter.cpp rotcoeff.cpp transcoeff.cpp xforms.cpp 	\
						lautil.cpp readutil.cpp molecule.cpp moldynamics.cpp system.cpp \
                        contact.cpp BDnam.cpp BDnam_main.cpp
BDNAM_OBJS = $(BDNAM_C:.cpp=.o)

### TARGETS ###

bdm: $(BDM_OBJS) Makefile
	rm -f $(BIN)/bdm
	$(CC) $(CFLAGS) $(PARAMS) $(BDM_OBJS)  -o bdm $(LINK)
	mv bdm $(BIN)

bdnam: $(BDNAM_OBJS) Makefile
	rm -f $(BIN)/bdnam_10
	$(CC) $(CFLAGS) $(PARAMS) $(BDNAM_OBJS) $(LINK) -o bdnam_10
	mv bdnam_10 $(BIN)

findContacts:  $(FC_OBJS)	 Makefile
	rm -f ../bin/findContacts
	$(CC)  $(CFLAGS) -D__FC  $(FC_OBJS) $(LIB) $(LIBS) -o ../bin/findContacts

pbsam: $(PB_OBJS) Makefile
	rm -f $(BIN)/pbsam_15
	$(CC) $(CFLAGS) $(PARAMS) $(PB_OBJS)  -o pbsam_15 $(LINK)
	mv pbsam_15	 $(BIN)

.cpp.o:
	$(CC) $(CFLAGS) $(PARAMS) -c $<

clean:
	rm -rf *.o 


