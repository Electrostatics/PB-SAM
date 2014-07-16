## COMPILE FLAGS ###

CC = g++
         
CFLAGS = -O3 -ffast-math -funroll-loops -fopenmp -m64 -I$(MKL_HOME)/include

PARAMS = -D__MKL -D__OMP -D__ACCURATE__  #-D__FEXT # -D__TOOLS -D__PBC  -D__POT

LINK =  -L$(MKL_HOME)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm

BIN     = ../bin

### OBJECT FILES ###

## FindContacts ##

FC_C =  contact.cpp findContacts.cpp
FC_OBJS      =      $(FC_C:.cpp=.o)

## Makesphere ##

MS_C =  spheres.cpp spheres_contact.cpp readutil.cpp \
					 makesphere.cpp expansion.cpp
MS_OBJS      =      $(MS_C:.cpp=.o)

## PQR ##

PQR_C = readutil.cpp expansion.cpp cfg2pqr.cpp
PQR_OBJS = $(PQR_C:.cpp=.o)

## PB ##
PB_C =  getmemory.cpp expansion.cpp triexpan.cpp expcenter.cpp  \
						rotcoeff.cpp transcoeff.cpp xforms.cpp lautil.cpp  \
						readutil.cpp molecule.cpp moldynamics.cpp 	\
						system.cpp selfp_makemat.cpp
PB_OBJS = $(PB_C:.cpp=.o)

## BDM ##
BDM_C = getmemory.cpp protein.cpp expansion.cpp triexpan.cpp expcenter.cpp  \
                		rotcoeff.cpp transcoeff.cpp xforms.cpp	\
						lautil.cpp readutil.cpp  \
                		molecule.cpp moldynamics.cpp system.cpp \
						contact.cpp BDmulti.cpp BDmulti_main.cpp
BDM_OBJS = $(BDM_C:.cpp=.o)

## NAFION ##
NAFION_C =    getmemory.cpp expansion.cpp triexpan.cpp 	\
						pemfc_tools.cpp expcenter.cpp rotcoeff.cpp \
                        transcoeff.cpp xforms.cpp lautil.cpp readutil.cpp 	\
						molecule.cpp moldynamics.cpp system.cpp \
						contact.cpp BDnafion.cpp BDnafion_main.cpp
NAFION_OBJS = $(NAFION_C:.cpp=.o)

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
	rm -f $(BIN)/bdnam
	$(CC) $(CFLAGS) $(PARAMS) $(BDNAM_OBJS) $(LINK) -o bdnam
	mv bdnam $(BIN)

findContacts:  $(FC_OBJS)	 Makefile
	rm -f ../bin/findContacts
	$(CC)  $(CFLAGS) findContacts.cpp $(LIB) $(LIBS) -o ../bin/findContacts

makesphere:  $(MS_OBJS) Makefile
		rm -f ../bin/makesphere
		$(CC)  $(CFLAGS) $(MS_OBJS) $(LIB) $(LIBS) -o ../bin/makesphere

nafion_nopbc: $(NAFION_OBJS) Makefile
	rm -f $(BIN)/nafion_nopbc
	$(CC) $(CFLAGS) $(PARAMS) $(NAFION_OBJS) $(LINK) -o nafion_nopbc
	mv nafion_nopbc $(BIN)

nafion: $(NAFION_OBJS) Makefile
	rm -f $(BIN)/nafion
	$(CC) $(CFLAGS) $(PARAMS) $(NAFION_OBJS) $(LINK) -o nafion
	mv nafion $(BIN)

pbsam: $(PB_OBJS) Makefile
	rm -f $(BIN)/pbsam_30
	$(CC) $(CFLAGS) $(PARAMS) $(PB_OBJS)  -o pbsam_30 $(LINK)
	mv pbsam_30	 $(BIN)

cfg2pqr: $(PQR_OBJS) Makefile
	rm -f cfg2pqr
	$(CC) $(CFLAGS) $(PARAMS) $(PQR_OBJS) $(LINK) -o cfg2pqr
	mv cfg2pqr $(BIN)

.cpp.o:
	$(CC) $(CFLAGS) $(PARAMS) -c $<

clean:
	rm -rf *.o 


