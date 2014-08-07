INCL	=   -I/usr/lib/gcc/x86_64-redhat-linux/4.4.6/include\
            -I/share/apps/intel/composer_xe_2013.2.146/mkl/include/
LIB	=   -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 \
            /share/apps/intel/composer_xe_2013.2.146/mkl/lib/intel64/libmkl_intel_lp64.so \
            /share/apps/intel/composer_xe_2013.2.146/mkl/lib/intel64/libmkl_gnu_thread.so \
            /share/apps/intel/composer_xe_2013.2.146/mkl/lib/intel64/libmkl_core.so 
#            /share/apps/intel/composer_xe_2013.2.146/compiler/lib/intel64/libiomp5.so 
         
#CFLAGS = -O3 -ffast-math -fopenmp -funroll-loops -fstrict-aliasing -D__LOWMEM -D__MKL  -D__OMP -D__ACCURATE__    #pbsam
CFLAGS = -O3 -ffast-math -fopenmp  -funroll-loops -fstrict-aliasing -D__LOWMEM -D__MKL -D__OMP  -D__ACCURATE__  #BD & etc

FLAGSTAT = -Wl --start-group
 
#CC      = /share/apps/gcc/bin/g++47
 CC      = /usr/bin/g++ 
LINK	= 

#####################################33
LIBS =  -lstdc++ -lm

BIN	= ../bin

PQR_C = readutil.cpp expansion.cpp cfg2pqr.cpp

PB_C =   getmemory.cpp pdb.cpp protein.cpp expansion.cpp triexpan.cpp expcenter.cpp rotcoeff.cpp transcoeff.cpp\
	 xforms.cpp lautil.cpp readutil.cpp molecule.cpp moldynamics.cpp  system.cpp\
	mytest.cpp 

BDM_C =   getmemory.cpp pdb.cpp protein.cpp expansion.cpp triexpan.cpp expcenter.cpp rotcoeff.cpp transcoeff.cpp\
	 xforms.cpp lautil.cpp readutil.cpp molecule.cpp moldynamics.cpp system.cpp\
	contact.cpp BDmulti.cpp BDmulti_main.cpp


COMPPORE_C =   getmemory.cpp pdb.cpp protein.cpp expansion.cpp triexpan.cpp expcenter.cpp rotcoeff.cpp transcoeff.cpp\
	 xforms.cpp lautil.cpp readutil.cpp molecule.cpp moldynamics.cpp system.cpp\
	contact.cpp pore.cpp computeporepotential.cpp
BDNAM_C =   getmemory.cpp pdb.cpp protein.cpp expansion.cpp triexpan.cpp expcenter.cpp rotcoeff.cpp transcoeff.cpp\
	 xforms.cpp lautil.cpp readutil.cpp molecule.cpp moldynamics.cpp system.cpp\
	contact.cpp BDnam.cpp BDnam_main.cpp

PQR_OBJS      =       $(PQR_C:.cpp=.o) 

PB_OBJS      =       $(PB_C:.cpp=.o) 

BDNAM_OBJS      =       $(BDNAM_C:.cpp=.o)

BDM_OBJS      =       $(BDM_C:.cpp=.o)
COMPPORE_OBJS =  $(COMPPORE_C:.cpp=.o)

bdm:  $(BDM_OBJS) Makefile
	rm -f $(BIN)/bdm
	$(CC)  $(CFLAGS) $(BDM_OBJS) $(LIB) $(LIBS) -o bdm $(LINK)
	mv bdm $(BIN)

bdnam:  $(BDNAM_OBJS) Makefile
	rm -f $(BIN)/bdnam
	$(CC)  $(CFLAGS) $(BDNAM_OBJS) $(LIB) $(LIBS) -o bdnam $(LINK)
	mv bdnam $(BIN)

comppore:  $(COMPPORE_OBJS) Makefile
	rm -f $(BIN)/comppore
	$(CC)  $(CFLAGS) $(COMPPORE_OBJS) $(LIB) $(LIBS) -o comppore $(LINK)
	mv comppore $(BIN)

pbsam:  $(PB_OBJS) Makefile
	rm -f $(BIN)/pbsam
	$(CC)    $(CFLAGS)  $(PB_OBJS) $(LIB) $(LIBS)   -o pbsam $(LINK)
	mv pbsam $(BIN)

cfg2pqr:$(PQR_OBJS) Makefile
	rm -f cfg2pqr
	$(CC)  $(CFLAGS) $(PQR_OBJS) $(LIB) $(LIBS) -o cfg2pqr $(LINK)
	mv cfg2pqr $(BIN)

.cpp.o:
	$(CC) $(CFLAGS)   $(INCL)   $(LIBS) -c $<


clean:
	rm *.o 

# DO NOT DELETE THIS LINE -- make depend depends on it.

