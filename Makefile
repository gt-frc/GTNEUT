# Multi-platform Makefile for the GTNEUT neutral transport code. 
# This Makefile has been written assuming use of the GNU Make. It may
# or may not work with other implementations of the Make utility.
# Written by John Mandrekas, GIT, 08/15/96, to replace the old Makefile.
# Currently supports CRAY, SUN and HP systems.

# June 10, 2003, Sparse matrix version (UMFPACK)

# Set the desired system:

# HP   : General Atomics HP
# SUN  : SUN SparcStation 2 running SUNOS
# SOL  : SUN ULTRA 10 running Solaris 7
# CRAY : NERSC CRAY

#SYS = PC
SYS  = PORTLAND
ifeq ($(SYS), SOL)
	F = .f
	O = .o
	E = 
	FF = f95
	LD = f95
	FFLAGS = -c -f77 -ftrap=%none -dalign -native -xarch=v8plusa -O4
	LDFLAGS = -f77 -ftrap=%none -dalign -native -xarch=v8plusa -O4
	LIBS = -L/usr/local/lib -lumfpack -xlic_lib=sunperf
endif
ifeq ($(SYS), PC)
	F = .f
	O = .o
	E = 
	FF = g77
	LD = g77
	FFLAGS = -O -c
	LDFLAGS = -O 
	LIBS = -L/usr/local/lib -llapack_sun4 -lblas_sun4
endif
ifeq ($(SYS), PORTLAND)
	F = .f
	O = .o
	E = 
	FF = pgf90
	LD = pgf90
	FFLAGS = -O -c 
	LDFLAGS = -O 
	LIBS = -L/usr/local/lib -L/c/source/PGI/pgi/linux86/7.2-2/lib /c/source/umfpack-2.2.1/lib/libumfpack.a /usr/lib/gcc/i386-redhat-linux/3.4.3/libg2c.a -llapack -lblas 
endif
#/usr/lib/libf2c.a -llapack -lblas 
ifeq ($(SYS), HP)
        F = .f
        O = .o
        E =
        FF = fort77
        LD = fort77
        FFLAGS = -c +O4 +DA1.1 +DS1.1 -K +U77
        LDFLAGS = -L/d/hp/lib -L/opt/fortran/lib
        LIBS = -llapack -lblas
endif                                                
ifeq ($(SYS), CRAY)
	F = .f
	O = .o
	E = 
	FF = f90
	LD = segldr
	FFLAGS = -ev -c
	LDFLAGS = 
	LIBS =
endif

SOURCES= main.f		\
         calctransm.f	\
	 transmcoeff.f	\
	 rectinp.f      \
	 checkinp.f 	\
	 calcrparms.f	\
	 calcRect.f 	\
	 tij.f		\
	 qgauss.f	\
	 calcmfp.f	\
	 svjanev.f	\
	 degasread.f    \
	 svdegas.f	\
	 calcrefln.f	\
	 reflect.f	\
	 escape.f	\
	 setup.f	\
	 solvers.f	\
	 output.f	\
	 calcxswms.f	\
	 wmsdata.f	\
	 simpson.f	\
	 bickley.f	\
	 ndata.f	\
	 zstop.f        \
	 fem.f          \
	 pbalance.f

OBJ = $(SOURCES:$F=$O)

xneut$E : $(OBJ)
	echo 'Makefile for GTNEUT     20081010 tbt'
	echo 'SYS = ' $(SYS)
	echo ' '
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

%$O : %$F
	$(FF) $(FFLAGS) $<

# Include file dependencies:

main.o : \
	 consts.inc \
	 neutGlob.inc \
	 comiou.inc 

rectinp.o :	\
	neutGlob.inc \
	comiou.inc

calcRect.o : 	\
	 locGeom.inc \
         neutGlob.inc

calcmfp.o : \
	 consts.inc	\
	 locGeom.inc	\
	 neutGlob.inc

calcrparms.o : \
	 consts.inc	\
	 locGeom.inc	\
	 neutGlob.inc

calctransm.o : \
	 neutGlob.inc

transmcoeff.o : \
	 consts.inc	\
	 locGeom.inc	\
	 neutGlob.inc

calcxswms.o : \
	 wmsdata.inc

checkinp.o : \
	 consts.inc \
	 neutGlob.inc \
	 comiou.inc

degasread.o : degdat.inc \
		comiou.inc

svdegas.o   : degdat.inc

calcrefln.o : neutGlob.inc

escape.o : \
	 consts.inc \
	 locGeom.inc \
	 neutGlob.inc \

output.o : \
	 consts.inc \
	 neutGlob.inc \
	 comiou.inc 

setup.o : \
	 consts.inc \
	 neutGlob.inc \


solvers.o : \
	 consts.inc \
	 neutGlob.inc

tij.o : \
	 locGeom.inc \

wmsdata.o : \
	 wmsdata.inc \

ndata.o  :	\
	   comiou.inc \
	   consts.inc

pbalance.o: \
	   neutGlob.inc

fem.o: \
	 neutGlob.inc \
	 esc.inc

clean:
	rm -f xneut$E $(OBJ)

