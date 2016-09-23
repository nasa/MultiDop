# Where to install things
PREFIX=/Users/tjlang/Documents/SvrWx/NASA_STORM/MultiDop/MultiDop/src

# Usual places to put things under $PREFIX
BINDIR=${PREFIX}/bin
MANDIR=${PREFIX}/share/man

RM = rm -fr
CP = cp -p -f

CC = gcc -std=c99
INCLUDES = -I/opt/local/include
CFLAGS_GENERIC = -Wall -Wmissing-prototypes ${INCLUDES}
CFLAGS_DEBUG = -g ${CFLAGS_GENERIC}
CFLAGS_OPT = -O3 ${CFLAGS_GENERIC}
CFLAGS = ${CFLAGS_OPT}

FC = gfortran
FFLAGS_GENERIC = -Wall
FFLAGS_DEBUG = -g ${FFLAGS_GENERIC}
FFLAGS_OPT = -O3 ${FFLAGS_GENERIC}
FFLAGS = ${FFLAGS_OPT}

LIB_DIRS = -L/opt/local/lib
#EFENCE_LIBS = -lefence
LIBS = ${EFENCE_LIBS} -lnetcdf -lm

EXECS = DDA

all : ${EXECS}

install : ${EXECS}
	mkdir -p ${BINDIR}
	${CP} ${EXECS} ${BINDIR}
	mkdir -p ${MANDIR}/man3
	${CP} ../man/man3/*.3 ${MANDIR}/man3

.SUFFIXES : .o

.c.o :
	${CC} ${CFLAGS} -c $<

.f.o :
	${FC} ${FFLAGS} -c $<

DDA_OBJ = DDA.o geog_proj.o geog_lib.o leise_filt.o coverage.o nnetcdf.o alloc.o

DDA : ${DDA_OBJ}
	${FC} -o DDA ${DDA_OBJ} ${LIB_DIRS} ${LIBS}

debug :
	make clean
	make CFLAGS="${CFLAGS_DEBUG}" FFLAGS="${FFLAGS_DEBUG}" ${DDA_OBJ}
	${CC} -o DDA ${DDA_OBJ} ${LIB_DIRS} ${LIBS}

DDA.o : DDA.c DDA.h

leise_filt.o : leise_filt.f

coverage.o : coverage.c DDA.h coverage.h

geog_proj.o : geog_proj.c geog_proj.h

geog_lib.o : geog_lib.c geog_lib.h

nnetcdf.o : nnetcdf.c nnetcdf.h

alloc.o : alloc.c alloc.h

clean :
	${RM} *.o ${EXECS} *.dSYM leise_example
