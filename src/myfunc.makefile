# Makefile for creating the library for predzinc program, including mypro.cpp,
# mypro.cpp and subsetc.pp 
CC         = g++
CFLAGS     = -c -Wall -O3
PREDZINC   = ../
LDFLAGS    =
MODE       = shared
AR         = ar
ARFLAGS    = rc
LIBS       = -lc
LIB_PATH   = $(PREDZINC)/lib
INCLUDE    = $(PREDZINC)/src
SRC        = myfunc.cpp mypro.cpp subset.cpp
OBJ        = $(SRC:.cpp=.o)
SONAME0    = libmyfunc.so
SONAME     = libmyfunc.so.1
SOLIBNAME  = libmyfunc.so.1.0.1
ALIBNAME   = libmyfunc.so.a
CP         = /bin/cp -f
LN         = /bin/ln -sf
RM         = /bin/rm -f

$(SOLIBNAME): $(OBJ)
	$(CC) -shared -Wl,-soname,$(SONAME) -o $@ $(OBJ) $(LIBS)

.cpp.o:
	$(CC) -fPIC $(CFLAGS) -I$(INCLUDE) $< -o $@

install:
	$(CP) $(SOLIBNAME) $(LIB_PATH)
	cd ${LIB_PATH} ; $(LN) $(SOLIBNAME) ${SONAME}; $(LN) $(SONAME) $(SONAME0)

clean:
	$(RM)  $(OBJ)  $(SOLIBNAME)
