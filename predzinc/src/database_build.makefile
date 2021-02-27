# Makefile for database_build.cpp
CC         = g++
CFLAGS     = -Wall -O3  
PREDZINC   = ../
LIBS       = -lm -lmyfunc
LIB_PATH   = $(PREDZINC)/lib
INCLUDE    = $(PREDZINC)/src
SRC        = database_build.cpp
OBJ        = $(SRC:.cpp=.o) 
EXE        = database_build
RM         = /bin/rm -f
CP         = /bin/cp -f
$(EXE): $(OBJ)
	$(CC) $(CFLAGS)  $(OBJ) -L$(LIB_PATH) -o $(EXE)  $(LIBS)
#compile and assemble C++/C source files into object files
# -c flag tells the compiler to create only OBJ files
$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c -I$(INCLUDE) $(SRC) 
install:
	$(CP)  $(EXE)  $(PREDZINC)/bin/
clean:
	$(RM)  $(OBJ) $(EXE)
