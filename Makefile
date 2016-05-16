#
# Makefile for program Adaptive Fast Multipole Poisson Boltzmann (AFMPB) solver
#

C = icc 
FORT = ifort 
CFLAG = -std=c99
FFLAG = 
LFLAG = -nofor_main -lcilkrts

SRC = ./src
OBJ = ./build

SRCC = $(wildcard $(SRC)/*.c)
SRCF = $(wildcard $(SRC)/*.f)

OBJC = $(patsubst $(SRC)/%.c, $(OBJ)/%.o, $(SRCC))
OBJF = $(patsubst $(SRC)/%.f, $(OBJ)/%.o, $(SRCF))

EXEC = afmpb

all: $(EXEC)

$(EXEC): $(OBJC) $(OBJF)
	$(FORT) $(LFLAG) -o $(EXEC) $(OBJC) $(OBJF)

$(OBJ)/%.o: $(SRC)/%.c
	$(C) $(CFLAG) -c -I./include -o $@ $<

$(OBJ)/%.o: $(SRC)/%.f
	$(FORT) $(FFLAG) -c -o $@ $<

clean:
	rm -rf $(OBJ)/*.o
	rm -rf $(EXEC)
	rm -rf $(SRC)/*~ ./include/*~
	rm -rf *~

