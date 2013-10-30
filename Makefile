EXE = DMRG
OBJ = ./build
FC = gfortran
TOUCH = touch
RM = rm
LIB = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/usr/lib/
CFLAGS = -fdefault-real-8 -fimplicit-none -ffree-line-length-none
LFLAGS = -framework vecLib 

all: $(EXE)
$(EXE): $(OBJ)/TENSOR.o $(OBJ)/MATHIO.o $(OBJ)/$(EXE).o
	@$(FC) -o $(EXE) $(OBJ)/TENSOR.o $(OBJ)/MATHIO.o $(OBJ)/$(EXE).o libarpack.a  -L $(LIB) $(LFLAGS) -I $(OBJ)
$(OBJ)/TENSOR.o: TENSOR.f90
	@$(FC) -c TENSOR.f90 -o $(OBJ)/TENSOR.o -L $(LIB) $(CFLAGS) -J $(OBJ)
$(OBJ)/MATHIO.o: MATHIO.f90
	@$(FC) -c MATHIO.f90 -o $(OBJ)/MATHIO.o -L $(LIB) $(CFLAGS) -J $(OBJ)
$(OBJ)/$(EXE).o: $(EXE).f90
	@$(FC) -c $(EXE).f90 -o $(OBJ)/$(EXE).o -L $(LIB) $(CFLAGS) -J $(OBJ)
clean: 
	$(RM) build/*
	$(RM) $(EXE)
