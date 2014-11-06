.SUFFIXES: 
.SUFFIXES: .f90 .o

FC = ifort
FFLAGS = -CB
INCLUDE = 
LIBS = 

EXE = WHAM.x

OBJS = precision_m.o constant.o lib.o T_WHAM.o T_WHAM_caller.o

all:	${EXE}


$(EXE):$(OBJS)
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(LIBS)

%.o %mod: %.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $<

clean:
	/bin/rm -f $(EXE) $(OBJS) *.mod
