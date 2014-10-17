FC=ifort
FFLAGS="-CB"
EXE="T_WHAM.x"
${FC} ${FFLAGS} -c fp_kind.f90
${FC} ${FFLAGS} -c constant.f90
${FC} ${FFLAGS} -c T_WHAM.f90
${FC} ${FFLAGS} -c T_WHAM_caller.f90
${FC} ${FFLAGS} -c lib.f90
${FC} ${FFLAGS} fp_kind.o constant.o T_WHAM.o T_WHAM_caller.o lib.o -o ${EXE}
rm *.o *.mod
