EXE = a.out

pwd = ./
debug 		= $(pwd)debug.o
module_mpi	= $(pwd)module_mpi_topology.o $(pwd)module_global.o $(pwd)module_mpi_subdomain.o
solver		= $(pwd)module_mpi_thermal.o $(pwd)module_mpi_momentum.o $(pwd)module_mpi_pressure.o 
analysis    = $(pwd)module_mpi_Post.o 
main 		= $(pwd)main.o 

OBJS   	= $(debug) $(module_mpi) $(solver) $(analysis) $(main)

# FLAG 		= -O0 -g -fPIC -fcheck=all --warn-all
# FLAG 		= -O0 -g -fPIC -check all -warn all -WB #ifort
# FLAG 		= -O0 -g -fPIC -check all -WB #ifort
FLAG 		= -O3 -fPIC -xMIC-AVX512#ifort

PaScaLINC = -I/${HOME}/ptdma/PaScaL_TDMA/include
PaScaLLIB = -L/${HOME}/ptdma/PaScaL_TDMA/lib -lpascal_tdma

FFTWLIB = -L/${HOME}/fft/lib		# For CSE XeonPhi
FFTWINC = -I/${HOME}/fft/include	# For CSE XeonPhi
# FFTWLIB = -L/usr/local/lib		# For CSE XeonPhi
# FFTWINC = -I/usr/local/include	# For CSE XeonPhi

FFTWFLAGS 	=  -lfftw3 

.SUFFIXES: .f90 .o
.f90.o:
	mpiifort -c ${FLAG} $(FFTWINC) $(FFTWLIB) $*.f90 ${PaScaLINC} $(FFTWFLAGS)

all:${OBJS}
	mpiifort ${FLAG} $(FFTWINC) $(FFTWLIB) ${OBJS} -o ${EXE} ${PaScaLINC} ${PaScaLLIB} $(FFTWFLAGS)

rm:
	clear;rm *.mod;rm *.o;rm *.out;rm *.plt *.PLT; rm fort.*;rm *.txt;rm NOB.e*