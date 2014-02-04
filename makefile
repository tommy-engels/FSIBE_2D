## Makefile for 2D DNS
## Build on Babel and on Duke

MOD_FILES = share_vars gif_util nrtype nrutil Interpolation mouvement CompressMatrixCSR FieldExport PerformanceMeasurement BeamForces 
SUB_FILES = cal_vis cofdx cofdxdx cofdxdy NavierStokes\
cofdy cofdydy dealiase_mask mean_velocity\
init_fields params poisson integrate_position BeamIO\
save_fields time_step init_beam create_mask
# mkl_lapack
PROG_FILE = dns


#----------------------------------------------------------------------------------------------------------------------
ifeq ($(CONF),duke)
	FPAR = -openmp -parallel -lpthread
 	FOPTS = -O3 -xW -fpp -r8 -static 
#  	FOPTS = -O3 -xW -fpp -r8 -static -fpe0 -warn -traceback -debug extended
	# -fpe0 -warn -traceback -debug extended
        #-CB -debug -traceback
	# -ftrapuv

	vpath %.f90 /Softs/intel/mkl/10.0.011/include

	MKL_BASE = /Softs/intel/mkl/10.1.0.015
	MKL_PATH = $(MKL_BASE)/lib/em64t
	
	IFACE_LIB=$(MKL_PATH)/libmkl_solver_lp64.a $(MKL_PATH)/libmkl_intel_lp64.a
	CORE_LIB=$(MKL_PATH)/libmkl_core.a


	THREADING_LIB=$(MKL_PATH)/libmkl_intel_thread.a
	LGUIDE = $(MKL_PATH)/libguide.a 



	MKL_LIBS=$(IFACE_LIB) -Wl,--start-group $(THREADING_LIB) $(CORE_LIB) -Wl,--end-group $(LGUIDE)

	FFT_LINK = -L$(MKL_PATH) $(MKL_LIBS)

	COF_FILE = cof_mkl 
	SUPPORT_FILE = mkl_dfti mkl_lapack mkl_pardiso

	FF = ifort
endif
#----------------------------------------------------------------------------------------------------------------------
ifeq ($(CONF),meso)
	# ----------------------------------------
	# to make on mesocentre:
	# module load intel/11.1.080
	# unset MAKEFLAGS
	# make meso
	# ----------------------------------------

	# the following helps finding mkl_dfti.f90
	vpath %.f90 /LOGINSTO/softs/intel/11.1.080/Compiler/11.1/080/mkl/include
	
	FPAR = -openmp -parallel -lpthread
 	FOPTS = -O3 -xW -fpp -r8 -static
#  	FOPTS = -O3 -xW -fpp -r8 -g -static -check all -warn -traceback -debug extended 

	MKL_BASE = /LOGINSTO/softs/intel/11.1.080/Compiler/11.1/080/mkl/
	MKL_PATH = $(MKL_BASE)/lib/em64t
	
	IFACE_LIB=$(MKL_PATH)/libmkl_solver_lp64.a $(MKL_PATH)/libmkl_intel_lp64.a
	CORE_LIB=$(MKL_PATH)/libmkl_core.a

	THREADING_LIB=$(MKL_PATH)/libmkl_intel_thread.a
	
	LGUIDE = /LOGINSTO/softs/intel/11.1.080/Compiler/11.1/080/lib/intel64/libguide.a

	MKL_LIBS=$(IFACE_LIB) -Wl,--start-group $(THREADING_LIB) $(CORE_LIB) -Wl,--end-group $(LGUIDE)

	FFT_LINK = -L$(MKL_PATH) $(MKL_LIBS)

	COF_FILE = cof_mkl 
	SUPPORT_FILE = mkl_dfti mkl_lapack mkl_pardiso SolidSolver

	FF = ifort
endif
#----------------------------------------------------------------------------------------------------------------------

COMMON = $(addsuffix .o ,$(MOD_FILES) $(SUPPORT_FILE) $(COF_FILE) $(SUB_FILES))
RES = $(addsuffix .res ,$(PROG_FILE))


#-------------------------------------------------------------------------------

help:
	@echo Usage: make {duke/clean/help}
duke:
	rm -f dns.out share_vars.o share_vars.res
	$(MAKE) build CONF=duke
meso:
	rm -f dns.out share_vars.o share_vars.res
	$(MAKE) build CONF=meso
build:
	$(MAKE) $(COMMON) $(RES)	
clean:
	rm -f dns.out *.o *.mod

#-------------------------------------------------------------------------------

$(COMMON): %.o: %.f90
	$(FF) $(FOPTS) $(FPAR) -c $< -o $@

$(RES): %.res: %.f90
	$(FF) $(FOPTS) $(COMMON) $< $(FPAR) $(FFT_LINK) -lm -o $*.out
	
	

#-------------------------------------------------------------------------------
