# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

DEBUG := t
DDEBUG := $(subst T,t,$(DEBUG))



OBJ_DIR = obj
MOD_DIR = obj
SRC_DIR = src
MAIN_DIR = app

#=================================================================================
# Some variables for the external libraries
#=================================================================================
MAIN_path   := $(shell pwd)

#=================================================================================
#=================================================================================
# Compiler?
#Possible values: ifort, ifx, gfortran (default), nagfor 
 FC := gfortran
#FC := ifort
#FC := nagfor
#
# Optimize? Empty: default Optimization; 0: No Optimization; 1 Optimization
OPT := 1
## OpenMP? Empty: default with OpenMP; 0: No OpenMP; 1 with OpenMP
OMP := 1
## Lapack/blas/mkl? Empty: default with Lapack; 0: without Lapack; 1 with Lapack
LAPACK := 1
## Arpack? Empty/0: without Arpack; 1: with Arpack
ARPACK := 0
## force the default integer (without kind) during the compillation.
## default 4: , INT=8 (for kind=8)
INT := 4
## change the real kind
## default real64: , possibilities, real32, real64, real128
RKIND := real64
# For some compilers (like lfortran), real128 (quadruple precision) is not implemented
# WITHRK16 = 1 (0) compilation with (without) real128
WITHRK16 := 0
# branch of the external libraries (main, dev)
BRANCH := main
# how to clean (recursively (1) or not (0)) the external libraries (*_loc)
RECCLEAN := 1
#=================================================================================
ifeq ($(FC),)
  override FC := gfortran
endif
FFC := $(FC)

ifeq ($(OPT),)
  override OPT := 1
endif
ifneq ($(OPT),$(filter $(OPT),0 1))
  $(info *********** OPT (optimisation):        $(OPT))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible option values)
endif
OOPT := $(OPT)

ifeq ($(OMP),)
  override OMP := 1
endif
ifneq ($(OMP),$(filter $(OMP),0 1))
  $(info *********** OMP (openmp):        $(OMP))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible option values)
endif
OOMP := $(OMP)

ifeq ($(LAPACK),)
  override LAPACK := 1
endif
ifneq ($(LAPACK),$(filter $(LAPACK),0 1))
  $(info *********** LAPACK:        $(LAPACK))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible option values)
endif
LLAPACK := $(LAPACK)

ifeq ($(ARPACK),)
  override ARPACK := 0
endif
# We cannot use ARPACK without lapack
ifeq ($(LAPACK),0)
  override ARPACK := 0
endif
ifneq ($(ARPACK),$(filter $(ARPACK),0 1))
  $(info *********** ARPACK:        $(ARPACK))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible option values)
endif
AARPACK := $(ARPACK)

# ARPACK-NG is installed locally in EXT_Lib when ARPACK=1.  No Homebrew,
# MacPorts or system ARPACK path is probed.  ARPACK_LIB is a project-local
# archive, built on the current OS/architecture by scripts/install-arpack.
ARPACK_VERSION ?= 3.9.1
ARPACK_ROOT := $(MAIN_path)/EXT_Lib/arpack-ng_loc
ARPACK_PREFIX := $(ARPACK_ROOT)/install
ARPACK_LIB := $(ARPACK_PREFIX)/lib/libarpack.a
ARPACK_LIBDIR := $(ARPACK_PREFIX)/lib
ARPACK_INCLUDEDIR := $(ARPACK_PREFIX)/include
ARPACK_TARGET :=
ifeq ($(ARPACK),1)
  ARPACK_TARGET := $(ARPACK_LIB)
endif

# Optional prefixes are supplied by the user (Homebrew is deliberately not
# probed here).  This also supports an explicit installation in EXT_Lib or
# any other prefix on Linux, macOS, BSD, etc.  ARPACK paths are only used when
# ARPACK=1, so ARPACK=0 cannot accidentally search or link ARPACK.
ifeq ($(strip $(ARPACK_LIBDIR)),)
  ifneq ($(strip $(ARPACK_PREFIX)),)
    ARPACK_LIBDIR := $(ARPACK_PREFIX)/lib
  endif
endif
ifeq ($(strip $(ARPACK_INCLUDEDIR)),)
  ifneq ($(strip $(ARPACK_PREFIX)),)
    ARPACK_INCLUDEDIR := $(ARPACK_PREFIX)/include
  endif
endif
ifeq ($(strip $(LAPACK_LIBDIR)),)
  ifneq ($(strip $(BLAS_PREFIX)),)
    LAPACK_LIBDIR := $(BLAS_PREFIX)/lib
  endif
endif

ifneq ($(INT),$(filter $(INT),4 8))
  $(info *********** INT (change default integer):        $(INT))
  $(info Possible values: 4, 8)
  $(error ERROR: Incompatible option values)
endif

ifneq ($(RKIND),$(filter $(RKIND),real32 real64 real128))
  $(info *********** RKIND (select the real kind):        $(RKIND))
  $(info Possible values (case sensitive): real32 real64 real128)
  $(error ERROR: Incompatible option values)
endif

WWITHRK16 := $(WITHRK16)
ifneq ($(WITHRK16),$(filter $(WITHRK16),0 1))
  $(info *********** WITHRK16 (compilation with real128):        $(WITHRK16))
  $(info Possible values: 0, 1)
  $(error ERROR: Incompatible option values)
endif
ifeq ($(RKIND),real128)
  ifeq ($(WWITHRK16),0)
    $(info "Incompatible options:")
    $(info ***********RKIND:        $(RKIND))
    $(info ***********WITHRK16:     $(WITHRK16))
    $(error ERROR: Incompatible RKIND and WITHRK16 option values)
  endif
endif
export RKIND WITHRK16 INT  LAPACK   ARPACK   FC  OPT  OMP BRANCH
export      WWITHRK16     LLAPACK  AARPACK  FFC OOPT OOMP
export ARPACK_LIBDIR ARPACK_INCLUDEDIR LAPACK_LIBDIR ARPACK_VERSION

#=================================================================================
# extension for the library (.a), objects and modules directory
#=================================================================================
ext_obj    :=_$(FC)_opt$(OPT)_omp$(OMP)_lapack$(LAPACK)_int$(INT)_$(RKIND)
#=================================================================================
# External Libraries : Tnum-Tana FOR_EVRT QuantumModelLib EVRT_dnSVM AD_dnSVM ConstPhys nDindex QDUtilLib
#=================================================================================
EXTLIB_LIST := Tnum-Tana FOR_EVRT QuantumModelLib EVRT_dnSVM AD_dnSVM ConstPhys nDindex QDUtilLib
ifneq ($(EXTLIB_LIST),)
  ifeq ($(ExtLibDIR),)
    ExtLibDIR := $(MAIN_path)/EXT_Lib
  endif
  OK := $(shell if test -d $(ExtLibDIR); then echo "ok"; fi)
  #$(info ***********OK:       $(OK))
  ifeq ($(OK),)
    $(info WARNING: $(ExtLibDIR) does not exist!)
    $(shell mkdir $(ExtLibDIR) ) 
  endif
  export ExtLibDIR
  EXTLib_DIR  := $(addprefix $(ExtLibDIR)/, $(EXTLIB_LIST))
  $(info ***********EXTLib_DIR:       $(EXTLib_DIR))
  EXTMod      := $(addsuffix /OBJ/obj$(ext_obj), $(EXTLib_DIR))
  $(info ***********EXTOBJDIR:       $(EXTMod))
  EXTLibOBJ      := $(addsuffix /*.o, $(EXTMod))
  $(info ***********EXTLibOBJ:       $(EXTLibOBJ))

  EXTMod      := $(addprefix -I,$(EXTMod))
  EXTLib := $(shell for LLIB in $(EXTLIB_LIST) ; do echo $(ExtLibDIR)/$$LLIB"/lib"$$LLIB""$(ext_obj)".a" ; done)
endif
# We use the full TnunTana library
EXTLib := $(MAIN_path)/EXT_Lib/Tnum-Tana/libTnum-TanaFull$(ext_obj).a
# Tnum-Tana calls the QuantumModelLib C/Fortran interface.  It is a
# separate archive and must be present even when ARPACK is disabled.
QMLLib := $(MAIN_path)/EXT_Lib/QuantumModelLib/libQuantumModelLibFull$(ext_obj).a
# special files
EXT_SRC  := calc_f2_f1Q.f90 Sub_X_TO_Q_ana.f90 Calc_Tab_dnQflex.f90 read_para.f90
EXT_OBJ  := $(EXT_SRC:.f90=.o)
EXT_OBJ  += sub_system.o
ifeq ($(DDEBUG),t)
  $(info ***********EXT_SRC:         $(EXT_SRC))
  $(info ***********EXT_OBJ files:   $(EXT_OBJ))
endif
EXT_OBJ  := $(addprefix $(OBJ_DIR)/,$(EXT_OBJ))

# =====================================================================================
# External library paths
# =====================================================================================
FFAGS +=  $(EXTMod) 

LIBS := $(EXTLib) $(QMLLib)
HOST_OS := $(shell uname -s)
## Apple ld does not always rescan all members of the monolithic Tnum-Tana
## archive.  Load every member on Darwin so the QDUtil symbols embedded in
## that archive are available.  Do not link QDUtilLib separately: Tnum-Tana
## already contains it, and doing so would create duplicate symbols.
ifeq ($(HOST_OS),Darwin)
  # Do not allow an inherited FORCE_LOAD=0 to reintroduce unresolved symbols
  # from the other members of the monolithic Tnum-Tana archive.
  override FORCE_LOAD := 1
else
  FORCE_LOAD ?= 0
endif
ifeq ($(HOST_OS),Darwin)
  ifeq ($(FORCE_LOAD),1)
    LDFLAGS := -Wl,-force_load,$(EXTLib) $(QMLLib)
  else
    LDFLAGS := $(LIBS)
  endif
else
  LDFLAGS := $(LIBS)
endif

# On Darwin, -force_load extracts every member of the Tnum-Tana archive.
# The five helper sources below are also present in that archive, so keeping
# their application objects in the link creates duplicate symbols (ld: ...
# duplicate symbols). They are not needed separately when force_load is active.
ifeq ($(HOST_OS),Darwin)
  ifeq ($(FORCE_LOAD),1)
    EXT_OBJ :=
  endif
endif
ifeq ($(ARPACK),1)
  ifneq ($(strip $(ARPACK_LIBDIR)),)
    LDFLAGS += -L$(ARPACK_LIBDIR)
  endif
  LDFLAGS += $(ARPACK_LIB)
endif
ifeq ($(LAPACK),1)
  ifneq ($(strip $(LAPACK_LIBDIR)),)
    LDFLAGS += -L$(LAPACK_LIBDIR)
  endif
  # Keep LAPACK/BLAS after ARPACK: this is required for static libarpack.a.
  LDFLAGS += -llapack -lblas
endif


# =====================================================================================
FFLAGS := -cpp -Wall -Wextra -Wimplicit-interface \
          -g -Og \
          -fopenmp \
          -fcheck=all -fbacktrace \
          -finit-real=nan -finit-integer=-999999 \
          -ffree-line-length-none \
          -Warray-temporaries \
          -J$(MOD_DIR) $(EXTMod)
ifeq ($(ARPACK),1)
  FFLAGS += -D__ARPACK=1
  ifneq ($(strip $(ARPACK_INCLUDEDIR)),)
    FFLAGS += -I$(ARPACK_INCLUDEDIR)
  endif
endif
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
$(info ***********************************************************************)
$(info ***********Library name:    $(LIB_NAME))
$(info ***********Library version: $(LIB_VERSION))
$(info ***********OS:              $(OS))
$(info ***********COMPILER:        $(FC))
$(info ***********COMPILER_VER:    $(FC_VER))
$(info ***********OPTIMIZATION:    $(OPT))
$(info ***********OpenMP:          $(OMP))
$(info ***********INT:             $(INT))
$(info ***********RKIND:           $(RKIND))
$(info ***********WITHRK16:        $(WWITHRK16))
$(info ***********LAPACK:          $(LAPACK))
$(info ***********FFLAGS:          $(FFLAGS))
$(info ***********FLIB:            $(FLIB))
$(info ***********ext_obj:         $(ext_obj))
$(info )
$(info ***********ExtLibDIR:       $(ExtLibDIR))
$(info ***********EXTLib:          $(EXTLib))
$(info ***********EXTMod:          $(EXTMod))
$(info ***********LDFLAGS:         $(LDFLAGS))
$(info ***********************************************************************)
# =====================================================================================
# Main program
# =====================================================================================
MAIN = TD_SCHROD
TARGET = $(MAIN).x

# =====================================================================================
# Source and object files
# =====================================================================================

LIBSRC = InputStaging_m.f90 NUFFT_m.f90 SplitOperator_m.f90 Molec_m.f90 NDindex_m.f90 NDindex_TestSupport_m.f90 poly0rtho_m.f90 \
         SmolyakSparseGridMetadata_m.f90 SmolyakSparseBasisAlgebra_m.f90 Basis_m.f90 psi_m.f90 Checkpoint_m.f90 \
         Ana_psi_m.f90 lanczos_m.f90 Op_m.f90 param_WP0_m.f90 \
         Auto_corr_m.f90 Propa_m.f90 Hagedorn_m.f90 \
         sub_propa_m.f90 Sub_Vp_m.f90 \
         HagedornLinalg_m.f90 HagedornQuadrature_m.f90 HagedornOverlap_m.f90 \
         NDindex_hagedorn_m.f90 HagedornBasis_m.f90 \
         HagedornHamiltonian_m.f90 Initial_GWP_Basis_m.f90 TISE_m.f90 HamiltonianSpectrum_m.f90 \
         HagedornDiracFrenkel_m.f90

# Object files
OBJ0 = $(LIBSRC:.f90=.o)
OBJ  = $(addprefix $(OBJ_DIR)/,$(OBJ0))

$(info [Build info] Object files: $(OBJ))

# =====================================================================================
# Targets
# =====================================================================================

.PHONY: all lib clean distclean info

all: $(TARGET)   # ← plus de dépendance vers dependencies.mk

# Compatibility target used by the external-library build instructions.
# The executable links the rebuilt external archives and all local modules.
lib: $(TARGET)

# =====================================================================================
# Link executable
# =====================================================================================

$(TARGET): $(OBJ_DIR)/$(MAIN).o $(OBJ) $(EXT_OBJ) | $(LIBS) $(ARPACK_TARGET)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# =====================================================================================
# Compile main program
# =====================================================================================

$(OBJ_DIR)/$(MAIN).o: $(MAIN_DIR)/$(MAIN).f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# =====================================================================================
# Generic Fortran compilation rule
# =====================================================================================

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# These Tnum helper sources are part of the application link, but they live
# in the downloaded Tnum-Tana tree rather than in src/.
$(OBJ_DIR)/calc_f2_f1Q.o: $(MAIN_path)/EXT_Lib/Tnum-Tana/SRC/sub_operator_T/calc_f2_f1Q.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/Sub_X_TO_Q_ana.o: $(MAIN_path)/EXT_Lib/Tnum-Tana/SRC/sub_operator_T/Sub_X_TO_Q_ana.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/Calc_Tab_dnQflex.o: $(MAIN_path)/EXT_Lib/Tnum-Tana/SRC/sub_operator_T/Calc_Tab_dnQflex.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/read_para.o: $(MAIN_path)/EXT_Lib/Tnum-Tana/SRC/sub_pot/read_para.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/sub_system.o: $(MAIN_path)/EXT_Lib/Tnum-Tana/SRC/sub_pot/sub_system.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# ============================================================
# Dépendances manuelles entre modules (format $(OBJ_DIR)/...)
# ============================================================

$(TARGET):                    | $(OBJ_DIR) $(EXTLib) $(QMLLib) $(ARPACK_TARGET)
$(OBJ):                       | $(OBJ_DIR) $(EXTLib)

$(OBJ_DIR)/Initial_GWP_Basis_m.o: \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/param_WP0_m.o

$(OBJ_DIR)/HamiltonianSpectrum_m.o: \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/Ana_psi_m.o \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/psi_m.o

$(OBJ_DIR)/TISE_m.o: \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Basis_m.o

$(OBJ_DIR)/HamiltonianSpectrum_m.o: $(OBJ_DIR)/TISE_m.o
$(OBJ_DIR)/$(MAIN).o: $(OBJ_DIR)/TISE_m.o
$(OBJ_DIR)/sub_propa_m.o: $(OBJ_DIR)/TISE_m.o

$(OBJ_DIR)/HagedornDiracFrenkel_m.o: \
    $(OBJ_DIR)/HagedornQuadrature_m.o \
    $(OBJ_DIR)/HagedornLinalg_m.o \
    $(OBJ_DIR)/NDindex_hagedorn_m.o \
    $(OBJ_DIR)/HagedornBasis_m.o

$(OBJ_DIR)/Ana_psi_m.o: \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/NDindex_m.o \
    $(OBJ_DIR)/psi_m.o

$(OBJ_DIR)/Auto_corr_m.o: \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/Hagedorn_m.o \
    $(OBJ_DIR)/sub_propa_m.o

$(OBJ_DIR)/SplitOperator_m.o: \
    $(OBJ_DIR)/NUFFT_m.o

$(OBJ_DIR)/Basis_m.o: \
    $(OBJ_DIR)/NUFFT_m.o \
    $(OBJ_DIR)/NDindex_m.o \
    $(OBJ_DIR)/poly0rtho_m.o \
    $(OBJ_DIR)/SmolyakSparseGridMetadata_m.o

$(OBJ_DIR)/Hagedorn_m.o: \
    $(OBJ_DIR)/HagedornBasis_m.o \
    $(OBJ_DIR)/HagedornLinalg_m.o \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/Ana_psi_m.o \
    $(OBJ_DIR)/poly0rtho_m.o \
    $(OBJ_DIR)/sub_propa_m.o \
    $(OBJ_DIR)/HagedornOverlap_m.o

$(OBJ_DIR)/HagedornBasis_m.o: \
    $(OBJ_DIR)/HagedornLinalg_m.o \
    $(OBJ_DIR)/NDindex_hagedorn_m.o

$(OBJ_DIR)/HagedornHamiltonian_m.o: \
    $(OBJ_DIR)/HagedornLinalg_m.o \
    $(OBJ_DIR)/HagedornBasis_m.o

$(OBJ_DIR)/HagedornOverlap_m.o: \
    $(OBJ_DIR)/HagedornLinalg_m.o \
    $(OBJ_DIR)/HagedornBasis_m.o \
    $(OBJ_DIR)/NDindex_hagedorn_m.o \
    $(OBJ_DIR)/poly0rtho_m.o

$(OBJ_DIR)/lanczos_m.o: \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Basis_m.o

$(OBJ_DIR)/Op_m.o: \
    $(OBJ_DIR)/SmolyakSparseBasisAlgebra_m.o \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/SmolyakSparseGridMetadata_m.o \
    $(OBJ_DIR)/Ana_psi_m.o \
    $(OBJ_DIR)/Molec_m.o \
    $(OBJ_DIR)/psi_m.o

$(OBJ_DIR)/Propa_m.o: \
    $(OBJ_DIR)/Checkpoint_m.o \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/Ana_psi_m.o \
    $(OBJ_DIR)/lanczos_m.o \
    $(OBJ_DIR)/Auto_corr_m.o \
    $(OBJ_DIR)/Hagedorn_m.o \
    $(OBJ_DIR)/sub_propa_m.o \
    $(OBJ_DIR)/Basis_m.o

$(OBJ_DIR)/Checkpoint_m.o: \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Basis_m.o

$(OBJ_DIR)/psi_m.o: \
    $(OBJ_DIR)/param_WP0_m.o \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/SmolyakSparseBasisAlgebra_m.o \
    $(OBJ_DIR)/SmolyakSparseGridMetadata_m.o \
    $(OBJ_DIR)/NDindex_m.o

$(OBJ_DIR)/NDindex_TestSupport_m.o: \
    $(OBJ_DIR)/NDindex_m.o

$(OBJ_DIR)/SmolyakSparseGridMetadata_m.o: \
    $(OBJ_DIR)/NDindex_m.o

$(OBJ_DIR)/SmolyakSparseBasisAlgebra_m.o: \
    $(OBJ_DIR)/SmolyakSparseGridMetadata_m.o

$(OBJ_DIR)/sub_propa_m.o: \
    $(OBJ_DIR)/SmolyakSparseBasisAlgebra_m.o \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/Ana_psi_m.o \
    $(OBJ_DIR)/lanczos_m.o \
    $(OBJ_DIR)/Sub_Vp_m.o \
    $(OBJ_DIR)/Basis_m.o

$(OBJ_DIR)/Sub_Vp_m.o: \
    $(OBJ_DIR)/SmolyakSparseBasisAlgebra_m.o \
    $(OBJ_DIR)/NDindex_m.o \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/poly0rtho_m.o

$(OBJ_DIR)/TD_SCHROD.o: \
    $(OBJ_DIR)/InputStaging_m.o \
    $(OBJ_DIR)/Initial_GWP_Basis_m.o \
    $(OBJ_DIR)/HamiltonianSpectrum_m.o \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/param_WP0_m.o \
    $(OBJ_DIR)/Propa_m.o \
    $(OBJ_DIR)/Ana_psi_m.o \
    $(OBJ_DIR)/lanczos_m.o \
    $(OBJ_DIR)/Sub_Vp_m.o

# =====================================================================================
# Create object directory
# =====================================================================================

$(OBJ_DIR):
	mkdir -p $@

#===============================================
#===============================================
#============= External libraries  =============
#===============================================
.PHONY: getlib
getlib: $(EXTLib)
#
$(EXTLib):
	$(MAKE) -C $(ExtLibDIR) -f $(MAIN_path)/scripts/makefile-extlib LIBA=$@ TD_SCHROD_ROOT=$(MAIN_path)

# QuantumModelLib is called directly by Tnum-Tana and must be built with the
# same compiler/options as the other external archive. It was present in
# LIBS/QMLLib but had no make rule, so a clean build stopped before compiling
# the application.
$(QMLLib):
	$(MAKE) -C $(ExtLibDIR) -f $(MAIN_path)/scripts/makefile-extlib LIBA=$@ TD_SCHROD_ROOT=$(MAIN_path)

.PHONY: arpack
arpack: $(ARPACK_TARGET)

ifeq ($(ARPACK),1)
$(ARPACK_LIB): $(MAIN_path)/scripts/install-arpack
	@chmod +x $(MAIN_path)/scripts/install-arpack
	$(MAIN_path)/scripts/install-arpack \
	  "$(ARPACK_ROOT)" "$(ARPACK_VERSION)" "$(ARPACK_PREFIX)" "$(FC)"
endif
# =====================================================================================
# Cleaning
# =====================================================================================
.PHONY: clean distclean cleanall cleanlocextlib
clean:
	rm -rf $(OBJ_DIR)/*.o
	rm -rf $(MOD_DIR)/*.mod

distclean: clean
	rm -f $(TARGET)


cleanall: distclean
	if [ "$(EXTLIB_LIST)" != "" -a "$(RECCLEAN)" = "1" ]; then ./scripts/cleanExtLib cleanall "$(ExtLibDIR)" "$(EXTLIB_LIST)" 0; fi
	@echo "  done remove the *.a libraries and the OBJ directory for "$(LIB_NAME)
cleanlocextlib: distclean
	if [ "$(EXTLIB_LIST)" != "" -a "$(RECCLEAN)" = "1" ] ; then ./scripts/cleanExtLib cleanlocextlib "$(ExtLibDIR)" "$(EXTLIB_LIST)" 0; fi
	@echo "  done remove all local library directories (..._loc) for "$(LIB_NAME)
# =====================================================================================
# Information
# =====================================================================================

info:
	@echo "Build Configuration:"
	@echo " Compiler : $(FC)"
	@echo " Target   : $(TARGET)"
	@echo " Flags    : $(FFLAGS)"
	@echo " Modules  : $(MOD_DIR)"
