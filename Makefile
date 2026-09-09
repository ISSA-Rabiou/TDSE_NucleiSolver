MAKEFLAGS += --no-builtin-rules --no-builtin-variables

FC := gfortran
#FC = ifor

OBJ_DIR = obj
MOD_DIR = obj
SRC_DIR = src
MAIN_DIR = app


FFLAGS := -cpp -Wall -Wextra -Wimplicit-interface \
          -g -Og \
          -fopenmp \
          -fcheck=all -fbacktrace \
          -finit-real=nan -finit-integer=-999999 \
          -ffree-line-length-none \
          -Warray-temporaries \
          -J$(MOD_DIR)

# =====================================================================================
# External library paths
# =====================================================================================
QML_DIR = EXT_Lib/QuantumModelLib
QDU_DIR = $(QML_DIR)/Ext_Lib/QDUtilLib
MOD_EXT = $(QDU_DIR)/OBJ/obj_gfortran_opt1_omp1_lapack1_int4_real64

TNUM_ROOT = EXT_Lib/Tnum-Tana
TNUM_LIB  = $(TNUM_ROOT)/libTnum-TanaFull.a
EVRT_LIB  = $(TNUM_ROOT)/Ext_Lib/EVRT_dnSVM_loc/libEVRT_dnSVM.a

MOD_TNUM  = $(TNUM_ROOT)/OBJ/obj_gfortran_opt1_omp1_lapack1_int4_real64
MOD_CONST = $(TNUM_ROOT)/Ext_Lib/ConstPhys/OBJ/obj_gfortran_opt1_omp1_lapack1_int4_real64
MOD_EVRT  = $(TNUM_ROOT)/Ext_Lib/EVRT_dnSVM_loc/OBJ/obj_gfortran_opt1_omp1_lapack1_int4_real64

FFLAGS += -I$(MOD_EXT) -I$(MOD_TNUM) -I$(MOD_CONST) -I$(MOD_EVRT)

LIBS = -L$(QML_DIR) -lQMLFull \
       -L$(QDU_DIR) -lQDUtilLib \
       -L$(TNUM_ROOT) -lTnum-TanaFull \
       -L$(TNUM_ROOT)/Ext_Lib/EVRT_dnSVM_loc -lEVRT_dnSVM \
       -lblas -llapack

LDFLAGS = $(LIBS)

# =====================================================================================
# Main program
# =====================================================================================
MAIN = TD_SCHROD
TARGET = $(MAIN).x

# =====================================================================================
# Source and object files
# =====================================================================================

LIBSRC = Molec_m.f90 NDindex_m.f90 poly0rtho_m.f90 \
         SmolyakSparseGridMetadata_m.f90 Basis_m.f90 psi_m.f90 \
         Ana_psi_m.f90 lanczos_m.f90 Op_m.f90 param_WP0_m.f90 \
         Auto_corr_m.f90 Propa_m.f90 Hagedorn_m.f90 \
         sub_propa_m.f90 Sub_Vp_m.f90 \
         HagedornLinalg_m.f90 HagedornOverlap_m.f90 \
         NDindex_hagedorn_m.f90 HagedornBasis_m.f90 \
         HagedornHamiltonian_m.f90

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

$(TARGET): $(OBJ_DIR)/$(MAIN).o $(OBJ) \
           $(QML_DIR)/libQMLFull.a \
           $(QDU_DIR)/libQDUtilLib.a \
           $(TNUM_LIB) \
           $(EVRT_LIB)
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

# ============================================================
# Dépendances manuelles entre modules (format $(OBJ_DIR)/...)
# ============================================================

$(OBJ_DIR)/Ana_psi_m.o: \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/NDindex_m.o \
    $(OBJ_DIR)/psi_m.o

$(OBJ_DIR)/Auto_corr_m.o: \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/Hagedorn_m.o \
    $(OBJ_DIR)/sub_propa_m.o

$(OBJ_DIR)/Basis_m.o: \
    $(OBJ_DIR)/NDindex_m.o \
    $(OBJ_DIR)/poly0rtho_m.o \
    $(OBJ_DIR)/SmolyakSparseGridMetadata_m.o

$(OBJ_DIR)/Hagedorn_m.o: \
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

$(OBJ_DIR)/HagedornDiracFrenkel_m.o: \
    $(OBJ_DIR)/HagedornLinalg_m.o \
    $(OBJ_DIR)/NDindex_hagedorn_m.o \
    $(OBJ_DIR)/HagedornBasis_m.o

$(OBJ_DIR)/HagedornHamiltonian_m.o: \
    $(OBJ_DIR)/HagedornLinalg_m.o \
    $(OBJ_DIR)/HagedornBasis_m.o

$(OBJ_DIR)/HagedornOverlap_m.o: \
    $(OBJ_DIR)/HagedornBasis_m.o \
    $(OBJ_DIR)/NDindex_hagedorn_m.o \
    $(OBJ_DIR)/poly0rtho_m.o

$(OBJ_DIR)/lanczos_m.o: \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Basis_m.o

$(OBJ_DIR)/Op_m.o: \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/SmolyakSparseGridMetadata_m.o \
    $(OBJ_DIR)/Ana_psi_m.o \
    $(OBJ_DIR)/Molec_m.o \
    $(OBJ_DIR)/psi_m.o

$(OBJ_DIR)/Propa_m.o: \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/Ana_psi_m.o \
    $(OBJ_DIR)/lanczos_m.o \
    $(OBJ_DIR)/Auto_corr_m.o \
    $(OBJ_DIR)/Hagedorn_m.o \
    $(OBJ_DIR)/sub_propa_m.o \
    $(OBJ_DIR)/Basis_m.o

$(OBJ_DIR)/psi_m.o: \
    $(OBJ_DIR)/param_WP0_m.o \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/SmolyakSparseGridMetadata_m.o \
    $(OBJ_DIR)/NDindex_m.o

$(OBJ_DIR)/SmolyakSparseGridMetadata_m.o: \
    $(OBJ_DIR)/NDindex_m.o

$(OBJ_DIR)/sub_propa_m.o: \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/Ana_psi_m.o \
    $(OBJ_DIR)/lanczos_m.o \
    $(OBJ_DIR)/Sub_Vp_m.o \
    $(OBJ_DIR)/Basis_m.o

$(OBJ_DIR)/Sub_Vp_m.o: \
    $(OBJ_DIR)/NDindex_m.o \
    $(OBJ_DIR)/Op_m.o \
    $(OBJ_DIR)/Basis_m.o \
    $(OBJ_DIR)/psi_m.o \
    $(OBJ_DIR)/poly0rtho_m.o

$(OBJ_DIR)/TD_SCHROD.o: \
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

# =====================================================================================
# External library verification
# =====================================================================================

$(QML_DIR)/libQMLFull.a:
	@test -e $@ || { echo "Error: Library $@ not found!"; exit 1; }

$(QDU_DIR)/libQDUtilLib.a:
	@test -e $@ || { echo "Error: Library $@ not found!"; exit 1; }

$(TNUM_LIB):
	@test -e $@ || { echo "Error: Library $@ not found!"; exit 1; }

$(EVRT_LIB):
	@test -e $@ || { echo "Error: Library $@ not found!"; exit 1; }

# =====================================================================================
# Cleaning
# =====================================================================================

clean:
	rm -rf $(OBJ_DIR)/*.o
	rm -rf $(MOD_DIR)/*.mod

distclean: clean
	rm -f $(TARGET)

# =====================================================================================
# Information
# =====================================================================================

info:
	@echo "Build Configuration:"
	@echo " Compiler : $(FC)"
	@echo " Target   : $(TARGET)"
	@echo " Flags    : $(FFLAGS)"
	@echo " Modules  : $(MOD_DIR)"