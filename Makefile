# =====================================================================================
# Compiler configuration and flags
# =====================================================================================

# Fortran compiler (gfortran by default, option for ifort commented)
FC = gfortran
#FC = ifort

# Compilation flags (optimized for production)
#FFLAGS = -Wall -Wextra -O3 -fopenmp -J$(MOD_DIR)
#FFLAGS = -g -O0 -Wall -Wextra -fcheck=all -fbacktrace -fopenmp -J$(MOD_DIR)
#FFLAGS = -Wall -Wextra -Wimplicit-interface -Wno-unused -fcheck=all -fbacktrace -fopenmp -J$(MOD_DIR)
 
# Alternative flags for debugging (uncomment if needed)
#FFLAGS = -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan -fopenmp -J$(MOD_DIR)
#FFLAGS= -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan -fopenmp -J$(MOD_DIR)
#FFLAGS= -Wall -Wextra -Wimplicit-interface -fPIC -O3 -march=native -ffast-math -funroll-loops -fopenmp -J$(MOD_DIR)
#FFLAGS= -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=5 -g -fcheck=all -fbacktrace -Og -g -fcheck=bounds -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan -fopenmp -J$(MOD_DIR)

FFLAGS = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=5 \
         -g -fcheck=all -fbacktrace -Og -fcheck=bounds -fwhole-file \
         -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan \
         -fopenmp -Warray-temporaries -J$(MOD_DIR)

OBJ_DIR=obj

# =====================================================================================
# Directory structure
# =====================================================================================
OBJ_DIR = obj
MOD_DIR = obj
SRC_DIR = src
MAIN_DIR = app

# =====================================================================================
# External library paths
# =====================================================================================
QML_DIR = EXT_Lib/QuantumModelLib
QDU_DIR = $(QML_DIR)/Ext_Lib/QDUtilLib
MOD_EXT = $(QDU_DIR)/OBJ/obj_gfortran_opt1_omp1_lapack1_int4

FFLAGS += -I$(MOD_EXT)

# =====================================================================================
# Main program configuration
# =====================================================================================
MAIN = TD_SCHROD
TARGET = $(MAIN).x

# =====================================================================================
# Source and object files
# =====================================================================================
LIBSRC = Molec_m.f90 NDindex_m.f90 poly0rtho_m.f90 Basis_m.f90 psi_m.f90 \
         Ana_psi_m.f90 lanczos_m.f90 Op_m.f90 param_WP0_m.f90 Auto_corr_m.f90 \
         Propa_m.f90 Hagedorn_m.f90 sub_propa_m.f90 Sub_Vp_m.f90 Sub_Vp_m.f90

# Generate object file list
OBJ0 = $(LIBSRC:.f90=.o)
OBJ = $(addprefix $(OBJ_DIR)/, $(OBJ0))

# Display build information (debug)
$(info [Build info] Object files: $(OBJ))

# =====================================================================================
# External libraries
# =====================================================================================
LIB = $(QML_DIR)/libQMLibFull_gfortran_opt1_omp1_lapack1_int4.a \
      $(QDU_DIR)/libQD_gfortran_opt1_omp1_lapack1_int4.a

LDFLAGS = -llapack -lblas

# =====================================================================================
# Main build rules
# =====================================================================================
.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ_DIR)/$(MAIN).o $(OBJ) $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/$(MAIN).o: $(MAIN_DIR)/$(MAIN).f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -o $@ -c $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -o $@ -c $<

# =====================================================================================
# External library verification
# =====================================================================================
$(QML_DIR)/libQMLibFull_gfortran_opt1_omp1_lapack1_int4.a:
	@test -e $@ || { echo "Error: Library $@ not found!"; exit 1; }

$(QDU_DIR)/libQD_gfortran_opt1_omp1_lapack1_int4.a:
	@test -e $@ || { echo "Error: Library $@ not found!"; exit 1; }

# =====================================================================================
# Create necessary directories
# =====================================================================================
$(OBJ_DIR):
	@mkdir -p $@
 
# =====================================================================================
# Module dependencies
# =====================================================================================
$(OBJ_DIR)/$(MAIN).o: $(OBJ)

$(OBJ_DIR)/Basis_m.o: $(OBJ_DIR)/NDindex_m.o $(OBJ_DIR)/poly0rtho_m.o
$(OBJ_DIR)/Hagedorn_m.o: $(OBJ_DIR)/sub_propa_m.o $(OBJ_DIR)/Op_m.o \
                        $(OBJ_DIR)/psi_m.o $(OBJ_DIR)/Basis_m.o \
                        $(OBJ_DIR)/poly0rtho_m.o $(OBJ_DIR)/Ana_psi_m.o

$(OBJ_DIR)/Sub_Vp_m.o: $(OBJ_DIR)/Op_m.o $(OBJ_DIR)/psi_m.o \
                      $(OBJ_DIR)/Basis_m.o $(OBJ_DIR)/NDindex_m.o

$(OBJ_DIR)/psi_m.o: $(OBJ_DIR)/Basis_m.o $(OBJ_DIR)/param_WP0_m.o
$(OBJ_DIR)/Ana_psi_m.o: $(OBJ_DIR)/psi_m.o
$(OBJ_DIR)/Auto_corr_m.o: $(OBJ_DIR)/psi_m.o $(OBJ_DIR)/Hagedorn_m.o \
                         $(OBJ_DIR)/Basis_m.o
$(OBJ_DIR)/Op_m.o: $(OBJ_DIR)/Molec_m.o $(OBJ_DIR)/Basis_m.o $(OBJ_DIR)/psi_m.o
$(OBJ_DIR)/lanczos_m.o: $(OBJ_DIR)/Op_m.o $(OBJ_DIR)/psi_m.o
$(OBJ_DIR)/Propa_m.o: $(OBJ_DIR)/sub_propa_m.o $(OBJ_DIR)/Auto_corr_m.o \
                     $(OBJ_DIR)/lanczos_m.o $(OBJ_DIR)/Op_m.o \
                     $(OBJ_DIR)/Ana_psi_m.o $(OBJ_DIR)/psi_m.o \
                     $(OBJ_DIR)/Hagedorn_m.o $(OBJ_DIR)/Basis_m.o
$(OBJ_DIR)/sub_propa_m.o: $(OBJ_DIR)/Sub_Vp_m.o $(OBJ_DIR)/lanczos_m.o \
                         $(OBJ_DIR)/Op_m.o $(OBJ_DIR)/psi_m.o \
                         $(OBJ_DIR)/Ana_psi_m.o $(OBJ_DIR)/Basis_m.o

# -------------------------------
# 8. Utility Targets
# -------------------------------
clean:
	rm -rf $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod

distclean: clean
	rm -rf $(BIN_DIR)/*.x build

info:
	@echo "Build Configuration:"
	@echo "  Compiler: $(FC)"
	@echo "  Build Type: $(BUILD)"
	@echo "  Flags: $(FFLAGS)"
	@echo "  Target: $(TARGET)"

# -------------------------------
# 9. Help System
# -------------------------------
help:
	@echo "Fortran Project Build System"
	@echo "Targets:"
	@echo "  all       - Build main executable (default)"
	@echo "  debug     - Build with debug flags (BUILD=debug)"
	@echo "  clean     - Remove object files"
	@echo "  distclean - Remove all generated files"
	@echo "  info      - Show build configuration"
	@echo ""
	@echo "Environment variables:"
	@echo "  FC=ifort   - Use Intel Fortran compiler"
	@echo "  BUILD=debug - Debug configuration"

# ================================================================
# End of Makefile
# ================================================================