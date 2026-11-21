MODULE SmolyakSparseBasisAlgebra_m
  USE QDUtil_m
  USE SmolyakSparseGridMetadata_m, ONLY: SmolyakSparseGridDescriptor_t, UnpackSmolyakTermsToFullBasis
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: DotSmolyakBasis
CONTAINS

  FUNCTION DotSmolyakBasis(a, b, Descriptor) RESULT(value)
    COMPLEX(KIND=Rkind), INTENT(IN) :: a(:), b(:)
    TYPE(SmolyakSparseGridDescriptor_t), INTENT(IN) :: Descriptor
    COMPLEX(KIND=Rkind) :: value
    COMPLEX(KIND=Rkind), ALLOCATABLE :: fa(:), fb(:)
    ALLOCATE(fa(Descriptor%nb*MAX(1,Descriptor%nsurf)))
    ALLOCATE(fb(SIZE(fa)))
    CALL UnpackSmolyakTermsToFullBasis(fa,a,Descriptor)
    CALL UnpackSmolyakTermsToFullBasis(fb,b,Descriptor)
    value=DOT_PRODUCT(fa,fb)
    DEALLOCATE(fa,fb)
  END FUNCTION DotSmolyakBasis


END MODULE SmolyakSparseBasisAlgebra_m
