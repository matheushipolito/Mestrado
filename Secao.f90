MODULE Secao
    IMPLICIT NONE
    PUBLIC
    ! =======================================================================================
    TYPE Type_Secao 
        PRIVATE
        INTEGER*4   ::  ID
        REAL*8      ::  Area    
    END TYPE Type_Secao
    ! =======================================================================================
CONTAINS
    ! =======================================================================================
        SUBROUTINE Ler_Secao(Self,unidade)
            IMPLICIT NONE
            TYPE(Type_Secao), INTENT(INOUT)   ::    Self
            INTEGER*8                                 ::  unidade
            INTEGER*4   ::  ID
            REAL*8      ::   Area
            READ(unidade,*)   ID, Area
            CALL    SET_Secao(Self, ID, Area)
        END SUBROUTINE Ler_Secao
  ! =======================================================================================
        SUBROUTINE SET_Secao(Self, ID, Area)
            IMPLICIT NONE
            TYPE(Type_Secao), INTENT(INOUT)   ::   Self
            INTEGER*4   ::   ID
            REAL*8      ::   Area
            Self%Id = ID
            Self%Area = Area       
        END SUBROUTINE SET_Secao
  ! =======================================================================================
        FUNCTION GET_Secao_Area(Self)    RESULT(Area)
        IMPLICIT NONE
        TYPE(Type_Secao), INTENT(IN)         ::   Self
        REAL*8                               ::   Area
        Area = Self%Area
        END FUNCTION GET_Secao_Area
  ! =======================================================================================
        FUNCTION GET_Secao_ID(Self)    RESULT(ID)
        IMPLICIT NONE
        TYPE(Type_Secao), INTENT(IN)         ::   Self
        REAL*8                               ::   ID
        ID = Self%ID
        END FUNCTION GET_Secao_ID
END MODULE Secao