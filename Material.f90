MODULE Material
IMPLICIT NONE
    !===================================================================================================
    TYPE Type_Material
        PRIVATE
        INTEGER*4                              :: ID
        REAL*8                                 :: Young
        REAL*8                                 :: Densidade
    END TYPE Type_Material
    !===============================================================================================
    CONTAINS
    !===============================================================================================
    SUBROUTINE Ler_Material(Self,unidade)                                           !Self é o tipo de material
    IMPLICIT NONE
    TYPE(Type_Material),    INTENT(INOUT)      :: Self
    INTEGER*4                                  :: ID
    REAL*8                                     :: Young
    REAL*8                                     :: Densidade
    INTEGER*8                                 ::  unidade
    READ(unidade,*) ID, Young, Densidade
    CALL SET_Material(Self, ID, Young, Densidade)
    END SUBROUTINE Ler_Material
    !================================================================================================
    SUBROUTINE SET_Material(Self, ID, Young, Densidade)
    TYPE(Type_Material),    INTENT(INOUT)      :: Self
    INTEGER*4                                  :: ID
    REAL*8                                     :: Young
    REAL*8                                     :: Densidade
    Self%ID=ID
    Self%Young=Young
    Self%Densidade=Densidade
    END SUBROUTINE SET_Material
    ! ========================================================================================
    FUNCTION GET_material_ID(Self)          RESULT(ID)
    TYPE(Type_material), INTENT(IN)     ::  Self
    INTEGER*4                           ::  ID
    ID = Self%ID
    END FUNCTION GET_material_ID
    ! ================================================================================
    FUNCTION GET_material_Young(Self)       RESULT(Young)
    TYPE(Type_material), INTENT(IN)     ::  Self
    REAL*8                              ::  Young
    Young = Self%Young
    END FUNCTION GET_material_young
    ! =================================================================================
    FUNCTION GET_material_Densidade(Self)   RESULT(Densidade)
    TYPE(Type_material), INTENT(IN)     ::  Self
    REAL*8                              ::  Densidade
    Densidade =  Self%Densidade
    END FUNCTION GET_material_Densidade
    ! =================================================================================

    END MODULE Material
    
    
