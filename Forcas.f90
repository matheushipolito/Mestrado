MODULE Forcas
    IMPLICIT NONE
    PUBLIC
    TYPE Type_Forcas
        INTEGER*4                            :: ID_Node
        REAL*8                               :: F_X
        REAL*8                               :: F_Y
        REAL*8                               :: F_Z 
    END TYPE Type_Forcas
    !===============================================================================================
    CONTAINS
    !===============================================================================================   
    
    SUBROUTINE Ler_Forcas(Self,unidade)                                           !Self são as forças
    TYPE(Type_Forcas),        INTENT(INOUT)   :: Self
    INTEGER*4                                 :: ID_Node
    REAL*8                                    :: F_X
    REAL*8                                    :: F_Y
    REAL*8                                    :: F_Z   
    INTEGER*8                                 ::  unidade
    READ(unidade,*) ID_Node, F_X, F_Y, F_Z
    CALL SET_Forcas(Self, ID_Node, F_X, F_Y, F_Z)
    END SUBROUTINE Ler_Forcas
    
    !================================================================================================
    SUBROUTINE SET_Forcas(Self, ID_Node,F_X, F_Y, F_Z)
    TYPE(Type_Forcas),        INTENT(INOUT)   :: Self
    INTEGER*4                                 :: ID_Node
    REAL*8                                    :: F_X
    REAL*8                                    :: F_Y
    REAL*8                                    :: F_Z
    Self%ID_Node=ID_Node
    Self%F_X=F_X
    Self%F_Y=F_Y
    Self%F_Z=F_Z
    END SUBROUTINE SET_Forcas
    !================================================================================================
    FUNCTION GET_Fext(Self,i)   RESULT(j)
    TYPE(Type_Forcas),    INTENT(IN)          ::  Self
    INTEGER*4                                 :: i
    REAL*8                                    :: j
    REAL*8,      DIMENSION(3)                 :: k
    k(1)=Self%F_X
    k(2)=Self%F_Y
    k(3)=Self%F_Z
    j=k(i)
    END FUNCTION GET_Fext
    !================================================================================================
    FUNCTION GET_Fext_ID(Self)   RESULT(j)
    TYPE(Type_Forcas),    INTENT(IN)          :: Self
    INTEGER*4                                 :: j
    j=Self%ID_Node
    END FUNCTION GET_Fext_ID
    
END MODULE Forcas