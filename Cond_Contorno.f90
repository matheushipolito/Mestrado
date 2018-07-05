MODULE Cond_Contorno
    IMPLICIT NONE
    PUBLIC 
    TYPE Type_Cond_Contorno
        INTEGER*4                            :: ID_Node
        INTEGER*4                            :: R_X
        REAL*8                               :: D_X
        INTEGER*4                            :: R_Y
        REAL*8                               :: D_Y
        INTEGER*4                            :: R_Z
        REAL*8                               :: D_Z
    END TYPE Type_Cond_Contorno
    !===============================================================================================
    CONTAINS
    !===============================================================================================   
    
    SUBROUTINE Ler_Cond_Contorno(Self,unidade)                            !Self são as cond de contorno
    TYPE(Type_Cond_Contorno), INTENT(INOUT)   :: Self
    INTEGER*4                                 :: ID_Node
    INTEGER*4                                 :: R_X
    REAL*8                                    :: D_X
    INTEGER*4                                 :: R_Y
    REAL*8                                    :: D_Y
    INTEGER*4                                 :: R_Z
    REAL*8                                    :: D_Z
    INTEGER*8                                 ::  unidade
    READ(unidade,*) ID_Node, R_X, D_X, R_Y, D_Y, R_Z, D_Z
    CALL SET_Cond_Contorno(Self, ID_Node, R_X, D_X, R_Y, D_Y, R_Z, D_Z)
    END SUBROUTINE Ler_Cond_Contorno
    
    !================================================================================================
    SUBROUTINE SET_Cond_Contorno(Self, ID_Node, R_X, D_X, R_Y, D_Y, R_Z, D_Z)
    TYPE(Type_Cond_Contorno), INTENT(INOUT)   :: Self
    INTEGER*4                                 :: ID_Node
    INTEGER*4                                 :: R_X
    REAL*8                                    :: D_X
    INTEGER*4                                 :: R_Y
    REAL*8                                    :: D_Y
    INTEGER*4                                 :: R_Z
    REAL*8                                    :: D_Z
    Self%ID_Node=ID_Node
    Self%R_X=R_X
    Self%D_X=D_X
    Self%R_Y=R_Y
    Self%D_Y=D_Y
    Self%R_Z=R_Z
    Self%D_Z=D_Z
    END SUBROUTINE SET_Cond_Contorno
    !================================================================================================
    FUNCTION GET_Cond_Cont(Self,i)   RESULT(j)
    TYPE(Type_Cond_Contorno),    INTENT(IN)   ::  Self
    INTEGER*4                                 :: i
    INTEGER*4                                 :: j
    INTEGER*4,   DIMENSION(3)                 :: k
    k(1)=Self%R_X
    k(2)=Self%R_Y
    k(3)=Self%R_Z
    j=k(i)
    END FUNCTION GET_Cond_Cont
    !================================================================================================
    FUNCTION GET_Cond_Cont_desl(Self)   RESULT(k)
    TYPE(Type_Cond_Contorno),    INTENT(IN)   ::  Self
    INTEGER*4                                 :: i
    INTEGER*4                                 :: j
    REAL*8,   DIMENSION(3)                 :: k
    k(1)=Self%D_X
    k(2)=Self%D_Y
    k(3)=Self%D_Z

    END FUNCTION GET_Cond_Cont_desl
    
END MODULE Cond_Contorno