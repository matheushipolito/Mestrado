MODULE Element
    USE Material
    USE Secao  
    USE Nodes
    IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------
    TYPE Type_Element
       
        INTEGER*4                             ::  Element_ID               ! ID número do elemento
        INTEGER*4                             ::  Material_ID              ! ID do material atribuido ao Element
        INTEGER*4                             ::  Secao_ID                 ! ID da seção atribuido ao Element
        INTEGER*4                             ::  Nodes_by_Element = 2     ! Número de nós por Element. Default: 2
        REAL*8                                ::  L0                       ! Comprimento inicial, L0
        REAL*8                                ::  Lf                       ! Comprimento atual, L
        REAL*8                                ::  Alpha                    ! Cosseno giratório em X
        REAL*8                                ::  Beta                     ! Cosseno giratório em Y
        REAL*8                                ::  Gama                     ! Cosseno giratório em Z
        REAL*8                                ::  Green                    ! Deformação de Green
        REAL*8                                ::  Massa                    ! Massa do elemento
        REAL*8                                ::  Stress                   ! Tensão no elemento
        
        
        INTEGER*4, DIMENSION(2)               ::  Incidencia
        REAL*8,    DIMENSION(:,:) , POINTER   ::  M                        ! Matriz de massa do elemento
        REAL*8,    DIMENSION(:,:) , POINTER   ::  H                        ! Matriz hessiana do elemento
        REAL*8,    DIMENSION(:,:) , POINTER   ::  D                        ! Matriz hessiana do elemento quarta ordem
        REAL*8,    DIMENSION(:,:) , POINTER   ::  B                        ! Matriz tensão deformação do elemento
        REAL*8,    DIMENSION(:,:) , POINTER   ::  X                        ! Matriz deslocamento específico do elemento
        REAL*8,    DIMENSION(3)               ::  Primeira_Derivada_Node_1 ! Derivada 1 do Nó 1 (local) do elemento 
        REAL*8,    DIMENSION(3)               ::  Primeira_Derivada_Node_2 ! Derivada 1 do Nó 2 (local) do elemento 
        REAL*8,    DIMENSION(6)               ::  Fi                       ! Força interna do elemento
        
        
        TYPE(Type_Nodes),POINTER              ::  Node1                     ! Apontador dos nós do elemento
        TYPE(Type_Nodes),POINTER              ::  Node2                     ! Apontador dos nós do elemento
        TYPE(Type_Material)      ,POINTER     ::  Material                 ! Apontador do material
        TYPE(Type_Secao)         ,POINTER     ::  Secao                    ! Apontador da seção
        
        
        
    END TYPE Type_Element
    
!----------------------------------------------------------------------------------------------------------
    CONTAINS
!----------------------------------------------------------------------------------------------------------
    
    SUBROUTINE Ler_Element(Self,unidade)              ! A subrotina atribui as propiedades do material e da seção ao Element
    IMPLICIT NONE
    TYPE(Type_Element),    INTENT(INOUT)    ::     Self
    INTEGER*4                               ::     Element_ID
    INTEGER*4                               ::     Material_ID
    INTEGER*4                               ::     Secao_ID
    INTEGER*8                                 ::  unidade
    READ(unidade,*)  Element_ID, Material_ID, Secao_ID
    CALL SET_Element(Self, Element_ID, Material_ID, Secao_ID)
    RETURN
    END SUBROUTINE Ler_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_RIGIDEZ_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(INOUT)    ::    Self
    REAL*8, DIMENSION(6,6)                  ::    H
    REAL*8                                  ::    L0
    REAL*8                                  ::    Area, E_Y
    Area = GET_Secao_Area(Self%Secao)
    E_Y = GET_material_Young(Self%Material)
    L0 = Self%L0  
    
    Self%H(1,1) =  Self%Alpha**2.0d0
    Self%H(1,2) =  Self%Alpha*Self%Beta 
    Self%H(1,3) =  Self%Alpha*Self%Gama 
    Self%H(1,4) = -Self%H(1,1) 
    Self%H(1,5) = -Self%H(1,2) 
    Self%H(1,6) = -Self%H(1,3)
    Self%H(2,1) =  Self%H(1,2)
    Self%H(2,2) =  Self%Beta**2.0d0 
    Self%H(2,3) =  Self%Beta*Self%Gama
    Self%H(2,4) = -Self%H(1,2)
    Self%H(2,5) = -Self%H(2,2)
    Self%H(2,6) = -Self%H(2,2)
    Self%H(3,1) =  Self%Alpha*Self%Gama 
    Self%H(3,2) =  Self%H(2,3) 
    Self%H(3,3) =  Self%Gama**2.0d0 
    Self%H(3,4) = -Self%H(1,3)
    Self%H(3,5) = -Self%H(1,3)
    Self%H(3,6) = -Self%H(3,3)
    Self%H(4,1) = -Self%H(1,1) 
    Self%H(4,2) = -Self%H(1,2)
    Self%H(4,3) = -Self%H(1,3)
    Self%H(4,4) =  Self%H(1,1)
    Self%H(4,5) =  Self%H(1,2)
    Self%H(4,6) =  Self%H(1,3) 
    Self%H(5,1) = -Self%H(1,2)
    Self%H(5,2) = -Self%H(2,2)
    Self%H(5,3) = -Self%H(2,3)
    Self%H(5,4) =  Self%H(1,2)
    Self%H(5,5) =  Self%H(2,2)
    Self%H(5,6) =  Self%H(2,3)
    Self%H(6,1) =  Self%H(1,6)
    Self%H(6,2) =  Self%H(2,6)
    Self%H(6,3) =  Self%H(3,6)
    Self%H(6,4) =  Self%H(4,6)
    Self%H(6,5) =  Self%H(5,6)
    Self%H(6,6) =  Self%H(3,3)
    
    Self%H = Self%H*E_y*Area/L0
    END SUBROUTINE SET_RIGIDEZ_Element
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_M_MASSA_Element(Self)   RESULT(M) 
    IMPLICIT NONE
    TYPE(Type_Element)       ::    Self
    REAL*8,     DIMENSION(6,6)          ::    M
    
    M=Self%M
    END FUNCTION GET_M_MASSA_Element
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_RIGIDEZ_Element(Self)   RESULT(H) 
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(IN)       ::    Self
    REAL*8,    DIMENSION(6,6)              ::    H
    H=Self%H
    END FUNCTION GET_RIGIDEZ_Element
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_Fi_Element(Self)   RESULT(Fi) 
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(IN)        ::    Self
    REAL*8,    DIMENSION(6)                 ::    Fi
    Fi=Self%Fi
    END FUNCTION GET_Fi_Element
!----------------------------------------------------------------------------------------------------------
 SUBROUTINE SET_Peso_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(INOUT)    ::     Self
    REAL*8                                 ::     Area
    REAL*8                                 ::     Lf
    REAL*8                                 ::     Densidade
    
    Area = GET_Secao_Area(Self%Secao)
    Densidade = GET_material_densidade(Self%Material)
    Lf = Self%Lf 
    
    Self%Massa = Densidade * Lf * Area
    
    END SUBROUTINE SET_Peso_Element
 !----------------------------------------------------------------------------------------------------------
    FUNCTION GET_Peso_Element(Self)    RESULT(Massa)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(IN)          ::   Self
    REAL*8                                  ::   Massa
    Massa = Self%Massa
    END FUNCTION GET_Peso_Element
 !----------------------------------------------------------------------------------------------------------
    FUNCTION GET_Stress_Element(Self)    RESULT(Stress)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(IN)          ::   Self
    REAL*8                                  ::   Stress
    Stress = Self%Stress
    END FUNCTION GET_Stress_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Alpha_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(INOUT)    ::     Self
    REAL*8                                  ::     x_i
    REAL*8                                  ::     x_f
    x_i = GET_x0(Self%Node1)
    x_f = GET_x0(Self%Node2)
    Self%Alpha = (x_f-x_i)/Self%L0
    END SUBROUTINE SET_Alpha_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Beta_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(INOUT)    ::     Self
    REAL*8                                  ::     y_i
    REAL*8                                  ::     y_f
    y_i = GET_y0(Self%Node1)
    y_f = GET_y0(Self%Node2)
    Self%Beta = (y_f-y_i) / Self%L0
    END SUBROUTINE SET_Beta_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Gama_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(INOUT)    ::     Self
    REAL*8                                  ::     z_i
    REAL*8                                  ::     z_f
    z_i = GET_z0(Self%Node1)
    z_f = GET_z0(Self%Node2)
    Self%Gama = (z_f-z_i)/Self%L0
    END SUBROUTINE SET_Gama_Element
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE Atualizar_Alpha_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(INOUT)    ::     Self
    REAL*8                                  ::     x_i
    REAL*8                                  ::     x_f
    x_i = GET_xf(Self%Node1)
    x_f = GET_xf(Self%Node2)
    Self%Alpha = (x_f-x_i)/Self%Lf
    END SUBROUTINE Atualizar_Alpha_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE Atualizar_Beta_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(INOUT)    ::     Self
    REAL*8                                  ::     y_i
    REAL*8                                  ::     y_f
    y_i = GET_yf(Self%Node1)
    y_f = GET_yf(Self%Node2)
    Self%Beta = (y_f-y_i) / Self%Lf
    END SUBROUTINE Atualizar_Beta_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE Atualizar_Gama_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(INOUT)    ::     Self
    REAL*8                                  ::     z_i
    REAL*8                                  ::     z_f
    z_i = GET_zf(Self%Node1)
    z_f = GET_zf(Self%Node2)
    Self%Gama = (z_f-z_i)/Self%Lf
    END SUBROUTINE Atualizar_Gama_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Green_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(INOUT)     ::     Self
    REAL*8                                  ::     z_i
    Self%Green = 5.0d-1 * (Self%Lf**2.0d0 / Self%L0**2.0d0 - 1.0d0)
    END SUBROUTINE SET_Green_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_L0_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element)  ,   INTENT(INOUT)    :: Self
    REAL*8                                    :: x_i, y_i, z_i, x_f, y_f, z_f
    x_i = GET_x0(Self%Node1)
    y_i = GET_y0(Self%Node1)
    z_i = GET_z0(Self%Node1)
    x_f = GET_x0(Self%Node2)
    y_f = GET_y0(Self%Node2)
    z_f = GET_z0(Self%Node2)
    
    Self%L0 = dsqrt((x_f-x_i)**2.0d0+(y_f-y_i)**2.0d0+(z_f-z_i)**2.0d0)
    
    END SUBROUTINE SET_L0_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Lf_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element)  ,   INTENT(INOUT)    :: Self
    REAL*8                                    :: x_i, y_i, z_i, x_f, y_f, z_f
    x_i = Get_xf(Self%Node1)
    y_i = Get_yf(Self%Node1)
    z_i = Get_zf(Self%Node1)
    x_f = Get_xf(Self%Node2)
    y_f = Get_yf(Self%Node2)
    z_f = Get_zf(Self%Node2)
    
    Self%Lf = dsqrt((x_f-x_i)**2.0d0+(y_f-y_i)**2.0d0+(z_f-z_i)**2.0d0)
    
    END SUBROUTINE SET_Lf_Element
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_L0_Element(Self)    RESULT(L0)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(IN)         ::   Self
    REAL*8                                  ::   L0
    L0 = Self%L0
    END FUNCTION GET_L0_Element
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_Alpha_Element(Self)    RESULT(Alpha)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(IN)         ::   Self
    REAL*8                                  ::   Alpha
    Alpha = Self%Alpha
    END FUNCTION GET_Alpha_Element
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_Beta_Element(Self)    RESULT(Beta)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(IN)         ::   Self
    REAL*8                                  ::   Beta
    Beta = Self%Beta
    END FUNCTION GET_Beta_Element
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_Gama_Element(Self)    RESULT(Gama)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(IN)         ::   Self
    REAL*8                                  ::   Gama
    Gama = Self%Gama
    END FUNCTION GET_Gama_Element
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_Element_ID(Self)    RESULT(ID)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(IN)         ::   Self
    INTEGER*8                              ::   ID
    ID = Self%Element_ID
    END FUNCTION GET_Element_ID
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Element(Self, Element_ID, Material_ID, Secao_ID)
    IMPLICIT NONE
    TYPE(Type_Element),    INTENT(INOUT)   ::     Self
    INTEGER*4,              INTENT(IN)      ::     Element_ID
    INTEGER*4,              INTENT(IN)      ::     Material_ID
    INTEGER*4,              INTENT(IN)      ::     Secao_ID
    Self%Element_ID = Element_ID
    Self%Material_ID = Material_ID
    Self%Secao_ID = Secao_ID
    END SUBROUTINE SET_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Nodes_by_Element(Self, Nodes_by_Element)
    IMPLICIT NONE
    TYPE(Type_Element),    INTENT(INOUT)   ::    Self
    INTEGER*4,              INTENT(IN)      ::    Nodes_by_Element
    Self%Nodes_by_Element = Nodes_by_Element
    END SUBROUTINE SET_Nodes_by_Element
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_Element_material_ID(Self)    RESULT(Material_ID)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(IN)          ::   Self
    INTEGER*4                               ::   Material_ID
    Material_ID = Self%Material_ID
    END FUNCTION GET_Element_material_ID
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_Element_Secao_ID(Self)    RESULT(Secao_ID)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(IN)         ::   Self
    INTEGER*4                               ::   Secao_ID
    Secao_ID = Self%Secao_ID
    END FUNCTION GET_Element_Secao_ID
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Material_do_Element(Self, Material)                  
    IMPLICIT NONE
    TYPE(Type_Element)  ,   INTENT(INOUT)  ::   Self
    TYPE(Type_Material)  ,   POINTER        ::   Material
    Self%Material => Material
    END SUBROUTINE SET_Material_do_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Secao_do_Element(Self, Secao)                  
    IMPLICIT NONE
    TYPE(Type_Element)  ,   INTENT(INOUT)  ::   Self
    TYPE(Type_Secao)     ,   POINTER        ::   Secao
    Self%Secao => Secao
    END SUBROUTINE SET_Secao_do_Element 
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE Ler_Incidencia(Self,unidade)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(INOUT)      ::  Self
    INTEGER*4                               ::  ID
    INTEGER*4,  DIMENSION(2)                ::  Incidencia
    INTEGER*4                               ::  i
    INTEGER*4                               ::  unidade
    READ(unidade,*) ID, (Incidencia(i), i = 1, Self%Nodes_by_Element)
    CALL SET_Incidencia(Self, incidencia)
    END SUBROUTINE Ler_Incidencia
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Incidencia(Self, Incidencia)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(INOUT)      ::  Self
    INTEGER*4, DIMENSION(2)                :: Incidencia
    Self%Incidencia = Incidencia
    END SUBROUTINE SET_Incidencia
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_Incidencia(Self,i)    RESULT(Incidencia)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(IN)          ::   Self
    INTEGER*4                               ::   Incidencia
    INTEGER*4                               ::   i
    Incidencia = Self%Incidencia(i)
    END FUNCTION GET_Incidencia
!----------------------------------------------------------------------------------------------------------
    FUNCTION GET_Nodes_by_Element(Self)    RESULT(Nodes_by_Element)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(IN)         ::   Self
    INTEGER*4                               ::   Nodes_by_Element
    Nodes_by_Element = Self%Nodes_by_Element
    END FUNCTION GET_Nodes_by_Element
!!----------------------------------------------------------------------------------------------------------
    !SUBROUTINE ALLOCATE_Nodes_of_Element(Self)
    !TYPE(Type_Element),      INTENT(INOUT)      ::   Self
    !ALLOCATE(Self%Node(Self%Nodes_by_Element))
    !END SUBROUTINE ALLOCATE_Nodes_of_Element
!!----------------------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_M_Element(Self)
    TYPE(Type_Element),      INTENT(INOUT)      ::   Self
    ALLOCATE(Self%M(6,6))
    END SUBROUTINE ALLOCATE_M_Element
    !!----------------------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_H_Element(Self)
    TYPE(Type_Element),      INTENT(INOUT)      ::   Self
    ALLOCATE(Self%H(6,6))
    END SUBROUTINE ALLOCATE_H_Element
    !!----------------------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_B_Element(Self)
    TYPE(Type_Element),      INTENT(INOUT)      ::   Self
    ALLOCATE(Self%B(6,6))
    END SUBROUTINE ALLOCATE_B_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Nodes_of_Element(Self, Node1, node2)
    IMPLICIT NONE
    TYPE(Type_Element), INTENT(INOUT)       ::   Self
    TYPE(Type_Nodes),   POINTER             ::   Node1
    TYPE(Type_Nodes),   POINTER             ::   Node2
    !INTEGER*4                               ::   index
    Self%Node1 => Node1
    Self%Node2 => Node2
    END SUBROUTINE SET_Nodes_of_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Element_Fi(Self,delta_pos)
    IMPLICIT NONE
    TYPE(Type_Element),INTENT(INOUT)        ::  Self
    REAL*8, DIMENSION(6)                    ::  delta_pos

    Self%Fi(1) = - Self%Alpha*delta_pos(1)
    Self%Fi(2) = - Self%Beta*delta_pos(2)
    Self%Fi(3) = - Self%Gama*delta_pos(3)
    Self%Fi(4) = Self%Alpha*delta_pos(4)
    Self%Fi(5) = Self%Beta*delta_pos(5)
    Self%Fi(6) = Self%Gama*delta_pos(6)
    Self%Fi = Self%Fi*GET_material_Young(Self%Material)/Self%L0
    
    !Self%Stress = Self%Fi(1) + Self%Fi(2) + Self%Fi(3) + Self%Fi(4) + Self%Fi(5) + Self%Fi(6)

    END SUBROUTINE SET_Element_Fi
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Stress_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(INOUT)    ::     Self
    REAL*8                                 ::     L0
    REAL*8                                 ::     Lf
    REAL*8                                 ::     Young
    
    Young = GET_material_young(Self%Material)
    L0 = Self%L0 
    Lf = Self%Lf
    
    
    Self%Stress = (L0-Lf)*Young / L0
    !Self%Stress = Self%Fi(1) + Self%Fi(2) + Self%Fi(3) + Self%Fi(4) + Self%Fi(5) + Self%Fi(6)
    END SUBROUTINE SET_Stress_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Stress_Nao_Lin_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),   INTENT(INOUT)     ::  Self
    REAL*8                                  ::  Young
    REAL*8                                  ::  L0
    REAL*8                                  ::  Lf
    REAL*8                                  ::  Green
    REAL*8                                  ::  Area
    
    Young = GET_material_Young(Self%Material)
    L0 = Self%L0
    Lf = Self%Lf
    Green = Self%Green
    Area = GET_Secao_Area(Self%Secao)
    Self%Fi = 0.0d0
    
    Self%Stress =  Young* ((Lf**2.0d0 - L0**2.0d0)/(2.0d0*(L0**2.0d0)))
    
    END SUBROUTINE SET_Stress_Nao_Lin_Element
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Element_Fi_N_Linear(Self)
    IMPLICIT NONE
    TYPE(Type_Element),INTENT(INOUT)        ::  Self
    REAL*8                                  ::  Young
    REAL*8                                  ::  L0
    REAL*8                                  ::  Green
    REAL*8                                  ::  Area
    Young = GET_material_Young(Self%Material)
    L0 = Self%L0
    Green = Self%Green
    Area = GET_Secao_Area(Self%Secao)
    Self%Fi = 0.0d0
    
    Self%Fi(1) = Self%Primeira_Derivada_Node_1(1)
    Self%Fi(2) = Self%Primeira_Derivada_Node_1(2)
    Self%Fi(3) = Self%Primeira_Derivada_Node_1(3)
    Self%Fi(4) = Self%Primeira_Derivada_Node_2(1)
    Self%Fi(5) = Self%Primeira_Derivada_Node_2(2)
    Self%Fi(6) = Self%Primeira_Derivada_Node_2(3)
    Self%Fi = Self%Fi*Young*Green*Area*L0

    END SUBROUTINE SET_Element_Fi_N_Linear
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE SET_Primeira_Derivada(Self)
    TYPE(Type_Element), INTENT(INOUT)      ::   Self
    INTEGER*4                               ::   i
    INTEGER*4                               ::   j
            self%Primeira_Derivada_Node_1(1) = (GET_xf(Self%Node1) - GET_xf(Self%Node2)) /(SELF%L0**2.0d0)
            self%Primeira_Derivada_Node_1(2) = (GET_yf(Self%Node1) - GET_yf(Self%Node2)) /(SELF%L0**2.0d0)
            self%Primeira_Derivada_Node_1(3) = (GET_zf(Self%Node1) - GET_zf(Self%Node2)) /(SELF%L0**2.0d0)
            self%Primeira_Derivada_Node_2(1) = (GET_xf(Self%Node2) - GET_xf(Self%Node1)) /(SELF%L0**2.0d0)
            self%Primeira_Derivada_Node_2(2) = (GET_yf(Self%Node2) - GET_yf(Self%Node1)) /(SELF%L0**2.0d0)
            self%Primeira_Derivada_Node_2(3) = (GET_zf(Self%Node2) - GET_zf(Self%Node1)) /(SELF%L0**2.0d0)
    END SUBROUTINE SET_Primeira_Derivada
!----------------------------------------------------------------------------------------------------------    
    SUBROUTINE SET_HESSIAN_ELEMENT(Self)
    IMPLICIT NONE
    TYPE(Type_Element),  INTENT(INOUT)      ::  Self
    REAL*8                                  ::  Young
    
    Self%H = 0.0d0
    Self%H(1,1) = Self%Primeira_Derivada_Node_1(1) * Self%Primeira_Derivada_Node_1(1) + Self%Green / (Self%L0**2.0d0)
    Self%H(2,1) = Self%Primeira_Derivada_Node_1(2) * Self%Primeira_Derivada_Node_1(1)
    Self%H(3,1) = Self%Primeira_Derivada_Node_1(3) * Self%Primeira_Derivada_Node_1(1)
    Self%H(4,1) = Self%Primeira_Derivada_Node_2(1) * Self%Primeira_Derivada_Node_1(1) - Self%Green / (Self%L0**2.0d0)
    Self%H(5,1) = Self%Primeira_Derivada_Node_2(2) * Self%Primeira_Derivada_Node_1(1)
    Self%H(6,1) = Self%Primeira_Derivada_Node_2(3) * Self%Primeira_Derivada_Node_1(1)
    Self%H(1,2) = Self%Primeira_Derivada_Node_1(1) * Self%Primeira_Derivada_Node_1(2)
    Self%H(2,2) = Self%Primeira_Derivada_Node_1(2) * Self%Primeira_Derivada_Node_1(2) + Self%Green / (Self%L0**2.0d0)
    Self%H(3,2) = Self%Primeira_Derivada_Node_1(3) * Self%Primeira_Derivada_Node_1(2)
    Self%H(4,2) = Self%Primeira_Derivada_Node_2(1) * Self%Primeira_Derivada_Node_1(2)
    Self%H(5,2) = Self%Primeira_Derivada_Node_2(2) * Self%Primeira_Derivada_Node_1(2) - Self%Green / (Self%L0**2.0d0)
    Self%H(6,2) = Self%Primeira_Derivada_Node_2(3) * Self%Primeira_Derivada_Node_1(2)
    Self%H(1,3) = Self%Primeira_Derivada_Node_1(1) * Self%Primeira_Derivada_Node_1(3)
    Self%H(2,3) = Self%Primeira_Derivada_Node_1(2) * Self%Primeira_Derivada_Node_1(3)
    Self%H(3,3) = Self%Primeira_Derivada_Node_1(3) * Self%Primeira_Derivada_Node_1(3) + Self%Green / (Self%L0**2.0d0)
    Self%H(4,3) = Self%Primeira_Derivada_Node_2(1) * Self%Primeira_Derivada_Node_1(3)
    Self%H(5,3) = Self%Primeira_Derivada_Node_2(2) * Self%Primeira_Derivada_Node_1(3)
    Self%H(6,3) = Self%Primeira_Derivada_Node_2(3) * Self%Primeira_Derivada_Node_1(3) - Self%Green / (Self%L0**2.0d0)
    Self%H(1,4) = Self%Primeira_Derivada_Node_1(1) * Self%Primeira_Derivada_Node_2(1) - Self%Green / (Self%L0**2.0d0)
    Self%H(2,4) = Self%Primeira_Derivada_Node_1(2) * Self%Primeira_Derivada_Node_2(1)
    Self%H(3,4) = Self%Primeira_Derivada_Node_1(3) * Self%Primeira_Derivada_Node_2(1)
    Self%H(4,4) = Self%Primeira_Derivada_Node_2(1) * Self%Primeira_Derivada_Node_2(1) + Self%Green / (Self%L0**2.0d0)
    Self%H(5,4) = Self%Primeira_Derivada_Node_2(2) * Self%Primeira_Derivada_Node_2(1)
    Self%H(6,4) = Self%Primeira_Derivada_Node_2(3) * Self%Primeira_Derivada_Node_2(1)
    Self%H(1,5) = Self%Primeira_Derivada_Node_1(1) * Self%Primeira_Derivada_Node_2(2)
    Self%H(2,5) = Self%Primeira_Derivada_Node_1(2) * Self%Primeira_Derivada_Node_2(2) - Self%Green / (Self%L0**2.0d0)
    Self%H(3,5) = Self%Primeira_Derivada_Node_1(3) * Self%Primeira_Derivada_Node_2(2)
    Self%H(4,5) = Self%Primeira_Derivada_Node_2(1) * Self%Primeira_Derivada_Node_2(2)
    Self%H(5,5) = Self%Primeira_Derivada_Node_2(2) * Self%Primeira_Derivada_Node_2(2) + Self%Green / (Self%L0**2.0d0)
    Self%H(6,5) = Self%Primeira_Derivada_Node_2(3) * Self%Primeira_Derivada_Node_2(2)
    Self%H(1,6) = Self%Primeira_Derivada_Node_1(1) * Self%Primeira_Derivada_Node_2(3)
    Self%H(2,6) = Self%Primeira_Derivada_Node_1(2) * Self%Primeira_Derivada_Node_2(3)
    Self%H(3,6) = Self%Primeira_Derivada_Node_1(3) * Self%Primeira_Derivada_Node_2(3) - Self%Green / (Self%L0**2.0d0)
    Self%H(4,6) = Self%Primeira_Derivada_Node_2(1) * Self%Primeira_Derivada_Node_2(3)
    Self%H(5,6) = Self%Primeira_Derivada_Node_2(2) * Self%Primeira_Derivada_Node_2(3)
    Self%H(6,6) = Self%Primeira_Derivada_Node_2(3) * Self%Primeira_Derivada_Node_2(3) + Self%Green / (Self%L0**2.0d0)
    
    Young = GET_Material_Young(Self%Material)
    Self%H = Self%H * GET_SECAO_AREA(Self%Secao) * Self%L0 * Young
    END SUBROUTINE SET_HESSIAN_ELEMENT
!----------------------------------------------------------------------------------------------------------    
    SUBROUTINE SET_M_MASSA_ELEMENT_LUMPED(Self)
    IMPLICIT NONE
    TYPE(Type_Element),  INTENT(INOUT)      ::  Self
    REAL*8                                  ::  Densidade
    REAL*8                                  ::  Area
    REAL*8                                  ::  L0
    REAL*8                                  ::  Constante
    INTEGER*4                               ::   j
    
    
    Densidade = GET_Material_Densidade(Self%Material)
    Area = GET_SECAO_AREA(Self%Secao)
    L0 = Self%L0
    Constante = Densidade*Area*L0
    Self%M = 0.0d0
      
    DO j = 1,6
        Self%M(J,J) = (Constante/2.0d0)
    END DO
    
    END SUBROUTINE SET_M_MASSA_ELEMENT_LUMPED
!----------------------------------------------------------------------------------------------------------    
    SUBROUTINE SET_M_MASSA_ELEMENT(Self)
    IMPLICIT NONE
    TYPE(Type_Element),  INTENT(INOUT)      ::  Self
    REAL*8                                  ::  Densidade
    REAL*8                                  ::  Area
    REAL*8                                  ::  L0
    REAL*8                                  ::  Constante
    INTEGER*4                               ::   j
    
    
    Densidade = GET_Material_Densidade(Self%Material)
    Area = GET_SECAO_AREA(Self%Secao)
    L0 = Self%L0
    Constante = Densidade*Area*L0
    Self%M = 0.0d0
      
    DO j = 1,6
        Self%M(J,J) = (Constante/2.0d0)
    END DO
    
    END SUBROUTINE SET_M_MASSA_ELEMENT
!----------------------------------------------------------------------------------------------------------    
    SUBROUTINE SET_Matriz_ten_def_ELEMENT(Self)
    IMPLICIT NONE
    TYPE(Type_Element),  INTENT(INOUT)      ::  Self
    REAL*8                                  ::  Young
    REAL*8                                  ::  Alpha
    REAL*8                                  ::  Beta
    REAL*8                                  ::  Gama
    REAL*8                                  ::  Area
    
    
    Alpha = Self%Alpha
    Beta = Self%Beta
    Gama = Self%Gama
    Young = GET_Material_Young(Self%Material)
    Area = GET_secao_Area(Self%Secao)
    
    IF (allocated(Self%D) .eqv. .false.) THEN
        ALLOCATE(Self%D(6,6))
    ELSE
        DEALLOCATE(Self%D)
        ALLOCATE(Self%D(6,6))
    ENDIF
    
    
    Self%D(1,1) = Alpha**4
    Self%D(1,2) = (Alpha**2)*(Beta**2)
    Self%D(1,3) = (Alpha**2)*(Gama**2)
    Self%D(1,4) = (Alpha**2)*(Beta*Gama)
    Self%D(1,5) = (Alpha**3)*Gama
    Self%D(1,6) = (Alpha**3)*Beta
    Self%D(2,2) = Beta**4
    Self%D(2,3) = (Beta**2)*(Gama**2)
    Self%D(3,3) = Gama**4
    Self%D(2,4) = (Beta**3)*Gama
    Self%D(3,4) = (Gama**3)*Beta
    Self%D(4,4) = (Beta**2)*(Gama**2)
    Self%D(2,5) = (Beta**2)*Alpha*Gama
    Self%D(3,5) = (Gama**3)*Alpha
    Self%D(4,5) = Alpha*Beta*Gama
    Self%D(5,5) = (Alpha**2)*(Gama**2)
    Self%D(2,6) = (Beta**3)*Alpha
    Self%D(3,6) = (Gama**2)*Alpha*Beta
    Self%D(4,6) = (Beta**2)*Alpha*Gama
    Self%D(5,6) = (Alpha**2)*Beta*Gama
    Self%D(6,6) = (Alpha**2)*(Beta**2)
    
    Self%D(2,1) = Self%D(1,2)
    Self%D(3,1) = Self%D(1,3)
    Self%D(4,1) = Self%D(1,4)
    Self%D(5,1) = Self%D(1,5)
    Self%D(6,1) = Self%D(1,6)
    Self%D(3,2) = Self%D(2,3)
    Self%D(4,2) = Self%D(2,4)
    Self%D(5,2) = Self%D(2,5)
    Self%D(6,2) = Self%D(2,6)
    Self%D(4,3) = Self%D(3,4)
    Self%D(5,3) = Self%D(3,5)
    Self%D(6,3) = Self%D(3,6)
    Self%D(5,4) = Self%D(4,5)
    Self%D(6,4) = Self%D(4,6)
    Self%D(6,5) = Self%D(5,6)
    
    Self%D = Self%D * Young * Self%L0 * Area
    END SUBROUTINE SET_Matriz_ten_def_ELEMENT
!----------------------------------------------------------------------------------------------------------    
    SUBROUTINE SET_Matriz_def_des_ELEMENT(Self)
    IMPLICIT NONE
    TYPE(Type_Element),  INTENT(INOUT)      ::  Self
    REAL*8                                  ::  Young
    REAL*8                                  ::  Alpha
    REAL*8                                  ::  Beta
    REAL*8                                  ::  Gama
    
    Alpha = Self%Alpha
    Beta = Self%Beta
    Gama = Self%Gama
    
    Self%B = 0.0d0
    Self%B(1,1) = -Alpha
    Self%B(1,4) = Alpha
    
    Self%B(2,2) = -Beta
    Self%B(2,5) = Beta
    
    Self%B(3,3) = -Gama
    Self%B(3,6) = Gama
    
    Self%B(4,1) = -Beta
    Self%B(4,2) = -Alpha
    Self%B(4,4) = Beta
    Self%B(4,5) = Alpha
    
    Self%B(5,2) = -Gama
    Self%B(5,3) = -Beta
    Self%B(5,5) = Gama
    Self%B(5,6) = Beta
    
    Self%B(6,1) = -Gama
    Self%B(6,3) = -Alpha
    Self%B(6,4) = Gama
    Self%B(6,6) = Alpha
    
    
    Self%B = Self%B / Self%L0
    END SUBROUTINE SET_Matriz_def_des_ELEMENT
    
!----------------------------------------------------------------------------------------------------------    
    SUBROUTINE SET_Matriz_des_esp_Element(Self)
    USE MA67_INTERFACE
    IMPLICIT NONE
    TYPE(Type_Element),  INTENT(INOUT)      ::  Self
    INTEGER*4                               ::  i
    INTEGER*4                               ::  j
    INTEGER*4                               ::  count
    INTEGER*4,  ALLOCATABLE                 ::  i_sparce(:), j_sparce(:)
    REAL*8, ALLOCATABLE                     ::  H_sparce(:)
    REAL*8, ALLOCATABLE                     ::  g_Vetor(:)
    REAL*8,  ALLOCATABLE                    ::  Hessiana(:,:)
    
    ALLOCATE (g_Vetor(6), Hessiana(6,6))
    
    
    Self%X = 0.0d0
    
    
    
    Hessiana = Self%H
    count = 0
    DO i = 1, 6
        DO j = i, 6
            IF ( abs(Hessiana(i, j)) > 1.0d-20)   count = count + 1
        END DO
    END DO
    IF (allocated(i_sparce) .eqv. .false.) THEN
        ALLOCATE(i_sparce(count),j_sparce(count),H_sparce(count))
    ELSE
        DEALLOCATE(i_sparce,j_sparce,H_sparce)
        ALLOCATE(i_sparce(count),j_sparce(count),H_sparce(count))
    ENDIF
    count = 0
    DO i = 1, 6
        DO j = i, 6
            IF ( abs(Hessiana(i, j)) > 1.0d-20) THEN
                count = count + 1
                i_sparce(count) = i
                j_sparce(count) = j
                H_sparce(count) = Hessiana(i,j)
            END IF
        END DO
    END DO
    DO i=1, 6
        g_Vetor  = Self%M(i,:)
        CALL MA67SOLVER(count,6,i_sparce, j_sparce, H_sparce, g_Vetor, Self%X(i,:)) 
    END DO
        
        
        
    END SUBROUTINE SET_Matriz_des_esp_Element
!----------------------------------------------------------------------------------------------------------    
    SUBROUTINE SET_Matriz_rig_homo_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Element),  INTENT(INOUT)      ::  Self
    REAL*8                                  ::  Area
    REAL*8                                  ::  L0
    REAL*8                                  ::  Constante
    INTEGER*4                               ::  i
    INTEGER*4                               ::  j
    REAL*8                                  ::  Btrans(6,6)
    
    Area = GET_SECAO_AREA(Self%Secao)
    L0 = Self%L0
    Constante = Area*L0
    DO i=1,6
        DO j=1,6
            Btrans(i,j)=Self%B(j,i)
        END DO
    END DO
    Self%M = MATMUL(Btrans,Self%D)
    Self%H = MATMUL(Self%M,Self%B)
    
    Self%M = Self%M*Constante
    Self%H = Self%H*Constante
    END SUBROUTINE SET_Matriz_rig_homo_Element

END MODULE Element