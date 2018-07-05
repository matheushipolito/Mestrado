MODULE Nodes

IMPLICIT NONE
 
!==========================================================================================================
TYPE Type_Nodes
    Private
    Integer*4                            ::   ID_Node                   !Número do nó
    Real*8                               ::   x0 = 0.0d0                !Coordenada do nó em X
    Real*8                               ::   y0 = 0.0d0                !Coordenada do nó em Y
    Real*8                               ::   z0 = 0.0d0                !Coordenada do nó em Z  
    Real*8                               ::   xf = 0.0d0                !Coordenada do nó em X
    Real*8                               ::   yf = 0.0d0                !Coordenada do nó em Y
    Real*8                               ::   zf = 0.0d0                !Coordenada do nó em Z 
    INTEGER*4                            ::   Element_por_node   
    INTEGER*4                            ::   Num_DOF                   !Número de graus de liberdade

   
    INTEGER*4,   DIMENSION(3)            ::   Node_DOF                  !Graus de liberdade global   
    
END TYPE Type_Nodes

    CONTAINS
    !=====================================================================================================
    !Ler dados dos nós no arquivo
    SUBROUTINE Ler_Nodes(Self,unidade)                             !Self é o Type_Nodes
    IMPLICIT NONE
    TYPE(Type_Nodes), INTENT(INOUT)      ::   Self    
    Integer*4                            ::   ID_Node      !Número do nó
    Real*8                               ::   x0, y0, z0   !Coordenadas
    INTEGER*8                            ::  unidade
    Read(unidade,*)ID_Node, x0, y0, z0                           !Leitura dos dados  
    CALL SET_Nodes(Self, ID_Node, x0, y0, z0)
    END SUBROUTINE Ler_Nodes
    !=====================================================================================================
    !Armazenar dados dos nós Self
    SUBROUTINE SET_Nodes(Self, ID_Node, x0, y0, z0)        !Armazenar ID e coordenadas no Self
    TYPE(Type_nodes), INTENT(INOUT)      ::  Self          !Tipo de objeto: nó
    INTEGER*4, INTENT(IN)                ::  ID_Node       !ID do nó
    REAL*8, INTENT(IN)                   ::  x0, y0, z0    !Coordenadas do nó na posição inicial
    Self%ID_Node=ID_Node
    Self%x0=x0
    Self%y0=y0
    Self%z0=z0
    Self%xf=x0
    Self%yf=y0
    Self%zf=z0
    END SUBROUTINE SET_Nodes
    !=================================================================================================
    SUBROUTINE  SET_Element_por_node(Self, Element_por_node)! Armacenar o número de elementos que chegan ao nó
    IMPLICIT NONE
    TYPE(Type_nodes), INTENT(INOUT)   ::  Self
    INTEGER*4                         ::  Element_por_node
    Self%Element_por_node = Element_por_node
    END SUBROUTINE SET_Element_por_node
    !==================================================================================================
    SUBROUTINE  SET_Num_DOF(Self, Num_DOF)
    IMPLICIT NONE
    TYPE(Type_nodes), INTENT(INOUT)   ::  Self
    INTEGER*4                         ::  Num_DOF          ! Graus de liberdade por nó
    Self%Num_DOF = Num_DOF
    END SUBROUTINE SET_Num_DOF
    !=================================================================================================
    SUBROUTINE  SET_Node_DOF(Self, DOF)
    IMPLICIT NONE
    TYPE(Type_nodes), INTENT(INOUT)          :: Self
    INTEGER*4       , DIMENSION(3)           :: Node_DOF
    INTEGER*4       , INTENT(INOUT)          :: DOF
    INTEGER*4                                :: i
    DO i=1, Self%Num_DOF
        DOF = DOF + 1
        Self%Node_DOF(i) = DOF
    END DO
    END SUBROUTINE SET_Node_DOF
    !=================================================================================================
    FUNCTION GET_Node_DOF(Self, Local_DOF)            RESULT(Node_DOF)
    IMPLICIT NONE
    TYPE(Type_nodes), INTENT(IN)   ::  Self
    INTEGER*4                      ::  Local_DOF
    INTEGER*4                      ::  Node_DOF
    Node_DOF = Self%Node_DOF(Local_DOF)
    END FUNCTION GET_Node_DOF
    ! =================================================================================================
    FUNCTION GET_x0(Self)       RESULT(x0)                                      
    IMPLICIT NONE
    TYPE(Type_nodes),          INTENT(IN)     ::     Self
    REAL*8   ::  x0
    x0 = Self%x0
    END FUNCTION GET_x0
    !========================================================================================
    FUNCTION GET_y0(Self)       RESULT(y0)                                      
    IMPLICIT NONE
    TYPE(Type_nodes),          INTENT(IN)     ::     Self
    REAL*8   ::  y0
    y0 = Self%y0
    END FUNCTION GET_y0
    !========================================================================================
    FUNCTION GET_z0(Self)       RESULT(z0)                                      
    IMPLICIT NONE
    TYPE(Type_nodes),          INTENT(IN)     ::     Self
    REAL*8   ::  z0
    z0 = Self%z0
    END FUNCTION GET_z0
    !========================================================================================
    FUNCTION GET_xf(Self)       RESULT(xf)                                      
    IMPLICIT NONE
    TYPE(Type_nodes),          INTENT(IN)     ::     Self
    REAL*8   ::  xf
    xf = Self%xf
    END FUNCTION GET_xf
    !========================================================================================
    FUNCTION GET_yf(Self)       RESULT(yf)                                      
    IMPLICIT NONE
    TYPE(Type_nodes),          INTENT(IN)     ::     Self
    REAL*8   ::  yf
    yf = Self%yf
    END FUNCTION GET_yf
    !========================================================================================
    FUNCTION GET_zf(Self)       RESULT(zf)                                      
    IMPLICIT NONE
    TYPE(Type_nodes),          INTENT(IN)       ::     Self
    REAL*8   ::  zf
    zf = Self%zf
    END FUNCTION GET_zf
    !========================================================================================
    SUBROUTINE SET_Node_pos_f(Self, xf, yf, zf) 
    TYPE(Type_nodes)                           ::  Self
    REAL*8                                     ::  xf, yf, zf
    Self%xf=xf
    Self%yf=yf
    Self%zf=zf
    END SUBROUTINE SET_Node_pos_f
    
END MODULE Nodes