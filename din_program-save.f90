MODULE Din_Program
USE Input_File
USE Nodes
USE Material
USE Secao
USE Element
USE Forcas
USE Cond_Contorno
USE Leitura
USE Output_File
IMPLICIT NONE
PUBLIC

!=============================================================================================================================
TYPE Type_Problem
    PRIVATE
    
    CHARACTER*50                                ::  Tipo_Estrutura
    
    INTEGER*4                                   ::  N_Nodes
    INTEGER*4                                   ::  N_Materiais
    INTEGER*4                                   ::  N_Secao
    INTEGER*4                                   ::  N_Element
    INTEGER*4                                   ::  N_Forcas
    INTEGER*4                                   ::  N_Cond_Contorno
    INTEGER*4                                   ::  N_DOF
    REAL*8                                      ::  Norm_dpos                   ! Norma R2 do delta_pos
    
    
    TYPE(Type_Nodes)         ,  ALLOCATABLE     ::  Node(:)                  
    TYPE(Type_Material)      ,  ALLOCATABLE     ::  Material(:)
    TYPE(Type_Secao)         ,  ALLOCATABLE     ::  Secao(:)
    TYPE(Type_Forcas)        ,  ALLOCATABLE     ::  Forcas(:) 
    TYPE(Type_Cond_Contorno) ,  ALLOCATABLE     ::  Cond_Contorno(:) 
    TYPE(Type_Element)       ,  ALLOCATABLE     ::  Element(:)
    
    INTEGER*4                                   ::  Nodes_by_Element            ! Número de nós por Element. Default: 2
    REAL*8                                      ::  Dt                          ! Incremento de Tempo
    REAL*8                                      ::  N_Dt                        ! Número de passos
    REAL*8                                      ::  Relogio                     ! Tempo atual
    REAL*8                                      ::  Beta                        ! Beta do Newton Raphson
    REAL*8                                      ::  Gama
    REAL*8, POINTER,    DIMENSION(:)            ::  Finer(:)
    
    REAL*8, POINTER,    DIMENSION(:)            ::  Fi_Global                   ! Vetor das forças Internas da Estrutura
    REAL*8, POINTER,    DIMENSION(:)            ::  Fext_Global                 ! Vetor das forças Externas da Estrutura
    REAL*8, POINTER,    DIMENSION(:,:)          ::  M                           ! Matriz de Massa da Estrutura
    REAL*8, POINTER,    DIMENSION(:,:)          ::  H                           ! Matriz Hessiana da Estrutura
    REAL*8, POINTER,    DIMENSION(:)            ::  F                           ! Vetor desbalanceado
    REAL*8, POINTER,    DIMENSION(:)            ::  v_0                          ! Velocidade inicial
    REAL*8, POINTER,    DIMENSION(:)            ::  ac_0                         ! initial acceleration
    REAL*8, POINTER,    DIMENSION(:)            ::  pos_f_p                      ! predictor position
    REAL*8, POINTER,    DIMENSION(:)            ::  v_p                         ! predictor velocity
    REAL*8, POINTER,    DIMENSION(:)            ::  ac_p                        ! predictor acceleration
    REAL*8, POINTER,    DIMENSION(:)            ::  pos_f_c                      ! corrector position
    REAL*8, POINTER,    DIMENSION(:)            ::  v_c                         ! corrector velocity
    REAL*8, POINTER,    DIMENSION(:)            ::  ac_c                        ! corrector acceleration
      
    REAL*8, POINTER,  DIMENSION(:)              ::  pos_0                        ! Posição inicial
    REAL*8, ALLOCATABLE                         ::  pos_f(:)                     ! Posição atual
    REAL*8, ALLOCATABLE                         ::  Delta_pos(:)                ! Delta pos: pos_f = pos_pred + Delta_pos
    
    CHARACTER*50                                ::  Nome_Arquivo                !Nome do arquivo de saída
     
    
END TYPE Type_problem

    CONTAINS
    ! ==============================================================================================
    ! ----------------------------------- ESTÁTICO LINEAR ------------------------------------------
    ! ==============================================================================================
    
    SUBROUTINE Trelica_Est_Lin(Self)
    IMPLICIT NONE
    TYPE (Type_Problem),     INTENT(INOUT)     :: Self

    CALL Preprocessing(Self)
    CALL CALC_L0_Element(Self)                               !Tamanho inicial do elemento
    CALL CALC_COSSENOS_GIR(Self)                             !Cossenos Giratórios de cada elemento
    CALL CALC_RIGIDEZ_ELEMENT(Self)                          !Matriz Hessiana de cada elemento
    CALL CALC_HESSIAN_GLOBAL(Self)                           !Matriz Hessiana da estrutura
    CALL CALC_Fext(Self)                                     !Calcular as forças externas
    CALL Aplicar_Cond_Cont_Hessian(Self)                     !Aplicar condições de contorno na matriz Hessiana da estrutura
    Self%F = Self%Fext_Global
    CALL CALC_Delta_pos_Solver(Self)                         !Cálcula os deslocamentos
    CALL Write_Step(Self,(0))
    CALL Atualizar_Lf(Self)                                  !Atualiza a posição final dos nós
    CALL Atualizar_pos_f(Self)
    CALL Atualizar_Node_pos_f(Self)
    
    CALL Write_Step(Self,(1))
    
    END SUBROUTINE Trelica_Est_Lin

    ! ==============================================================================================
    ! -------------------------- NEWTON RAPHSON ESTÁTICO NÂO LINEAR --------------------------------
    ! ==============================================================================================  
    SUBROUTINE Trelica_Est_Nao_Lin(Self)
    IMPLICIT NONE
    TYPE (Type_Problem),     INTENT(INOUT)      ::  Self
    INTEGER*4                                   ::  i
    INTEGER*4                                   ::  j
    INTEGER*4                                   ::  Max_Iter
    INTEGER*4                                   ::  Stop_Iteration
    INTEGER*4                                   ::  Stop_Iteration2
    INTEGER*4                                   ::  Stop_Iteration3
    REAL*8                                      ::  Tolerance
    REAL*8                                      ::  Var

    
    Tolerance = 1.0d-8
    Self%Relogio = 0.0d0
    MAX_Iter = 150
    CALL Preprocessing_Nao_Lin(Self)
    CALL CALC_L0_Element(Self)                                      !Tamanho inicial do elemento
    CALL Write_Step(Self,(0))
    DO i=1,Self%N_Dt
        
        Self%Relogio = Self%Relogio + Self%Dt
        Stop_Iteration = 0
        Stop_Iteration2 = 0
        Stop_Iteration3 = 0
        Var = 1.0d0
        CALL CALC_Fext(Self) 
        
        
        DO WHILE (( Var >= Tolerance) .AND. (Stop_Iteration<= MAX_Iter))
           
            Stop_Iteration = Stop_Iteration + 1
            Stop_Iteration2 = Stop_Iteration2 + 1
            Stop_Iteration3 = Stop_Iteration3 + 1
            CALL Atualizar_Lf(Self)
            CALL Atualizar_Green(Self)
            CALL Atualizar_Primeira_Derivada(Self)                          !Calcular primeira derivada de E
            CALL CALC_HESSIAN_ELEMENT(Self)                                 !Calcular a Hessiana de cada elemento
            CALL CALC_HESSIAN_GLOBAL(Self)                                  !Calcular a Hessiana da estrutura

            CALL CALC_Fi_Element_N_Linear(Self)
            CALL CALC_Fi_Global(Self)
            CALL CALC_g_desbalanceamento(Self)
            CALL Aplicar_Cond_Cont_Hessian(Self)                            !Aplicar condições de contorno na matriz Hessiana da estrutura
            CALL Aplicar_Cond_Cont_desbal(Self)
            CALL CALC_Delta_pos_Solver(Self)
            CALL Atualizar_pos_f(Self)
            CALL Atualizar_Node_pos_f(Self)
            
            
            
            IF (Stop_Iteration2.GE.10) THEN
                Print *,Var
                Stop_Iteration2 = 0
            END IF
            IF (Stop_Iteration.GE.2) THEN
                Var = Norm(Self,Self%Delta_pos)
            END IF
        END DO
        
    END DO
    
    CALL Write_Step(Self,(1))
    
    END SUBROUTINE Trelica_Est_Nao_Lin
    ! ==============================================================================================
    ! ------------------------ NEWTON RAPHSON DINÂMICO NÂO LINEAR ----------------------------------
    ! ==============================================================================================  
    SUBROUTINE Trelica_Din_Nao_Lin(Self)
    IMPLICIT NONE
    TYPE (Type_Problem),     INTENT(INOUT)      ::  Self
    INTEGER*4                                   ::  Max_Iter
    INTEGER*4                                   ::  i
    INTEGER*4                                   ::  Stop_Iteration
    REAL*8                                      ::  Tolerance
    REAL*8                                      ::  Var

    Tolerance = 1.0d-8
    Self%Relogio = 0.0d0
    MAX_Iter = 150
    
    
    CALL Preprocessing_Din_Nao_Lin(Self)
    CALL CALC_L0_Element(Self)                                      !Tamanho inicial do elemento
    CALL CALC_M_MASSA_ELEMENT(Self)                                 !Calcular a Matriz de Massa de cada elemento
    CALL CALC_M_MASSA_GLOBAL(Self)                                  !Calcular a Matriz de Massa da estrutura
    CALL Write_Step(Self,(0))
    
    
    DO i=1,Self%N_Dt
        
        Self%Relogio = Self%Relogio + Self%Dt
        Stop_Iteration = 0
        Var = 1.0d0
        CALL CALC_Fext(Self)  
        
        DO WHILE (( Var >= Tolerance) .AND. (Stop_Iteration<= MAX_Iter))
    
            Stop_Iteration = Stop_Iteration + 1
            CALL Atualizar_Lf(Self)
            CALL Atualizar_Green(Self)
            CALL Atualizar_Primeira_Derivada(Self)                      !Calcular primeira derivada de E
            CALL CALC_HESSIAN_ELEMENT(Self)                             !Calcular a Hessiana de cada elemento
            CALL CALC_HESSIAN_GLOBAL(Self)                              !Calcular a Hessiana da estrutura
            CALL CALC_HESSIAN_GLOBAL_Dyn(Self)                          !Calcular a Hessiana da estrutura para o sistema dinâmico
            
            CALL Aplicar_Cond_Cont_Hessian(Self)                        !Aplicar condições de contorno na matriz Hessiana da estrutura
            CALL Aplicar_Cond_Cont_Massa(Self)                          !Aplicar condições de contorno na matriz Hessiana da estrutura
            CALL CALC_Fi_Dyn_Global(Self)
            CALL CALC_g_desbalanceamento(Self)
            CALL CALC_Delta_pos_Solver(Self)                            !Cálcula os deslocamentos
            CALL Atualizar_pos_f(Self)                                  !Atualiza a posição final
            CALL Atualizar_Node_pos_f(Self)
            
            IF (Stop_Iteration.GE.2) THEN
                Var = Norm(Self,Self%Delta_pos)
            END IF
        END DO
        CALL CALC_Corretivos(Self)                                      !Calcula a velocidade e acelerações corretivas
        CALL CALC_Preditivos(Self)                                      !Calcula a velocidade e acelerações preditivas
        CALL Write_Step(Self,(i))
        Print *,i
    END DO
    
    
    END SUBROUTINE Trelica_Din_Nao_Lin
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Preprocessing(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),   INTENT(INOUT)             ::  Self
    CALL write_dates(Self)
    CALL SET_pos_0_Global(Self)
    CALL CALC_node_DOF(Self)
    CALL Inicial_pos_f_Global(Self)
    Self%Norm_Dpos = Inicial_Norm(Self)
    CALL Inicial_Fi_Global(Self)   
        
    END SUBROUTINE Preprocessing
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Preprocessing_Nao_Lin(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),   INTENT(INOUT)             ::  Self
    
    CALL write_dates(Self)
    CALL SET_pos_0_Global(Self)

    CALL CALC_node_DOF(Self)
    CALL CALC_L0_Element(Self)
    CALL Inicial_pos_f_Global(Self)
    Self%Norm_Dpos = Inicial_Norm(Self)

    CALL Inicial_Fi_Global(Self)  
    CALL ALLOCATE_g_Desbalanceamento(self)
    CALL Inicial_Primeira_Derivada(Self)
        
    END SUBROUTINE Preprocessing_Nao_Lin
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Preprocessing_Din_Nao_Lin(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),   INTENT(INOUT)             ::  Self
    INTEGER*4                                       ::  i
    CALL write_dates(Self)
    CALL SET_pos_0_Global(Self)
    !CALL Inic_pos_0(Self)
    CALL Inicial_pos_f_Global(Self)
    Self%Norm_Dpos = Inicial_Norm(Self)
    CALL Inicial_Fi_Global(Self)  
    CALL ALLOCATE_g_Desbalanceamento(self)
    CALL Inicial_Dyn(Self)
    CALL CALC_Inicial_Preditivos(Self)
    CALL Inicial_Primeira_Derivada(Self)
        
    END SUBROUTINE Preprocessing_Din_Nao_Lin
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Corretivos(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    CALL CALC_Ac_c(Self)
    CALL CALC_V_c(Self)
    END SUBROUTINE CALC_Corretivos
    !----------------------------------------------------------------------------------------------++++++++++++
    SUBROUTINE CALC_V_c(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    Self%V_c = Self%Gama * (Self%pos_f - Self%pos_0)/(Self%Beta * Self%Dt) + Self%V_p - Self%Ac_c*Self%Dt*Self%Gama  
    END SUBROUTINE CALC_V_c
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Ac_c(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    Self%Ac_c = ((Self%pos_f - Self%pos_0)/(Self%Beta * Self%Dt**2.0d0)) - Self%Ac_p
    END SUBROUTINE CALC_Ac_c
    !----------------------------------------------------------------------------------------------++++++++++++
    SUBROUTINE CALC_Preditivos(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    CALL CALC_Ac_p(Self)
    CALL CALC_V_p(Self)
    END SUBROUTINE CALC_Preditivos
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_V_p(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    Self%V_p = Self%V_c + Self%Ac_c * Self%Dt * (1.0d0 - Self%Gama)
    END SUBROUTINE CALC_V_p
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Ac_p(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    Self%Ac_p = (Self%pos_f - Self%pos_0)/(Self%Beta * Self%Dt**2.0d0) + Self%V_c/(Self%Beta * Self%Dt) + Self%Ac_p * (( 1.0d0 / 2.0d0*Self%Beta) - 1)
    END SUBROUTINE CALC_Ac_p
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Inicial_Preditivos(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    CALL CALC_Inicial_Ac_p(Self)
    CALL CALC_Inicial_V_p(Self)
    END SUBROUTINE CALC_Inicial_Preditivos
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Inicial_Ac_p(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    Self%Ac_p = (Self%pos_f - Self%pos_0)/(Self%Beta * Self%Dt**2.0d0) + Self%V_0/(Self%Beta * Self%Dt) + Self%Ac_0 * (( 1.0d0 / 2.0d0*Self%Beta) - 1)
    END SUBROUTINE CALC_Inicial_Ac_p
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Inicial_V_p(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    Self%V_p = Self%V_0 + Self%Ac_0 * Self%Dt * (1.0d0 - Self%Gama)
    END SUBROUTINE CALC_Inicial_V_p
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Inicial_Dyn(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    ALLOCATE (Self%V_0(3*Self%N_Nodes))
    ALLOCATE (Self%Ac_0(3*Self%N_Nodes))
    ALLOCATE (Self%pos_f_p(3*Self%N_Nodes))
    ALLOCATE (Self%V_p(3*Self%N_Nodes))
    ALLOCATE (Self%Ac_p(3*Self%N_Nodes))
    ALLOCATE (Self%pos_f_c(3*Self%N_Nodes))
    ALLOCATE (Self%V_c(3*Self%N_Nodes))
    ALLOCATE (Self%Ac_c(3*Self%N_Nodes))
    ALLOCATE (Self%Finer(3*Self%N_Nodes))
    Self%V_0 = 0.0d0
    Self%Ac_0 = 0.0d0
    Self%pos_f_C = 0.0d0
    Self%V_p = 0.0d0
    Self%Ac_p = 0.0d0
    Self%pos_f_c = 0.0d0
    Self%V_c = 0.0d0
    Self%Ac_c = 0.0d0
    Self%Beta = 0.5d0
    Self%Gama = 2.5d-1
    Self%Finer = 0.0d0
    
    END SUBROUTINE Inicial_Dyn
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_HESSIAN_GLOBAL_Dyn(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    
    Self%H = Self%H + Self%M / ( Self%Beta * Self%Dt**2.0d0 )
    
    END SUBROUTINE CALC_HESSIAN_GLOBAL_Dyn
    !----------------------------------------------------------------------------------------------   
    SUBROUTINE CALC_Fi_Element(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*4                              :: i
    REAL*8,                 DIMENSION(6)   :: deslocamentos
    INTEGER*4                              :: I1,I2
    DO i=1,Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2)
        deslocamentos(1) = Get_xf(Self%Node(I1)) - Get_x0(Self%Node(I1))
        deslocamentos(2) = Get_yf(Self%Node(I1)) - Get_y0(Self%Node(I1))
        deslocamentos(3) = Get_zf(Self%Node(I1)) - Get_z0(Self%Node(I1))
        deslocamentos(4) = Get_xf(Self%Node(I2)) - Get_x0(Self%Node(I2))
        deslocamentos(5) = Get_yf(Self%Node(I2)) - Get_y0(Self%Node(I2))
        deslocamentos(6) = Get_zf(Self%Node(I2)) - Get_z0(Self%Node(I2))
        CALL SET_Element_Fi(Self%Element(i),deslocamentos)
    END DO
    END SUBROUTINE CALC_Fi_Element
    !----------------------------------------------------------------------------------------------   
    SUBROUTINE CALC_Fi_Element_N_Linear(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*4                              :: i
    DO i=1,Self%N_Element
        CALL SET_Element_Fi_N_Linear(Self%Element(i))
    END DO
    END SUBROUTINE CALC_Fi_Element_N_Linear
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_M_MASSA_ELEMENT(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*4                              :: i
    DO i=1,Self%N_Element
        CALL SET_M_MASSA_Element(Self%Element(i))   
    END DO
    END SUBROUTINE CALC_M_MASSA_ELEMENT
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_HESSIAN_ELEMENT(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*4                              :: i
    DO i=1,Self%N_Element
        CALL SET_HESSIAN_Element(Self%Element(i))   
    END DO
    END SUBROUTINE CALC_HESSIAN_ELEMENT
    !----------------------------------------------------------------------------------------------    
    SUBROUTINE CALC_g_desbalanceamento(SELF)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)    ::     SELF
    INTEGER*4  :: i
    REAL*8, DIMENSION(3*self%N_Nodes)     :: Fint
    REAL*8, DIMENSION(3*self%N_Nodes)     :: Fext
    Self%F = 0.0d0
    DO i = 1 ,  3*self%N_Nodes
        Fint = self%Fi_Global
        Fext = self%Fext_Global
        Self%F(i) = Fint(i) - Fext(i)
    END DO
    END SUBROUTINE CALC_g_desbalanceamento
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Read_File(Self, Nome_Arquivo)
    IMPLICIT NONE
    TYPE(Type_Problem) , INTENT(INOUT)     :: Self                       !Representa o tipo de problema estrutural (treliça)
    CHARACTER*50       , INTENT(IN)        :: Nome_Arquivo               !INPUT File_Name
    CHARACTER*50                           :: buscar_palavra             !Palavra Auxiliar pra buscar
    INTEGER*4                              :: i                          !Índice para loops
    INTEGER*4          , PARAMETER         :: Nodes_por_Element = 2      ! Número de nós por Element
    
    CALL Open_File(Nome_Arquivo, 1)                                      !Subroutine: "Input_File"

    buscar_palavra = 'Tipo_Estrutura'                                    !Busca "Tipo_Estrutura" em 'Input_File'
    CALL Buscar(buscar_palavra)                                          !Subroutine: "Leitura"
    READ(1,*) Self%Tipo_Estrutura                                        !Lê as linhas embaixo de buscar palavra
    
    buscar_palavra = 'Coordenadas'
    Self%N_Nodes = counter(buscar_palavra)                               !Conta o número de nós e joga no N_Nodes dentro do Type_Problem
    ALLOCATE( Self%Node( Self%N_Nodes ) )                                !Aloca o Type_Nodes
    CALL Buscar(buscar_palavra)                                          !Volta o cursor para a primeira posição do nó
    DO i = 1, Self%N_Nodes
        CALL Ler_Nodes (Self%Node (i))                                   !Leitura das coordenadas iniciais dos nós
    END DO
    
    buscar_palavra = 'Material'
    Self%N_Materiais = counter(buscar_palavra)                           !Conta o número de materiais e joga no N_Materiais dentro do Type_Problem
    ALLOCATE( Self%Material( Self%N_Materiais ) )                        !Aloca o Type_Material
    CALL Buscar(buscar_palavra)                                          !Volta o cursor para a primeira posição do material
    DO i = 1, Self%N_Materiais
        CALL Ler_Material (Self%Material (i))                            !Leitura dos dados dos materiais
    END DO
    
    buscar_palavra = 'Secao'
    Self%N_Secao = counter(buscar_palavra)                          
    ALLOCATE( Self%Secao( Self%N_Secao ) )                       
    CALL Buscar(buscar_palavra)                                      
    DO i = 1, Self%N_Secao
        CALL Ler_Secao (Self%Secao (i))                          
    END DO
    
    
    buscar_palavra = 'Atribuir_Propriedades'
    Self%N_Element = counter(buscar_palavra)
    ALLOCATE( Self%Element( Self%N_Element ) )
    CALL Buscar(buscar_palavra)
    DO i = 1, Self%N_Element
        CALL Ler_Element(Self%Element(i))
        CALL SET_Nodes_by_Element(Self%Element(i), Nodes_por_Element)
    END DO
    CALL CALC_Material_do_Element(Self)
    CALL CALC_Secao_do_Element(Self)
    
    
    buscar_palavra = 'Incidencia'
    CALL Buscar(buscar_palavra)
    DO i = 1, Self%N_Element
        CALL Ler_Incidencia( Self%Element(i) )                               ! Leitura da incidencia por elemento, nos locais por elemento
    END DO
    

    buscar_palavra = 'Forcas'
    Self%N_Forcas = counter(buscar_palavra)                          
    ALLOCATE( Self%Forcas( Self%N_Forcas ) )                       
    CALL Buscar(buscar_palavra)                                      
    DO i = 1, Self%N_Forcas
        CALL Ler_Forcas (Self%Forcas (i))                          
    END DO
    
    
    buscar_palavra = 'Cond_Contorno'
    Self%N_Cond_Contorno = counter(buscar_palavra)                          
    ALLOCATE( Self%Cond_Contorno( Self%N_Cond_Contorno ) )                       
    CALL Buscar(buscar_palavra)                                      
    DO i = 1, Self%N_Cond_Contorno
        CALL Ler_Cond_Contorno (Self%Cond_Contorno (i))                          
    END DO
    
    buscar_palavra = 'Tempo'
    Self%N_Dt = counter(buscar_palavra)                                            
    CALL Buscar(buscar_palavra)                                      
    CALL Ler_Tempo(Self)                         
    
    
    
    CALL CALC_Element_por_node(Self)                                    ! Número de elementos que chegan ao nó
    CALL CALC_Num_DOF(Self)
    
    CALL CALC_node_dof(Self)
    DO i=1, Self%N_Element
        CALL ALLOCATE_Nodes_of_Element(Self%Element(i))
        CALL ALLOCATE_M_Element(Self%Element(i))
        CALL ALLOCATE_H_Element(Self%Element(i))
    END DO
    CALL CALC_Nodes_of_Element(Self)

    CALL CALC_N_DOF(Self)
    
    !inicializar outras variaveis
    CALL ALLOCATE_Fext_Global(Self)
    CALL ALLOCATE_Fi_Global(Self)
    CALL ALLOCATE_g_desbalanceamento(Self)
    CALL ALLOCATE_Global_Hessian(Self)
    CALL ALLOCATE_M_Global(Self)
    CALL ALLOCATE_delta_pos_Global(Self)

    
    


    
    
    END SUBROUTINE Read_File
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Ler_Tempo(Self)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)    ::   Self    
    Real*8                               ::   Dt, N_Dt  
    Read(1,*)Dt, N_Dt                                      !Leitura dos dados
    Call Set_Tempo(Self, Dt, N_Dt)
    END SUBROUTINE Ler_Tempo
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Set_Tempo(Self, Dt, N_Dt)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    REAL*8,             INTENT(INOUT)   ::  Dt  
    REAL*8,             INTENT(INOUT)   ::  N_Dt 
    Self%Dt=Dt
    Self%N_Dt=N_Dt
    END SUBROUTINE Set_Tempo
    !----------------------------------------------------------------------------------------------    
    SUBROUTINE Atualizar_Lf(Self)
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*4                              :: i
    DO i=1,Self%N_Element
        CALL SET_Lf_Element(Self%Element(i))
    END DO
    END SUBROUTINE Atualizar_Lf
    !----------------------------------------------------------------------------------------------    
    SUBROUTINE Atualizar_Green(Self)
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*4                              :: i
    DO i=1,Self%N_Element
        CALL SET_Green_Element(Self%Element(i))
    END DO
    END SUBROUTINE Atualizar_Green
    !---------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Fext(Self)
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*4                              :: i
    INTEGER*4                              :: j
    Self%Fext_Global = 0.0d0
    DO i=1,Self%N_Nodes
        DO j=1,3
            Self%Fext_Global((3*(i-1))+j)= Self%Fext_Global((3*(i-1))+j) + GET_Fext(Self%Forcas(i),j)
        END DO
    END DO
    END SUBROUTINE CALC_Fext
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Aplicar_Cond_Cont_Massa(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)     :: Self
    INTEGER*4                              :: i
    INTEGER*4                              :: j
    INTEGER*4                              :: k
    INTEGER*4                              :: cc

    DO i=1,Self%N_Nodes
        DO j=1,3
            cc=GET_Cond_Cont(Self%Cond_Contorno(i),j)
            IF (cc==1) THEN
                DO k=1,Self%N_DOF
                    Self%M(k,((3*(i-1))+j))=0.0d0
                    Self%M(((3*(i-1))+j),k)=0.0d0
                    Self%M(((3*(i-1))+j),((3*(i-1))+j))=1.0d0
                ENDDO
            ENDIF
        ENDDO
      ENDDO    
    END SUBROUTINE Aplicar_Cond_Cont_Massa
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Aplicar_Cond_Cont_Hessian(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)     :: Self
    INTEGER*4                              :: i
    INTEGER*4                              :: j
    INTEGER*4                              :: k
    INTEGER*4                              :: cc

    DO i=1,Self%N_Nodes
        DO j=1,3
            cc=GET_Cond_Cont(Self%Cond_Contorno(i),j)
            IF (cc==1) THEN
                DO k=1,Self%N_DOF
                    Self%H(k,((3*(i-1))+j))=0.0d0
                    Self%H(((3*(i-1))+j),k)=0.0d0
                ENDDO
                Self%H(((3*(i-1))+j),((3*(i-1))+j))=1.0d0
            ENDIF
        ENDDO
      ENDDO    
    END SUBROUTINE Aplicar_Cond_Cont_Hessian
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Aplicar_Cond_Cont_desbal(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)     :: Self
    INTEGER*4                              :: i
    INTEGER*4                              :: j
    INTEGER*4                              :: k
    INTEGER*4                              :: cc

    DO i=1,(Self%N_Nodes)
        DO j=1,3
            cc=GET_Cond_Cont(Self%Cond_Contorno(i),j)
            IF (cc==1) THEN
                DO k=1,Self%N_DOF
                    Self%F((3*(i-1))+j)=0.0d0
                ENDDO
            ENDIF
        ENDDO
      ENDDO    
    END SUBROUTINE Aplicar_Cond_Cont_desbal
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_M_MASSA_GLOBAL(Self)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)      :: Self
    INTEGER*4                              :: i
    INTEGER*4                              :: j
    INTEGER*4                              :: k
    INTEGER*4                              :: dof1
    INTEGER*4                              :: dof2
    REAL*8,             DIMENSION(6,6)     :: M_elemento 
    Self%M=0.0d0
    DO i=1,Self%N_Element
        dof1=GET_Incidencia(Self%Element(i),1)
        dof2=GET_Incidencia(Self%Element(i),2)
        dof1=3*dof1
        dof2=3*dof2
        M_elemento=GET_M_MASSA_Element(Self%Element(i))
        DO j=0,2
            DO k=0,2
            Self%M(dof1-j,dof1-k) = Self%M(dof1-j,dof1-k) + M_elemento(3-j,3-k)
            Self%M(dof1-j,dof2-k) = Self%M(dof1-j,dof2-k) + M_elemento(3-j,6-k)
            Self%M(dof2-j,dof2-k) = Self%M(dof2-j,dof2-k) + M_elemento(6-j,6-k)
            Self%M(dof2-j,dof1-k) = Self%M(dof2-j,dof1-k) + M_elemento(6-j,3-k)  
            END DO
        END DO
    END DO   
    END SUBROUTINE CALC_M_MASSA_GLOBAL
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_HESSIAN_GLOBAL(Self)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)      :: Self
    INTEGER*4                              :: i
    INTEGER*4                              :: j
    INTEGER*4                              :: k
    INTEGER*4                              :: dof1
    INTEGER*4                              :: dof2
    REAL*8,             DIMENSION(6,6)     :: H_elemento
    Self%H=0.0d0
    DO i=1,Self%N_Element
        dof1=GET_Incidencia(Self%Element(i),1)
        dof2=GET_Incidencia(Self%Element(i),2)
        dof1=3*dof1
        dof2=3*dof2
        H_elemento=GET_RIGIDEZ_Element(Self%Element(i))
        DO j=0,2
            DO k=0,2
            Self%H(dof1-j,dof1-k) = Self%H(dof1-j,dof1-k) + H_elemento(3-j,3-k)
            Self%H(dof1-j,dof2-k) = Self%H(dof1-j,dof2-k) + H_elemento(3-j,6-k)
            Self%H(dof2-j,dof2-k) = Self%H(dof2-j,dof2-k) + H_elemento(6-j,6-k)
            Self%H(dof2-j,dof1-k) = Self%H(dof2-j,dof1-k) + H_elemento(6-j,3-k)  
            END DO
        END DO
    END DO   
    END SUBROUTINE CALC_HESSIAN_GLOBAL
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_RIGIDEZ_ELEMENT(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)     :: Self
    INTEGER*4                              :: i
    DO i=1, Self%N_Element
        CALL SET_RIGIDEZ_Element(Self%Element(i))
    END DO
    END SUBROUTINE CALC_RIGIDEZ_ELEMENT
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_N_DOF(Self)
    TYPE(Type_Problem),    INTENT(INOUT)   :: Self
    INTEGER*4                              :: N_DOF
    N_DOF=3*Self%N_Nodes
    Self%N_DOF=N_DOF
    END SUBROUTINE CALC_N_DOF
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_COSSENOS_GIR(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)     :: Self
    INTEGER*4                              ::  i
    DO i=1, Self%N_Element
        CALL SET_Alpha_Element(Self%Element(i))
        CALL SET_Beta_Element(Self%Element(i))
        CALL SET_Gama_Element(Self%Element(i))
    END DO
    END SUBROUTINE CALC_COSSENOS_GIR
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_L0_ELEMENT(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)     :: Self
    INTEGER*4                              ::  i
    DO i=1, Self%N_Element
        CALL SET_L0_Element(Self%Element(i))
    END DO
    END SUBROUTINE CALC_L0_ELEMENT
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Inicial_Primeira_Derivada(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)      ::     Self
    INTEGER*4                               ::     i
    DO i = 1,Self%N_Element
        CALL SET_Primeira_Derivada(Self%Element(i))
    END DO
    END SUBROUTINE Inicial_Primeira_Derivada
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Atualizar_Primeira_Derivada(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)      ::     Self
    INTEGER*4                               ::     i
    DO i = 1,Self%N_Element
        CALL SET_Primeira_Derivada(Self%Element(i))
    END DO
    END SUBROUTINE Atualizar_Primeira_Derivada
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Inicial_Fi_Global(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)      ::     Self
    INTEGER*4                               ::     i
    Self%Fi_Global = 0.0
    END SUBROUTINE Inicial_Fi_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Fi_Dyn_Global(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)      ::     Self
    Self%Fi_Global = MATMUL(Self%H,(Self%pos_f-Self%pos_0)) + MATMUL(Self%M,(Self%Ac_p))
    END SUBROUTINE CALC_Fi_Dyn_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Calc_Fi_Global(Self)
    TYPE(Type_Problem),     INTENT(INOUT)   :: Self
    INTEGER*4                               :: i
    INTEGER*4                               :: j
    INTEGER*4                               :: dof1
    INTEGER*4                               :: dof2
    REAL*8,     DIMENSION(6)                :: Fi_element
    REAL*8,     DIMENSION(3*Self%N_Nodes)   :: Fi_global
    Fi_Global = 0.0d0
    DO i=1,Self%N_Element
        dof1=GET_Incidencia(Self%Element(i),1)
        dof2=GET_Incidencia(Self%Element(i),2)
        dof1=3*dof1
        dof2=3*dof2
        Fi_element = GET_Fi_Element (Self%Element(i))
        DO j=0,2
            Fi_Global(dof1-j) = Fi_Global(dof1-j) + Fi_element(3-j)
            Fi_Global(dof2-j) = Fi_Global(dof2-j) + Fi_element(6-j)
        END DO
    END DO    
    Self%Fi_Global = Fi_global
    END SUBROUTINE CALC_Fi_Global
    !----------------------------------------------------------------------------------------------
    FUNCTION Inicial_Norm(Self)      RESULT(Den_Norm)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)      ::     Self
    REAL*8     ::   Den_Norm
    INTEGER*4  ::   i
    Den_Norm = 0.00d0
    DO i=1, 3*Self%N_Nodes
        Den_Norm = Den_Norm + Self%pos_0(i)**2
    END DO
    END FUNCTION Inicial_Norm
    !----------------------------------------------------------------------------------------------
    FUNCTION Norm(Self,F)                RESULT(Norm_atual)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)      ::  Self
    REAL*8                                  ::  Norm_atual
    REAL*8                                  ::  Norm_base
    REAL*8,       DIMENSION(3*Self%N_Nodes) ::  F
    INTEGER*4                               ::  i

    Norm_atual = 0.00d0
    Norm_base = 0.00d0
    DO i=1, 3*Self%N_Nodes
        Norm_base = Norm_base + Self%pos_0(i)**2
    END DO
    DO i=1, 3*Self%N_Nodes
        Norm_atual = Norm_atual + F(i)**2
    END DO
    Norm_atual = DSQRT(Norm_atual)
    
    END FUNCTION Norm
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Inicial_pos_f_Global(Self)
    IMPLICIT NONE
    TYPE(type_problem), INTENT(INOUT)    ::  Self
    INTEGER*4  :: i, j,gdof(3)
    Self%pos_f = Self%pos_0                                           
    DO i=1, Self%N_nodes
        gdof(1) = GET_Node_DOF(Self%node(i), 1)
        gdof(2) = GET_Node_DOF(Self%node(i), 2)
        gdof(3) = GET_Node_DOF(Self%node(i), 3)
        CALL SET_Node_pos_f(Self%node(i), Self%pos_f(Gdof(1)), Self%pos_f(Gdof(2)), Self%pos_f(Gdof(3))) 
    END DO
    END SUBROUTINE Inicial_pos_f_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Atualizar_node_pos_f(Self)
    IMPLICIT NONE
    TYPE(type_problem), INTENT(INOUT)    ::  Self
    INTEGER*4  :: i, j,gdof(3)                                           
    DO i=1, Self%N_nodes
        gdof(1) = GET_Node_DOF(Self%node(i), 1)
        gdof(2) = GET_Node_DOF(Self%node(i), 2)
        gdof(3) = GET_Node_DOF(Self%node(i), 3)
        CALL SET_Node_pos_f(Self%node(i), Self%pos_f(Gdof(1)), Self%pos_f(Gdof(2)), Self%pos_f(Gdof(3))) 
    END DO
    CALL CALC_Nodes_of_Element(Self)
    END SUBROUTINE Atualizar_node_pos_f
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_delta_pos_Global(Self)
    IMPLICIT NONE
    TYPE(type_problem), INTENT(INOUT)    ::  Self
    ALLOCATE(Self%Delta_pos( 3*Self%N_Nodes ))
    Self%Delta_pos = 0.0d0
    END SUBROUTINE ALLOCATE_delta_pos_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE SET_pos_0_Global(Self)
    IMPLICIT NONE
    TYPE(type_problem), INTENT(INOUT)    ::  Self
    INTEGER*4  :: i
    CALL ALLOCATE_pos_0_Global(Self)
    Self%pos_0 = 0.00d0
    DO i = 1, Self%N_Nodes
        Self%pos_0(GET_Node_DOF(Self%Node(i),1)) = GET_x0( Self%Node(i))
        Self%pos_0(GET_Node_DOF(Self%Node(i),2)) = GET_y0( Self%Node(i))
        Self%pos_0(GET_Node_DOF(Self%Node(i),3)) = GET_z0( Self%Node(i))
    END DO
    END SUBROUTINE SET_pos_0_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_pos_0_Global(Self)
    IMPLICIT NONE
    TYPE(type_problem), INTENT(INOUT)    ::  Self
    ALLOCATE(Self%pos_0(3*Self%N_Nodes))
    END SUBROUTINE ALLOCATE_pos_0_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Write_Dates(Self)
    IMPLICIT NONE
    TYPE(type_problem)  , INTENT(INOUT)             ::  Self
    CHARACTER*50                                    ::  Nome_Arquivo
    Nome_Arquivo='Viga_em_balanco.csv'
    !CALL temp_open_file(Nome_Arquivo)
    Self%Nome_Arquivo = Nome_Arquivo
    Close(3)
    END SUBROUTINE Write_Dates
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Material_do_Element(Self)
    IMPLICIT NONE
    TYPE(Type_problem)  , TARGET,   INTENT(INOUT)   ::  Self
    TYPE(Type_Material) , POINTER                   ::  Material
    INTEGER*4                                       ::  Material_ID
    INTEGER*4                                       ::  i
    DO i = 1, Self%N_Element
        Material_ID = GET_Element_Material_ID(Self%Element(i))
        Material => Self%Material(Material_ID)
        CALL SET_Material_do_Element(Self%Element(i), Material)
    END DO
    END SUBROUTINE CALC_Material_do_Element
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Secao_do_Element(Self)                         
    IMPLICIT NONE
    TYPE(Type_Problem)  , TARGET,   INTENT(INOUT)   ::  Self
    TYPE(Type_Secao)    , POINTER                   ::  Secao
    INTEGER*4                                       ::  Secao_ID
    INTEGER*4                                       ::  i
    DO i = 1, Self%N_Element
        Secao_ID = GET_Element_Secao_ID(Self%Element(i))
        Secao => Self%Secao(Secao_ID)
        CALL SET_Secao_do_Element(Self%Element(i), Secao)
    END DO
    END SUBROUTINE CALC_Secao_do_Element
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Element_por_node(Self)                                         ! Subrotina AUX para calcular os elementos por nó
    IMPLICIT NONE
    TYPE(type_problem), TARGET, INTENT(INOUT)       ::  Self
    INTEGER*4,          DIMENSION(Self%N_nodes)     ::  Element_por_node          ! Número de elementos que chegan ao nó
    INTEGER*4                                       ::  i                          ! Contador para os Element
    INTEGER*4                                       ::  j                          ! Contador, nós locais no Element
    INTEGER*4                                       ::  Global_node                ! Número de nós
    INTEGER*4                                       ::  Nodes_por_Element = 2      ! Número de nós por Element, DEFAULT =2
    Element_por_node = 0
    do1: DO i = 1, Self%N_Element
        do2: DO j = 1, Nodes_por_Element
            Global_node = GET_Incidencia( Self%Element(i), j )                    ! node id: given the element and local node in the element
            Element_por_node(Global_node) = Element_por_node( global_node ) + 1  ! add 1 in element_by_node(global_node) each time that node 'global_node' is repeated
        END DO do2
    END DO do1
    do3: DO i = 1, Self%N_nodes
        CALL SET_Element_por_node( Self%Node(i), Element_por_node(i) )            ! store in Self%node(i) the elements by node
    END DO do3
    END SUBROUTINE CALC_Element_por_node
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Num_DOF(Self)
    IMPLICIT NONE
    TYPE(type_problem), INTENT(INOUT)   ::  Self
    INTEGER*4                           ::  Num_DOF                   !Numero de graus de liberdade por Nó
    INTEGER*4                           ::  N_DOF
    INTEGER*4                           ::  i
    DO i = 1, Self%N_Nodes
        Num_DOF = 3
        N_DOF =  Num_DOF
        CALL SET_Num_DOF(Self%Node(i), N_DOF)
    END DO
    END SUBROUTINE CALC_Num_DOF
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_node_DOF(Self)                                     ! assign the numbering of each degree of freedom
    IMPLICIT NONE
    TYPE(type_problem), INTENT(INOUT)   ::  Self
    INTEGER*4                           ::  contar_DOF
    INTEGER*4                           ::  i
    contar_DOF = 0                                                     ! counter that accumulates de number of dof
    do1: DO i = 1, Self%N_Nodes
        !CALL ALLOCATE_Node_DOF(Self%Node(i))
        CALL SET_Node_DOF(Self%node(i), contar_DOF)
    END DO do1
    END SUBROUTINE CALC_node_DOF
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Nodes_of_Element(Self)                 ! Associar os "Pointers" dos nós ao Element
    IMPLICIT NONE
    TYPE(type_problem),TARGET, INTENT(INOUT)   ::  Self
    TYPE(Type_Nodes), POINTER                   ::  Node
    INTEGER*4   ::  Node_ID
    INTEGER*4   ::  Nodes_by_Element
    INTEGER*4   ::  i
    INTEGER*4   ::  j
    DO i = 1, Self%N_Element
        Nodes_by_Element = GET_Nodes_by_Element(Self%Element(i))
        !CALL ALLOCATE_Nodes_of_Element(Self%Element(i))
        DO j = 1, Nodes_by_Element
            Node_ID = GET_incidencia( Self%Element(i), j )
            Node => Self%Node(node_ID)
            CALL SET_Nodes_of_Element(Self%Element(i), Node, j)
        END DO
    END DO
    RETURN
    END SUBROUTINE CALC_Nodes_of_Element
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_Fi_Global(Self)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    ALLOCATE(Self%Fi_Global(3*Self%N_Nodes))
    END SUBROUTINE ALLOCATE_Fi_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_Fext_Global(Self)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    ALLOCATE(Self%Fext_Global(3*Self%N_Nodes))
    END SUBROUTINE ALLOCATE_Fext_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_g_Desbalanceamento(Self)
    IMPLICIT NONE
    TYPE(Type_problem), INTENT(INOUT)    ::  Self
    ALLOCATE(Self%F(3*Self%N_Nodes))
    END SUBROUTINE ALLOCATE_g_Desbalanceamento
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_Global_Hessian(Self)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    ALLOCATE(Self%H(3*Self%N_Nodes,3*Self%N_Nodes))
    END SUBROUTINE ALLOCATE_Global_Hessian
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_M_Global(Self)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    ALLOCATE(Self%M(3*Self%N_Nodes,3*Self%N_Nodes))
    END SUBROUTINE ALLOCATE_M_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Atualizar_pos_f(Self)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    Self%pos_f = Self%pos_f + Self%Delta_pos
    END SUBROUTINE Atualizar_pos_f
    !----------------------------------------------------------------------------------------------  
    SUBROUTINE CALC_Delta_pos_SOLVER(Self)
    USE MA67_INTERFACE
    IMPLICIT NONE
    TYPE(Type_Problem),       INTENT(INOUT)    ::    Self
    INTEGER*4 :: count, i,j, DOF
    INTEGER*4,  ALLOCATABLE :: i_sparce(:), j_sparce(:)
    REAL*8, ALLOCATABLE  ::  H_sparce(:)
    REAL*8, ALLOCATABLE, DIMENSION(:)  ::  g_Vetor
    REAL*8,  ALLOCATABLE, DIMENSION(:,:)    ::  Hessiana
    DOF=Self%N_DOF
    ALLOCATE (g_Vetor(DOF), Hessiana(DOF,DOF))

    
    g_Vetor  = Self%F
    Hessiana = Self%H
    count = 0
    DO i = 1, DOF
        DO j = i, DOF
            IF ( Hessiana(i, j) > 1.0d-15)   count = count + 1
        END DO
    END DO
    IF (allocated(i_sparce) .eqv. .false.) THEN
        ALLOCATE(i_sparce(count),j_sparce(count),H_sparce(count))
    ELSE
        DEALLOCATE(i_sparce,j_sparce,H_sparce)
        ALLOCATE(i_sparce(count),j_sparce(count),H_sparce(count))
    ENDIF
    count = 0
    DO i = 1, DOF
        DO j = i, DOF
            IF ( Hessiana(i, j) > 1.0d-15) THEN
                count = count + 1
                i_sparce(count) = i
                j_sparce(count) = j
                H_sparce(count) = Hessiana(i,j)
            END IF
        END DO
    END DO
    CALL MA67SOLVER(count,DOF,i_sparce, j_sparce, H_sparce, g_Vetor, Self%Delta_pos)
    
    END SUBROUTINE CALC_Delta_pos_SOLVER
    !----------------------------------------------------------------------------------------------   
    SUBROUTINE Write_Header2(Self)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*2                         ::  C1
    CHARACTER*2                         ::  C2
    CHARACTER*2                         ::  C3
    CHARACTER*2                         ::  C4
    CHARACTER*2                         ::  C5
    CHARACTER*2                         ::  C6
    CHARACTER*2                         ::  C7
    CHARACTER*2                         ::  C8
    C1 = 'PX'
    C2 = 'PY'
    C3 = 'PZ'
    C4 = 'DX'
    C5 = 'DY'
    C6 = 'DZ'
    C7 = 'DD'
    C8 = 'Fi'
    WRITE(4,1) C1, C2, C3, C4, C5, C6, C7, C8               
1   FORMAT(5X,A2,',',A2,',',A2,',',A2,',',A2,',',A2,',',A2',',A2)
    END SUBROUTINE Write_Header2
    !----------------------------------------------------------------------------------------------   
    SUBROUTINE Write_Step(Self,t)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*50                        ::  Nome_Arquivo
    CHARACTER*50                        ::  t1
    CHARACTER*8                         ::  fmt
    INTEGER*4                           ::  t
    INTEGER*4                           ::  i
    REAL*8                              ::  X
    REAL*8                              ::  Y
    REAL*8                              ::  Z
    REAL*8                              ::  DX
    REAL*8                              ::  DY
    REAL*8                              ::  DZ
    REAL*8                              ::  D
    REAL*8                              ::  Fi
    
    fmt = '(I5.5)'
    WRITE(t1,fmt) t
    Nome_Arquivo = (TRIM(Self%Nome_Arquivo)//'.'//TRIM(t1))
    CALL temp_open_file2(Nome_Arquivo)
    CALL Write_Header2(Self)
    DO i=1,Self%N_Nodes
        t = 3*i
        X = Self%pos_f(t-2)
        Y = Self%pos_f(t-1)
        Z = Self%pos_f(t)
        DX = X - Self%pos_0(t-2) 
        DY = Y - Self%pos_0(t-1)
        DZ = Z - Self%pos_0(t)
        D = DSQRT(DX**2 + DY**2 + DZ**2)
        Fi = DSQRT(Self%Fi_Global(t-2)**2 + Self%Fi_Global(t-1)**2 + Self%Fi_Global(t)**2)
        Write(4,10) X, Y, Z, DX, DY, DZ, D, Fi
10      FORMAT(5X,ES15.5e3,',',ES15.5e3,',',ES15.5e3,',',ES15.5e3,',',ES15.5e3,',',ES15.5e3,',',ES15.5e3',',ES15.5e3)
    END DO      
    
    CALL temp_close_file2
    END SUBROUTINE Write_Step
    

    

END MODULE Din_Program