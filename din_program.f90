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
    
    
    character*90                                ::  Tipo_Estrutura
    
    INTEGER*8                                   ::  Unidade
    INTEGER*8                                   ::  N_Nodes
    INTEGER*8                                   ::  N_Materiais
    INTEGER*8                                   ::  N_Secao
    INTEGER*8                                   ::  N_Element
    INTEGER*8                                   ::  N_Forcas
    INTEGER*8                                   ::  N_Cond_Contorno
    INTEGER*8                                   ::  N_DOF
    INTEGER*8                                   ::  N_dof_res
    INTEGER*8                                   ::  N_Plano_x0
    INTEGER*8                                   ::  N_Plano_y0
    INTEGER*8                                   ::  N_Plano_z0
    INTEGER*8                                   ::  N_Plano_xf
    INTEGER*8                                   ::  N_Plano_yf
    INTEGER*8                                   ::  N_Plano_zf
    REAL*8                                      ::  Norm_dpos                   ! Norma R2 do delta_pos
    REAL*8                                      ::  L(3)
    
    
    
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
    REAL*8                                      ::  Alpha
    REAL*8                                      ::  Beta                        ! Beta do Newton Raphson
    REAL*8                                      ::  Gama
    REAL*8                                      ::  Lambda
    REAL*8                                      ::  F_Reacao
    REAL*8, ALLOCATABLE                         ::  Finer(:)
    
    REAL*8, ALLOCATABLE                         ::  Fi_Global(:)                ! Vetor das forças Internas da Estrutura
    REAL*8, ALLOCATABLE                         ::  Fext_Global(:)              ! Vetor das forças Externas da Estrutura
    REAL*8, ALLOCATABLE                         ::  M(:,:)                      ! Matriz de Massa da Estrutura
    REAL*8, ALLOCATABLE                         ::  H(:,:)                      ! Matriz Hessiana da Estrutura
    REAL*8, ALLOCATABLE                         ::  MR(:,:)                     ! Matriz de Massa da Estrutura com restrições de deslocamento
    REAL*8, ALLOCATABLE                         ::  HR(:,:)                     ! Matriz Hessiana da Estrutura com restrições de deslocamento
    REAL*8, ALLOCATABLE                         ::  F(:)                        ! Vetor desbalanceado
    REAL*8, ALLOCATABLE                         ::  Reac(:,:)
    REAL*8, ALLOCATABLE                         ::  Reacao(:)  
    REAL*8, ALLOCATABLE                         ::  v_0(:)                      ! Velocidade inicial
    REAL*8, ALLOCATABLE                         ::  ac_0(:)                     ! initial acceleration
    REAL*8, ALLOCATABLE                         ::  pos_f_p(:)                  ! predictor position
    REAL*8, ALLOCATABLE                         ::  v_p(:)                      ! predictor velocity
    REAL*8, ALLOCATABLE                         ::  ac_p(:)                     ! predictor acceleration
    REAL*8, ALLOCATABLE                         ::  pos_f_c(:)                  ! corrector position
    REAL*8, ALLOCATABLE                         ::  v_c(:)                      ! corrector velocity
    REAL*8, ALLOCATABLE                         ::  ac_c(:)                     ! corrector acceleration
      
    REAL*8, ALLOCATABLE                         ::  pos_0(:)                    ! Posição inicial
    REAL*8, ALLOCATABLE                         ::  pos_f(:)                    ! Posição atual
    REAL*8, ALLOCATABLE                         ::  Delta_pos(:)                ! Delta pos: pos_f = pos_pred + Delta_pos
    
    CHARACTER*90                                ::  Nome_Arquivo                ! Nome do arquivo de saída
    
    INTEGER*8, ALLOCATABLE                      ::  Plano_x0(:)                 ! Nós no plano de deslocamento
    INTEGER*8, ALLOCATABLE                      ::  Plano_y0(:)                 ! Nós no plano de deslocamento
    INTEGER*8, ALLOCATABLE                      ::  Plano_z0(:)                 ! Nós no plano de deslocamento
    INTEGER*8, ALLOCATABLE                      ::  Plano_xf(:)                 ! Nós no plano de deslocamento
    INTEGER*8, ALLOCATABLE                      ::  Plano_yf(:)                 ! Nós no plano de deslocamento
    INTEGER*8, ALLOCATABLE                      ::  Plano_zf(:)                 ! Nós no plano de deslocamento
    REAL*8, ALLOCATABLE                         ::  T_rest(:,:)                 ! Matriz T de restição de deslocamento
    REAL*8, ALLOCATABLE                         ::  T_rest_tran(:,:)                 ! Matriz T de restição de deslocamento transposta
     
    
END TYPE Type_problem

    CONTAINS
    ! ==============================================================================================
    ! ----------------------------------- ESTÁTICO LINEAR ------------------------------------------
    ! ==============================================================================================
    
    SUBROUTINE Trelica_Est_Lin(Self,Output)
    IMPLICIT NONE
    TYPE (Type_Problem),     INTENT(INOUT)     :: Self
    character*90                               :: Output
    
    CALL Preprocessing(Self,Output)
    CALL CALC_L0_Element(Self)                               !Tamanho inicial do elemento
    CALL CALC_COSSENOS_GIR(Self)                             !Cossenos Giratórios de cada elemento
    CALL CALC_RIGIDEZ_ELEMENT(Self)                          !Matriz Hessiana de cada elemento
    CALL CALC_HESSIAN_GLOBAL(Self)                           !Matriz Hessiana da estrutura
    CALL CALC_Fext(Self)                                     !Calcular as forças externas
    CALL Aplicar_Cond_Cont_Hessian(Self)                     !Aplicar condições de contorno na matriz Hessiana da estrutura
    Self%F = -Self%Fext_Global
    CALL CALC_Delta_pos_Solver(Self)                         !Cálcula os deslocamentos
    CALL Write_VTK(Self,(0))
    CALL Atualizar_pos_f(Self)
    CALL Atualizar_Node_pos_f(Self)                          !Atualiza a posição final dos nós
    CALL Atualizar_Lf(Self)  
    CALL CALC_Fi_Element(Self)
    
    CALL Write_VTK(Self,(1))
    
    END SUBROUTINE Trelica_Est_Lin
    ! ==============================================================================================
    ! --------------------------- NEWTON RAPHSON DINÂMICO LINEAR -----------------------------------
    ! ==============================================================================================  
    SUBROUTINE Trelica_Din_Lin(Self,Output)
    IMPLICIT NONE
    TYPE (Type_Problem),     INTENT(INOUT)      ::  Self
    INTEGER*4                                   ::  Max_Iter
    INTEGER*4                                   ::  i
    INTEGER*4                                   ::  Stop_Iteration
    REAL*8                                      ::  Tolerance
    REAL*8                                      ::  Var
    character*90   ,     INTENT(IN)             ::  Output
    
    Tolerance = 1.0d-12
    Self%Relogio = 0.0d0
    MAX_Iter = 150
    
    CALL Preprocessing_Din_Lin(Self,Output)
    CALL Write_VTK(Self,(0))
    CALL CALC_L0_Element(Self)                                      !Tamanho inicial do elemento
    CALL CALC_COSSENOS_GIR(Self)                                    !Cossenos Giratórios de cada elemento
    CALL CALC_RIGIDEZ_ELEMENT(Self)                                 !Matriz Hessiana de cada elemento
    CALL CALC_HESSIAN_GLOBAL(Self)                                  !Calcular a Hessiana da estrutura
    
    CALL CALC_M_MASSA_ELEMENT_LUMPED(Self)                          !Calcular a Matriz de Massa de cada elemento
    
    CALL CALC_Fext(Self)


    
    DO i=1,Self%N_Dt
        
        Self%Relogio = Self%Relogio + Self%Dt
        Stop_Iteration = 0
        Var = 1.0d0  
        
        DO WHILE (( Var >= Tolerance) .AND. (Stop_Iteration<= MAX_Iter))
    
            Stop_Iteration = Stop_Iteration + 1
            CALL CALC_M_MASSA_GLOBAL(Self)                              !Calcular a Matriz de Massa da estrutura
            CALL CALC_HESSIAN_GLOBAL(Self)                              !Calcular a Hessiana da estrutura
            CALL CALC_HESSIAN_GLOBAL_Dyn(Self)                          !Calcular a Hessiana da estrutura para o sistema dinâmico
            
            CALL CALC_Fi_Element(Self)
            CALL CALC_Fi_Global(Self)
            CALL CALC_Fi_Dyn_Global(Self)
            CALL CALC_g_desbalanceamento(Self)
            CALL Aplicar_Cond_Cont_Hessian(Self)                        !Aplicar condições de contorno na matriz Hessiana da estrutura
            CALL Aplicar_Cond_Cont_Massa(Self)                          !Aplicar condições de contorno na matriz Hessiana da estrutura
            CALL Aplicar_Cond_Cont_Desbal(Self)                         !Aplicar condições de contorno na matriz Hessiana da estrutura
            CALL CALC_Delta_pos_Solver(Self)                            !Cálcula os deslocamentos
            CALL Atualizar_pos_f(Self)                                  !Atualiza a posição final
            CALL Atualizar_Node_pos_f(Self)
            
            IF (Stop_Iteration.GE.2) THEN
                Var = Norm(Self)
            END IF
        END DO
        CALL CALC_Corretivos(Self)                                      !Calcula a velocidade e acelerações corretivas
        CALL CALC_Preditivos(Self)                                      !Calcula a velocidade e acelerações preditivas
        CALL Write_VTK(Self,(i))
        
        Print *,'Tempo calculado em segundos:',' ', Self%Relogio
        Print *,'Numero de passos de tempo:',' ',i,'/',Self%N_Dt


    END DO
    
    
    END SUBROUTINE Trelica_Din_Lin
    ! ==============================================================================================
    ! -------------------------- NEWTON RAPHSON ESTÁTICO NÂO LINEAR --------------------------------
    ! ==============================================================================================  
    SUBROUTINE Trelica_Est_Nao_Lin(Self,Output,k)
    IMPLICIT NONE
    TYPE (Type_Problem),     INTENT(INOUT)      ::  Self
    INTEGER*4                                   ::  i
    INTEGER*4                                   ::  k
    INTEGER*4                                   ::  j
    INTEGER*4                                   ::  Dir
    INTEGER*4                                   ::  Max_Iter
    INTEGER*4                                   ::  Stop_Iteration
    INTEGER*4                                   ::  Stop_Iteration2
    REAL*8                                      ::  Tolerance
    REAL*8                                      ::  Var
    character*90                                ::  Output

    
    
    
        Tolerance = 1.0d-13
        Self%Relogio = 0.0d0
        MAX_Iter = 500000
        CALL Preprocessing_Nao_Lin(Self,Output)
        CALL CALC_L0_Element(Self)                                      !Tamanho inicial do elemento
        CALL CALC_COSSENOS_GIR(Self)
        CALL CALC_RIGIDEZ_ELEMENT(Self)   
        CALL Planos_de_deslocamento(Self)
!DO j=1, 1
!        Dir = j
        !CALL Create_T_rest(Self,Dir)
        CALL Write_VTK(Self,(0))
        CALL CALC_Fext(Self)
        !CALL Create_rest_Fext(Self)
    
        DO i=1,Self%N_Dt
        
            Self%Relogio = Self%Relogio + Self%Dt
            Stop_Iteration = 0
            Stop_Iteration2 = 0
            Var = 1.0d0
        
            DO WHILE (( Var >= Tolerance) .AND. (Stop_Iteration<= MAX_Iter))
           
                Stop_Iteration = Stop_Iteration + 1
                Stop_Iteration2 = Stop_Iteration2 + 1
                CALL Atualizar_Lf(Self)
                CALL Atualizar_Green(Self)
                CALL Atualizar_Primeira_Derivada(Self)                          !Calcular primeira derivada de E
                CALL CALC_HESSIAN_ELEMENT(Self)                                 !Calcular a Hessiana de cada elemento
                CALL CALC_HESSIAN_GLOBAL(Self)                                  !Calcular a Hessiana da estrutura
                CALL CALC_Fi_Element_N_Linear(Self)
                CALL CALC_Fi_Global(Self)
                !CALL Create_rest_Fi(Self)
                CALL CALC_g_desbalanceamento(Self)
                !CALL Create_rest_H(Self,Dir)
                
                CALL Aplicar_Cond_Cont_Hessian(Self)                            !Aplicar condições de contorno na matriz Hessiana da estrutura
                CALL Aplicar_Cond_Cont_desbal(Self)
                CALL CALC_Delta_pos_Solver(Self)
                !CALL Corrigir_delta_pos(Self)
                CALL Atualizar_pos_f(Self)
                CALL Atualizar_Node_pos_f(Self)
            
                Var = Norm(Self)
                IF ((Stop_iteration.GE.2000).AND.(Var.GE.1.0d-3)) THEN
                    Self%pos_f = Self%pos_0
                    CALL Atualizar_Node_pos_f(Self)
                    GO TO 123
                END IF
            
            END DO
        
        
        END DO
        

    !END DO
        
    CALL Atualizar_Lf(Self)
    CALL Atualizar_Green(Self)
    CALL Atualizar_Primeira_Derivada(Self)                          !Calcular primeira derivada de E
    CALL CALC_HESSIAN_ELEMENT(Self)                                 !Calcular a Hessiana de cada elemento
    CALL CALC_HESSIAN_GLOBAL(Self)                                  !Calcular a Hessiana da estrutura
    CALL CALC_Reacoes(Self)    
        
    IF (Stop_iteration.GE.MAX_iter) THEN
        Self%pos_f = Self%pos_0
        CALL Atualizar_Node_pos_f(Self)
    ELSE
        CALL Write_VTK(Self,(k))
    END IF
123 CONTINUE
    
    END SUBROUTINE Trelica_Est_Nao_Lin
    ! ==============================================================================================
    ! ------------------------ NEWTON RAPHSON DINÂMICO NÂO LINEAR ----------------------------------
    ! ==============================================================================================  
    SUBROUTINE Trelica_Din_Nao_Lin(Self,Output)
    IMPLICIT NONE
    TYPE (Type_Problem),     INTENT(INOUT)      ::  Self
    INTEGER*4                                   ::  Max_Iter
    INTEGER*4                                   ::  i
    INTEGER*4                                   ::  Stop_Iteration
    REAL*8                                      ::  Tolerance
    REAL*8                                      ::  Var
    REAL*8                                      ::  X
    REAL*8                                      ::  Y
    REAL*8                                      ::  Z
    character*90   ,     INTENT(IN)             ::  Output
    
    Tolerance = 1.0d-15
    Self%Relogio = 0.0d0
    MAX_Iter = 100000
    
    CALL Preprocessing_Din_Nao_Lin(Self,Output)
    CALL CALC_L0_Element(Self)                                      !Tamanho inicial do elemento
    CALL CALC_M_MASSA_ELEMENT_LUMPED(Self)                          !Calcular a Matriz de Massa de cada elemento
    CALL CALC_M_MASSA_GLOBAL(Self)                                  !Calcular a Matriz de Massa da estrutura
    CALL Write_VTK(Self,(0))


    
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
            CALL CALC_Fi_Element_N_Linear(Self)
            CALL CALC_Fi_Global(Self)
            CALL CALC_Fi_Dyn_Global(Self)
            CALL CALC_g_desbalanceamento(Self)
            CALL Aplicar_Cond_Cont_Hessian(Self)                        !Aplicar condições de contorno na matriz Hessiana da estrutura
            CALL Aplicar_Cond_Cont_Massa(Self)                          !Aplicar condições de contorno na matriz Hessiana da estrutura
            CALL Aplicar_Cond_Cont_Desbal(Self)                         !Aplicar condições de contorno na matriz Hessiana da estrutura
            CALL CALC_Delta_pos_Solver(Self)                            !Cálcula os deslocamentos
            CALL Atualizar_pos_f(Self)                                  !Atualiza a posição final
            CALL Atualizar_Node_pos_f(Self)
            
            IF (Stop_Iteration.GE.2) THEN
                Var = Norm(Self)
            END IF
        END DO
        CALL CALC_Corretivos(Self)                                      !Calcula a velocidade e acelerações corretivas
        CALL CALC_Preditivos(Self)                                      !Calcula a velocidade e acelerações preditivas
        CALL Write_VTK(Self,(i))
        
        Print *,'Tempo calculado em segundos:',' ', Self%Relogio
        Print *,'Numero de passos de tempo:',' ',i,'/',Self%N_Dt
        !Y = Self%v_c(5)
        !X = Self%Relogio
        !Z = Self%pos_f(5)
        !
        !WRITE(5,*) Y,X,Z

    END DO
    
    
    END SUBROUTINE Trelica_Din_Nao_Lin

    !----------------------------------------------------------------------------------------------
    SUBROUTINE Preprocessing(Self,Output)
    IMPLICIT NONE
    TYPE(Type_Problem),   INTENT(INOUT)             ::  Self
    character*90 , INTENT(IN)                       :: Output
    CALL write_dates(Self,Output)
    CALL SET_pos_0_Global(Self)
    CALL CALC_node_DOF(Self)
    CALL Inicial_pos_f_Global(Self)
    Self%Norm_Dpos = Inicial_Norm(Self)
    CALL Inicial_Fi_Global(Self)   
        
    END SUBROUTINE Preprocessing
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Preprocessing_Din_Lin(Self,Output)
    IMPLICIT NONE
    TYPE(Type_Problem),   INTENT(INOUT)             ::  Self
    INTEGER*4                                       ::  i
    character*90                                    :: Output
    CALL write_dates(Self,Output)
    CALL SET_pos_0_Global(Self)
    CALL CALC_node_DOF(Self)
    CALL Inicial_pos_f_Global(Self)
    Self%Norm_Dpos = Inicial_Norm(Self)
    CALL Inicial_Fi_Global(Self)  
    CALL ALLOCATE_g_Desbalanceamento(self)
    CALL Inicial_Dyn(Self)
    CALL CALC_Inicial_Preditivos(Self)
        
    END SUBROUTINE Preprocessing_Din_Lin
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Preprocessing_Nao_Lin(Self,Output)
    IMPLICIT NONE
    TYPE(Type_Problem),   INTENT(INOUT)             ::  Self
    character*90                                    :: Output
    
    CALL write_dates(Self,Output)
    CALL SET_pos_0_Global(Self)

    CALL CALC_node_DOF(Self)
    CALL CALC_L0_Element(Self)
    CALL Inicial_pos_f_Global(Self)
    Self%Norm_Dpos = Inicial_Norm(Self)
    Self%n_Dt = 1.0d0
    CALL Inicial_Fi_Global(Self)  
    CALL ALLOCATE_g_Desbalanceamento(self)
    CALL Inicial_Primeira_Derivada(Self)
        
    END SUBROUTINE Preprocessing_Nao_Lin
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Preprocessing_Din_Nao_Lin(Self,Output)
    IMPLICIT NONE
    TYPE(Type_Problem),   INTENT(INOUT)             ::  Self
    INTEGER*4                                       ::  i
    character*90                                    :: Output
    CALL write_dates(Self,Output)
    CALL SET_pos_0_Global(Self)
    CALL Inicial_pos_f_Global(Self)
    Self%Norm_Dpos = Inicial_Norm(Self)
    CALL Inicial_Fi_Global(Self)  
    CALL ALLOCATE_g_Desbalanceamento(self)
    CALL Inicial_Dyn(Self)
    CALL CALC_Inicial_Preditivos(Self)
    CALL Inicial_Primeira_Derivada(Self)
        
    END SUBROUTINE Preprocessing_Din_Nao_Lin
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Inicial_pos_f_Global(Self)
    IMPLICIT NONE
    TYPE(type_problem), INTENT(INOUT)    ::  Self
    INTEGER*4  :: i, j,gdof(3)
    REAL*8 :: desloc(3)
    CALL ALLOCATE_pos_f_Global(Self)
    Self%pos_f = Self%pos_0    
    Desloc = 0.0d0
    DO i=1, Self%N_nodes
        gdof(1) = GET_Node_DOF(Self%node(i), 1)
        gdof(2) = GET_Node_DOF(Self%node(i), 2)
        gdof(3) = GET_Node_DOF(Self%node(i), 3)
        Desloc = GET_Cond_Cont_desl(Self%Cond_contorno(i))
        Self%pos_f(Gdof(1)) = Self%pos_f(Gdof(1)) + Desloc(1)
        Self%pos_f(Gdof(2)) = Self%pos_f(Gdof(2)) + Desloc(2)
        Self%pos_f(Gdof(3)) = Self%pos_f(Gdof(3)) + Desloc(3)
        
        CALL SET_Node_pos_f(Self%node(i), Self%pos_f(Gdof(1)), Self%pos_f(Gdof(2)), Self%pos_f(Gdof(3))) 
    END DO
    END SUBROUTINE Inicial_pos_f_Global
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
    Self%V_c = Self%Ac_c*Self%Dt*Self%Gama + Self%V_p
    END SUBROUTINE CALC_V_c
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Ac_c(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    Self%Ac_c = ((Self%pos_f)/(Self%Beta * Self%Dt**2.0d0)) - Self%Ac_p
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
    Self%Ac_p = Self%pos_f/(Self%Beta * Self%Dt**2.0d0) + Self%V_c/(Self%Beta * Self%Dt) + Self%Ac_c * ( 1.0d0 / (2.0d0*Self%Beta) - 1.0d0)
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
    Self%Ac_p = Self%pos_f/(Self%Beta * Self%Dt**2.0d0) + Self%V_0/(Self%Beta * Self%Dt) + Self%Ac_0 * (( 1.0d0 / 2.0d0*Self%Beta) - 1.0d0)
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
    Self%Beta = 2.5d-1
    Self%Gama = 5.0d-1
    Self%Finer = 0.0d0
    Self%Lambda = 0.0d0
    
    END SUBROUTINE Inicial_Dyn
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_HESSIAN_GLOBAL_Dyn(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*2                               :: i
    INTEGER*2                               :: j

    do i = 1 , Self%N_DOF
        do j = 1 , Self%N_DOF
            Self%H(i,j) = Self%H(i,j) + (Self%M(i,j)/( Self%Beta * Self%Dt**2.0d0 ))
        end do
    end do
    
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
    SUBROUTINE CALC_M_MASSA_ELEMENT_LUMPED(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*4                              :: i
    DO i=1,Self%N_Element
        CALL SET_M_MASSA_Element_LUMPED(Self%Element(i))   
    END DO
    END SUBROUTINE CALC_M_MASSA_ELEMENT_LUMPED
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
    TYPE(Type_Problem),  INTENT(INOUT)      ::  SELF
    Self%F = 0.0d0
    Self%F = self%Fi_Global - self%Fext_Global
    END SUBROUTINE CALC_g_desbalanceamento
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Read_File(Self,Nome_Arquivo,unidade)
    IMPLICIT NONE
    TYPE(Type_Problem) , INTENT(INOUT)     :: Self                       !Representa o tipo de problema estrutural (treliça)
    character*90       , INTENT(IN)        :: Nome_Arquivo               !INPUT File_Name
    character*90                           :: buscar_palavra             !Palavra Auxiliar pra buscar
    INTEGER*4                              :: i                          !Índice para loops
    INTEGER*4                              :: Nodes_por_Element = 2      ! Número de nós por Element
    INTEGER*2                              :: unidade
    
    Self%unidade = unidade
    
    CALL Open_File(Nome_Arquivo,unidade)                                      !Subroutine: "Input_File"

    buscar_palavra = 'Tipo_Estrutura'                                    !Busca "Tipo_Estrutura" em 'Input_File'
    CALL Buscar(buscar_palavra,unidade)                                          !Subroutine: "Leitura"
    READ(unidade,*) Self%Tipo_Estrutura                                        !Lê as linhas embaixo de buscar palavra
    
    buscar_palavra = 'Coordenadas'
    Self%N_Nodes = counter(buscar_palavra,unidade)                               !Conta o número de nós e joga no N_Nodes dentro do Type_Problem
    ALLOCATE( Self%Node( Self%N_Nodes ) )                                !Aloca o Type_Nodes
    CALL Buscar(buscar_palavra,unidade)                                          !Volta o cursor para a primeira posição do nó
    DO i = 1, Self%N_Nodes
        CALL Ler_Nodes(Self%Node(i),unidade)                                   !Leitura das coordenadas iniciais dos nós
    END DO
    
    buscar_palavra = 'Material'
    Self%N_Materiais = counter(buscar_palavra,unidade)                           !Conta o número de materiais e joga no N_Materiais dentro do Type_Problem
    ALLOCATE( Self%Material( Self%N_Materiais ) )                        !Aloca o Type_Material
    CALL Buscar(buscar_palavra,unidade)                                          !Volta o cursor para a primeira posição do material
    DO i = 1, Self%N_Materiais
        CALL Ler_Material (Self%Material(i),unidade)                            !Leitura dos dados dos materiais
    END DO
    
    buscar_palavra = 'Secao'
    Self%N_Secao = counter(buscar_palavra,unidade)                          
    ALLOCATE( Self%Secao( Self%N_Secao ) )                       
    CALL Buscar(buscar_palavra, unidade)                                      
    DO i = 1, Self%N_Secao
        CALL Ler_Secao (Self%Secao (i),unidade)                          
    END DO
    
    
    buscar_palavra = 'Atribuir_Propriedades'
    Self%N_Element = counter(buscar_palavra,unidade)
    ALLOCATE( Self%Element( Self%N_Element ) )
    CALL Buscar(buscar_palavra, unidade)
    DO i = 1, Self%N_Element
        CALL Ler_Element(Self%Element(i),unidade)
        CALL SET_Nodes_by_Element(Self%Element(i), Nodes_por_Element)
    END DO
    CALL CALC_Material_do_Element(Self)
    CALL CALC_Secao_do_Element(Self)
    
    
    buscar_palavra = 'Incidencia'
    CALL Buscar(buscar_palavra, unidade)
    DO i = 1, Self%N_Element
        CALL Ler_Incidencia( Self%Element(i),unidade)                               ! Leitura da incidencia por elemento, nos locais por elemento
    END DO
    

    buscar_palavra = 'Forcas'
    Self%N_Forcas = counter(buscar_palavra,unidade)                          
    ALLOCATE( Self%Forcas( Self%N_Forcas ) )                       
    CALL Buscar(buscar_palavra, unidade)                                      
    DO i = 1, Self%N_Forcas
        CALL Ler_Forcas (Self%Forcas (i),unidade)                          
    END DO
    
    
    buscar_palavra = 'Cond_Contorno'
    Self%N_Cond_Contorno = counter(buscar_palavra,unidade)                          
    ALLOCATE( Self%Cond_Contorno( Self%N_Cond_Contorno ) )                       
    CALL Buscar(buscar_palavra, unidade)                                      
    DO i = 1, Self%N_Cond_Contorno
        CALL Ler_Cond_Contorno (Self%Cond_Contorno (i),unidade)                          
    END DO
    
    buscar_palavra = 'Tempo'
    Self%N_Dt = counter(buscar_palavra,unidade)                                            
    CALL Buscar(buscar_palavra, unidade)                                      
    CALL Ler_Tempo(Self,unidade)                         
    
    
    
    !CALL CALC_Element_por_node(Self)                                    ! Número de elementos que chegan ao nó
    CALL CALC_Num_DOF(Self)
    
    CALL CALC_node_dof(Self)
    DO i=1, Self%N_Element
        !CALL ALLOCATE_Nodes_of_Element(Self%Element(i))
        CALL ALLOCATE_M_Element(Self%Element(i))
        CALL ALLOCATE_H_Element(Self%Element(i))
        CALL ALLOCATE_B_Element(Self%Element(i))
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
    SUBROUTINE Ler_Tempo(Self,unidade)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)    ::   Self    
    Real*8                               ::   Dt, N_Dt  
    INTEGER*8                                 ::  unidade
    Read(unidade,*)Dt, N_Dt                                      !Leitura dos dados
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
    SUBROUTINE CALC_peso_element(Self)
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*4                              :: i
    DO i=1,Self%N_Element
        CALL SET_peso_Element(Self%Element(i))
    END DO
    END SUBROUTINE CALC_peso_element
    !----------------------------------------------------------------------------------------------    
    SUBROUTINE CALC_Stress_element(Self)
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*4                              :: i
    DO i=1,Self%N_Element
        CALL SET_Stress_Element(Self%Element(i))
    END DO
    END SUBROUTINE CALC_Stress_element
    !----------------------------------------------------------------------------------------------    
    SUBROUTINE CALC_Stress_Nao_Lin_element(Self)
    TYPE(Type_Problem),     INTENT(INOUT)  :: Self
    INTEGER*4                              :: i
    DO i=1,Self%N_Element
        CALL SET_Stress_Nao_Lin_Element(Self%Element(i))
    END DO
    END SUBROUTINE CALC_Stress_Nao_Lin_element
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
    INTEGER*4                              :: DOF
    INTEGER*4                              :: Node_ID
    Self%Fext_Global = 0.0d0
    DO i=1,Self%N_Forcas
        Node_ID = Get_Fext_ID(Self%Forcas(i))
        DOF = GET_Node_DOF(Self%Node(Node_ID), 1) - 1
        DO j=1,3
            Self%Fext_Global(DOF+j)= GET_Fext(Self%Forcas(i),j)
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
    TYPE(Type_Problem), INTENT(INOUT)       :: Self
    INTEGER*4                               :: i
    INTEGER*4                               :: j
    INTEGER*4                               :: k
    INTEGER*4                               :: dof1
    INTEGER*4                               :: dof2
    REAL*8, ALLOCATABLE                     :: M_elemento(:,:) 
    ALLOCATE (M_elemento(6,6))
    Self%M = 0.0d0
    M_elemento = 0.0d0
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
    SUBROUTINE Atualizar_COSSENOS_GIR(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)     :: Self
    INTEGER*4                              ::  i
    DO i=1, Self%N_Element
        CALL Atualizar_Alpha_Element(Self%Element(i))
        CALL Atualizar_Beta_Element(Self%Element(i))
        CALL Atualizar_Gama_Element(Self%Element(i))
    END DO
    END SUBROUTINE Atualizar_COSSENOS_GIR
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
    Self%Fi_Global = 0.0d0
    END SUBROUTINE Inicial_Fi_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Fi_Dyn_Global(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)  ::  SELF
    REAL*8, ALLOCATABLE                 ::  F_iner(:)
    REAL*8, ALLOCATABLE                 ::  M(:,:)
    INTEGER*4 ::  i
    INTEGER*4 ::  j
    REAL*8    ::  Lambda
    REAL*8    ::  Gamma
    REAL*8    ::  Beta
    REAL*8    ::  dt
    REAL*8    ::  Node_ID
    REAL*8    ::  DOF
    ALLOCATE(F_iner(Self%N_DOF), M(Self%N_DOF,Self%N_DOF))
    Gamma = self%Gama
    Beta  = self%Beta
    Lambda = self%Lambda
    dt = self%dt
    F_iner = 0.0d0
    F_iner = (( 1 / (Beta*dt**2.0d0) ) * self%pos_f)
    F_iner = F_iner + ( (Gamma*Lambda) / (Beta*dt) ) * self%pos_f 
    F_iner = F_iner - ( 1.0d0 + (Gamma*dt*Lambda) )*self%ac_p + (Lambda*Self%v_p)
    
    ! Matriz de massa completa
    DO i = 1, self%N_DOF
        DO j = i, self%N_DOF
            M(j,i) = self%M(i,j)
        END DO
    END DO
    
    self%Fi_Global = self%Fi_Global + MATMUL(M, F_iner)
    
    DO i = 1, Self%N_DOF
        IF (abs(Self%Fi_Global(i)) < 1.0d-20) THEN
        Self%Fi_Global(i) = 0.0d0
        END IF
    END DO
    
    RETURN
    END SUBROUTINE CALC_Fi_Dyn_Global
!---------------------------------------------------------------------------------------------------------------
    SUBROUTINE Calc_Fi_Global(Self)
    TYPE(Type_Problem),     INTENT(INOUT)   :: Self
    INTEGER*4                               :: i
    INTEGER*4                               :: j
    INTEGER*4                               :: dof1
    INTEGER*4                               :: dof2
    REAL*8,     DIMENSION(6)                :: Fi_element
    REAL*8, ALLOCATABLE                     :: Fi_global(:)
    ALLOCATE(Fi_global(Self%N_DOF))
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
        Den_Norm = Den_Norm + Self%pos_0(i)**2.0d0
    END DO
    END FUNCTION Inicial_Norm
    !----------------------------------------------------------------------------------------------

    FUNCTION Norm(SELF)      RESULT(Norm_dy)
    IMPLICIT NONE
    TYPE(Type_Problem),  INTENT(INOUT)      ::     SELF
    INTEGER*4   ::  i
    INTEGER*4   ::  DOF
    INTEGER*4   ::  r_x1
    INTEGER*4   ::  r_x2
    INTEGER*4   ::  r_x3
    REAL*8      ::  Norm_dy , Den_norm_dy

    DO i=1, self%N_Nodes
        r_x1    = GET_Cond_Cont(self%Cond_Contorno(i),1)
        IF ( r_x1 == 1 )  THEN
            DOF = GET_Node_DOF(self%node(i), 1)
            self%Delta_pos(DOF) = 0.0d0
        END IF
        r_x2    = GET_Cond_Cont(self%Cond_Contorno(i),2)
        IF ( r_x2 == 1 )  THEN
            DOF  = GET_Node_DOF(self%node(i), 2)
            self%Delta_pos(DOF) = 0.0d0
        END IF
        r_x3    = GET_Cond_Cont(self%Cond_Contorno(i),3)
        IF ( r_x3 == 1 )  THEN
            DOF  = GET_Node_DOF(self%node(i), 3)
            self%Delta_pos(DOF) = 0.0d0
        END IF
    END DO

    Norm_dy = 0.00d0
    DO i=1, 3*self%N_Nodes
        Norm_dy = Norm_dy + (self%Delta_pos(i))**2
    END DO
    Den_Norm_dy = 0.00d0
    DO i=1, 3*self%N_Nodes
        Den_Norm_dy = Den_Norm_dy + self%pos_0(i)**2
    END DO
    Norm_dy=DSQRT(Norm_dy)/DSQRT(Den_Norm_dy)
    !WRITE(*,'(f25.14)') Norm_dy
    END FUNCTION Norm
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
    SUBROUTINE ALLOCATE_pos_f_Global(Self)
    IMPLICIT NONE
    TYPE(type_problem), INTENT(INOUT)    ::  Self
    ALLOCATE(Self%pos_f(3*Self%N_Nodes))
    END SUBROUTINE ALLOCATE_pos_f_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_pos_0_Global(Self)
    IMPLICIT NONE
    TYPE(type_problem), INTENT(INOUT)    ::  Self
    ALLOCATE(Self%pos_0(3*Self%N_Nodes))
    END SUBROUTINE ALLOCATE_pos_0_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE Write_Dates(Self,Output)
    IMPLICIT NONE
    TYPE(type_problem)  , INTENT(INOUT)             ::  Self
    character*90                                    ::  Output
    Self%Nome_Arquivo = Output
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
    INTEGER*4, ALLOCATABLE                          ::  Element_por_node(:)        ! Número de elementos que chegan ao nó
    INTEGER*4                                       ::  i                          ! Contador para os Element
    INTEGER*4                                       ::  j                          ! Contador, nós locais no Element
    INTEGER*4                                       ::  Global_node                ! Número de nós
    INTEGER*4                                       ::  Nodes_por_Element = 2      ! Número de nós por Element, DEFAULT =2
    ALLOCATE(Element_por_node(Self%N_Nodes))
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
    TYPE(Type_Nodes), POINTER                   ::  Node2
    INTEGER*4   ::  Node_ID
    INTEGER*4   ::  Node_ID2
    INTEGER*4   ::  Nodes_by_Element
    INTEGER*4   ::  i
    INTEGER*4   ::  j
    DO i = 1, Self%N_Element
        Nodes_by_Element = GET_Nodes_by_Element(Self%Element(i))
        !CALL ALLOCATE_Nodes_of_Element(Self%Element(i))
        Node_ID = GET_incidencia( Self%Element(i), 1 )
        Node_ID2 = GET_incidencia( Self%Element(i), 2 )
        Node => Self%Node(node_ID)
        Node2 => Self%Node(node_id2)
        CALL SET_Nodes_of_Element(Self%Element(i), Node, Node2)
    END DO
    RETURN
    END SUBROUTINE CALC_Nodes_of_Element
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_Fi_Global(Self)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    IF (allocated(Self%Fi_Global) .eqv. .false.) THEN
        ALLOCATE(Self%Fi_Global(Self%N_DOF))
    ELSE
        DEALLOCATE(Self%Fi_Global)
        ALLOCATE(Self%Fi_Global(Self%N_DOF))
    ENDIF
    END SUBROUTINE ALLOCATE_Fi_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_Fext_Global(Self)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    IF (allocated(Self%Fext_Global) .eqv. .false.) THEN
        ALLOCATE(Self%Fext_Global(Self%N_DOF))
    ELSE
        DEALLOCATE(Self%Fext_Global)
        ALLOCATE(Self%Fext_Global(Self%N_DOF))
    ENDIF
    END SUBROUTINE ALLOCATE_Fext_Global
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_g_Desbalanceamento(Self)
    IMPLICIT NONE
    TYPE(Type_problem), INTENT(INOUT)    ::  Self
    IF (allocated(Self%F) .eqv. .false.) THEN
        ALLOCATE(Self%F(Self%N_DOF))
    ELSE
        DEALLOCATE(Self%F)
        ALLOCATE(Self%F(Self%N_DOF))
    ENDIF
    END SUBROUTINE ALLOCATE_g_Desbalanceamento
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_Global_Hessian(Self)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    IF (allocated(Self%H) .eqv. .false.) THEN
        ALLOCATE(Self%H(Self%N_DOF,Self%N_DOF))
    ELSE
        DEALLOCATE(Self%H)
        ALLOCATE(Self%H(Self%N_DOF,Self%N_DOF))
    ENDIF
    END SUBROUTINE ALLOCATE_Global_Hessian
    !----------------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_M_Global(Self)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    IF (allocated(Self%M) .eqv. .false.) THEN
        ALLOCATE(Self%M(Self%N_DOF,Self%N_DOF))
    ELSE
        DEALLOCATE(Self%M)
        ALLOCATE(Self%M(Self%N_DOF,Self%N_DOF))
    ENDIF
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
    INTEGER*4,  ALLOCATABLE     ::  i_sparce(:), j_sparce(:)
    REAL*8, ALLOCATABLE         ::  H_sparce(:)
    REAL*8, ALLOCATABLE         ::  g_Vetor(:)
    REAL*8,  ALLOCATABLE        ::  Hessiana(:,:)
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
    Self%Delta_pos = -1 * Self%Delta_pos
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
    SUBROUTINE Write_CSV(Self,t)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*90                        ::  Nome_Arquivo
    character*90                        ::  t1
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
    Nome_Arquivo = (TRIM(Self%Nome_Arquivo)//'.csv.'//TRIM(t1))
    CALL temp_open_file2(Nome_Arquivo)
    CALL Write_Header2(Self)
    DO i=1,Self%N_Nodes
        t = 3*i
        X = Self%pos_f(t-2.0d0)
        Y = Self%pos_f(t-1.0d0)
        Z = Self%pos_f(t)
        DX = X - Self%pos_0(t-2.0d0) 
        DY = Y - Self%pos_0(t-1.0d0)
        DZ = Z - Self%pos_0(t)
        D = DSQRT(DX**2.0d0 + DY**2.0d0 + DZ**2.0d0)
        Fi = DSQRT(Self%Fi_Global(t-2.0d0)**2.0d0 + Self%Fi_Global(t-1.0d0)**2.0d0 + Self%Fi_Global(t)**2.0d0)
        Write(4,10) X, Y, Z, DX, DY, DZ, D, Fi
10      FORMAT(5X,ES15.5e3,',',ES15.5e3,',',ES15.5e3,',',ES15.5e3,',',ES15.5e3,',',ES15.5e3,',',ES15.5e3',',ES15.5e3)
    END DO      
    
    CALL temp_close_file2
    END SUBROUTINE Write_CSV
    !----------------------------------------------------------------------------------------------   
    SUBROUTINE Write_Header(Self)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*26                        ::  C1
    character*90                        ::  C2
    CHARACTER*5                         ::  C3
    CHARACTER*25                        ::  C4
    CHARACTER*2                         ::  C5
    CHARACTER*2                         ::  C6
    CHARACTER*2                         ::  C7
    CHARACTER*2                         ::  C8
    C1 = '# vtk DataFile Version 3.0'
    C2 = Self%Nome_Arquivo
    C3 = 'ASCII'
    C4 = 'DATASET UNSTRUCTURED_GRID'
    WRITE(4,2) C1
2   FORMAT(A26)
    WRITE(4,3) C2
3   FORMAT(A50)
    WRITE(4,4) C3
4   FORMAT(A5)
    WRITE(4,*)
    WRITE(4,5) C4
5   FORMAT(A25)
    END SUBROUTINE Write_Header
    !----------------------------------------------------------------------------------------------   
    SUBROUTINE Write_VTK(Self,t)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*90                        ::  Nome_Arquivo
    character*90                        ::  t1
    CHARACTER*8                         ::  fmt
    INTEGER*4                           ::  t
    INTEGER*4                           ::  i
    INTEGER*4                           ::  ID
    INTEGER*4                           ::  Element2
    REAL*8                              ::  X
    REAL*8                              ::  Y
    REAL*8                              ::  Z
    REAL*8                              ::  ZY
    REAL*8                              ::  Fi(6)
    INTEGER*2                           ::  Ninc
    INTEGER*2                           ::  I1
    INTEGER*2                           ::  I2
    INTEGER*2                           ::  Xi
    INTEGER*2                           ::  Yi
    INTEGER*2                           ::  Zi
    
    
    fmt = '(I6.6)'
    WRITE(t1,fmt) t
    Nome_Arquivo = (TRIM(Self%Nome_Arquivo)//TRIM(t1)//'.vtk')
    CALL temp_open_file2(Nome_Arquivo)
    
    
    CALL Write_Header(Self)
    
    
    Write(4,6) Self%N_Nodes
6   FORMAT('POINTS ',I5.5,' double')
    
    
    
    
    DO i=1,Self%N_Nodes
        X = Get_xf(Self%Node(i)) + 10
        Y = Get_yf(Self%Node(i)) + 10
        Z = Get_zf(Self%Node(i)) + 10
        Write(4,10) X, Y, Z
10      FORMAT(ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
    END DO
    Write(4,*) 
    
    IF (Self%N_Element.LE.10500) THEN
        Ninc = 3 * Self%N_Element
    
    
        Write(4,20) Self%N_Element, Ninc
    20  FORMAT('CELLS ',I5.5,' ',I5.5)
    
        DO i=1,Self%N_Element
            I1 = GET_Incidencia(Self%Element(i),1) - 1
            I2 = GET_Incidencia(Self%Element(i),2) - 1
            Write(4,30) I1, I2 
    30      FORMAT('2 ',I5.5,' ',I5.5)
        END DO
        Write(4,*) 
    
        Write(4,40) Self%N_Element
    40  FORMAT('CELL_TYPES ',I5.5)
        DO i=1,Self%N_Element
            Write(4,50)
    50      FORMAT('3')
        END DO
        
    ELSE IF (Self%N_Element.LE.100000) THEN
        Element2 = 0
        DO i=1,Self%N_Element
            CALL CALC_COSSENOS_GIR(Self)
            X = Get_Alpha_Element(Self%Element(i))
            Y = Get_Beta_Element(Self%Element(i))
            Z = Get_Gama_Element(Self%Element(i))
            X=abs(X)
            Y=abs(Y)
            Z=abs(Z)
            IF ((X.EQ.1).OR.(Y.EQ.1).OR.(Z.EQ.1)) THEN
                Element2 = Element2 + 1
            END IF
        END DO
        
        Ninc = 3 * Element2
        Write(4,20) Element2, Ninc
        
        DO i=1,Self%N_Element
            X = Get_Alpha_Element(Self%Element(i))
            Y = Get_Beta_Element(Self%Element(i))
            Z = Get_Gama_Element(Self%Element(i))
            X=abs(X)
            Y=abs(Y)
            Z=abs(Z)
            IF ((X.EQ.1).OR.(Y.EQ.1).OR.(Z.EQ.1)) THEN
                I1 = GET_Incidencia(Self%Element(i),1) - 1
                I2 = GET_Incidencia(Self%Element(i),2) - 1
                Write(4,30) I1, I2 
            END IF
            
        END DO

        Write(4,*) 
    
        Write(4,40) Element2
        
        
        DO i=1,Element2
            Write(4,50)
        END DO
    ELSE 
        Element2 = Self%N_Plano_x0 + Self%N_Plano_y0 + Self%N_Plano_z0 + Self%N_Plano_xf + Self%N_Plano_yf + Self%N_Plano_zf
        
        Ninc = 3 * Element2
        Write(4,20) Element2, Ninc
        
        DO i=1,Self%N_Plano_x0
            ID = Self%Plano_x0(i)
            I1 = GET_Incidencia(Self%Element(ID),1) - 1
            I2 = GET_Incidencia(Self%Element(ID),2) - 1
            Write(4,30) I1, I2 
        END DO
        DO i=1,Self%N_Plano_y0
            ID = Self%Plano_y0(i)
            I1 = GET_Incidencia(Self%Element(ID),1) - 1
            I2 = GET_Incidencia(Self%Element(ID),2) - 1
            Write(4,30) I1, I2 
        END DO
        DO i=1,Self%N_Plano_z0
            ID = Self%Plano_z0(i)
            I1 = GET_Incidencia(Self%Element(ID),1) - 1
            I2 = GET_Incidencia(Self%Element(ID),2) - 1
            Write(4,30) I1, I2 
        END DO
        DO i=1,Self%N_Plano_xf
            ID = Self%Plano_xf(i)
            I1 = GET_Incidencia(Self%Element(ID),1) - 1
            I2 = GET_Incidencia(Self%Element(ID),2) - 1
            Write(4,30) I1, I2 
        END DO
        DO i=1,Self%N_Plano_yf
            ID = Self%Plano_yf(i)
            I1 = GET_Incidencia(Self%Element(ID),1) - 1
            I2 = GET_Incidencia(Self%Element(ID),2) - 1
            Write(4,30) I1, I2 
        END DO
        DO i=1,Self%N_Plano_zf
            ID = Self%Plano_zf(i)
            I1 = GET_Incidencia(Self%Element(ID),1) - 1
            I2 = GET_Incidencia(Self%Element(ID),2) - 1
            Write(4,30) I1, I2 
        END DO

        Write(4,*) 
    
        Write(4,40) Element2
        
        
        DO i=1,Element2
            Write(4,50)
        END DO
    END IF
    
    Write(4,*)
    Write(4,60) Self%N_Nodes
60  FORMAT('POINT_DATA ',I5.5)
    
    Write(4,70)
70  FORMAT('SCALARS deslocamentos_em_X double')
    
    Write(4,80)
80  FORMAT('LOOKUP_TABLE default')
    
    DO i=1,Self%N_Nodes
        X = Get_xf(Self%Node(i)) - Get_x0(Self%Node(i))
        WRITE(4,90) X
90      FORMAT(F22.16)
    END DO


    Write(4,102)
102 FORMAT('SCALARS deslocamentos_em_Y double')
    
    Write(4,80)
    
    DO i=1,Self%N_Nodes
        Y = Get_yf(Self%Node(i)) - Get_y0(Self%Node(i))
        WRITE(4,90) abs(Y)
    END DO
    
    
       Write(4,110)
110 FORMAT('SCALARS deslocamentos_em_Z double')
    
    Write(4,80)
    
    DO i=1,Self%N_Nodes
        Z = Get_zf(Self%Node(i)) - Get_z0(Self%Node(i))
        WRITE(4,90) abs(Z)
    END DO
    
        Write(4,115)
115 FORMAT('SCALARS deslocamentos_em_ZY double')
    
    Write(4,80)
    
    DO i=1,Self%N_Nodes
        ZY = SQRT((Get_zf(Self%Node(i)) - Get_z0(Self%Node(i)))**2 + (Get_yf(Self%Node(i)) - Get_y0(Self%Node(i)))**2)
        IF (ZY.LE.1.0d-14) THEN
            ZY = 0.0d0
        END IF
        WRITE(4,90) ZY
    END DO
    Write(4,116)
116 FORMAT('SCALARS deslocamentos_em_XY double')
    
    Write(4,80)
    
    DO i=1,Self%N_Nodes
        ZY = SQRT((Get_xf(Self%Node(i)) - Get_x0(Self%Node(i)))**2.0d0 + (Get_yf(Self%Node(i)) - Get_y0(Self%Node(i)))**2.0d0)
        IF (ZY.LE.1.0d-14) THEN
            ZY = 0.0d0
        END IF
        WRITE(4,90) ZY
    END DO
    Write(4,117)
117 FORMAT('SCALARS deslocamentos_em_XZ double')
    
    Write(4,80)
    
    DO i=1,Self%N_Nodes
        ZY = SQRT((Get_xf(Self%Node(i)) - Get_x0(Self%Node(i)))**2.0d0 + (Get_zf(Self%Node(i)) - Get_z0(Self%Node(i)))**2.0d0)
        IF (ZY.LE.1.0d-14) THEN
            ZY = 0.0d0
        END IF
        WRITE(4,90) ZY
    END DO
    
    
    
    Write(4,120)
120 FORMAT('VECTORS deslocamentos double')
    
    
    DO i=1,Self%N_Nodes
        X = Get_xf(Self%Node(i)) - Get_x0(Self%Node(i))
        Y = Get_yf(Self%Node(i)) - Get_y0(Self%Node(i))
        Z = Get_zf(Self%Node(i)) - Get_z0(Self%Node(i))
        WRITE(4,130) X, Y,Z
130     FORMAT(F22.16,' ',F22.16,' ',F22.16)
    END DO
 

    CALL temp_close_file2
    END SUBROUTINE Write_VTK
    !----------------------------------------------------------------------------------------------   
    SUBROUTINE CALC_L0_Cell(Self)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    INTEGER*4                           ::  i
    REAL*8                              ::  L(3)
    REAL*8                              ::  X
    REAL*8                              ::  Y
    REAL*8                              ::  Z
    
    L=0.0d0
    DO i=1, Self%N_Nodes
        X = Get_x0(Self%Node(i))
        IF (X.GT.L(1)) THEN
            L(1) = X
        END IF
        Y = Get_y0(Self%Node(i))
        IF (Y.GT.L(2)) THEN
            L(2) = Y
        END IF
        Z = Get_z0(Self%Node(i))
        IF (Z.GT.L(3)) THEN
            L(3) = Z
        END IF
    END DO
    
    Self%L = L 
     
    END SUBROUTINE CALC_L0_Cell

!!---------------------------------------------------------------------------------------------------------
    SUBROUTINE Planos_de_deslocamento(Self)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    INTEGER*8                           ::  i
    INTEGER*8                           ::  j
    INTEGER*8                           ::  k
    REAL*8                              ::  X1
    REAL*8                              ::  X2
    REAL*8                              ::  Y1
    REAL*8                              ::  Y2
    REAL*8                              ::  Z1
    REAL*8                              ::  Z2
    INTEGER*8                           ::  count_x
    INTEGER*8                           ::  count_y
    INTEGER*8                           ::  count_z
    INTEGER*8                           ::  count_x0
    INTEGER*8                           ::  count_y0
    INTEGER*8                           ::  count_z0
    INTEGER*8, ALLOCATABLE              ::  Plano_x0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_y0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_z0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_xf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_yf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_zf(:)
    REAL*8,    ALLOCATABLE              ::  Nodes(:,:)
    
    
    CALL CALC_L0_Cell(Self)
    
    count_x = 0
    count_y = 0
    count_z = 0
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    DO i=1, Self%N_Nodes
        X1 = Get_x0(Self%Node(i)) 
        Y1 = Get_y0(Self%Node(i)) 
        Z1 = Get_z0(Self%Node(i)) 
        X2 = Self%L(1)-X1
        Y2 = Self%L(2)-Y1
        Z2 = Self%L(3)-Z1
        IF (X2.EQ.0.0d0) THEN
            count_x = count_x + 1
        ELSE IF (X1.EQ.0.0d0) THEN
            count_x0 = count_x0 + 1
        END IF

        IF (Y2.EQ.0.0d0) THEN
            count_y = count_y + 1
        ELSE IF (Y1.EQ.0.0d0) THEN
            count_y0 = count_y0 + 1
        END IF
        
        IF (Z2.EQ.0.0d0) THEN
            count_z = count_z + 1
        ELSE IF (Z1.EQ.0.0d0) THEN
            count_z0 = count_z0 + 1
        END IF
        
    END DO
    
    
    ALLOCATE (Plano_x0(count_x0), Plano_y0(count_y0), Plano_z0(count_z0), Plano_xf(count_x), Plano_yf(count_y), Plano_zf(count_z))
    
    
    IF (Allocated(Self%Plano_x0).eqv. .false.) THEN
        ALLOCATE (Self%Plano_x0(count_x0), Self%Plano_y0(count_y0), Self%Plano_z0(count_z0), Self%Plano_xf(count_x), Self%Plano_yf(count_y), Self%Plano_zf(count_z))
    END IF
    
    
    count_x = 0
    count_y = 0
    count_z = 0
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    DO i=1, Self%N_Nodes
        X1 = Get_x0(Self%Node(i)) 
        Y1 = Get_y0(Self%Node(i)) 
        Z1 = Get_z0(Self%Node(i)) 
        X2 = Self%L(1)-X1
        Y2 = Self%L(2)-Y1
        Z2 = Self%L(3)-Z1
        
        IF (X2.EQ.0.0d0) THEN
            count_x = count_x + 1
            Plano_xf(count_x) = i
        ELSE IF (X1.EQ.0.0d0) THEN
            count_x0 = count_x0 + 1
            Plano_x0(count_x0) = i
        END IF
        
        
        IF (Y2.EQ.0.0d0) THEN
            count_y = count_y + 1
            Plano_yf(count_y) = i
        ELSE IF (Y1.EQ.0.0d0) THEN
            count_y0 = count_y0 + 1
            Plano_y0(count_y0) = i
        END IF
        
        IF (Z2.EQ.0.0d0) THEN
            count_z = count_z + 1
            Plano_zf(count_z) = i
        ELSE IF (Z1.EQ.0.0d0) THEN
            count_z0 = count_z0 + 1
            Plano_z0(count_z0) = i
        END IF
        
        
    END DO
    
    
    Self%Plano_x0 = Plano_x0
    Self%Plano_y0 = Plano_y0
    Self%Plano_z0 = Plano_z0
    Self%Plano_xf = Plano_xf
    Self%Plano_yf = Plano_yf
    Self%Plano_zf = Plano_zf
    
    Self%N_Plano_x0 = count_x0
    Self%N_Plano_y0 = count_y0
    Self%N_Plano_z0 = count_z0
    Self%N_Plano_xf = count_x
    Self%N_Plano_yf = count_y
    Self%N_Plano_zf = count_z
    
    
    END SUBROUTINE Planos_de_deslocamento
!------------------------------------------------------------------------------------------
    SUBROUTINE Create_T_rest(Self,Dir)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    INTEGER*4                           ::  N_sup(3)
    INTEGER*4                           ::  N_inf(3)
    INTEGER*4                           ::  Dir
    INTEGER*4                           ::  i
    INTEGER*4                           ::  j
    INTEGER*4                           ::  k
    INTEGER*4                           ::  N_DOF
    INTEGER*4                           ::  N_dof_res
    
    
    N_DOF = Self%N_DOF
    N_inf(1) = Self%N_Plano_x0
    N_inf(2) = Self%N_Plano_y0
    N_inf(3) = Self%N_Plano_z0
    N_sup(1) = Self%N_Plano_xf
    N_sup(2) = Self%N_Plano_yf
    N_sup(3) = Self%N_Plano_zf
    
    
    
    
    IF (ALLOCATED(Self%T_rest) .eqv. .TRUE.) THEN
        DEALLOCATE (Self%T_rest)
        DEALLOCATE (Self%T_rest_tran)
    END IF
    
    ALLOCATE (Self%T_rest(N_DOF,N_DOF))
    ALLOCATE (Self%T_rest_tran(N_DOF,N_DOF))
    
    
    
    Self%T_rest= 0
    Self%T_rest_tran = 0
    
    DO i=1, N_DOF
        Self%T_rest(i,i) = 1
    END DO
    
    
    IF (Dir.EQ.1) THEN
        k = 3 * Self%Plano_xf(1) - 2
        DO i=2, N_sup(Dir)
            j = (3 * Self%Plano_xf(i) - 2)
            Self%T_rest(j,j) = 0
            Self%T_rest(k,j) = 1
        END DO
        
        k = 3 * Self%Plano_x0(1) - 2
        DO i=2, N_inf(Dir)
            j = (3 * Self%Plano_x0(i) - 2)
            Self%T_rest(j,j) = 0
            Self%T_rest(k,j) = 1
        END DO
        
    ELSE IF (Dir.EQ.2) THEN
        
        k = 3 * Self%Plano_yf(1) - 1
        DO i=2, N_sup(Dir)
            j = (3 * Self%Plano_yf(i) - 1)
            Self%T_rest(j,j) = 0
            Self%T_rest(k,j) = 1
        END DO
        
        k = 3 * Self%Plano_y0(1) - 1
        DO i=2, N_inf(Dir)
            j = (3 * Self%Plano_y0(i) - 1)
            Self%T_rest(j,j) = 0
            Self%T_rest(k,j) = 1
        END DO
        
    ELSE IF (Dir.EQ.3) THEN
        
        k = 3 * Self%Plano_zf(1)
        DO i=2, N_sup(Dir)
            j = 3 * Self%Plano_zf(i)
            Self%T_rest(j,j) = 0
            Self%T_rest(k,j) = 1
        END DO
        
        k = 3 * Self%Plano_z0(1)
        DO i=2, N_inf(Dir)
            j = 3 * Self%Plano_z0(i)
            Self%T_rest(j,j) = 0
            Self%T_rest(k,j) = 1
        END DO
        
    END IF
            
            
    DO i=1, N_DOF
        DO j=1, N_DOF
            Self%T_rest_tran(i,j) = Self%T_rest(j,i)
        END DO
    END DO        
    
    
    
    
    END SUBROUTINE Create_T_rest
!------------------------------------------------------------------------------------------
    SUBROUTINE Create_rest_H(Self,Dir)
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    INTEGER*4                           ::  Dir
    INTEGER*4                           ::  i
    INTEGER*4                           ::  j
    INTEGER*4                           ::  N_DOF
    INTEGER*4                           ::  N_sup(3)
    INTEGER*4                           ::  N_inf(3)
    
    N_DOF = Self%N_DOF
    N_inf(1) = Self%N_Plano_x0
    N_inf(2) = Self%N_Plano_y0
    N_inf(3) = Self%N_Plano_z0
    N_sup(1) = Self%N_Plano_xf
    N_sup(2) = Self%N_Plano_yf
    N_sup(3) = Self%N_Plano_zf
    
    
    
    

    
    Self%H = MATMUL(Self%T_rest,Self%H)
    Self%H = MATMUL(Self%H,Self%T_rest)
    
    IF (Dir.EQ.1) THEN
    
        DO i=2, N_sup(Dir)
            j = (3 * Self%Plano_xf(i) - 2)
            Self%H(j,j) = 1
        END DO
        
        DO i=2, N_inf(Dir)
            j = (3 * Self%Plano_x0(i) - 2)
            Self%H(j,j) = 1
        END DO
        
    ELSE IF (Dir.EQ.2) THEN
        
        DO i=2, N_sup(Dir)
            j = (3 * Self%Plano_yf(i) - 1)
            Self%H(j,j) = 1
        END DO
        
        DO i=2, N_inf(Dir)
            j = (3 * Self%Plano_y0(i) - 1)
            Self%H(j,j) = 1
        END DO
        
    ELSE IF (Dir.EQ.3) THEN
        
        DO i=2, N_sup(Dir)
            j = 3 * Self%Plano_zf(i)
            Self%H(j,j) = 1
        END DO
        
        DO i=2, N_inf(Dir)
            j = 3 * Self%Plano_z0(i)
            Self%H(j,j) = 1
        END DO
        
    END IF
    
    END SUBROUTINE Create_rest_H
    
    
!------------------------------------------------------------------------------------------
    SUBROUTINE Create_rest_M(Self,Dir)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    INTEGER*4                           ::  Dir
    INTEGER*4                           ::  i
    INTEGER*4                           ::  j
    INTEGER*4                           ::  N_DOF
    INTEGER*4                           ::  N_sup(3)
    INTEGER*4                           ::  N_inf(3)
    REAL*8, ALLOCATABLE                 ::  MR_temp(:,:)
    
    N_DOF = Self%N_DOF
    N_inf(1) = Self%N_Plano_x0
    N_inf(2) = Self%N_Plano_y0
    N_inf(3) = Self%N_Plano_z0
    N_sup(1) = Self%N_Plano_xf
    N_sup(2) = Self%N_Plano_yf
    N_sup(3) = Self%N_Plano_zf
    
    
    
    
    ALLOCATE (MR_temp(N_DOF,N_DOF))
    
    MR_Temp = MATMUL(Self%T_rest,Self%M)
    Self%M = MATMUL(MR_Temp,Self%T_rest)
    
    IF (Dir.EQ.1) THEN

        DO i=2, N_sup(Dir)
            j = (3 * Self%Plano_xf(i) - 2)
            Self%M(j,j) = 1
        END DO
        
        DO i=2, N_inf(Dir)
            j = (3 * Self%Plano_x0(i) - 2)
            Self%M(j,j) = 1
        END DO
        
    ELSE IF (Dir.EQ.2) THEN
        
        DO i=2, N_sup(Dir)
            j = (3 * Self%Plano_yf(i) - 1)
            Self%M(j,j) = 1
        END DO
        
        DO i=2, N_inf(Dir)
            j = (3 * Self%Plano_y0(i) - 1)
            Self%M(j,j) = 1
        END DO
        
    ELSE IF (Dir.EQ.3) THEN
        
        DO i=2, N_sup(Dir)
            j = 3 * Self%Plano_zf(i)
            Self%M(j,j) = 1
        END DO
        
        DO i=2, N_inf(Dir)
            j = 3 * Self%Plano_z0(i)
            Self%M(j,j) = 1
        END DO
        
    END IF
    
    END SUBROUTINE Create_rest_M
!------------------------------------------------------------------------------------------------
    SUBROUTINE Create_rest_Fext(Self)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    
    Self%Fext_global = MATMUL(Self%T_rest,Self%Fext_global)
    
    
    END SUBROUTINE 
!------------------------------------------------------------------------------------------------
    SUBROUTINE Create_rest_Fi(Self)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    
    Self%Fi_global = MATMUL(Self%T_rest,Self%Fi_global)
    
    
    END SUBROUTINE 
!------------------------------------------------------------------------------------------
    SUBROUTINE Corrigir_delta_pos(Self)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    INTEGER*4                           ::  i
    INTEGER*4                           ::  j
    INTEGER*4                           ::  k
    INTEGER*4                           ::  N_DOF
    
    N_DOF = Self%N_DOF
    
    DO i=1, N_DOF
        IF (Self%T_rest_tran(i,i).EQ.0.0d0) THEN
            j=0
            DO WHILE (k.NE.1)
                j = j + 1
                k = Self%T_rest_tran(i,j)
            END DO
            Self%delta_pos(i) = Self%delta_pos(j)
            k = 0
        END IF
    END DO
    
    END SUBROUTINE Corrigir_delta_pos
!------------------------------------------------------------------------------------------------
    SUBROUTINE Cell_Repeater_x(Self,Output,Nx,Ny,Nz)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*90                        ::  Output
    INTEGER*8                           ::  Nx
    INTEGER*8                           ::  Ny
    INTEGER*8                           ::  Nz
    REAL*8                              ::  desloc
    INTEGER*8                           ::  i
    INTEGER*8                           ::  j
    INTEGER*8                           ::  k
    INTEGER*8                           ::  m
    INTEGER*8                           ::  n
    INTEGER*8                           ::  I1
    INTEGER*8                           ::  I2
    REAL*8                              ::  ID
    INTEGER*8                           ::  N_Cell
    INTEGER*8                           ::  element_total
    INTEGER*8                           ::  node_total
    REAL*8                              ::  X1
    REAL*8                              ::  X2
    REAL*8                              ::  Y1
    REAL*8                              ::  Y2
    REAL*8                              ::  Z1
    REAL*8                              ::  Z2
    REAL*8                              ::  dX1
    REAL*8                              ::  dX2
    REAL*8                              ::  dY1
    REAL*8                              ::  dY2
    REAL*8                              ::  dZ1
    REAL*8                              ::  dZ2
    REAL*8                              ::  Area
    REAL*8                              ::  Young
    REAL*8                              ::  Densidade
    REAL*8                              ::  Max_node(3)
    INTEGER*8                           ::  count_x
    INTEGER*8                           ::  count_y
    INTEGER*8                           ::  count_z
    INTEGER*8                           ::  count_x0
    INTEGER*8                           ::  count_y0
    INTEGER*8                           ::  count_z0
    INTEGER*8                           ::  count_aresta_x
    INTEGER*8                           ::  count_quina_x
    INTEGER*8                           ::  count_aresta_y
    INTEGER*8                           ::  count_quina_y
    INTEGER*8                           ::  count_aresta_z
    INTEGER*8                           ::  count_quina_z
    INTEGER*8                           ::  count_aresta_xy
    INTEGER*8                           ::  count_aresta_xz
    INTEGER*8                           ::  count_aresta_yx
    INTEGER*8                           ::  count_aresta_yz
    INTEGER*8                           ::  count_aresta_zx
    INTEGER*8                           ::  count_aresta_zy
    INTEGER*8                           ::  count_cell
    INTEGER*8, ALLOCATABLE              ::  Plano_x0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_y0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_z0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_xf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_yf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_zf(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_x(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_y(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_z(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_xy(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_xz(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_yx(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_yz(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_zx(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_zy(:)
    INTEGER*8, ALLOCATABLE              ::  quina_x(:)
    INTEGER*8, ALLOCATABLE              ::  quina_y(:)
    INTEGER*8, ALLOCATABLE              ::  quina_z(:)
    INTEGER*8, ALLOCATABLE              ::  Cell(:)
    REAL*8,    ALLOCATABLE              ::  Nodes(:,:)
    REAL*8,    ALLOCATABLE              ::  Nodes2(:,:)
    REAL*8,    ALLOCATABLE              ::  Elements(:,:)
    REAL*8,    ALLOCATABLE              ::  P(:,:,:,:)
    
    
    CALL CALC_L0_Cell(Self)
    
    
    !Colocando um ponto base para cada célula
    ALLOCATE (P(Nx,Ny,Nz,3))
    DO i=1, (Nx)
        DO j=1, (Ny)
            DO k=1, (Nz)
                P(i,j,k,1) = (i-1)*Self%L(1)
                P(i,j,k,2) = (j-1)*Self%L(2)
                P(i,j,k,3) = (k-1)*Self%L(3)
            END DO
        END DO
    END DO
    
    N_Cell = (i-1)*(j-1)*(k-1) !número de células
    
    Node_total = 2*N_Cell*Self%N_Nodes
    
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    count_cell = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ( ((X1.EQ.0) .AND. (X2.EQ.0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2)) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) ) THEN
            count_x0 = count_x0 + 1
        ELSE IF ( ((Y1.EQ.0) .AND. (Y2.EQ.0) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) )) THEN
            count_y0 = count_y0 + 1
        ELSE IF ( ((Z1.EQ.0) .AND. (Z2.EQ.0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) ) THEN
            count_z0 = count_z0 + 1
        ELSE IF ( ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_cell = count_cell + 1
        END IF
    END DO
    
    count_x = 0
    count_y = 0
    count_z = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ((X1.EQ.Self%L(1)) .AND. (X2.EQ.Self%L(1))) THEN
            count_x = count_x + 1
        ELSE IF ((Y1.EQ.Self%L(2)) .AND. (Y2.EQ.Self%L(2))) THEN
            count_y = count_y + 1
        ELSE IF ((Z1.EQ.Self%L(3)) .AND. (Z2.EQ.Self%L(3))) THEN
            count_z = count_z + 1
        END IF
    END DO
    
    !Salvando ID de cada elemento no seu vetor
    ALLOCATE (Plano_x0(count_x0), Plano_y0(count_y0), Plano_z0(count_z0), Cell(count_cell))
    ALLOCATE (Plano_xf(count_x), Plano_yf(count_y), Plano_zf(count_z))
    
    
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    count_cell = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ( ((X1.EQ.0) .AND. (X2.EQ.0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_x0 = count_x0 + 1
            Plano_x0(count_x0) = GET_Element_ID(Self%Element(i))
        ELSE IF ( ((Y1.EQ.0) .AND. (Y2.EQ.0) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) )) THEN
            count_y0 = count_y0 + 1
            Plano_y0(count_y0) = GET_Element_ID(Self%Element(i))
        ELSE IF ( ((Z1.EQ.0) .AND. (Z2.EQ.0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) ) THEN
            count_z0 = count_z0 + 1
            Plano_z0(count_z0) = GET_Element_ID(Self%Element(i))
        ELSE IF ( ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_cell = count_cell + 1
            Cell(count_cell) = GET_Element_ID(Self%Element(i))
        END IF
    END DO
    

    
    
    

    count_x = 0
    count_y = 0
    count_z = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ((X1.EQ.Self%L(1)) .AND. (X2.EQ.Self%L(1))) THEN
            count_x = count_x + 1
            Plano_xf(count_x) = GET_Element_ID(Self%Element(i))
        ELSE IF ((Y1.EQ.Self%L(2)) .AND. (Y2.EQ.Self%L(2))) THEN
            count_y = count_y + 1
            Plano_yf(count_y) = GET_Element_ID(Self%Element(i))
        ELSE IF ((Z1.EQ.Self%L(3)) .AND. (Z2.EQ.Self%L(3))) THEN
            count_z = count_z + 1
            Plano_zf(count_z) = GET_Element_ID(Self%Element(i))
        END IF
    END DO
    
    
    
    
    
    
    !Número total de elementos depois da repetição
    Element_total = N_cell*(count_cell + count_x0 + count_y0 + count_z0) + (Ny*Nz)*count_x + (Nx*Nz)*count_y + (Nx*Ny)*count_z
    ALLOCATE(Elements(Element_total,11))

    !   1 - ID elemento novo
    !   2 - X nó 1
    !   3 - Y nó 1
    !   4 - Z nó 1
    !   5 - X nó 2
    !   6 - Y nó 2
    !   7 - Z nó 2
    !   8 - ID Seção
    !   9 - ID Material
    !   10 - Incidência nova 1
    !   11 - Incidência nova 2
    
    Elements = 0.0d0
    
    
    !repetindo elementos
    n = 1
    DO i=1, Nx
        DO j=1, Ny
            DO k=1, Nz                
                
                    DO m=1, count_x0    !Repetindo plano x0
                        ID = Plano_x0(m)
                    
                        Elements(n,1) = n
                    
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                    
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                        n = n + 1
                    END DO
                IF (i.EQ.Nx) THEN
                    DO m=1, count_x
                        ID = Plano_xf(m)
                        
                        Elements(n,1) = n
                        
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                        
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                        n = n + 1
                    END DO
                END IF
                    
        
                    DO m=1, count_y0    !Repetindo plano y0
                        ID = Plano_y0(m)
                    
                        Elements(n,1) = n
                    
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                    
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))

                        n = n + 1
                    END DO
                IF (j.EQ.Ny) THEN
                    DO m=1, count_y        
                        ID = Plano_yf(m)
                        Elements(n,1) = n
                        
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                    
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))

                        n = n + 1
                    END DO
                END IF
                
                
                
                    DO m=1, count_z0    !Repetindo plano z0
                        ID = Plano_z0(m)
                    
                        Elements(n,1) = n
                    
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                    
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))

                        n = n + 1
                    END DO
                IF (k.EQ.Nz) THEN
                    DO m=1, count_z
                        ID = Plano_zf(m)
                        Elements(n,1) = n
                        
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                        
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))

                        n = n + 1
                    END DO
                END IF
                
                
                
                
                DO m=1, count_cell    !Repetindo interiores
                    ID = Cell(m)
                    
                    Elements(n,1) = n
                    
                    I1=GET_Incidencia(Self%Element(ID),1)
                    I2=GET_Incidencia(Self%Element(ID),2) 
                    X1 = Get_x0(Self%Node(I1))
                    X2 = Get_x0(Self%Node(I2))
                    Y1 = Get_y0(Self%Node(I1))
                    Y2 = Get_y0(Self%Node(I2))
                    Z1 = Get_z0(Self%Node(I1))
                    Z2 = Get_z0(Self%Node(I2))
                    
                    Elements(n,2) = P(i,j,k,1) + X1
                    Elements(n,3) = P(i,j,k,2) + Y1
                    Elements(n,4) = P(i,j,k,3) + Z1
                    
                    Elements(n,5) = P(i,j,k,1) + X2
                    Elements(n,6) = P(i,j,k,2) + Y2
                    Elements(n,7) = P(i,j,k,3) + Z2
                    
                    Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
                    Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                    n = n + 1
                END DO
    
            END DO
        END DO
    END DO

    
    
    
    ALLOCATE (Nodes2(node_total,7))
    Nodes2 = 0.0d0

    n = n-1
    j = 1
    i = 1

    dX1=Self%L(1)/50
    dY1=Self%L(2)/50
    dZ1=Self%L(3)/50
        
    DO j=1, element_total
            
            
        X1 = Elements(j,2) 
        Y1 = Elements(j,3)
        Z1 = Elements(j,4)

        DO k=1, (i-1)       !Teste se nó já existe
            X2 = abs(Nodes2(k,2) - X1)
            Y2 = abs(Nodes2(k,3) - Y1)
            Z2 = abs(Nodes2(k,4) - Z1)
            IF ((X2.LE.dX1) .AND. (Y2.LE.dY1) .AND. (Z2.LE.dZ1) ) THEN
                Elements(j,10) = k
                GO TO 115
            END IF
        END DO
        Nodes2(i,1) = i
        Nodes2(i,2) = X1
        Nodes2(i,3) = Y1
        Nodes2(i,4) = Z1
        Elements(j,10) = i
        i = i + 1
            
115     CONTINUE
            
        X1 = Elements(j,5)
        Y1 = Elements(j,6)
        Z1 = Elements(j,7)

        DO k=1, (i-1)       !Teste se nó já existe
            X2 = abs(Nodes2(k,2) - X1)
            Y2 = abs(Nodes2(k,3) - Y1)
            Z2 = abs(Nodes2(k,4) - Z1)    
            IF ((X2.LE.dX1) .AND. (Y2.LE.dY1) .AND. (Z2.LE.dZ1) ) THEN
                Elements(j,11) = k
                GO TO 116
            END IF
        END DO
        Nodes2(i,1) = i
        Nodes2(i,2) = X1
        Nodes2(i,3) = Y1
        Nodes2(i,4) = Z1
        Elements(j,11) = i
        i = i + 1
116     CONTINUE
        
    END DO
    !Print *,i,'/',element_total
    node_total = i - 1
    ALLOCATE (Nodes((i-1),7))
    DO j=1, node_total
        Nodes(j,:) = Nodes2(j,:)
    END DO
    DEALLOCATE (Nodes2)
    
    

   
    !Criando forças
    !Definindo o final em X, Y, Z
    Max_node(:) = P(Nx,Ny,Nz,:) + Self%L(:)
    
    !Contando quantos nós estão em cada final
    count_x = 0
    count_y = 0
    count_z = 0
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    count_aresta_x = 0
    count_quina_x = 0
    count_aresta_y = 0
    count_quina_y = 0
    count_aresta_z = 0
    count_quina_z = 0
    count_aresta_xy = 0
    count_aresta_xz = 0
    count_aresta_yx = 0
    count_aresta_yz = 0
    count_aresta_zx = 0
    count_aresta_zy = 0
    DO i=1, node_total
        X1 = Nodes(i,2) 
        Y1 = Nodes(i,3)
        Z1 = Nodes(i,4)
        X2 = Max_node(1)-X1
        Y2 = Max_node(2)-Y1
        Z2 = Max_node(3)-Z1
        IF (X2.EQ.0.0d0) THEN
            !aresta da força cisalhante
            IF (Y2.EQ.0.0d0) THEN
                count_aresta_xy = count_aresta_xy + 1
            END IF
            IF (Z2.EQ.0.0d0) THEN
                count_aresta_xz = count_aresta_xz + 1
            END IF
            !aresta da força de tração/compressão
            IF (( (Y2.LE.0.0d0).OR.(Y1.EQ.0.0d0)) .AND. ((Z2.EQ.0.0d0).OR.(Z1.EQ.0.0d0)) ) THEN
                count_quina_x = count_quina_x + 1                
            ELSE IF ((Y2.EQ.0.0d0) .OR. (Y1.EQ.0.0d0)) THEN
                count_aresta_x = count_aresta_x + 1
            ELSE IF ((Z2.EQ.0.0d0) .OR. (Z1.EQ.0.0d0)) THEN
                count_aresta_x = count_aresta_x + 1
            ELSE
                count_x = count_x + 1
            END IF
        ELSE IF (X1.EQ.0.0d0) THEN
            count_x0 = count_x0 + 1
        END IF

        IF (Y2.EQ.0.0d0) THEN
            !aresta da força cisalhante
            IF (X2.EQ.0.0d0) THEN
                count_aresta_yx = count_aresta_yx + 1
            END IF
            IF (Z2.EQ.0.0d0) THEN
                count_aresta_yz = count_aresta_yz + 1
            END IF
            !aresta da força de tração/compressão
            IF ( ((X2.EQ.0.0d0).OR.(X1.EQ.0.0d0)) .AND. ((Z2.EQ.0.0d0).OR.(Z1.EQ.0.0d0)) ) THEN
                count_quina_y = count_quina_y + 1                
            ELSE IF ((X2.EQ.0.0d0) .OR. (X1.EQ.0.0d0)) THEN
                count_aresta_y = count_aresta_y + 1
            ELSE IF ((Z2.EQ.0.0d0) .OR. (Z1.EQ.0.0d0)) THEN
                count_aresta_y = count_aresta_y + 1
            ELSE
                count_y = count_y + 1
            END IF
        ELSE IF (Y1.EQ.0.0d0) THEN
            count_y0 = count_y0 + 1
        END IF
        
        IF (Z2.EQ.0.0d0) THEN
            !aresta da força cisalhante
            IF (X2.EQ.0.0d0) THEN
                count_aresta_zx = count_aresta_zx + 1
            END IF
            IF (Y2.EQ.0.0d0) THEN
                count_aresta_zy = count_aresta_zy + 1
            END IF
            !aresta da força de tração/compressão
            IF ( ((Y2.EQ.0.0d0).OR.(Y1.EQ.0.0d0)) .AND. ((X2.EQ.0.0d0).OR.(X1.EQ.0.0d0)) ) THEN
                count_quina_z = count_quina_z + 1                
            ELSE IF ((Y2.EQ.0.0d0) .OR. (Y1.EQ.0.0d0)) THEN
                count_aresta_z = count_aresta_z + 1
            ELSE IF ((X2.EQ.0.0d0) .OR. (X1.EQ.0.0d0)) THEN
                count_aresta_z = count_aresta_z + 1
            ELSE
                count_z = count_z + 1
            END IF
        ELSE IF (Z1.EQ.0.0d0) THEN
            count_z0 = count_z0 + 1
        END IF
        
    END DO
    
    DEALLOCATE (Plano_x0, Plano_y0, Plano_z0, Plano_xf, Plano_yf, Plano_zf)
    ALLOCATE (Plano_x0(count_x0), Plano_y0(count_y0), Plano_z0(count_z0), Plano_xf(count_x), Plano_yf(count_y), Plano_zf(count_z))
    ALLOCATE (aresta_x(count_aresta_x), quina_x(count_quina_x))
    ALLOCATE (aresta_y(count_aresta_y), quina_y(count_quina_y))
    ALLOCATE (aresta_z(count_aresta_z), quina_z(count_quina_z))
    ALLOCATE (Aresta_xy(count_aresta_xy), Aresta_xz(count_aresta_xz))
    ALLOCATE (Aresta_yx(count_aresta_yx), Aresta_yz(count_aresta_yz))
    ALLOCATE (Aresta_zx(count_aresta_zx), Aresta_zy(count_aresta_zy))

    count_x = 0
    count_y = 0
    count_z = 0
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    count_aresta_x = 0
    count_quina_x = 0
    count_aresta_y = 0
    count_quina_y = 0
    count_aresta_z = 0
    count_quina_z = 0
    count_aresta_xy = 0
    count_aresta_xz = 0
    count_aresta_yx = 0
    count_aresta_yz = 0
    count_aresta_zx = 0
    count_aresta_zy = 0
    DO i=1, node_total
        X1 = Nodes(i,2) 
        Y1 = Nodes(i,3)
        Z1 = Nodes(i,4)
        X2 = Max_node(1)-X1
        Y2 = Max_node(2)-Y1
        Z2 = Max_node(3)-Z1
        
        IF (X2.EQ.0) THEN
            !aresta da força cisalhante
            IF (Y2.EQ.0) THEN
                count_aresta_xy = count_aresta_xy + 1
                Aresta_xy(count_aresta_xy) = i
            END IF
            IF (Z2.EQ.0) THEN
                count_aresta_xz = count_aresta_xz + 1
                Aresta_xz(count_aresta_xz) = i 
            END IF
            !aresta da força de tração/compressão
            IF ( ((Y2.EQ.0.0d0).OR.(Y1.EQ.0.0d0)) .AND. ((Z2.EQ.0.0d0).OR.(Z1.EQ.0.0d0)) ) THEN
                count_quina_x = count_quina_x + 1                
                quina_x(count_quina_x) = i
            ELSE IF ((Y2.EQ.0.0d0) .OR. (Y1.EQ.0.0d0)) THEN
                count_aresta_x = count_aresta_x + 1
                aresta_x(count_aresta_x) = i
            ELSE IF ((Z2.EQ.0.0d0) .OR. (Z1.EQ.0.0d0)) THEN
                count_aresta_x = count_aresta_x + 1
                aresta_x(count_aresta_x) = i
            ELSE
                count_x = count_x + 1
                Plano_xf(count_x) = i
            END IF
        ELSE IF (X1.EQ.0.0d0) THEN
            count_x0 = count_x0 + 1
            Plano_x0(count_x0) = i
        END IF
        
        
        IF (Y2.EQ.0.0d0) THEN
            !aresta da força cisalhante
            IF (X2.EQ.0.0d0) THEN
                count_aresta_yx = count_aresta_yx + 1
                Aresta_yx(count_aresta_yx) = i
            END IF
            IF (Z2.EQ.0.0d0) THEN
                count_aresta_yz = count_aresta_yz + 1
                Aresta_yz(count_aresta_yz) = i 
            END IF
            !aresta da força de tração/compressão
            IF ( ((X2.EQ.0.0d0).OR.(X1.EQ.0.0d0)) .AND. ((Z2.EQ.0.0d0).OR.(Z1.EQ.0.0d0)) ) THEN
                count_quina_y = count_quina_y + 1                
                quina_y(count_quina_y) = i
            ELSE IF ((X2.EQ.0.0d0) .OR. (X1.EQ.0.0d0)) THEN
                count_aresta_y = count_aresta_y + 1
                aresta_y(count_aresta_y) = i
            ELSE IF ((Z2.EQ.0.0d0) .OR. (Z1.EQ.0.0d0)) THEN
                count_aresta_y = count_aresta_y + 1
                aresta_y(count_aresta_y) = i
            ELSE
                count_y = count_y + 1
                Plano_yf(count_y) = i
            END IF
        ELSE IF (Y1.EQ.0.0d0) THEN
            count_y0 = count_y0 + 1
            Plano_y0(count_y0) = i
        END IF
        
        IF (Z2.EQ.0.0d0) THEN
            !aresta da força cisalhante
            IF (X2.EQ.0.0d0) THEN
                count_aresta_zx = count_aresta_zx + 1
                Aresta_zx(count_aresta_zx) = i
            END IF
            IF (Y2.EQ.0.0d0) THEN
                count_aresta_zy = count_aresta_zy + 1
                Aresta_zy(count_aresta_zy) = i 
            END IF
            !aresta da força de tração/compressão
            IF ( ((X2.EQ.0.0d0).OR.(X1.EQ.0.0d0)) .AND. ((Y2.EQ.0.0d0).OR.(Y1.EQ.0.0d0)) ) THEN
                count_quina_z = count_quina_z + 1                
                quina_z(count_quina_z) = i
            ELSE IF ((X2.EQ.0.0d0) .OR. (X1.EQ.0.0d0)) THEN
                count_aresta_z = count_aresta_z + 1
                aresta_z(count_aresta_z) = i
            ELSE IF ((Y2.EQ.0.0d0) .OR. (Y1.EQ.0.0d0)) THEN
                count_aresta_z = count_aresta_z + 1
                aresta_z(count_aresta_z) = i
            ELSE
                count_z = count_z + 1
                Plano_zf(count_z) = i
            END IF
        ELSE IF (Z1.EQ.0.0d0) THEN
            count_z0 = count_z0 + 1
            Plano_z0(count_z0) = i
        END IF
        
        
    END DO
    
    
    
        !=================================================================
    
    
    

      
    
    
    !Gravando no arquivo de saída
    
    
    CALL Temp_open_file4(Output,Self%unidade) !Unit 5
    
    
    WRITE(Self%unidade,*) 'Tipo_Estrutura'
    WRITE(Self%unidade,*) '-------------------------------------'
    WRITE(Self%unidade,*) 'TRUSS OU TRUSS'
    WRITE(Self%unidade,*) '-------------------------------------'
    WRITE(Self%unidade,*) 'TRUSS'
    WRITE(Self%unidade,*) '====================================='
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*) 'Coordenadas'
    WRITE(Self%unidade,*) '-------------------------------------'
    WRITE(Self%unidade,*) 'ID		x1		   y1		  z1'
    WRITE(Self%unidade,*) '-------------------------------------'
    DO i=1, node_total
        nodes(i,2) = nodes(i,2)
        nodes(i,3) = nodes(i,3)
        nodes(i,4) = nodes(i,4)
        WRITE(Self%unidade,100) Nodes(i,1), Nodes(i,2), Nodes(i,3), Nodes(i,4)
100     FORMAT(F7.0,' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
    END DO
    WRITE(Self%unidade,*) '====================================='
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*) 'Incidencia'
    WRITE(Self%unidade,*) '-------------------------------------'
    WRITE(Self%unidade,*) 'ID elem      Nó 1         Nó 2'
    WRITE(Self%unidade,*) '-------------------------------------'
    DO i=1, element_total
        WRITE(Self%unidade,150) Elements(i,1), Elements(i,10), Elements(i,11)
150     FORMAT(F7.0,' ',F7.0,' ',F7.0)
    END DO
    WRITE(Self%unidade,*) '====================================='
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*) 'Atribuir_Propriedades'
    WRITE(Self%unidade,*) '--------------------------------------'
    WRITE(Self%unidade,*) 'ID ELEM     Material  Seção   '
    WRITE(Self%unidade,*) '--------------------------------------'
    DO i=1, element_total
        WRITE(Self%unidade,150) Elements(i,1), Elements(i,9), Elements(i,8)
    END DO
    DEALLOCATE (Elements)
    WRITE(Self%unidade,*) '====================================='
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*) 'Material'
    WRITE(Self%unidade,*) '--------------------------------------------'
    WRITE(Self%unidade,*) 'ID  Mód. de elasticidade		  Densidade   '
    WRITE(Self%unidade,*) '--------------------------------------------'
    DO i=1, Self%N_Materiais
        ID = GET_material_ID(Self%Material(i))
        Young = GET_material_Young(Self%Material(i))
        Densidade = GET_material_Densidade(Self%Material(i))
        WRITE(Self%unidade,200) ID, Young, Densidade
200     FORMAT(F7.0,' ',ES15.7e3,' ',ES15.7e3)
    END DO
    WRITE(Self%unidade,*) '====================================='
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*) 'Secao'
    WRITE(Self%unidade,*) '-------------------------------------'
    WRITE(Self%unidade,*) 'ID			 Área em m²'
    WRITE(Self%unidade,*) '-------------------------------------'
    DO i=1, Self%N_Secao
        ID = GET_secao_ID(Self%Secao(i))
        Area = GET_secao_Area(Self%Secao(i))
        WRITE(Self%unidade,250) ID, Area
250     FORMAT(F7.0,' ',ES15.7e3)
    END DO
    WRITE(Self%unidade,*) '====================================='
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*) 'Cond_Contorno'
    WRITE(Self%unidade,*) '------------------------------------'
    WRITE(Self%unidade,*) 'ID Nó      X    r_x    Y    r_y     Z     r_z (1=restrito, 0=livre)'
    WRITE(Self%unidade,*) '------------------------------------'
    
    Nodes(:,2) = 0.0d0
    Nodes(:,3) = 0.0d0
    Nodes(:,4) = 0.0d0
    Nodes(:,5) = 0.0d0 !r_x
    Nodes(:,6) = 0.0d0 !r_y
    Nodes(:,7) = 0.0d0 !r_z
    DO i=1, count_x0
        ID = Plano_x0(i)
        nodes(ID,2) = 1
        !nodes(ID,3) = 1
        !nodes(ID,4) = 1
    END DO
    
    
    desloc = Self%L(1)*0.15d0
    
    
    DO i=1, count_x
        ID = Plano_xf(i)
        nodes(ID,2) = 1
        !nodes(ID,3) = 1
        !nodes(ID,4) = 1
        nodes(ID,5) = desloc
    END DO
    DO i=1, count_y0
        ID = Plano_y0(i)
        nodes(ID,3) = 1
    END DO
    DO i=1, count_z0
        ID = Plano_z0(i)
        nodes(ID,4) = 1
    END DO
    DO i=1, count_aresta_x
        ID = aresta_x(i)
        nodes(ID,2) = 1
        !nodes(ID,3) = 1
        !nodes(ID,4) = 1
        nodes(ID,5) = desloc
    END DO
    DO i=1, 4
        ID = quina_x(i)
        nodes(ID,2) = 1
        !nodes(ID,3) = 1
        !nodes(ID,4) = 1
        nodes(ID,5) = desloc
    END DO
    nodes(1,2) = 1
    nodes(1,3) = 1
    nodes(1,4) = 1
    
    
    
    DO i=1,node_total
        WRITE(Self%unidade,300) Nodes(i,1), Nodes(i,2), Nodes(i,5), Nodes(i,3), Nodes(i,6), Nodes(i,4), Nodes(i,7)
300     FORMAT(F7.0,'  ',F7.0,' ',ES15.7e3,' ',F7.0,' ',ES15.7e3,' ',F7.0,' ',ES15.7e3)
!        WRITE(Self%unidade,300) Nodes(i,1), Nodes(i,2), Nodes(i,3), Nodes(i,4)
!300     FORMAT(F7.0,'  ',F7.0,' ',F7.0,' ',F7.0)
    END DO
    WRITE(Self%unidade,*) '====================================='
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*) 
    WRITE(Self%unidade,*) 'Forcas'
    WRITE(Self%unidade,*) '------------------------------------'
    WRITE(Self%unidade,*) 'ID Nó      X        Y        Z'
    WRITE(Self%unidade,*) '------------------------------------'
        X1 = GET_Fext(Self%Forcas(1),1)
        Y1 = GET_Fext(Self%Forcas(1),2)
        Z1 = GET_Fext(Self%Forcas(1),3)
    DO i=1, count_x
        ID = plano_xf(i)
        WRITE(Self%unidade,350) Nodes(ID,1), X1, Y1, Z1
350     FORMAT(F7.0,' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
    END DO
    DO i=1, count_aresta_x
        ID = aresta_x(i)
        X2 = X1/2
        Y2 = Y1/2
        Z2 = Z1/2
        WRITE(Self%unidade,350) Nodes(ID,1), X2, Y2, Z2
    END DO
    DO i=1, 4
        ID = quina_x(i)
        X2 = X1/4
        Y2 = Y1/4
        Z2 = Z1/4
        WRITE(Self%unidade,350) Nodes(ID,1), X2, Y2, Z2
    END DO
    WRITE(Self%unidade,*) '====================================='
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*) 
    WRITE(Self%unidade,*) 'Tempo '
    WRITE(Self%unidade,*) '--------------------------------------------------------'
    WRITE(Self%unidade,*) 'Dt (incremento)				N_Dt (Número de incrementos)'
    WRITE(Self%unidade,*) '--------------------------------------------------------'
	WRITE(Self%unidade,*) Self%DT, Self%N_Dt
400 FORMAT(F7.0,' ',ES15.7e3)
    WRITE(Self%unidade,*) '========================================================'
    WRITE(Self%unidade,*)
    WRITE(Self%unidade,*) 'Fim_do_Arquivo'
    Close(Self%unidade)
    
    

    
    END SUBROUTINE Cell_Repeater_x
    !----------------------------------------------------------------------------------------------   
    SUBROUTINE Cell_Creator(Self,Output,Nx,Ny,Nz,gen_size)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*90                        ::  Output
    INTEGER*8                           ::  Nx
    INTEGER*8                           ::  Ny
    INTEGER*8                           ::  Nz
    INTEGER*8                           ::  i
    INTEGER*8                           ::  j
    INTEGER*8                           ::  k
    INTEGER*8                           ::  m
    INTEGER*8                           ::  n
    INTEGER*8                           ::  t
    INTEGER*8                           ::  I1
    INTEGER*8                           ::  I2
    REAL*8                              ::  ID
    INTEGER*8                           ::  N_Cell
    INTEGER*8                           ::  element_total
    INTEGER*8                           ::  node_total
    REAL*8                              ::  X1
    REAL*8                              ::  X2
    REAL*8                              ::  Y1
    REAL*8                              ::  Y2
    REAL*8                              ::  Z1
    REAL*8                              ::  Z2
    REAL*8                              ::  dY2
    REAL*8                              ::  dX2
    REAL*8                              ::  dZ2
    REAL*8                              ::  dY1
    REAL*8                              ::  dX1
    REAL*8                              ::  dZ1
    REAL*8                              ::  Area
    REAL*8                              ::  Young
    REAL*8                              ::  Densidade
    REAL*8                              ::  Max_node(3)
    INTEGER*8                           ::  count_x
    INTEGER*8                           ::  count_y
    INTEGER*8                           ::  count_z
    INTEGER*8                           ::  count_x0
    INTEGER*8                           ::  count_y0
    INTEGER*8                           ::  count_z0
    INTEGER*8                           ::  count_aresta_x
    INTEGER*8                           ::  count_quina
    INTEGER*8                           ::  count_cell
    INTEGER*8, ALLOCATABLE              ::  Plano_x0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_y0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_z0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_xf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_yf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_zf(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_x(:)
    INTEGER*8, ALLOCATABLE              ::  quina(:)
    INTEGER*8, ALLOCATABLE              ::  Cell(:)
    REAL*8,    ALLOCATABLE              ::  Nodes(:,:)
    REAL*8,    ALLOCATABLE              ::  Nodes2(:,:)
    REAL*8,    ALLOCATABLE              ::  Elements(:,:)
    REAL*8,    ALLOCATABLE              ::  P(:,:,:,:)
    INTEGER*8, ALLOCATABLE              ::  gen_map(:) 
    INTEGER*8                           ::  gen_size
    
    
    
    CALL CALC_L0_Cell(Self)
    
    
    !Colocando um ponto base para cada célula
    ALLOCATE (P(Nx,Ny,Nz,3))
    DO i=1, (Nx)
        DO j=1, (Ny)
            DO k=1, (Nz)
                P(i,j,k,1) = (i-1)*Self%L(1)
                P(i,j,k,2) = (j-1)*Self%L(2)
                P(i,j,k,3) = (k-1)*Self%L(3)
            END DO
        END DO
    END DO
    
    N_Cell = Nx*Ny*Nz !número de células
    
    Node_total = 2*N_Cell*Self%N_Nodes
    
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    count_cell = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ( ((X1.EQ.0.0d0) .AND. (X2.EQ.0.0d0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2)) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) ) THEN
            count_x0 = count_x0 + 1
        ELSE IF ( ((Y1.EQ.0.0d0) .AND. (Y2.EQ.0.0d0) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) )) THEN
            count_y0 = count_y0 + 1
        ELSE IF ( ((Z1.EQ.0.0d0) .AND. (Z2.EQ.0.0d0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) ) THEN
            count_z0 = count_z0 + 1
        ELSE IF ( ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_cell = count_cell + 1
        END IF
    END DO
    
    count_x = 0
    count_y = 0
    count_z = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ((X1.EQ.Self%L(1)) .AND. (X2.EQ.Self%L(1))) THEN
            count_x = count_x + 1
        ELSE IF ((Y1.EQ.Self%L(2)) .AND. (Y2.EQ.Self%L(2))) THEN
            count_y = count_y + 1
        ELSE IF ((Z1.EQ.Self%L(3)) .AND. (Z2.EQ.Self%L(3))) THEN
            count_z = count_z + 1
        END IF
    END DO
    
    !Salvando ID de cada elemento no seu vetor
    ALLOCATE (Plano_x0(count_x0), Plano_y0(count_y0), Plano_z0(count_z0), Cell(count_cell))
    ALLOCATE (Plano_xf(count_x), Plano_yf(count_y), Plano_zf(count_z))
    
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    count_cell = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ( ((X1.EQ.0.0d0) .AND. (X2.EQ.0.0d0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_x0 = count_x0 + 1
            Plano_x0(count_x0) = GET_Element_ID(Self%Element(i))
        ELSE IF ( ((Y1.EQ.0.0d0) .AND. (Y2.EQ.0.0d0) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) )) THEN
            count_y0 = count_y0 + 1
            Plano_y0(count_y0) = GET_Element_ID(Self%Element(i))
        ELSE IF ( ((Z1.EQ.0.0d0) .AND. (Z2.EQ.0.0d0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) ) THEN
            count_z0 = count_z0 + 1
            Plano_z0(count_z0) = GET_Element_ID(Self%Element(i))
        ELSE IF ( ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_cell = count_cell + 1
            Cell(count_cell) = GET_Element_ID(Self%Element(i))
        END IF
    END DO
    

    
    
    

    count_x = 0
    count_y = 0
    count_z = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ((X1.EQ.Self%L(1)) .AND. (X2.EQ.Self%L(1))) THEN
            count_x = count_x + 1
            Plano_xf(count_x) = GET_Element_ID(Self%Element(i))
        ELSE IF ((Y1.EQ.Self%L(2)) .AND. (Y2.EQ.Self%L(2))) THEN
            count_y = count_y + 1
            Plano_yf(count_y) = GET_Element_ID(Self%Element(i))
        ELSE IF ((Z1.EQ.Self%L(3)) .AND. (Z2.EQ.Self%L(3))) THEN
            count_z = count_z + 1
            Plano_zf(count_z) = GET_Element_ID(Self%Element(i))
        END IF
    END DO
    
    
    
    
    
    
    !Número total de elementos depois da repetição
    Element_total = N_cell*(count_cell + count_x0 + count_y0 + count_z0) + (Ny*Nz)*count_x + (Nx*Nz)*count_y + (Nx*Ny)*count_z
    ALLOCATE(Elements(Element_total,11))
    
    !   1 - ID elemento novo
    !   2 - X nó 1
    !   3 - Y nó 1
    !   4 - Z nó 1
    !   5 - X nó 2
    !   6 - Y nó 2
    !   7 - Z nó 2
    !   8 - ID Seção
    !   9 - ID Material
    !   10 - Incidência nova 1
    !   11 - Incidência nova 2
    
    Elements = 0.0d0
    
    !Tamanho do gene e mapeamento genético
    Gen_size = element_total-((Ny*Nz)*count_x + (Nx*Nz)*count_y + (Nx*Ny)*count_z)
    ALLOCATE (Gen_map(element_total))
    
    
    !repetindo elementos
    n = 1
    t = 1
    DO i=1, Nx
        DO j=1, Ny
            DO k=1, Nz                
                
                    DO m=1, count_x0    !Repetindo plano x0
                        ID = Plano_x0(m)
                    
                        Elements(n,1) = n
                    
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                    
                        Elements(n,8) = t   !Cell creator
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                        Gen_map(n) = 0
                        n = n + 1
                        t = t + 1
                    END DO
                IF (i.EQ.Nx) THEN
                    DO m=1, count_x
                        ID = Plano_xf(m)
                        
                        Elements(n,1) = n
                        
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                        
                        Elements(n,8) = 0
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                        Gen_map(n) = 1
                        n = n + 1
                    END DO
                END IF
                    
        
                    DO m=1, count_y0    !Repetindo plano y0
                        ID = Plano_y0(m)
                    
                        Elements(n,1) = n
                    
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                    
                        Elements(n,8) = t !Cell creator
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                        Gen_map(n) = 0
                        t = t + 1
                        n = n + 1
                    END DO
                IF (j.EQ.Ny) THEN
                    DO m=1, count_y        
                        ID = Plano_yf(m)
                        Elements(n,1) = n
                        
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                    
                        Elements(n,8) = 0
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                        Gen_map(n) = 1
                        n = n + 1
                    END DO
                END IF
                
                
                
                    DO m=1, count_z0    !Repetindo plano z0
                        ID = Plano_z0(m)
                    
                        Elements(n,1) = n
                    
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                    
                        Elements(n,8) = t !Cell creator
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                        Gen_map(n) = 0
                        t = t + 1
                        n = n + 1
                    END DO
                IF (k.EQ.Nz) THEN
                    DO m=1, count_z
                        ID = Plano_zf(m)
                        Elements(n,1) = n
                        
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                        
                        Elements(n,8) = 0
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                        Gen_map(n) = 1
                        n = n + 1
                    END DO
                END IF
                
                
                
                
                DO m=1, count_cell    !Repetindo interiores
                    ID = Cell(m)
                    
                    Elements(n,1) = n
                    
                    I1=GET_Incidencia(Self%Element(ID),1)
                    I2=GET_Incidencia(Self%Element(ID),2) 
                    X1 = Get_x0(Self%Node(I1))
                    X2 = Get_x0(Self%Node(I2))
                    Y1 = Get_y0(Self%Node(I1))
                    Y2 = Get_y0(Self%Node(I2))
                    Z1 = Get_z0(Self%Node(I1))
                    Z2 = Get_z0(Self%Node(I2))
                    
                    Elements(n,2) = P(i,j,k,1) + X1
                    Elements(n,3) = P(i,j,k,2) + Y1
                    Elements(n,4) = P(i,j,k,3) + Z1
                    
                    Elements(n,5) = P(i,j,k,1) + X2
                    Elements(n,6) = P(i,j,k,2) + Y2
                    Elements(n,7) = P(i,j,k,3) + Z2
                    
                    Elements(n,8) = t !Cell creator
                    
                    Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                    Gen_map(n) = 0
                    t = t + 1
                    n = n + 1
                END DO
    
            END DO
        END DO
    END DO
    
    
    
    
    
    
    
    
    
    ALLOCATE (Nodes2(node_total,4))
    Nodes2 = 0.0d0

    n = n-1
    j = 1
    i = 1
    
    dX1=Self%L(1)/10
    dY1=Self%L(2)/10
    dZ1=Self%L(3)/10
        
        
    DO j=1, element_total
            
            
        X1 = Elements(j,2) 
        Y1 = Elements(j,3)
        Z1 = Elements(j,4)

        DO k=1, (i-1)       !Teste se nó já existe
            X2 = abs(Nodes2(k,2) - X1)
            Y2 = abs(Nodes2(k,3) - Y1)
            Z2 = abs(Nodes2(k,4) - Z1)
            IF ((X2.LE.dX1) .AND. (Y2.LE.dY1) .AND. (Z2.LE.dZ1) ) THEN
                Elements(j,10) = k
                GO TO 15
            END IF
        END DO
        Nodes2(i,1) = i
        Nodes2(i,2) = X1
        Nodes2(i,3) = Y1
        Nodes2(i,4) = Z1
        Elements(j,10) = i
        i = i + 1
            
15      CONTINUE
            
        X1 = Elements(j,5)
        Y1 = Elements(j,6)
        Z1 = Elements(j,7)

        DO k=1, (i-1)       !Teste se nó já existe
            X2 = abs(Nodes2(k,2) - X1)
            Y2 = abs(Nodes2(k,3) - Y1)
            Z2 = abs(Nodes2(k,4) - Z1)    
            IF ((X2.LE.dX1) .AND. (Y2.LE.dY1) .AND. (Z2.LE.dZ1) ) THEN
                Elements(j,11) = k
                GO TO 16
            END IF
        END DO
        Nodes2(i,1) = i
        Nodes2(i,2) = X1
        Nodes2(i,3) = Y1
        Nodes2(i,4) = Z1
        Elements(j,11) = i
        i = i + 1
16      CONTINUE
        Print *,i,'/',element_total
    END DO
    
    node_total = i - 1
    ALLOCATE (Nodes(i,4))
    DO j=1, node_total
        Nodes(j,:) = Nodes2(j,:)
    END DO
    DEALLOCATE (Nodes2)
    
   
    !Criando forças
    !Definindo o final em X, Y, Z
    Max_node(:) = P(Nx,Ny,Nz,:) + Self%L(:)
    
    
    !Contando quantos nós estão em cada final
    count_x = 0
    count_y = 0
    count_z = 0
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    count_aresta_x = 0
    count_quina = 0
    DO i=1, node_total
        X1 = Nodes(i,2) 
        Y1 = Nodes(i,3)
        Z1 = Nodes(i,4)
        
        IF (Max_node(1).EQ.X1) THEN
            IF ( ((Max_node(2).EQ.Y1).OR.(Y1.EQ.0.0d0)) .AND. ((Max_node(3).EQ.Z1).OR.(Z1.EQ.0.0d0)) ) THEN
                count_quina = count_quina + 1                
            ELSE IF ((Max_node(2).EQ.Y1) .OR. (Y1.EQ.0.0d0)) THEN
                count_aresta_x = count_aresta_x + 1
            ELSE IF ((Max_node(3).EQ.Z1) .OR. (Z1.EQ.0.0d0)) THEN
                count_aresta_x = count_aresta_x + 1
            ELSE
                count_x = count_x + 1
            END IF
        ELSE IF (X1.EQ.0.0d0) THEN
            count_x0 = count_x0 + 1
        END IF
        
        IF (Max_node(2).EQ.Y1) THEN
            count_y = count_y +1
        ELSE IF (Y1.EQ.0.0d0) THEN
            count_y0 = count_y0 + 1
        END IF
        
        IF (Max_node(3).EQ.Z1) THEN
            count_z = count_z +1
        ELSE IF (Y1.EQ.0.0d0) THEN
            count_z0 = count_z0 + 1
        END IF
        
    END DO
    
    DEALLOCATE (Plano_x0, Plano_y0, Plano_z0, Plano_xf, Plano_yf, Plano_zf)
    ALLOCATE (Plano_x0(count_x0), Plano_y0(count_y0), Plano_z0(count_z0), Plano_xf(count_x), Plano_yf(count_y), Plano_zf(count_z))
    ALLOCATE (aresta_x(count_aresta_x), quina(count_quina))

    count_x = 0
    count_y = 0
    count_z = 0
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    count_aresta_x = 0
    count_quina = 0
    DO i=1, node_total
        X1 = Nodes(i,2) 
        Y1 = Nodes(i,3)
        Z1 = Nodes(i,4)
        
        IF (Max_node(1).EQ.X1) THEN
            IF ( ((Max_node(2).EQ.Y1).OR.(Y1.EQ.0.0d0)) .AND. ((Max_node(3).EQ.Z1).OR.(Z1.EQ.0.0d0)) ) THEN
                count_quina = count_quina + 1                
                quina(count_quina) = i
            ELSE IF ((Max_node(2).EQ.Y1) .OR. (Y1.EQ.0.0d0)) THEN
                count_aresta_x = count_aresta_x + 1
                aresta_x(count_aresta_x) = i
            ELSE IF ((Max_node(3).EQ.Z1) .OR. (Z1.EQ.0.0d0)) THEN
                count_aresta_x = count_aresta_x + 1
                aresta_x(count_aresta_x) = i
            ELSE
                count_x = count_x + 1
                Plano_xf(count_x) = i
            END IF
        ELSE IF (X1.EQ.0.0d0) THEN
            count_x0 = count_x0 + 1
            Plano_x0(count_x0) = i
        END IF
        
        IF (Max_node(2).EQ.Y1) THEN
            count_y = count_y +1
            Plano_yf(count_y) = i
        ELSE IF (Y1.EQ.0.0d0) THEN
            count_y0 = count_y0 + 1
            Plano_y0(count_y0) = i
        END IF
        
        IF (Max_node(3).EQ.Z1) THEN
            count_z = count_z +1
            Plano_zf(count_z) = i
        ELSE IF (Y1.EQ.0.0d0) THEN
            count_z0 = count_z0 + 1
            Plano_z0(count_z0) = i
        END IF
        
        
    END DO
    
    
    !corrigindo a simetria de fronteiras da célula final
    

       
    DO i=1, element_total
        j = gen_map(i)
        IF (j.EQ.1) THEN
            X1 = Elements(i,2)
            X2 = Elements(i,5)
            Y1 = Elements(i,3)
            Y2 = Elements(i,6)
            Z1 = Elements(i,4)
            Z2 = Elements(i,7)
            IF ((X1.EQ.max_node(1)) .AND. (X2.EQ.max_node(1))) THEN
                DO k=1, element_total
                    dX1 = abs(Elements(k,2)) !igual a zero
                    dX2 = abs(Elements(k,5)) !igual a zero
                    dY1 = abs(Elements(k,3) - Y1)
                    dY2 = abs(Elements(k,6) - Y2)
                    dZ1 = abs(Elements(K,4) - Z1)
                    dZ2 = abs(Elements(k,7) - Z2)
                    IF ((dX1.EQ.0.0d0).AND.(dX2.EQ.0.0d0).AND.(dY1.EQ.0.0d0).AND.(dY2.EQ.0.0d0).AND.(dZ1.EQ.0.0d0).AND.(dZ2.EQ.0.0d0)) THEN
                        Elements(i,8) =  Elements(k,8)
                        GO TO 37
                    ELSE
                        dX2 = abs(Elements(k,2)) !igual a zero
                        dX1 = abs(Elements(k,5)) !igual a zero
                        dY2 = abs(Elements(k,3) - Y2)
                        dY1 = abs(Elements(k,6) - Y1)
                        dZ2 = abs(Elements(K,4) - Z2)
                        dZ1 = abs(Elements(k,7) - Z1)
                        IF ((dX1.EQ.0.0d0).AND.(dX2.EQ.0.0d0).AND.(dY1.EQ.0.0d0).AND.(dY2.EQ.0.0d0).AND.(dZ1.EQ.0.0d0).AND.(dZ2.EQ.0.0d0)) THEN
                            Elements(i,8) =  Elements(k,8)
                            GO TO 37
                        END IF
                    END IF
                END DO
                
                
            ELSE IF ((Y1.EQ.max_node(2)) .AND. (Y2.EQ.max_node(2))) THEN
                DO k=1, element_total
                    dX1 = abs(Elements(k,2) - X1) 
                    dX2 = abs(Elements(k,5) - X2) 
                    dY1 = abs(Elements(k,3)) !igual a zero
                    dY2 = abs(Elements(k,6)) !igual a zero
                    dZ1 = abs(Elements(K,4) - Z1)
                    dZ2 = abs(Elements(k,7) - Z2)
                    IF ((dX1.EQ.0.0d0).AND.(dX2.EQ.0.0d0).AND.(dY1.EQ.0.0d0).AND.(dY2.EQ.0.0d0).AND.(dZ1.EQ.0.0d0).AND.(dZ2.EQ.0.0d0)) THEN
                        Elements(i,8) = Elements(k,8)
                        GO TO 37
                    ELSE
                        dX2 = abs(Elements(k,2) - X2) 
                        dX1 = abs(Elements(k,5) - X1) 
                        dY2 = abs(Elements(k,3)) !igual a zero
                        dY1 = abs(Elements(k,6)) !igual a zero
                        dZ2 = abs(Elements(K,4) - Z2)
                        dZ1 = abs(Elements(k,7) - Z1)
                        IF ((dX1.EQ.0.0d0).AND.(dX2.EQ.0.0d0).AND.(dY1.EQ.0.0d0).AND.(dY2.EQ.0.0d0).AND.(dZ1.EQ.0.0d0).AND.(dZ2.EQ.0.0d0)) THEN
                            Elements(i,8) = Elements(k,8)
                            GO TO 37
                        END IF
                    END IF
                END DO    
            ELSE IF ((Z1.EQ.max_node(3)) .AND. (Z2.EQ.max_node(3))) THEN
                DO k=1, element_total
                    dX1 = abs(Elements(k,2) - X1) 
                    dX2 = abs(Elements(k,5) - X2) 
                    dY1 = abs(Elements(k,3) - Y1)
                    dY2 = abs(Elements(k,6) - Y2)
                    dZ1 = abs(Elements(K,4)) !igual a zero
                    dZ2 = abs(Elements(k,7)) !igual a zero
                    IF ((dX1.EQ.0.0d0).AND.(dX2.EQ.0.0d0).AND.(dY1.EQ.0.0d0).AND.(dY2.EQ.0.0d0).AND.(dZ1.EQ.0.0d0).AND.(dZ2.EQ.0.0d0)) THEN
                        Elements(i,8) = Elements(k,8)
                        GO TO 37
                    ELSE
                        dX2 = abs(Elements(k,2) - X2) 
                        dX1 = abs(Elements(k,5) - X1) 
                        dY2 = abs(Elements(k,3) - Y2)
                        dY1 = abs(Elements(k,6) - Y1)
                        dZ2 = abs(Elements(K,4)) !igual a zero
                        dZ1 = abs(Elements(k,7)) !igual a zero
                        IF ((dX1.EQ.0.0d0).AND.(dX2.EQ.0.0d0).AND.(dY1.EQ.0.0d0).AND.(dY2.EQ.0.0d0).AND.(dZ1.EQ.0.0d0).AND.(dZ2.EQ.0.0d0)) THEN
                            Elements(i,8) = Elements(k,8)
                            GO TO 37
                        END IF
                    END IF
                END DO  
            END IF
        END IF
37      CONTINUE
    END DO

    
    
    !Gravando no arquivo de saída
    
    CALL Temp_open_file3(Output) !Unit 5
    
    
    WRITE(5,*) 'Tipo_Estrutura'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS OU TRUSS'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS'
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Coordenadas'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID		x1		   y1		  z1'
    WRITE(5,*) '-------------------------------------'
    DO i=1, node_total
        WRITE(5,100) Nodes(i,1), Nodes(i,2), Nodes(i,3), Nodes(i,4)
100     FORMAT(F7.0,' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Incidencia'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID elem      Nó 1         Nó 2'
    WRITE(5,*) '-------------------------------------'
    DO i=1, element_total
        WRITE(5,150) Elements(i,1), Elements(i,10), Elements(i,11)
150     FORMAT(F7.0,' ',F7.0,' ',F7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Atribuir_Propriedades'
    WRITE(5,*) '--------------------------------------'
    WRITE(5,*) 'ID ELEM     Material  Seção   '
    WRITE(5,*) '--------------------------------------'
    DO i=1, element_total
        WRITE(5,155) Elements(i,1), Elements(i,8)
155     FORMAT(F7.0,'  ','1.',' ',F7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Material'
    WRITE(5,*) '--------------------------------------------'
    WRITE(5,*) 'ID  Mód. de elasticidade		  Densidade   '
    WRITE(5,*) '--------------------------------------------'
    DO i=1, Self%N_Materiais
        ID = GET_material_ID(Self%Material(i))
        Young = GET_material_Young(Self%Material(i))
        Densidade = GET_material_Densidade(Self%Material(i))
        WRITE(5,200) ID, Young, Densidade
200     FORMAT(F7.0,' ',ES15.7e3,' ',ES15.7e3)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Secao'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID			 Área em m²'
    WRITE(5,*) '-------------------------------------'
    DO i=1, gen_size              !Cell creator
        WRITE(5,250) i
250     FORMAT(I7.0,'     0')
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Cond_Contorno'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z   (1=restrito, 0=livre)'
    WRITE(5,*) '------------------------------------'
    
    Nodes(:,2) = 0
    Nodes(:,3) = 0
    Nodes(:,4) = 0
    DO i=1, count_x0
        ID = Plano_x0(i)
        nodes(ID,2) = 1
        nodes(ID,3) = 1
        nodes(ID,4) = 1
    END DO
    
    DO i=1,node_total
        WRITE(5,300) Nodes(i,1), Nodes(i,2), Nodes(i,3), Nodes(i,4)
300     FORMAT(F7.0,' ',F7.0,' ',F7.0,' ',F7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Forcas'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z'
    WRITE(5,*) '------------------------------------'
        X1 = GET_Fext(Self%Forcas(1),1)
        Y1 = GET_Fext(Self%Forcas(1),2)
        Z1 = GET_Fext(Self%Forcas(1),3)
        WRITE(5,350)  X1, Y1, Z1
350     FORMAT('  1.','  ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Tempo '
    WRITE(5,*) '--------------------------------------------------------'
    WRITE(5,*) 'Dt (incremento)				N_Dt (Número de incrementos)'
    WRITE(5,*) '--------------------------------------------------------'
	WRITE(5,*) Self%DT, Self%N_Dt
400 FORMAT(ES15.7e3,' ',ES15.7e3)
    WRITE(5,*) '========================================================'
    WRITE(5,*)
    WRITE(5,*) 'Fim_do_Arquivo'
    DEALLOCATE (Elements)
    Close(5)
    
    
    END SUBROUTINE Cell_Creator    

    
    
!----------------------------------------------------------------------------------------------   
    SUBROUTINE Conversor_de_gene(Self,Output,Gen)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*90                        ::  Output
    INTEGER*8                           ::  i
    INTEGER*8                           ::  j
    INTEGER*8                           ::  k
    INTEGER*8                           ::  m
    INTEGER*8                           ::  n
    INTEGER*8                           ::  gen_size
    INTEGER*8                           ::  ID
    INTEGER*8                           ::  I1
    INTEGER*8                           ::  I2
    INTEGER*8                           ::  element_total
    INTEGER*8                           ::  node_total
    REAL*8                              ::  X1
    REAL*8                              ::  X2
    REAL*8                              ::  Y1
    REAL*8                              ::  Y2
    REAL*8                              ::  Z1
    REAL*8                              ::  Z2
    REAL*8                              ::  Area
    REAL*8                              ::  Young
    REAL*8                              ::  Densidade
    REAL*8                              ::  Max_node(3)
    INTEGER*8                           ::  count_x
    INTEGER*8                           ::  count_y
    INTEGER*8                           ::  count_z
    INTEGER*8                           ::  count_x0
    INTEGER*8                           ::  count_y0
    INTEGER*8                           ::  count_z0
    INTEGER*8                           ::  count_aresta_x
    INTEGER*8                           ::  count_quina
    INTEGER*8                           ::  count_cell
    INTEGER*8, ALLOCATABLE              ::  Plano_x0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_y0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_z0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_xf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_yf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_zf(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_x(:)
    INTEGER*8, ALLOCATABLE              ::  quina(:)
    INTEGER*8, ALLOCATABLE              ::  Cell(:)
    REAL*8,    ALLOCATABLE              ::  Nodes(:,:)
    REAL*8,    ALLOCATABLE              ::  Nodes2(:,:)
    REAL*8,    ALLOCATABLE              ::  Elements(:,:)
    REAL*8                              ::  Gen(:)

    
    gen_size = SIZE(Gen)
    
    !Gravando no arquivo de saída
    
    CALL Temp_open_file3(Output) !Unit 5
100 FORMAT(I7.0,' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
150 FORMAT(I7.0,' ',I7.0,' ',I7.0)
200 FORMAT(I7.0,' ',ES15.7e3,' ',ES15.7e3)
350 FORMAT(I7.0,' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
400 FORMAT(ES15.7e3,' ',ES15.7e3)
    
    WRITE(5,*) 'Tipo_Estrutura'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS OU TRUSS'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS'
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Coordenadas'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID		x1		   y1		  z1'
    WRITE(5,*) '-------------------------------------'
    DO i=1, Self%N_Nodes
        X1 = Get_x0(Self%Node(i))
        Y1 = Get_y0(Self%Node(i))
        Z1 = Get_z0(Self%Node(i))
        WRITE(5,100) i, X1, Y1, Z1
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Incidencia'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID elem      Nó 1         Nó 2'
    WRITE(5,*) '-------------------------------------'
    DO i=1, Self%N_Element
        I1 = Get_incidencia(self%Element(i),1)
        I2 = Get_incidencia(self%Element(i),2)
        WRITE(5,150) i, I1, I2
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Atribuir_Propriedades'
    WRITE(5,*) '--------------------------------------'
    WRITE(5,*) 'ID ELEM     Material  Seção   '
    WRITE(5,*) '--------------------------------------'
    DO i=1, Self%N_Element
        ID = GET_Element_Secao_ID(Self%Element(i))
        WRITE(5,155) i, ID   
155     FORMAT(I7.0,'    ','1','  ',I7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Material'
    WRITE(5,*) '--------------------------------------------'
    WRITE(5,*) 'ID  Mód. de elasticidade		  Densidade   '
    WRITE(5,*) '--------------------------------------------'
    DO i=1, Self%N_Materiais
        ID = GET_material_ID(Self%Material(i))
        Young = GET_material_Young(Self%Material(i))
        Densidade = GET_material_Densidade(Self%Material(i))
        WRITE(5,200) ID, Young, Densidade
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Secao'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID			 Área em m²'
    WRITE(5,*) '-------------------------------------'
    DO i=1, gen_size              !Gen modificator
        ID = i
        Area = Gen(i)
        WRITE(5,250) ID, Area
250     FORMAT(I7.0,' ',ES15.7e3)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Cond_Contorno'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z   (1=restrito, 0=livre)'
    WRITE(5,*) '------------------------------------'
125 FORMAT(I7,' ',I1,' ','0',' ',I1,' ','0',' ',I1,' ','0')
    DO i=1, Self%N_Cond_Contorno
        j = GET_Cond_Cont(Self%Cond_Contorno(i),1)
        k = GET_Cond_Cont(Self%Cond_Contorno(i),2)
        m = GET_Cond_Cont(Self%Cond_Contorno(i),3)
        WRITE(5,125) i, j, k, m
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Forcas'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z'
    WRITE(5,*) '------------------------------------'
    DO i=1,Self%N_Forcas
        j = GET_Fext_ID(Self%Forcas(i))
        X1 = GET_Fext(Self%Forcas(i),1)
        Y1 = GET_Fext(Self%Forcas(i),2)
        Z1 = GET_Fext(Self%Forcas(i),3)
        WRITE(5,350) j, X1, Y1, Z1
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Tempo '
    WRITE(5,*) '--------------------------------------------------------'
    WRITE(5,*) 'Dt (incremento)				N_Dt (Número de incrementos)'
    WRITE(5,*) '--------------------------------------------------------'
	WRITE(5,*) Self%DT, Self%N_Dt
    WRITE(5,*) '========================================================'
    WRITE(5,*)
    WRITE(5,*) 'Fim_do_Arquivo'
    Close(5)
    
    
    END SUBROUTINE Conversor_de_gene 
    
    
!----------------------------------------------------------------------------------------------   
    SUBROUTINE Conversor_de_gene_discreto(Self,Output,Gen)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*90                        ::  Output
    INTEGER*8                           ::  i
    INTEGER*8                           ::  j
    INTEGER*8                           ::  k
    INTEGER*8                           ::  m
    INTEGER*8                           ::  n
    INTEGER*8                           ::  gen_size
    INTEGER*8                           ::  ID
    INTEGER*8                           ::  I1
    INTEGER*8                           ::  I2
    INTEGER*8                           ::  element_total
    INTEGER*8                           ::  node_total
    REAL*8                              ::  X1
    REAL*8                              ::  X2
    REAL*8                              ::  Y1
    REAL*8                              ::  Y2
    REAL*8                              ::  Z1
    REAL*8                              ::  Z2
    REAL*8                              ::  Area
    REAL*8                              ::  Young
    REAL*8                              ::  Densidade
    REAL*8                              ::  Max_node(3)
    INTEGER*8                           ::  count_x
    INTEGER*8                           ::  count_y
    INTEGER*8                           ::  count_z
    INTEGER*8                           ::  count_x0
    INTEGER*8                           ::  count_y0
    INTEGER*8                           ::  count_z0
    INTEGER*8                           ::  count_aresta_x
    INTEGER*8                           ::  count_quina
    INTEGER*8                           ::  count_cell
    INTEGER*8, ALLOCATABLE              ::  Plano_x0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_y0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_z0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_xf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_yf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_zf(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_x(:)
    INTEGER*8, ALLOCATABLE              ::  quina(:)
    INTEGER*8, ALLOCATABLE              ::  Cell(:)
    REAL*8,    ALLOCATABLE              ::  Nodes(:,:)
    REAL*8,    ALLOCATABLE              ::  Nodes2(:,:)
    REAL*8,    ALLOCATABLE              ::  Elements(:,:)
    REAL*8                              ::  Gen(:)

    
    gen_size = SIZE(Gen)
    
    !Gravando no arquivo de saída
    
    CALL Temp_open_file3(Output) !Unit 5
100 FORMAT(I7.0,' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
150 FORMAT(I7.0,' ',I7.0,' ',I7.0)
200 FORMAT(I7.0,' ',ES15.7e3,' ',ES15.7e3)
350 FORMAT(I7.0,' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
400 FORMAT(ES15.7e3,' ',ES15.7e3)
    
    WRITE(5,*) 'Tipo_Estrutura'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS OU TRUSS'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS'
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Coordenadas'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID		x1		   y1		  z1'
    WRITE(5,*) '-------------------------------------'
    DO i=1, Self%N_Nodes
        X1 = Get_x0(Self%Node(i))
        Y1 = Get_y0(Self%Node(i))
        Z1 = Get_z0(Self%Node(i))
        WRITE(5,100) i, X1, Y1, Z1
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Incidencia'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID elem      Nó 1         Nó 2'
    WRITE(5,*) '-------------------------------------'
    DO i=1, Self%N_Element
        I1 = Get_incidencia(self%Element(i),1)
        I2 = Get_incidencia(self%Element(i),2)
        WRITE(5,150) i, I1, I2
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Atribuir_Propriedades'
    WRITE(5,*) '--------------------------------------'
    WRITE(5,*) 'ID ELEM     Material  Seção   '
    WRITE(5,*) '--------------------------------------'
    DO i=1, Self%N_Element
        ID = GET_Element_Secao_ID(Self%Element(i))
        WRITE(5,155) i, Gen(i)
155     FORMAT(I7.0,'    ','1','  ',F7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Material'
    WRITE(5,*) '--------------------------------------------'
    WRITE(5,*) 'ID  Mód. de elasticidade		  Densidade   '
    WRITE(5,*) '--------------------------------------------'
    DO i=1, Self%N_Materiais
        ID = GET_material_ID(Self%Material(i))
        Young = GET_material_Young(Self%Material(i))
        Densidade = GET_material_Densidade(Self%Material(i))
        WRITE(5,200) ID, Young, Densidade
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Secao'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID			 Área em m²'
    WRITE(5,*) '-------------------------------------'
    DO i=1, Self%N_Secao             
        Area = GET_Secao_area(Self%Secao(i))
        WRITE(5,250) i, Area
250     FORMAT(I7.0,' ',ES15.7E1)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Cond_Contorno'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z   (1=restrito, 0=livre)'
    WRITE(5,*) '------------------------------------'
125 FORMAT(I7,' ',I1,' ',I1,' ',I1)
    DO i=1, Self%N_Cond_Contorno
        j = GET_Cond_Cont(Self%Cond_Contorno(i),1)
        k = GET_Cond_Cont(Self%Cond_Contorno(i),2)
        m = GET_Cond_Cont(Self%Cond_Contorno(i),3)
        WRITE(5,125) i, j, k, m
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Forcas'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z'
    WRITE(5,*) '------------------------------------'
    DO i=1,Self%N_Forcas
        j = GET_Fext_ID(Self%Forcas(i))
        X1 = GET_Fext(Self%Forcas(i),1)
        Y1 = GET_Fext(Self%Forcas(i),2)
        Z1 = GET_Fext(Self%Forcas(i),3)
        WRITE(5,350) j, X1, Y1, Z1
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Tempo '
    WRITE(5,*) '--------------------------------------------------------'
    WRITE(5,*) 'Dt (incremento)				N_Dt (Número de incrementos)'
    WRITE(5,*) '--------------------------------------------------------'
	WRITE(5,*) Self%DT, Self%N_Dt
    WRITE(5,*) '========================================================'
    WRITE(5,*)
    WRITE(5,*) 'Fim_do_Arquivo'
    Close(5)
    
    
    END SUBROUTINE Conversor_de_gene_discreto
!----------------------------------------------------------------------------------------------   
    SUBROUTINE Conversor_de_gene_discreto_25(Self,Output,Gen)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*90                        ::  Output
    INTEGER*8                           ::  i
    INTEGER*8                           ::  j
    INTEGER*8                           ::  k
    INTEGER*8                           ::  m
    INTEGER*8                           ::  n
    INTEGER*8                           ::  gen_size
    INTEGER*8                           ::  ID
    INTEGER*8                           ::  I1
    INTEGER*8                           ::  I2
    INTEGER*8                           ::  element_total
    INTEGER*8                           ::  node_total
    REAL*8                              ::  X1
    REAL*8                              ::  X2
    REAL*8                              ::  Y1
    REAL*8                              ::  Y2
    REAL*8                              ::  Z1
    REAL*8                              ::  Z2
    REAL*8                              ::  Area
    REAL*8                              ::  Young
    REAL*8                              ::  Densidade
    REAL*8                              ::  Max_node(3)
    INTEGER*8                           ::  count_x
    INTEGER*8                           ::  count_y
    INTEGER*8                           ::  count_z
    INTEGER*8                           ::  count_x0
    INTEGER*8                           ::  count_y0
    INTEGER*8                           ::  count_z0
    INTEGER*8                           ::  count_aresta_x
    INTEGER*8                           ::  count_quina
    INTEGER*8                           ::  count_cell
    INTEGER*8, ALLOCATABLE              ::  Plano_x0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_y0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_z0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_xf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_yf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_zf(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_x(:)
    INTEGER*8, ALLOCATABLE              ::  quina(:)
    INTEGER*8, ALLOCATABLE              ::  Cell(:)
    REAL*8,    ALLOCATABLE              ::  Nodes(:,:)
    REAL*8,    ALLOCATABLE              ::  Nodes2(:,:)
    REAL*8,    ALLOCATABLE              ::  Elements(:,:)
    REAL*8                              ::  Gen(:)

    
    gen_size = SIZE(Gen)
    
    !Gravando no arquivo de saída
    
    CALL Temp_open_file3(Output) !Unit 5
100 FORMAT(I7.0,' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
150 FORMAT(I7.0,' ',I7.0,' ',I7.0)
200 FORMAT(I7.0,' ',ES15.7e3,' ',ES15.7e3)
350 FORMAT(I7.0,' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
400 FORMAT(ES15.7e3,' ',ES15.7e3)
    
    WRITE(5,*) 'Tipo_Estrutura'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS OU TRUSS'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS'
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Coordenadas'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID		x1		   y1		  z1'
    WRITE(5,*) '-------------------------------------'
    DO i=1, Self%N_Nodes
        X1 = Get_x0(Self%Node(i))
        Y1 = Get_y0(Self%Node(i))
        Z1 = Get_z0(Self%Node(i))
        WRITE(5,100) i, X1, Y1, Z1
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Incidencia'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID elem      Nó 1         Nó 2'
    WRITE(5,*) '-------------------------------------'
    DO i=1, Self%N_Element
        I1 = Get_incidencia(self%Element(i),1)
        I2 = Get_incidencia(self%Element(i),2)
        WRITE(5,150) i, I1, I2
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Atribuir_Propriedades'
    WRITE(5,*) '--------------------------------------'
    WRITE(5,*) 'ID ELEM     Material  Seção   '
    WRITE(5,*) '--------------------------------------'
    DO i=1, Self%N_Element
        ID = GET_Element_Secao_ID(Self%Element(i))
        WRITE(5,155) i, ID   !Gen modificator
155     FORMAT(I7.0,'    ','1','  ',I7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Material'
    WRITE(5,*) '--------------------------------------------'
    WRITE(5,*) 'ID  Mód. de elasticidade		  Densidade   '
    WRITE(5,*) '--------------------------------------------'
    DO i=1, Self%N_Materiais
        ID = GET_material_ID(Self%Material(i))
        Young = GET_material_Young(Self%Material(i))
        Densidade = GET_material_Densidade(Self%Material(i))
        WRITE(5,200) ID, Young, Densidade
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Secao'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID			 Área em m²'
    WRITE(5,*) '-------------------------------------'
    DO i=1, gen_size              !Gen modificator
        ID = Gen(i)
        Area = GET_secao_area(Self%Secao(ID))
        WRITE(5,250) i, Area
250     FORMAT(I7.0,' ',F9.3)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Cond_Contorno'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z   (1=restrito, 0=livre)'
    WRITE(5,*) '------------------------------------'
125 FORMAT(I7,' ',I1,' ',I1,' ',I1)
    DO i=1, Self%N_Cond_Contorno
        j = GET_Cond_Cont(Self%Cond_Contorno(i),1)
        k = GET_Cond_Cont(Self%Cond_Contorno(i),2)
        m = GET_Cond_Cont(Self%Cond_Contorno(i),3)
        WRITE(5,125) i, j, k, m
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Forcas'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z'
    WRITE(5,*) '------------------------------------'
    DO i=1,Self%N_Forcas
        j = GET_Fext_ID(Self%Forcas(i))
        X1 = GET_Fext(Self%Forcas(i),1)
        Y1 = GET_Fext(Self%Forcas(i),2)
        Z1 = GET_Fext(Self%Forcas(i),3)
        WRITE(5,350) j, X1, Y1, Z1
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Tempo '
    WRITE(5,*) '--------------------------------------------------------'
    WRITE(5,*) 'Dt (incremento)				N_Dt (Número de incrementos)'
    WRITE(5,*) '--------------------------------------------------------'
	WRITE(5,*) Self%DT, Self%N_Dt
    WRITE(5,*) '========================================================'
    WRITE(5,*)
    WRITE(5,*) 'Fim_do_Arquivo'
    Close(5)
    
    
    END SUBROUTINE Conversor_de_gene_discreto_25
    !----------------------------------------------------------------------------------------------   
    SUBROUTINE Cell_creator_Mirror(Self,Output,gen_size)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*90                        ::  Output
    INTEGER*8                           ::  Nx
    INTEGER*8                           ::  Ny
    INTEGER*8                           ::  Nz
    INTEGER*8                           ::  i
    INTEGER*8                           ::  j
    INTEGER*8                           ::  k
    INTEGER*8                           ::  m
    INTEGER*8                           ::  n
    INTEGER*8                           ::  gen_size
    INTEGER*8                           ::  I1
    INTEGER*8                           ::  I2
    REAL*8                              ::  ID
    INTEGER*8                           ::  N_Cell
    INTEGER*8                           ::  element_total
    INTEGER*8                           ::  node_total
    REAL*8                              ::  X1
    REAL*8                              ::  X2
    REAL*8                              ::  Y1
    REAL*8                              ::  Y2
    REAL*8                              ::  Z1
    REAL*8                              ::  Z2
    REAL*8                              ::  dX1
    REAL*8                              ::  dX2
    REAL*8                              ::  dY1
    REAL*8                              ::  dY2
    REAL*8                              ::  dZ1
    REAL*8                              ::  dZ2
    REAL*8                              ::  Area
    REAL*8                              ::  Young
    REAL*8                              ::  Densidade
    REAL*8                              ::  Max_node(3)
    INTEGER*8                           ::  count_x
    INTEGER*8                           ::  count_y
    INTEGER*8                           ::  count_z
    INTEGER*8                           ::  count_x0
    INTEGER*8                           ::  count_y0
    INTEGER*8                           ::  count_z0
    INTEGER*8                           ::  count_aresta_x
    INTEGER*8                           ::  count_quina_x
    INTEGER*8                           ::  count_aresta_y
    INTEGER*8                           ::  count_quina_y
    INTEGER*8                           ::  count_aresta_z
    INTEGER*8                           ::  count_quina_z
    INTEGER*8                           ::  count_aresta_xy
    INTEGER*8                           ::  count_aresta_xz
    INTEGER*8                           ::  count_aresta_yx
    INTEGER*8                           ::  count_aresta_yz
    INTEGER*8                           ::  count_aresta_zx
    INTEGER*8                           ::  count_aresta_zy
    INTEGER*8                           ::  count_cell
    INTEGER*8, ALLOCATABLE              ::  Plano_x0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_y0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_z0(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_xf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_yf(:)
    INTEGER*8, ALLOCATABLE              ::  Plano_zf(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_x(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_y(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_z(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_xy(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_xz(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_yx(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_yz(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_zx(:)
    INTEGER*8, ALLOCATABLE              ::  aresta_zy(:)
    INTEGER*8, ALLOCATABLE              ::  quina_x(:)
    INTEGER*8, ALLOCATABLE              ::  quina_y(:)
    INTEGER*8, ALLOCATABLE              ::  quina_z(:)
    INTEGER*8, ALLOCATABLE              ::  Cell(:)
    REAL*8,    ALLOCATABLE              ::  Nodes(:,:)
    REAL*8,    ALLOCATABLE              ::  Nodes2(:,:)
    REAL*8,    ALLOCATABLE              ::  Elements(:,:)
    REAL*8,    ALLOCATABLE              ::  P(:,:,:,:)
    
    
    CALL CALC_L0_Cell(Self)
    Nx = 2
    Ny = 2
    Nz = 2
    
    !Colocando um ponto base para cada célula
    ALLOCATE (P(Nx,Ny,Nz,3))
    DO i=1, (Nx)
        DO j=1, (Ny)
            DO k=1, (Nz)
                IF (i.EQ.1) THEN
                    P(i,j,k,1) = (i-1)*Self%L(1)
                ELSE
                    P(i,j,k,1) = i*Self%L(1)
                END IF
                
                IF (j.EQ.1) THEN
                    P(i,j,k,2) = (j-1)*Self%L(2)
                ELSE
                    P(i,j,k,2) = j*Self%L(2)
                END IF
                
                IF (k.EQ.1) THEN
                    P(i,j,k,3) = (k-1)*Self%L(3)
                ELSE
                    P(i,j,k,3) = k*Self%L(3)
                END IF
            END DO
        END DO
    END DO
    
    N_Cell = (i-1)*(j-1)*(k-1) !número de células
    
    Node_total = 2*N_Cell*Self%N_Nodes
    
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    count_cell = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ( ((X1.EQ.0.0d0) .AND. (X2.EQ.0.0d0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2)) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) ) THEN
            count_x0 = count_x0 + 1
        ELSE IF ( ((Y1.EQ.0.0d0) .AND. (Y2.EQ.0.0d0)) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_y0 = count_y0 + 1
        ELSE IF ( ((Z1.EQ.0.0d0) .AND. (Z2.EQ.0.0d0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) ) THEN
            count_z0 = count_z0 + 1
        ELSE IF ( ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_cell = count_cell + 1
        END IF
    END DO
    
    count_x = 0
    count_y = 0
    count_z = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF (((X1.EQ.Self%L(1)) .AND. (X2.EQ.Self%L(1))) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2)) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) ) THEN
            count_x = count_x + 1
        END IF
        IF (((Y1.EQ.Self%L(2)) .AND. (Y2.EQ.Self%L(2))) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_y = count_y + 1
        END IF
        IF (((Z1.EQ.Self%L(3)) .AND. (Z2.EQ.Self%L(3))) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) ) THEN
            count_z = count_z + 1
        END IF
    END DO
    
    count_aresta_x = 0
    count_aresta_y = 0
    count_aresta_z = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ( ((Y1.EQ.Self%L(2)) .AND. (Y2.EQ.Self%L(2))) .AND. ((Z1.EQ.Self%L(3)) .AND. (Z2.EQ.Self%L(3))) ) THEN
            count_aresta_x = count_aresta_x + 1
        END IF
        IF ( ((X1.EQ.Self%L(1)) .AND. (X2.EQ.Self%L(1))) .AND. ((Z1.EQ.Self%L(3)) .AND. (Z2.EQ.Self%L(3))) ) THEN
            count_aresta_y = count_aresta_y + 1
        END IF
        IF ( ((Y1.EQ.Self%L(2)) .AND. (Y2.EQ.Self%L(2))) .AND. ((X1.EQ.Self%L(1)) .AND. (X2.EQ.Self%L(1))) ) THEN
            count_aresta_z = count_aresta_z + 1
        END IF
    END DO
    
    !Salvando ID de cada elemento no seu vetor
    ALLOCATE (Plano_x0(count_x0), Plano_y0(count_y0), Plano_z0(count_z0), Cell(count_cell))
    ALLOCATE (Plano_xf(count_x), Plano_yf(count_y), Plano_zf(count_z))
    ALLOCATE (Aresta_x(count_aresta_x), Aresta_y(count_aresta_y), Aresta_z(count_aresta_z))
    
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    count_cell = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ( ((X1.EQ.0.0d0) .AND. (X2.EQ.0.0d0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_x0 = count_x0 + 1
            Plano_x0(count_x0) = GET_Element_ID(Self%Element(i))
        ELSE IF ( ((Y1.EQ.0.0d0) .AND. (Y2.EQ.0.0d0)) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_y0 = count_y0 + 1
            Plano_y0(count_y0) = GET_Element_ID(Self%Element(i))
        ELSE IF ( ((Z1.EQ.0.0d0) .AND. (Z2.EQ.0.0d0)) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) ) THEN
            count_z0 = count_z0 + 1
            Plano_z0(count_z0) = GET_Element_ID(Self%Element(i))
        ELSE IF ( ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_cell = count_cell + 1
            Cell(count_cell) = GET_Element_ID(Self%Element(i))
        END IF
    END DO
    

    
    
    

    count_x = 0
    count_y = 0
    count_z = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF (((X1.EQ.Self%L(1)) .AND. (X2.EQ.Self%L(1))) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3))) ) THEN
            count_x = count_x + 1
            Plano_xf(count_x) = GET_Element_ID(Self%Element(i))
        END IF
        IF (((Y1.EQ.Self%L(2)) .AND. (Y2.EQ.Self%L(2))) .AND. ((X1.NE.Self%L(1)) .OR. (X2.NE.Self%L(1))) .AND. ((Z1.NE.Self%L(3)) .OR. (Z2.NE.Self%L(3)))) THEN
            count_y = count_y + 1
            Plano_yf(count_y) = GET_Element_ID(Self%Element(i))
        END IF
        IF (((Z1.EQ.Self%L(3)) .AND. (Z2.EQ.Self%L(3))) .AND. ((Y1.NE.Self%L(2)) .OR. (Y2.NE.Self%L(2))) .AND. ((X1.NE.Self%L(3)) .OR. (X2.NE.Self%L(3))) ) THEN
            count_z = count_z + 1
            Plano_zf(count_z) = GET_Element_ID(Self%Element(i))
        END IF
    END DO
    
    count_aresta_x = 0
    count_aresta_y = 0
    count_aresta_z = 0
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ( ((Y1.EQ.Self%L(2)) .AND. (Y2.EQ.Self%L(2))) .AND. ((Z1.EQ.Self%L(3)) .AND. (Z2.EQ.Self%L(3))) ) THEN
            count_aresta_x = count_aresta_x + 1
            aresta_x(count_aresta_x) = i
        END IF
        IF ( ((X1.EQ.Self%L(1)) .AND. (X2.EQ.Self%L(1))) .AND. ((Z1.EQ.Self%L(3)) .AND. (Z2.EQ.Self%L(3))) ) THEN
            count_aresta_y = count_aresta_y + 1
            aresta_y(count_aresta_y) = i
        END IF
        IF ( ((Y1.EQ.Self%L(2)) .AND. (Y2.EQ.Self%L(2))) .AND. ((X1.EQ.Self%L(1)) .AND. (X2.EQ.Self%L(1))) ) THEN
            count_aresta_z = count_aresta_z + 1
            aresta_z(count_aresta_z) = i
        END IF
    END DO
    
    
    
    
    !Número total de elementos depois da repetição
    Element_total = N_cell*(count_cell + count_x0 + count_y0 + count_z0) + 4*(count_x + count_y + count_z) + 2*(count_aresta_x + count_aresta_y + count_aresta_z)
    ALLOCATE(Elements(Element_total,11))

    !   1 - ID elemento novo
    !   2 - X nó 1
    !   3 - Y nó 1
    !   4 - Z nó 1
    !   5 - X nó 2
    !   6 - Y nó 2
    !   7 - Z nó 2
    !   8 - ID Seção
    !   9 - ID Material
    !   10 - Incidência nova 1
    !   11 - Incidência nova 2
    
    Elements = 0.0d0
    
    
    !repetindo elementos
    n = 1
    DO i=1, 2
        DO j=1, 2
            DO k=1, 2                
                        
                    DO m=1, count_x0    !Repetindo plano x0
                        ID = Plano_x0(m)
                        
                        Elements(n,1) = n
                        
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        IF (i.EQ.1) THEN
                            X1 = Get_x0(Self%Node(I1))
                            X2 = Get_x0(Self%Node(I2))
                        ELSE
                            X1 = -1 * (Get_x0(Self%Node(I1)))
                            X2 = -1 * (Get_x0(Self%Node(I2)))
                        END IF
                        
                        IF (j.EQ.1) THEN
                            Y1 = Get_y0(Self%Node(I1))
                            Y2 = Get_y0(Self%Node(I2))
                        ELSE
                            Y1 = -1 * (Get_y0(Self%Node(I1)))
                            Y2 = -1 * (Get_y0(Self%Node(I2)))
                        END IF
                        
                        IF (k.EQ.1) THEN
                            Z1 = Get_z0(Self%Node(I1))
                            Z2 = Get_z0(Self%Node(I2))
                        ELSE
                            Z1 = -1 * (Get_z0(Self%Node(I1)))
                            Z2 = -1 * (Get_z0(Self%Node(I2)))
                        END IF
                        
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                        
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                        
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                        
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                        n = n + 1
                    END DO
                IF (i.EQ.1) THEN
                    DO m=1, count_x
                        ID = Plano_xf(m)
                        
                        Elements(n,1) = n
                        
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        IF (j.EQ.1) THEN
                            Y1 = Get_y0(Self%Node(I1))
                            Y2 = Get_y0(Self%Node(I2))
                        ELSE
                            Y1 = -1 * (Get_y0(Self%Node(I1)))
                            Y2 = -1 * (Get_y0(Self%Node(I2)))
                        END IF
                        
                        IF (k.EQ.1) THEN
                            Z1 = Get_z0(Self%Node(I1))
                            Z2 = Get_z0(Self%Node(I2))
                        ELSE
                            Z1 = -1 * (Get_z0(Self%Node(I1)))
                            Z2 = -1 * (Get_z0(Self%Node(I2)))
                        END IF
                    
                        Elements(n,2) = self%L(1)
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = self%L(1)
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                        
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                        n = n + 1
                    END DO
                END IF
                    
        
                    DO m=1, count_y0    !Repetindo plano y0
                        ID = Plano_y0(m)
                    
                        Elements(n,1) = n
                    
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        IF (i.EQ.1) THEN
                            X1 = Get_x0(Self%Node(I1))
                            X2 = Get_x0(Self%Node(I2))
                        ELSE
                            X1 = -1 * (Get_x0(Self%Node(I1)))
                            X2 = -1 * (Get_x0(Self%Node(I2)))
                        END IF
                        
                        IF (j.EQ.1) THEN
                            Y1 = Get_y0(Self%Node(I1))
                            Y2 = Get_y0(Self%Node(I2))
                        ELSE
                            Y1 = -1 * (Get_y0(Self%Node(I1)))
                            Y2 = -1 * (Get_y0(Self%Node(I2)))
                        END IF
                        
                        IF (k.EQ.1) THEN
                            Z1 = Get_z0(Self%Node(I1))
                            Z2 = Get_z0(Self%Node(I2))
                        ELSE
                            Z1 = -1 * (Get_z0(Self%Node(I1)))
                            Z2 = -1 * (Get_z0(Self%Node(I2)))
                        END IF
                        
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                    
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))

                        n = n + 1
                    END DO
                IF (j.EQ.1) THEN
                    DO m=1, count_y        
                        ID = Plano_yf(m)
                        Elements(n,1) = n
                        
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        IF (i.EQ.1) THEN
                            X1 = Get_x0(Self%Node(I1))
                            X2 = Get_x0(Self%Node(I2))
                        ELSE
                            X1 = -1 * (Get_x0(Self%Node(I1)))
                            X2 = -1 * (Get_x0(Self%Node(I2)))
                        END IF
                                                
                        IF (k.EQ.1) THEN
                            Z1 = Get_z0(Self%Node(I1))
                            Z2 = Get_z0(Self%Node(I2))
                        ELSE
                            Z1 = -1 * (Get_z0(Self%Node(I1)))
                            Z2 = -1 * (Get_z0(Self%Node(I2)))
                        END IF
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = Self%L(2)
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = Self%L(2)
                        Elements(n,7) = P(i,j,k,3) + Z2
                    
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))

                        n = n + 1
                    END DO
                END IF
                
                
                
                    DO m=1, count_z0    !Repetindo plano z0
                        ID = Plano_z0(m)
                    
                        Elements(n,1) = n
                    
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        IF (i.EQ.1) THEN
                            X1 = Get_x0(Self%Node(I1))
                            X2 = Get_x0(Self%Node(I2))
                        ELSE
                            X1 = -1 * (Get_x0(Self%Node(I1)))
                            X2 = -1 * (Get_x0(Self%Node(I2)))
                        END IF
                        
                        IF (j.EQ.1) THEN
                            Y1 = Get_y0(Self%Node(I1))
                            Y2 = Get_y0(Self%Node(I2))
                        ELSE
                            Y1 = -1 * (Get_y0(Self%Node(I1)))
                            Y2 = -1 * (Get_y0(Self%Node(I2)))
                        END IF
                        
                        IF (k.EQ.1) THEN
                            Z1 = Get_z0(Self%Node(I1))
                            Z2 = Get_z0(Self%Node(I2))
                        ELSE
                            Z1 = -1 * (Get_z0(Self%Node(I1)))
                            Z2 = -1 * (Get_z0(Self%Node(I2)))
                        END IF
                        
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = P(i,j,k,3) + Z1
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = P(i,j,k,3) + Z2
                    
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))

                        n = n + 1
                    END DO
                IF (k.EQ.1) THEN
                    DO m=1, count_z
                        ID = Plano_zf(m)
                        Elements(n,1) = n
                        
                        I1=GET_Incidencia(Self%Element(ID),1)
                        I2=GET_Incidencia(Self%Element(ID),2) 
                        IF (i.EQ.1) THEN
                            X1 = Get_x0(Self%Node(I1))
                            X2 = Get_x0(Self%Node(I2))
                        ELSE
                            X1 = -1 * (Get_x0(Self%Node(I1)))
                            X2 = -1 * (Get_x0(Self%Node(I2)))
                        END IF
                        
                        IF (j.EQ.1) THEN
                            Y1 = Get_y0(Self%Node(I1))
                            Y2 = Get_y0(Self%Node(I2))
                        ELSE
                            Y1 = -1 * (Get_y0(Self%Node(I1)))
                            Y2 = -1 * (Get_y0(Self%Node(I2)))
                        END IF
                        
                    
                        Elements(n,2) = P(i,j,k,1) + X1
                        Elements(n,3) = P(i,j,k,2) + Y1
                        Elements(n,4) = Self%L(3)
                    
                        Elements(n,5) = P(i,j,k,1) + X2
                        Elements(n,6) = P(i,j,k,2) + Y2
                        Elements(n,7) = Self%L(3)
                        
                        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
                        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))

                        n = n + 1
                    END DO
                END IF
                
                
                
                
                DO m=1, count_cell    !Repetindo interiores
                    ID = Cell(m)
                    
                    Elements(n,1) = n
                    
                    I1=GET_Incidencia(Self%Element(ID),1)
                    I2=GET_Incidencia(Self%Element(ID),2) 
                    IF (i.EQ.1) THEN
                        X1 = Get_x0(Self%Node(I1))
                        X2 = Get_x0(Self%Node(I2))
                    ELSE
                        X1 = -1 * (Get_x0(Self%Node(I1)))
                        X2 = -1 * (Get_x0(Self%Node(I2)))
                    END IF
                        
                    IF (j.EQ.1) THEN
                        Y1 = Get_y0(Self%Node(I1))
                        Y2 = Get_y0(Self%Node(I2))
                    ELSE
                        Y1 = -1 * (Get_y0(Self%Node(I1)))
                        Y2 = -1 * (Get_y0(Self%Node(I2)))
                    END IF
                        
                    IF (k.EQ.1) THEN
                        Z1 = Get_z0(Self%Node(I1))
                        Z2 = Get_z0(Self%Node(I2))
                    ELSE
                        Z1 = -1 * (Get_z0(Self%Node(I1)))
                        Z2 = -1 * (Get_z0(Self%Node(I2)))
                    END IF
                        
                    
                    Elements(n,2) = P(i,j,k,1) + X1
                    Elements(n,3) = P(i,j,k,2) + Y1
                    Elements(n,4) = P(i,j,k,3) + Z1
                    
                    Elements(n,5) = P(i,j,k,1) + X2
                    Elements(n,6) = P(i,j,k,2) + Y2
                    Elements(n,7) = P(i,j,k,3) + Z2
                    
                    Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
                    Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
                    n = n + 1
                END DO
    
            END DO
        END DO
    END DO
    
    !Repetindo arestas
    
    !Aresta x
    DO m=1, count_aresta_x
        ID = aresta_x(m)
                    
        Elements(n,1) = n
                    
        I1=GET_Incidencia(Self%Element(ID),1)
        I2=GET_Incidencia(Self%Element(ID),2) 
        X1 = (Get_x0(Self%Node(I1)))
        X2 = (Get_x0(Self%Node(I2)))
        
        Elements(n,2) = X1
        Elements(n,3) = Self%L(2)
        Elements(n,4) = Self%L(3)
                    
        Elements(n,5) = X2
        Elements(n,6) = Self%L(2)
        Elements(n,7) = Self%L(3)
                    
        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
        
        n = n + 1
        
        Elements(n,1) = n
                    
        I1=GET_Incidencia(Self%Element(ID),1)
        I2=GET_Incidencia(Self%Element(ID),2) 
        X1 = -1 * (Get_x0(Self%Node(I1)))
        X2 = -1 * (Get_x0(Self%Node(I2)))
        
        Elements(n,2) = 2*Self%L(1) + X1
        Elements(n,3) = Self%L(2)
        Elements(n,4) = Self%L(3)
                    
        Elements(n,5) = 2*Self%L(1) + X2
        Elements(n,6) = Self%L(2)
        Elements(n,7) = Self%L(3)
                    
        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
        
        n = n + 1
    END DO
    
    !Aresta y
    DO m=1, count_aresta_y
        ID = aresta_y(m)
                    
        Elements(n,1) = n
                    
        I1=GET_Incidencia(Self%Element(ID),1)
        I2=GET_Incidencia(Self%Element(ID),2) 
        Y1 = (Get_y0(Self%Node(I1)))
        Y2 = (Get_y0(Self%Node(I2)))
        
        Elements(n,2) = Self%L(1)
        Elements(n,3) = Y1
        Elements(n,4) = Self%L(3)
                    
        Elements(n,5) = Self%L(1)
        Elements(n,6) = Y2
        Elements(n,7) = Self%L(3)
                    
        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
        
        n = n + 1
        
        Elements(n,1) = n
                    
        I1=GET_Incidencia(Self%Element(ID),1)
        I2=GET_Incidencia(Self%Element(ID),2) 
        Y1 = -1 * (Get_y0(Self%Node(I1)))
        Y2 = -1 * (Get_y0(Self%Node(I2)))
        
        Elements(n,2) = Self%L(1)
        Elements(n,3) = 2*Self%L(2) + Y1
        Elements(n,4) = Self%L(3)
                    
        Elements(n,5) = Self%L(1)
        Elements(n,6) = 2*Self%L(2) + Y2
        Elements(n,7) = Self%L(3)
                    
        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
        
        n = n + 1
    END DO
    
    !Aresta z
    DO m=1, count_aresta_z
        ID = aresta_z(m)
                    
        Elements(n,1) = n
                    
        I1=GET_Incidencia(Self%Element(ID),1)
        I2=GET_Incidencia(Self%Element(ID),2) 
        Z1 = (Get_z0(Self%Node(I1)))
        Z2 = (Get_z0(Self%Node(I2)))
        
        Elements(n,2) = Self%L(1)
        Elements(n,3) = Self%L(2)
        Elements(n,4) = Z1
                    
        Elements(n,5) = Self%L(1)
        Elements(n,6) = Self%L(2)
        Elements(n,7) = Z2
                    
        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
        
        n = n + 1
        
        
        Elements(n,1) = n
                    
        I1=GET_Incidencia(Self%Element(ID),1)
        I2=GET_Incidencia(Self%Element(ID),2) 
        Z1 = -1 * (Get_z0(Self%Node(I1)))
        Z2 = -1 * (Get_z0(Self%Node(I2)))
        
        Elements(n,2) = Self%L(1)
        Elements(n,3) = Self%L(2)
        Elements(n,4) = 2*Self%L(3) + Z1
                    
        Elements(n,5) = Self%L(1)
        Elements(n,6) = Self%L(2)
        Elements(n,7) = 2*Self%L(3) + Z2
                    
        Elements(n,8) = GET_Element_Secao_ID(Self%Element(ID))
                    
        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
        
        n = n + 1
    END DO
    
    
    
    ALLOCATE (Nodes2(node_total,4))
    Nodes2 = 0.0d0

    n = n-1
    j = 1
    i = 1

    dX1=Self%L(1)/10    
    dY1=Self%L(2)/10
    dZ1=Self%L(3)/10
    
    
    DO j=1, element_total
            
            
        X1 = Elements(j,2) 
        Y1 = Elements(j,3)
        Z1 = Elements(j,4)

        DO k=1, (i-1)       !Teste se nó já existe
            X2 = abs(Nodes2(k,2) - X1)
            Y2 = abs(Nodes2(k,3) - Y1)
            Z2 = abs(Nodes2(k,4) - Z1)
            IF ((X2.LE.dX1) .AND. (Y2.LE.dY1) .AND. (Z2.LE.dZ1) ) THEN
                Elements(j,10) = k
                GO TO 115
            END IF
        END DO
        Nodes2(i,1) = i
        Nodes2(i,2) = X1
        Nodes2(i,3) = Y1
        Nodes2(i,4) = Z1
        Elements(j,10) = i
        i = i + 1
            
115     CONTINUE
            
        X1 = Elements(j,5)
        Y1 = Elements(j,6)
        Z1 = Elements(j,7)

        DO k=1, (i-1)       !Teste se nó já existe
            X2 = abs(Nodes2(k,2) - X1)
            Y2 = abs(Nodes2(k,3) - Y1)
            Z2 = abs(Nodes2(k,4) - Z1)    
            IF ((X2.LE.dX1) .AND. (Y2.LE.dY1) .AND. (Z2.LE.dZ1) ) THEN
                Elements(j,11) = k
                GO TO 116
            END IF
        END DO
        Nodes2(i,1) = i
        Nodes2(i,2) = X1
        Nodes2(i,3) = Y1
        Nodes2(i,4) = Z1
        Elements(j,11) = i
        i = i + 1
116     CONTINUE
        Print *,i,'/',element_total
    END DO
    
    node_total = i - 1
    ALLOCATE (Nodes(i,4))
    DO j=1, node_total
        Nodes(j,:) = Nodes2(j,:)
    END DO
    DEALLOCATE (Nodes2)
    
    

   
    !Criando forças
    !Definindo o final em X, Y, Z
    Max_node(:) = P(Nx,Ny,Nz,:) + Self%L(:)
    
    !Contando quantos nós estão em cada final
    count_x = 0
    count_y = 0
    count_z = 0
    count_x0 = 0
    count_y0 = 0
    count_z0 = 0
    count_aresta_x = 0
    count_quina_x = 0
    count_aresta_y = 0
    count_quina_y = 0
    count_aresta_z = 0
    count_quina_z = 0
    count_aresta_xy = 0
    count_aresta_xz = 0
    count_aresta_yx = 0
    count_aresta_yz = 0
    count_aresta_zx = 0
    count_aresta_zy = 0
    DO i=1, node_total
        X1 = Nodes(i,2) 
        Y1 = Nodes(i,3)
        Z1 = Nodes(i,4)
        X2 = Max_node(1)-X1
        Y2 = Max_node(2)-Y1
        Z2 = Max_node(3)-Z1
        IF (X2.EQ.0.0d0) THEN
            !aresta da força cisalhante
            IF (Y2.EQ.0.0d0) THEN
                count_aresta_xy = count_aresta_xy + 1
            ELSE IF (Z2.EQ.0.0d0) THEN
                count_aresta_xz = count_aresta_xz + 1
            END IF
            !aresta da força de tração/compressão
            IF (( (Y2.EQ.0.0d0).OR.(Y1.EQ.0.0d0)) .AND. ((Z2.EQ.0.0d0).OR.(Z1.EQ.0.0d0)) ) THEN
                count_quina_x = count_quina_x + 1                
            ELSE IF ((Y2.EQ.0.0d0) .OR. (Y1.EQ.0.0d0)) THEN
                count_aresta_x = count_aresta_x + 1
            ELSE IF ((Z2.EQ.0.0d0) .OR. (Z1.EQ.0.0d0)) THEN
                count_aresta_x = count_aresta_x + 1
            ELSE
                count_x = count_x + 1
            END IF
        ELSE IF (X1.EQ.0.0d0) THEN
            count_x0 = count_x0 + 1
        END IF

        IF (Y2.EQ.0.0d0) THEN
            !aresta da força cisalhante
            IF (X2.EQ.0.0d0) THEN
                count_aresta_yx = count_aresta_yx + 1
            ELSE IF (Z2.EQ.0.0d0) THEN
                count_aresta_yz = count_aresta_yz + 1
            END IF
            !aresta da força de tração/compressão
            IF ( ((X2.EQ.0.0d0).OR.(X1.EQ.0.0d0)) .AND. ((Z2.EQ.0.0d0).OR.(Z1.EQ.0.0d0)) ) THEN
                count_quina_y = count_quina_y + 1                
            ELSE IF ((X2.EQ.0.0d0) .OR. (X1.EQ.0.0d0)) THEN
                count_aresta_y = count_aresta_y + 1
            ELSE IF ((Z2.EQ.0.0d0) .OR. (Z1.EQ.0.0d0)) THEN
                count_aresta_y = count_aresta_y + 1
            ELSE
                count_y = count_y + 1
            END IF
        ELSE IF (Y1.EQ.0.0d0) THEN
            count_y0 = count_y0 + 1
        END IF
        
        IF (Z2.EQ.0.0d0) THEN
            !aresta da força cisalhante
            IF (X2.EQ.0.0d0) THEN
                count_aresta_zx = count_aresta_zx + 1
            ELSE IF (Y2.EQ.0.0d0) THEN
                count_aresta_zy = count_aresta_zy + 1
            END IF
            !aresta da força de tração/compressão
            IF ( ((Y2.EQ.0.0d0).OR.(Y1.EQ.0.0d0)) .AND. ((X2.EQ.0.0d0).OR.(X1.EQ.0.0d0)) ) THEN
                count_quina_z = count_quina_z + 1                
            ELSE IF ((Y2.EQ.0.0d0) .OR. (Y1.EQ.0.0d0)) THEN
                count_aresta_z = count_aresta_z + 1
            ELSE IF ((X2.EQ.0.0d0) .OR. (X1.EQ.0.0d0)) THEN
                count_aresta_z = count_aresta_z + 1
            ELSE
                count_z = count_z + 1
            END IF
        ELSE IF (Z1.EQ.0.0d0) THEN
            count_z0 = count_z0 + 1
        END IF
        
    END DO
    
    DEALLOCATE (Plano_x0, Plano_y0, Plano_z0, Plano_xf, Plano_yf, Plano_zf, aresta_x, aresta_y, aresta_z)
    
    
    
    !=================================================================
    
    
   !Gravando no arquivo de saída
    
    CALL Temp_open_file3(Output) !Unit 5
    
100 FORMAT(I7.0,' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
150 FORMAT(I7.0,' ',I7.0,' ',I7.0)

350 FORMAT('1',' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
400 FORMAT(ES15.7e3,' ',ES15.7e3)
    WRITE(5,*) 'Tipo_Estrutura'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS OU TRUSS'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS'
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Coordenadas'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID		x1		   y1		  z1'
    WRITE(5,*) '-------------------------------------'
    DO i=1, node_total
        WRITE(5,100) i, Nodes(i,2), Nodes(i,3), Nodes(i,4)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Atribuir_Propriedades'
    WRITE(5,*) '--------------------------------------'
    WRITE(5,*) 'ID ELEM     Material  Seção   '
    WRITE(5,*) '--------------------------------------'
    DO i=1, element_total
        WRITE(5,155) i, Elements(i,8)
155     FORMAT(I7.0,'    1   ',F7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Incidencia'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID elem      Nó 1         Nó 2'
    WRITE(5,*) '-------------------------------------'
    DO i=1, element_total
        WRITE(5,160) Elements(i,1), Elements(i,10), Elements(i,11)
160     FORMAT(F7.0,' ',F7.0,' ',F7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Secao'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID			 Área em m²'
    WRITE(5,*) '-------------------------------------'
    DO i=1, self%N_Secao              !Gen modificator
        ID = GET_secao_ID(Self%Secao(i))
        Area = GET_secao_Area(Self%Secao(i))
        WRITE(5,250) ID, Area
250     FORMAT(F7.0,' ',F7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Cond_Contorno'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z   (1=restrito, 0=livre)'
    WRITE(5,*) '------------------------------------'
    DO i=1, node_total
        WRITE(5,300) i
300     FORMAT(I7.0,' ','0',' ','0',' ','0',' ','0',' ','0',' ','0')
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Material'
    WRITE(5,*) '--------------------------------------------'
    WRITE(5,*) 'ID  Mód. de elasticidade		  Densidade   '
    WRITE(5,*) '--------------------------------------------'
    DO i=1, Self%N_Materiais
        ID = GET_material_ID(Self%Material(i))
        Young = GET_material_Young(Self%Material(i))
        Densidade = GET_material_Densidade(Self%Material(i))
        WRITE(5,200) ID, Young, Densidade
200     FORMAT(F7.0,' ',ES15.7e3,' ',ES15.7e3)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Forcas'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z'
    WRITE(5,*) '------------------------------------'
        X1 = GET_Fext(Self%Forcas(1),1)
        Y1 = GET_Fext(Self%Forcas(1),2)
        Z1 = GET_Fext(Self%Forcas(1),3)
        WRITE(5,350)  X1, Y1, Z1
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Tempo '
    WRITE(5,*) '--------------------------------------------------------'
    WRITE(5,*) 'Dt (incremento)				N_Dt (Número de incrementos)'
    WRITE(5,*) '--------------------------------------------------------'
	WRITE(5,400) Self%DT, Self%N_Dt
    WRITE(5,*) '========================================================'
    WRITE(5,*)
    WRITE(5,*) 'Fim_do_Arquivo'
    DEALLOCATE (Elements)
    Close(5)
    gen_size = self%N_secao
    END SUBROUTINE Cell_creator_Mirror
    !----------------------------------------------------------------------------------------------   
    SUBROUTINE Cell_creator_Isotropic(Self,Output,gen_size)
    IMPLICIT NONE
    TYPE(Type_Problem), INTENT(INOUT)   ::  Self
    CHARACTER*90                        ::  Output
    INTEGER*8                           ::  Nx
    INTEGER*8                           ::  Ny
    INTEGER*8                           ::  Nz
    INTEGER*8                           ::  i
    INTEGER*8                           ::  j
    INTEGER*8                           ::  k
    INTEGER*8                           ::  m
    INTEGER*8                           ::  n
    INTEGER*8                           ::  gen_size
    INTEGER*8                           ::  I1
    INTEGER*8                           ::  I2
    REAL*8                              ::  ID
    INTEGER*8                           ::  N_Cell
    INTEGER*8                           ::  element_total
    INTEGER*8                           ::  node_total
    REAL*8                              ::  X1
    REAL*8                              ::  X2
    REAL*8                              ::  Y1
    REAL*8                              ::  Y2
    REAL*8                              ::  Z1
    REAL*8                              ::  Z2
    REAL*8                              ::  dX1
    REAL*8                              ::  dX2
    REAL*8                              ::  dY1
    REAL*8                              ::  dY2
    REAL*8                              ::  dZ1
    REAL*8                              ::  dZ2
    REAL*8                              ::  Area
    REAL*8                              ::  Young
    REAL*8                              ::  Densidade

    INTEGER*8                           ::  count_xyz

    INTEGER*8                           ::  count_cell

    INTEGER*8, ALLOCATABLE              ::  Eixo_xyz(:)
    INTEGER*8, ALLOCATABLE              ::  Cell(:)
    
    REAL*8,    ALLOCATABLE              ::  Nodes(:,:)
    REAL*8,    ALLOCATABLE              ::  Nodes2(:,:)
    REAL*8,    ALLOCATABLE              ::  Elements(:,:)

    
    
    CALL CALC_L0_Cell(Self)
    
    
    
    
    N_Cell = 3 !número de células
    
    count_xyz = 0!Número de elementos na diagonal principal 
    count_cell = 0!Número de elementos a serem rotacionados
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ( ((X1.EQ.Y1) .AND. (Y1.EQ.Z1)) .AND. ((X2.EQ.Y2) .AND. (Y2.EQ.Z2) ) ) THEN
            count_xyz = count_xyz + 1
        ELSE IF ((X1.NE.Z1) .OR. (X2.NE.Z2)) THEN
            count_cell = count_cell + 1
        END IF
    END DO
    
    
    
    !Salvando ID de cada elemento no seu vetor
    ALLOCATE (Eixo_xyz(count_xyz), Cell(count_cell))
    
    
    count_xyz = 0!Número de elementos na diagonal principal 
    count_cell = 0!Número de elementos a serem rotacionados
    DO i=1, Self%N_Element
        I1=GET_Incidencia(Self%Element(i),1)
        I2=GET_Incidencia(Self%Element(i),2) 
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        IF ( ((X1.EQ.Y1) .AND. (Y1.EQ.Z1)) .AND. ((X2.EQ.Y2) .AND. (Y2.EQ.Z2) ) ) THEN
            count_xyz = count_xyz + 1
            Eixo_xyz(count_xyz) = GET_Element_ID(Self%Element(i))
        ELSE IF ((X1.NE.Z1) .OR. (X2.NE.Z2)) THEN
            count_cell = count_cell + 1
            Cell(count_cell) = GET_Element_ID(Self%Element(i))
        END IF
    END DO
    
    
    
    !Número total de elementos depois da repetição
    Element_total = N_cell*count_cell + count_xyz
    ALLOCATE(Elements(Element_total,11))

    !   1 - ID elemento novo
    !   2 - X nó 1
    !   3 - Y nó 1
    !   4 - Z nó 1
    !   5 - X nó 2
    !   6 - Y nó 2
    !   7 - Z nó 2
    !   8 - ID Seção
    !   9 - ID Material
    !   10 - Incidência nova 1
    !   11 - Incidência nova 2
    
    Elements = 0.0d0
    
    
    !repetindo elementos
    n = 1
    
    DO i=1, count_xyz    !Inserindo eixo principal
        ID = Eixo_xyz(i)
                        
        Elements(n,1) = n
                        
        I1=GET_Incidencia(Self%Element(ID),1)
        I2=GET_Incidencia(Self%Element(ID),2) 
        
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))      
        
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
        
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
        
                        
        Elements(n,2) = X1
        Elements(n,3) = Y1
        Elements(n,4) = Z1
                        
        Elements(n,5) = X2
        Elements(n,6) = Y2
        Elements(n,7) = Z2
                        
        Elements(n,8) = n
                                           
        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
        n = n + 1
    END DO
    

    DO j=1, count_cell    !Repetindo primeira célula
        ID = Cell(j)
                    
        Elements(n,1) = n
                    
        I1=GET_Incidencia(Self%Element(ID),1)
        I2=GET_Incidencia(Self%Element(ID),2) 
                    
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
                    
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
                    
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
                    
                        
                    
        Elements(n,2) = X1
        Elements(n,3) = Y1
        Elements(n,4) = Z1
                    
        Elements(n,5) = X2
        Elements(n,6) = Y2
        Elements(n,7) = Z2
                    
        Elements(n,8) = n
                    
        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
        n = n + 1
    END DO
    
    DO j=1, count_cell    !Repetindo segunda célula
        ID = Cell(j)
                    
        Elements(n,1) = n
                    
        I1=GET_Incidencia(Self%Element(ID),1)
        I2=GET_Incidencia(Self%Element(ID),2) 
                    
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
                    
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
                    
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
                    
                        
                    
        Elements(n,2) = Z1
        Elements(n,3) = X1
        Elements(n,4) = Y1
                    
        Elements(n,5) = Z2
        Elements(n,6) = X2
        Elements(n,7) = Y2
                    
        Elements(n,8) = j+count_xyz
                    
        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
        n = n + 1
    END DO
    
    DO j=1, count_cell    !Repetindo terceira célula
        ID = Cell(j)
                    
        Elements(n,1) = n
                    
        I1=GET_Incidencia(Self%Element(ID),1)
        I2=GET_Incidencia(Self%Element(ID),2) 
                    
        X1 = Get_x0(Self%Node(I1))
        X2 = Get_x0(Self%Node(I2))
                    
        Y1 = Get_y0(Self%Node(I1))
        Y2 = Get_y0(Self%Node(I2))
                    
        Z1 = Get_z0(Self%Node(I1))
        Z2 = Get_z0(Self%Node(I2))
                    
                        
                    
        Elements(n,2) = Y1
        Elements(n,3) = Z1
        Elements(n,4) = X1
                    
        Elements(n,5) = Y2
        Elements(n,6) = Z2
        Elements(n,7) = X2
                    
        Elements(n,8) = j+count_xyz
                    
        Elements(n,9) = GET_Element_Material_ID(Self%Element(ID))
        n = n + 1
    END DO

    
    
    node_total = Self%N_Nodes * N_Cell
    ALLOCATE (Nodes2(node_total,4))
    Nodes2 = 0.0d0

    n = n-1
    j = 1
    i = 1

    dX1=Self%L(1)/10
    dY1=Self%L(2)/10
    dZ1=Self%L(3)/10
        
    DO j=1, element_total
            
            
        X1 = Elements(j,2) 
        Y1 = Elements(j,3)
        Z1 = Elements(j,4)

        DO k=1, (i-1)       !Teste se nó já existe
            X2 = abs(Nodes2(k,2) - X1)
            Y2 = abs(Nodes2(k,3) - Y1)
            Z2 = abs(Nodes2(k,4) - Z1)
            IF ((X2.LE.dX1) .AND. (Y2.LE.dY1) .AND. (Z2.LE.dZ1) ) THEN
                Elements(j,10) = k
                GO TO 115
            END IF
        END DO
        Nodes2(i,1) = i
        Nodes2(i,2) = X1
        Nodes2(i,3) = Y1
        Nodes2(i,4) = Z1
        Elements(j,10) = i
        i = i + 1
            
115     CONTINUE
            
        X1 = Elements(j,5)
        Y1 = Elements(j,6)
        Z1 = Elements(j,7)

        DO k=1, (i-1)       !Teste se nó já existe
            X2 = abs(Nodes2(k,2) - X1)
            Y2 = abs(Nodes2(k,3) - Y1)
            Z2 = abs(Nodes2(k,4) - Z1)    
            IF ((X2.LE.dX1) .AND. (Y2.LE.dY1) .AND. (Z2.LE.dZ1) ) THEN
                Elements(j,11) = k
                GO TO 116
            END IF
        END DO
        Nodes2(i,1) = i
        Nodes2(i,2) = X1
        Nodes2(i,3) = Y1
        Nodes2(i,4) = Z1
        Elements(j,11) = i
        i = i + 1
116     CONTINUE
        Print *,i,'/',element_total
    END DO
    
    node_total = i - 1
    ALLOCATE (Nodes(i,4))
    DO j=1, node_total
        Nodes(j,:) = Nodes2(j,:)
    END DO
    DEALLOCATE (Nodes2)
    
    
    !=================================================================
    
    
   !Gravando no arquivo de saída
    
    CALL Temp_open_file3(Output) !Unit 5
    
100 FORMAT(I7.0,' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
150 FORMAT(I7.0,' ',I7.0,' ',I7.0)

350 FORMAT('1',' ',ES15.7e3,' ',ES15.7e3,' ',ES15.7e3)
400 FORMAT(ES15.7e3,' ',ES15.7e3)
    WRITE(5,*) 'Tipo_Estrutura'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS OU TRUSS'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'TRUSS'
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Coordenadas'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID		x1		   y1		  z1'
    WRITE(5,*) '-------------------------------------'
    DO i=1, node_total
        WRITE(5,100) i, Nodes(i,2), Nodes(i,3), Nodes(i,4)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Atribuir_Propriedades'
    WRITE(5,*) '--------------------------------------'
    WRITE(5,*) 'ID ELEM     Material  Seção   '
    WRITE(5,*) '--------------------------------------'
    DO i=1, element_total
        WRITE(5,155) i, Elements(i,8)
155     FORMAT(I7.0,'    1   ',F7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Incidencia'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID elem      Nó 1         Nó 2'
    WRITE(5,*) '-------------------------------------'
    DO i=1, element_total
        WRITE(5,160) Elements(i,1), Elements(i,10), Elements(i,11)
160     FORMAT(F7.0,' ',F7.0,' ',F7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Secao'
    WRITE(5,*) '-------------------------------------'
    WRITE(5,*) 'ID			 Área em m²'
    WRITE(5,*) '-------------------------------------'
    gen_size = count_cell + count_xyz
    DO i=1, gen_size            
        ID = i
        Area = 1.0d0
        WRITE(5,250) ID, Area
250     FORMAT(F7.0,' ',F7.0)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Cond_Contorno'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z   (1=restrito, 0=livre)'
    WRITE(5,*) '------------------------------------'
    DO i=1, node_total
        WRITE(5,300) i
300     FORMAT(I7.0,' ','0  0.0d0',' ','0  0.0d0',' ','0  0.0d0')
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Material'
    WRITE(5,*) '--------------------------------------------'
    WRITE(5,*) 'ID  Mód. de elasticidade		  Densidade   '
    WRITE(5,*) '--------------------------------------------'
    DO i=1, Self%N_Materiais
        ID = GET_material_ID(Self%Material(i))
        Young = GET_material_Young(Self%Material(i))
        Densidade = GET_material_Densidade(Self%Material(i))
        WRITE(5,200) ID, Young, Densidade
200     FORMAT(F7.0,' ',ES15.7e3,' ',ES15.7e3)
    END DO
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Forcas'
    WRITE(5,*) '------------------------------------'
    WRITE(5,*) 'ID Nó      X        Y        Z'
    WRITE(5,*) '------------------------------------'
        X1 = GET_Fext(Self%Forcas(1),1)
        Y1 = GET_Fext(Self%Forcas(1),2)
        Z1 = GET_Fext(Self%Forcas(1),3)
        WRITE(5,350)  X1, Y1, Z1
    WRITE(5,*) '====================================='
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*)   
    WRITE(5,*) 'Tempo '
    WRITE(5,*) '--------------------------------------------------------'
    WRITE(5,*) 'Dt (incremento)				N_Dt (Número de incrementos)'
    WRITE(5,*) '--------------------------------------------------------'
	WRITE(5,400) Self%DT, Self%N_Dt
    WRITE(5,*) '========================================================'
    WRITE(5,*)
    WRITE(5,*) 'Fim_do_Arquivo'
    DEALLOCATE (Elements)
    Close(5)
    gen_size = count_cell + count_xyz
    END SUBROUTINE Cell_creator_isotropic
    !----------------------------------------------------------------------------------------------
    SUBROUTINE CALC_Reacoes(Self)
    IMPLICIT NONE
    TYPE(Type_Problem),     INTENT(INOUT)   ::  Self
    INTEGER*4                               ::  i
    INTEGER*4                               ::  j
    INTEGER*4                               ::  N_DOF
    INTEGER*4                               ::  DOF

    
    N_DOF = Self%N_DOF
    i = Self%N_Plano_xf
    
    ALLOCATE (Self%Reac(i,N_DOF))
    
    
    
    DO i=1, Self%N_Plano_xf
        DOF = 3*Self%Plano_xf(i)-2
        DO j=1, Self%N_DOF
            Self%Reac(i,j) = Self%H(DOF,j)
        END DO
    END DO
    
    !DEALLOCATE (Self%H)
    
    Self%Delta_pos(:) = Self%pos_f(:) - Self%pos_0(:)
    ALLOCATE (Self%Reacao(Self%N_Plano_xf))
    
    
    Self%Reacao = MATMUL(Self%Reac,Self%Delta_pos)
    Self%F_Reacao = 0.0d0
    DO i=1, Self%N_Plano_xf
        Self%F_Reacao = Self%F_Reacao + Self%Reacao(i)
    END DO
    
    END SUBROUTINE CALC_Reacoes
    !----------------------------------------------------------------------------------------------
    
END MODULE Din_Program