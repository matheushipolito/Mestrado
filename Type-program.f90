
    PROGRAM Type_Program
    USE Din_Program
    USE Gene
    IMPLICIT NONE
    INTEGER*4           ::       i
    INTEGER*4           ::       j
    INTEGER*4           ::       k
    INTEGER*4           ::       Ninhada
    INTEGER*4           ::       invasores
    INTEGER*4           ::       ger
    INTEGER*4           ::       max_ger
    INTEGER*4           ::       n_pai
    INTEGER*4           ::       pai(10)
    REAL*8, ALLOCATABLE ::       pai_rate(:)
    REAL*8              ::       pai_res(10)
    INTEGER*8           ::       pop_size
    INTEGER*8           ::       max_unit
    CHARACTER*90        ::       File       !entrada do repetidor
    CHARACTER*90        ::       File1     
    CHARACTER*8         ::       fmt 
    CHARACTER*90        ::       t1
    CHARACTER*90        ::       File2  
    CHARACTER*90        ::       File3
    CHARACTER*90        ::       base
    CHARACTER*90, ALLOCATABLE        ::       Output(:)     !Saída dos passos
    CHARACTER*90        ::       saida(10)   
    CHARACTER*90, ALLOCATABLE ::       Individuos(:)  
    INTEGER*8           ::       gen_size
    REAL*8, ALLOCATABLE ::       Gen(:,:)
    REAL*8              ::       Valores(14)
    REAL*8, ALLOCATABLE ::       Pontos(:)
    REAL*8, ALLOCATABLE ::       Sopa_genetica(:)
    REAL*8, ALLOCATABLE ::       Encubadora(:,:)
    REAL*8              ::       Resultado
    REAL*8              ::       Resultado2
    REAL*8              ::       Res_anterior
    REAL*8              ::       evolucao
    REAL*8              ::       evolucao2
    REAL*8              ::       mut_rate 
    REAL*8              ::       virus_rate
    INTEGER*4           ::       virus_effect
    REAL*8              ::       min
    REAL*8              ::       max
    REAL*8              ::       Poisson
    INTEGER*4           ::       discreto
    TYPE(Type_Problem)  ::       TRUSS
    TYPE(Type_Problem)  ::       TRUSS2
    
    Base ='isotropico.txt'                                  !Estrutura inicial
    File='Cubo.txt'                                         !Estrutura revolucionada
    File1='espelhado.txt'                                   !Celula base
    File2='Celula_tronco'                                   !Celulas modificada geneticamente
    File3='repetido.txt'                                    !Celulas modificada geneticamente repetidas
    Output='C:\\Users\\Matheus Cenoura\\Desktop\\pais\\saida'
    
    
    saida(1)='C:\\Users\\Matheus Cenoura\\Desktop\\pais\\Pai-1.txt'
    saida(2)='C:\\Users\\Matheus Cenoura\\Desktop\\pais\\Pai-2.txt'
    saida(3)='C:\\Users\\Matheus Cenoura\\Desktop\\pais\\Pai-3.txt' 
    saida(4)='C:\\Users\\Matheus Cenoura\\Desktop\\pais\\Pai-4.txt'
    saida(5)='C:\\Users\\Matheus Cenoura\\Desktop\\pais\\Pai-5.txt'
    saida(6)='C:\\Users\\Matheus Cenoura\\Desktop\\pais\\Pai-6.txt'
    saida(7)='C:\\Users\\Matheus Cenoura\\Desktop\\pais\\Pai-7.txt'
    saida(8)='C:\\Users\\Matheus Cenoura\\Desktop\\pais\\Pai-8.txt'
    saida(9)='C:\\Users\\Matheus Cenoura\\Desktop\\pais\\Pai-9.txt'
    saida(10)='C:\\Users\\Matheus Cenoura\\Desktop\\pais\\Pai-10.txt'

    
    
    
    !========================================================================================
    !=======Definições da otimização genética================================================
    !========================================================================================
    pop_size = 30        !Tamanho da população em cada geração
    ger = 0             !Colocar 1 se deseja continuar uma otimização anterior, 0 se deseja iniciar uma nova  
    max_ger = 100        !Máximo de gerações
    mut_rate = 0.02      !Taxa de mutação aleatória
    virus_rate = 0.0d0  !Taxa de ação do vírus na população
    virus_effect = 0.0d0   !Taxa em porcento da redução máxima causada pelo vírus no gene (se for discreto, a redução será sempre 1)
    min = 1.0d0             !Menor valor do gene
    max = 3.0d1            !Maior valor do gene
    n_pai = 5           !Número de pais que passam genes para a geração seguinte
    invasores = 0       !Número de individuos aleatórios inseridos em cada populacao
    discreto = 0        !1 - discreto, 0 - continuo
    Poisson = 4.5d-1
    
    ALLOCATE(pai_rate(10))
    
    Pai_rate(1) = 50    !Número em porcentagem de influência nos filhos
    Pai_rate(2) = 20    !Número em porcentagem de influência nos filhos
    Pai_rate(3) = 15    !Número em porcentagem de influência nos filhos
    Pai_rate(4) = 10     !Número em porcentagem de influência nos filhos
    Pai_rate(5) = 5     !Número em porcentagem de influência nos filhos
    Pai_rate(6) = 0     !Número em porcentagem de influência nos filhos
    Pai_rate(7) = 0     !Número em porcentagem de influência nos filhos
    Pai_rate(8) = 0     !Número em porcentagem de influência nos filhos
    Pai_rate(9) = 0     !Número em porcentagem de influência nos filhos
    Pai_rate(10) = 0    !Número em porcentagem de influência nos filhos
    !========================================================================================
    Valores(1) = 1.9635d-007
    Valores(2) = 2.8274d-007
    Valores(3) = 3.8484d-007
    Valores(4) = 5.0264d-007
    Valores(5) = 6.3617d-007
    Valores(6) = 7.8540d-007
    Valores(7) = 9.5033d-007
    Valores(8) = 1.13097d-006
    Valores(9) = 1.32732d-006
    Valores(10) = 1.53938d-006
    Valores(11) = 1.76714d-006
    Valores(12) = 2.02732d-006
    Valores(13) = 2.53938d-006
    Valores(14) = 3.00014d-006
    
    !========================================================================================
    !=======Criador de célula================================================================
    !========================================================================================
    CALL Read_File(TRUSS,base,1)
    CALL Cell_creator_isotropic(TRUSS,File,gen_size)
    !CALL Trelica_est_nao_Lin(TRUSS,File3)
    Close(1)
    CALL Read_File(TRUSS2,File,1)
    !CALL Cell_Creator(TRUSS,File1,2,2,2,gen_size)!TRUSS, Saída, repetições em X, Y, Z
    CALL Cell_Creator_mirror(TRUSS2,File1,gen_size)!TRUSS, Saída
    !CALL Read_File(TRUSS2,File1)
    !CALL Trelica_est_nao_Lin(TRUSS2,Output)
    Close(1)
    !========================================================================================
    ALLOCATE (Gen(pop_size,gen_size),Encubadora(100,gen_size))
    ALLOCATE (Pontos(pop_size),Sopa_genetica(pop_size),Individuos(pop_size),Output(pop_size))
    !========================================================================================
    !=======Criando genes iniciais aleatórios================================================
    !========================================================================================   
     
    CALL RANDOM_NUMBER(Gen)
    
    !Colocando limites superior e inferior em cada gene

    Gen = Gen * (max-min) + min

    IF (discreto.EQ.1) THEN
        Gen = INT(Gen) + 1
        
        DO i=1, pop_size
            DO j=1, gen_size
                k = Gen(i,j)
                Gen(i,j) = Valores(k)
            END DO
        END DO
        
    END IF
    
    
    !========================================================================================
    !=====Iniciando loop das gerações========================================================
    !========================================================================================
    Ninhada = pop_size - n_pai - invasores   !Tamanho da população de filhos em cada geração
    evolucao = 10
    Res_anterior = 1
    DO WHILE ((ger.LE.max_ger))!.AND.(evolucao.GE.1.0d-2))
        ger = ger + 1
    !========================================================================================
    !=====Lendo pais da geração anterior=====================================================
    !========================================================================================
    IF (ger.NE.1) THEN
        DO i=1,n_pai
            CALL Open_File(saida(i), 3)
            CALL Ler_pai(3, Gen(i,:),gen_size)
        END DO
    
        CALL Open_File(saida(1), 3)
        CALL Ler_pai_res(3, Res_anterior)
        Close(3)
    
        
    !========================================================================================
    !=======Criando prole====================================================================
    !========================================================================================
        k=1
        DO i=1, n_pai
            IF (pai_rate(i).NE.0) THEN
                DO j=1, pai_rate(i)
                    Encubadora(k,:) = Gen(i,:)
                    k = k+1
                END DO
            END IF
        END DO
        
        DO i=1, Ninhada
            DO j=1, gen_size 
                CALL RANDOM_NUMBER(evolucao)
                k = evolucao*100
                IF (k.EQ.0) THEN
                    k=1
                END IF
                Gen((n_pai+i),j) = Encubadora(k,j)
            END DO
        END DO
        
    
    !========================================================================================
    !=======Criando mutação==================================================================
    !========================================================================================
    IF (ger.LE.(max_ger-(max_ger/10))) THEN
        DO i=1, Ninhada
            DO j=1, gen_size
                CALL RANDOM_NUMBER(evolucao)
                IF (evolucao.LE.mut_rate) THEN
                    CALL RANDOM_NUMBER(evolucao)
                    IF (discreto.EQ.1) THEN
                        evolucao = INT(evolucao * max + 1)
                        Gen((n_pai+i),j) = Valores(evolucao)
                    ELSE
                        Gen((n_pai+i),j) = evolucao * (max-min) + min
                    END IF
                END IF
            END DO
        END DO
    END IF
    END IF
    !========================================================================================
    !=======Simulação da vida================================================================
    !========================================================================================
    
    IF (ger.NE.1) THEN
        DO i=1, n_pai
            CALL Open_File(saida(i), 3)
            CALL Ler_pai_res(3, Pontos(i))
            Close(3)
        END DO
        j=n_pai+1
    ELSE
        j=1
    END IF
    
      
    DO i=1,pop_size  
        fmt = '(I6.6)'
        WRITE(t1,fmt) i
        Individuos(i) = ('espelhado'//TRIM(t1)//'.txt')
        Output(i) = ('Resultado_ind'//TRIM(t1)//'_')
        
    END DO

    
    CALL GENE_CONVERTER(File1,Gen,discreto,Individuos)
    
        !$OMP PARALLEL
        !$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(i)
        DO i=j,pop_size
            !Print*,'geracao: ',ger,' / individuo: ',i
            !Print*,' '
            CALL Life_Simulator(Individuos(i),resultado,Output(i),i,Poisson)
            !Resultado é o valor da função objetivo do gene
            Pontos(i) = resultado
            IF (resultado.NE.1.0d12) THEN
                resultado2=1.0d0/resultado
                Print*,'Individuo',i,' viveu: ',resultado2,'anos'
            ELSE
                Print*,'Individuo',i,': morreu no nascimento'
            END IF
            !Print*,' '
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
   
    !========================================================================================
    !=======Pegando os indivíduos mais aptos=================================================
    !========================================================================================
        pai=0
        pai_res=9990
        Sopa_genetica = Pontos
        CALL KB06AD(Sopa_genetica,pop_size)
        
        
        
        
        DO j=1, n_pai
            DO i=1, pop_size
                IF (Sopa_genetica(pop_size+1-j).EQ.Pontos(i)) THEN!minimizar
                !IF (Sopa_genetica(j).EQ.Pontos(i)) THEN !maximizar
                    IF (j.EQ.1) THEN
                        Pai(1) = i
                        pai_res(1) = Pontos(i)
                    ELSE IF (Pai_res(j-1).NE.Pontos(i)) THEN
                        Pai(j) = i
                        pai_res(j) = Pontos(i)
                    END IF
                END IF
            END DO
        END DO
        
        resultado2 = 1.0d0/Pai_res(1)
        Print*,'geracao: ',ger,' chegou ao fim, o individuo mais velho viveu',resultado2,'anos'
        
        
        DO i=1,n_pai
            IF (pai_res(i).NE.9990) THEN
                CALL Open_File(saida(i), 3)
                CALL Imprimir_pai(3,Gen(Pai(i),:),gen_size,pai_res(i))
                Close(3)
            END IF
        END DO
        Close(1)
    
        
        
        IF (ger.NE.1) THEN
            evolucao = Res_anterior - pai_res(1)
        END IF
        
    END DO
    
    
    
    
    end program Type_Program

