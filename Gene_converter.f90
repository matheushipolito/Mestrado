MODULE Gene
USE Input_File
USE Nodes
USE Material
USE Secao
USE Element
USE Forcas
USE Cond_Contorno
USE Leitura
USE Output_File
USE Din_Program
IMPLICIT NONE
PUBLIC


CONTAINS
!=============================================================================================================================
SUBROUTINE GENE_CONVERTER(File1,Gen,discreto,File2)
    IMPLICIT NONE
    TYPE(Type_Problem)  ::      CELLBASE
    TYPE(Type_Problem)  ::      CELL
    CHARACTER*90        ::      File1       
    CHARACTER*90        ::      File2(:)      
    REAL*8              ::      Gen(:,:)
    INTEGER*4           ::      discreto
    INTEGER*4           ::      i
    INTEGER*4           ::      j
    
    j=size(File2)

    CALL Read_File(CELLBASE,File1,1)

    DO i=1, j
            CALL Conversor_de_gene(CELLBASE,File2(i),Gen(i,:))!TRUSS, Saída, gene
        Close(1)
    END DO
    
END SUBROUTINE GENE_CONVERTER
!=============================================================================================================================
SUBROUTINE Life_Simulator(File,F_obj,Output,i,Meta)
    IMPLICIT NONE
    TYPE(Type_Problem)  ::       TRUSS 
    TYPE(Type_Problem)  ::       TRUSS2
    CHARACTER*90        ::       File 
    CHARACTER*90        ::       File2
    CHARACTER*90        ::       Output
    REAL*8              ::       Resultado(4)
    REAL*8              ::       F_obj
    REAL*8              ::       dX
    REAL*8              ::       dY
    REAL*8              ::       dZ
    REAL*8              ::       F
    REAL*8              ::       Forca
    REAL*8              ::       Poisson
    INTEGER*4           ::       i
    CHARACTER*8         ::       fmt 
    CHARACTER*90        ::       t1
    REAL*8              ::       Meta
    
    fmt = '(I6.6)'
    WRITE(t1,fmt) i
    File2 = ('temp'//TRIM(t1)//'.txt')
    
    CALL Read_File(TRUSS,File,i)
    CALL Cell_repeater_x(TRUSS,File2,1,1,1)
    Close(i)
    CALL Read_File(TRUSS2,File2,i)
    Close(i)
    CALL Trelica_est_nao_Lin(TRUSS2,Output,i)
    CALL Calc_resultados(TRUSS2,dX,dY,dZ,F)
    dX = abs(dX)    
    dY = abs(dY)
    dZ = dZ
    Close(i)
    
    Poisson = dY/dX

    F_obj = abs(Poisson - Meta)/Meta

    

    IF (dY.EQ.0.0d0) THEN
        F_obj = 1.0d12
    END IF
    

END SUBROUTINE Life_Simulator
!=============================================================================================================================
SUBROUTINE Calc_Resultados(Self,mX,mY,mZ,F)
    IMPLICIT NONE
    TYPE(Type_Problem)  ::      Self     
    INTEGER*4           ::      i
    INTEGER*4           ::      j
    INTEGER*4           ::      ID
    INTEGER*8           ::      N_Plano
    REAL*8              ::      X1
    REAL*8              ::      Y1
    REAL*8              ::      Z1
    REAL*8              ::      mX
    REAL*8              ::      mY
    REAL*8              ::      mZ
    REAL*8              ::      F
    REAL*8, ALLOCATABLE ::      X(:)
    REAL*8, ALLOCATABLE ::      Y(:)
    REAL*8, ALLOCATABLE ::      Z(:)
    CALL Planos_de_deslocamento(Self)
    ALLOCATE (X(Self%N_Plano_xf), Y(Self%N_Plano_yf), Z(Self%N_Plano_zf))
    X = 0
    Y = 0
    Z = 0
    DO i=1, Self%N_Plano_xf
        ID = Self%Plano_xf(i)
        X(i) = Get_xf(Self%Node(ID)) - Get_x0(Self%Node(ID))
    END DO
    
    DO i=1, Self%N_Plano_yf
        ID = Self%Plano_yf(i)
        Y(i) = Get_yf(Self%Node(ID)) - Get_y0(Self%Node(ID))
    END DO
    
    DO i=1, Self%N_Plano_zf
        ID = Self%Plano_zf(i)
        Z(i) = Get_zf(Self%Node(ID)) - Get_z0(Self%Node(ID))
    END DO
    
    CALL KB06AD(X,Self%N_Plano_xf)
    CALL KB06AD(Y,Self%N_Plano_yf)
    CALL KB06AD(Z,Self%N_Plano_zf)
    
    
    mX = 0
    mY = 0
    mZ = 0
    
    j = Self%N_Plano_xf/3
    DO i=1,j
        mX = mX + X(i)
    END DO
    
    mX = mX/j
    
    j = Self%N_Plano_yf/3
    DO i=1,j
        mY = mY + Y(i)
    END DO
    
    mY = mY/j
    
    j = Self%N_Plano_zf/3
    DO i=1,j
        mZ = mZ + Z(i)
    END DO
    
    mZ = mZ/j
    
    F = Self%F_Reacao
    
END SUBROUTINE Calc_Resultados
!=============================================================================================================================
SUBROUTINE Imprimir_pai(unit,Gen,N,res)
    IMPLICIT NONE
    REAL*8              ::      Gen(:)           
    INTEGER*4           ::      i
    INTEGER*4           ::      unit
    INTEGER*8           ::      N
    REAL*8              ::      res
    
    WRITE(unit,100) Res
100 FORMAT(F30.15)
    DO i=1, N
        WRITE(unit,110) Gen(i)
110     FORMAT(F22.13)
    END DO
END SUBROUTINE Imprimir_pai
!=============================================================================================================================
SUBROUTINE Ler_pai(unit,Gen,N)
    IMPLICIT NONE
    REAL*8              ::      Gen(:)           
    INTEGER*4           ::      i
    INTEGER*4           ::      unit
    INTEGER*8           ::      N
    
    READ(unit,*)
    DO i=1, N
        READ(unit,*) Gen(i)
    END DO
END SUBROUTINE Ler_pai
!=============================================================================================================================
SUBROUTINE Ler_pai_res(unit,res)
    IMPLICIT NONE
    REAL*8              ::      res          
    INTEGER*4           ::      unit

    READ(unit,*) res
END SUBROUTINE Ler_pai_res


END MODULE Gene