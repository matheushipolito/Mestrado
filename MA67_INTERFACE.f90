MODULE MA67_INTERFACE
  IMPLICIT NONE
  !'IN ORDER TO USE THIS MODULE, FILES MA67.for, HLS_Dependencies_01.for and Dependent_BLAS_package.for are needed!'
  
  CONTAINS
  
SUBROUTINE MA67SOLVER(NE,N,ROWS, COLUMNS, VALUES, VECTOR, SOLUTION)
  IMPLICIT NONE
  !INPUT AND OUTPUT VARIABLES
  INTEGER, INTENT(IN)::NE,N ! NE - NON-NULL TERMS IN THE SPARSE SIMETRIC MATRIX; N - MATRIX ORDER (NxN)
  INTEGER, DIMENSION(NE), INTENT(IN)::ROWS, COLUMNS !VECTORS CONTAINING ROW AND COLUMN INDEXES
  REAL(8), DIMENSION(NE),INTENT(IN)::VALUES !VECTOR CONTAINING MATRIX TERMS VALUES 
  REAL(8), DIMENSION(N),INTENT(IN)::VECTOR !VECTOR OF THE EQUATION SYSTEM
  REAL(8), DIMENSION(N),INTENT(OUT)::SOLUTION !VECTOR WHERE THE SOLUTION WILL BE RETURNED
  !MA67 COMPUTATION VARIABLES  
  INTEGER::LFACT, LIFACT, ICNTL(10), INFO(24)
  INTEGER, DIMENSION(:),ALLOCATABLE::IFACT,IW
  REAL(8), DIMENSION(:),ALLOCATABLE::FACT,W
  REAL(8)::CNTL(5),RINFO(4)
  ! ALLOCATING VARIABLES FOR THE MA67 COMPUTATIONS
  LFACT = 99*(NE*2+5*N+4) !IF NECESSARY, CHANGE 2*(...) FOR 5* OR HIGHER!
  LIFACT = 99*(NE*2+12*N+10) !IF NECESSARY, CHANGE 3*(...) FOR 5* OR HIGHER!
  ALLOCATE(FACT(LFACT),IFACT(LIFACT),W(N),IW(2*N+2))
  !SET DEFAULT PARAMETERS
  CALL MA67ID(CNTL,ICNTL)
  !FACTORIZE MATRIX
  CALL MA67AD(N,NE,VALUES,ROWS,COLUMNS,LFACT,FACT,LIFACT,IFACT,CNTL,ICNTL,RINFO,INFO)
  !SOLVE THE EQUATIONS
  SOLUTION = VECTOR
  CALL MA67CD(N,LFACT,FACT,LIFACT,IFACT,W,SOLUTION,IW,ICNTL)
END SUBROUTINE MA67SOLVER


END MODULE MA67_INTERFACE
