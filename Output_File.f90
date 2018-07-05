MODULE Output_File

    IMPLICIT NONE
    PUBLIC
    CONTAINS
     
!------------------------------------------------------------------------------------------------------   
    SUBROUTINE  temp_open_file(file_name)
    character*90, INTENT(IN) :: file_name 
    LOGICAL exist_file
    INQUIRE(FILE = TRIM(file_name), EXIST = exist_file)
    IF(exist_file)THEN
        OPEN(UNIT = 3, FILE = file_name, STATUS = 'OLD' )
    ELSE
        OPEN(UNIT = 3, FILE = file_name, STATUS = 'UNKNOWN' )
    END IF
    END SUBROUTINE temp_open_file
!------------------------------------------------------------------------------------------------------        
    SUBROUTINE temp_close_file
        CLOSE(3)
    END SUBROUTINE temp_close_file
!------------------------------------------------------------------------------------------------------
    SUBROUTINE  temp_open_file2(file_name)
    CHARACTER*90, INTENT(IN) :: file_name 
    LOGICAL exist_file
    INQUIRE(FILE = TRIM(file_name), EXIST = exist_file)
    OPEN(UNIT = 4, FILE = TRIM(file_name) )
    END SUBROUTINE temp_open_file2
!------------------------------------------------------------------------------------------------------
    SUBROUTINE  temp_open_file3(file_name)
    character*90, INTENT(IN) :: file_name 
    LOGICAL exist_file
    INQUIRE(FILE = TRIM(file_name), EXIST = exist_file)
    IF(exist_file)THEN
        OPEN(UNIT = 5, FILE = file_name, STATUS = 'OLD' )
    ELSE
        OPEN(UNIT = 5, FILE = file_name, STATUS = 'UNKNOWN' )
    END IF

    END SUBROUTINE temp_open_file3
!------------------------------------------------------------------------------------------------------
    SUBROUTINE  temp_open_file4(file_name,unidade)
    character*90, INTENT(IN) :: file_name 
    INTEGER*4                   :: unidade
    LOGICAL exist_file
    INQUIRE(FILE = TRIM(file_name), EXIST = exist_file)
    IF(exist_file)THEN
        OPEN(UNIT = unidade, FILE = file_name, STATUS = 'OLD' )
    ELSE
        OPEN(UNIT = unidade, FILE = file_name, STATUS = 'UNKNOWN' )
    END IF

    END SUBROUTINE temp_open_file4
!------------------------------------------------------------------------------------------------------
    SUBROUTINE temp_close_file2
    CLOSE(4)
    END SUBROUTINE temp_close_file2
        
        
END MODULE Output_File