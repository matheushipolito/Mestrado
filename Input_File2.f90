MODULE Input_File

IMPLICIT NONE
PUBLIC
    CONTAINS
    SUBROUTINE Open_File(File_Name, Unit_Number)
    character*90     ,       INTENT(IN)  ::   File_Name                                 ! Variavel de entrada só
    INTEGER*4        ,       INTENT(IN)  ::   Unit_Number
    LOGICAL exist_file
    INQUIRE(FILE = TRIM(file_name), EXIST = exist_file)
    IF(exist_file)THEN
        OPEN(Unit_Number, FILE = file_name, STATUS = 'OLD' )
    ELSE
        OPEN(Unit_Number, FILE = file_name, STATUS = 'UNKNOWN' )
    END IF
    END SUBROUTINE  Open_File

    SUBROUTINE Close_File(Unit_Number)
    INTEGER*4        ,       INTENT(IN)  ::   Unit_Number
    CLOSE(Unit_Number)                                                                  ! Fechar o arquivo
    END SUBROUTINE Close_File

END MODULE Input_File