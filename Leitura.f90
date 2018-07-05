MODULE Leitura

   IMPLICIT NONE 
   PUBLIC
   CONTAINS
   !========================================================================================
        SUBROUTINE Buscar(buscar_palavra,unidade)
             IMPLICIT NONE
             character*90   ::    buscar_palavra, ler_palavra
             INTEGER*4      ::      unidade
             REWIND(unidade)                                       !Volta o cursor para o início do arquivo
             ler_palavra = ''
             
             DO WHILE(TRIM(buscar_palavra)/=TRIM(ler_palavra))      !Remover espaços vazios localizados na estremidade e no final
                  READ(unidade,*) ler_palavra
                  IF(TRIM(ler_palavra)=='Fim_do_Arquivo') RETURN    !Chegou no arquivo de dados e não encontrou a palavra
             END DO
             IF(TRIM(buscar_palavra)==TRIM(ler_palavra)) THEN       !Pular as 3 linhas da barra superior do bloco de dados
                       READ(unidade,*)
                       READ(unidade,*)
                       READ(unidade,*)
             END IF
             RETURN                                                 !
        END SUBROUTINE Buscar
   !========================================================================================
        FUNCTION Counter(buscar_palavra,unidade)      RESULT(Numero_contar) !Contar o número de itens no arquivo de dados
                IMPLICIT NONE
                character*90        ::          buscar_palavra      !Palavra a ser encontrada
                CHARACTER           ::          Fin_char            !Variável para finalizar a contagem
                INTEGER*4           ::          Numero_contar       !Contador
                INTEGER*4           ::          unidade
                CALL Buscar(buscar_palavra,unidade)                         !Coloca o cursor na posição do primeiro item
                Numero_contar = -1                                  !-1 para descontar o loop final
                Fin_char = ''                                       
                DO WHILE(Fin_char /= '=')                           !O bloco de dados termina em "="
                      READ(unidade,*) Fin_char
                      Numero_contar = numero_contar + 1
                END DO 
                RETURN                                              !
        END FUNCTION Counter
   !========================================================================================   
   END MODULE Leitura