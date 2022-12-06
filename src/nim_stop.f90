!-----------------------------------------------------------------------
!     file nim_stop.f:                                                  
!     this is an abbreviated version of what appears in nimrod.         
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      SUBROUTINE nim_stop(message) 
      USE local 
      USE io 
      IMPLICIT NONE 
                                                                        
      CHARACTER(*), INTENT(IN) :: message 
!-----------------------------------------------------------------------
!     write completion message.                                         
!-----------------------------------------------------------------------
      WRITE(nim_wr,'(2a)') 'NIM_STOP => ', TRIM(message) 
      IF (ofname/='none') THEN
        WRITE(out_unit,'(2a)') 'NIM_STOP => ', TRIM(message) 
        CLOSE(UNIT=out_unit) 
      ENDIF

      STOP 
      END SUBROUTINE nim_stop 
