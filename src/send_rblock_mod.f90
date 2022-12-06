!-----------------------------------------------------------------------
!     $Id: send_rblock_mod.f90 4620 2015-12-03 19:18:40Z jking $
!     routines for sending portion of rblock dependent variables between
!     fluid and closure nodes                                           
!-----------------------------------------------------------------------
!     code organization.                                                
!-----------------------------------------------------------------------
!     1.  send_fluid_rblock                                             
!     2.  send_close_rblock                                             
!     3.  slagr_quad                                                    
!     4.  rlagr_quad                                                    
!-----------------------------------------------------------------------
      MODULE send_rblock_mod 
      USE local 
      USE rblock_type_mod 
      USE input 
      USE mpi_nim 
      USE pardata 
      USE time 
      IMPLICIT NONE 
      PUBLIC :: send_fluid_rblock, send_close_rblock 
      PRIVATE :: slagr_quad, rlagr_quad 
      INTEGER(i4) :: sproc,rproc 
      INTEGER(i4) :: status(mpi_status_size) 
      INTEGER(i4) :: notify 
                                                                        
      CONTAINS 
                                                                        
      SUBROUTINE send_fluid_rblock(rbsend,rbrecv) 
      IMPLICIT NONE 
      TYPE (rblock_type), INTENT(IN) :: rbsend 
      TYPE (rblock_type), INTENT(INOUT) :: rbrecv 
      INTEGER(i4) :: ierr 
                                                                        
      IF (global_node == sproc) then 
                                                                        
        CALL mpi_recv(notify,1,mpi_nim_int,rproc,0,                     &
     &       mpi_comm_world,status,ierr)                                
                                                                        
        CALL slagr_quad(rbsend%be) 
        CALL slagr_quad(rbsend%ve) 
        CALL slagr_quad(rbsend%nd) 
        CALL slagr_quad(rbsend%tele) 
        CALL slagr_quad(rbsend%tion) 
                                                                        
      ELSE IF (global_node == rproc) THEN 
                                                                        
        CALL mpi_send(notify,1,mpi_nim_int,sproc,0,                     &
     &       mpi_comm_world,ierr)                                       
                                                                        
        CALL rlagr_quad(rbrecv%be) 
        CALL rlagr_quad(rbrecv%ve) 
        CALL rlagr_quad(rbrecv%nd) 
        CALL rlagr_quad(rbrecv%tele) 
        CALL rlagr_quad(rbrecv%tion) 
                                                                        
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE send_fluid_rblock 
                                                                        
      SUBROUTINE send_close_rblock(rbsend,rbrecv) 
      IMPLICIT NONE 
      TYPE (rblock_type), INTENT(IN) :: rbsend 
      TYPE (rblock_type), INTENT(INOUT) :: rbrecv 
      INTEGER(i4) :: ierr 
                                                                        
      IF (global_node == sproc) then 
                                                                        
        CALL mpi_recv(notify,1,mpi_nim_int,rproc,0,                     &
     &       mpi_comm_world,status,ierr)                                
                                                                        
        CALL slagr_quad(rbsend%qpe) 
        CALL slagr_quad(rbsend%qpi) 
        CALL slagr_quad(rbsend%ppe) 
        CALL slagr_quad(rbsend%ppi) 
        CALL slagr_quad(rbsend%phot) 
                                                                        
        CALL mpi_send(time_cel,1,mpi_nim_real,rproc,0,                  &
     &       mpi_comm_world,ierr)                                       
                                                                        
      ELSE IF (global_node == rproc) THEN 
                                                                        
        CALL mpi_send(notify,1,mpi_nim_int,sproc,0,                     &
     &       mpi_comm_world,ierr)                                       
                                                                        
        CALL rlagr_quad(rbrecv%qpe) 
        CALL rlagr_quad(rbrecv%qpi) 
        CALL rlagr_quad(rbrecv%ppe) 
        CALL rlagr_quad(rbrecv%ppi) 
        CALL rlagr_quad(rbrecv%phot) 
                                                                        
        CALL mpi_recv(time_cel,1,mpi_nim_real,sproc,0,                  &
     &       mpi_comm_world,status,ierr)                                
                                                                        
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE send_close_rblock 
                                                                        
      SUBROUTINE slagr_quad(laq) 
      IMPLICIT NONE 
      INTEGER(i4) :: ierr 
                                                                        
      TYPE(lagr_quad_type), INTENT(IN) :: laq 
                                                                        
      CALL mpi_send(laq%fs,SIZE(laq%fs),                                &
     &     mpi_nim_comp,rproc,0,mpi_comm_world,ierr)                    
      IF (laq%n_side>0) THEN 
        CALL mpi_send(laq%fsh,SIZE(laq%fsh),                            &
     &       mpi_nim_comp,rproc,0,mpi_comm_world,ierr)                  
        CALL mpi_send(laq%fsv,SIZE(laq%fsv),                            &
     &       mpi_nim_comp,rproc,0,mpi_comm_world,ierr)                  
      ENDIF 
      IF (laq%n_int>0) THEN 
        CALL mpi_send(laq%fsi,SIZE(laq%fsi),                            &
     &       mpi_nim_comp,rproc,0,mpi_comm_world,ierr)                  
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE slagr_quad 
                                                                        
      SUBROUTINE rlagr_quad(laq) 
      IMPLICIT NONE 
                                                                        
      TYPE(lagr_quad_type), INTENT(INOUT) :: laq 
      INTEGER(i4) :: ierr 
                                                                        
      CALL mpi_recv(laq%fs,SIZE(laq%fs),                                &
     &     mpi_nim_comp,sproc,0,mpi_comm_world,status,ierr)             
      IF (laq%n_side>0) THEN 
        CALL mpi_recv(laq%fsh,SIZE(laq%fsh),                            &
     &       mpi_nim_comp,sproc,0,mpi_comm_world,status,ierr)           
        CALL mpi_recv(laq%fsv,SIZE(laq%fsv),                            &
     &       mpi_nim_comp,sproc,0,mpi_comm_world,status,ierr)           
      ENDIF 
      IF (laq%n_int>0) THEN 
        CALL mpi_recv(laq%fsi,SIZE(laq%fsi),                            &
     &       mpi_nim_comp,sproc,0,mpi_comm_world,status,ierr)           
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE rlagr_quad 
                                                                        
      END MODULE send_rblock_mod 
