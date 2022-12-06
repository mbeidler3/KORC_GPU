!-----------------------------------------------------------------------
!     $Id: edge_type_mod.f90 6393 2019-02-12 16:34:01Z jking $
!     module containing edge_type and vertex_type definitions.          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     declarations.                                                     
!-----------------------------------------------------------------------
      MODULE edge_type_mod 
      USE local 
      IMPLICIT NONE 
                                                                        
      TYPE :: vertex_type 
        INTEGER(i4), DIMENSION(2) :: intxy 
        INTEGER(i4) :: nimage 
        INTEGER(i4), DIMENSION(:), POINTER :: order 
        INTEGER(i4), DIMENSION(:,:), POINTER :: ptr,ptr2 
        INTEGER(i4) :: applV0=0
        REAL(r8), DIMENSION(:), POINTER :: seam_in,seam_out,tang,norm,  &
     &            seam_save, grdn
        REAL(r8), DIMENSION(:,:), POINTER :: seam_hold 
        COMPLEX(r8), DIMENSION(:), POINTER :: seam_cin,seam_cout,       &
     &            seam_csave                                            
        COMPLEX(r8), DIMENSION(:,:), POINTER :: seam_chold 
        REAL(r8) :: ave_factor,ave_factor_pre,rgeom,zgeom 
      END TYPE vertex_type 
                                                                        
      TYPE :: segment_type 
        INTEGER(i4), DIMENSION(2) :: intxyn,intxyp,intxys,ptr 
        INTEGER(i4) :: applV0=0
        LOGICAL :: h_side 
        REAL(r8), DIMENSION(:), POINTER :: seam_in,seam_out,            &
     &            seam_save, grdn
        REAL(r8), DIMENSION(:,:), POINTER :: tang,norm 
        REAL(r8), DIMENSION(:,:,:), POINTER :: seam_mat_in,seam_mat_out,&
     &            seam_mat_save                                         
        COMPLEX(r8), DIMENSION(:), POINTER :: seam_cin,                 &
     &            seam_cout,seam_csave                                  
        COMPLEX(r8), DIMENSION(:,:,:), POINTER :: seam_mat_cin,         &
     &            seam_mat_cout,seam_mat_csave                          
        REAL(r8) :: ave_factor,ave_factor_pre 
      END TYPE segment_type 
                                                                        
      TYPE :: edge_type 
        CHARACTER(64) :: name 
        INTEGER(i4) :: id 
        INTEGER(i4) :: nvert,npt
        INTEGER(i4) :: nexvert
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ive
        TYPE(vertex_type), DIMENSION(:), POINTER :: vertex 
        TYPE(segment_type), DIMENSION(:), POINTER :: segment 
        CHARACTER(1), DIMENSION(:), POINTER :: vtype 
        LOGICAL, DIMENSION(:), POINTER :: expoint 
        LOGICAL, DIMENSION(:), POINTER :: excorner 
        LOGICAL, DIMENSION(:), POINTER :: r0point 
      END TYPE edge_type 
!-----------------------------------------------------------------------
!     close module                                                      
!-----------------------------------------------------------------------
      END MODULE edge_type_mod
