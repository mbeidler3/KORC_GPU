!-----------------------------------------------------------------------
!     $Id: computation_pointers.f90 5161 2017-02-21 23:52:27Z tbechtel $
!     module containing pointers for use during finite element          
!     computations.  the rhs, sln, and cell_rhs pointers are            
!     associated with (otherwise) unnamed memory locations.             
!-----------------------------------------------------------------------
      MODULE computation_pointers 
      USE vector_type_mod 
      IMPLICIT NONE 
                                                                        
!-----------------------------------------------------------------------
!     working array structures.                                         
!-----------------------------------------------------------------------
      TYPE(vector_type), DIMENSION(:), POINTER :: rhs,cell_rhs,sln,     &
     &                   vectr,lump_mass,lump_summed                    
      TYPE(cvector_type), DIMENSION(:), POINTER :: crhs,cell_crhs,cvecn,&
     &                    cvecnn                                        
      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: csln,cvectr,ctemp
!-----------------------------------------------------------------------
!     close module                                                      
!-----------------------------------------------------------------------
      END MODULE computation_pointers 
