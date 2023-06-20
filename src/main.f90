program main

use interp_gpu
use params_gpu
use pusher_gpu

implicit none

REAL(rp),ALLOCATABLE,DIMENSION(:) :: X_X,X_Y,X_Z
REAL(rp),ALLOCATABLE,DIMENSION(:) :: V_X,V_Y,V_Z
REAL(rp),ALLOCATABLE,DIMENSION(:) :: B_X,B_Y,B_Z
REAL(rp),ALLOCATABLE,DIMENSION(:) :: E_X,E_Y,E_Z
REAL(rp),DIMENSION(20) :: XF=0._rp,YF=0._rp
REAL(rp),DIMENSION(20,20) :: BF_X=0._rp,BF_Y=0._rp,BF_Z=0._rp
REAL(rp),DIMENSION(20,20) :: EF_X=0._rp,EF_Y=0._rp,EF_Z=0._rp
#ifdef PSPLINE
TYPE(KORC_2D_FIELDS_INTERPOLANT)      :: bfield_2d
!! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
!! the magnetic field.
TYPE(KORC_2D_FIELDS_INTERPOLANT)   :: efield_2d
!! An instance of KORC_2D_FIELDS_INTERPOLANT for interpolating
!! the electric field.
#endif PSPLINE
NAMELIST /input_parameters/ nRE,simulation_time,field_type

#ifdef ACC 
#ifdef PSPLINE
  !$acc routine (interp_fields) seq
#endif 
#endif

!! initialize system_clock
CALL system_clock(count_rate=cr)
rate = REAL(cr)

CALL system_clock(c1)

!! find input/output file
argn = command_argument_count()

call get_command_argument(1,path_to_inputs)
call get_command_argument(2,path_to_outputs)

!! open log and data output files
OPEN(UNIT=data_write, &
FILE=TRIM(path_to_outputs)//"data.korc", &
STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')

OPEN(UNIT=output_write, &
FILE=TRIM(path_to_outputs)//"output.korc", &
STATUS='UNKNOWN',FORM='FORMATTED',POSITION='REWIND')

!! set defaults for inputs, open input file, read from input file

nRE=1
simulation_time=0._rp
field_type='UNIFORM'

OPEN(UNIT=default_unit_open,FILE=TRIM(path_to_inputs), &
STATUS='OLD',FORM='formatted',POSITION='REWIND')
read(default_unit_open,nml=input_parameters,IOSTAT=read_stat)
close(default_unit_open)

!! Allocate solution fields

ALLOCATE(X_X(nRE))
ALLOCATE(X_Y(nRE))
ALLOCATE(X_Z(nRE))
ALLOCATE(V_X(nRE))
ALLOCATE(V_Y(nRE))
ALLOCATE(V_Z(nRE))
ALLOCATE(B_X(nRE))
ALLOCATE(B_Y(nRE))
ALLOCATE(B_Z(nRE))
ALLOCATE(E_X(nRE))
ALLOCATE(E_Y(nRE))
ALLOCATE(E_Z(nRE))
ALLOCATE(rnd1(nRE))
ALLOCATE(gam(nRE))

!! Initialize fields
write(output_write,'("* * * * * * * * * Fields * * * * * * * * *")')
if (field_type.eq.'UNIFORM') then
   
   write(output_write,'("* * * * USING UNIFORM MAGNETIC FIELD * * * *")')

else if (field_type.eq.'PSPLINE') then

   do it=1,size(XF)
      XF(it)=-0.1+(it-1)*0.2/(size(XF)-1)
      do jt=1,size(YF)
         BF_X(it,jt)=0._rp
         BF_Y(it,jt)=0._rp
         BF_Z(it,jt)=1._rp

         EF_X(it,jt)=0._rp
         EF_Y(it,jt)=0._rp
         EF_Z(it,jt)=0._rp
      enddo
   enddo
   YF=XF
endif

!! Set kinetic energy and pitch, random gyrophase, and then velocity
Eo=10E6_rp
eta0=90._rp*C_PI/180._rp
gam0=1._rp+(Eo*C_E/(C_ME*C_C**2))
v0=C_C*sqrt(1._rp-1/gam0**2)

!call RANDOM_NUMBER(rnd1)
!chi0=2*C_PI*rnd1(1)
chi0=0

v_norm=C_C
if (field_type.eq.'UNIFORM') then
   B_norm=1._rp
else if (field_type.eq.'PSPLINE') then
   B_norm=BF_Z(1,1)
end if
t_norm=C_ME/(C_E*B_norm)
x_norm=V_norm*t_norm

#ifdef ACC
!$acc  parallel loop
#endif ACC
do pp=1,nRE
   !initialize particle fields
   B_X(pp)=0._rp/b_norm
   B_Y(pp)=0._rp/b_norm
   B_Z(pp)=1._rp/b_norm

   E_X(pp)=0._rp/(b_norm*v_norm)
   E_Y(pp)=0._rp/(b_norm*v_norm)
   E_Z(pp)=0._rp/(b_norm*v_norm)

   !! Initialize location
   X_X(pp)=0._rp/x_norm
   X_Y(pp)=0._rp/x_norm
   X_Z(pp)=0._rp/x_norm

   V_X(pp)=v0*sin(eta0)*cos(chi0)/v_norm
   V_Y(pp)=v0*sin(eta0)*sin(chi0)/v_norm
   V_Z(pp)=v0*cos(eta0)/v_norm
   gam(pp)=gam0
end do
#ifdef ACC  
  !$acc end parallel loop
#endif ACC

#ifdef PSPLINE
if (field_type.eq.'PSPLINE') then
   XF=XF/x_norm
   YF=YF/x_norm

   BF_X=BF_X/b_norm
   BF_Y=BF_Y/b_norm
   BF_Z=BF_Z/b_norm

   EF_X=EF_X/(b_norm*v_norm)
   EF_Y=EF_Y/(b_norm*v_norm)
   EF_Z=EF_Z/(b_norm*v_norm)

   call initialize_interpolants(XF,YF,BF_X,BF_Y,BF_Z,EF_X,EF_Y,EF_Z, &
     bfield_2d,efield_2d)

#ifdef ACC 
  !$acc  parallel loop &
  !$acc  copyin(bfield_2d,efield_2d)
  !$acc& private(X_X_loop,X_Y_loop,X_Z_loop,V_X_loop, &
  !$acc& V_Y_loop,V_Z_loop,gam_loop)
#endif ACC
   do pp=1,nRE

      X_X_loop=X_X(pp)
      X_Y_loop=X_Y(pp)
   
      B_X_loop=B_X(pp)
      B_Y_loop=B_Y(pp)
      B_Z_loop=B_Z(pp)
   
      E_X_loop=E_X(pp)
      E_Y_loop=E_Y(pp)
      E_Z_loop=E_Z(pp)

      call interp_fields(X_X_loop,X_Y_loop,B_X_loop,B_Y_loop,B_Z_loop, &
         E_X_loop,E_Y_loop,E_Z_loop,bfield_2d,efield_2d)

      X_X(pp)=X_X_loop
      X_Y(pp)=X_Y_loop

      B_X(pp)=B_X_loop
      B_Y(pp)=B_Y_loop
      B_Z(pp)=B_Z_loop

      E_X(pp)=E_X_loop
      E_Y(pp)=E_Y_loop
      E_Z(pp)=E_Z_loop
   enddo
#ifdef ACC  
  !$acc end parallel loop
#endif ACC

endif
#endif PSPLINE

write(output_write,'("B: ",E17.10,E17.10,E17.10)') &
   B_X(1)*b_norm,B_Y(1)*b_norm,B_Z(1)*b_norm
write(output_write,'("E: ",E17.10,E17.10,E17.10)') &
   E_X(1)*(b_norm*v_norm),E_Y(1)*(b_norm*v_norm),E_Z(1)*(b_norm*v_norm)

write(output_write,*) '* * * * * * * * * Initial Conditions * * * * * * * * *'

write(output_write,'("Number of electrons: ",I16)') nRE
write(output_write,'("gam0,eta0,chi0: ",E17.10,E17.10,E17.10)') gam0,eta0,chi0
write(output_write,'("X0: ",E17.10,E17.10,E17.10)') X_X(1)*x_norm,X_Y(1)*x_norm,X_Z(1)*x_norm
write(output_write,'("V0: ",E17.10,E17.10,E17.10)') V_X(1)*v_norm,V_Y(1)*v_norm,V_Z(1)*v_norm
write(data_write,'("X0: ",E17.10,E17.10,E17.10)') X_X(1)*x_norm,X_Y(1)*x_norm,X_Z(1)*x_norm
write(data_write,'("V0: ",E17.10,E17.10,E17.10)') V_X(1)*v_norm,V_Y(1)*v_norm,V_Z(1)*v_norm
flush(output_write)

!! Set timestep to resolve relativistic gyrofrequency, simulation time and
!! number of time steps
dt = 0.01_rp*(2.0_rp*C_PI/(C_E*B_Z(1)/( gam0*C_ME )))/t_norm

simulation_time=simulation_time/t_norm
t_steps=ceiling(simulation_time/dt)
if (t_steps.eq.0) then
   dt=0._rp
else
   dt=simulation_time/real(t_steps)
endif

write(output_write,*) '* * * * * * * * * Timings * * * * * * * * *'

write(output_write,'("simulation time: ",E17.10)') simulation_time*t_norm
write(output_write,'("dt: ",E17.10)') dt*t_norm
write(output_write,'("t_steps: ",I16)') t_steps
flush(output_write)

!! Particle push

!! Allocating work arrays

call system_clock(c2)

write(output_write,'("Setup time: ",E17.10)') (c2-c1)/rate

write(output_write,'("* * * * * * * * * Begin Orbits * * * * * * * * *")')

#ifdef PSPLINE
call FO_push(nRE,dt,t_steps,field_type,x_norm,v_norm,X_X,X_Y,X_Z,V_X,V_Y,V_Z, &
  gam,B_X,B_Y,B_Z,E_X,E_Y,E_Z,bfield_2d,efield_2d)
#else
call FO_push(nRE,dt,t_steps,field_type,x_norm,v_norm,X_X,X_Y,X_Z,V_X,V_Y,V_Z, &
  gam,B_X,B_Y,B_Z,E_X,E_Y,E_Z)
#endif PSPLINE

write(output_write,'("* * * * * * * * * Final Conditions * * * * * * * * *")')
write(output_write,'("V: ",E17.10,E17.10,E17.10)') V_X(1)*v_norm,V_Y(1)*v_norm,V_Z(1)*v_norm
write(output_write,'("X: ",E17.10,E17.10,E17.10)') X_X(1)*x_norm,X_Y(1)*x_norm,X_Z(1)*x_norm


call system_clock(c1)

write(output_write,'("Pusher time: ",E17.10)') (c1-c2)/rate

! * * * FINALIZING SIMULATION * * *

#ifdef PSPLINE
if (field_type.eq.'PSPLINE') then
   call finalize_interpolants(bfield_2d,efield_2d)
endif
#endif PSPLINE

write(output_write,'("KORC ran suppessfully!")')
close(output_write)
close(data_write)

end program main