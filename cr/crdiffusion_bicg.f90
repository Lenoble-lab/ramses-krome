subroutine crdiffusion_cg (ilevel,Nsub)
  use amr_commons
  use hydro_commons
  use cooling_module,ONLY:kB,mH,clight
  !  use radiation_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !===========================================================================
  ! Iterative solver with Conjugate Gradient method and Newton Raphson
  ! to solve A x = b
  !   r1      : stored in unew(i,1)
  !   p1      : stored in unew(i,2)
  !   Diag(A) : stored in unew(i,4)
  !   Ap1     : stored in unew(i,3)
  !  x1(i)    : stored in uold(i,8+nrad+ntp+ncr)(i.e. new thermal energy at time n+1)
  !  b1(n)    : stored in unew(i,8+nrad+ntp+ncr)(i.e. thermal energy at time n)
  ! x1(i-1)   : stored in unew(i,7)
  !===========================================================================
  integer,intent(IN)::ilevel,Nsub
  complex*16 :: final_sum
  real(dp)::error,error_ini,epsilon,error_nr,error_nrm1,error_nrm2,error_nrm3
  real(dp)::error_nr_loc,error_nr_all,error_cg_loc,error_cg_all
  real(dp)::alpha_cg,beta_cg,Cv,told,tnew,wdt,rho,dt_exp
  real(dp)::r2_old,r2,pAp,rhs_norm1
  real(dp)::density,kpara_ana,kperp_ana,norm_Eth,temp
  integer :: i,info,ind,iter,iskip,itermax,icpu,igroup,igrp
  integer :: this,iter_nr,nx_loc,nsub_imp,isub
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::rhs,lhs,scale_kappa
  real(dp)::dx,dx_loc,surf_loc,vol_loc,scale
  real(dp)::min_ener,min_ener_all,max_ener,max_ener_all
  real(dp)::dener
  logical::exist_leaf_cell=.true.

  real(dp)::omega,scale_omega,Ti_new,Ti_old,T_new
  real(dp)::n_electron,n_ion
  integer ::ind_cr,ivar
  if(verbose)write(*,111)
  if(numbtot(1,ilevel)==0)return

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kappa=scale_d*scale_l*scale_v**3
  scale_omega = (scale_d*scale_v**2*scale_l**3)*(gamma-1.0d0)/(dtnew(ilevel)*kB)
  !  call compute_Tau_equi_vol(Tau_ei_vol)

  ! Rescaling factors
  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  surf_loc = dx_loc**(ndim-1)
  vol_loc  = dx_loc**ndim

  allocate(liste_ind (1:twotondim*active(ilevel)%ngrid))

  nb_ind = 0

  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        if(son(active(ilevel)%igrid(i)+iskip) == 0)then
           nb_ind = nb_ind+1 
           liste_ind(nb_ind) = active(ilevel)%igrid(i)+iskip
        end if
     end do
  end do

!!$  do ind=1,twotondim
!!$     iskip=ncoarse+(ind-1)*ngridmax
!!$     do icpu=1,ncpu
!!$        do i=1,reception(icpu,ilevel)%ngrid
!!$           heat_flux(reception(icpu,ilevel)%igrid(i)+iskip)=0.0d0
!!$        end do
!!$        do i=1,reception(icpu,ilevel-1)%ngrid
!!$           heat_flux(reception(icpu,ilevel-1)%igrid(i)+iskip)=0.0d0
!!$        end do
!!$     end do
!!$  end do

  if (nb_ind == 0)then
     print*,'No leaf-cell - myid=',myid
     exist_leaf_cell=.false.
  end if

  do i=1,nb_ind
     this = liste_ind(i)
     unew(this,1:nvar)=0.0d0
!!$     unew(this,1:8+nrad+ntp)=0.0d0
!!$     unew(this,ind_Cve)=0.d0
!!$     unew(this,ind_eint)=0.d0
!!$     unew(this,ind_Teold)=0.0d0
!!$     divu(this)=0.0d0
!!$     enew(this)=0.0d0
  end do

  !===================================================================
  ! Compute gas temperature stored in unew(i,ind_Teold) and in unew(i,ind_Tenew)
  ! 1: Etot     -> Etot-Ecr
  ! 2: Etot-Ecr -> Etot
  !===================================================================
  call cmp_energy_cr(1)

  do ivar=1,ncr
     ind_cr = 8+ivar
     if(twotemp)ind_cr = 9+ivar

  !===================================================================
  ! Begin of subcycles....
  !===================================================================
  dt_exp = dtnew(ilevel)
  dt_imp = dtnew(ilevel)

  dener=0.0d0
  max_ener=0.0d0
  min_ener=1.0d30
  error_cg_loc=0.0d0
  do i=1,nb_ind
     this = liste_ind(i)
     max_ener=max(max_ener, unew(liste_ind(i),ind_cr))
     min_ener=min(min_ener, unew(liste_ind(i),ind_cr))
     error_cg_loc = error_cg_loc + unew(liste_ind(i),ind_cr)
  end do

  ! Compute maximum error
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(max_ener,max_ener_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  max_ener=max_ener_all
  call MPI_ALLREDUCE(min_ener,min_ener_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  min_ener=min_ener_all
  call MPI_ALLREDUCE(error_cg_loc,error_cg_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  error_cg_loc=error_cg_all
#endif
  dener=max_ener/min_ener

  nsub_imp=1
!!$  if(maxval(dener) .gt. 1.d4)then
!!$     nsub_imp=int(maxval(dener))
!!$     dt_imp=dt_imp/real(nsub_imp)
!!$  endif
  if(myid==1)print*,'ilevel',ilevel,'MAXIMUM OF DENER=',dener,'NSUB_IMP=',nsub_imp,error_cg_loc

  do isub=1,nsub_imp

!!$     !===================================================
!!$     ! Compute CR diffusion coefficient :
!!$
!!$     ! Room left for Spitzer relaxation term between ion and electron in enew(indcell(i)) 
!!$     !===================================================
!!$     do i=1,nb_ind
!!$        this = liste_ind(i)
!!$        temp = unew(this,ind_Tenew)
!!$        ! Compute Rosseland opacity (Compute kappa*rho)
!!$        divu(this)= Dcr / scale_kappa       
!!$     end do

     ! Update boundaries
!!$     call make_virtual_fine_dp(unew(1,2),ilevel)
!!$     call make_virtual_fine_dp(unew(1,3),ilevel)
!!$     call make_virtual_fine_dp(unew(1,1),ilevel)
!!$     call make_virtual_fine_dp(unew(1,4),ilevel)
     call make_virtual_fine_dp(unew(1,ind_cr),ilevel)

     call make_boundary_crdiffusion(ilevel,ind_cr)

     !===================================================
     ! Compute r1 = b1 - A1x1 and store it into unew(i,1)
     ! Also set p1 = r1 and store it into unew(i,2)
     !===================================================
     call cmp_matrixA_crdiffusion (ilevel,Nsub, 1,ind_cr)

     !        call make_virtual_reverse_dp(unew(1,1),ilevel)
     call make_virtual_fine_dp(unew(1,1),ilevel)

     !        call make_virtual_reverse_dp(unew(1,2),ilevel)
     call make_virtual_fine_dp(unew(1,2),ilevel)

     !===================================
     ! Compute right-hand side norm (max)
     !===================================
     call dot_product(unew(:,1),unew(:,1),rhs_norm1,final_sum)

     !===================================================
     ! Compute Preconditionner M=1/diag(A) and store it in unew(i,4)
     !===================================================
     !ben
     call make_boundary_crdiffusion(ilevel,ind_cr)
     call cmp_matrixA_crdiffusion (ilevel,Nsub, 3,ind_cr)

     !        call make_virtual_reverse_dp(unew(1,4),ilevel)
     call make_virtual_fine_dp(unew(1,4),ilevel)
     !====================================
     ! MAIN ITERATION LOOP
     !====================================     

     iter=0; itermax=5000

     error_ini=sqrt(rhs_norm1)

     !     error_ini=sqrt(real(final_sum))
     error=error_ini
     error_cg_loc=1.0d0

     do while(error/error_ini>epsilon_diff .and.iter<itermax)
        !     do while(error_cg_loc .gt. epsilon)

        iter=iter+1

        !============================================
        ! Compute z = Mr and store it into unew(i,3)
        !============================================

        do i=1,nb_ind
           this = liste_ind(i)
           unew(this,3) = unew(this,1) * unew(this,4)
        end do

        !====================================
        ! Compute scalar r.z
        !====================================

        call dot_product(unew(:,1),unew(:,3),r2,final_sum)
        r2=r2!real(final_sum)
        !        r2=real(final_sum)

        !====================================
        ! Compute beta factor
        !====================================

        if(iter==1)then
           beta_cg = 0.0d0
        else
           beta_cg = r2/r2_old
        end if
        r2_old=r2

        !====================================
        ! Recurrence on p = z + beta*p
        !====================================
        call cX_plus_Y_to_Z (beta_cg,unew(:,2),unew(:,3),unew(:,2))
        ! Update boundaries
        call make_boundary_crdiffusion(ilevel,ind_cr)
        call make_virtual_fine_dp(unew(1,2),ilevel)

        !==============================================
        ! Compute q1 = Ap1 and store it into unew(i,3)
        !==============================================
        call cmp_matrixA_crdiffusion (ilevel,Nsub, 2,ind_cr)

        !        call make_virtual_reverse_dp(unew(1,3),ilevel)
        call make_virtual_fine_dp(unew(1,3),ilevel)

        !====================================
        ! Compute p.Ap scalar product
        !====================================

        call dot_product(unew(:,2),unew(:,3),pAp,final_sum)
        pap = pap!real(final_sum) !DDP
        !        pap = real(final_sum) !DDP

        !====================================
        ! Compute alpha factor
        !====================================
        alpha_cg = r2/pAp

        !====================================
        ! Recurrence on x = x + alpha*p
        !====================================
        error_cg_loc=0.0d0
        do i=1,nb_ind
           this = liste_ind(i)
           !        unew(liste_ind(i),8+igroup) = max(unew(liste_ind(i),8+igroup),eray_min/scale_E0)
           error_cg_loc=max(error_cg_loc, abs((alpha_cg*unew(liste_ind(i),2))/unew(liste_ind(i),ind_cr)))
           !error_cg_loc=error_cg_loc+  unew(liste_ind(i),ind_cr)*(0.5**ilevel)**3
        end do
        ! Compute maximum error
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(error_cg_loc,error_cg_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        error_cg_loc=error_cg_all
#endif

        call cX_plus_Y_to_Z (alpha_cg,unew(:,2),unew(:,ind_cr),unew(:,ind_cr))

        !====================================
        ! Recurrence on r (unew(i,1))
        !====================================
        call cX_plus_Y_to_Z (- alpha_cg ,unew(:,3),unew(:,1),unew(:,1))

        !===================================
        ! Compute right-hand side norm (max)
        !===================================
        call dot_product(unew(:,1),unew(:,1),rhs_norm1,final_sum)

        error=SQRT(rhs_norm1)
        !        error=SQRT(real(final_sum))
        !        error = error_cg_loc

        if(verbose)write(*,112)iter,error,error/error_ini,error_cg_loc,r2,pap

     end do
     ! End main iteration loop

     if(iter >= itermax)then
        if(myid==1)write(*,*)'Diffusion fail to converge...'
     end if

!!$     !====================================
!!$     ! Copie des flux
!!$     !====================================
!!$     call cmp_matrixA (ilevel,Nsub, 4,igroup)

     if(myid==1) write(*,*)' CG :',iter, ' error=',error/error_ini,error_ini,error_cg_loc,error
     niter=niter+iter     

  end do
  !ENd loop over subcycles
end do

!=============================
! Update energy value
!=============================
call cmp_energy_cr(2)

! Update boundaries
call make_virtual_fine_dp(uold(1,5),ilevel)
do ivar=1,ncr
  call make_virtual_fine_dp(uold(1,ind_cr),ilevel)
end do
!  call upload_fine(ilevel)
!  call upload_fine(ilevel-1)


111 format('   Entering conduction_cg')
112 format('   ==> Step=',i5,' Error=',5(1pe10.3,1x))!,e23.15)
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))

deallocate(liste_ind)

contains

  !###########################################################
  !###########################################################

subroutine dot_product(fact1,fact2,pdtvect,local_sum) !                pdtvect = sum(fact1*fact2)
 implicit none
 real(dp),dimension(1:ncoarse+twotondim*ngridmax),intent(IN)::fact1,fact2
 real(dp),intent(OUT)::pdtvect
 complex*16,intent(OUT)::local_sum

 real(dp)::pdtvect_all
 complex*16 ::global_sum
 integer::this

 pdtvect=0.0d0
 local_sum = cmplx(0.0d0,0.0d0)
 global_sum = cmplx(0.0d0,0.0d0)

 do i=1,nb_ind
    this = liste_ind(i)
    !call DDPDD (cmplx(fact1(this)*fact2(this), 0.0,dp), local_sum, 1, itype)
    pdtvect = pdtvect + fact1(this)*fact2(this)
 end do

 ! Compute global norms
#ifndef WITHOUTMPI
 call MPI_ALLREDUCE(pdtvect,pdtvect_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
 pdtvect   = pdtvect_all
 !call MPI_ALLREDUCE(local_sum,global_sum,1,MPI_COMPLEX,MPI_SUMDD,MPI_COMM_WORLD,info)
 !local_sum = global_sum
 local_sum = pdtvect
#endif

end subroutine dot_product

!###########################################################
!###########################################################
subroutine cX_plus_Y_to_Z (cste,vectX,vectY,vectZ)! vectZ = cste*vectX+vectY
 implicit none
 real(dp),dimension(1:ncoarse+twotondim*ngridmax),intent(IN)::vectX,vectY
 real(dp),intent(IN)::cste
 real(dp),dimension(1:ncoarse+twotondim*ngridmax),intent(OUT)::vectZ


 do i=1,nb_ind
    vectZ(liste_ind(i)) = vectY(liste_ind(i)) + cste*vectX(liste_ind(i)) 
 end do

end subroutine cX_plus_Y_to_Z


end subroutine crdiffusion_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_matrixA_crdiffusion (ilevel,Nsub,compute,igroup)
  !------------------------------------------------------------------
  ! This routine computes the matrix A to vect_in and create vect_out
  ! compute = 1 : residu           	return B - Ax
  ! compute = 2 : Produit                 return  A.p
  ! compute = 3 : Preconditionner         return diag(A)
  ! compute = 4 : Compute flux in rad_flux
  !------------------------------------------------------------------

use amr_commons
use hydro_commons
implicit none

integer,intent(IN)::compute,ilevel,Nsub,igroup

integer ,dimension(1:nvector,0:2*ndim),save::  igridn
integer ,dimension(1:nvector),save ::          ind_cell , ind_grid

integer :: i,idim,ind,igrid,ngrid,ncache,iskip

if(numbtot(1,ilevel)==0)return

! Loop over active grids by vector sweeps
ncache=active(ilevel)%ngrid
do igrid=1,ncache,nvector
  ngrid=MIN(nvector,ncache-igrid+1)
  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
  end do
  call crdifffine1(ind_grid,ngrid,ilevel,compute,igroup)
end do

end subroutine cmp_matrixA_crdiffusion
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine make_boundary_crdiffusion(ilevel,ind_cr)
use amr_commons
use hydro_commons
use cooling_module,ONLY:kB,mH,clight
!  use radiation_parameters
implicit none
! -------------------------------------------------------------------
! This routine set up boundary conditions for fine levels.
! -------------------------------------------------------------------
integer,intent(IN)::ilevel,ind_cr
integer::ibound,boundary_dir,idim,inbor
integer::i,ncache,ivar,igrid,ngrid,ind,iperp1,iperp2,iplane,icell
integer::iskip,iskip_ref,nx_loc,ix,iy,iz,igrp,gdim
integer,dimension(1:8)::ind_ref,alt
integer,dimension(1:2,1:4)::ind0
integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

real(dp)::dx,dx_loc,scale,switch
real(dp)::rosseland_ana
real(dp),dimension(1:3)::skip_loc,gs
real(dp),dimension(1:twotondim,1:3)::xc
real(dp),dimension(1:nvector,1:ndim),save::xx
real(dp),dimension(1:nvector,1:nvar+3),save::uu
real(dp)::dd,t2,usquare,emag,erad_loc,eps,ekin,Cv,rho
real(dp)::scale_kappa
real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
real(dp)::a_R,radiation_Source,kpara_ana

if(.not. simple_boundary)return

! Conversion factor from user units to cgs units
call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

scale_kappa=scale_d*scale_l*scale_v**3

! Mesh size at level ilevel
dx=0.5D0**ilevel

! Rescaling factors
nx_loc=(icoarse_max-icoarse_min+1)
skip_loc=(/0.0d0,0.0d0,0.0d0/)
if(ndim>0)skip_loc(1)=dble(icoarse_min)
if(ndim>1)skip_loc(2)=dble(jcoarse_min)
if(ndim>2)skip_loc(3)=dble(kcoarse_min)
scale=boxlen/dble(nx_loc)
dx_loc=dx*scale

! Set position of cell centers relative to grid center
do ind=1,twotondim
  iz=(ind-1)/4
  iy=(ind-1-4*iz)/2
  ix=(ind-1-2*iy-4*iz)
  if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
  if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
  if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
end do

! Loop over boundaries
do ibound=1,nboundary
   ! Compute direction of reference neighbors
  boundary_dir=boundary_type(ibound)-10*(boundary_type(ibound)/10)
  if(boundary_dir==1)inbor=2

  if(boundary_dir==2)inbor=1
  if(boundary_dir==3)inbor=4
  if(boundary_dir==4)inbor=3
  if(boundary_dir==5)inbor=6
  if(boundary_dir==6)inbor=5

  ! Compute index of reference cells
  ! Zero flux
  if(boundary_type(ibound)== 1)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
  if(boundary_type(ibound)== 2)ind_ref(1:8)=(/2,1,4,3,6,5,8,7/)
  if(boundary_type(ibound)== 3)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
  if(boundary_type(ibound)== 4)ind_ref(1:8)=(/3,4,1,2,7,8,5,6/)
  if(boundary_type(ibound)== 5)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
  if(boundary_type(ibound)== 6)ind_ref(1:8)=(/5,6,7,8,1,2,3,4/)
  ! Zero flux
  if(boundary_type(ibound)==11)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
  if(boundary_type(ibound)==12)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
  if(boundary_type(ibound)==13)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
  if(boundary_type(ibound)==14)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
  if(boundary_type(ibound)==15)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
  if(boundary_type(ibound)==16)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)
  ! Imposed boundary
  if(boundary_type(ibound)==21)ind_ref(1:8)=(/1,1,3,3,5,5,7,7/)
  if(boundary_type(ibound)==22)ind_ref(1:8)=(/2,2,4,4,6,6,8,8/)
  if(boundary_type(ibound)==23)ind_ref(1:8)=(/1,2,1,2,5,6,5,6/)
  if(boundary_type(ibound)==24)ind_ref(1:8)=(/3,4,3,4,7,8,7,8/)
  if(boundary_type(ibound)==25)ind_ref(1:8)=(/1,2,3,4,1,2,3,4/)
  if(boundary_type(ibound)==26)ind_ref(1:8)=(/5,6,7,8,5,6,7,8/)
  ! For magnetic field, we have only free boundary amd imposed boundary
  if(boundary_dir==1)alt(1:8)=-(/2,1,2,1,2,1,2,1/)
  if(boundary_dir==2)alt(1:8)=+(/1,2,1,2,1,2,1,2/)
  if(boundary_dir==3)alt(1:8)=-(/1,1,2,2,1,1,2,2/)
  if(boundary_dir==4)alt(1:8)=+(/2,2,1,1,2,2,1,1/)
  if(boundary_dir==5)alt(1:8)=-(/1,1,1,1,2,2,2,2/)
  if(boundary_dir==6)alt(1:8)=+(/2,2,2,2,1,1,1,1/)

  ! Velocity sign switch for reflexive boundary conditions
  gs=(/1,1,1/)
  if(boundary_type(ibound)==1.or.boundary_type(ibound)==2)gs(1)=-1
  if(boundary_type(ibound)==3.or.boundary_type(ibound)==4)gs(2)=-1
  if(boundary_type(ibound)==5.or.boundary_type(ibound)==6)gs(3)=-1

  ! Direction of gravity vector for hydrostatic equilibrium
  if(boundary_dir==1.or.boundary_dir==2)gdim=1
  if(boundary_dir==3.or.boundary_dir==4)gdim=2
  if(boundary_dir==5.or.boundary_dir==6)gdim=3

  if(boundary_dir==1)ind0(1:2,1:4)=RESHAPE((/2,4,6,8,1,3,5,7/),SHAPE=(/2, 4/))
  if(boundary_dir==2)ind0(1:2,1:4)=RESHAPE((/1,3,5,7,2,4,6,8/),SHAPE=(/2, 4/))
  if(boundary_dir==3)ind0(1:2,1:4)=RESHAPE((/3,4,7,8,1,2,5,6/),SHAPE=(/2, 4/))
  if(boundary_dir==4)ind0(1:2,1:4)=RESHAPE((/1,2,5,6,3,4,7,8/),SHAPE=(/2, 4/))
  if(boundary_dir==5)ind0(1:2,1:4)=RESHAPE((/5,6,7,8,1,2,3,4/),SHAPE=(/2, 4/))
  if(boundary_dir==6)ind0(1:2,1:4)=RESHAPE((/1,2,3,4,5,6,7,8/),SHAPE=(/2, 4/))

  if(boundary_dir==1)then
     iperp1=6; iperp2=nvar+1
  endif
  if(boundary_dir==2)then
     iperp1=nvar+1; iperp2=6
  endif
  if(boundary_dir==3)then
     iperp1=7; iperp2=nvar+2
  endif
  if(boundary_dir==4)then
     iperp1=nvar+2; iperp2=7
  endif
  if(boundary_dir==5)then
     iperp1=8; iperp2=nvar+3
  endif
  if(boundary_dir==6)then
     iperp1=nvar+3; iperp2=8
  endif

  ! Loop over grids by vector sweeps
  ncache=boundary(ibound,ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=boundary(ibound,ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring reference grid
     do i=1,ngrid
        ind_grid_ref(i)=son(nbor(ind_grid(i),inbor))
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather neighboring reference cell
        iskip_ref=ncoarse+(ind_ref(ind)-1)*ngridmax
        do i=1,ngrid
           ind_cell_ref(i)=iskip_ref+ind_grid_ref(i)
        end do

        ! Zero flux boundary conditions
        if((boundary_type(ibound)/10).ne.2)then

           ! Gather reference hydro variables
           do ivar=1,nvar+3
              do i=1,ngrid
                 uu(i,ivar)=uold(ind_cell_ref(i),ivar)
              end do
           end do
           do i=1,ngrid
              uold(ind_cell(i),ind_cr)=uold(ind_cell_ref(i),ind_cr)
              unew(ind_cell(i),2)     = unew(ind_cell_ref(i),2)
           end do

           ! Imposed boundary conditions
        else

           ! Compute cell center in code units and rescale position from code units to user units
           do idim=1,ndim
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0)then
                    xx(i,idim)=(xg(ind_grid(i),idim)+xc(ind,idim)-skip_loc(idim))*scale
                 end if
              end do
           end do

           call boundana(xx,uu,dx_loc,ibound,ngrid)

           ! Scatter variables
           do i=1,ngrid 
              if(son(ind_cell(i)) == 0)then
                 unew(ind_cell(i),2)=  0.0d0
                 uold(ind_cell(i),ind_cr)=  uu(i,ind_cr)
              end if
           end do
        end if

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end do
! End loop over boundaries


111 format('   Entering make_boundary_conduction for level ',I2)

end subroutine make_boundary_crdiffusion
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine cmp_energy_cr(Etype)
use amr_commons
use hydro_commons
use cooling_module,ONLY:kB,mH,clight
!  use radiation_parameters
implicit none
integer,intent(in) :: Etype ! Etype=1 : beginning ; Etype=2 : end
integer ::i,idim,this,mvar,igroup
real(dp)::usquare,Cv,eps,rho,ecr
real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

! EOS
real(dp) :: dd,ee
integer ::ind_cr

! Conversion factor from user units to cgs units
call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

ind_cr=8
if(twotemp)ind_cr=9

do i=1,nb_ind
  this = liste_ind(i)
  rho   = uold(this,1)

  ! Compute total cosmic ray energy
  ecr=0.0d0
  do igroup=1,ncr
     ecr = ecr + uold(this,ind_cr+igroup)
  end do

  if(Etype==1)then
     unew(this,5)=uold(this,5)-ecr ! Store old total energy minus CR energy
     do igroup=1,ncr
        unew(this,ind_cr+igroup) = uold(this,ind_cr+igroup)
!        print*,unew(this,ind_cr+igroup),'1'
     end do

  elseif(Etype==2)then
     ! update total energy
     uold(this,5) = unew(this,5)
     do igroup=1,ncr
        uold(this,ind_cr+igroup) = unew(this,ind_cr+igroup)
        uold(this,5) = uold(this,5) + unew(this,ind_cr+igroup)
!        print*,unew(this,ind_cr+igroup),'2'
     end do

  end if
end do



end subroutine cmp_energy_cr

!################################################################
!################################################################
!################################################################ 
!################################################################
