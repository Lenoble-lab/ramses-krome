subroutine crdiffusion_cg (ilevel)
  use amr_commons
  use hydro_commons
  use constants,ONLY:kB,mH 
  !  use radiation_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !===========================================================================
  ! Y. Dubois (05/15/2020): STS with RKL2 (adapted from Meyer+12)
  ! Y. Dubois (05/14/2020): Automatic switch to use an explcit or an implicit
  ! solver depending on the ratio of time steps
  ! Y. Dubois (06/2019)   : BiCGSTAB with slope limiter
  ! (adapted from CG from Dubois & Commerçon (2016))
  ! Iterative solver with Preconditionned Bi-Conjugate Gradient Stabilized 
  ! method to solve A x = b, where b is the vector of olf variables, x is 
  ! the vector of new (searched) variables, and A is the matrix containing the
  ! the diffusion coefficients. M is the preconditionning matrix, which is set
  ! to the identity matrix at the moment (i.e. the BiCGSTAB is not
  ! preconditionned).
  ! The corresponding algorithm is:
  ! 
  ! r_0=b-A*x_0
  ! rho_0=alpha=omega_0=1
  ! v_0=p_0=0
  ! i=i+1
  ! do i=1,2,3...
  !    rho_i=(r_0,r_{i-1})
  !    beta=(rho_i/rho_{i-1})*(alpha/omega_{i-1})
  !    p_i=r_{i-1}+beta*(p_i-omega_{i-1}*v_{i-1})
  !    y=M^-1*p_i
  !    v_i=A*y
  !    alpha=rho_i/(r_zero,v_i)
  !    s=r_{i-1}-alpha*v_{i}
  !    if(||s|| < tol) then
  !       x_i=x_{i-1}+alpha*p_i
  !       exit loop
  !    z=M^-1*s
  !    t=A*z
  !    omega_i=(t,s)/(t,t)
  !    r_i=s-omega_i*t
  !    x_i=x_{i-1}+alpha*y+omega_i*z
  !    if(|r_i|/|r_0| < tol) exit loop
  ! enddo
  !
  !===========================================================================
  !   ri      : stored in unew(i,1)
  !   pi      : stored in unew(i,2)
  !   y       : stored in unew(i,3)
  !   M^-1    : stored in unew(i,4) (code structure imposes diag preconditionner)
  !   x(i)    : stored in unew(i,8+nrad+ntp+ncr)(new energy)
  !   vi      : stored in unew(i,6)
  !   s       : stored in unew(i,7)
  !   z       : stored in unew(i,8)
  !   t       : stored in unew(i,nvar+1)
  !   rzero   : stored in unew(i,nvar+2)
  !   M^-1*t  : stored in unew(i,nvar+3)
  !  Do NOT use unew(i,5): placeholder for total energy minus cr (cmp_energy)
  !===========================================================================
  integer,intent(IN)::ilevel
  complex*16 :: final_sum
  real(dp)::error,error_ini,error_cg_loc,error_cg_all
  real(dp)::alpha_cg,beta_cg,dt_exp,error_s_ini
  real(dp)::r2,pAp,rhs_norm1,rhs_norm2
  integer ::i,info,ind,iter,iskip,itermax
  integer ::this,nx_loc,nsub_imp,isub
  real(dp)::dx,dx_loc,surf_loc,scale,vol_loc
  real(dp)::min_ener,min_ener_all,max_ener,max_ener_all
  real(dp)::dener,max_loc,max_loc_all
  logical::exist_leaf_cell=.true.
  integer ::ind_cr,ivar,irad
  logical::continue_loop
  real(dp),dimension(2)::omega_cg,rho_cg
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_kappa
  real(dp)::dt_cond,ratio_timestep
  integer::s_rkl2

  if(verbose)write(*,111)
  if(numbtot(1,ilevel)==0)return

  ! Rescaling factors
  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  surf_loc = dx_loc**(ndim-1)
  vol_loc  = dx_loc**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kappa = scale_l**2/scale_t

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

  if (nb_ind == 0)then
     !print*,'No leaf-cell - myid=',myid
     exist_leaf_cell=.false.
  end if

  do i=1,nb_ind
     this = liste_ind(i)
     unew(this,1:nvar+3)=0.0d0
     divu(this)=0.0d0
  end do

  !===================================================================
  ! Compute gas temperature stored in unew(i,ind_cr+1)
  ! ind_cr=9  (if twotemp=.false.)
  ! ind_cr=10 (if twotemp=.true.)
  ! 1: Etot     -> Etot-Ecr
  ! 2: Etot-Ecr -> Etot
  !===================================================================
  call cmp_energy_cr(1)

  do irad=1,ncr
     ind_cr = 8+irad
     if(twotemp)ind_cr = 9+irad

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
     min_ener=min(min_ener, MAX(unew(liste_ind(i),ind_cr),smallcr))
     error_cg_loc = error_cg_loc + unew(liste_ind(i),ind_cr)
     divu(this)= unew(liste_ind(i),ind_cr)
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

  !===================================================
  ! Compute Diffusion explicit time step : only works
  ! with constant diffusion coef, not streaming or 
  ! varying efficiencies diff coef which are computed
  ! later in crdiffusion_fine.f90
  !===================================================
  if(explicit_diff.or.autoswitch_diff.or.sts_diff)then
     dt_cond=0.5d0*dx_loc**2/(Dcr/scale_kappa)/dble(ndim)
     ratio_timestep=dt_cond/dtnew(ilevel)
     if(ratio_timestep<1.0)then
        nsub_imp=int(1d0/ratio_timestep)+1
        s_rkl2=int(0.5d0*(-1d0+sqrt(1d0+4d0*(4d0/ratio_timestep+1d0))))+1
     else
        nsub_imp=1
        s_rkl2=1
     endif
  else
     nsub_imp=1
  endif
  if(myid==1)write(*,'(A,1pe12.3,A,1pe12.3,A,I5,A,1pe12.3)')'dt_diff:',dt_cond,',ratio_dt:',ratio_timestep,',nstep_diff:',nsub_imp,',dt_imp:',dt_imp/dble(nsub_imp)
  if(autoswitch_diff)then
     if(nsub_imp.le.nmax_explicit)then
        explicit_diff=.true.
     else
        explicit_diff=.false. !just to be sure (but it's redundant)
        nsub_imp=1
     endif
  endif
  dt_imp=dt_imp/real(nsub_imp)
  if(autoswitch_diff)then
     if(myid==1)then
        if(explicit_diff)then
           write(*,*)'Doing ',nsub_imp,' subcycled explicit steps for CR diffusion'
        else
           write(*,*)'Doing ',nsub_imp,' subcycled implicit steps for CR diffusion'
        endif
     endif
  endif
  if(sts_diff.and.myid==1)then
     write(*,*)'It should be doing ',s_rkl2,' STS sub-time steps for half a RKL-2 step'
!!$     write(*,*)'mu =',mu_rkl2(2),mu_rkl2(3)
!!$     write(*,*)'mut=',mut_rkl2(1,3),mut_rkl2(2,3),mut_rkl2(3,3)
!!$     write(*,*)'nu =',nu_rkl2(2),nu_rkl2(3)
!!$     write(*,*)'gt =',gt_rkl2(2,3),gt_rkl2(3,3)
     if(s_rkl2>3)then
!!$     dt_imp=0.5d0*dt_exp ! Choose to use the hydro time step (s is set up according to this choice)
!!$                         ! And half that because this is Runge-Kutta 2
     endif
  endif

  if(myid==1)write(*,'(A,I4,A,e12.3,A,I4,E12.3)') '(CRdiff) ilevel',ilevel,' MAXIMUM OF DENER=',dener,' NSUB_IMP=',nsub_imp,error_cg_loc

  do isub=1,nsub_imp

     do ivar=1,nvar+3
        call make_virtual_fine_dp(unew(1,ivar),ilevel)
     enddo
     call make_virtual_fine_dp(divu(1),ilevel)

     ! Update physical boundaries
     call make_boundary_crdiffusion(ilevel,ind_cr)
     if(explicit_diff)then

        if(sts_diff.and.s_rkl2>3)then
           !Tricky part...
           call cmp_matrixA_crdiffusion(ilevel,0,ind_cr)


        else
           call cmp_matrixA_crdiffusion(ilevel,0,ind_cr)
           call cmp_energy_cr(2)
           ! Update boundaries
           call make_virtual_fine_dp(uold(1,5     ),ilevel)
           call make_virtual_fine_dp(uold(1,ind_cr),ilevel)
        endif

     else

     !===================================================
     ! Compute r1 = b1 - A1x1 and store it into unew(i,1)
     ! Also set rzero = r1    and store it into unew(i,nvar+2)
     !===================================================
     call cmp_matrixA_crdiffusion (ilevel,1,ind_cr)
     ! Are not involved into matrix multiplication...
!!$     call make_virtual_fine_dp(unew(1,1),ilevel)
!!$     call make_virtual_fine_dp(unew(1,nvar+2),ilevel)

     !===================================
     ! Compute right-hand side norm (max)
     !===================================
     call dot_product(unew(:,1),unew(:,1),rhs_norm1,final_sum)
     if(debug_bicg.and.myid==1)write(*,*)'r_0',sqrt(rhs_norm1)

     !===================================================
     ! Compute Preconditionner M=1/diag(A) and store it in unew(i,4)
     !===================================================
     call make_boundary_crdiffusion(ilevel,ind_cr)
     call cmp_matrixA_crdiffusion (ilevel,3,ind_cr)
     call make_virtual_fine_dp(unew(1,4),ilevel)

     !====================================
     ! MAIN ITERATION LOOP
     !====================================     

     iter=0; itermax=5000

     error_ini=sqrt(rhs_norm1)

     !     error_ini=sqrt(real(final_sum))
     error=1.d2*error_ini
     error_cg_loc=1.0d0

     rho_cg  (2)=1.0d0
     alpha_cg   =1.0d0
     omega_cg(2)=1.0d0
     continue_loop=.true.
if(dener>1)then
!!$     do while((error/error_ini>epsilon_diff .or. error_cg_loc>epsilon_diff) .and.iter<itermax)
     do while(continue_loop)

        iter=iter+1

        !=======================================================
        ! Update old state with previous new state values
        !=======================================================
        rho_cg  (1)=rho_cg  (2)
        omega_cg(1)=omega_cg(2)

        !=======================================================
        ! rho_i = (r_0,r_{i-1})
        !=======================================================
        call dot_product(unew(:,1),unew(:,nvar+2),r2,final_sum)
        rho_cg(2)=r2
        if(debug_bicg.and.myid==1)write(*,*)'================================================='
        if(debug_bicg.and.myid==1)write(*,'(A15,es10.2)')'rho_i',r2
        
        !=======================================================
        ! beta=(rho_i/rho_{i-1})*(alpha/omega_{i-1})
        !=======================================================
        beta_cg=rho_cg(2)/rho_cg(1)*alpha_cg/omega_cg(1)
        if(debug_bicg.and.myid==1)write(*,'(A15,es10.2)')'beta_cg',beta_cg
        if(debug_bicg.and.myid==1)write(*,'(A15,es10.2)')'omega_{i-1}',omega_cg(1)

        !=======================================================
        ! p_i = r_{i-1}+beta*(p_{i-1}-omega_{i-1}*v_{i-1})
        !=======================================================
!!$        do i=1,nb_ind
!!$           this=liste_ind(i)
!!$           write(*,'(A15,es10.2)')'p_{i-1}(this)',unew(this,2)
!!$        enddo
!!$        do i=1,nb_ind
!!$           this=liste_ind(i)
!!$           write(*,'(A15,es10.2)')'v_{i-1}(this)',unew(this,6)
!!$        enddo
!!$        do i=1,nb_ind
!!$           this=liste_ind(i)
!!$           write(*,'(A15,es10.2)')'r_{i-1}(this)',unew(this,1)
!!$        enddo

        call cX_plus_Y_to_Z(beta_cg,unew(:,2)-omega_cg(1)*unew(:,6),unew(:,1),unew(:,2))
!!$        call cX_plus_Y_to_Z( beta_cg            ,unew(:,2),unew(:,1),unew(:,2))
!!$        call cX_plus_Y_to_Z(-beta_cg*omega_cg(1),unew(:,6),unew(:,2),unew(:,2))

        call make_boundary_crdiffusion(ilevel,ind_cr)
        call make_virtual_fine_dp(unew(1,2),ilevel) ! When M will be non-diagonal (one day) it will be important
!!$        do i=1,nb_ind
!!$           this=liste_ind(i)
!!$           write(*,'(A15,es10.2)')'p(this)',unew(this,2)
!!$        enddo

        !=======================================================
        ! Compute y = M^-1*p_i and store it into unew(i,3)
        ! (here it works because diagonal matrix M)
        !=======================================================
        do i=1,nb_ind
           this = liste_ind(i)
           unew(this,3) = unew(this,4) * unew(this,2)
        end do
        call make_virtual_fine_dp(unew(1,3),ilevel)
        call make_boundary_crdiffusion(ilevel,ind_cr) !YD: deal with boundaries!!! (update boundaries unew(3))

        !=======================================================
        ! Compute v_i=Ay and store it into unew(i,6)
        !=======================================================
        call cmp_matrixA_crdiffusion(ilevel,2,ind_cr)
        call make_virtual_fine_dp(unew(1,6),ilevel) 
!!$        do i=1,nb_ind
!!$           this=liste_ind(i)
!!$           write(*,'(A15,es10.2)')'v(this)',unew(this,6)
!!$        enddo
       
        !=======================================================
        ! Compute (rzero,v_i)
        !=======================================================
        call dot_product(unew(:,nvar+2),unew(:,6),r2,final_sum)
        if(debug_bicg.and.myid==1)write(*,'(A15,es10.2)')'(rzero,v_i)',r2
        !=======================================================
        ! Compute alpha=rho_i/(rzero,v_i)
        !=======================================================
        alpha_cg=rho_cg(2)/r2
        if(debug_bicg.and.myid==1)write(*,'(A15,es10.2)')'alpha_cg',alpha_cg

        !=======================================================
        ! s=r_{i-1}-alpha*v_i
        ! unew(i,1): r_{i-1}
        ! unew(i,6): v_i
        ! unew(i,7): s
        !=======================================================
        call cX_plus_Y_to_Z(-alpha_cg,unew(:,6),unew(:,1),unew(:,7))
        call make_virtual_fine_dp(unew(1,7),ilevel) ! When M will be non-diagonal (one day) it will be important
!!$        do i=1,nb_ind
!!$           this=liste_ind(i)
!!$           write(*,'(A15,es10.2)')'s(this)',unew(this,7)
!!$        enddo
        

        ! This test can be done (if done properly and have to check here)
        ! but unecessary CPU-wise for N iter steps, 2N calculations done
        ! and this check might save 1 calculation
        !=======================================================
        ! Compute (s,s)
        !=======================================================
        call dot_product(unew(:,7),unew(:,7),rhs_norm1,final_sum)
        error=SQRT(rhs_norm1)
        if(iter==1)error_s_ini=error
        if(verbose)write(*,112)iter,error,error/error_ini,error_cg_loc,r2,pap
        if(debug_bicg.and.myid==1)write(*,'(A15,2es10.2)')'|s|',error,error/error_s_ini

        if(error/error_ini .le. epsilon_diff)then

           continue_loop=.false.
           !=======================================================
           ! Compute x_i=x_{i-1}+alpha*p_i
           ! unew(i,2     ): p_i
           ! unew(i,ind_cr): x_{i-1}->x_i
           !=======================================================
           call cX_plus_Y_to_Z(alpha_cg,unew(:,2),unew(:,ind_cr),unew(:,ind_cr))           
           
        else
           
        !=======================================================
        ! Compute z = M^-1*s and store it into unew(i,8)
        ! (here it works because diagonal matrix M)
        !=======================================================
        do i=1,nb_ind
           this = liste_ind(i)
           unew(this,8) = unew(this,4) * unew(this,7)
        end do
        call make_virtual_fine_dp(unew(1,8),ilevel)
        call make_boundary_crdiffusion(ilevel,ind_cr) !YD: deal with boundaries!!! (update boundaries unew(8))
        
        !=======================================================
        ! Compute t=Az and store it into unew(i,nvar+1)
        !=======================================================
!!$        call make_boundary_crdiffusion(ilevel,ind_cr) !YD: deal with boundaries!!!
        call cmp_matrixA_crdiffusion(ilevel,4,ind_cr)
        call make_virtual_fine_dp(unew(1,nvar+1),ilevel)
        
        !=======================================================
        ! Compute tprime=M^-1*t and store it into unew(i,nvar+3)
        !=======================================================
        do i=1,nb_ind
           this = liste_ind(i)
           unew(this,nvar+3) = unew(this,4) * unew(this,nvar+1)
        end do
        call make_virtual_fine_dp(unew(1,nvar+3),ilevel)

        !=======================================================
        ! Compute omega_i=(M^-1*t,M^-1*s)/(M^-1*t,M^-1*t)
        !=======================================================
        call dot_product(unew(:,nvar+3),unew(:,8     ),rhs_norm1,final_sum)
        call dot_product(unew(:,nvar+3),unew(:,nvar+3),rhs_norm2,final_sum)
        if(rhs_norm2.ne.0.0d0)then
           omega_cg(2)=rhs_norm1/rhs_norm2
        else
           ! WARNING: this might just help the code to produce NaN more easily...
           ! Do something different
           omega_cg(2)=0.0d0
        endif
!!$           do i=1,nb_ind
!!$              this=liste_ind(i)
!!$              write(*,'(A15,es10.2)')'t(this)',unew(this,nvar+1)
!!$           enddo
        if(debug_bicg.and.myid==1)write(*,'(A15,es10.2)')'(t,s)',rhs_norm1
        if(debug_bicg.and.myid==1)write(*,'(A15,es10.2)')'(t,t)',rhs_norm2
        if(debug_bicg.and.myid==1)write(*,'(A15,es10.2)')'omega',omega_cg(2)
        
        !=======================================================
        ! Compute x_i=x_{i-1}+alpha*y+omega_i*z
        !=======================================================
        ! unew(i,ind_cr): x_{i-1}->x_i
        ! unew(i,3): y
        ! unew(i,8): z
        !=======================================================
        call cX_plus_Y_to_Z(alpha_cg   ,unew(:,3),unew(:,ind_cr),unew(:,ind_cr))
        call cX_plus_Y_to_Z(omega_cg(2),unew(:,8),unew(:,ind_cr),unew(:,ind_cr))
        ! call make_boundary_crdiffusion(ilevel,ind_cr)
        call make_virtual_fine_dp(unew(1,ind_cr),ilevel)
        
        !=======================================================
        ! Compute r_i=s-omega_i*t
        !=======================================================
        ! unew(i,1     ): r_i
        ! unew(i,7     ): s
        ! unew(i,nvar+1): t
        !=======================================================
        call cX_plus_Y_to_Z(-omega_cg(2),unew(:,nvar+1),unew(:,7),unew(:,1))
        ! call make_boundary_crdiffusion(ilevel,ind_cr)
        
        !=======================================================
        ! Compute (r_i,r_i)
        !=======================================================
        call dot_product(unew(:,1),unew(:,1),rhs_norm1,final_sum)
        error=SQRT(rhs_norm1)
        if(debug_bicg.and.myid==1)write(*,'(A15,2es10.2,I5)')'|r_i|,iter',error,error/error_ini,iter
        if(verbose)write(*,112)iter,error,error/error_ini,error_cg_loc,r2,pap
        
!!$        if(error/error_ini .le. epsilon_diff)continue_loop=.false.
        
        max_loc=0.0d0
        do i=1,nb_ind
           this = liste_ind(i)
           max_loc=MAX(max_loc,ABS((alpha_cg*unew(this,3)+omega_cg(2)*unew(this,8))/uold(this,ind_cr)))
        enddo
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(max_loc,max_loc_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
        max_loc=max_loc_all
#endif
        if(debug_bicg.and.myid==1)write(*,'(A15,es10.2,I5)')'max_loc',max_loc
        if(   (error/error_ini .le. epsilon_diff).and.&
             &(max_loc         .le. epsilon_diff))continue_loop=.false.
        

        call dot_product(unew(:,1),unew(:,nvar+2),rhs_norm1,final_sum)
        if(abs(rhs_norm1).lt.epsilon_restartbicg)then
           !=======================================================
           ! Restart rzero if ``accuracy'' below some threshold
           ! unew(i,1     ): r_i
           ! unew(i,nvar+2): rzero
           ! unew(i,2     ): p_i
           !=======================================================           
           unew(:,nvar+2)=unew(:,1)
           unew(:,     2)=unew(:,1)
           call make_virtual_fine_dp(unew(1,2),ilevel)
           if(debug_bicg.and.myid==1)write(*,*)'Restarting BiCGSTAB'
        endif

        endif

        if(iter.ge.itermax)continue_loop=.false.

     end do
     ! End main iteration loop

     if(iter >= itermax)then
        if(myid==1)write(*,*)'Diffusion fail to converge...'
     end if

     if(myid==1) write(*,'(A,i4,A,4e11.3)')' BiCGSTAB :',iter, ' error=',error/error_ini,error_ini,error_cg_loc,error
     niter=niter+iter     
  end if
  end if
  end do
  !End loop over subcycles
end do

!=============================
! Update energy value
!=============================
if(.not.explicit_diff)then
   call cmp_energy_cr(2)
   ! Update boundaries
   call make_virtual_fine_dp(uold(1,5),ilevel)
   do irad=1,ncr
      call make_virtual_fine_dp(uold(1,ind_cr),ilevel)
   end do
endif
!  call upload_fine(ilevel)
!  call upload_fine(ilevel-1)

! Reset explicit_diff to false to allow for the next step to decide
if(autoswitch_diff)explicit_diff=.false.


111 format('   Entering conduction_cg')
112 format('   ==> Step=',i5,' Error=',5(1pe10.3,1x))!,e23.15)

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
 integer::this

 pdtvect=0.0d0
 local_sum = cmplx(0.0,0.0)

 do i=1,nb_ind
    this = liste_ind(i)
!!$    write(*,'(A,3e10.2)')'pdtvect=',pdtvect,fact1(this),fact2(this)
    pdtvect = pdtvect + fact1(this)*fact2(this)
    
 end do

 ! Compute global norms
#ifndef WITHOUTMPI
 call MPI_ALLREDUCE(pdtvect,pdtvect_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
 pdtvect   = pdtvect_all
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

!################################################################
function a_rkl2(j)
  implicit none
  integer::j
  real(dp)::a_rkl2
  
  a_rkl2=1d0-b_rkl2(j)

end function a_rkl2
!################################################################
function b_rkl2(j)
  implicit none
  integer::j
  real(dp)::b_rkl2
  real(dp)::x
  
  if(j==0.or.j==1.or.j==2)then
     b_rkl2=1d0/3d0
  else
     x=dble(j)
     b_rkl2=(x*x+x-2d0)/(2d0*x*(x+1d0))
  endif

end function b_rkl2
!################################################################
function mu_rkl2(j)
  implicit none
  integer::j
  real(dp)::mu_rkl2
  real(dp)::x
  
  x=dble(j)
  mu_rkl2=(2d0*x-1d0)/x * b_rkl2(j)/b_rkl2(j-1)

end function mu_rkl2
!################################################################
function nu_rkl2(j)
  implicit none
  integer::j
  real(dp)::nu_rkl2
  real(dp)::x
  
  x=dble(j)
  nu_rkl2=(1d0-x)/x * b_rkl2(j)/b_rkl2(j-2)

end function nu_rkl2
!################################################################
function mut_rkl2(j,s)
  implicit none
  integer::j,s
  real(dp)::mut_rkl2
  real(dp)::x,y
  
  y=dble(s)
  if(j==1)then
     mut_rkl2=4d0/(3d0*(y*y+y-2d0))
  else
     x=dble(j)
     mut_rkl2=4d0*(2d0*x-1d0)/(x*(y*y+y-2d0)) * b_rkl2(j)/b_rkl2(j-1)
  endif

end function mut_rkl2
!################################################################
function gt_rkl2(j,s)
  implicit none
  integer::j,s
  real(dp)::gt_rkl2
  real(dp)::x

  x=dble(j)
  gt_rkl2=-a_rkl2(j-1)*mut_rkl2(j,s)

end function gt_rkl2
!################################################################
end subroutine crdiffusion_cg

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_matrixA_crdiffusion (ilevel,compute,igroup)

  !------------------------------------------------------------------
  ! This routine computes the matrix A to vect_in and create vect_out
  ! compute = 1 : residu             	return B - Ax
  ! compute = 2 : Produit                 return  A.p
  ! compute = 3 : Preconditionner         return diag(A)
  ! compute = 4 : Compute flux in rad_flux
  !------------------------------------------------------------------
use amr_commons
use hydro_commons
implicit none
integer,intent(in)::compute,ilevel,igroup
integer ,dimension(1:nvector),save :: ind_grid
integer :: i,igrid,ngrid,ncache

if(numbtot(1,ilevel)==0)return
if(verbose)write(*,111)ilevel,compute

! Loop over active grids by vector sweeps
ncache=active(ilevel)%ngrid
do igrid=1,ncache,nvector
  ngrid=MIN(nvector,ncache-igrid+1)
  do i=1,ngrid
     ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
  end do
  call crdifffine1(ind_grid,ngrid,ilevel,compute,igroup)
end do

111 format('   Entering cmp_matrixA for level ',i2,i2)

end subroutine cmp_matrixA_crdiffusion
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine make_boundary_crdiffusion(ilevel,ind_cr)
use amr_commons
use hydro_commons
use constants,ONLY:kB,mH
!  use radiation_parameters
implicit none
! -------------------------------------------------------------------
! This routine set up boundary conditions for fine levels.
! -------------------------------------------------------------------
integer,intent(in)::ilevel,ind_cr
integer::ibound,boundary_dir,idim,inbor
integer::i,ncache,ivar,igrid,ngrid,ind,iperp1,iperp2
integer::iskip,iskip_ref,nx_loc,ix,iy,iz,gdim
integer,dimension(1:8)::ind_ref,alt
integer,dimension(1:2,1:4)::ind0
integer,dimension(1:nvector),save::ind_grid,ind_grid_ref
integer,dimension(1:nvector),save::ind_cell,ind_cell_ref

real(dp)::dx,dx_loc,scale
real(dp),dimension(1:3)::skip_loc,gs
real(dp),dimension(1:twotondim,1:3)::xc
real(dp),dimension(1:nvector,1:ndim),save::xx
real(dp),dimension(1:nvector,1:nvar+3),save::uu

if(.not. simple_boundary)return

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
              unew(ind_cell(i),3)     = unew(ind_cell_ref(i),3) !YD: for y(=M^-1*p)
              unew(ind_cell(i),8)     = unew(ind_cell_ref(i),8) !YD: for z(=M^-1*s)
              divu(ind_cell(i))=divu(ind_cell_ref(i))
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
                 divu(ind_cell(i))=unew(ind_cell(i),ind_cr)
              end if
           end do
        end if

     end do
     ! End loop over cells

  end do
  ! End loop over grids

end do
! End loop over boundaries

end subroutine make_boundary_crdiffusion
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine cmp_energy_cr(Etype)
use amr_commons
use hydro_commons
!  use radiation_parameters
implicit none
integer,intent(in) :: Etype ! Etype=1 : beginning ; Etype=2 : end
integer ::i,this,igroup,ind_cr
real(dp)::rho,ecr,TCR,coef
real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

! Conversion factor from user units to cgs units
call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

ind_cr=8
if(twotemp)ind_cr=9

if(fix_temp_diff.and.Etype==2)then
   do i=1,nb_ind
      this = liste_ind(i)
      do igroup=1,ncr
         unew(this,ind_cr+igroup)=MAX(unew(this,ind_cr+igroup),smallcr)
      end do
   enddo
endif

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
!!$        print*,unew(this,ind_cr+igroup),'1'
     end do

  elseif(Etype==2)then
     ! update total energy
     uold(this,5) = unew(this,5)
     do igroup=1,ncr
        uold(this,ind_cr+igroup) = unew(this,ind_cr+igroup)
        uold(this,5) = uold(this,5) + unew(this,ind_cr+igroup)
!!$        print*,unew(this,ind_cr+igroup),'2'
     end do
     if(TCRmax.gt.0.0d0)then
        do igroup=1,ncr
           uold(this,5) = uold(this,5) - uold(this,ind_cr+igroup)
           coef=(gamma_rad(igroup)-1d0)/rho*scale_T2
           TCR=coef*uold(this,ind_cr+igroup)
           uold(this,ind_cr+igroup) = MIN(TCR,TCRmax)/coef
           uold(this,5) = uold(this,5) + uold(this,ind_cr+igroup)
        end do
     endif
     if(TCRmin.gt.0.0d0)then
        do igroup=1,ncr
           uold(this,5) = uold(this,5) - uold(this,ind_cr+igroup)
           coef=(gamma_rad(igroup)-1d0)/rho*scale_T2
           TCR=coef*uold(this,ind_cr+igroup)
           uold(this,ind_cr+igroup) = MAX(TCR,TCRmin)/coef
           uold(this,5) = uold(this,5) + uold(this,ind_cr+igroup)
           if (TCRmin.gt.TCR)print*,'TCR gt TCRmin',TCRmin/coef
        end do
     endif

  end if
end do



end subroutine cmp_energy_cr

!################################################################
!################################################################
!################################################################ 
!################################################################
