subroutine conduction_cg (ilevel,itemp)
  use amr_commons
  use hydro_commons
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
  !  x1(i)    : stored in uold(i,8+nrad+ntp)(i.e. new thermal energy at time n+1)
  !  b1(n)    : stored in unew(i,8+nrad+ntp)(i.e. thermal energy at time n)
  ! x1(i-1)   : stored in unew(i,7)
  !  Te^n     : stored in unew(i,nvar+3)
  !  Te(i)    : stored in unew(i,5)
  !  T^n      : stored in unew(i,nvar+2)
  !===========================================================================
  integer,intent(IN)::ilevel,itemp
  complex*16 :: final_sum
  real(dp)::error,error_ini,epsilon
  real(dp)::error_cg_loc,error_cg_all
  real(dp)::alpha_cg,beta_cg,Cv
  real(dp)::r2_old,r2,pAp,rhs_norm1
  real(dp)::kpara_ana,temp
  integer :: i,info,ind,iter,iskip,itermax
  integer :: this,nsub_imp
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::scale_kappa
  real(dp)::min_ener,min_ener_all,max_ener,max_ener_all
  real(dp)::dener,ratio
  logical::exist_leaf_cell=.true.

  real(dp)::Tau_ei_vol,Te_new,Te_old,Cvi,Cve
  real(dp)::scale_tau,Ti_new,Ti_old,T_new
  real(dp)::n_electron,n_ion,cfl,cfl_loc,dt_cond,temp_loc
  real(dp)::tt1,tt2,tt2b,deltatt_cond,tt1iter,tt2end,deltatt_end
  integer ::ind_eint,ind_Tenew,ind_Teold,ind_Cve
  if(verbose)write(*,111)
  if(numbtot(1,ilevel)==0)return

#ifndef WITHOUTMPI
  dtc12=0.0d0;dtc23=0.0d0;dtc34=0.0d0;dtcfath=0.0d0;dtcneig=0.0d0
  deltatt_cond =0.0d0
  tt1=MPI_WTIME()
#endif

  ind_eint = nvar+2
  ind_Tenew = 5 
  ind_Cve = nvar+1
  ind_Teold = nvar+3

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kappa=scale_d*scale_l*scale_v**3
  scale_tau=1d0/(scale_t*scale_nH) ! Time in code units
  call compute_Tau_equi_vol(Tau_ei_vol)

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
     !write(*,*)'No leaf-cell - myid=',myid
     exist_leaf_cell=.false.
  end if

  do i=1,nb_ind
     this = liste_ind(i)
     unew(this,1:nvar)=0.0d0
     divu(this)=0.0d0
     enew(this)=0.0d0
  end do

  ! Set constants
  epsilon = epsilon_diff

  !===================================================================
  ! Compute gas temperature stored in unew(i,ind_Teold) and in unew(i,ind_Tenew)
  ! 1: Etot -> Tp
  ! 2: Tp   -> Etot
  !===================================================================
  call cmp_energy_conduction(1,itemp)

  !===================================================================
  ! Begin of subcycles....
  !===================================================================
  dt_imp = dtnew(ilevel)

  dener=0.0d0
  max_ener=0.0d0
  min_ener=1.0d30
  error_cg_loc=0.0d0
  do i=1,nb_ind
     this = liste_ind(i)
     max_ener=max(max_ener, unew(liste_ind(i),ind_Tenew))
     min_ener=min(min_ener, unew(liste_ind(i),ind_Tenew))
     error_cg_loc = error_cg_loc + unew(liste_ind(i),ind_Tenew)
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
  dt_imp=dt_imp/real(nsub_imp)

!!$  if(dener .gt. 1.d4)then
!!$     nsub_imp=MIN(int(dener**(2./5.)),10)
!!$     dt_imp=dt_imp/real(nsub_imp)
!!$  endif

  if(myid==1)then
     write(*,*)'=============================================='
     if(itemp==1)write(*,'(A,I4,A,e10.3,A,e10.3,A,e10.3)')'(Conduc E) ilevel',ilevel,' MAXIMUM OF DENER=',dener,' err_cg_loc=',error_cg_loc,' dt=',dt_imp
     if(itemp==2)write(*,'(A,I4,A,e10.3,A,e10.3,A,e10.3)')'(Conduc I) ilevel',ilevel,' MAXIMUM OF DENER=',dener,' err_cg_loc=',error_cg_loc,' dt=',dt_imp
  endif


  !===================================================
  ! Compute thermal coefficient :
  ! Spitzer conductivity : stored in divu(indcell(i))
  ! Room left for Spitzer relaxation term between ion and electron in enew(indcell(i)) 
  !===================================================
  cfl_loc=0.0d0
  temp_loc=0.0d0
  do i=1,nb_ind
     this = liste_ind(i)
     temp = unew(this,ind_Tenew)
     divu(this)= kpara_ana(temp,itemp)/scale_kappa
     if(.not.testcase)then
        dt_cond=kpara_ana(temp,itemp)/(uold(this,1)*scale_d*1.38d-16/(mu_electron*1.66d-24)/(gamma-1.0d0))/(scale_l**2/scale_t)
        if(itemp==1)then
           dt_cond=dt_cond/mu_electron
        else
           dt_cond=dt_cond/mu_ion
        endif
     else
        dt_cond=kpara_ana(temp,itemp)/uold(this,1)
     endif
     cfl_loc =max(cfl_loc,dt_cond)
     temp_loc=max(temp,temp_loc)
  end do
  ! Compute maximum cfl
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(cfl_loc,cfl,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  cfl_loc=cfl
  call MPI_ALLREDUCE(temp_loc,temp,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  temp_loc=temp
#endif

  if(myid==1)write(*,*)'CFL=',0.5*(boxlen*0.5**ilevel)**2/(2.0d0*cfl_loc)/dt_imp,temp_loc

  ! Update virtual boundaries
  call make_virtual_fine_dp(unew(1,2),ilevel)
  call make_virtual_fine_dp(unew(1,3),ilevel)
  call make_virtual_fine_dp(unew(1,1),ilevel)
  call make_virtual_fine_dp(unew(1,4),ilevel)
  call make_virtual_fine_dp(divu(1),ilevel)
  call make_virtual_fine_dp(unew(1,ind_Tenew),ilevel)
  call make_virtual_fine_dp(unew(1,ind_Teold),ilevel)
  call make_virtual_fine_dp(unew(1,ind_Cve),ilevel)
  call make_virtual_fine_dp(unew(1,ind_eint),ilevel)
  ! Update physical boundaries
  call make_boundary_conduction(ilevel,itemp)

  if(semi_implicit)then 
     ! Do the first explicit step for the component
     ! tangential to the face
     call cmp_matrixA_cond (ilevel,0,itemp)
     do i=1,nb_ind
        Cv = unew(liste_ind(i),ind_Cve)
        ! Only works with 1T
        unew(liste_ind(i),ind_Tenew)=unew(liste_ind(i),ind_eint)/Cv
        if(fix_temp_diff)then
           unew(liste_ind(i),ind_Tenew)=MAX(unew(liste_ind(i),ind_Tenew),Tfloor)
        endif
        if(twotemp)then
           if(itemp==1)then
              Cve = Cv
              Cvi = Cv*(mu_electron/mu_gas-1.0d0)
              Te_old = unew(liste_ind(i),ind_Teold)
              Te_new = unew(liste_ind(i),ind_Tenew)
              Ti_old = (unew(liste_ind(i),ind_eint)-uold(liste_ind(i),9))/Cvi
              n_ion = uold(liste_ind(i),1)/mu_ion
              n_electron = uold(liste_ind(i),1)/mu_electron                 
              Ti_new = Ti_old
              unew(liste_ind(i),ind_Teold)=unew(liste_ind(i),ind_Tenew)
              uold(liste_ind(i),9)=unew(liste_ind(i),ind_Tenew)*Cv
              unew(liste_ind(i),ind_eint)=Cvi*Ti_new+Cve*Te_new
           else
              Cvi = Cv
              Cve = Cv*(mu_ion/mu_gas-1.0d0)
              Ti_old = unew(liste_ind(i),ind_Teold)
              Ti_new = unew(liste_ind(i),ind_Tenew)
              Te_old = uold(liste_ind(i),9)/Cve
              unew(liste_ind(i),ind_Teold)=unew(liste_ind(i),ind_Tenew)
              unew(liste_ind(i),ind_eint)=Cvi*Ti_new+Cve*Te_old
           endif
        else
           T_new = unew(liste_ind(i),ind_Tenew)
           unew(liste_ind(i),ind_eint)=Cv*T_new
        end if
     end do
     call cmp_energy_conduction(2,itemp)

     ! Update virtual boundaries
     call make_virtual_fine_dp(unew(1,2),ilevel)
     call make_virtual_fine_dp(unew(1,3),ilevel)
     call make_virtual_fine_dp(unew(1,1),ilevel)
     call make_virtual_fine_dp(unew(1,4),ilevel)
     call make_virtual_fine_dp(divu(1),ilevel)
     call make_virtual_fine_dp(unew(1,ind_Tenew),ilevel)
     call make_virtual_fine_dp(unew(1,ind_Teold),ilevel)
     call make_virtual_fine_dp(unew(1,ind_Cve),ilevel)
     call make_virtual_fine_dp(unew(1,ind_eint),ilevel)
     call make_virtual_fine_dp(uold(1,5),ilevel)
     if(twotemp)call make_virtual_fine_dp(uold(1,9),ilevel)
     ! Update physical boundaries
     call make_boundary_conduction(ilevel,itemp)

     call cmp_energy_conduction(1,itemp)
  endif

  !===================================================
  ! Compute r1 = b1 - A1x1 and store it into unew(i,1)
  ! Also set p1 = r1 and store it into unew(i,2)
  !===================================================
  call cmp_matrixA_cond (ilevel,1,itemp)

  call make_virtual_fine_dp(unew(1,1),ilevel)
  call make_virtual_fine_dp(unew(1,2),ilevel)

  !===================================
  ! Compute right-hand side norm (max)
  !===================================
  call dot_product(unew(:,1),unew(:,1),rhs_norm1,final_sum)

  !===================================================
  ! Compute Preconditionner M=1/diag(A) and store it in unew(i,4)
  !===================================================
  !ben
  call make_boundary_conduction(ilevel,itemp)
  call cmp_matrixA_cond (ilevel,3,itemp)

  call make_virtual_fine_dp(unew(1,4),ilevel)
  !====================================
  ! MAIN ITERATION LOOP
  !====================================     

  iter=0; itermax=5000

  error_ini=sqrt(rhs_norm1)

  !     error_ini=sqrt(real(final_sum))
  error=error_ini
  error_cg_loc=1.0d0

  !     write(*,*)error,error_ini
  if(error_ini.ne.0.0d0)then
     ratio=error/error_ini
  else
     ratio=epsilon
  endif
  !ratio=error/error_ini ! bypass previous conditional statement

  !do while(error/error_ini>epsilon .and.iter<itermax)
#ifndef WITHOUTMPI
  tt1iter=MPI_WTIME()
#endif
  do while((ratio>epsilon .or. error_cg_loc>epsilon) .and.iter<itermax)

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
     r2=r2
     !r2=real(final_sum)

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
     call make_boundary_conduction(ilevel,itemp)
     call make_virtual_fine_dp(unew(1,2),ilevel)

     !==============================================
     ! Compute q1 = Ap1 and store it into unew(i,3)
     !==============================================
!!$#ifndef WITHOUTMPI
!!$        call MPI_BARRIER(MPI_COMM_WORLD,info)
!!$        tt1b=MPI_WTIME()
!!$#endif
     call cmp_matrixA_cond (ilevel,2,itemp)
!!$#ifndef WITHOUTMPI
!!$        call MPI_BARRIER(MPI_COMM_WORLD,info)
!!$        tt2b=MPI_WTIME()
!!$        deltatt_cond=deltatt_cond+(tt2b-tt1b)
!!$#endif

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
        error_cg_loc=max(error_cg_loc, abs((alpha_cg*unew(liste_ind(i),2))/unew(liste_ind(i),ind_Teold)))
     end do

     ! Compute maximum error
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(error_cg_loc,error_cg_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     error_cg_loc=error_cg_all
#endif

     call cX_plus_Y_to_Z (alpha_cg,unew(:,2),unew(:,ind_Tenew),unew(:,ind_Tenew))

     !====================================
     ! Recurrence on r (unew(i,1))
     !====================================

     call cX_plus_Y_to_Z (- alpha_cg ,unew(:,3),unew(:,1),unew(:,1))

     !===================================
     ! Compute right-hand side norm (max)
     !===================================
     call dot_product(unew(:,1),unew(:,1),rhs_norm1,final_sum)

     error=SQRT(rhs_norm1)
     ratio=error/error_ini

     !        error=SQRT(real(final_sum))
     !        error = error_cg_loc
     if(verbose) write(*,112)iter,error,error/error_ini,error_cg_loc,r2,pap
!!$        if(myid==1) write(*,112)iter,error,error/error_ini,error_cg_loc,r2,pap

#ifndef WITHOUTMPI
     tt2end=MPI_WTIME()
     deltatt_end=deltatt_end+(tt2end-tt2b)
#endif
  end do
  ! End main iteration loop

  if(iter >= itermax)then
     if(myid==1)write(*,*)'Conduction fail to converge...'
  end if

  if(myid==1) write(*,'(A,i5,A,4e10.3)')' CG :',iter, ' error=',ratio,error_ini,error_cg_loc,error
  niter=niter+iter

  !====================================
  ! Update gas temperature
  !====================================
  do i=1,nb_ind
     Cv = unew(liste_ind(i),ind_Cve)
     if(fix_temp_diff)then
        unew(liste_ind(i),ind_Tenew)=MAX(unew(liste_ind(i),ind_Tenew),Tfloor)
     endif
     if(twotemp)then
        if(itemp==1)then
           Cve = Cv
           Cvi = Cv*(mu_electron/mu_gas-1.0d0)
           Te_old = unew(liste_ind(i),ind_Teold)
           Te_new = unew(liste_ind(i),ind_Tenew)
           Ti_old = (unew(liste_ind(i),ind_eint)-uold(liste_ind(i),9))/Cvi
           Ti_new = Ti_old

           !update electron temperature and electron internal energy
           unew(liste_ind(i),ind_Teold)=unew(liste_ind(i),ind_Tenew)
           uold(liste_ind(i),9)=unew(liste_ind(i),ind_Tenew)*Cv
           !update internal energy
           unew(liste_ind(i),ind_eint)=Cvi*Ti_new+Cve*Te_new

        else

           Cvi = Cv
           Cve = Cv*(mu_ion/mu_gas-1.0d0)
           Ti_old = unew(liste_ind(i),ind_Teold)
           Ti_new = unew(liste_ind(i),ind_Tenew)
           Te_old = uold(liste_ind(i),9)/Cve
           !update electron temperature and electron internal energy
           unew(liste_ind(i),ind_Teold)=unew(liste_ind(i),ind_Tenew)
           !update internal energy
           unew(liste_ind(i),ind_eint)=Cvi*Ti_new+Cve*Te_old

        endif


     else
        T_new = unew(liste_ind(i),ind_Tenew)
        unew(liste_ind(i),ind_eint)=Cv*T_new
     end if

  end do

  !=============================
  ! Update energy value
  !=============================
  call cmp_energy_conduction(2,itemp)

#ifndef WITHOUTMPI
  tt2=MPI_WTIME()
  if(myid==1)write(*,*)'Time elapsed in conduction:',tt2-tt1
!!$  if(myid==1)write(*,*)'Time elapsed in conduction (iteration):',tt2-tt1iter
!!$  if(myid==1)write(*,*)'Time elapsed in cmpMatrix:',deltatt_cond
  call MPI_BARRIER(MPI_COMM_WORLD,info)
!!$  if(nb_ind>0)write(*,'(I3,A,5e12.5)')myid,',Time=',dtc12,dtc23,dtc34,dtcneig,dtcfath
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

  ! Update boundaries
  call make_virtual_fine_dp(uold(1,5),ilevel)
  if(twotemp)call make_virtual_fine_dp(uold(1,9),ilevel)

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
    complex*16 ::global_sum
    integer::this

    pdtvect=0.0d0
    local_sum = cmplx(0.0,0.0)
    global_sum = cmplx(0.0,0.0)

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

end subroutine conduction_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_matrixA_cond (ilevel,compute,itemp)
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
  integer,intent(IN)::compute,ilevel,itemp
  integer ,dimension(1:nvector) :: ind_grid
  integer :: i,igrid,ngrid,ncache

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call condfine1(ind_grid,ngrid,ilevel,compute,itemp)
  end do
111 format('   Entering cmp_matrixA_cond')

end subroutine cmp_matrixA_cond
!################################################################
!################################################################
!################################################################ 
!################################################################
function kpara_ana(Tp,itemp)
  use amr_commons,ONLY:dp
  use hydro_commons,only:fudge_kspitzer,testcase!,mu_electron,Tfloor
  implicit none
  real(dp)::Tp,kpara_ana
  integer::itemp

  if(.not.testcase)then
     if(itemp==1)then
        kpara_ana=fudge_kspitzer*MAX(Tp,0.0d0)**2.5
     else
        kpara_ana=fudge_kspitzer*MAX(Tp,0.0d0)**2.5 *0.023337031 ! (me/mp)^0.5=0.023337031
     endif
  else
     kpara_ana=1.0d0
  endif
  !kpara_ana=0.0d0
  !kpara_ana=MAX(kpara_ana,1.0d0)
  !write(*,*)'kpar=',kpara_ana
  
  return 
end function kpara_ana
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine compute_Tau_equi_vol(Tau_ei_vol)
  use amr_commons,ONLY:dp
  use hydro_commons,only:coulomb_log
  use constants,ONLY:kB,mH
  implicit none
  real(dp)::Z_i,electron_charge,electron_mass
  real(dp),intent(OUT)::Tau_ei_vol

  ! Compute coupling time between ion and electron
  Z_i = 1d0 ! For H: mH/ZH^2=mp/1 For He: mHe/ZHe^2=4mp/(2ZH)^2=mp/1
            ! (does not work for heavier elements)
  electron_charge = 4.803204d-10
  electron_mass = 9.10938291d-28
  Tau_ei_vol = 3.0d0*electron_mass*mH/(8.0d0*sqrt(2.0d0*acos(-1.0d0))*(Z_i**2)*(electron_charge**4) &
       & * coulomb_log) *((kB/electron_mass)**(1.5))

  ! The value of tau_ei_vol is in seconds per K^3/2 and per L^3

  return 
end subroutine compute_Tau_equi_vol
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine make_boundary_conduction(ilevel,itemp)
  use amr_commons
  use hydro_commons
  use constants,ONLY:kB,mH
  implicit none
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer,intent(IN)::ilevel,itemp
  integer::ibound,boundary_dir,idim,inbor
  integer::i,ncache,ivar,igrid,ngrid,ind,iperp1,iperp2
  integer::iskip,iskip_ref,nx_loc,ix,iy,iz,gdim
  integer,dimension(1:8)::ind_ref,alt
  integer,dimension(1:2,1:4)::ind0
  integer,dimension(1:nvector)::ind_grid,ind_grid_ref
  integer,dimension(1:nvector)::ind_cell,ind_cell_ref

  real(dp)::dx,dx_loc,scale,switch
  real(dp),dimension(1:3)::skip_loc,gs
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(dp),dimension(1:nvector,1:nvar+3)::uu
  real(dp)::dd,usquare,emag,erad_loc,eps,ekin,rho
  real(dp)::scale_kappa
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::kpara_ana,Cv1,Cv1_electron,Cv1_ion
  integer ::ind_eint,ind_Tenew,ind_Teold,ind_Cve
#if NENER>0
  integer::igroup
#endif

  ind_eint = nvar+2
  ind_Tenew = 5 
  ind_Cve = nvar+1
  ind_Teold = nvar+3

  if(.not. simple_boundary)return
  if(verbose)write(*,111)
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  scale_kappa=scale_d*scale_l*scale_v**3
  if(.not.testcase)then
     Cv1         =kB/(mu_gas*mH*(gamma-1.0d0))/scale_v**2
     if(twotemp)then
        Cv1_electron=kB/(mu_electron*mH*(gamma_rad(1)-1.0d0))/scale_v**2
        Cv1_ion     =kB/(mu_ion     *mH*(gamma       -1.0d0))/scale_v**2
     else
        Cv1_electron=kB/(mu_gas     *mH*(gamma       -1.0d0))/scale_v**2 ! Not used in that case
     endif
  else
     if(twotemp)then
        Cv1         =2.0d0
        Cv1_electron=1.0d0
     else
        Cv1         =1.0d0 ! For Heaviside test
        Cv1_electron=Cv1
     endif
  endif
!!$  write(*,*)'=============='
!!$  write(*,*)Cv1_electron
!!$  write(*,*)kB
!!$  write(*,*)mu_electron
!!$  write(*,*)mH
!!$  write(*,*)gamma_rad(1)
!!$  write(*,*)scale_v
!!$  write(*,*)'=============='

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
              ! Remove magnetic energy
              do i=1,ngrid
                 emag=0.125d0*((uu(i,6)+uu(i,nvar+1))**2+ &
                      &        (uu(i,7)+uu(i,nvar+2))**2+ &
                      &        (uu(i,8)+uu(i,nvar+3))**2)
                 uu(i,5)=uu(i,5)-emag
              end do

              ! Scatter to boundary region
              do ivar=1,nvar+3
                 switch=1
                 if(ivar>1.and.ivar<=4)switch=gs(ivar-1)
                 if(ivar.ne.(5+gdim).and.ivar.ne.(nvar+gdim))then
                    do i=1,ngrid
                       uold(ind_cell(i),ivar)=uu(i,ivar)*switch
                    end do
                 endif
                 if(ivar==(5+gdim))then
                    do i=1,ngrid
                       uold(ind_cell(i),5+gdim)=uu(i,5+gdim)+(uu(i,nvar+gdim)-uu(i,5+gdim))*alt(ind)
                    end do
                 endif
                 if(ivar==(nvar+gdim))then
                    do i=1,ngrid
                       uold(ind_cell(i),nvar+gdim)=uu(i,nvar+gdim)+(uu(i,nvar+gdim)-uu(i,5+gdim))*alt(ind)
                    end do
                 endif
              end do
              ! Add magnetic energy
              do i=1,ngrid
                 emag=0.125d0*((uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))**2+ &
                      &        (uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))**2+ &
                      &        (uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))**2)
                 uold(ind_cell(i),5)=uold(ind_cell(i),5)+emag
                 unew(ind_cell(i),2)=unew(ind_cell_ref(i),2)

                 unew(ind_cell(i),ind_Tenew)     = unew(ind_cell_ref(i),ind_Tenew)
                 unew(ind_cell(i),ind_Teold)     = unew(ind_cell_ref(i),ind_Teold)
                 unew(ind_cell(i),ind_eint )     = unew(ind_cell_ref(i),ind_eint )
                 unew(ind_cell(i),ind_Cve  )     = unew(ind_cell_ref(i),ind_Cve  )
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
                 dd=max(uu(i,1),smallr)

                 usquare=0.0d0
                 do idim=1,ndim
                    usquare=usquare+(uu(i,idim+1)/uu(i,1))**2
                 end do
                 ! Compute total magnetic energy
                 emag = 0.0d0
                 do ivar=1,3
                    emag = emag + 0.125d0*(uu(i,5+ivar) &
                         &  +uu(i,nvar+ivar))**2
                 end do

                 ! Compute total non-thermal energy
                 erad_loc = 0.0d0
#if NENER>0
                 if(twotemp)then
#if NENER>2
                    do igroup=2,nener
                       erad_loc = erad_loc + uu(i,8+igroup)
                    enddo
#endif
                 else
                    do igroup=1,nener
                       erad_loc = erad_loc + uu(i,8+igroup)
                    enddo
                 end if
#endif

                 rho   = uu(i,1)
                 ekin  = rho*usquare/2.0d0
                 eps   = (uu(i,5)-ekin-emag-erad_loc)

                 if(.not.twotemp)then                 
                    unew(ind_cell(i),nvar+3) = eps/(Cv1_electron*rho)
                 else
                    if(itemp==1)then
                       unew(ind_cell(i),nvar+3) = uu(i,nvar-ncr)/(Cv1_electron*rho)
                    else
                       unew(ind_cell(i),nvar+3) = (eps-uu(i,nvar-ncr))/(Cv1_ion*rho)
                    endif
                 endif

                 divu(ind_cell(i))=kpara_ana(unew(ind_cell(i),nvar+3),itemp)/scale_kappa

                 unew(ind_cell(i),2)=  0.0d0
                 do ivar=1,nvar
                    uold(ind_cell(i),ivar)=uu(i,ivar)
                 enddo

              end if
           end do
           end if

        end do
        ! End loop over cells

     end do
     ! End loop over grids

  end do
  ! End loop over boundaries

if(verbose)write(*,112)
111 format('   Entering make_boundary_conduction for level ',I2)
112 format('   Leaving make_boundary_conduction for level ',I2)

end subroutine make_boundary_conduction
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine cmp_energy_conduction(Etype,itemp)
  use amr_commons
  use hydro_commons
  use constants,ONLY:kB,mH
  implicit none
  integer,intent(in) :: Etype ! Etype=1 : beginning ; Etype=2 : end
  integer,intent(in) :: itemp
  integer ::i,idim,this,mvar
  real(dp)::usquare,Cv,eps,ekin,emag,rho,erad_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::Tp_loc
  real(dp)::eps_electron,eps_ion
  real(dp)::Cv_electron,Cv_ion,Cv_loc
  real(dp)::Cv1_electron,Cv1_ion,Cv1
#if NENER>0
  integer::igroup
#endif
  ! EOS
  integer ::ind_eint,ind_Tenew,ind_Teold,ind_Cve

  ind_eint = nvar+2
  ind_Tenew = 5 
  ind_Cve = nvar+1
  ind_Teold = nvar+3

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(.not.testcase)then
     Cv1         =kB/(mu_gas*mH*(gamma-1.0d0))/scale_v**2
     Cv1_ion     =kB/(mu_ion*mH*(gamma-1.0d0))/scale_v**2
     if(twotemp)then
        Cv1_electron=kB/(mu_electron*mH*(gamma_rad(1)-1.0d0))/scale_v**2
     else
        Cv1_electron=kB/(mu_electron*mH*(gamma       -1.0d0))/scale_v**2
     endif
  else
     if(twotemp)then
        Cv1         =2.0d0
        Cv1_ion     =2.0d0
        Cv1_electron=1.0d0
     else
        Cv1         =1.0d0
        Cv1_ion     =1.0d0
        Cv1_electron=Cv1
     endif
  endif
  
  do i=1,nb_ind
     this = liste_ind(i)
     rho   = uold(this,1)

     Cv = rho*Cv1 
     if(twotemp)then
        Cv_electron = rho*Cv1_electron
        Cv_ion      = rho*Cv1_ion
     endif

     ! Compute total kinetic energy
     usquare=0.0d0
     do idim=1,ndim
        usquare=usquare+(uold(this,idim+1)/uold(this,1))**2
     end do
     ekin  = rho*usquare*0.5d0

     ! Compute total magnetic energy
     emag = 0.0d0
     do mvar=1,3
        emag = emag + 0.125d0*(uold(this,5+mvar)+uold(this,nvar+mvar))**2
     end do

     ! Compute total non-thermal energy
     erad_loc = 0.0d0
#if NENER>0
     if(twotemp)then
#if NENER>1
        do igroup=2,nener
           erad_loc = erad_loc + uold(this,8+igroup)
        enddo
#endif
     else
        do igroup=1,nener
           erad_loc = erad_loc + uold(this,8+igroup)
        enddo
     end if
#endif
     
     if(Etype==1)then
        if(.not.twotemp)then
           eps    = uold(this,5)-ekin-emag-erad_loc
           Tp_loc = eps/Cv
           Cv_loc = Cv
        else
           eps_electron = uold(this,9)
           eps_ion      = uold(this,5)-ekin-emag-erad_loc-eps_electron
           eps          = uold(this,5)-ekin-emag-erad_loc
           Tp_loc = eps_electron/Cv_electron
           Cv_loc = Cv_electron
        end if

        if(itemp==1)then
           unew(this,ind_Teold) = Tp_loc ! Store electron temperature
           unew(this,ind_Tenew) = unew(this,ind_Teold)
           unew(this,ind_Cve  ) = Cv_loc ! Store Cve = rho*kB/(mu_electron*mH*(gamma-1)) in code units
        else
           unew(this,ind_Teold) = eps_ion/Cv_ion ! Store ion temperature
           unew(this,ind_Tenew) = unew(this,ind_Teold)
           unew(this,ind_Cve  ) = Cv_ion ! Store Cvi = rho*kB/(mu_ion*mH*(gamma-1)) in code units
        endif

        ! Store gas total internal energy
        unew(this,ind_eint) = eps

     else if(Etype==2)then
        eps = unew(this,ind_eint)
        uold(this,5) = eps + ekin + emag + erad_loc
     end if
  end do

end subroutine cmp_energy_conduction
!################################################################
!################################################################
!################################################################ 
!################################################################
