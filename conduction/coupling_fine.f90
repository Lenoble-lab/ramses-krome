subroutine coupling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute coupling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for coupling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call couplfine1(ind_grid,ngrid,ilevel)
  end do

  ! Update boundaries
  call make_virtual_fine_dp(uold(1,5),ilevel)
  call make_virtual_fine_dp(uold(1,nvar),ilevel)

111 format('   Entering coupling_fine for level',i2)

end subroutine coupling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine couplfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8)::Tau_ei_vol,Cv1_electron,Cv1_ion,Cv_electron,Cv_ion
  real(kind=8)::Te,Ti,ekkscalar,emagscalar,errscalar,n_ion,n_electron
  real(kind=8)::Tau_ei_reduced,Tau_ie_reduced,Tau_ei,Tau_ie,dt_loc,halfdt,onehalfdt
  real(kind=8)::oneoverte05,oneoverte15,aa,bb1,bb2
  real(kind=8)::f1,f2,df1d1,df2d1,df1d2,df2d2,oneoverden
  real(kind=8)::dTi,dTe,Ti2,Te2,error,deltaTe,deltaTi
  real(kind=8)::dtsum,dTirel,dTerel,dTrelmax,diffT,dt_loc2,Te0,Ti0
  integer::ninc
  integer::niterimp
  logical::oknonzero
#if NENER>0
  integer::irad
#endif
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do
     if(nleaf.eq.0)cycle

     call  compute_Tau_equi_vol(Tau_ei_vol)
     !scale_tau=1d0/scale_t ! Time in code units
     Cv1_electron=kB/(mu_electron*mH*(gamma_rad(1)-1.0d0))/scale_v**2
     Cv1_ion     =kB/(mu_ion     *mH*(gamma       -1.0d0))/scale_v**2
     do i=1,nleaf
        Cv_electron=uold(ind_leaf(i),1)*Cv1_electron
        Cv_ion     =uold(ind_leaf(i),1)*Cv1_ion
        Te=uold(ind_leaf(i),9)/Cv_electron
        ! Compute ion temperature
        Ti=uold(ind_leaf(i),5)
        ekkscalar=0.0d0
        do idim=1,3
           ekkscalar=ekkscalar+0.5*uold(ind_leaf(i),idim+1)**2/MAX(uold(ind_leaf(i),1),smallr)
        end do
        emagscalar=0.0d0
        do idim=1,3
           emagscalar=emagscalar+0.125d0*(uold(ind_leaf(i),idim+5)+uold(ind_leaf(i),idim+nvar))**2
        end do
        errscalar=0.0d0
#if NENER>0
        do irad=1,nener
           errscalar=errscalar+uold(ind_leaf(i),8+irad)
        end do
#endif
        Ti=(Ti-ekkscalar-emagscalar-errscalar)/Cv_ion
        Te0=Te;Ti0=Ti
        !write(*,'(A,I10,5es10.3)')'before',ind_leaf(i),Te*Cv_electron,Ti*Cv_ion,ekkscalar,emagscalar,errscalar
        n_ion=uold(ind_leaf(i),1)/mu_ion*scale_d/mH
        n_electron=n_ion*mu_ion/mu_electron
        
        dt_loc=dtnew(ilevel)
        dtsum=0.0d0
        ninc=0
        oknonzero=.true.
        do while(dtsum .lt. dtnew(ilevel))
           ninc=ninc+1
           if(oknonzero)then
              Tau_ei_reduced = Tau_ei_vol*MAX(Te,Tfloor)**(1.5)/n_ion/scale_t
              Tau_ie_reduced = Tau_ei_reduced*n_ion/n_electron
              Tau_ei = Tau_ei_vol/n_ion     /scale_t
              Tau_ie = Tau_ei_vol/n_electron/scale_t
              
              diffT=abs(Te-Ti)
              dt_loc2=dtnew(ilevel)-dtsum
              dTi=diffT*dt_loc2/Tau_ie_reduced
              dTe=diffT*dt_loc2/Tau_ei_reduced
              dTirel=dTi/Ti;dTerel=dTe/Te;dTrelmax=MAX(dTirel,dTerel)
              if(dTrelmax > RelVar)then
                 dt_loc=dt_loc2/(dTrelmax/RelVar)
              else
                 dt_loc=dt_loc2
              endif
              if(dtsum+dt_loc>dtnew(ilevel))dt_loc=dtnew(ilevel)-dtsum
           else
              ! if zero or negative values of temperature have been found
              ! try with half reduced time step
              dt_loc=dt_loc/2d0
           endif
              
           halfdt   =0.5d0*dt_loc
           onehalfdt=1.5d0*dt_loc

           Ti2=Ti
           Te2=Te
           niterimp=0
           error=1.0
           deltaTi=0.0d0;deltaTe=0.0d0
           oknonzero=.true.
           do while(error > 1d-4 .and. oknonzero)
!!$              if(ind_leaf(i)==2012316.and.myid==10)then
!!$                 write(*,*)'Te =',Te
!!$                 write(*,*)'Ti =',Ti
!!$                 write(*,*)'Te2=',Te2
!!$                 write(*,*)'Ti2=',Ti2
!!$                 write(*,*)'dTe=',dTe
!!$                 write(*,*)'dTi=',dTi
!!$                 write(*,*)'dTerel=',dTerel
!!$                 write(*,*)'dTirel=',dTirel
!!$                 write(*,*)'dTe (guess1)=',diffT*dt_loc2/Tau_ei_reduced
!!$                 write(*,*)'dTi (guess1)=',diffT*dt_loc2/Tau_ie_reduced
!!$                 write(*,*)'dt_loc =',dt_loc
!!$                 write(*,*)'dt_loc2=',dt_loc2
!!$                 write(*,*)'ind_leaf(i)=',ind_leaf(i),myid
!!$                 write(*,*)'niter=',niterimp
!!$              endif
              oneoverte05 = Te2**(-0.5)
              oneoverte15 = Te2**(-1.5)
              
              aa =dt_loc*(oneoverte05-Ti2*oneoverte15)
              bb1=dt_loc* oneoverte15
              bb2=halfdt*oneoverte15-onehalfdt*Ti2/Te2**2.5
              
              f1= aa-Tau_ie*(Ti2-Ti)
              f2=-aa-Tau_ei*(Te2-Te)
              
              df1d1=-bb1-tau_ie
              df2d1= bb1
              df1d2=-bb2
              df2d2= bb2-tau_ei
              oneoverden=1d0/(df1d1*df2d2-df2d1*df1d2)
              dTi=(-f1*df2d2+f2*df1d2)*oneoverden
              dTe=(-f2*df1d1+f1*df2d1)*oneoverden
              
              Ti2=Ti2+dTi
              Te2=Te2+dTe

              error=MAX(abs(dTe/Te2),abs(dTi/Ti2))
              niterimp=niterimp+1
              if(Ti2.le.0d0.or. Te2.le.0d0)then
                 oknonzero=.false.
              endif
           enddo
           if(oknonzero)then
              Ti=Ti2
              Te=Te2
              dtsum=dtsum+dt_loc
           endif
        enddo
        if(abs(Ti-Te).gt.abs(Ti2-Te2))then
           write(*,*)'Temperature difference has increased in coupling_fine, stopping...'
           write(*,*)'Ti=',Ti0,Ti2,myid
           write(*,*)'Te=',Te0,Te2,myid
           stop
        endif
        uold(ind_leaf(i),9)=Te2*Cv_electron
        errscalar=0.0d0
#if NENER>0
        do irad=1,nener
           errscalar=errscalar+uold(ind_leaf(i),8+irad)
        enddo
#endif
        uold(ind_leaf(i),5)=Ti2*Cv_ion+ekkscalar+emagscalar+errscalar

     enddo

  end do
  ! End loop over cells

end subroutine couplfine1



