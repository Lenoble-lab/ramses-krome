!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
! This subroutine (conduction_fine) is unused
!###########################################################
subroutine conduction_fine(ilevel,compute,itemp)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel,itemp
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the thermal conduction scheme.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,igrid,ncache,ngrid,compute
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call condfine1(ind_grid,ngrid,ilevel,compute,itemp)
  end do

111 format('   Entering conduction_fine for level ',i2)

end subroutine conduction_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine condfine1(ind_grid,ncache,ilevel,compute,itemp)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use cooling_module
  implicit none
  integer,intent(IN)::ilevel,ncache,compute,itemp
  integer,dimension(1:nvector),intent(IN)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! thermal conduction solver that compute energy flux. This flux is zeroed at 
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated 
  ! and stored in array unew(:), both at the current level and at the 
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     )::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       )::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         )::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         )::ind1
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3)::u2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::facdx
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux

  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::residu_loc=0.0d0

  integer,dimension(1:nvector)::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist
  logical,dimension(1:nvector)::okleaf
  real(dp),dimension(1:nvector):: residu

  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale
  real(dp)::kpara_ana,facdxamr

  real(dp)::Tau_ei_vol
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_tau,scale_kappa
  real(dp)::CviIntoCve,Cv1,Cv1_electron,Cv1_ion,ekin,emag,erad
  integer ::ind_eint,ind_Tenew,ind_Teold,ind_Cve

  real(dp)::d,u,v,w,bb
  real(dp)::oneoverd,oneoverCv1,oneoverCv1e,oneoverCv1i
#if NENER>0
  integer::irad
#endif
  
#ifndef WITHOUTMPI
  include 'mpif.h'
  real(dp)::tc1,tc2,tc3,tc4,tcneig1,tcneig2,tcfath1,tcfath2
  tc1=MPI_WTIME()
#endif

  ind_eint = nvar+2
  ind_Tenew = 5 
  ind_Cve = nvar+1
  ind_Teold = nvar+3
  CviIntoCve=mu_electron/mu_gas-1.0d0 ! 1.0d0
  if(interpol_type_cond.gt.0)then
     facdxamr=1.0d0
  else
     facdxamr=1.5d0
  endif

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_tau=1d0/(scale_t*scale_nH) ! Time in code units
  scale_kappa=scale_d*scale_l*scale_v**3

  call compute_Tau_equi_vol(Tau_ei_vol)
  if(.not.testcase)then
     Cv1         =kB/(mu_gas*mH*(gamma-1.0d0))/scale_v**2
     oneoverCv1=1d0/Cv1
     if(twotemp)then
        Cv1_electron=kB/(mu_electron*mH*(gamma_rad(1)-1.0d0))/scale_v**2
        Cv1_ion     =kB/(mu_ion     *mH*(gamma       -1.0d0))/scale_v**2        
        oneoverCv1e=1d0/Cv1_electron
        oneoverCv1i=1d0/Cv1_ion
     else
        Cv1_electron=kB/(mu_electron*mH*(gamma       -1.0d0))/scale_v**2 ! Not used in that case
        oneoverCv1e=1d0/Cv1_electron
     endif
  else
     if(twotemp)then
        Cv1         =2.0d0
        Cv1_electron=1.0d0
     else
        Cv1         =1.0d0
        Cv1_electron=Cv1
     endif
     oneoverCv1 =1d0/Cv1
     oneoverCv1e=1d0/Cv1_electron
  endif

  residu     = 0.0d0
  residu_loc = 0.0d0

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
 
 ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     
     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
           nbuffer=nbuffer+1
           ind_nexist(nbuffer)=i
           ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do

     ! If not, interpolate variables from parent cells
     if(nbuffer>0)then
#ifndef WITHOUTMPI
        tcneig1=MPI_WTIME()
#endif
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar+3
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_hydro_cond(u1,ind1,u2,nbuffer)
#ifndef WITHOUTMPI
        tcneig2=MPI_WTIME()
        dtcneig=dtcneig+tcneig2-tcneig1
#endif
     end if

#ifndef WITHOUTMPI
     tcfath1=MPI_WTIME()
#endif
     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do
        
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2
        
        okleaf=.false.
        do i=1,nexist
           if(son(ind_cell(i))==0)okleaf(i)=.true.
        enddo

        do i=1,nexist
           uloc(ind_exist(i),i3,j3,k3,1:nvar+3)=uold(ind_cell(i),1:nvar+3)
           facdx(ind_exist(i),i3,j3,k3)=1.0d0
        enddo

        ! neighbor cell at ilevel+1, put value to zero in vector because it is considered as a Dirichlet BC, i.e. in the RHS
        ! ivar=2
        if(compute==2)then 
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,2)=unew(ind_cell(i),2)
              if(.not.okleaf(i)) uloc(ind_exist(i),i3,j3,k3,2)=0.0d0
           enddo
        endif
        ! ivar=3&4
        if(.not.twotemp)then
           do i=1,nexist
              if(okleaf(i))then
                 uloc(ind_exist(i),i3,j3,k3,3)=unew(ind_cell(i),ind_Teold)
                 uloc(ind_exist(i),i3,j3,k3,4)=divu(ind_cell(i))
              else
                 d=uold(ind_cell(i),1)
                 u=uold(ind_cell(i),2)
                 v=uold(ind_cell(i),3)
                 w=uold(ind_cell(i),4)
                 oneoverd=1d0/d
                 ekin = 0.5d0*(u*u+v*v+w*w)*oneoverd
                 emag=0.0d0
                 do idim=1,3
                    bb=uold(ind_cell(i),5+idim)+uold(ind_cell(i),nvar+idim)
                    emag=emag+bb*bb
                 end do
                 emag=emag*0.125d0
                 erad=0.0d0
                 ! Neighbor more refined, compute temperature from total energy
#if NENER>0
                 do irad=1,nener
                    erad=erad+uold(ind_cell(i),8+irad)
                 end do
#endif
                 uloc(ind_exist(i),i3,j3,k3,3)=(uold(ind_cell(i),5)-ekin-emag-erad)*oneoverd*oneoverCv1
                 ! Neighbor more refined, compute conductivity from temperature
                 uloc(ind_exist(i),i3,j3,k3,4)=kpara_ana(uloc(ind_exist(i),i3,j3,k3,3),itemp)/scale_kappa
              endif
           enddo
        else
           if(itemp==1)then
              do i=1,nexist
                 if(okleaf(i))then
                    uloc(ind_exist(i),i3,j3,k3,3)=unew(ind_cell(i),ind_Teold)
                    uloc(ind_exist(i),i3,j3,k3,4)=divu(ind_cell(i))
                 else
                    ! Neighbor more refined, compute conductivity from temperature
                    uloc(ind_exist(i),i3,j3,k3,3)=uold(ind_cell(i),9)/uold(ind_cell(i),1)*oneoverCv1e
                    uloc(ind_exist(i),i3,j3,k3,4)=kpara_ana(uloc(ind_exist(i),i3,j3,k3,3),itemp)/scale_kappa
                 endif
              enddo
           else
              do i=1,nexist
                 if(okleaf(i))then
                    uloc(ind_exist(i),i3,j3,k3,3)=unew(ind_cell(i),ind_Teold)
                    uloc(ind_exist(i),i3,j3,k3,4)=divu(ind_cell(i))
                 else
                    d=uold(ind_cell(i),1)
                    u=uold(ind_cell(i),2)
                    v=uold(ind_cell(i),3)
                    w=uold(ind_cell(i),4)
                    oneoverd=1d0/d
                    ekin = 0.5d0*(u*u+v*v+w*w)*oneoverd
                    emag=0.0d0
                    do idim=1,3
                       bb=uold(ind_cell(i),5+idim)+uold(ind_cell(i),nvar+idim)
                       emag=emag+bb*bb
                    end do
                    emag=emag*0.125d0
                    erad=0.0d0
                    ! Neighbor more refined, compute temperature from total energy
#if NENER>1
                    do irad=2,nener
                       erad=erad+uold(ind_cell(i),8+irad)
                    end do
#endif
                    uloc(ind_exist(i),i3,j3,k3,3)=(uold(ind_cell(i),5)-ekin-emag-erad-uold(ind_cell(i),9)) &
                         & *oneoverd*oneoverCv1i
                    ! Neighbor more refined, compute conductivity from temperature
                    uloc(ind_exist(i),i3,j3,k3,4)=kpara_ana(uloc(ind_exist(i),i3,j3,k3,3),itemp)/scale_kappa
                 endif
              enddo
           endif
        endif

        do i=1,nbuffer
           facdx(ind_nexist(i),i3,j3,k3)=facdxamr
        enddo
        do ivar=1,nvar+3
           do i=1,nbuffer
              uloc (ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           enddo
        enddo
        do i=1,nbuffer
           d=u2(i,ind_son,1)
           u=u2(i,ind_son,2)
           v=u2(i,ind_son,3)
           w=u2(i,ind_son,4)
           oneoverd=1d0/d
           ekin = 0.5d0*(u*u+v*v+w*w)*oneoverd
           emag=0.0d0
           do idim=1,3
              bb=u2(i,ind_son,5+idim)+u2(i,ind_son,nvar+idim)
              emag=emag+bb*bb
           end do
           emag=emag*0.125d0
           erad=0.0d0
           if(.not.twotemp)then
#if NENER>0
              do irad=1,nener
                 erad=erad+u2(i,ind_son,8+irad)
              end do
#endif               
              uloc(ind_nexist(i),i3,j3,k3,3)=(u2(i,ind_son,5)-ekin-emag-erad)*oneoverd*oneoverCv1
           else
#if NENER>1
              do irad=2,nener
                 erad=erad+u2(i,ind_son,8+irad)
              end do
#endif               
              if(itemp==1)then
                 uloc(ind_nexist(i),i3,j3,k3,3)=u2(i,ind_son,9)*oneoverd*oneoverCv1e
              else
                 uloc(ind_nexist(i),i3,j3,k3,3)=(u2(i,ind_son,5)-ekin-emag-erad-u2(i,ind_son,9)) &
                      & *oneoverd*oneoverCv1i
              endif
           end if
           ! Compute conductivity from temperature
           uloc(ind_nexist(i),i3,j3,k3,4)=kpara_ana(uloc(ind_nexist(i),i3,j3,k3,3),itemp)/scale_kappa
           ! neighbor cell at ilevel-1, put value to zero in vector because it is considered as a Dirichlet BC, i.e. in the RHS
           if(compute==2)uloc(ind_nexist(i),i3,j3,k3,2)=0.0
        enddo
        
     end do
     end do
     end do
     ! End loop over cells
#ifndef WITHOUTMPI
     tcfath2=MPI_WTIME()
     dtcfath=dtcfath+tcfath2-tcfath1
#endif

  end do
  end do
  end do
  ! End loop over neighboring grids

#ifndef WITHOUTMPI
  tc2=MPI_WTIME()
  dtc12=dtc12+tc2-tc1
#endif
  !-----------------------------------------------
  ! Compute energy flux due to thermal conduction
  !-----------------------------------------------
  if(semi_implicit)then
     call cond_split_semi(uloc,flux,dx,dt_imp,ncache,compute,facdx)
  else
     call cond_split(uloc,flux,dx,dt_imp,ncache,compute,facdx,itemp)
  endif
#ifndef WITHOUTMPI
  tc3=MPI_WTIME()
  dtc23=dtc23+tc3-tc2
#endif

  !-----------------------------------------------------
  ! update at level ilevel
  !-----------------------------------------------------
  i0=0; j0=0; k0=0
  if(idim==1)i0=1
  if(idim==2)j0=1
  if(idim==3)k0=1
  do k2=k2min,k2max
     do j2=j2min,j2max
        do i2=i2min,i2max
           ind_son=1+i2+2*j2+4*k2
           iskip=ncoarse+(ind_son-1)*ngridmax
           do i=1,ncache
              ind_cell(i)=iskip+ind_grid(i)
           end do
           i3=1+i2
           j3=1+j2
           k3=1+k2
           ! Update conservative variables new state vector
           do i=1,ncache
              if(son(ind_cell(i))==0)then        
                 
                 if(compute==1)then
                    residu_loc(i,i3,j3,k3)= 0.0d0
                 else if(compute==2)then ! compute Ap 
                    residu_loc(i,i3,j3,k3)=(unew(ind_cell(i),ind_Cve))*unew(ind_cell(i),2)
                 else if(compute==3)then
                    residu_loc(i,i3,j3,k3)=(unew(ind_cell(i),ind_Cve))
                 else if(compute==0)then
                    residu_loc(i,i3,j3,k3)=unew(ind_cell(i),ind_Cve)*unew(ind_cell(i),ind_Teold)
                 end if
              
              endif
           end do
        end do
     end do
  end do

  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        ! Update conservative variables new state vector
        do i=1,ncache
           if(son(ind_cell(i))==0)then
              residu_loc(i,i3   ,j3   ,k3   )=residu_loc(i,i3   ,j3   ,k3   )+ &
                   & (flux(i,i3   ,j3   ,k3   ,5,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,5,idim))
           endif
        end do
     end do
     end do
     end do
  end do

  i0=0; j0=0; k0=0
  if(idim==1)i0=1
  if(idim==2)j0=1
  if(idim==3)k0=1
  do k2=k2min,k2max
     do j2=j2min,j2max
        do i2=i2min,i2max
           ind_son=1+i2+2*j2+4*k2
           iskip=ncoarse+(ind_son-1)*ngridmax
           do i=1,ncache
              ind_cell(i)=iskip+ind_grid(i)
           end do
           i3=1+i2
           j3=1+j2
           k3=1+k2
           ! Update conservative variables new state vector
           do i=1,ncache
              if(son(ind_cell(i))==0)then
                 
                 if(compute==1)then
                    unew(ind_cell(i),1) = -residu_loc(i,i3   ,j3   ,k3  ) ! r0
                    unew(ind_cell(i),2) = -residu_loc(i,i3   ,j3   ,k3  ) ! p0
                 else if(compute==2)then ! compute Ap 
                    unew(ind_cell(i),3) = residu_loc(i,i3   ,j3   ,k3   )
                 else if(compute==3)then ! Compute 1/A
                    unew(ind_cell(i),4) = 1.0d0/residu_loc(i,i3   ,j3   ,k3   )
                 else if(compute==0)then
                    unew(ind_cell(i),ind_eint) = residu_loc(i,i3   ,j3   ,k3   )
                 end if
              endif
           end do
        end do
     end do
  end do
#ifndef WITHOUTMPI
  tc4=MPI_WTIME()
  dtc34=dtc34+tc4-tc3
#endif
end subroutine condfine1
