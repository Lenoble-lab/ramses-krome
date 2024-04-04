!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine crdiff_fine(ilevel,compute)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the thermal conduction scheme.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,igrid,ncache,ngrid,compute,igroup
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
     call crdifffine1(ind_grid,ngrid,ilevel,compute,igroup)
  end do

111 format('   Entering conduction_fine for level ',i2)

end subroutine crdiff_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine crdifffine1(ind_grid,ncache,ilevel,compute,igroup)
  use amr_commons
  use hydro_commons
  use poisson_commons
!  use radiation_parameters,ONLY:dt_imp
  use cooling_module
  implicit none
  integer::ilevel,ncache,compute,igroup
  integer,dimension(1:nvector)::ind_grid
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
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  integer,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::indpass
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::upass
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::facdx
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux

  !real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::xloc=0.0d0
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2),save::residu_loc=0.0d0

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist
  real(dp),dimension(1:nvector),save:: residu

  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale,oneontwotondim,length,twodx,twodx2,bnorm2,DCRmax_code

  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_kappa
  real(dp)::vx,vy,vz,dens,bnorm,bx,by,bz,vnorm,va,Ma,Dpara,kperp,kpar,deCR
  real(dp)::deCRx,deCRy,deCRz,bdeCR
  integer::ind_cr1,l

!!$  scale_tau=1d0/(scale_t/scale_l**3) ! Time in code units

  ind_cr1=9
  if(twotemp)ind_cr1=10

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kappa = scale_l**2/scale_t
  kpar=Dcr/scale_kappa
  kpar=1.d10*flinj*boxlen*scale_l/scale_kappa ! Linj*c/3
  DCRmax_code=DCRmax/scale_kappa

  oneontwotondim = 1.d0/dble(twotondim)

  residu     = 0.0d0
  residu_loc = 0.0d0

  uloc =0.0d0
  upass=0.0d0

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  twodx =2.0d0*dx
  twodx2=twodx**2
 
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
     end if

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
        
        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nexist
              indpass(ind_exist(i),i3,j3,k3)=ind_cell(i) !! YD: test
              uloc (ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
              upass(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
              facdx(ind_exist(i),i3,j3,k3)=1.0d0

              if(compute==2 .and. ivar==1)then 
                 upass(ind_exist(i),i3,j3,k3,ivar)=unew(ind_cell(i),3) ! y of PBiCGSTAB
                 ! neighbor cell at ilevel+1, put value to zero in vector
                 ! because it is considered as a Dirichlet BC, i.e. in the RHS
                 if(son(ind_cell(i))>0) upass(ind_exist(i),i3,j3,k3,ivar)=0.0d0
              endif
              if(compute==4 .and. ivar==2)then 
                 upass(ind_exist(i),i3,j3,k3,ivar)=unew(ind_cell(i),8)  ! z of PBiCGSTAB
                 ! neighbor cell at ilevel+1, put value to zero in vector
                 ! because it is considered as a Dirichlet BC, i.e. in the RHS
                 if(son(ind_cell(i))>0) upass(ind_exist(i),i3,j3,k3,ivar)=0.0d0
              endif

           end do
           do i=1,nbuffer
              indpass(ind_nexist(i),i3,j3,k3)=ind_son !! YD: test
              uloc (ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
              upass(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
              ! neighbor cell at ilevel+1, put value to zero in vector
              ! because it is considered as a Dirichlet BC, i.e. in the RHS
              if(compute==2 .and. ivar==1)upass(ind_nexist(i),i3,j3,k3,ivar)=0.0
              if(compute==4 .and. ivar==2)upass(ind_nexist(i),i3,j3,k3,ivar)=0.0
           end do
        end do
        
        do i=1,nbuffer
           if(interpol_type_cond.gt.0)then
              facdx(ind_nexist(i),i3,j3,k3)=1.0d0
           else
              facdx(ind_nexist(i),i3,j3,k3)=1.5d0
           endif
        end do

        if(streaming_diffusion)then
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,5)=divu(ind_cell(i))              
           enddo
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,5)=uloc(ind_nexist(i),i3,j3,k3,igroup)
           enddo
        endif

        if(alfven_diff_coeff)then
           ! Compute Alfvenic Mach number = V/V_A and store it into uloc(:,3)
           do i=1,nexist
              dens   = uold(ind_cell(i),1)
              vx     = uold(ind_cell(i),2)/dens
              vy     = uold(ind_cell(i),3)/dens
              vz     = uold(ind_cell(i),4)/dens
              vnorm  = (vx**2+vy**2+vz**2)**0.5
              bx     = 0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
              by     = 0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
              bz     = 0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
              bnorm  = (bx**2+by**2+bz**2)**0.5
              
              va = bnorm/(dens**0.5)
              Ma = vnorm/va
              !==========================================================
              ! Ma>1 : lambda = Linj/Ma^3, isotrop, D=c*lambda/3
              ! Ma<1 : lambda_perp = lambda_para*Ma^4 and lambda_para=Linj/Ma^2
              !==========================================================
              if(Ma > 1)then
                 Dpara = kpar/(Ma**3)
                 kperp = 1.0d0
              else
                 Dpara = kpar*(Ma**2)
                 kperp=Ma**4
              end if
              upass(ind_exist(i),i3,j3,k3,3)=Dpara
              upass(ind_exist(i),i3,j3,k3,4)=kperp ! Dperp=kperp*Dpara
           end do
           do i=1,nbuffer
              dens= u2(i,ind_son,1)
              vx = u2(i,ind_son,2)/dens
              vy = u2(i,ind_son,3)/dens
              vz = u2(i,ind_son,4)/dens
              vnorm  = (vx**2+vy**2+vz**2)**0.5
              bx = 0.5*(u2(i,ind_son,6)+u2(i,ind_son,nvar+1))
              by = 0.5*(u2(i,ind_son,7)+u2(i,ind_son,nvar+2))
              bz = 0.5*(u2(i,ind_son,8)+u2(i,ind_son,nvar+3))
              bnorm = (bx**2+by**2+bz**2)**0.5
              va = bnorm/dens**0.5
              Ma = vnorm/va
              !==========================================================
              ! Ma>1 : lambda = Linj/Ma^3, isotrop, D=c*lambda/3
              ! Ma<1 : lambda_perp = lambda_para*Ma^4 and lambda_para=Linj/Ma^2
              !==========================================================
              if(Ma > 1)then
                 Dpara = kpar/Ma**3
                 kperp = 1.0d0
              else
                 Dpara = kpar*Ma**2
                 kperp=Ma**4
              end if
              
              upass(ind_nexist(i),i3,j3,k3,3)=Dpara
              upass(ind_nexist(i),i3,j3,k3,4)=kperp ! Dperp=kperp*Dpara
           end do
        endif

     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids


  if(streaming_diffusion)then

#if NDIM==1
     ! Yo: It does not work with NCR>1
     j1=1;k1=1
     do l=1,ncache
        do i1=0,3
           deCR=uloc(l,i1+1,j1,k1,ind_cr1)-uloc(l,i1-1,j1,k1,ind_cr1)
           dens   = uloc(l,i1,j1,k1,1)
           bx     = 0.5d0*(uloc(l,i1,j1,k1,6)+uloc(l,i1,j1,k1,nvar+1))
           by     = 0.5d0*(uloc(l,i1,j1,k1,7)+uloc(l,i1,j1,k1,nvar+2))
           bz     = 0.5d0*(uloc(l,i1,j1,k1,8)+uloc(l,i1,j1,k1,nvar+3))
           bnorm2  = bx**2+by**2+bz**2
           va = sqrt(bnorm2/dens)*fudge_streamboost
           bdeCR= abs(bx*deCR) / twodx / sqrt(bnorm2)

           deCR  =MAX(bdeCR,va*uloc(l,i1,j1,k1,ind_cr1)/DCRmax_code)
           length=MIN(uloc(l,i1,j1,k1,ind_cr1)/deCR,nlength_str*dx)
           upass(l,i1,j1,k1,3)=va*length*gamma_rad(1)
           upass(l,i1,j1,k1,4)=0.0d0
        enddo
     enddo
#endif     
#if NDIM==2
     ! Yo: It does not work with NCR>1
     k1=1
     do l=1,ncache
        do j1=0,3
        do i1=0,3
!!$           deCRx=uloc(l,i1+1,j1,k1,5)-uloc(l,i1-1,j1,k1,5)
!!$           deCRy=uloc(l,i1,j1+1,k1,5)-uloc(l,i1,j1-1,k1,5)
           deCRx=uloc(l,i1+1,j1,k1,ind_cr1)-uloc(l,i1-1,j1,k1,ind_cr1)
           deCRy=uloc(l,i1,j1+1,k1,ind_cr1)-uloc(l,i1,j1-1,k1,ind_cr1)
           dens   = uloc(l,i1,j1,k1,1)
           bx     = 0.5d0*(uloc(l,i1,j1,k1,6)+uloc(l,i1,j1,k1,nvar+1))
           by     = 0.5d0*(uloc(l,i1,j1,k1,7)+uloc(l,i1,j1,k1,nvar+2))
           bz     = 0.5d0*(uloc(l,i1,j1,k1,8)+uloc(l,i1,j1,k1,nvar+3))
           bnorm2  = bx**2+by**2+bz**2
           va = sqrt(bnorm2/dens)*fudge_streamboost
           bdeCR= abs(bx*deCRx + by*deCRy) / twodx / sqrt(bnorm2)
!!$           bdeCR= sqrt(deCRx*deCRx + deCRy*deCRy) / twodx

           deCR  =MAX(bdeCR,va*uloc(l,i1,j1,k1,ind_cr1)/DCRmax_code)
           length=MIN(uloc(l,i1,j1,k1,ind_cr1)/deCR,nlength_str*dx)
!!$           deCR  =MAX(bdeCR,va*uloc(l,i1,j1,k1,5)/DCRmax_code)
!!$           length=MIN(uloc(l,i1,j1,k1,5)/deCR,nlength_str*dx)
           upass(l,i1,j1,k1,3)=va*length*gamma_rad(1)
           upass(l,i1,j1,k1,4)=0.0d0
        enddo
        enddo
     enddo
#endif     
#if NDIM==3
     ! Yo: It does not work with NCR>1
     do l=1,ncache
        do k1=0,3
        do j1=0,3
        do i1=0,3
           deCRx=uloc(l,i1+1,j1  ,k1  ,ind_cr1)-uloc(l,i1-1,j1  ,k1  ,ind_cr1)
           deCRy=uloc(l,i1  ,j1+1,k1  ,ind_cr1)-uloc(l,i1  ,j1-1,k1  ,ind_cr1)
           deCRz=uloc(l,i1  ,j1  ,k1+1,ind_cr1)-uloc(l,i1  ,j1  ,k1-1,ind_cr1)
           dens   = uloc(l,i1,j1,k1,1)
           bx     = 0.5d0*(uloc(l,i1,j1,k1,6)+uloc(l,i1,j1,k1,nvar+1))
           by     = 0.5d0*(uloc(l,i1,j1,k1,7)+uloc(l,i1,j1,k1,nvar+2))
           bz     = 0.5d0*(uloc(l,i1,j1,k1,8)+uloc(l,i1,j1,k1,nvar+3))
           bnorm2  = bx**2+by**2+bz**2
           va = sqrt(bnorm2/dens)*fudge_streamboost
           bdeCR= abs(bx*deCRx + by*deCRy + bz*deCRz) / twodx / sqrt(bnorm2)
!!$           bdeCR= sqrt(deCRx**2 + deCRy**2 + deCRz**2) / twodx

           deCR  =MAX(bdeCR,va*uloc(l,i1,j1,k1,ind_cr1)/DCRmax_code)
           length=MIN(uloc(l,i1,j1,k1,ind_cr1)/deCR,nlength_str*dx)
           upass(l,i1,j1,k1,3)=va*length*gamma_rad(1)
           upass(l,i1,j1,k1,4)=0.0d0
        enddo
        enddo
        enddo
     enddo
#endif     

  endif

  !-----------------------------------------------
  ! Compute energy flux due to thermal conduction
  !-----------------------------------------------
  call crdiff_split(upass,flux,dx,dt_imp,ncache,compute,facdx,igroup)

  !-----------------------------------------------------
  ! update at level ilevel
  !-----------------------------------------------------
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
                 if(compute==0)then
                    residu_loc(i,i3   ,j3   ,k3   )=uold(ind_cell(i),igroup)
                 else if(compute==1)then
                    residu_loc(i,i3   ,j3   ,k3   )=0.0d0
                 else if(compute==2)then
                    residu_loc(i,i3   ,j3   ,k3   )=unew(ind_cell(i),3)
                 else if(compute==3)then
                    residu_loc(i,i3   ,j3   ,k3   )=1.0d0
                 else if(compute==4)then
                    residu_loc(i,i3   ,j3   ,k3   )=unew(ind_cell(i),8)
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
              if(compute==0)then
                 residu_loc(i,i3   ,j3   ,k3   )=residu_loc(i,i3   ,j3   ,k3   )+ &
                      & (flux(i,i3+i0,j3+j0,k3+k0,5,idim) &
                      & -flux(i,i3   ,j3   ,k3   ,5,idim))
              else
                 residu_loc(i,i3   ,j3   ,k3   )=residu_loc(i,i3   ,j3   ,k3   )+ &
                      & (flux(i,i3   ,j3   ,k3   ,5,idim) &
                      & -flux(i,i3+i0,j3+j0,k3+k0,5,idim))
              endif
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
                 if(compute==0)then
                    unew(ind_cell(i),igroup) = residu_loc(i,i3   ,j3   ,k3  ) ! Update the CR pressure
                 else if(compute==1)then
                    unew(ind_cell(i),1     ) = -residu_loc(i,i3   ,j3   ,k3  ) ! r_0 of PBiCGSTAB
                    unew(ind_cell(i),nvar+2) = -residu_loc(i,i3   ,j3   ,k3  ) ! rzero of PBiCGSTAB
                 else if(compute==2)then
                    unew(ind_cell(i),6     ) = residu_loc(i,i3   ,j3   ,k3   ) ! vi of PBiCGSTAB
                 else if(compute==3)then ! Compute 1/A
                    unew(ind_cell(i),4     ) = 1.0d0!/residu(i)
                 else if(compute==4)then
                    unew(ind_cell(i),nvar+1) = residu_loc(i,i3   ,j3   ,k3   ) ! t of PBiCGSTAB
                 end if
              endif
           end do
        end do
     end do
  end do


end subroutine crdifffine1
