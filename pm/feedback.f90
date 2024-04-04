!################################################################
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine thermal_feedback(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info2,dummy_io
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and
  ! the metal mass dumped in the gas by stars (SNII, SNIa, winds).
  ! This routine is called every fine time step.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ivar
  integer::ig,ip,npart1,npart2,icpu,ilun,idim
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  character(LEN=80)::filename,filedir,fileloc,filedirini
  character(LEN=5)::nchar,ncharcpu
  logical::file_exist
  integer,parameter::tag=1120

  if(sf_log_properties) then
     call title(ifout-1,nchar)
     if(IOGROUPSIZEREP>0) then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        filedirini='output_'//TRIM(nchar)//'/'
        filedir='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/'
     else
        filedir='output_'//TRIM(nchar)//'/'
     endif
     filename=TRIM(filedir)//'stars_'//TRIM(nchar)//'.out'
     ilun=myid+103
     call title(myid,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     inquire(file=fileloc,exist=file_exist)
     if(.not.file_exist) then
        open(ilun, file=fileloc, form='formatted')
        write(ilun,'(A24)',advance='no') '# event id  ilevel  mp  '
        do idim=1,ndim
           write(ilun,'(A2,I1,A2)',advance='no') 'xp',idim,'  '
        enddo
        do idim=1,ndim
           write(ilun,'(A2,I1,A2)',advance='no') 'vp',idim,'  '
        enddo
        do ivar=1,nvar
           if(ivar.ge.10) then
              write(ilun,'(A1,I2,A2)',advance='no') 'u',ivar,'  '
           else
              write(ilun,'(A1,I1,A2)',advance='no') 'u',ivar,'  '
           endif
        enddo
        write(ilun,'(A5)',advance='no') 'tag  '
        write(ilun,'(A1)') ' '
     else
        open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
     endif
  endif


  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Gather star particles only.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if ( is_star(typep(ipart)) ) then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        ! Gather star particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only star particles
              if ( is_star(typep(ipart)) ) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

  if(sf_log_properties) close(ilun)

111 format('   Entering thermal_feedback for level ',I2)

end subroutine thermal_feedback
#endif
!################################################################
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine feedbk(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use random
  use constants, only: M_sun, Myr2sec, pc2cm
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine feedback. Each stellar particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! unew.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ivar,ilun
  real(kind=8)::RandNum
  real(dp)::SN_BOOST,mstar,dx_min,vol_min
  real(dp)::t0,ESN,mejecta,zloss,e,uvar
  real(dp)::ERAD,RAD_BOOST,tauIR,msne_min,mstar_max,eta_sn2
  real(dp)::delta_x,tau_factor,rad_factor
  real(dp)::dx,dx_loc,scale,birth_time,current_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,ekinetic,dteff
  real(dp),dimension(1:nvector),save::vol_loc
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
#if NENER>0
  integer::irad
#endif

  if(trim(feedback_model).eq.'agertz') then
     call agertz_feedbk(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
     return
  endif

  if(sf_log_properties) ilun=myid+103
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Minimum star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif
  msne_min=mass_sne_min*M_sun/(scale_d*scale_l**3)
  mstar_max=mass_star_max*M_sun/(scale_d*scale_l**3)

  ! Compute stochastic boost to account for target GMC mass
  SN_BOOST=MAX(mass_gmc*M_sun/(scale_d*scale_l**3)/mstar,1d0)

  ! Massive star lifetime from Myr to code units
  if(use_proper_time)then
     t0=t_sne*Myr2sec/(scale_t/aexp**2)
     current_time=texp
  else
     t0=t_sne*Myr2sec/scale_t
     current_time=t
  endif

  ! Type II supernova specific energy from cgs to code units
  ESN=1d51/(10d0*M_sun)/scale_v**2

  ! Life time radiation specific energy from cgs to code units
  ERAD=1d53/(10d0*M_sun)/scale_v**2

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=int(x(j,idim))
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
        vol_loc(j)=vol_loc(j)*2**ndim ! ilevel-1 cell volume
     end if
  end do

  ! Compute individual time steps
  do j=1,np
     dteff(j)=dtnew(levelp(ind_part(j)))
  end do

  if(use_proper_time)then
     do j=1,np
        dteff(j)=dteff(j)*aexp**2
     end do
  endif

  ! Reset ejected mass, metallicity, thermal energy
  do j=1,np
     mloss(j)=0d0
     mzloss(j)=0d0
     ethermal(j)=0d0
  end do

  ! Compute stellar mass loss and thermal feedback due to supernovae
  if(f_w==0)then
     do j=1,np
        birth_time=tp(ind_part(j))
        ! Make sure that we don't count feedback twice
        if(birth_time.lt.(current_time-t0).and.birth_time.ge.(current_time-t0-dteff(j)))then
           eta_sn2   = eta_sn
           if(sf_imf)then
              if(mp(ind_part(j)).le.mstar_max)then
                 if(mp(ind_part(j)).ge.msne_min) eta_sn2 = eta_ssn
                 if(mp(ind_part(j)).lt.msne_min) eta_sn2 = 0
              endif
           endif
           ! Stellar mass loss
           mejecta=eta_sn2*mp(ind_part(j))
           mloss(j)=mloss(j)+mejecta/vol_loc(j)
           ! Thermal energy
           ethermal(j)=ethermal(j)+mejecta*ESN/vol_loc(j)
           ! Metallicity
           if(metal)then
              zloss=yield+(1d0-yield)*zp(ind_part(j))
              mzloss(j)=mzloss(j)+mejecta*zloss/vol_loc(j)
           endif
           ! Reduce star particle mass
           mp(ind_part(j))=mp(ind_part(j))-mejecta
           ! Boost SNII energy and depopulate accordingly
           if(SN_BOOST>1d0)then
              call ranf(localseed,RandNum)
              if(RandNum<1d0/SN_BOOST)then
                 mloss(j)=SN_BOOST*mloss(j)
                 mzloss(j)=SN_BOOST*mzloss(j)
                 ethermal(j)=SN_BOOST*ethermal(j)
              else
                 mloss(j)=0d0
                 mzloss(j)=0d0
                 ethermal(j)=0d0
              endif
           endif
           if(sf_log_properties) then
              write(ilun,'(I10)',advance='no') 1
              write(ilun,'(2I10,E24.12)',advance='no') idp(ind_part(j)),ilevel,mp(ind_part(j))
              do idim=1,ndim
                 write(ilun,'(E24.12)',advance='no') xp(ind_part(j),idim)
              enddo
              do idim=1,ndim
                 write(ilun,'(E24.12)',advance='no') vp(ind_part(j),idim)
              enddo
              write(ilun,'(E24.12)',advance='no') unew(indp(j),1)
              do ivar=2,nvar
                 if(ivar.eq.ndim+2)then
                    e=0.0d0
                    do idim=1,ndim
                       e=e+0.5d0*unew(ind_cell(i),idim+1)**2/max(unew(ind_cell(i),1),smallr)
                    enddo
#if NENER>0
                    do irad=0,nener-1
                       e=e+unew(ind_cell(i),inener+irad)
                    enddo
#endif
#ifdef SOLVERmhd
                    do idim=1,ndim
                       e=e+0.125d0*(unew(ind_cell(i),idim+ndim+2)+unew(ind_cell(i),idim+nvar))**2
                    enddo
#endif
                    ! Temperature
                    uvar=(gamma-1.0d0)*(unew(ind_cell(i),ndim+2)-e)*scale_T2
                 else
                    uvar=unew(indp(j),ivar)
                 endif
                 write(ilun,'(E24.12)',advance='no') uvar/unew(indp(j),1)
              enddo
              write(ilun,'(I10)',advance='no') typep(ind_part(i))%tag
              write(ilun,'(A1)') ' '
           endif
        endif
     end do
  endif

  ! Update hydro variables due to feedback

  ! For IR radiation trapping,
  ! we use a fixed length to estimate the column density of gas
  delta_x=200d0*pc2cm
  if(metal)then
     tau_factor=kappa_IR*delta_x*scale_d/0.02d0
  else
     tau_factor=kappa_IR*delta_x*scale_d*z_ave
  endif
  rad_factor=ERAD/ESN

  do j=1,np

     ! Infrared photon trapping boost
     if(metal)then
        tauIR=tau_factor*max(uold(indp(j),imetal),smallr)
     else
        tauIR=tau_factor*max(uold(indp(j),1),smallr)
     endif
     if(uold(indp(j),1)*scale_nH > 10.)then
        RAD_BOOST=rad_factor*(1d0-exp(-tauIR))
     else
        RAD_BOOST=0
     endif

     ! Specific kinetic energy of the star
     ekinetic(j)=0.5d0*(vp(ind_part(j),1)**2 &
          &            +vp(ind_part(j),2)**2 &
          &            +vp(ind_part(j),3)**2)

     ! Update hydro variable in NGP cell
     unew(indp(j),1)=unew(indp(j),1)+mloss(j)
     unew(indp(j),2)=unew(indp(j),2)+mloss(j)*vp(ind_part(j),1)
     unew(indp(j),3)=unew(indp(j),3)+mloss(j)*vp(ind_part(j),2)
     unew(indp(j),4)=unew(indp(j),4)+mloss(j)*vp(ind_part(j),3)
     unew(indp(j),5)=unew(indp(j),5)+mloss(j)*ekinetic(j)+ &
          & ethermal(j)*(1d0+RAD_BOOST)
#if NENER>0
     ! Cosmic rays injection!!
     unew(indp(j),inener)=unew(indp(j),inener)+ethermal(j)*(1d0+RAD_BOOST)*fecr
#endif                  
  end do

  ! Add metals
  if(metal)then
     do j=1,np
        unew(indp(j),imetal)=unew(indp(j),imetal)+mzloss(j)
     end do
  endif

  ! Add delayed cooling switch variable
  if(delayed_cooling)then
     do j=1,np
        unew(indp(j),idelay)=unew(indp(j),idelay)+mloss(j)
     end do
  endif

end subroutine feedbk
#endif
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kinetic_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  use constants, only:Myr2sec
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,sSN_all,ZSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all
#endif
  !----------------------------------------------------------------------
  ! This subroutine compute the kinetic feedback due to SNII and
  ! imolement this using exploding GMC particles.
  ! This routine is called only at coarse time step.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::ip,icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_loc,nSN_tot,iSN,ilevel,ivar
  integer,dimension(1:ncpu)::nSN_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0
  real(dp)::current_time
  real(dp)::scale,dx_min,vol_min,mstar
  integer::nx_loc
  integer,dimension(:),allocatable::ind_part,ind_grid
  logical,dimension(:),allocatable::ok_free
  integer,dimension(:),allocatable::indSN
  real(dp),dimension(:),allocatable::mSN,sSN,ZSN,m_gas,vol_gas,ekBlast
  real(dp),dimension(:,:),allocatable::xSN,vSN,u_gas,dq

  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)'Entering make_sn'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Minimum star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif


  ! Lifetime of Giant Molecular Clouds from Myr to code units
  ! Massive star lifetime from Myr to code units
  if(use_proper_time)then
     t0=t_sne*Myr2sec/(scale_t/aexp**2)
     current_time=texp
  else
     t0=t_sne*Myr2sec/scale_t
     current_time=t
  endif

  !------------------------------------------------------
  ! Gather GMC particles eligible for disruption
  !------------------------------------------------------
  nSN_loc=0
  ! Loop over levels
  do icpu=1,ncpu
  ! Loop over cpus
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count old enough GMC particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if ( is_debris(typep(ipart)) .and. tp(ipart).lt.(current_time-t0) ) then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
            end do
        endif
        nSN_loc=nSN_loc+npart2   ! Add SNe to the total
        igrid=next(igrid)   ! Go to next grid
     end do
  end do
  ! End loop over levels
  nSN_icpu=0
  nSN_icpu(myid)=nSN_loc
#ifndef WITHOUTMPI
  ! Give an array of number of SN on each cpu available to all cpus
  call MPI_ALLREDUCE(nSN_icpu,nSN_icpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_icpu=nSN_icpu_all
#endif

  nSN_tot=sum(nSN_icpu(1:ncpu))

  if (nSN_tot .eq. 0) return

  if(myid==1)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of GMC to explode=',nSN_tot
     write(*,*)'-----------------------------------------------'
  endif

  ! Allocate arrays for the position and the mass of the SN
  allocate(xSN(1:nSN_tot,1:3),vSN(1:nSN_tot,1:3))
  allocate(mSN(1:nSN_tot),sSN(1:nSN_tot),ZSN(1:nSN_tot))
  xSN=0; vSN=0; mSN=0; sSN=0; ZSN=0
  ! Allocate arrays for particles index and parent grid
  if(nSN_loc>0)then
     allocate(ind_part(1:nSN_loc),ind_grid(1:nSN_loc),ok_free(1:nSN_loc))
  endif

  !------------------------------------------------------
  ! Store position and mass of the GMC into the SN array
  !------------------------------------------------------
  if(myid==1)then
     iSN=0
  else
     iSN=sum(nSN_icpu(1:myid-1))
  endif
  ! Loop over levels
  ip=0
  do icpu=1,ncpu
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        ! Count old enough star particles that have not exploded
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if ( is_debris(typep(ipart)) .and. tp(ipart).lt.(current_time-t0) ) then
                 iSN=iSN+1
                 xSN(iSN,1)=xp(ipart,1)
                 xSN(iSN,2)=xp(ipart,2)
                 xSN(iSN,3)=xp(ipart,3)
                 vSN(iSN,1)=vp(ipart,1)
                 vSN(iSN,2)=vp(ipart,2)
                 vSN(iSN,3)=vp(ipart,3)
                 mSN(iSN)=mp(ipart)
                 sSN(iSN)=dble(-idp(ipart))*mstar
                 if(metal)ZSN(iSN)=zp(ipart)
                 ip=ip+1
                 ind_grid(ip)=igrid
                 ind_part(ip)=ipart
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        igrid=next(igrid)   ! Go to next grid
     end do
  end do
  ! End loop over levels

  ! Remove GMC particle
  if(nSN_loc>0)then
     ok_free=.true.
     call remove_list(ind_part,ind_grid,ok_free,nSN_loc)
     call add_free_cond(ind_part,ok_free,nSN_loc)
     deallocate(ind_part,ind_grid,ok_free)
  endif

#ifndef WITHOUTMPI
  allocate(xSN_all(1:nSN_tot,1:3),vSN_all(1:nSN_tot,1:3),mSN_all(1:nSN_tot),sSN_all(1:nSN_tot),ZSN_all(1:nSN_tot))
  call MPI_ALLREDUCE(xSN,xSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vSN,vSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mSN,mSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(sSN,sSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ZSN,ZSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  xSN=xSN_all
  vSN=vSN_all
  mSN=mSN_all
  sSN=sSN_all
  ZSN=ZSN_all
  deallocate(xSN_all,vSN_all,mSN_all,sSN_all,ZSN_all)
#endif

  nSN=nSN_tot
  allocate(m_gas(1:nSN),u_gas(1:nSN,1:3),vol_gas(1:nSN),dq(1:nSN,1:3),ekBlast(1:nSN))
  allocate(indSN(1:nSN))

  ! Compute the grid discretization effects
  call average_SN(xSN,vol_gas,dq,ekBlast,indSN,nSN)

  ! Modify hydro quantities to account for a Sedov blast wave
  call Sedov_blast(xSN,vSN,mSN,sSN,ZSN,indSN,vol_gas,dq,ekBlast,nSN)

  deallocate(xSN,vSN,mSN,sSN,ZSN,indSN,m_gas,u_gas,vol_gas,dq,ekBlast)

  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo

end subroutine kinetic_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_SN(xSN,vol_gas,dq,ekBlast,ind_blast,nSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  use constants, only: pc2cm
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,iSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_SN,u,v,w,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nSN)::ind_blast
  real(dp),dimension(1:nSN)::vol_gas,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,dq,u2Blast
#ifndef WITHOUTMPI
  real(dp),dimension(1:nSN)::vol_gas_all,ekBlast_all
  real(dp),dimension(1:nSN,1:3)::dq_all,u2Blast_all
#endif
  logical ,dimension(1:nvector),save::ok

  if(nSN==0)return
  if(verbose)write(*,*)'Entering average_SN'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(2.0d0*dx_min*scale_l/aexp,rbubble*pc2cm)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Initialize the averaged variables
  vol_gas=0; dq=0; u2Blast=0; ekBlast=0; ind_blast=-1

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities)
     dx=0.5D0**ilevel
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                    if(dr_SN.lt.rmax2)then
                       vol_gas(iSN)=vol_gas(iSN)+vol_loc
                       ! Take account for grid effects on the conservation of the
                       ! normalized linear momentum
                       u=dxx/rmax
                       v=dyy/rmax
                       w=dzz/rmax
                       ! Add the local normalized linear momentum to the total linear
                       ! momentum of the blast wave (should be zero with no grid effect)
                       dq(iSN,1)=dq(iSN,1)+u*vol_loc
                       dq(iSN,2)=dq(iSN,2)+v*vol_loc
                       dq(iSN,3)=dq(iSN,3)+w*vol_loc
                       u2Blast(iSN,1)=u2Blast(iSN,1)+u*u*vol_loc
                       u2Blast(iSN,2)=u2Blast(iSN,2)+v*v*vol_loc
                       u2Blast(iSN,3)=u2Blast(iSN,3)+w*w*vol_loc
                    endif
                    if(dr_cell.le.dx_loc/2.0)then
                       ind_blast(iSN)=ind_cell(i)
                       ekBlast  (iSN)=vol_loc
                    endif
                 end do
              endif
           end do

        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(vol_gas,vol_gas_all,nSN  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dq     ,dq_all     ,nSN*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(u2Blast,u2Blast_all,nSN*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ekBlast,ekBlast_all,nSN  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas=vol_gas_all
  dq     =dq_all
  u2Blast=u2Blast_all
  ekBlast=ekBlast_all
#endif
  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        dq(iSN,1)=dq(iSN,1)/vol_gas(iSN)
        dq(iSN,2)=dq(iSN,2)/vol_gas(iSN)
        dq(iSN,3)=dq(iSN,3)/vol_gas(iSN)
        u2Blast(iSN,1)=u2Blast(iSN,1)/vol_gas(iSN)
        u2Blast(iSN,2)=u2Blast(iSN,2)/vol_gas(iSN)
        u2Blast(iSN,3)=u2Blast(iSN,3)/vol_gas(iSN)
        u2=u2Blast(iSN,1)-dq(iSN,1)**2
        v2=u2Blast(iSN,2)-dq(iSN,2)**2
        w2=u2Blast(iSN,3)-dq(iSN,3)**2
        ekBlast(iSN)=max(0.5d0*(u2+v2+w2),0.0d0)
     endif
  end do

  if(verbose)write(*,*)'Exiting average_SN'

end subroutine average_SN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine Sedov_blast(xSN,vSN,mSN,sSN,ZSN,indSN,vol_gas,dq,ekBlast,nSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  use constants, only: M_sun, pc2cm
  use mpi_mod
  implicit none
  !------------------------------------------------------------------------
  ! This routine merges SN using the FOF algorithm.
  !------------------------------------------------------------------------
  integer::ilevel,iSN,nSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_SN,u,v,w,ESN,mstar,eta_sn2,msne_min,mstar_max
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax,vol_min
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nSN)::mSN,sSN,ZSN,p_gas,d_gas,d_metal,vol_gas,uSedov,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,vSN,dq
  integer ,dimension(1:nSN)::indSN
  logical ,dimension(1:nvector),save::ok

  if(nSN==0)return
  if(verbose)write(*,*)'Entering Sedov_blast'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax
  vol_min=dx_min**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(2.0d0*dx_min*scale_l/aexp,rbubble*pc2cm)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Minimum star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif
  msne_min=mass_sne_min*M_sun/(scale_d*scale_l**3)
  mstar_max=mass_star_max*M_sun/(scale_d*scale_l**3)
  ! Supernova specific energy from cgs to code units
  ESN=(1d51/(10d0*M_sun))/scale_v**2

  do iSN=1,nSN
     eta_sn2    = eta_sn
     if(sf_imf)then
        if(mSN(iSN).le.mstar_max)then
           if(mSN(iSN).ge.msne_min) eta_sn2 = eta_ssn
           if(mSN(iSN).lt.msne_min) eta_sn2 = 0
        endif
     endif
     if(vol_gas(iSN)>0d0)then
        d_gas(iSN)=mSN(iSN)/vol_gas(iSN)
        if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/vol_gas(iSN)
        if(ekBlast(iSN)==0d0)then
           p_gas(iSN)=eta_sn2*sSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=0d0
        else
           p_gas(iSN)=(1d0-f_ek)*eta_sn2*sSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=sqrt(f_ek*eta_sn2*sSN(iSN)*ESN/mSN(iSN)/ekBlast(iSN))
        endif
     else
        d_gas(iSN)=mSN(iSN)/ekBlast(iSN)
        p_gas(iSN)=eta_sn2*sSN(iSN)*ESN/ekBlast(iSN)
        if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/ekBlast(iSN)
     endif
  end do

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities)
     dx=0.5D0**ilevel
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    if(dr_SN.lt.rmax2)then
                       ! Compute the mass density in the cell
                       uold(ind_cell(i),1)=uold(ind_cell(i),1)+d_gas(iSN)
                       ! Compute the metal density in the cell
                       if(metal)uold(ind_cell(i),imetal)=uold(ind_cell(i),imetal)+d_metal(iSN)
                       ! Velocity at a given dr_SN linearly interpolated between zero and uSedov
                       u=uSedov(iSN)*(dxx/rmax-dq(iSN,1))+vSN(iSN,1)
                       v=uSedov(iSN)*(dyy/rmax-dq(iSN,2))+vSN(iSN,2)
                       w=uSedov(iSN)*(dzz/rmax-dq(iSN,3))+vSN(iSN,3)
                       ! Add each momentum component of the blast wave to the gas
                       uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas(iSN)*u
                       uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas(iSN)*v
                       uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas(iSN)*w
                       ! Finally update the total energy of the gas
                       uold(ind_cell(i),5)=uold(ind_cell(i),5)+0.5d0*d_gas(iSN)*(u*u+v*v+w*w)+p_gas(iSN)
                    endif
                 end do
              endif
           end do

        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  do iSN=1,nSN
     if(vol_gas(iSN)==0d0)then
        u=vSN(iSN,1)
        v=vSN(iSN,2)
        w=vSN(iSN,3)
        if(indSN(iSN)>0)then
           uold(indSN(iSN),1)=uold(indSN(iSN),1)+d_gas(iSN)
           uold(indSN(iSN),2)=uold(indSN(iSN),2)+d_gas(iSN)*u
           uold(indSN(iSN),3)=uold(indSN(iSN),3)+d_gas(iSN)*v
           uold(indSN(iSN),4)=uold(indSN(iSN),4)+d_gas(iSN)*w
           uold(indSN(iSN),5)=uold(indSN(iSN),5)+d_gas(iSN)*0.5d0*(u*u+v*v+w*w)+p_gas(iSN)
           if(metal)uold(indSN(iSN),imetal)=uold(indSN(iSN),imetal)+d_metal(iSN)
        endif
     endif
  end do

  if(verbose)write(*,*)'Exiting Sedov_blast'

end subroutine Sedov_blast
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM==3
subroutine agertz_feedbk(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use random
  use constants, only: M_sun, Myr2sec, pc2cm
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine feedback. Each stellar particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! unew.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,ilun,iii,ii,jj,kk
  real(kind=8)::RandNum
  real(dp)::dx_min,vol_min
  real(dp)::t0,ESN,mejecta,zloss,time_simu
  real(dp)::tauIR
  real(dp)::dx,dx_loc,scale,birth_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::L1,Cr,KappaIR,KappaIR_0
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,dteff
  real(dp),dimension(1:nvector)::mlossFe,mlossO
  real(dp),dimension(1:nvector)::ptot
  real(dp),dimension(1:nvector),save::vol_loc,dx_loc2
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::MBmin,MBmax,Mremnant,SNIametals,numII,ENSN,SNIaFe,SNIaO
  real(dp)::mstarmin,mstarmax,NumSNIa
  real(dp)::masslossIa,ESNIa
  real(dp)::masslossFW,masslossW,metalsW
  real(dp)::max_massloss,m_not_lost
  real(dp)::theint,t1,t2
  real(dp)::Mwindmin,Mwindmax,mett
  real(dp)::vol_loc2
  real(dp):: NSNIa,MIMF,SNIIFe,SNIIO,SNIIej,NSNII,MIMFChabrier,IMFChabrier
  external NSNIa,MIMF,SNIIFe,SNIIO,SNIIej,NSNII,MIMFChabrier,IMFChabrier
  real(dp)::eta1,eta2,alpha1,alpha2,beta,Cr1,Cr2,mumax,eps_cl,tcl,Mclmin,Mclmax
  real(dp)::alpha,mtrans,Zscale,Zgas,tcut
  integer::iicell,iskip
  integer,dimension(1:nvector,1:2,1:2,1:2)::indcube2  !Oscar
  real(dp)::xcont,ycont,zcont,contr
  real(dp)::Lum,p10  
  real(dp)::tt,tekin,vkick,maxPrad,maxcl,maxTau,maxMet
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  real(dp),dimension(1:nvector)::mcl,Prad,agecl,Lumcl
  integer ,dimension(1:nvector)::ibinning
  integer::irad,icenter
  real(dp)::vmax,momx,momy,momz,pST,n0,imfboost
  integer::indpmax,iradmax
  integer,dimension(:),allocatable::ind
  integer ,dimension(1:nvector)::indrad
  real(dp)::momxt,momyt,momzt,momtot,momtott,pmiss,clustermass,clusterage,meanmass,meanmassM
  real(dp)::SNmin,SNmax,Tmaxfb,rsf,cellsize
  real(dp)::maxadv,vxold,vyold,vzold,vxnew,vynew,vznew,Emax
  integer::indd
  real(dp)::rhofb,scale_m,twind,numresidual
!  real(dp),dimension(1:nvector,1:5)::workarray

  if(sf_log_properties) ilun=myid+103
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=1.98892d33/scale_d/scale_l/scale_l/scale_l !Msun into grams and then internal units

!  If necessary, initialize random number generator                                           
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  dx_loc2(1:nvector)=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim   !new approach
  vol_loc2=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim


! ---------------------------- SN feedback parameters ----------------------------
  ! Type Ia supernova parameters
  MBmin=3.0d0 !Observ. parameter of binaries going SNIa
  MBmax=16.0d0 !Observ. parameter of binaries going SNIa
  Mremnant=1.4d0   !Chandrasekhar mass, nothing remains
  SNIaFe=0.63d0
  SNIaO=0.13d0
  SNIametals=SNIaFe+SNIaO  
  !Mass in Msun of enriched material (Raiteri et al. 1996) 0.63 Msun Fe + 0.13Msun O (e.g. Theilemann 1986)

  ! Scale to internal units 
  SNIametals=SNIametals*scale_m
  SNIaFe=SNIaFe*scale_m
  SNIaO=SNIaO*scale_m
  Mremnant=Mremnant*scale_m

  ! AGB wind parameters
  Mwindmin=0.5d0
  Mwindmax=8.0d0
  twind=1.0d7 !'fast winds'

  ! Massive star lifetime from Myr to code units
  t0=t_sne*1d6*(365.*24.*3600.)/scale_t

!  SNenergy=1.d51 !erg/SN
  ! Type II supernova specific energy from cgs to code units
  ESN=E_SNII/(10.*2d33)/scale_v/scale_v  !energy per 10 Msun in internal units 
  ENSN=E_SNII/scale_d/scale_l/scale_l/scale_l/scale_v/scale_v   !energy in internal units
  p10=12.0d0*3000.0d5*scale_m/scale_v !12 Msun at 3000 km/s per SNII event, calibrated to get same \dot{p} as SB99 
  ESNIa=ENSN  !Oscar: assume same energy release as SNII

  SNmin=8.0d0
  SNmax=80.0d0   !Limits for SNII events

  ! Radiation pressure
  eta1=2.0d0 !See paper  
  eta2=eta_rap
  alpha2=0.0d0   
  alpha1=0.4d0  
  beta=1.7d0  !spectrum slope
  mumax=1.0d0  
  eps_cl=0.2d0 !Cluster formation efficiency, observed  
  tcl=6.0d6 !assumed clump lifetime in years 
  tcut=3.0d6 !3.0d6 !From constant Lumonisity to powerlaw                 
  mtrans=3.0d4*scale_m  
  Cr2=2.5d0*3.08568d18/scale_l !2.5 parsec in internal units
  Cr1=Cr2/(mtrans)**0.4 
  imfboost=0.3143d0/0.224468d0 !Kroupa to Chabrier
  KappaIR_0=5.0d0*scale_d*scale_l !g-1 cm2 into internal units. Depending on dust mix, it can be as high as 30 in IR
  Mclmin=100.0*scale_m  !100 Msol, internal units

  ! Radiation pressure, bolometric
  L1=imfboost*scale_t*1.9d-7/scale_v !3.0  !Specific Lbol/c: Lbol/1d6 Msun/c in CGS converted into internal. Units are cm/s^2

!-------------------------------------------------------------------------------------

  vmax=vmaxFB
  Tmaxfb=T2max
  maxadv=maxadvfb*1d5/scale_v !into internal units

!--- Oscar: move age caluclation here to only do book-keeping if necessary?

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
        x(j,idim)=x(j,idim)/dx
        id(j,idim)=int(x(j,idim))   ! NGP at level ilevel
        igd(j,idim)=id(j,idim)/2    ! Compute parent grids
     end do
  end do
 
  ok(1:np)=.true.
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
     ok(j)=ok(j).and.igrid(j)>0    ! Check if particles are entirely in level ilevel
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  indpmax=1
  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
      else
        indp(j)= nbors_father_cells(ind_grid_part(j),kg(j))  
        vol_loc(j)=vol_loc(j)*2**ndim 
        dx_loc2(j)=dx_loc2(j)*2.0    
     end if
    ! max address 
     indpmax=max(indpmax,indp(j)) 
  end do

  allocate(ind(1:indpmax))
  ind=0.0
  indrad=0.0

! Loop over cells 
  indd=0
  do kk=1,2
     do jj=1,2
        do ii=1,2
           indd=indd+1  !counter from 1,twotondim
           iskip=ncoarse+(indd-1)*ngridmax
           do j=1,np
              ind_cell(j)=iskip+igrid(j)        !Oscar: Oct indices. do feedback on these
              indcube2(j,ii,jj,kk)=ind_cell(j)  !Oscar change back
           end do
        enddo
     enddo
  enddo

  !-------------------------------------------------------------

  ! Compute individual time steps
  do j=1,np
     dteff(j)=dtnew(levelp(ind_part(j)))   !Oscar: also new (c.f. old versions)
     if(use_proper_time)then
        dteff(j)=dteff(j)*aexp**2  !dteff is in conformal
     endif
  end do
  
  ! Reset ejected mass, metallicity, thermal energy
  mloss(:)=0.d0
  mzloss(:)=0.d0
  ethermal(:)=0.d0
  Prad(:)=0.0d0
  ptot(:)=0.0d0
  mlossFe(:)=0.0 !new
  mlossO(:)=0.0  !new
  m_not_lost=0.0 !nex (Maxime)
  
  if(cosmo) then
     ! Find neighboring expansion factors                                 
     i=1
     do while(aexp_frw(i)>aexp.and.i<n_frw)
        i=i+1
     end do
     ! Interploate time
     time_simu=t_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
          & t_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))
  endif
  
  ! --------- Prad binning---------
  maxPrad=0.0
  maxcl=0.0
  maxTau=0.0
  maxMet=0.0
  irad=0
  mcl=0.0
  ibinning=0
  agecl=0.0
  Lumcl=0.0
  ! -------------------------------
  !   workarray=0.0
  
  ! Compute feedback
  do j=1,np
     mejecta=0.0
     numII=0.0
     zloss=0.0

     if(cosmo) then
        if(use_proper_time)then    !Correct as tp is in proper, and dteff is proper. Now get to yrs
           call getAgeGyr(tp(ind_part(j)), t1)          !  End-of-dt age [Gyrs]
           call getAgeGyr(tp(ind_part(j))-dteff(j), t2) !  End-of-dt age [Gyrs]
           t1=t1*1.0d9
           t2=t2*1.0d9
        else
           ! Compute star age in years,conformal time to years
           iii=1
           do while(tau_frw(iii)>tp(ind_part(j)).and.iii<n_frw)
              iii=iii+1
           end do
           
           t2=t_frw(iii)*(tp(ind_part(j))-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &                       
                & t_frw(iii-1)*(tp(ind_part(j))-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
           t2=(time_simu-t2)/(h0*1.d5/3.08d24)/(365.*24.*3600.)                        !Units of years; particle age 
           
           t1=t_frw(iii)*(tp(ind_part(j))+dteff(j)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                & t_frw(iii-1)*(tp(ind_part(j))+dteff(j)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
           t1=(time_simu-t1)/(h0*1.d5/3.08d24)/(365.*24.*3600.)                        !Units of years; particle age - dt  
           
           birth_time=tp(ind_part(j))
           t1=max(t1,0.0d0) 
           t2=max(t2,0.0d0)            
        endif
     else
        t2=t-tp(ind_part(j))           !Age at t
        t2=t2*scale_t/365./24./60./60. !For non-cosmo units, in years
        t1=t-tp(ind_part(j))-dteff(j)  !Age at t-dt, in years
        t1=t1*scale_t/365./24./60./60.
        birth_time=tp(ind_part(j))     !internal units
        t1=max(t1,0.0d0)
     endif
        
     !--------------------------- Get dying stars ----------------------
     !mett=zp(ind_part(j))
     
     !  ---------- metals is star and host cell
     ! imetal=Fe
     ! imetal+1=O
     mett=1.06d0*zp(ind_part(j))
     Zgas=1.06d0*unew(indp(j),imetal)/max(unew(indp(j),1),smallr)/0.02d0
     if(is_oxygen) then
        mett=mett+2.09d0*zp_ox(ind_part(j))!Asplund, used for total Z assuming solar mixture
        Zgas=Zgas+2.09d0*unew(indp(j),imetal+1)/ &
             & max(unew(indp(j),1),smallr)/0.02d0 !Average, solar mix in gas (Asplund)
     endif
     Zgas=max(Zgas,0.0)
     
     if (mett<0.0004) then ! Raiteri et al. 1996 limits
        mett=0.0004
     endif
     if(mett>0.05) then
        mett=0.05
     endif
     if(t1.gt.0.0)then
        call agemass(t1,mett,mstarmax)  !Get stellar masses that exit main sequence during dt
        call agemass(t2,mett,mstarmin)  
     else
        mstarmin=120.d0
        mstarmax=120.d0
     endif
     if(mstarmin.le.SNmax.and.mstarmax.ge.SNmax) then
        mstarmax=SNmax-1d-2
     endif
     if(mstarmin.le.SNmin.and.mstarmax.ge.SNmin) then
        mstarmin=SNmin+1d-2
     endif
     if(mstarmin.le.SNmin.and.mstarmax.ge.SNmax) then !if timestep is >40 Myr
        mstarmin=SNmin+1d-2
        mstarmax=SNmax-1d-2
     endif
     meanmass=(mstarmax+mstarmin)/2.0d0
     !---------------------------
     !--------------------------- Supernovae
     !---------------------------
     if(supernovae) then
        !---------------------------
        !--------------------------- Type II events
        !---------------------------
        if(mstarmax.le.SNmax.and.mstarmin.ge.SNmin) then 
           call SNIInum(mstarmin,mstarmax,theint)  !no numerical integral
           numII=mp0(ind_part(j))*theint/scale_m 
           call ranf(localseed,RandNum) 
           numresidual=numII-int(numII)  !int <1 --> 0
           numII=int(numII)
           if(RandNum<numresidual) then
              numII=numII+1       ! Another star from residual
           endif
           
           if(numII>0.0) then
              
              Zscale=(max(Zgas,0.01d0))**(-0.2)                    !---------- Thornton et al. 1998.
              n0=max(unew(indp(j),1),smallr)*scale_nH              !---------- mH/cc   
              pST=4.0d5*numII**(0.941)*n0**(-0.1176)*Zscale        !---------- km/sST momentum Blondin et al (1998). 
              pST=pST*1.0d5*scale_m/scale_v                        !---------- Into g cm/s and then into code units
              !--------------------------------------------------------------- Check cooling radius
              rsf=30.*numII**0.29*(n0**(-0.43))*((Zgas+0.01)**(-0.18)) !------ Shell formation radius (Hopkins et al. 2013, Coiffi, et al.)
              cellsize=dx_loc2(j)*scale_l/3.08d18                  !---------- local resolution element
              
              if(cellsize.lt.rsf/Nrcool)then                       !----------- Kim & Ostriker criterion. Momentum: 1/3, energy: 1/10. 
                 momST=.false.
              else
                 momST=.true.
              endif
              
              rhofb=1.d-4 !hardcoded (..and very low)
              if(n0<rhofb)momST=.false.                            !--- Optional momentum limiter in very diffuse gas
              
              if(momST) then
                 ptot(j)=ptot(j)+pST                               !------- Post Sedov Taylor momentum
              else
                 ptot(j)=ptot(j)+p10*numII                         !------- Conservative initial bblastwave momentum
              endif
              
              ! -------  Stellar mass loss and SNII massloading
              max_massloss=mp(ind_part(j))-0.1*mp0(ind_part(j))
              mejecta=min(numII*(0.7682*meanmass**1.056)*scale_m,max_massloss)       !------- Woosley Weaver 1995, Raiteri 1996 
              m_not_lost=m_not_lost+max(numII*(0.7682*meanmass**1.056)*scale_m-max_massloss,0.0)
              mloss(j)=mloss(j)+mejecta/vol_loc(j)
              ethermal(j)=ethermal(j)+numII*ENSN/vol_loc(j)        !------- ENSN per SNII event (cooled away in low res cases)
              
              ! --- Metallicity (here only Fe and O, Raiteri 1996)
              meanmassM=min(meanmass,40.0)                         !------ yields from very massive stars assumed to be same as
              mlossFe(j)=mlossFe(j)+numII*0.375d0*exp(-17.94d0/meanmassM)*scale_m/vol_loc(j)  !------ Wosley & Heger (2007)
              mlossO(j)=mlossO(j)+numII*27.66d0*exp(-51.81d0/meanmassM)*scale_m/vol_loc(j)    !------ Wosley & Heger (2007)

              ! --- Reduce star particle mass
              mp(ind_part(j))=mp(ind_part(j))-mejecta
            !   if(mp(ind_part(j)).le.0.0)then
            !      write(*,*) "Error, mp=0 from ejecta, type II"              !------ Oscar: at mp< few 100 Msun, stochastic SN can trigger this. To be solved...
            !   endif
           endif
        endif
        
        !---------------------------
        !--------------------------- Type Ia
        !---------------------------
        if(mstarmin.gt.MBmin/2.0.and.mstarmax.lt.MBmax/2.0) then  !----- Mass limits in solarmasses (see Raiteri et al. 1996)
           call SNIanum(mstarmax,mstarmin,theint)  
           NumSNIa=theint*mp0(ind_part(j))/scale_m
           
           ! Do random sampling of type Ia SNe
           call ranf(localseed,RandNum)
           numresidual=NumSNIa-int(NumSNIa)              !----- int <1 --> 0
           NumSNIa=int(NumSNIa)
           if(RandNum<numresidual) then
              NumSNIa=NumSNIa+1      
           endif
           
           if(NumSNIa>0.0) then         
              ! -- SNIA mass loss
              max_massloss=mp(ind_part(j))-0.1*mp0(ind_part(j))
              masslossIa=min(NumSNIa*Mremnant,max_massloss)
              m_not_lost=m_not_lost+max(NumSNIa*Mremnant-max_massloss,0.0)
              mloss(j)=mloss(j)+masslossIa/vol_loc(j)
              
              ! ------ Momentum, same approch as for typeII
              if(numII<=0.)then                           !---------- if cell is already check for ST momentum in Type II event, we can save so time
                 Zscale=(max(Zgas,0.01d0))**(-0.2)        !---------- Thornton et al. 1998. Zgas computed above assuming solar mixture
                 n0=max(unew(indp(j),1),smallr)*scale_nH 
                 pST=4.0d5*NumSNIa**(0.941)*n0**(-0.1176)*Zscale 
                 pST=pST*1.0d5*scale_m/scale_v      
                 rsf=30.*NumSNIa**0.29*(n0**(-0.43))*((Zgas+0.01)**(-0.18))
                 cellsize=dx_loc2(j)*scale_l/3.08d18    
                 if(cellsize.lt.rsf/Nrcool)then 
                    momST=.false.
                 else
                    momST=.true.
                 endif
              endif
              
              if(momST) then
                 ptot(j)=ptot(j)+pST  !ST momentum
              else
                 ptot(j)=ptot(j)+p10*NumSNIa  !Conservative initial momentum
              endif
              
              ! -- SNIA thermal energy
              ethermal(j)=ethermal(j)+NumSNIa*ESNIa/vol_loc(j)  
              
              if(metal)then
                 mlossFe(j)=mlossFe(j)+NumSNIa*SNIaFe/vol_loc(j)
                 mlossO(j)=mlossO(j)+NumSNIa*SNIaO/vol_loc(j)
              endif

              ! --- Reduce star particle mass
              mp(ind_part(j))=mp(ind_part(j))-masslossIa
            !   if(mp(ind_part(j)).le.0.0)then
            !      write(*,*) "Error, mp=0 from ejecta, type Ia"               !------ Oscar: at mp< few 100 Msun, stochastic SN can trigger this. To be solved...
            !   endif
           endif
        endif
     endif
     
     !---------------------------
     !--------------------------- Winds
     !---------------------------      
     if(winds) then 
        !---------------------------
        !--------------------------- High mass stars, fast winds
        !---------------------------
        if(t2.ge.0.0.and.t1.le.twind) then
           call fm_w(t1,t2,mett,theint)                                      !------- Get fraction os stellar mass lost in winds
           max_massloss=mp(ind_part(j))-0.1*mp0(ind_part(j))
           masslossFW=min(mp0(ind_part(j))*theint,max_massloss)
           m_not_lost=m_not_lost+max(mp0(ind_part(j))*theint-max_massloss,0.0)
           mloss(j)=mloss(j)+masslossFW/vol_loc(j)  
           mlossFe(j)=mlossFe(j)+zp(ind_part(j))*masslossFW/vol_loc(j)
           if(is_oxygen) then !No proper yields, just metals from stars (to be included in the future)
              mlossO(j)=mlossO(j)+zp_ox(ind_part(j))*masslossFW/vol_loc(j)
           endif

           ! --- Reduce star particle mass
           mp(ind_part(j))=mp(ind_part(j))-masslossFW
         !   if(mp(ind_part(j)).le.0.0)then
         !      write(*,*) "Error, mp=0 from ejecta, fast winds"               !------ Oscar: at mp< few 100 Msun, stochastic SN can trigger this. To be solved...
         !   endif
           
           call p_w(t1,t2,mett,theint)                                       !------ Get specific momentum
           ptot(j)=ptot(j)+mp0(ind_part(j))*theint/scale_v                   !------ scale by vol_loc later
           
           call E_w(t1,t2,mett,theint)                                       !------ Get energy
           ethermal(j)=ethermal(j)+mp0(ind_part(j))*theint/vol_loc(j)/scale_v/scale_v
        endif
        !---------------------------
        !--------------------------- Low mass stars (Kalirai et al. 2008)
        !---------------------------
        if(mstarmin<=Mwindmax.and.mstarmax>=Mwindmin.and.mstarmin<=mstarmax) then     
           call AGBmassloss(mstarmin,mstarmax,theint)                        !------ Integration limits are from m_i to m_(i-1)
           max_massloss=mp(ind_part(j))-0.1*mp0(ind_part(j))
           masslossW=min(mp0(ind_part(j))*theint,max_massloss)
           m_not_lost=m_not_lost+max(mp0(ind_part(j))*theint-max_massloss,0.0)
           mloss(j)=mloss(j)+masslossW/vol_loc(j)
           if(metal)then
              metalsW=zp(ind_part(j))                                        !------ Same as host
              mlossFe(j)=mlossFe(j)+masslossW*metalsW/vol_loc(j)
              if(is_oxygen) then
                 metalsW=zp_ox(ind_part(j))                                  !------ Same as host
                 mlossO(j)=mlossO(j)+masslossW*metalsW/vol_loc(j)
              endif
           endif

           ! --- Reduce star particle mass
           mp(ind_part(j))=mp(ind_part(j))-masslossW
         !   if(mp(ind_part(j)).le.0.0)then
         !      write(*,*) "Error, mp=0 from ejecta, slow winds"              !------ Oscar: at mp< few 100 Msun, stochastic SN can trigger this. To be solved...
         !   endif
        endif
     endif
     !---------------------------
     !--------------------------- Radiation pressure, not to be used when full RT is used
     !---------------------------
     if(radpressure) then
        if(t2.gt.0.0.and.t1.lt.tcl) then                      ! -- Bin young star particles get "cluster mass" in cell. Use tcl~1-10Myr
           if(t2<=tcut) then                                  ! -- t2 is current age
              Lum=L1*mp0(ind_part(j))       
           else        
              Lum=(L1*(t2/tcut)**(-1.25d0))*mp0(ind_part(j)) 
           endif
           if(t2>40.0d6) then !No more
              Lum=0.0
           endif
           if(ind(indp(j)).gt.0) then                           ! -- Already a particle at location
              icenter=ind(indp(j)) 
              mcl(icenter)=mcl(icenter)+mp0(ind_part(j))
              agecl(icenter)=agecl(icenter)+t2*mp0(ind_part(j)) ! -- get average age of stars in cluster
              Lumcl(icenter)=Lumcl(icenter)+Lum*(t2-t1)*365.*24.*3600./scale_t 
           else
              irad=irad+1
              ind(indp(j))=irad !Associate entry (cell index) to Prad entry
              mcl(irad)=mcl(irad)+mp0(ind_part(j))
              agecl(irad)=agecl(irad)+t2*mp0(ind_part(j)) !to get average age of stars in cluster                               
              Lumcl(irad)=Lumcl(irad)+Lum*(t2-t1)*365.*24.*3600./scale_t
              indrad(irad)=j !We need book-keeping to get back to indcube.
           endif
        endif
     endif
  enddo
  if(m_not_lost.gt.0.0) then
     write(*,*) "Mass not lost due to 10% threshold: ", m_not_lost, " Msun"
  endif
  !---------------------------
  !--------------------------- Radiation pressure magnitude (abbove step is finding mass in young stars)
  !---------------------------
  iradmax=irad
  if(radpressure) then
     if(iradmax.gt.0) then
        do i=1,iradmax                       ! -- over cells now
           if(Lumcl(i).gt.0.0) then          ! -- only do for cells with actual young stars in them
              agecl(i)=agecl(i)/mcl(i)       ! -- Normalize cluster age
              iicell=indp(indrad(i))         ! -- get particle "j"
              
              if(metalscaling) then          ! -- for dust opacity
                 Zgas=1.06d0*unew(iicell,imetal)/ &
                      & max(unew(iicell,1),smallr)/0.02d0
                 if(is_oxygen) then
                    Zgas=Zgas+2.09d0*unew(iicell,imetal+1)/&
                         & max(unew(iicell,1),smallr)/0.02d0 !Average, solar mix in gas (Asplund)
                 endif
                 Zgas=max(Zgas,0.01) 
              else
                 Zgas=1.0
              endif
              Mclmax=mumax*mcl(i)  
              KappaIR=KappaIR_0*Zgas !Scaled by Z/Z_sun to get dust-to-gas ratio dependency. Current SN ejecta is included
              if(mcl(i)<=mtrans) then 
                 Cr=Cr1        
                 alpha=alpha1     
              else        
                 Cr=Cr2    
                 alpha=alpha2     
              endif
              if(tau_IR.ge.0) then
                 tauIR=tau_IR !set in paramter file
              else
                 tauIR=KappaIR*((1.-eps_cl)*(2.-beta)/(2.*acos(-1.)*Cr**2)/(3.-beta-2.*alpha))
                 tauIR=tauIR*(mumax/eps_cl)**(1.-2.*alpha)*(1.-(Mclmin/Mclmax)**(3.-2.*alpha-beta))
                 tauIR=tauIR/(1.-(Mclmin/Mclmax)**(2.-beta))
                 tauIR=tauIR*mcl(i)**(1.-2.*alpha)   !Correct, as we multiply by mp below
              endif
              
              tauIR=min(tauIR,50.0d0)  !Prad limiter
              
              if(agecl(i).lt.tcl)then !If clump is still intact
                 Prad(i)=(eta1+eta2*tauIR)*Lumcl(i)  !Lumcl=L1*mcl, we don't need mcl here!! dteff is accounted for above!
              else
                 Prad(i)=(eta1+eta2*KappaIR*(unew(iicell,1))*dx_loc)*Lumcl(i)  !Lumcl=L1*mcl*dteff
              endif
              
              if(ysc_stats) then  ! Get YSC stats printed to screen (change to file)
                 if(agecl(i).lt.10.0d6)then 
                    clustermass=mcl(i)*scale_d*scale_l*scale_l*scale_l/2.0d33 !Young cluster mass
                    clusterage=agecl(i)/1.0d6  !In Myr
                    write(*,"(a,6e14.5,i9)") "YSO:",t,clusterage,clustermass,unew(iicell,1)*scale_nH,tauIR,Zgas,ilevel
                 endif
              endif
           endif
        enddo
     endif
  endif
  
  !----------- Inject feedback ----------------
  
  do j=1,np
     if(mloss(j)>0.or.ethermal(j)>0.or.ptot(j)>0.) then  ! -- only enter if star actually injects something
        if(ok(j)) then                                 ! -- for particles in oct
           do ii=1,2                                      ! -- Do feedback over 3x3 cube
              do jj=1,2
                 do kk=1,2
                    iicell=indcube2(j,ii,jj,kk)
                    
                    if(iicell.gt.0) then
                       !-------------- Do momentum and mass-loss under assumption of no Etherm increase
                       tt=unew(iicell,ndim+2)  
                       tekin=0.0d0
                       do idim=1,ndim
                          tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr) ! -- Kinetic E
                       end do
                       tt=tt-tekin  !Etherm                                                                  
                       pmiss=0.0
                       
                       !------- Geomtrical factors and cell index ---------------------------------------
                       xcont=-1.0+2.0*(ii-1)
                       ycont=-1.0+2.0*(jj-1)
                       zcont=-1.0+2.0*(kk-1)
                       contr=8.0*(xcont**2+ycont**2+zcont**2)**0.5  ! -- Each cell get 1/8th of momentum
                       
                       !----- Return ejected mass -------------------------------------------------------
                       
                       if(mloss(j)>0.) then
                          unew(iicell,1)=unew(iicell,1)+(mloss(j))/8.0                     ! -- Spread over 8 cells
                          unew(iicell,2)=unew(iicell,2)+(mloss(j))*vp(ind_part(j),1)/8.0   ! -- Ejected mass momentum due to stars motion
                          unew(iicell,3)=unew(iicell,3)+(mloss(j))*vp(ind_part(j),2)/8.0
                          unew(iicell,4)=unew(iicell,4)+(mloss(j))*vp(ind_part(j),3)/8.0
                       endif
                       
                       if(momentum) then
                          if(ptot(j)>0.) then
                             !------ Momentum from winds and SNe -----------------
                             vkick=scale_v*ptot(j)/8.d0/1.d5/max(unew(iicell,1),smallr)/vol_loc(j) !Velocity for exactly ptot/8/mas\
                             
                             vxold=unew(iicell,2)/max(unew(iicell,1),smallr) !extra safety
                             vyold=unew(iicell,3)/max(unew(iicell,1),smallr)
                             vzold=unew(iicell,4)/max(unew(iicell,1),smallr)
                             
                             if(vkick.gt.vmax) then  !Limit momentum in this way (for stabbility)
                                momx=xcont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density, factor of 8 is to cancel contr. Vmax is the correct velocity
                                momy=ycont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density                       
                                momz=zcont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density                   
                                ! missing -------
                                momtot=(momx**2+momy**2+momz**2)**0.5
                                momxt=xcont*ptot(j)/contr/vol_loc(j)
                                momyt=ycont*ptot(j)/contr/vol_loc(j)
                                momzt=zcont*ptot(j)/contr/vol_loc(j)
                                momtott=(momxt**2+momyt**2+momzt**2)**0.5
                                pmiss=pmiss+(momtott-momtot)
                             else
                                momx=xcont*ptot(j)/contr/vol_loc(j)
                                momy=ycont*ptot(j)/contr/vol_loc(j)
                                momz=zcont*ptot(j)/contr/vol_loc(j)
                             endif
                             unew(iicell,2)=unew(iicell,2)+momx 
                             unew(iicell,3)=unew(iicell,3)+momy 
                             unew(iicell,4)=unew(iicell,4)+momz 
                          endif
                          ! ------------------------------
                          if(fbsafety) then                       
                             vxnew=unew(iicell,2)/max(unew(iicell,1),smallr)  !extra safety
                             vynew=unew(iicell,3)/max(unew(iicell,1),smallr)
                             vznew=unew(iicell,4)/max(unew(iicell,1),smallr)
                             
                             if(abs(vxnew).gt.maxadv) then
                                unew(iicell,2)=sign(maxadv,vxnew)*unew(iicell,1) 
                             endif
                             if(abs(vynew).gt.maxadv) then
                                unew(iicell,3)=sign(maxadv,vynew)*unew(iicell,1)
                             endif
                             if(abs(vznew).gt.maxadv) then
                                unew(iicell,4)=sign(maxadv,vznew)*unew(iicell,1)
                             endif
                          endif
                          
                          ! ----- All momentum is now added, calculate new Ekin and update Etot. This is consistent with SNe mass loading
                          tekin=0.0d0
                          do idim=1,ndim
                             tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)            ! -- Kin E 
                          enddo
                          unew(iicell,ndim+2)=tt+tekin    
                       endif
                       !tt is old Etherm, which should not change due to momentum additions
                       if(energy) then
                          if(ethermal(j)>0.) then
                             !------ Update thermal energy -----------------------------------------------------
                             if(iicell==indp(j))then  !"central cell"
                                
                                Emax=tekin+Tmaxfb*unew(iicell,1)/scale_T2/(gamma-1.0d0)
                                unew(iicell,ndim+2)=unew(iicell,ndim+2)+ethermal(j)
                                
                                !       if(pnonthermal) then  !adding missed momentum to cell via pressure
                                !          unew(iicell,ndim+2)=unew(iicell,ndim+2)+pmiss/6.0/dx_loc2(j)/dteff(j)
                                !       endif
                                
                                unew(iicell,ndim+2)=min(unew(iicell,ndim+2),Emax)  !Oscar, safety
                                
                             endif
                          endif
                       endif
                       ! --------- Finally Return metals -----
                       if(mlossFe(j)>0.or.mlossO(j)>0) then
                          unew(iicell,imetal)=unew(iicell,imetal)+mlossFe(j)/8.0   !8 cell stencil, neighbours must have same volume (oct)
                          if(is_oxygen) then !8 cell stencil, neighbours must have same volume (oct)
                             unew(iicell,imetal+1)=unew(iicell,imetal+1)+mlossO(j)/8.0
                          endif
                       endif
                    endif
                    if(iicell.le.0) then
                       write(*,*) "Error in oct, iicell=0",iicell
                    endif
                 enddo
              enddo
           enddo
        else   ! if not OK ---------------------------------------  DRIFTERS
           
           iicell=indp(j)
           
           !-------------- Do momentum and mass-loss under assumption of no Etherm increase
           tt=unew(iicell,ndim+2)  !Tot E 
           tekin=0.0d0
           do idim=1,ndim
              !         tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E                       
              tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E
           end do
           tt=tt-tekin  !Etherm                                                                  
           pmiss=0.0
           !------- Geomtrical factors and cell index ---------------------------------------------
           xcont=-1.0+2.0*(ii-1) !rand this
           ycont=-1.0+2.0*(jj-1)
           zcont=-1.0+2.0*(kk-1)
           contr=(xcont**2+ycont**2+zcont**2)**0.5  !Each cell get 1/8th of momentum
           
           !----- Return ejected mass ------------------------------------------------------------- 
           if(mloss(j)>0.) then      
              unew(iicell,1)=unew(iicell,1)+(mloss(j))  !Spread over 1 cells
              unew(iicell,2)=unew(iicell,2)+(mloss(j))*vp(ind_part(j),1)  !-------- Ejected mass momentum due to stars motion
              unew(iicell,3)=unew(iicell,3)+(mloss(j))*vp(ind_part(j),2)
              unew(iicell,4)=unew(iicell,4)+(mloss(j))*vp(ind_part(j),3)
           endif
           !------ Momentum
           
           ! Future: add momentum as random kick
           
           ! ----- All momentum is now added, calculate new Ekin and update Etot. This is consistent with SNe mass loading
           tekin=0.0d0
           do idim=1,ndim
              !           tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E                          
              tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E                          \
              
           enddo
           unew(iicell,ndim+2)=tt+tekin    !tt is old Etherm, which should not change due to momentum additions
           
           if(energy) then
              !------ Update thermal energy -----------------------------------------------------
              
              Emax=tekin+Tmaxfb*unew(iicell,1)/scale_T2/(gamma-1.0d0)
              unew(iicell,ndim+2)=unew(iicell,ndim+2)+ethermal(j)
              
              !if(pnonthermal) then !momentum done in nonthermal way
              !   unew(iicell,ndim+2)=unew(iicell,ndim+2)+ptot(j)/6.0/dx_loc2(j)/dteff(j)
              !endif
              
              unew(iicell,ndim+2)=min(unew(iicell,ndim+2),Emax)  !Oscar, safety
           endif
           
           ! --------- Finally Return metals -----
           if(mlossFe(j)>0.or.mlossO(j)>0) then
              unew(iicell,imetal)=unew(iicell,imetal)+mlossFe(j)      !1 cell for drifters
              if(is_oxygen) then
                 unew(iicell,imetal+1)=unew(iicell,imetal+1)+mlossO(j)
              endif
           endif
        endif
     endif
  enddo
  
  !-------------------------- Radiation pressure ------------------
  if(radpressure) then !Put this with the rest later
     if(iradmax.gt.0) then 
        do i=1,iradmax                   
           if(Prad(i)>0) then
              !We now the number of entries in Prad array, not particles    
              icenter=indrad(i)                          ! -- Particle index   
              if(ok(icenter)) then                       ! -- Only do for particles on local oct
                 do ii=1,2  !j is now picked
                    do jj=1,2
                       do kk=1,2
                          iicell=indcube2(icenter,ii,jj,kk) 
                          
                          ! -- Kicks are applied assuming constant etherm
                          tt=unew(iicell,ndim+2)  !Tot E                                                                           
                          tekin=0.0d0
                          do idim=1,ndim
                             tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E
                          end do
                          tt=tt-tekin  !Etherm     
                          
                          !------- Geomtrical factors and cell index ---------------------------------------------
                          xcont=-1.0+2.0*(ii-1)
                          ycont=-1.0+2.0*(jj-1)
                          zcont=-1.0+2.0*(kk-1)
                          contr=8.0*(xcont**2+ycont**2+zcont**2)**0.5  !Each cell get 1/8th of momentum
                          
                          !               vkick=scale_v*Prad(i)/8.d0/1.d5/max(unew(iicell,1),smallr)/vol_loc(icenter)
                          
                          vkick=scale_v*Prad(i)/8.d0/1.d5/max(unew(iicell,1),smallr)/vol_loc(icenter)
                          
                          if(vkick.gt.vmax) then  !Limit momentum
                             momx=xcont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density, factor of 8 is to cancel contr
                             momy=ycont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density                       
                             momz=zcont*8.*unew(iicell,1)*vmax*1.d5/scale_v/contr  !mom density                   
                          else
                             momx=xcont*Prad(i)/contr/vol_loc(icenter)
                             momy=ycont*Prad(i)/contr/vol_loc(icenter)
                             momz=zcont*Prad(i)/contr/vol_loc(icenter)
                          endif
                          unew(iicell,2)=unew(iicell,2)+momx 
                          unew(iicell,3)=unew(iicell,3)+momy 
                          unew(iicell,4)=unew(iicell,4)+momz 
                          
                          if(fbsafety) then
                             
                             vxnew=unew(iicell,2)/max(unew(iicell,1),smallr)
                             vynew=unew(iicell,3)/max(unew(iicell,1),smallr)
                             vznew=unew(iicell,4)/max(unew(iicell,1),smallr)
                             
                             if(abs(vxnew).gt.maxadv) then
                                unew(iicell,2)=sign(maxadv,vxnew)*unew(iicell,1)
                             endif
                             if(abs(vynew).gt.maxadv) then
                                unew(iicell,3)=sign(maxadv,vynew)*unew(iicell,1)
                             endif
                             if(abs(vznew).gt.maxadv) then
                                unew(iicell,4)=sign(maxadv,vznew)*unew(iicell,1)
                             endif
                          endif
                          
                          ! ----- All momentum is now added, calculate new Ekin and update Etot
                          tekin=0.0d0
                          do idim=1,ndim
                             tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr)  !Kin E     
                          enddo
                          unew(iicell,ndim+2)=tt+tekin  
                          !--------------------------------- Oscar: add pnonthemal for missing momentum
                       enddo
                    enddo
                 enddo
              else
                 iicell=indp(icenter)         
                 ! -- Kicks are applied assuming constant etherm
                 if(pnonthermal) then  !adding missed momentum to cell via pressure
                    tt=unew(iicell,ndim+2)  !Tot E
                    tekin=0.0d0
                    do idim=1,ndim
                       tekin=tekin+0.5*unew(iicell,idim+1)**2/max(unew(iicell,1),smallr) !Kin E
                    end do
                    tt=tt-tekin  !Etherm + Erad
                    
                    ! nonthermal
                    unew(iicell,ndim+2)=unew(iicell,ndim+2)+Prad(i)/6.0/dx_loc2(icenter)/dteff(icenter)
                    
                    ! -- limiter
                    Emax=tekin+Tmaxfb*unew(iicell,1)/scale_T2/(gamma-1.0d0)
                    unew(iicell,ndim+2)=min(unew(iicell,ndim+2),Emax)  !Oscar, safety 
                    
                 endif
              endif
           endif
        enddo
     endif
  endif
  
  deallocate(ind) 

end subroutine agertz_feedbk
#endif
!####################################################################
!####################################################################
!####################################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!---------------------------------------
subroutine SNIInum(m1,m2,NSNII)
!---------------------------------------
  implicit none
  REAL(kind=8),intent(out) :: NSNII
  REAL(kind=8), intent(in) :: m1,m2
  REAL(kind=8):: A,ind

 ! A=0.2244557d0  !K01
  A=0.31491d0 !Chabrier 2003
  ind=-2.3d0

  NSNII=(-A/1.3)*(m2**(-1.3)-m1**(-1.3))

END subroutine SNIInum
!###########################################################
!---------------------------------------
subroutine SNIanum(m1,m2,NSNIa)
!---------------------------------------
  implicit none
  REAL(kind=8),intent(out) :: NSNIa
  REAL(kind=8), intent(in) :: m1,m2
  REAL(kind=8):: A,Ap,N1,N2,SNIafrac

 ! A=0.2244557d0  !K01
  SNIafrac=0.16d0 !Bergh & McClure 1994, rate of SN per century in a MW type galaxy
  A=0.31491d0 !Chabrier 2003
  Ap=SNIafrac*A

  N1=(-Ap*m1**2/3.3)*((2.*m1)**(-3.3)-(m1+8.)**(-3.3))  ! (eq 14 in Agertz et al. 2013)
  N2=(-Ap*m2**2/3.3)*((2.*m2)**(-3.3)-(m2+8.)**(-3.3))

  NSNIa=(m2-m1)*(N1+N2)/2.  !Trapez. For relevant timesteps, error is epsilon

END subroutine SNIanum

!###########################################################
!---------------------------------------
subroutine AGBmassloss(m1,m2,AGB)
!---------------------------------------
  implicit none
  REAL(kind=8),intent(out) :: AGB
  REAL(kind=8), intent(in) :: m1,m2
  REAL(kind=8):: A,N1,N2

 ! A=0.2244557d0  !K01
  A=0.31491d0 !Chabrier 2003

  N1=(m1**(-0.3))*(0.3031/m1-2.97)  !Agertz et al. 2013
  N2=(m2**(-0.3))*(0.3031/m2-2.97)
  AGB=A*(N2-N1)

END subroutine AGBmassloss
!###########################################################
!------------------------------------------------
SUBROUTINE fm_w(t_1,t_2,smet,fmw)
!-----------------------------------------------
  implicit none
  real(kind=8),intent(in)::t_1,t_2,smet
  real(kind=8),intent(out)::fmw
  real(kind=8)::a,b,ts,metalscale,imfboost

  imfboost=0.3143d0/0.224468d0  !K01 to Chabrier
  
  !--- Fitting parameters --- 
  a=0.024357d0*imfboost
  b=0.000460697d0  
  ts=1.0d7
  metalscale=a*log(smet/b+1.d0)
  fmw=0.0d0

  if(t_2.le.ts) then
     fmw=metalscale*(t_2-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.lt.ts.and.t_2.gt.ts) then
     fmw=metalscale*(ts-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.ge.ts) then
     fmw=0.0
  endif

end SUBROUTINE fm_w

!###########################################################
!###########################################################
!###########################################################
!###########################################################
!------------------------------------------------
SUBROUTINE p_w(t_1,t_2,smet,momW)
!-----------------------------------------------
  implicit none
  real(kind=8),intent(in)::t_1,t_2,smet
  real(kind=8),intent(out)::momW
  real(kind=8)::a,b,c,ts,metalscale,imfboost

  imfboost=0.31430400d0/0.224468d0  !K01 to Chabrier
  !--- Fitting parameters --- 
  a=imfboost*1.8d46/1.0d6/2.d33 !Scale to per gram
  b=0.00961529d0
  c=0.363086d0
  ts=6.5d6 
  momW=0.0d0

  metalscale=a*(smet/b)**c

  if(t_2.le.ts) then
     momW=metalscale*(t_2-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.lt.ts.and.t_2.gt.ts) then
     momW=metalscale*(ts-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.ge.ts) then
     momW=0.0
  endif

end SUBROUTINE p_w

!###########################################################
!###########################################################
!###########################################################
!###########################################################
!------------------------------------------------
SUBROUTINE E_w(t_1,t_2,smet,EW)
!-----------------------------------------------
  implicit none
  real(kind=8),intent(in)::t_1,t_2,smet
  real(kind=8),intent(out)::EW
  real(kind=8)::a,b,c,ts,metalscale,imfboost

  imfboost=0.31430400d0/0.224468d0  !K01 to Chabrier

  !--- Fitting parameters --- 
  a=imfboost*1.9d54/1.0d6/2.d33 !Scale to per gram
  b=0.0101565d0
  c=0.41017d0
  ts=6.5d6 
  EW=0.0d0

  metalscale=a*(smet/b)**c

  if(t_2.le.ts) then
     EW=metalscale*(t_2-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.lt.ts.and.t_2.gt.ts) then
     EW=metalscale*(ts-t_1)/ts  !Linear mass loss, m(t2)-m(t1)
  endif

  if(t_1.ge.ts) then
     EW=0.0
  endif

end SUBROUTINE E_w

!###########################################################
!###########################################################
!###########################################################
!----------------------------------------------
SUBROUTINE agemass(time,met,mass)
!---------------------------------------
  implicit none
  real*8::a0,a1,a2,a,b,c,zzz
  real(kind=8),intent(in)::time, met
  real(kind=8),intent(out)::mass
  
  !      IMPLICIT REAL*8 (A-H,L-Z)
  !
  !     Following Raiteri et al.
  !
  !     Masses: 0.6-120.0 M_sun and Z: 7e-5 to 3e-2
  !
  if(met.lt.7.0d-5) then
     zzz=7.0d-5
  else
     zzz=met
  endif
  if(met.gt.3.0d-2) then
     zzz=3.0d-2
  else
     zzz=met
  endif
  
  a0=10.13+0.07547*log10(zzz)-0.008084*(log10(zzz))**2 
  a1=-4.424-0.7939*log10(zzz)-0.1187*(log10(zzz))**2
  a2=1.262+0.3385*log10(zzz)+0.05417*(log10(zzz))**2
  
  c=(-log10(time)+a0)
  b=a1
  a=a2
  
  if(b*b-4.*a*c.ge.0.0) then 
     mass=-b-sqrt(b*b-4.0*a*c)
     mass=mass/(2.0*a)
     mass=10.0**mass
  else
     mass=120.0
  endif
  
END SUBROUTINE agemass

!###########################################################
!###########################################################
!###########################################################
subroutine mechanical_feedback_fine(ilevel,icount)
  use pm_commons
  use amr_commons
  use mechanical_commons
#ifdef RT
  use rt_parameters,only:group_egy
  use SED_module,only:nSEDgroups,inp_SED_table
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine computes the energy liberated from supernova II ,
  ! and inject momentum and energy to the surroundings of the young stars.
  ! This routine is called every fine time step.
  ! ind_pos_cell: position of the cell within an oct
  ! m8,mz8,nph8: temporary variable necessary for an oct to add up 
  !              the mass, metal, etc. on a cell by cell basis
  ! mejecta: mass of the ejecta from SNe
  ! Zejecta: metallicity of the ejecta (not yield)
  ! mZSNe : total metal mass from SNe in each cell
  ! mchSNe: total mass of each chemical element in each cell
  ! nphSNe: total production rate of ionising radiation in each cell.
  !         This is necessary to estimate the Stromgren sphere
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::npart1,npart2,icpu,icount,idim,ip
  integer::ind,ind_son,ind_cell,ilevel,iskip,info
  integer::nSNc,nSNc_mpi
  integer,dimension(1:nvector),save::ind_grid,ind_pos_cell
  real(dp)::nsnII_star,mass0,mass_t,nsnII_tot,nsnII_mpi
  real(dp)::tyoung,current_time,dteff
  real(dp)::skip_loc(1:3),scale,dx,dx_loc,vol_loc,x0(1:3)
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_msun
  real(dp),dimension(1:twotondim),save::m8, mz8, n8, nph8 ! SNe
  real(dp),dimension(1:twotondim,1:3),save:: p8  ! SNe
  real(dp),dimension(1:nvector),save::nSNe, mSNe, mZSNe, nphSNe
  real(dp),dimension(1:nvector,1:3),save::pSNe
  real(dp),dimension(1:nvector,1:nchem),save::mchSNe
#if NCHEM>0  
  integer::ich
  real(dp),dimension(1:twotondim,1:nchem),save::mch8 ! SNe
#endif
  real(dp)::mejecta,Zejecta,mfrac_snII,M_SN_var,snII_freq
  real(dp),parameter::msun2g=2d33
  real(dp),parameter::myr2s=3.15576000d+13
  real(dp)::ttsta
  logical::ok,done_star
#ifdef RT
  real(dp),allocatable,dimension(:)::L_star
  real(dp)::Z_star,age,L_star_ion
  integer::igroup
  Z_star=z_ave
  allocate(L_star(1:nSEDgroups))
#endif
 
#if NDIM==3

  if(icount==2) return
  if(.not.hydro) return
  if(ndim.ne.3)  return
  if(numbtot(1,ilevel)==0)return
  if(nstar_tot==0)return

#ifndef WITHOUTMPI
  if(myid.eq.1) ttsta=MPI_WTIME(info)
#endif 
  nSNc=0; nsnII_tot=0d0

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g

  ! Mesh spacing in that level
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc) 

  ! To filter out old particles and compute individual time steps
  ! NB: the time step is always the coarser level time step, since no feedback for icount=2
  if (ilevel==levelmin)then
     dteff = dtnew(ilevel)
  else
     dteff = dtnew(ilevel-1)
  endif

  if (use_proper_time)then
     tyoung = t_sne*myr2s/(scale_t/aexp**2) 
     current_time=texp
     dteff = dteff*aexp**2
  else
     tyoung = t_sne*myr2s/scale_t 
     current_time=t
  endif
  tyoung = current_time - tyoung

  ncomm_SN=0  ! important to initialize; number of communications (not SNe)
#ifndef WITHOUTMPI
  xSN_comm=0d0;ploadSN_comm=0d0;mSN_comm=0d0
  mloadSN_comm=0d0;mZloadSN_comm=0d0;iSN_comm=0
  floadSN_comm=0d0;eloadSN_comm=0d0
#endif

  ! Type II Supernova frequency per Msun; irrelevant for BPASS_v2
  snII_freq = eta_sn / M_SNII

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ip=0

     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
 
        ! Count star particles
        if(npart1>0)then
           do idim=1,ndim
              x0(idim) = xg(igrid,idim) -dx -skip_loc(idim) 
           end do
 
           ipart=headp(igrid)

           m8=0d0;mz8=0d0;p8=0d0;n8=0d0;nph8=0d0
#if NCHEM>0  
           mch8=0d0
#endif
           
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
#ifdef POP3
              if(pop3 .and. (zp(ipart).lt.Zcrit_pop3) ) then
                 ok=.false. ! this will be done elsewhere
                 done_star=.false.
              else
#endif
                 nsnII_star=0d0

                 ! active star particles?
                 ok = is_star(typep(ipart)) .and. idp(ipart).lt.0

                 if(ok)then
                    ! initial mass
                    mass0 = mp0(ipart)*scale_msun
                    mass_t = mp (ipart)*scale_msun

                    ok=.false.
                    if(mechanical_bpass)then
                       call get_number_of_sn2_bpass (tp(ipart), dteff, zp(ipart), idp(ipart),&
                                & mass0, mass_t, nsnII_star, done_star)
                       if(nsnII_star>0)ok=.true.
                    else if(sn2_real_delay)then
                       if(tp(ipart).ge.tyoung)then  ! if younger than t_sne
                          call get_number_of_sn2  (tp(ipart), dteff, zp(ipart), idp(ipart),&
                                & mass0, mass_t, nsnII_star, done_star)
                          if(nsnII_star>0)ok=.true.
                       endif
                    else ! single SN event
                       if(tp(ipart).le.tyoung)then   ! if older than t_sne
                          ok=.true. 
                          ! number of sn doesn't have to be an integer
                          nsnII_star = mass0*eta_sn/M_SNII
                       endif
                    endif
                 endif
#ifdef POP3
              endif
#endif

              if(ok)then
                 ! Find the cell index to get the position of it
                 ind_son=1
                 do idim=1,ndim
                    ind = int((xp(ipart,idim)/scale - x0(idim))/dx)
                    ind_son=ind_son+ind*2**(idim-1)
                 end do 
                 iskip=ncoarse+(ind_son-1)*ngridmax
                 ind_cell=iskip+igrid
                 if(son(ind_cell)==0)then  ! leaf cell

                    !----------------------------------
                    ! For Type II explosions
                    !----------------------------------
                    M_SN_var = M_SNII
                    if(metal) then
                       if(variable_yield_SNII)then
                          call SNII_yield (zp(ipart), mfrac_snII, Zejecta, Zejecta_chem_II)
                          ! Adjust M_SNII mass not to double-count the mass loss from massive stars
                          M_SN_var = mass_loss_boost * mfrac_snII / snII_freq ! ex) 0.1 / 0.01 = 10 Msun
                       else
                          Zejecta = zp(ipart)+(1d0-zp(ipart))*yield
                       endif
                    endif

                    ! total ejecta mass in code units
                    mejecta = M_SN_var/scale_msun*nsnII_star

                    ! number of SNII
                    n8(ind_son) = n8(ind_son) + nsnII_star
                    ! mass return from SNe
                    m8 (ind_son)  = m8(ind_son) + mejecta
                    ! momentum from the original star, not the one generated by SNe
                    p8 (ind_son,1) = p8(ind_son,1) + mejecta*vp(ipart,1)
                    p8 (ind_son,2) = p8(ind_son,2) + mejecta*vp(ipart,2)
                    p8 (ind_son,3) = p8(ind_son,3) + mejecta*vp(ipart,3)
                    ! metal mass return from SNe including the newly synthesised one
                    if(metal)then
                       mz8(ind_son) = mz8(ind_son) + mejecta*Zejecta
                    endif
#if NCHEM>0  
                    do ich=1,nchem
                       mch8(ind_son,ich) = mch8(ind_son,ich) + mejecta*Zejecta_chem_II(ich)
                    end do 
#endif
                    ! subtract the mass return
                    mp(ipart)=mp(ipart)-mejecta
                  
                    ! mark if we are done with this particle
                    if(sn2_real_delay) then
                       if(done_star)then ! only if all SNe exploded
                          idp(ipart)=abs(idp(ipart)) 
                       endif 
                    else
                       idp(ipart)=abs(idp(ipart))
                    endif
                  
#ifdef RT
                    ! Enhanced momentum due to pre-processing of the ISM due to radiation
                    if(rt.and.mechanical_geen) then
                       ! Let's count the total number of ionising photons per sec 
                       call getAgeGyr(tp(ipart), age)
                       if(metal) Z_star=zp(ipart)
                       Z_star=max(Z_star,10.d-5)

                       ! compute the number of ionising photons from SED
                       call inp_SED_table(age, Z_star, 1, .false., L_star) ! L_star = [# s-1 Msun-1]
                       L_star_ion = 0d0
                       do igroup=1,nSEDgroups
                          if(group_egy(igroup).ge.13.6) L_star_ion = L_star_ion + L_star(igroup)
                       end do
                       nph8 (ind_son)=nph8(ind_son) + mass0*L_star_ion ! [# s-1]
                    endif
#endif

                 endif
              endif
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
    
           do ind=1,twotondim
              if (abs(n8(ind))>0d0)then
                 ip=ip+1
                 ind_grid(ip)=igrid
                 ind_pos_cell(ip)=ind
     
                 ! collect information 
                 nSNe(ip)=n8(ind)
                 mSNe(ip)=m8(ind)
                 mZSNe(ip)=mz8(ind)
                 pSNe(ip,1)=p8(ind,1)
                 pSNe(ip,2)=p8(ind,2)
                 pSNe(ip,3)=p8(ind,3)
                 nphSNe(ip)=nph8(ind)  ! mechanical_geen
#if NCHEM>0  
                 do ich=1,nchem
                    mchSNe(ip,ich)=mch8(ind,ich)
                 end do
#endif   
                 ! statistics
                 nSNc=nSNc+1
                 nsnII_tot = nsnII_tot + nsnII_star
 
                 if(ip==nvector)then
                    call mech_fine(ind_grid,ind_pos_cell,ip,ilevel,nSNe,mSNe,pSNe,mZSNe,nphSNe,mchSNe)
                    ip=0
                 endif 
              endif
           enddo
    
    
        end if
        igrid=next(igrid)   ! Go to next grid
     end do ! End loop over grids
    
     if (ip>0) then
        call mech_fine(ind_grid,ind_pos_cell,ip,ilevel,nSNe,mSNe,pSNe,mZSNe,nphSNe,mchSNe)
        ip=0
     endif

  end do ! End loop over cpus


#ifndef WITHOUTMPI
  nSNc_mpi=0; nsnII_mpi=0d0
  ! Deal with the stars around the boundary of each cpu (need MPI)
  call mech_fine_mpi(ilevel)
  call MPI_ALLREDUCE(nSNc,nSNc_mpi,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nsnII_tot,nsnII_mpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  nSNc = nSNc_mpi
  nsnII_tot = nSNII_mpi
!  if(myid.eq.1.and.nSNc>0.and.log_mfb) then
!     ttend=MPI_WTIME(info)
!     write(*,*) '--------------------------------------'
!     write(*,*) 'Time in mechanical_fine [sec]:', sngl(ttend-ttsta), nSNc, sngl(nsnII_tot) 
!     write(*,*) '--------------------------------------'
!  endif
#endif

#ifdef RT
   deallocate(L_star)
#endif

#endif
end subroutine mechanical_feedback_fine
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine mech_fine(ind_grid,ind_pos_cell,np,ilevel,nSN,mSN,pSN,mZSN,nphSN,mchSN)
  use amr_commons
  use pm_commons
  use hydro_commons
  use mechanical_commons
  implicit none
  integer::np,ilevel ! actually the number of cells
  integer,dimension(1:nvector)::ind_grid,ind_pos_cell
  real(dp),dimension(1:nvector)::nSN,mSN,mZSN,nphSN
  real(dp),dimension(1:nvector)::mloadSN,mZloadSN,eloadSN
  real(dp),dimension(1:nvector,1:3)::pSN,ploadSN
  real(dp),dimension(1:nvector,1:nchem)::mchSN
#if NCHEM>0  
  real(dp),dimension(1:nvector,1:nchem)::mchloadSN
  real(dp),dimension(1:nchem)::chload,z_ch
#endif
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine mechanical_feedback_fine 
  !-----------------------------------------------------------------------
  integer::i,j,nwco,nwco_here,idim,icell,igrid,ista,iend,ilevel2
  integer::ind_cell,ii
  real(dp)::d,u,v,w,e,z,eth,ekk,emag,Tk,d0,u0,v0,w0
  real(dp)::dx,dx_loc,scale,vol_loc,nH_cen,fleftSN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_msun,msun2g=2d33,scale_erg
  real(dp)::skip_loc(1:3),Tk0,ekk0,eth0,etot0,T2min
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  ! Grid based arrays
  real(dp),dimension(1:ndim,1:nvector),save::xc2
  real(dp),dimension(1:nvector,1:nSNnei), save::p_solid,ek_solid
  real(dp)::d_nei,Z_nei,Z_neisol,dm_ejecta,vol_nei
  real(dp)::mload,vload,Zload=0d0,f_esn2
  real(dp)::num_sn,nH_nei,Zdepen=1d0,f_w_cell,f_w_crit
  real(dp)::m_cen
  real(dp)::pvar(1:nvarMHD)
  ! For stars affecting across the boundary of a cpu
  integer, dimension(1:nSNnei),save::icpuSNnei
  logical,dimension(1:nvector,1:nSNnei),save ::snowplough
  real(dp),dimension(1:nvector)::rStrom ! in pc
  real(dp)::dx_loc_pc,psn_tr,chi_tr,psn_thor98,psn_geen15,fthor
  real(dp)::km2cm=1d5,M_SN_var,boost_geen_ad=0d0,p_hydro,vload_rad,f_wrt_snow
  ! chemical abundance
  ! fractional abundances ; for ionisation fraction and ref, etc
  real(dp),dimension(1:NVAR),save::fractions ! not compatible with delayed cooling
  integer::i_fractions
#if NENER>0
  integer::irad
#endif
#if NCHEM>0
  integer::ich
#endif
  
  ! starting index for passive variables except for imetal and chem
  i_fractions = imetal+nchem+1

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g
  scale_erg  = scale_l**3*scale_d*scale_v**2

  ! Mesh variables
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)
  dx_loc_pc = dx_loc*scale_l/3.08d18

  ! Record position of each cell [0.0-1.0] regardless of boxlen
  xc2=0d0
  do i=1,np
     do idim=1,ndim
        xc2(idim,i)=xg(ind_grid(i),idim)-skip_loc(idim)+xc(ind_pos_cell(i),idim)
     end do 
  end do

  !======================================================================
  ! Determine p_solid before redistributing mass 
  !   (momentum along some solid angle or cell) 
  ! - This way is desirable when two adjacent SNe explode simulataenously.
  ! - if the neighboring cell does not belong to myid, this will be done 
  !      in mech_fine_mpi
  !======================================================================
  p_solid=0d0;ek_solid=0d0;snowplough=.false.

  do i=1,np
     ind_cell = ncoarse+ind_grid(i)+(ind_pos_cell(i)-1)*ngridmax

     ! redistribute the mass/metals to the central cell
     call get_icell_from_pos (xc2(1:3,i), ilevel+1, igrid, icell,ilevel2)

     ! Sanity Check
     if((cpu_map(father(igrid)).ne.myid).or.&
        (ilevel.ne.ilevel2).or.&
        (ind_cell.ne.icell))     then 
        print *,'>>> fatal error in mech_fine'
        print *, cpu_map(father(igrid)),myid
        print *, ilevel, ilevel2
        print *, ind_cell, icell
        stop 
     endif
 
     num_sn    = nSN(i) * (E_SNII/1d51) ! doesn't have to be an integer
     M_SN_var = mSN(i)*scale_msun / num_sn
     nH_cen    = uold(icell,1)*scale_nH
     m_cen     = uold(icell,1)*vol_loc*scale_msun

     emag=0.0d0
#ifdef SOLVERmhd
     do idim=1,ndim
        emag=emag+0.125d0*(uold(icell,ndim+2+idim)+uold(icell,nvar+idim))**2
     enddo
#endif
     d   = uold(icell,1)
     u   = uold(icell,2)/d
     v   = uold(icell,3)/d
     w   = uold(icell,4)/d
     e   = uold(icell,5)
     e   = e-0.5d0*d*(u**2+v**2+w**2) - emag
#if NENER>0
     do irad=0,nener-1
        e = e - uold(icell,inener+irad) 
     end do
#endif
     Tk  = e/d*scale_T2*(gamma-1.0)
     if(Tk<0)then
        print *,'TKERR : mech fbk (pre-call): TK<0', TK,icell
        print *,'nH [H/cc]= ',d*scale_nH
        print *,'u  [km/s]= ',u*scale_v/1d5
        print *,'v  [km/s]= ',v*scale_v/1d5
        print *,'w  [km/s]= ',w*scale_v/1d5
!        stop
     endif

     ! For stability - T<0 seems to happen if sharp contact discontinuities occurs 
     ! due to mechanical feedback
     T2min=T2_star*(d*scale_nH/n_star)**(g_star-1.0)
     if(Tk<T2min)then
        e=T2min*d/scale_T2/(gamma-1.0)*1.2195
        uold(icell,5) = 0.5d0*d*(u**2+v**2+w**2) + e + emag
#if NENER>0
        do irad=0,nener-1
           uold(icell,5) = uold(icell,5) + uold(icell,inener+irad)
        end do
#endif
     endif
    
     if(metal)then
        z   = uold(icell,imetal)/d
     else
        z   = z_ave*0.02
     endif

     !==========================================
     ! estimate Stromgren sphere (relevant to RHD simulations only)
     ! (mechanical_geen=.true)
     !==========================================
     if(mechanical_geen.and.rt) rStrom(i) = (3d0*nphSN(i)/4./3.141592/2.6d-13/nH_cen**2d0)**(1d0/3d0)/3.08d18 ! [pc] 

     if(log_mfb)then
398     format('MFB z= ',f9.5,' N= ',f7.1,' nH,T,Z= ',3(f7.3,1x),' dx= ',f9.5)
        write(*,398) 1./aexp-1,num_sn,log10(d*scale_nH),log10(Tk),log10(z/0.02),log10(dx_loc*scale_l/3.08d18)
     endif

     dm_ejecta = f_LOAD*mSN(i)/dble(nSNnei)  ! per solid angle
     mload     = f_LOAD*mSN(i) + uold(icell,1)*vol_loc*f_LOAD  ! total SN ejecta + host cell
     if(metal) Zload = (f_LOAD*mZSN(i) + uold(icell,imetal)*vol_loc*f_LOAD)/mload
#if NCHEM>0
     do ich=1,nchem
        chload(ich) = (f_LOAD*mchSN(i,ich) + uold(icell,ichem+ich-1)*vol_loc*f_LOAD)/mload
     end do
#endif
 
     do j=1,nSNnei
        call get_icell_from_pos (xc2(1:3,i)+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell,ilevel2)
        if(cpu_map(father(igrid)).eq.myid) then ! if belong to myid

           Z_nei = z_ave*0.02 ! For metal=.false. 
           if(ilevel>ilevel2)then ! touching level-1 cells
              d_nei     = unew(icell,1)
              if(metal) Z_nei = unew(icell,imetal)/d_nei
           else
              d_nei     = uold(icell,1)
              if(metal) Z_nei = uold(icell,imetal)/d_nei
           endif

           f_w_cell  = (mload/dble(nSNnei) + d_nei*vol_loc/8d0)/dm_ejecta - 1d0
           nH_nei    = d_nei*scale_nH
           Z_neisol  = max(0.01, Z_nei/0.02)
           Zdepen    = Z_neisol**(expZ_SN*2d0) 

           ! transition mass loading factor (momentum conserving phase)
           ! psn_tr = sqrt(2*chi_tr*Nsn*Esn*Msn*fe)
           ! chi_tr = (1+f_w_crit)
           psn_thor98   = A_SN * num_sn**(expE_SN) * nH_nei**(expN_SN) * Z_neisol**(expZ_SN)  !km/s Msun
           psn_tr       = psn_thor98 
           if(mechanical_geen)then 
              ! For snowplough phase, psn_tr will do the job
              psn_geen15 = A_SN_Geen * num_sn**(expE_SN)* Z_neisol**(expZ_SN)  !km/s Msun

              if(rt)then
                 fthor   = exp(-dx_loc_pc/rStrom(i))
                 psn_tr  = psn_thor98*fthor + psn_geen15*(1d0-fthor)
              else
                 psn_tr  = psn_geen15
              endif
              psn_tr = max(psn_tr, psn_thor98)

              ! For adiabatic phase
              ! psn_tr =  A_SN * (E51 * boost_geen)**expE_SN_boost * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
              !        =  p_hydro * boost_geen**expE_SN_boost
              p_hydro = A_SN * num_sn**(expE_SN_boost) * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
              boost_geen_ad = (psn_tr / p_hydro)**(1d0/expE_SN_boost)
              boost_geen_ad = max(boost_geen_ad-1d0,0.0)
           endif

           chi_tr   = psn_tr**2d0 / (2d0 * num_sn**2d0 * (E_SNII/msun2g/km2cm**2d0) * M_SN_var * f_ESN)
           !          (Msun*km/s)^2                 (Msun *km2/s2)              (Msun)
           f_w_crit = max(chi_tr-1d0, 0d0)

           !f_w_crit = (A_SN/1d4)**2d0/(f_ESN*M_SNII)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
           !f_w_crit = max(0d0,f_w_crit)
           vload_rad = dsqrt(2d0*f_ESN*E_SNII*(1d0+f_w_crit)/(M_SN_var*msun2g))/scale_v/(1d0+f_w_cell)/f_LOAD/f_PCAN
           if(f_w_cell.ge.f_w_crit)then ! radiative phase
              ! ptot = sqrt(2*chi_tr*Mejtot*(fe*Esntot))
              ! vload = ptot/(chi*Mejtot) = sqrt(2*chi_tr*fe*Esntot/Mejtot)/chi = sqrt(2*chi_tr*fe*Esn/Mej)/chi
              vload = vload_rad
              snowplough(i,j)=.true.
           else ! adiabatic phase
              ! ptot = sqrt(2*chi*Mejtot*(fe*Esntot))
              ! vload = ptot/(chi*Mejtot) = sqrt(2*fe*Esntot/chi/Mejtot) = sqrt(2*fe*Esn/chi/Mej)
              f_esn2 = 1d0-(1d0-f_ESN)*f_w_cell/f_w_crit ! to smoothly correct the adibatic to the radiative phase
              vload = dsqrt(2d0*f_esn2*E_SNII/(1d0+f_w_cell)/(M_SN_var*msun2g))/scale_v/f_LOAD
              if(mechanical_geen) then
                 !f_wrt_snow = (f_esn2 - f_ESN)/(1d0-f_ESN)
                 f_wrt_snow = 2d0-2d0/(1d0+exp(-f_w_cell/f_w_crit/0.3)) ! 0.3 is obtained by calibrating
                 vload = vload * dsqrt(1d0 + boost_geen_ad*f_wrt_snow)
                 ! NB. this sometimes give too much momentum because expE_SN_boost != expE_SN. A limiter is needed
              endif
              ! safety device: limit the maximum velocity so that it does not exceed p_{SN,final}
              if(vload>vload_rad) vload = vload_rad
              snowplough(i,j)=.false.
           endif
           p_solid(i,j)=(1d0+f_w_cell)*dm_ejecta*vload
           ek_solid(i,j)=ek_solid(i,j)+p_solid(i,j)*(vload*f_LOAD)/2d0 !ek=(m*v)*v/2, not (d*v)*v/2

           if(log_mfb_mega)then
             write(*,'(" MFBN nHcen=", f6.2," nHnei=", f6.2, " mej=", f6.2, " mcen=", f6.2, " vload=", f6.2, " lv2=",I3," fwcrit=", f6.2, " fwcell=", f6.2, " psol=",f6.2, " mload/48=",f6.2," mnei/8=",f6.2," mej/48=",f6.2)') &
            & log10(nH_cen),log10(nH_nei),log10(mSN(i)*scale_msun),log10(m_cen),log10(vload*scale_v/1d5),ilevel2-ilevel,&
            & log10(f_w_crit),log10(f_w_cell),log10(p_solid(i,j)*scale_msun*scale_v/1d5),log10(mload*scale_msun/48),&
            & log10(d_nei*vol_loc/8d0*scale_msun),log10(dm_ejecta*scale_msun)
            endif

        endif
        
     enddo ! loop over neighboring cells
  enddo ! loop over SN cells

  !-----------------------------------------
  ! Redistribute mass from the SN cell
  !-----------------------------------------
  do i=1,np
     icell = ncoarse+ind_grid(i)+(ind_pos_cell(i)-1)*ngridmax
     d     = uold(icell,1)
     u     = uold(icell,2)/d
     v     = uold(icell,3)/d
     w     = uold(icell,4)/d
     e     = uold(icell,5)
#if NENER>0
     do irad=0,nener-1
        e = e - uold(icell,inener+irad)
     enddo
#endif
     emag=0.0d0
#ifdef SOLVERmhd
     do idim=1,ndim
        emag=emag+0.125d0*(uold(icell,ndim+2+idim)+uold(icell,nvar+idim))**2
     enddo
#endif

     ekk   = 0.5*d*(u**2 + v**2 + w**2)
     eth   = e - ekk - emag  ! thermal pressure 

     ! ionisation fractions, ref, etc.
     do ii=i_fractions,nvar
        fractions(ii) = uold(icell,ii)/d
     end do

     mloadSN (i) = mSN (i)*f_LOAD + d*vol_loc*f_LOAD
     if(metal)then
        z = uold(icell,imetal)/d
        mZloadSN(i) = mZSN(i)*f_LOAD + d*z*vol_loc*f_LOAD
     endif

#if NCHEM>0
     do ich=1,nchem
        z_ch(ich) = uold(icell,ichem+ich-1)/d
        mchloadSN(i,ich) = mchSN(i,ich)*f_LOAD + d*z_ch(ich)*vol_loc*f_LOAD
     end do 
#endif
     ! original momentum by star + gas entrained from the SN cell
     ploadSN(i,1) = pSN(i,1)*f_LOAD + vol_loc*d*u*f_LOAD
     ploadSN(i,2) = pSN(i,2)*f_LOAD + vol_loc*d*v*f_LOAD
     ploadSN(i,3) = pSN(i,3)*f_LOAD + vol_loc*d*w*f_LOAD

     ! update the hydro variable
     fleftSN = 1d0 - f_LOAD
     uold(icell,1) = uold(icell,1)*fleftSN + mSN(i)  /vol_loc*f_LEFT 
     uold(icell,2) = uold(icell,2)*fleftSN + pSN(i,1)/vol_loc*f_LEFT  ! rho*v, not v
     uold(icell,3) = uold(icell,3)*fleftSN + pSN(i,2)/vol_loc*f_LEFT 
     uold(icell,4) = uold(icell,4)*fleftSN + pSN(i,3)/vol_loc*f_LEFT 
     if(metal)then
        uold(icell,imetal) = mZSN(i)/vol_loc*f_LEFT + d*z*fleftSN
     endif
#if NCHEM>0
     do ich=1,nchem
        uold(icell,ichem+ich-1) = mchSN(i,ich)/vol_loc*f_LEFT + d*z_ch(ich)*fleftSN
     end do
#endif
     do ii=i_fractions,nvar
        uold(icell,ii) = fractions(ii) * uold(icell,1)
     end do

     ! original kinetic energy of the gas entrained
     eloadSN(i) = ekk*vol_loc*f_LOAD 

     ! original thermal energy of the gas entrained (including the non-thermal part)
     eloadSN(i) = eloadSN(i) + eth*vol_loc*f_LOAD

     ! reduce total energy as we are distributing it to the neighbours
     !uold(icell,5) = uold(icell,5)*fleftSN 
     uold(icell,5) = uold(icell,5) - (ekk+eth)*f_LOAD

     ! add the contribution from the original kinetic energy of SN particle
     d = mSN(i)/vol_loc
     u = pSN(i,1)/mSN(i)
     v = pSN(i,2)/mSN(i)
     w = pSN(i,3)/mSN(i)
     uold(icell,5) = uold(icell,5) + 0.5d0*d*(u**2 + v**2 + w**2)*f_LEFT
#if NENER>0
     ! Adding cosmic ray energy!!!
     uold(icell,inener) = uold(icell,inener) + fecr*nSN(i)*E_SNII/scale_erg/vol_loc*f_LEFT
     uold(icell,5     ) = uold(icell,5     ) + fecr*nSN(i)*E_SNII/scale_erg/vol_loc*f_LEFT
#endif
    
     ! add the contribution from the original kinetic energy of SN to outflow
     eloadSN(i) = eloadSN(i) + 0.5d0*mSN(i)*(u**2 + v**2 + w**2)*f_LOAD

     ! update ek_solid     
     ek_solid(i,:) = ek_solid(i,:) + eloadSN(i)/dble(nSNnei)

  enddo  ! loop over SN cell


  !-------------------------------------------------------------
  ! Find and save stars affecting across the boundary of a cpu
  !-------------------------------------------------------------
  do i=1,np

     ind_cell = ncoarse+ind_grid(i)+(ind_pos_cell(i)-1)*ngridmax

     nwco=0; icpuSNnei=0
     do j=1,nSNnei
        call get_icell_from_pos (xc2(1:3,i)+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell,ilevel2)
 
        if(cpu_map(father(igrid)).ne.myid) then ! need mpi
           nwco=nwco+1
           icpuSNnei(nwco)=cpu_map(father(igrid))
        else  ! can be handled locally
           vol_nei = vol_loc*(2d0**ndim)**(ilevel-ilevel2)
           pvar(:) = 0d0 ! temporary primitive variable
           if(ilevel>ilevel2)then ! touching level-1 cells
              pvar(1:nvarMHD) = unew(icell,1:nvarMHD)
           else
              pvar(1:nvarMHD) = uold(icell,1:nvarMHD)
           endif 
           do ii=i_fractions,nvar ! fractional quantities that we don't want to change
              fractions(ii) = pvar(ii)/pvar(1)
           end do

           emag=0.0d0
#ifdef SOLVERmhd
           do idim=1,ndim
              emag=emag+0.125d0*(pvar(ndim+2+idim)+pvar(nvar+idim))**2
           enddo
#endif
           d0=pvar(1)
           u0=pvar(2)/d0
           v0=pvar(3)/d0
           w0=pvar(4)/d0
           ekk0=0.5d0*d0*(u0**2+v0**2+w0**2)
           eth0=pvar(5)-ekk0-emag
#if NENER>0
           do irad=0,nener-1
              eth0=eth0-pvar(inener+irad)
           end do
#endif
           ! For stability
           Tk0 =eth0/d0*scale_T2*(gamma-1.0)
           T2min=T2_star*(d0*scale_nH/n_star)**(g_star-1.0)
           if(Tk0<T2min)then
              eth0=T2min*d0/scale_T2/(gamma-1.0)*1.2195 !Joki: what's the magic number?
           endif 

           d= mloadSN(i  )/dble(nSNnei)/vol_nei
           u=(ploadSN(i,1)/dble(nSNnei)+p_solid(i,j)*vSNnei(1,j))/vol_nei/d
           v=(ploadSN(i,2)/dble(nSNnei)+p_solid(i,j)*vSNnei(2,j))/vol_nei/d
           w=(ploadSN(i,3)/dble(nSNnei)+p_solid(i,j)*vSNnei(3,j))/vol_nei/d
           pvar(1)=pvar(1)+d
           pvar(2)=pvar(2)+d*u
           pvar(3)=pvar(3)+d*v
           pvar(4)=pvar(4)+d*w
           !ekk_ej = 0.5*d*(u**2 + v**2 + w**2)   
           etot0  = eth0+ekk0+emag+ek_solid(i,j)/vol_nei ! additional energy from SNe+entrained gas

           ! the minimum thermal energy input floor
           d   = pvar(1)
           u   = pvar(2)/d
           v   = pvar(3)/d
           w   = pvar(4)/d
           ekk = 0.5*d*(u**2 + v**2 + w**2)
           pvar(5) = max(etot0, ekk+eth0+emag)

           ! sanity check
           Tk = (pvar(5)-ekk-emag)/d*scale_T2*(gamma-1)
           if(Tk<0)then
              print *,'TKERR: mech (post-call): Tk<0 =',Tk
              print *,'nH [H/cc]= ',d*scale_nH
              print *,'u  [km/s]= ',u*scale_v/1d5
              print *,'v  [km/s]= ',v*scale_v/1d5
              print *,'w  [km/s]= ',w*scale_v/1d5
              print *,'T0 [K]   = ',Tk0
              stop
           endif 

#if NENER>0
           ! Inject cosmic rays!!
           do irad=0,nener-1
              pvar(inener+irad) = pvar(inener+irad)+fecr*nSN(i)*E_SNII/scale_erg/dble(nSNnei)/vol_nei
              pvar(5) = pvar(5) + pvar(inener+irad)
           end do
#endif
           if(metal)then
               pvar(imetal)=pvar(imetal)+mzloadSN(i)/dble(nSNnei)/vol_nei
           end if
#if NCHEM>0
           do ich=1,nchem
              pvar(ichem+ich-1)=pvar(ichem+ich-1)+mchloadSN(i,ich)/dble(nSNnei)/vol_nei
           end do
#endif
           do ii=i_fractions,nvar
               pvar(ii)=fractions(ii)*pvar(1)
           end do

           ! update the hydro variable
           if(ilevel>ilevel2)then ! touching level-1 cells
              unew(icell,1:nvarMHD) = pvar(1:nvarMHD)
           else
              uold(icell,1:nvarMHD) = pvar(1:nvarMHD)
           endif 

        end if
     end do ! loop over 48 neighbors


#ifndef WITHOUTMPI 
     if(nwco>0)then  ! for SNs across different cpu
        if(nwco>1)then
           nwco_here=nwco
           ! remove redundant cpu list for this SN cell
           call redundant_non_1d(icpuSNnei(1:nwco), nwco_here, nwco)
        endif
        ista=ncomm_SN+1
        iend=ista+nwco-1

        if(iend>ncomm_max)then
           write(*,*) 'Error: increase ncomm_max in mechanical_fine.f90', ncomm_max, iend
           call clean_stop
        endif
        iSN_comm (ista:iend)=icpuSNnei(1:nwco)
        nSN_comm (ista:iend )=nSN(i)
        mSN_comm (ista:iend )=mSN(i)
        mloadSN_comm (ista:iend  )=mloadSN(i)
        xSN_comm (1,ista:iend)=xc2(1,i)
        xSN_comm (2,ista:iend)=xc2(2,i)
        xSN_comm (3,ista:iend)=xc2(3,i)
        ploadSN_comm (1,ista:iend)=ploadSN(i,1)
        ploadSN_comm (2,ista:iend)=ploadSN(i,2)
        ploadSN_comm (3,ista:iend)=ploadSN(i,3)
        floadSN_comm (ista:iend  )=f_LOAD
        eloadSN_comm (ista:iend  )=eloadSN(i)
        if(metal) mZloadSN_comm(ista:iend  )=mZloadSN(i)
#if NCHEM>0
        do ich=1,nchem
            mchloadSN_comm(ista:iend,ich)=mchloadSN(i,ich)
         end do
#endif
        if(mechanical_geen.and.rt) rSt_comm (ista:iend)=rStrom(i)
        ncomm_SN=ncomm_SN+nwco
     endif
#endif

  end do ! loop over SN cell

end subroutine mech_fine
#endif
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine mech_fine_mpi(ilevel)
  use amr_commons
  use mechanical_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ilevel
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::i,j,info,nSN_tot,icpu,ncpu_send,ncpu_recv,ncc
  integer::ncell_recv,ncell_send,cpu2send,cpu2recv,tag,np
  integer::isend_sta,irecv_sta,irecv_end
  real(dp),dimension(:,:),allocatable::SNsend,SNrecv,p_solid,ek_solid
  integer ,dimension(:),allocatable::list2recv,list2send
  integer, dimension(:),allocatable::reqrecv,reqsend
  integer, dimension(:,:),allocatable::statrecv,statsend
  ! SN variables
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::mloadSN_i,zloadSN_i,ploadSN_i(1:3),mSN_i,xSN_i(1:3),fload_i
  real(dp)::f_esn2,d_nei,Z_nei,Z_neisol,f_w_cell,f_w_crit,nH_nei
  real(dp)::num_sn,vload,Zdepen=1d0,vol_nei,dm_ejecta
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_msun,msun2g=2d33,scale_erg, pvar(1:nvarMHD),etot0
  real(dp)::skip_loc(1:3),d,u,v,w,ekk,emag,d0,u0,v0,w0,eth0,ekk0,Tk0,T2min!,ekk_ej
  integer::igrid,icell,ilevel2,ii
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  logical,allocatable,dimension(:,:)::snowplough 
#if NCHEM>0
  real(dp),dimension(1:nchem),save::chloadSN_i
  integer::ich
#endif
#if NENER>0
  integer::irad
#endif
  real(dp)::rSt_i,dx_loc_pc,psn_tr,chi_tr,psn_thor98,psn_geen15,fthor
  real(dp)::km2cm=1d5,M_SN_var,boost_geen_ad=0d0,p_hydro,vload_rad,f_wrt_snow
  ! fractional abundances ; for ionisation fraction and ref, etc
  real(dp),dimension(1:NVAR),save::fractions ! not compatible with delayed cooling
  integer::i_fractions
#ifdef SOLVERmhd
  integer::idim
#endif
  if(ndim.ne.3) return

  ! starting index for passive variables except for imetal and chem
  i_fractions = imetal+nchem+1

  !============================================================
  ! For MPI communication
  !============================================================
  ncpu_send=0;ncpu_recv=0

  ncomm_SN_cpu=0 
  ncomm_SN_mpi=0
  ncomm_SN_mpi(myid)=ncomm_SN
  ! compute the total number of communications needed
  call MPI_ALLREDUCE(ncomm_SN_mpi,ncomm_SN_cpu,ncpu,&
                   & MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_tot = sum(ncomm_SN_cpu)
  if(nSN_tot==0) return

  allocate(icpuSN_comm    (1:nSN_tot,1:2))
  allocate(icpuSN_comm_mpi(1:nSN_tot,1:2))

  ! index for mpi variable
  if(myid==1)then
     isend_sta = 0 
  else
     isend_sta = sum(ncomm_SN_cpu(1:myid-1)) 
  endif

  icpuSN_comm=0
  do i=1,ncomm_SN_cpu(myid)
     icpuSN_comm(isend_sta+i,1)=myid
     icpuSN_comm(isend_sta+i,2)=iSN_comm(i)
     ! iSN_comm:   local variable
     ! icpuSN_comm:  local (but extended) variable to be passed to a mpi variable
     ! icpuSN_comm_mpi: mpi variable
  end do

  ! share the list of communications
  icpuSN_comm_mpi=0
  call MPI_ALLREDUCE(icpuSN_comm,icpuSN_comm_mpi,nSN_tot*2,&
                   & MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

  ncell_send = ncomm_SN_cpu(myid)
  ncell_recv = count(icpuSN_comm_mpi(:,2).eq.myid, 1)

  ! check if myid needs to send anything
  if(ncell_send>0)then
     allocate( SNsend (1:nvarSN,1:ncell_send)) ! x(3),m,mz,p,pr
     allocate( list2send(1:ncell_send)    )
     list2send=0;  SNsend=0d0
     list2send=icpuSN_comm_mpi(isend_sta+1:isend_sta+ncell_send,2)
     ncpu_send=1
     if(ncell_send>1) call redundant_non_1d (list2send,ncell_send,ncpu_send)
     ! ncpu_send = No. of cpus to which myid should send info 
     allocate( reqsend  (1:ncpu_send) )
     allocate( statsend (1:MPI_STATUS_SIZE,1:ncpu_send) )
     reqsend=0; statsend=0
  endif

  ! check if myid needs to receive anything
  if(ncell_recv>0)then
     allocate( SNrecv (1:nvarSN,1:ncell_recv)) ! x(3),m,mz,p,pr
     allocate( list2recv(1:ncell_recv)    )
     list2recv=0;  SNrecv=0d0
     j=0
     do i=1,nSN_tot
        if(icpuSN_comm_mpi(i,2).eq.myid)then
           j=j+1 
           list2recv(j) = icpuSN_comm_mpi(i,1)
        endif
     end do
 
     ncc = j
     if(ncc.ne.ncell_recv)then ! sanity check
        write(*,*) 'Error in mech_fine_mpi: ncc != ncell_recv',ncc,ncell_recv,myid
        call clean_stop
     endif

     ncpu_recv=1 ! No. of cpus from which myid should receive info 
     if(j>1) call redundant_non_1d(list2recv,ncc,ncpu_recv)

     allocate( reqrecv  (1:ncpu_recv) )
     allocate( statrecv (1:MPI_STATUS_SIZE,1:ncpu_recv) )
     reqrecv=0; statrecv=0
  endif

  ! prepare one variable and send
  if(ncell_send>0)then
     do icpu=1,ncpu_send
        cpu2send = list2send(icpu)
        ncc=0 ! number of SN host cells that need communications with myid=cpu2send
        do i=1,ncell_send
           j=i+isend_sta
           if(icpuSN_comm_mpi(j,2).eq.cpu2send)then
              ncc=ncc+1
              SNsend(1:3,ncc)=xSN_comm (1:3,i)
              SNsend(4  ,ncc)=mSN_comm (    i)
              SNsend(5  ,ncc)=mloadSN_comm (i)
              SNsend(6:8,ncc)=ploadSN_comm (1:3,i)
              SNsend(9  ,ncc)=floadSN_comm (i)
              SNsend(10 ,ncc)=eloadSN_comm (i)
              SNsend(11 ,ncc)=nSN_comm     (i)
              if(metal)SNsend(12,ncc)=mZloadSN_comm(i)
              if(mechanical_geen.and.rt)SNsend(13,ncc)=rSt_comm(i)
#if NCHEM>0
              do ich=1,nchem
                 SNsend(13+ich,ncc)=mchloadSN_comm(i,ich)
              end do
#endif
           endif
        end do ! i

        tag = myid + cpu2send + ncc
        call MPI_ISEND (SNsend(1:nvarSN,1:ncc),ncc*nvarSN,MPI_DOUBLE_PRECISION, &
                      & cpu2send-1,tag,MPI_COMM_WORLD,reqsend(icpu),info) 
     end do ! icpu

  endif ! ncell_send>0


  ! receive one large variable
  if(ncell_recv>0)then
     irecv_sta=1
     do icpu=1,ncpu_recv
        cpu2recv = list2recv(icpu)
        ncc=0 ! number of SN host cells that need communications with cpu2recv
        do i=1,nSN_tot
           if((icpuSN_comm_mpi(i,1)==cpu2recv).and.&
             &(icpuSN_comm_mpi(i,2)==myid)      )then
              ncc=ncc+1
           endif
        end do
        irecv_end=irecv_sta+ncc-1
        tag = myid + cpu2recv + ncc
        
        call MPI_IRECV (SNrecv(1:nvarSN,irecv_sta:irecv_end),ncc*nvarSN,MPI_DOUBLE_PRECISION,&
                     & cpu2recv-1,tag,MPI_COMM_WORLD,reqrecv(icpu),info)

        irecv_sta=irecv_end+1
     end do ! icpu 

  endif ! ncell_recv >0

  if(ncpu_send>0)call MPI_WAITALL(ncpu_send,reqsend,statsend,info)
  if(ncpu_recv>0)call MPI_WAITALL(ncpu_recv,reqrecv,statrecv,info)


  !============================================================
  ! inject mass/metal/momentum
  !============================================================

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g
  scale_erg  = scale_l**3*scale_d*scale_v**2
  
  ! Mesh variables
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)
  dx_loc_pc = dx_loc*scale_l/3.08d18

  np = ncell_recv
  if(ncell_recv>0) then
     allocate(p_solid(1:np,1:nSNnei))
     allocate(ek_solid(1:np,1:nSNnei))
     allocate(snowplough(1:np,1:nSNnei))
     p_solid=0d0;ek_solid=0d0;snowplough=.false.
  endif

  ! Compute the momentum first before redistributing mass
  do i=1,np

     xSN_i(1:3) = SNrecv(1:3,i)
     mSN_i      = SNrecv(4,i)
     mloadSN_i  = SNrecv(5,i)
     fload_i    = SNrecv(9,i)
     dm_ejecta  = f_LOAD*mSN_i/dble(nSNnei)
     num_sn     = SNrecv(11,i) * (E_SNII/1d51) ! doesn't have to be an integer
     ek_solid(i,:) = SNrecv(10,i)/dble(nSNnei) ! kinetic energy of the gas mass entrained from the host cell + SN
     if(metal) ZloadSN_i = SNrecv(12,i)/SNrecv(5,i)
     if(mechanical_geen.and.rt) rSt_i = SNrecv(13,i) ! Stromgren sphere in pc 

     M_SN_var = mSN_i*scale_msun/num_sn

     do j=1,nSNnei
        call get_icell_from_pos (xSN_i+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell,ilevel2)
        if(cpu_map(father(igrid)).eq.myid) then ! if belong to myid
           Z_nei = z_ave*0.02 ! for metal=.false.
           if(ilevel>ilevel2)then ! touching level-1 cells
              d_nei     = unew(icell,1)
              if(metal) Z_nei = unew(icell,imetal)/d_nei
           else
              d_nei     = uold(icell,1)
              if(metal) Z_nei = uold(icell,imetal)/d_nei
           endif
           f_w_cell = (mloadSN_i/dble(nSNnei) + d_nei*vol_loc/8d0)/dm_ejecta - 1d0
           nH_nei   = d_nei*scale_nH
           Z_neisol = max(0.01,Z_nei/0.02)
           Zdepen   = Z_neisol**(expZ_SN*2d0) !From Thornton+(98)

           ! transition mass loading factor (momentum conserving phase)
           ! psn_tr = sqrt(2*chi_tr*Nsn*Esn*Msn*fe)
           ! chi_tr = (1+f_w_crit)
           psn_thor98   = A_SN * num_sn**(expE_SN) * nH_nei**(expN_SN) * Z_neisol**(expZ_SN)  !km/s Msun
           psn_tr       = psn_thor98 
           if(mechanical_geen)then
              ! For snowplough phase, psn_tr will do the job
              psn_geen15= A_SN_Geen * num_sn**(expE_SN) * Z_neisol**(expZ_SN)  !km/s Msun
              if(rt)then
                 fthor   = exp(-dx_loc_pc/rSt_i)
                 psn_tr  = psn_thor98*fthor + psn_geen15*(1d0-fthor)
              else
                 psn_tr  = psn_geen15
              endif
              psn_tr = max(psn_tr, psn_thor98)

              ! For adiabatic phase
              ! psn_tr =  A_SN * (E51 * boost_geen)**expE_SN_boost * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
              !        =  p_hydro * boost_geen**expE_SN_boost
              p_hydro = A_SN * num_sn**(expE_SN_boost) * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
              boost_geen_ad = (psn_tr / p_hydro)**(1d0/expE_SN_boost)
              boost_geen_ad = max(boost_geen_ad-1d0,0.0)
           endif
           chi_tr   = psn_tr**2d0 / (2d0 * num_sn**2d0 * (E_SNII/msun2g/km2cm**2d0) * M_SN_var * f_ESN)
           !          (Msun*km/s)^2                 (Msun *km2/s2)              (Msun)
           f_w_crit = max(chi_tr-1d0, 0d0)
 
           !f_w_crit = (A_SN/1d4)**2d0/(f_ESN*M_SNII)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
           !f_w_crit = max(0d0,f_w_crit)

           vload_rad = dsqrt(2d0*f_ESN*E_SNII*(1d0+f_w_crit)/(M_SN_var*msun2g))/(1d0+f_w_cell)/scale_v/f_LOAD/f_PCAN
           if(f_w_cell.ge.f_w_crit)then ! radiative phase
              vload = vload_rad
              snowplough(i,j)=.true.
           else ! adiabatic phase
              f_esn2 = 1d0-(1d0-f_ESN)*f_w_cell/f_w_crit
              vload = dsqrt(2d0*f_esn2*E_SNII/(1d0+f_w_cell)/(M_SN_var*msun2g))/scale_v/f_LOAD
              if(mechanical_geen) then
                 f_wrt_snow = 2d0-2d0/(1d0+exp(-f_w_cell/f_w_crit/0.30))
                 boost_geen_ad = boost_geen_ad * f_wrt_snow
                 vload = vload * dsqrt(1d0+boost_geen_ad) ! f_boost
                 ! NB. this sometimes give too much momentum because expE_SN_boost != expE_SN. A limiter is needed
              endif
              ! safety device: limit the maximum velocity so that it does not exceed p_{SN,final}
              if(vload>vload_rad) vload = vload_rad
              snowplough(i,j)=.false.
           endif
           p_solid(i,j)=(1d0+f_w_cell)*dm_ejecta*vload
           ek_solid(i,j)=ek_solid(i,j)+p_solid(i,j)*(vload*f_LOAD)/2d0 !ek
        endif
       
     enddo ! loop over neighboring cells

  end do



  ! Apply the SNe to this cpu domain
  do i=1,np

     xSN_i(1:3) = SNrecv(1:3,i)
     mSN_i      = SNrecv(4,i)
     mloadSN_i  = SNrecv(5,i)
     ploadSN_i(1:3)=SNrecv(6:8,i)
     if(metal) ZloadSN_i = SNrecv(12,i)/mloadSN_i
#if NCHEM>0
     do ich=1,nchem
        chloadSN_i(ich) = SNrecv(13+ich,i)/mloadSN_i
     end do
#endif
     
     do j=1,nSNnei
 
        call get_icell_from_pos (xSN_i(1:3)+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell, ilevel2)
        if(cpu_map(father(igrid))==myid)then

           vol_nei = vol_loc*(2d0**ndim)**(ilevel-ilevel2)
           if(ilevel>ilevel2)then ! touching level-1 cells
              pvar(1:nvarMHD) = unew(icell,1:nvarMHD)
           else
              pvar(1:nvarMHD) = uold(icell,1:nvarMHD)
           endif
           do ii=i_fractions,nvar
              fractions(ii)=pvar(ii)/pvar(1)
           enddo
     
           emag=0.0d0
#ifdef SOLVERmhd
           do idim=1,ndim
              emag=emag+0.125d0*(pvar(ndim+2+idim)+pvar(nvar+idim))**2
           enddo
#endif
           d0=pvar(1)
           u0=pvar(2)/d0
           v0=pvar(3)/d0
           w0=pvar(4)/d0
           ekk0=0.5d0*d0*(u0**2d0 + v0**2d0 + w0**2d0)
           eth0=pvar(5)-ekk0-emag
#if NENER>0
           do irad=0,nener-1
              eth0 = eth0 - pvar(inener+irad)
           end do
#endif
           ! For stability
           Tk0 = eth0/d0*scale_T2*(gamma-1)
           T2min = T2_star*(d0*scale_nH/n_star)**(g_star-1.0)
           if(Tk0<T2min)then
              eth0 = T2min*d0/scale_T2/(gamma-1.0)*1.2195
           endif 

           d= mloadSN_i   /dble(nSNnei)/vol_nei
           u=(ploadSN_i(1)/dble(nSNnei)+p_solid(i,j)*vSNnei(1,j))/vol_nei/d
           v=(ploadSN_i(2)/dble(nSNnei)+p_solid(i,j)*vSNnei(2,j))/vol_nei/d
           w=(ploadSN_i(3)/dble(nSNnei)+p_solid(i,j)*vSNnei(3,j))/vol_nei/d

           pvar(1)=pvar(1)+d
           pvar(2)=pvar(2)+d*u
           pvar(3)=pvar(3)+d*v
           pvar(4)=pvar(4)+d*w
           !ekk_ej = 0.5*d*(u**2 + v**2 + w**2)
           etot0  = eth0+ekk0+emag+ek_solid(i,j)/vol_nei  ! additional energy from SNe+entrained gas

           ! the minimum thermal energy input floor
           d   = pvar(1)
           u   = pvar(2)/d
           v   = pvar(3)/d
           w   = pvar(4)/d
           ekk = 0.5*d*(u**2 + v**2 + w**2)

           pvar(5) = max(etot0, ekk+eth0+emag)

#if NENER>0
           ! Inject cosmic rays!!
           do irad=0,nener-1
              pvar(inener+irad) = pvar(inener+irad)+fecr*SNrecv(11,i)*E_SNII/scale_erg/dble(nSNnei)/vol_nei
              pvar(5) = pvar(5) + pvar(inener+irad)
           end do
#endif

           if(metal)then
              pvar(imetal)=pvar(imetal)+mloadSN_i/dble(nSNnei)*ZloadSN_i/vol_nei
           end if
#if NCHEM>0
           do ich=1,nchem
              pvar(ichem+ich-1)=pvar(ichem+ich-1)+mloadSN_i/dble(nSNnei)*chloadSN_i(ich)/vol_nei
           end do
#endif
           do ii=i_fractions,nvar
              pvar(ii)=fractions(ii)*pvar(1)
           end do

           ! update the hydro variable
           if(ilevel>ilevel2)then ! touching level-1 cells
              unew(icell,1:nvarMHD) = pvar(1:nvarMHD)
           else
              uold(icell,1:nvarMHD) = pvar(1:nvarMHD)
           endif
 
        endif ! if this belongs to me


     end do ! loop over neighbors
  end do ! loop over SN cells


  deallocate(icpuSN_comm_mpi, icpuSN_comm)
  if(ncell_send>0) deallocate(list2send,SNsend,reqsend,statsend)
  if(ncell_recv>0) deallocate(list2recv,SNrecv,reqrecv,statrecv)
  if(ncell_recv>0) deallocate(p_solid,ek_solid,snowplough)

  ncomm_SN=nSN_tot
#endif

end subroutine mech_fine_mpi
#endif
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine get_number_of_sn2(birth_time,dteff,zp_star,id_star,mass0,mass_t,nsn,done_star)
  use amr_commons, ONLY:dp,M_SNII,eta_sn
  use random
  implicit none
  real(kind=dp)::birth_time,zp_star,mass0,mass_t,dteff ! birth_time in code, mass in Msun
  real(kind=dp)::nsn
  integer::nsn_tot,nsn_sofar,nsn_age2
  integer::i,localseed,id_star   ! necessary for the random number
  real(kind=dp)::age1,age2,xdum,ydum,logzpsun
  real(kind=dp)::co0,co1,co2
  ! fit to the cumulative number fraction for Kroupa IMF
  ! ~kimm/soft/starburst99/output_hires_new/snrate.pro 
  ! ~kimm/soft/starburst99/output_hires_new/fit.pro 
  ! ~kimm/soft/starburst99/output_hires_new/sc.pro 
  real(kind=dp),dimension(1:3)::coa=(/-2.677292E-01,1.392208E-01,-5.747332E-01/)
  real(kind=dp),dimension(1:3)::cob=(/4.208666E-02, 2.152643E-02, 7.893866E-02/)
  real(kind=dp),dimension(1:3)::coc=(/-8.612668E-02,-1.698731E-01,1.867337E-01/)
  real(kind=dp),external::ran1 
  logical:: done_star


  ! determine the SNrateCum curve for Zp
  logzpsun=max(min(zp_star,0.05),0.008) ! no extrapolation
  logzpsun=log10(logzpsun/0.02)

  co0 = coa(1)+coa(2)*logzpsun+coa(3)*logzpsun**2
  co1 = cob(1)+cob(2)*logzpsun+cob(3)*logzpsun**2
  co2 = coc(1)+coc(2)*logzpsun+coc(3)*logzpsun**2

  ! RateCum = co0+sqrt(co1*Myr+co2)

  ! get stellar age
  call getStarAgeGyr(birth_time+dteff, age1)
  call getStarAgeGyr(birth_time      , age2)

  ! convert Gyr -> Myr
  age1 = age1*1d3
  age2 = age2*1d3

  nsn=0d0; done_star=.false.
  if(age2.le.(-co2/co1))then
     return
  endif

  ! total number of SNe
  nsn_tot   = NINT(mass0*(eta_sn/M_SNII),kind=4)
  if(nsn_tot.eq.0)then
      write(*,*) 'Fatal error: please increase the mass of your star particle'
      stop
  endif

  ! number of SNe up to this point
  nsn_sofar = NINT((mass0-mass_t)/M_SNII,kind=4)
  if(nsn_sofar.ge.nsn_tot)then
     done_star=.true.
     return
  endif

  localseed = -abs(id_star) !make sure that the number is negative

  nsn_age2 = 0
  do i=1,nsn_tot
     xdum =  ran1(localseed)
     ! inverse function for y=co0+sqrt(co1*x+co2)
     ydum = ((xdum-co0)**2.-co2)/co1
     if(ydum.le.age2)then
        nsn_age2=nsn_age2+1
     endif
  end do

  nsn_age2 = min(nsn_age2,nsn_tot)  ! just in case...
  nsn = max(nsn_age2 - nsn_sofar, 0) 

  if(nsn_age2.ge.nsn_tot) done_star=.true.

end subroutine get_number_of_sn2
#endif
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine get_number_of_sn2_bpass(birth_time,dteff,zp_star,id_star,mass0,mass_t,nsn,done_star)
  use amr_commons, ONLY:dp,M_SNII
  use mechanical_commons
  implicit none
  real(kind=dp)::birth_time,zp_star,mass0,mass_t,dteff ! birth_time in code, mass in Msun
  real(kind=dp)::nsn
  logical:: done_star
  integer::nsn_tot,nsn_age2,nsn_sofar,isn,ix,iZ(1)
  integer::localseed,id_star   ! necessary for the random number
  real(kind=dp)::age1,age2,log_yr,yr_sn,xran,yran,logzp,max_yr_sn
  real(kind=dp),external::ran1 

  if(variable_yield_SNII)then
     print *, 'ERR: not ready for the bpass with variable M_SNII'
     stop
  endif

  ! get stellar age
  call getStarAgeGyr(birth_time+dteff, age1)
  call getStarAgeGyr(birth_time      , age2)

  ! convert Gyr -> yr
  age1 = age1*1d9
  age2 = age2*1d9

  ! determine the SNrateCum curve for Zp
  logzp=log10(max(zp_star,0.0001))

  ! find the metallicity bin
  iZ = minloc(abs(Zgrid_bpass2-logzp)) 
  nsn_tot = nint(mass0 * snr_bpass2_sum(iZ(1)))
  nsn_sofar = nint( (mass0-mass_t)/M_SNII)

  if(nsn_tot.eq.0)then
      write(*,*) 'Fatal error: please increase the mass of your star particle'
      stop
  endif

  localseed = -abs(id_star) !make sure that the number is negative

  ! random sampling based on the rejection method
  nsn=0d0; isn=0; done_star=.false.; max_yr_sn=0d0; nsn_age2=0
  do while (isn.lt.nsn_tot) 
     xran = ran1(localseed)
     log_yr = xran* (logmax_yr_bpass2 - logmin_yr_bpass2) + logmin_yr_bpass2
     ix = int((log_yr - logmin_yr_bpass2)/binsize_yr_bpass2) + 1
     ix = min(max(1,ix),nt_bpass2)
     yran = ran1(localseed)*snr_bpass2_max(iZ(1))
     if (yran.le.snr_bpass2(ix,iZ(1)))then
        isn=isn+1
        yr_sn = 10d0**log_yr
        if(yr_sn.le.age2) nsn_age2 = nsn_age2+1
        if(max_yr_sn.lt.yr_sn) max_yr_sn=yr_sn
     endif 
  end do 

  nsn_age2 = min(nsn_age2,nsn_tot)  ! just in case...
  nsn = max(nsn_age2 - nsn_sofar, 0) 

  if(max_yr_sn.lt.age2) done_star=.true.

end subroutine get_number_of_sn2_bpass
#endif
!################################################################
!################################################################
!################################################################
!################################################################
function ran1(idum)
   implicit none
   integer:: idum,IA,IM,IQ,IR,NTAB,NDIV
   real(kind=8):: ran1,AM,EPS,RNMX
   parameter(IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,&
            &NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   integer::j,k,iv(NTAB),iy
   save iv,iy
   data iv /NTAB*0/, iy /0/ 
   if (idum.le.0.or.iy.eq.0) then! initialize
      idum=max(-idum,1)
      do j=NTAB+8,1,-1
         k=idum/IQ
         idum=IA*(idum-k*IQ)-IR*k
         if (idum.lt.0) idum=idum+IM
         if (j.le.NTAB) iv(j)=idum
      end do
      iy=iv(1)
   end if
   k=idum/IQ
   idum=IA*(idum-k*IQ)-IR*k
   if (idum.lt.0) idum=idum+IM
   j=1+iy/NDIV
   iy=iv(j)
   iv(j)=idum
   ran1=min(AM*iy,RNMX)
   return
end function ran1
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getStarAgeGyr(birth_time,age_star)
   use amr_commons
   implicit none
   real(dp)::age_star,age_star_pt,birth_time

   if (use_proper_time)then
      call getAgeGyr    (birth_time, age_star)
   else
      call getProperTime(birth_time, age_star_pt)
      call getAgeGyr    (age_star_pt,age_star)
   endif

end subroutine getStarAgeGyr
!################################################################
!################################################################
!################################################################
subroutine redundant_non_1d(list,ndata,ndata2)
   implicit none
   integer,dimension(1:ndata)::list,list2
   integer,dimension(1:ndata)::ind
   integer::ndata,i,y1,ndata2

   ! sort first
   call heapsort(list,ind,ndata) 

   list(:)=list(ind)

   ! check the redundancy
   list2(:) = list(:)

   y1=list(1)
   ndata2=1
   do i=2,ndata
       if(list(i).ne.y1)then
          ndata2=ndata2+1
          list2(ndata2)=list(i)
          y1=list(i)
       endif
   end do

   list =list2

end subroutine redundant_non_1d
!################################################################
!################################################################
!################################################################
subroutine heapsort(ain,ind,n)
   implicit none
   integer::n
   integer,dimension(1:n)::ain,aout,ind
   integer::i,j,l,ir,idum,rra

   l=n/2+1
   ir=n
   do i=1,n
      aout(i)=ain(i)                        ! Copy input array to output array
      ind(i)=i                                   ! Generate initial idum array
   end do
   if(n.eq.1) return                            ! Special for only one record
10 continue
   if(l.gt.1)then
      l=l-1
      rra=aout(l)
      idum=ind(l)
   else
      rra=aout(ir)
      idum=ind(ir)
      aout(ir)=aout(1)
      ind(ir)=ind(1)
      ir=ir-1
      if(ir.eq.1)then
        aout(1)=rra
        ind(1)=idum
        return
      endif
    endif
    i=l
    j=l+l
20  if(j.le.ir)then
       if(j.lt.ir)then
          if(aout(j).lt.aout(j+1))j=j+1
       endif
       if(rra.lt.aout(j))then
          aout(i)=aout(j)
          ind(i)=ind(j)
          i=j
          j=j+j
       else
          j=ir+1
       endif
       go to 20
    endif
    aout(i)=rra
    ind(i)=idum
    go to 10
end subroutine heapsort
!################################################################
!################################################################
!################################################################
subroutine mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)
  use amr_commons
  implicit none
  integer::ilevel,ind,ix,iy,iz,nx_loc
  real(dp)::skip_loc(1:3),scale,dx,dx_loc,vol_loc
  real(dp),dimension(1:twotondim,1:ndim):: xc

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

end subroutine mesh_info
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine get_icell_from_pos (fpos,ilevel_max,ind_grid,ind_cell,ilevel_out)
  use amr_commons
  implicit none
  real(dp)::fpos(1:3)
  integer ::ind_grid,ind_cell
  !-------------------------------------------------------------------
  ! This routnies find the index of the leaf cell for a given position
  ! fpos: positional info, from [0.0-1.0] (scaled by scale)
  ! ilevel_max: maximum level you want to search
  ! ind_cell: index of the cell
  ! ind_grid: index of the grid that contains the cell
  ! ilevel_out: level of this cell
  ! You can check whether this grid belongs to this cpu
  !     by asking cpu_map(father(ind_grid))
  !-------------------------------------------------------------------
  integer::ilevel_max
  integer::ilevel_out,i,ind,idim
  real(dp)::dx,fpos2(1:3)
  real(dp)::skip_loc(1:3),x0(1:3)
  logical ::not_found

  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
 
  fpos2=fpos
  if (fpos2(1).gt.1d0) fpos2(1)=fpos2(1)-1d0
  if (fpos2(2).gt.1d0) fpos2(2)=fpos2(2)-1d0
  if (fpos2(3).gt.1d0) fpos2(3)=fpos2(3)-1d0
  if (fpos2(1).lt.0d0) fpos2(1)=fpos2(1)+1d0
  if (fpos2(2).lt.0d0) fpos2(2)=fpos2(2)+1d0
  if (fpos2(3).lt.0d0) fpos2(3)=fpos2(3)+1d0

  not_found=.true.
  ind_grid =1  ! this is level=1 grid
  ilevel_out=1
  do while (not_found)
     dx = 0.5D0**ilevel_out
     x0 = xg(ind_grid,1:ndim)-dx-skip_loc  !left corner of *this* grid [0-1], not cell (in ramses, grid is basically a struture containing 8 cells)

     ind  = 1 
     do idim = 1,ndim
        i = int( (fpos2(idim) - x0(idim))/dx)
        ind = ind + i*2**(idim-1)
     end do

     ind_cell = ind_grid+ncoarse+(ind-1)*ngridmax
!     write(*,'(2(I2,1x),2(I10,1x),3(f10.8,1x))') ilevel_out,ilevel_max,ind_grid,ind_cell,fpos2
     if(son(ind_cell)==0.or.ilevel_out==ilevel_max) return

     ind_grid=son(ind_cell)
     ilevel_out=ilevel_out+1
  end do

end subroutine get_icell_from_pos
#endif
!################################################################
!################################################################
!################################################################
#if NDIM==3
subroutine SNII_yield (zp_star, ej_m, ej_Z, ej_chem)
  use amr_commons, ONLY:dp
  use mechanical_commons, ONLY:nchem,chem_list
#if NCHEM>0
  use hydro_parameters, ONLY:ichem
#endif
  implicit none
  real(dp)::zp_star,ej_m,ej_Z,ej_chem(1:nchem)
!-----------------------------------------------------------------
! Notice that even if the name is 'yield', 
! the return value is actually a metallicity fraction for simplicity
! These numbers are based on the outputs from Starburst99 
!                   (i.e. essentially Woosley & Weaver 95)
!-----------------------------------------------------------------
  real(dp),dimension(1:5)::log_SNII_m, log_Zgrid, log_SNII_Z
  real(dp),dimension(1:5)::log_SNII_H,log_SNII_He,log_SNII_C,log_SNII_N,log_SNII_O
  real(dp),dimension(1:5)::log_SNII_Mg,log_SNII_Si,log_SNII_S,log_SNII_Fe
  real(dp)::log_Zstar,fz
  integer::nz_SN=5, izg
#if NCHEM>0  
  integer::ich
  character(len=2)::element_name
  real(dp),dimension(1:5)::dum1d
#endif
  
  ! These are the numbers calculated from Starburst99 (Kroupa with 50Msun cut-off)
  ! (check library/make_stellar_winds.pro)
  log_SNII_m = (/-0.85591807,-0.93501857,-0.96138483,-1.0083450,-1.0544419/)
  log_Zgrid = (/-3.3979400,-2.3979400,-2.0969100,-1.6989700,-1.3010300/)
  log_SNII_Z=(/-0.99530662,-0.98262223,-0.96673581,-0.94018599,-0.93181853/)
  log_SNII_H=(/-0.28525316, -0.28988675, -0.29588988, -0.30967822, -0.31088065/)
  log_SNII_He=(/-0.41974152, -0.41688929, -0.41331511, -0.40330181, -0.40426739/)
  log_SNII_C=(/-1.9739731, -1.9726015, -1.9692265, -1.9626259, -1.9635311/)
  log_SNII_N=(/-3.7647616, -3.1325173, -2.8259748, -2.4260355, -2.1260417/)
  log_SNII_O=(/-1.1596291, -1.1491227, -1.1344617, -1.1213435, -1.1202089/)
  log_SNII_Mg=(/-2.4897201, -2.4368979, -2.4043654, -2.3706062, -2.3682933/)
  log_SNII_Si=(/-2.2157073, -2.1895132, -2.1518758, -2.0431845, -2.0444829/)
  log_SNII_S=(/-2.5492508, -2.5045016, -2.4482936, -2.2964020, -2.2988656/)
  log_SNII_Fe=(/-2.0502141, -2.0702598, -2.1074876, -2.2126987, -2.3480877/)

  ! search for the metallicity index
  log_Zstar = log10(zp_star)
  call binary_search(log_Zgrid, log_Zstar, nz_SN, izg )

  fz  = (log_Zgrid(izg+1) - log_Zstar )/( log_Zgrid(izg+1) - log_Zgrid(izg) )
  ! no extraploation
  if (fz  < 0.0) fz  = 0.0
  if (fz  > 1.0) fz  = 1.0


  ej_m = log_SNII_m(izg)*fz + log_SNII_m(izg+1)*(1d0-fz)
  ej_m = 10d0**ej_m

  ej_Z = log_SNII_Z(izg)*fz + log_SNII_Z(izg+1)*(1d0-fz)
  ej_Z = 10d0**ej_Z

#if NCHEM>0  
  do ich=1,nchem
      element_name=chem_list(ich)
      select case (element_name)
         case ('H ')
            dum1d = log_SNII_H
         case ('He')
            dum1d = log_SNII_He
         case ('C ')
            dum1d = log_SNII_C
         case ('N ')
            dum1d = log_SNII_N
         case ('O ')
            dum1d = log_SNII_O
         case ('Mg')
            dum1d = log_SNII_Mg
         case ('Si')
            dum1d = log_SNII_Si
         case ('S ')
            dum1d = log_SNII_S
         case ('Fe')
            dum1d = log_SNII_Fe
         case default
            dum1d = 0d0
      end select
     ej_chem(ich) = dum1d(izg)*fz + dum1d(izg+1)*(1d0-fz)
     ej_chem(ich) = 10d0**ej_chem(ich)
  end do
#endif

end subroutine SNII_yield
#endif
!################################################################
!################################################################
!################################################################
!################################################################
subroutine binary_search(database,xtarget,ndata,i)
   use amr_commons,ONLY:dp
   implicit none 
   integer::i,j,k
   integer,intent(in)::ndata
   real(dp),intent(in)::database(1:ndata),xtarget

   i=1  
   j=ndata
   do   
     k=(i+j)/2
     if (xtarget<database(k)) then 
         j=k  
     else 
         i=k  
     end if
     if (i+1>=j) exit 
   end do

end subroutine binary_search


