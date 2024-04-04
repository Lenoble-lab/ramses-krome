!************************************************************************
SUBROUTINE rt_godunov_fine(ilevel, dt)

! This routine is a wrapper to the grid solver for radiative transfer.
! Small grids (2x2x2) are gathered from level ilevel and sent to the
! RT solver. On entry, RT variables are gathered from array rtuold.
! On exit, rtunew has been updated.
!------------------------------------------------------------------------
  use amr_commons
  use rt_hydro_commons
  implicit none
  integer::ilevel
  integer::i,igrid,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid
  real(dp)::dt
!------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return  ! # of grids at ilevel
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid  ! total # of grids at level ilevel
  do igrid=1,ncache,nvector    ! take steps of 500 grids up to ncache
     ngrid=MIN(nvector,ncache-igrid+1) ! # of grids in each sweep
     do i=1,ngrid              ! collect grid indices for one sweep
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call rt_godfine1(ind_grid,ngrid,ilevel, dt) ! solve for max 500 grids
  end do

  if(verbose)write(*,112)ilevel

111 format('   Entering rt_godunov_fine for level ',i2)
112 format('   Exiting rt_godunov_fine for level ',i2)

END SUBROUTINE rt_godunov_fine

!************************************************************************
SUBROUTINE rt_set_unew(ilevel)

! This routine sets array rtunew to its initial value rtuold before calling
! the rt scheme. rtunew is set to zero in virtual boundaries.
!------------------------------------------------------------------------
  use amr_commons
  use rt_hydro_commons
  implicit none
  integer::ilevel
  integer::i,ivar,ind,icpu,iskip
!------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set rtunew to rtuold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nrtvartrace
        do i=1,active(ilevel)%ngrid
           rtunew(active(ilevel)%igrid(i)+iskip,ivar) = rtuold(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do
  end do

  ! Set rtunew to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nrtvartrace
        do i=1,reception(icpu,ilevel)%ngrid
           rtunew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
        end do
     end do
  end do
  end do

111 format('   Entering rt_set_unew for level ',i2)

END SUBROUTINE rt_set_unew

!************************************************************************
SUBROUTINE rt_set_uold(ilevel)

! This routine sets array rtuold to its new value rtunew after the
! hydro step.
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters
  use rt_hydro_commons
  implicit none
  integer::ilevel
  integer::i,ivar,ind,iskip,ig,icell
  real(dp)::Npc,fred
!------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! If rt_smooth==true, then cooling_fine takes care of this, adding to
  ! the RT variables incrementally during the cooling-substeps.
  ! This is useful to keep down the cooling-substepping.

  ! We do have to update the tracers here
  ! Set rtuold to rtunew for myid cells
#ifdef NTRACEGROUPS
  if (ntracegroups.gt.0) then
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=nrtvar+1,nrtvartrace
           do i=1,active(ilevel)%ngrid
              rtuold(active(ilevel)%igrid(i)+iskip,ivar) = rtunew(active(ilevel)%igrid(i)+iskip,ivar)
           end do
        end do
     end do
  endif
#endif
  if(rt_smooth) return

  ! Set rtuold to rtunew for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nrtvar
        do i=1,active(ilevel)%ngrid
           rtuold(active(ilevel)%igrid(i)+iskip,ivar) = rtunew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
     end do

     ! Make a photon conservation fix (prevent beam induced crash)
     do ig=1,nGroups
        do i=1,active(ilevel)%ngrid
           icell=active(ilevel)%igrid(i)+iskip
           ! No negative photon densities:
           rtuold(icell,iGroups(ig)) = max(rtuold(icell,iGroups(ig)),smallNp)
           Npc=rtuold(icell,iGroups(ig))*rt_c(ilevel)
           ! Reduced flux, should always be .le. 1
           fred = sqrt(sum((rtuold(icell,iGroups(ig)+1:iGroups(ig)+ndim))**2))/Npc
           if(fred .gt. 1d0) then ! Too big so normalize flux to one
              rtuold(icell,iGroups(ig)+1:iGroups(ig)+ndim) &
                   = rtuold(icell,iGroups(ig)+1:iGroups(ig)+ndim)/fred
           endif
        end do
     end do
     ! End photon conservation fix

  end do

111 format('   Entering rt_set_uold for level ',i2)

END SUBROUTINE rt_set_uold

!*************************************************************************
SUBROUTINE rt_godfine1(ind_grid, ncache, ilevel, dt)

! This routine first gathers RT variables from neighboring grids
! to set initial conditions in a 6x6x6 grid. It interpolates from
! coarser level to missing grid variables. It then calls the solver
! that computes fluxes. These fluxes are zeroed at coarse-fine boundaries,
! since contribution from finer levels has already been taken into
! account. Conservative variables are updated and stored in array rtunew(:),
! both at the current level and at the coarser level if necessary.
!
! in ind_grid: Indexes of grids/octs to solve in
! in ncache:   Length of ind_grid (i.e. number of grids)
! in ilevel:   Level at which the grids are
! in dt:       Timestep length
!-------------------------------------------------------------------------
  use amr_commons
  use rt_hydro_commons
  use rt_flux_module
  use rt_parameters
  implicit none
  integer::ilevel,ncache
  real(dp)::dt
  integer   ,dimension(1:nvector)::ind_grid
  integer   ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer   ,dimension(1:nvector,1:twotondim  ),save::nbors_father_grids
  integer   ,dimension(1:nvector,0:twondim    ),save::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim  ,1:nrtvartrace),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nrtvartrace),save::u2

  ! 500*6*6*6*4:
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nrtvartrace),save::uloc
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nrtvar,1:ndim),save::flux
  logical,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok
  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer
  integer,dimension(1:nvector),save::ind_exist,ind_nexist

  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale,oneontwotondim

  logical,dimension(1:nvector),save:: rt_per_bnd=.false.
  integer::ind_nbor, tracegroup
  real(dp)::dx8,maxDist,rt_c_diff

#ifdef NTRACEGROUPS
  !Photon tracer stuff
  integer::hindex, hstart, hfin
  real(dp),dimension(ntraceGroups)::tracedn,traceup,tracecell
  real(dp)::cellphotdens,fluxup,fluxdn
#endif
!------------------------------------------------------------------------
  oneontwotondim = 1d0/dble(twotondim) ! 1/8 in 3D

  ! Mesh spacing
  nx_loc=icoarse_max-icoarse_min+1    ! =1
  scale=boxlen/dble(nx_loc)           ! length per coarse oct (=boxlen)
  dx=0.5D0**ilevel*scale              ! length per oct/grid at ilevel
  dx8=8.*dx                           ! Outflow boundary

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
  ! in 3D:
  !            min      max    tot
  ! -----------------------------------------
  ! i,j,k1      0        2      3
  ! i,j,k2      0        1      2
  ! i,j,k3      1        2      2

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  ! of grid (ilevel-1). ind_cell are indexes
  ! in rtuold.
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  ! ..and father cells of neighbor grids:
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  ! now for the parent cell (ind_cell(i)) of each grid i in cache:
  !
  ! nbors_father_cells contains indexes of all it's neighbor cells
  ! (ilevel-1), plus itself, total 3^ndim.
  !
  ! nbors_father_grids contains indexes of all the neighboring
  ! (and containing) grids of the father cell, total 2^ndim
  ! (in case interpolation is needed I guess)

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells (ilevel-1)
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max

     ! Check if neighbor has a grid (or is just a cell)

     ! # of neighbors (out of all the ncache cells) at k1, j1, i1
     ! that are only a cell (not a grid):
     nbuffer=0
     ! # of neighbors (out of all the ncache cells) at k1, j1, i1
     ! that have subgrids:
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
     !--------------------------------------------------------------------
     ! we should now have: nexist+nbuffer==ncache
     !
     ! ind_exist has size nexist and contains indexes of cells in cache
     ! that have a neighboring father grid at k1, j1, i1
     !
     ! ind_nexist has size nbuffer and contains indexes of cells in cache
     ! that only have a neighboring father cell at k1, j1, i1
     !
     ! ind_buffer has size nbuffer and contains tree indexes of these
     ! father nongrid-cells
     !--------------------------------------------------------------------

     if(nbuffer>0)then
        ! For those nongrid cells we interpolate variables from parent cells:
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nrtvartrace
              do i=1,nbuffer
                 u1(i,j,ivar)=rtuold(ibuffer_father(i,j),ivar)
              end do
           end do
        end do
        call rt_interpol_hydro(u1,u2,nbuffer)
     endif

     ! RT ouflow boundary begin-------------------------------------------
     if(rt_is_outflow_bound) then
        rt_per_bnd(:)=.false.           ! Flag indicates periodic boundary
        do i=1,ncache
           ind_nbor = son(nbors_father_cells(i,ind_father))
           if(ind_nbor .le. 0) then
              ! Neighbor is leaf.......need container grid to get position
              ind_nbor = (nbors_father_cells(i,ind_father)-ncoarse-1)    &
                   / ngridmax + 1
              ind_nbor = nbors_father_cells(i,ind_father)                &
                   - ncoarse - (ind_nbor-1)*ngridmax
           endif
           maxDist = MAXVAL(ABS( xg(ind_grid(i),:)- xg(ind_nbor,:) ))
           if( maxDist .gt. dx8 ) then
              rt_per_bnd(i)=.true.  ! Neighbor is across periodic boundary
           else
              rt_per_bnd(i)=.false.
           endif
        end do
     endif
     ! RT outflow boundary end---------------------------------------------

     ! Loop 2x2x2 cells within father cell and add them to stencil (uloc)
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do

        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2 ! From -1 to 4 over outer loop, but
        if(ndim>1)j3=1+2*(j1-1)+j2 ! only over 2x2x2 indexes in inner loop
        if(ndim>2)k3=1+2*(k1-1)+k2

        ! Gather hydro variables

        ! RT outflow boundary begin---------------------------------------
        if(rt_is_outflow_bound) then
           do i=1,nexist
              if (rt_per_bnd(ind_exist(i))) then
                 uloc(ind_exist(i),i3,j3,k3,1:nrtvartrace) =  0d0
              else              
                 do ivar=1,nrtvartrace
                    uloc(ind_exist(i),i3,j3,k3,ivar) =                   &
                                                  rtuold(ind_cell(i),ivar)
                    if(son(ind_cell(i))>0                                &
                      .and. mod(ivar,ndim+1)==1                          &
                      .and. ivar.le.nrtvar) then
                       ! VSLA: finer level and different light speed
                       uloc(ind_exist(i),i3,j3,k3,ivar)        &
                            = uloc(ind_exist(i),i3,j3,k3,ivar) &
                            * rt_c(ilevel+1)/rt_c(ilevel)
                    endif
                 end do
              endif
           end do
           do i=1,nbuffer
              if (rt_per_bnd(ind_nexist(i))) then
                 uloc(ind_nexist(i),i3,j3,k3,1:nrtvartrace) =  0d0
              else
                 do ivar=1,nrtvartrace
                    uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
                    if(mod(ivar,ndim+1)==1 .and. ivar.le.nrtvar) then
                       ! VSLA: coarser level and different light speed
                       uloc(ind_nexist(i),i3,j3,k3,ivar)        &
                            = uloc(ind_nexist(i),i3,j3,k3,ivar) &
                            * rt_c(ilevel-1)/rt_c(ilevel)
                    endif
                 end do
              endif
           end do
        else
           do ivar=1,nrtvartrace
              do i=1,nexist
                 uloc(ind_exist(i),i3,j3,k3,ivar)=rtuold(ind_cell(i),ivar)
                 if(son(ind_cell(i))>0 .and. mod(ivar,ndim+1)==1 .and. ivar.le.nrtvar) then
                    ! VSLA: finer level and different light speed
                    uloc(ind_exist(i),i3,j3,k3,ivar) &
                         = uloc(ind_exist(i),i3,j3,k3,ivar) &
                         * rt_c(ilevel+1)/rt_c(ilevel)
                 endif
              end do
              do i=1,nbuffer
                 uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
                 if(mod(ivar,ndim+1)==1 .and. ivar.le.nrtvar) then
                    ! VSLA: coarser level and different light speed
                    uloc(ind_nexist(i),i3,j3,k3,ivar)        &
                         = uloc(ind_nexist(i),i3,j3,k3,ivar) &
                         * rt_c(ilevel-1)/rt_c(ilevel)
                    
                 endif
              end do
           end do
        endif
        ! RT outflow boundary end-----------------------------------------

        ! Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do

     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

  !----------------------------------------------------------------------
  ! Compute fluxes of each photon group, using Eddington tensor
  !----------------------------------------------------------------------
  do i = 1,nGroups
     !We don't want to send photon tracers into this function
     call cmp_rt_faces(uloc(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nrtvar) &
              ,flux,dx,dx,dx,dt,iGroups(i),ncache,ilevel)
  end do

  !----------------------------------------------------------------------
  ! Reset flux along direction at refined interface, if not rt-subcycling
  !----------------------------------------------------------------------
  if (rt_nsubcycle == 1)then
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        !flux only has dimension nrtvar so we dont loop over nrtvartrace
        do ivar=1,nrtvar
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 flux(i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
     end do
     end do
     end do
  end do
  endif
  !--------------------------------------
  ! Conservative update at level ilevel
  !--------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max      ! all from 0 to 1
     do j2=j2min,j2max      ! => update 2x2x2 cells
     do i2=i2min,i2max      ! i.e. update one grid at level ilevel
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2   ! just because the flux indexes are (1:3), not (0:2)
        k3=1+k2
        ! Update conservative variables new state vector
        tracegroup = 0
        do ivar=1,nrtvar
           if (mod(ivar-1,ndim+1).eq.0) tracegroup = tracegroup + 1 !check if we are looking at photon number density
           do i=1,ncache
#ifdef NTRACEGROUPS
              if (ntraceGroups.gt.0) then
                 if (ivar.eq.iGroups(tracegroup)) then                !This means we are updating the number density
                    hindex = ((ivar - 1)/(nDim+1))*ntraceGroups      !Offset from the generic starting index for this photon group
                    hstart = (nGroups*(nDim+1)) + 1                  !The generic starting index of the photon tracers 
                    hstart = hstart+hindex                           !The starting index of the photon tracers for this photon group
                    hfin   = hstart+ntraceGroups-1                   !The ending index of the photon tracers for this photon group
                    cellphotdens = rtunew(ind_cell(i),ivar)          !Number density of photons in the cell before the update
                    fluxdn = flux(i,i3   ,j3   ,k3   ,ivar,idim)     !Photon flux across the left boundary
                    fluxup = flux(i,i3+i0,j3+j0,k3+k0,ivar,idim)     !Photon flux across the right boundary
                    tracedn = uloc(i,i3-i0,j3-j0,k3-k0,hstart:hfin)  !Tracer fractions in the cell to the left
                    traceup = uloc(i,i3+i0,j3+j0,k3+k0,hstart:hfin)  !Tracer fractions in the cell to the right
                    tracecell = rtunew(ind_cell(i),hstart:hfin)      !Tracer fractions in the 
                    call update_tracers(cellphotdens,fluxdn,fluxup,tracedn,traceup,tracecell) !Calculate the updated tracer fractions
                    rtunew(ind_cell(i),hstart:hfin) = tracecell      !Set the host cell value to the new tracer fractions
                 end if
              end if ! if(ntracegroups.gt.0)
#endif
              rtunew(ind_cell(i),ivar)=                     &
                   &  rtunew(ind_cell(i),ivar)              &
                   & +(flux(i,i3   ,j3   ,k3   ,ivar,idim)  &
                   & - flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
              !if (mod(ivar,ndim+1)==1) rtunew(ind_cell(i),ivar) = max(rtunew(ind_cell(i),ivar),smallNp) ! Light explosion problem here
           end do
        end do
     end do
     end do
     end do
  end do

  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  if (rt_nsubcycle == 1)then
     ! Loop over dimensions
     do idim=1,ndim
        i0=0; j0=0; k0=0
        if(idim==1)i0=1
        if(idim==2)j0=1
        if(idim==3)k0=1

        !----------------------
        ! Left flux at boundary
        !----------------------
        ! Check if grids sits near left boundary
        ! and gather neighbor father cells index
        nb_noneigh=0
        do i=1,ncache
           if (son(nbor(ind_grid(i),2*idim-1))==0) then
              nb_noneigh = nb_noneigh + 1
              ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
              ind_cell(nb_noneigh) = i
           end if
        end do
        ! Conservative update of new state variables
        tracegroup = 0
        do ivar=1,nrtvar
           ! For VSLA, when updating coarser level, rescale radiation flux
           ! to the expression which would be seen from there.
           rt_c_diff = rt_c(ilevel-1)/rt_c(ilevel)
           if (mod(ivar,ndim+1)==1) then
              rt_c_diff=1.d0
              tracegroup = tracegroup + 1
           end if
           ! Loop over boundary cells
           do k3=k3min,k3max-k0 ! 1 to 1 if dim=3, 1 to 2 otherwise
              do j3=j3min,j3max-j0 ! 1 to 1 if dim=2, 1 to 2 otherwise
                 do i3=i3min,i3max-i0 ! 1 to 1 if dim=1, 1 to 2 otherwise
                    do i=1,nb_noneigh
#ifdef NTRACEGROUPS
                       if (ntraceGroups.gt.0) then
                          if (ivar.eq.iGroups(tracegroup)) then 
                             hindex = ((ivar - 1)/(nDim+1))*ntraceGroups
                             hstart = (nGroups*(nDim+1)) + 1 
                             hstart = hstart+hindex
                             hfin   = hstart+ntraceGroups-1  
                             !Only update the fractions if photons are moving into the cell
                             if (flux(ind_cell(i),i3,j3,k3,ivar,idim).lt.0.0d0) then
                                cellphotdens = rtunew(ind_buffer(i),ivar)
                                fluxdn = 0.0d0                                                !In this case, there is no down flux so we set it to 0
                                fluxup = flux(ind_cell(i),i3,j3,k3,ivar,idim)*oneontwotondim  !Photon flux in the refined cell moving into the coarse cell
                                tracedn = 1.d0/ntracegroups                                   !There is no flux but set the tracer fractions to 1/ntracegroup (dummy argument)
                                traceup = uloc(ind_cell(i),i3,j3,k3,hstart:hfin)              !Tracer fractions in the refined cell 
                                tracecell = rtunew(ind_buffer(i),hstart:hfin)                 !Tracer fractions in the coarse cell   
                                call update_tracers(cellphotdens,fluxdn,fluxup,tracedn,traceup,tracecell,1,ilevel)!Calculate the updated tracer fractions
                                rtunew(ind_buffer(i),hstart:hfin) = tracecell                 !Set the host cell value to the new tracer fractions
                             endif
                          endif
                       endif
#endif 
                       rtunew(ind_buffer(i),ivar) =                      &
                           & rtunew(ind_buffer(i),ivar)                  &
                           & - flux(ind_cell(i),i3,j3,k3,ivar,idim)      &
                           & * oneontwotondim * rt_c_diff
                       !if (mod(ivar,ndim+1)==1) rtunew(ind_buffer(i),ivar) = max(rtunew(ind_buffer(i),ivar),smallNp) !harley edit for neg phot dens
                    end do
                 end do
              end do
           end do
        end do

        !-----------------------
        ! Right flux at boundary
        !-----------------------
        ! Check if grids sits near right boundary
        ! and gather neighbor father cells index
        nb_noneigh=0
        do i=1,ncache
           if (son(nbor(ind_grid(i),2*idim))==0) then
              nb_noneigh = nb_noneigh + 1
              ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
              ind_cell(nb_noneigh) = i
           end if
        end do
        ! Conservative update of new state variables
        tracegroup = 0
        do ivar=1,nrtvar
           ! Rescale for VSLA, as above
           rt_c_diff = rt_c(ilevel-1)/rt_c(ilevel)
           if (mod(ivar,ndim+1)==1) then
              rt_c_diff=1.d0
              tracegroup = tracegroup + 1
           end if
           ! Loop over boundary cells
           do k3=k3min+k0,k3max
              do j3=j3min+j0,j3max
                 do i3=i3min+i0,i3max
                    do i=1,nb_noneigh
#ifdef NTRACEGROUPS
                       if (ntraceGroups.gt.0) then
                          if (ivar.eq.iGroups(tracegroup)) then 
                             hindex = ((ivar - 1)/(nDim+1))*ntraceGroups
                             hstart = (nGroups*(nDim+1)) + 1 
                             hstart = hstart+hindex
                             hfin   = hstart+ntraceGroups-1
                             if (flux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim).gt.0.0d0) then
                                cellphotdens = rtunew(ind_buffer(i),ivar)
                                fluxdn = flux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim)*oneontwotondim   !Photon flux in the refined cell moving into the coarse cell
                                fluxup = 0.0d0                                                          !In this case, there is no up flux so we set it to 0
                                tracedn = uloc(ind_cell(i),i3,j3,k3,hstart:hfin)                        !Tracer fractions in the refined cell 
                                traceup = 1.d0/ntracegroups                                             !There is no flux but set the tracer fractions to 1/ntracegroup (dummy argument) 
                                tracecell = rtunew(ind_buffer(i),hstart:hfin)                           !Tracer fractions in the coarse cell  
                                call update_tracers(cellphotdens,fluxdn,fluxup,tracedn,traceup,tracecell,1,ilevel)!Calculate the updated tracer fractions
                                rtunew(ind_buffer(i),hstart:hfin) = tracecell                           !Set the host cell value to the new tracer fractions
                             endif
                          endif
                       endif
#endif
                       rtunew(ind_buffer(i),ivar) =                          &
                           & rtunew(ind_buffer(i),ivar)                      &
                           & + flux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim) &
                           & * oneontwotondim * rt_c_diff
                       !if (mod(ivar,ndim+1)==1) rtunew(ind_buffer(i),ivar) = max(rtunew(ind_buffer(i),ivar),smallNp) !harley edit for neg phot dens
                    end do
                 end do
              end do
           end do
        end do

     end do
     ! End loop over dimensions
  end if
  ! End if-clause for rt-subcycling
END SUBROUTINE rt_godfine1

#ifdef NTRACEGROUPS
SUBROUTINE update_tracers(Nphotcell,ndn,nup,fracdn,fracup,fraccell)
  use amr_commons
  use rt_parameters
  implicit none
  real(dp)::Nphotcell,ndn,nup,npcellnewtot,Np_MIN
  real(dp),dimension(1:ntraceGroups)::fracdn,fracup,fraccell,npdn,npup,npcell,npcellnewgroup
  integer::i

  !Np_MIN=1.d-13 !Photon density floor (same as rt_cooling_module.f90)
  Np_MIN=0.0d0

  !total number of photons in the cell after the update
  Nphotcell = max(Nphotcell,0.0d0)
  npcellnewtot = Nphotcell+(ndn-nup)

  !this is an array of the number of photons in each 'tracer group' from the down flux
  if (ndn.lt.0d0) then !photons are leaving the cell
     npdn = ndn*fraccell
  else !photons are moving into the cell
     npdn = ndn*fracdn
  end if

  !this is an array of the number of photons in each 'tracer group' from the up flux 
  !it is negative because if nup is positive, photons are moving out of the cell
  if (nup.gt.0d0) then !photons are leaving the cell
     npup = nup*fraccell 
  else !photons are moving into the cell
     npup = nup*fracup
  endif

  !this is an array of the number of photons in each 'tracer group' for the current cell 
  npcell = Nphotcell*fraccell


  !total number of photons in each group after the update
  npcellnewgroup = npcell+(npdn-npup)
  
  !Check to make sure that we dont get negative fractions
  !If we have, leave fractions as they were and don't update (this only happens with AMR)
  !Note that these get updated very quickly once actual photons move into the cell
  !Also note that this doesn't affect the photoionization rate calculation since Np=0
  if (npcellnewtot.le.Np_MIN) then
     !write(*,*) 'Error: Your number density has gone negative in rt_godunov_fine',npcellnewtot
     !fraccell = 1.d0/ntraceGroups
     fraccell = fraccell
  else
     fraccell = npcellnewgroup/npcellnewtot
  endif

  !Next, check to make sure that none of the fractions have gone too small
  do i=1,ntracegroups
     if (fraccell(i).lt.Np_MIN) then
        fraccell(i) = Np_MIN
     endif
  end do

  !Since we have adjusted the cells which have become small, we need to make sure 
  !that the sum is still 1
  fraccell = fraccell/sum(fraccell)

END SUBROUTINE update_tracers
#endif
