!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine
  use amr_commons
  use pm_commons
#ifdef SOLVERmhd
  use hydro_commons,only:napot
#endif
  implicit none
  !-------------------------------------------
  ! This routine builds the initial AMR grid
  !-------------------------------------------
  integer::ilevel

  if(myid==1)write(*,*)'Building initial AMR grid'
  init=.true.

#ifdef SOLVERmhd
  ! Generate the uniform grid of vector potentials
  if(napot>0)call generate_randomapot
#endif

  ! Base refinement
  do ilevel=1,levelmin
     call flag
     call refine
  end do

  ! Further refinements if necessary
  do ilevel=levelmin+1,nlevelmax
     if(initfile(levelmin).ne.' '.and.initfile(ilevel).eq.' ')exit
     if(hydro)call init_flow
#ifdef RT
     if(rt)call rt_init_flow
#endif
     if(ivar_refine==0)call init_refmap
     call flag
     call refine
     if(nremap>0)call load_balance
     if(numbtot(1,ilevel)==0)exit
  end do

  ! Final pass to initialize the flow
  init=.false.
  if(hydro)call init_flow
#ifdef RT
  if(rt)call rt_init_flow
#endif

end subroutine init_refine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_2
  !--------------------------------------------------------------
  ! This routine builds additional refinements to the
  ! the initial AMR grid for filetype ne 'grafic'
  !--------------------------------------------------------------
  use amr_commons
  use hydro_commons
#ifdef RT
  use rt_hydro_commons
#endif
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,i,ivar

  if(filetype.eq.'grafic')return

  do i=levelmin,nlevelmax+1

     call refine_coarse
     do ilevel=1,nlevelmax
        call build_comm(ilevel)
        call make_virtual_fine_int(cpu_map(1),ilevel)
        call refine_fine(ilevel)
        if(hydro)call init_flow_fine(ilevel)
#ifdef RT
        if(rt)call rt_init_flow_fine(ilevel)
#endif
     end do

     if(nremap>0)call load_balance

     do ilevel=levelmin,nlevelmax
        if(pic)call make_tree_fine(ilevel)
        if(poisson)call rho_fine(ilevel,2)
        if(pic)then
           call kill_tree_fine(ilevel)
           call virtual_tree_fine(ilevel)
        endif
     end do

     do ilevel=nlevelmax,levelmin,-1
        if(pic)call merge_tree_fine(ilevel)
        if(hydro)then
           call upload_fine(ilevel)
#ifdef SOLVERmhd
           do ivar=1,nvar+3
#else
           do ivar=1,nvar
#endif
              call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
           end do
#else
           end do
#endif
           if(simple_boundary)call make_boundary_hydro(ilevel)
        endif
#ifdef RT
        if(rt)then
           call rt_upload_fine(ilevel)
           do ivar=1,nrtvartrace
              call make_virtual_fine_rtdp(rtuold(1,ivar),ilevel)
           end do
           if(simple_boundary)call rt_make_boundary_hydro(ilevel)
        end if
#endif
     end do

     do ilevel=nlevelmax,1,-1
        call flag_fine(ilevel,2)
     end do
     call flag_coarse

  end do

#ifdef RT
  if(neq_chem .and. rt_is_init_xion) then
     if(myid==1) write(*,*) 'Initializing ionization states from T profile'
     do ilevel=nlevelmax,1,-1
        call rt_init_xion(ilevel)
        call upload_fine(ilevel)
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
        end do
#else
        end do
#endif
        if(simple_boundary)call make_boundary_hydro(ilevel)
     end do
  endif
#endif

#ifdef SOLVERmhd
  if(napot>0)deallocate(apotx,apoty,apotz)
#endif 
  
end subroutine init_refine_2

#ifdef SOLVERmhd
!################################################################
!################################################################
!################################################################
!################################################################
subroutine generate_randomapot
  use amr_commons
  use hydro_commons, ONLY: apotx,apoty,apotz,napot,napotsmall,ncellpot,seedpot
  implicit none
  integer::i,j,k,icelll,icellr,napot_loc,size
  integer, allocatable :: seed(:)

  if(verbose)write(*,*)'Entering generate_randomA'

  if(myid==1)write(*,*)'Asking for random A pot on a grid with ncells (1D)=',napot
  if(napot>napotsmall)then
     if(myid==1)write(*,*)'Changing for random A pot on a grid with ncells (1D)=',napotsmall
     if(myid==1)write(*,*)'replicated ',napot/napotsmall,' times'
     napot_loc=napotsmall
  else
     napot_loc=napot
  endif
  ncellpot=(napot_loc+1)*(napot_loc+1)*(napot_loc+1)

  allocate(apotx(ncellpot))
  if(myid==1)write(*,'(f8.6,A)')sizeof(apotx)*3d0/1d9,' GB of memory being used for Apot'
  allocate(apoty(ncellpot))
  allocate(apotz(ncellpot))

  call random_seed(size = size)
  allocate(seed(size))
  seed(:)=seedpot(1)
  call random_seed(put=seed)
  call random_number(apotx)
  call random_number(apoty)
  call random_number(apotz)

  ! Copy the left-most faces into the right-most faces (x-direction)
  do j=1,napot_loc+1
  do k=1,napot_loc+1
     icellr=icell_apot(napot_loc+1,j,k,napot_loc+1)
     icelll=icell_apot(1          ,j,k,napot_loc+1)
     apotx(icellr)=apotx(icelll)
     apoty(icellr)=apoty(icelll)
     apotz(icellr)=apotz(icelll)
  enddo
  enddo
  ! Copy the front-most faces into the back-most faces (y-direction)
  do i=1,napot_loc+1
  do k=1,napot_loc+1
     icellr=icell_apot(i,napot_loc+1,k,napot_loc+1)
     icelll=icell_apot(i,1          ,k,napot_loc+1)
     apotx(icellr)=apotx(icelll)
     apoty(icellr)=apoty(icelll)
     apotz(icellr)=apotz(icelll)
  enddo
  enddo
  ! Copy the bottom-most faces into the top-most faces (z-direction)
  do i=1,napot_loc+1
  do j=1,napot_loc+1
     icellr=icell_apot(i,j,napot_loc+1,napot_loc+1)
     icelll=icell_apot(i,j,1          ,napot_loc+1)
     apotx(icellr)=apotx(icelll)
     apoty(icellr)=apoty(icelll)
     apotz(icellr)=apotz(icelll)
  enddo
  enddo

  apotx=apotx/dble(napot_loc)
  apoty=apoty/dble(napot_loc)
  apotz=apotz/dble(napot_loc)

#if NDIM==2
  apotx=0.0d0
  apoty=0.0d0
#endif

contains

function icell_apot(i,j,k,nx)
  implicit none
  integer::i,j,k,nx
  integer::icell_apot

  icell_apot=i+(j-1)*nx+(k-1)*nx*nx
end function icell_apot



end subroutine generate_randomapot
#endif
