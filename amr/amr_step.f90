recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
#ifdef RT
  use rt_hydro_commons
  use SED_module
  use UV_module
  use coolrates_module, only: update_coolrates_tables
  use rt_cooling_module, only: update_UVrates
#endif
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::mpi_err
#endif
  integer::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first.                          !
  !-------------------------------------------------------------------!
  integer::i,idim,ivar
  logical::ok_defrag,output_now_all
  logical,save::first_step=.true.
#ifdef SOLVERmhd
  integer::itemp
#endif
#ifdef NTRACEGROUPS
  real(dp) :: adiff, adiffnew
  integer  :: minind_temp, jj
  logical  ::okhprop
  integer:: hnhs, ih
  character(len=256)::hhprops                             ! Filenames
#endif

  if(numbtot(1,ilevel)==0)return

  if(verbose)write(*,999)icount,ilevel

  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
                               call timer('refine','start')
  if(levelmin.lt.nlevelmax .and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then
     if(ilevel==levelmin.or.icount>1)then
        do i=ilevel,nlevelmax
           if(i>levelmin)then

              !--------------------------
              ! Build communicators
              !--------------------------
              call build_comm(i)

              !--------------------------
              ! Update boundaries
              !--------------------------
              call make_virtual_fine_int(cpu_map(1),i)
              if(hydro)then
#ifdef SOLVERmhd
                 do ivar=1,nvar+3
#else
                 do ivar=1,nvar
#endif
                    call make_virtual_fine_dp(uold(1,ivar),i)
#ifdef SOLVERmhd
                 end do
#else
                 end do
#endif
                 if(momentum_feedback>0)call make_virtual_fine_dp(pstarold(1),i)
                 if(simple_boundary)call make_boundary_hydro(i)
              end if
#ifdef RT
              if(rt)then
                 do ivar=1,nrtvartrace
                    call make_virtual_fine_rtdp(rtuold(1,ivar),i)
                 end do
                 if(simple_boundary)call rt_make_boundary_hydro(i)
              end if
#endif
              if(poisson)then
                 call make_virtual_fine_dp(phi(1),i)
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),i)
                 end do
                 if(simple_boundary)call make_boundary_force(i)
              end if
           end if

           !--------------------------
           ! Refine grids
           !--------------------------
           call refine_fine(i)
        end do
     end if
  end if

  !--------------------------
  ! Load balance
  !--------------------------
                               call timer('load balance','start')
  ok_defrag=.false.
  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin)then
        if(nremap>0)then
           ! Skip first load balance because it has been performed before file dump
           if(nrestart>0.and.first_step)then
              if(nrestart.eq.nrestart_quad) restart_remap=.true.
              if(restart_remap) then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
              first_step=.false.
           else
              if(MOD(nstep_coarse,nremap)==0)then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
           end if
        end if
     endif
  end if

  !-----------------
  ! Update sink cloud particle properties
  !-----------------
#if NDIM==3
                               call timer('sinks','start')
  if(sink)call update_cloud(ilevel)
#endif
  !-----------------
  ! Particle leakage
  !-----------------
                               call timer('particles','start')
  if(pic)call make_tree_fine(ilevel)

  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==levelmin)then

#ifdef WITHOUTMPI
     output_now_all = output_now
#else
     ! check if any of the processes received a signal for output
     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_ALLREDUCE(output_now,output_now_all,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,mpi_err)
#endif
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout).or.output_now_all.EQV..true.)then
                               call timer('io','start')
        if(.not.ok_defrag)then
           call defrag
        endif

        ! Run the clumpfinder, (produce output, don't keep arrays alive on output)
        ! CAREFUL: create_output is used to destinguish between the case where
        ! the clumpfinder is called from create_sink or directly from amr_step.
#if NDIM==3
        if(clumpfind .and. ndim==3) call clump_finder(.true.,.false.)
#endif

        call dump_all


        ! Dump lightcone
        if(lightcone .and. ndim==3) call output_cone()

        if (output_now_all.EQV..true.) then
          output_now=.false.
        endif

     endif

  endif

  !----------------------------
  ! Output frame to movie dump (without synced levels)
  !----------------------------
  if(movie) then
     if(imov.le.imovout)then
        if(aexp>=amovout(imov).or.t>=tmovout(imov))then
                               call timer('movie','start')
           call output_frame()
        endif
     endif
  end if

  !-----------------------------------------------------------
  ! Put here all stuffs that are done only at coarse time step
  !-----------------------------------------------------------
  if(ilevel==levelmin)then
     !----------------------------------------------------
     ! Kinetic feedback from giant molecular clouds
     !----------------------------------------------------
                               call timer('feedback','start')
     if(hydro.and.star.and.eta_sn>0.and.f_w>0 .and. trim(feedback_model).eq.'kinetic') &
          call kinetic_feedback
  endif

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
                               call timer('poisson','start')
     !save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)
                               call timer('rho','start')
     call rho_fine(ilevel,icount)
  endif

  !-------------------------------------------
  ! Sort particles between ilevel and ilevel+1
  !-------------------------------------------
  if(pic)then
     ! Remove particles to finer levels
                               call timer('particles','start')
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
     call virtual_tree_fine(ilevel)
  end if

  if(poisson)then ! necessary because we want make_virtual_fine_dp in line 275 
     ! Mechanical feedback from stars 
     if(hydro.and.star.and.trim(feedback_model).eq.'mechanical') then
                               call timer('feedback','start')
        call mechanical_feedback_fine(ilevel,icount)
        do ivar=1,nvar
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
     endif
  endif

  !---------------
  ! Gravity update
  !---------------
  if(poisson)then
                               call timer('poisson','start')

     ! Remove gravity source term with half time step and old force
     if(hydro)then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif

     ! Compute gravitational potential
     if(ilevel>levelmin)then
        if(ilevel .ge. cg_levelmin) then
           call phi_fine_cg(ilevel,icount)
        else
           call multigrid_fine(ilevel,icount)
        end if
     else
        call multigrid_fine(levelmin,icount)
     end if
     !when there is no old potential...
     if (nstep==0)call save_phi_old(ilevel)

     ! Compute gravitational acceleration
     call force_fine(ilevel,icount)

     ! Synchronize remaining particles for gravity
     if(pic)then
                               call timer('particles','start')
        if(static_dm.or.static_stars)then
           call synchro_fine_static(ilevel)
        else
           call synchro_fine(ilevel)
        end if
     end if

     if(hydro)then
                               call timer('poisson','start')

        ! Add gravity source term with half time step and new force
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

        ! Update boundaries
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
        if(simple_boundary)call make_boundary_hydro(ilevel)

        ! Compute Bondi-Hoyle accretion parameters
#if NDIM==3
                               call timer('sinks','start')
        if(sink)call collect_acczone_avg(ilevel)
#endif
     end if
  end if

#ifdef RT
  ! Turn on RT in case of rt_stars and first stars just created:
  ! Update photon packages according to star particles
                               call timer('SEDs','start')
  if(rt .and. rt_star) call update_star_RT_feedback(ilevel)
#endif

  !----------------------
  ! Compute new time step
  !----------------------
                               call timer('courant','start')
  call newdt_fine(ilevel)
#ifdef SOLVERmhd
  if(dt_max_user>0.0d0.and.ilevel==nlevelmax)dtnew(ilevel)=MIN(dtnew(ilevel),dt_max_user)
#endif
  if(ilevel>levelmin)then
     dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if

  ! Set unew equal to uold
                               call timer('hydro - set unew','start')
  if(hydro)call set_unew(ilevel)

#ifdef RT
  ! Set rtunew equal to rtuold
                               call timer('rt unew','start')
  if(rt)call rt_set_unew(ilevel)
#endif

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<nlevelmax)then
     if(numbtot(1,ilevel+1)>0)then
        if(nsubcycle(ilevel)==2)then
           call amr_step(ilevel+1,1)
           call amr_step(ilevel+1,2)
        else
           call amr_step(ilevel+1,1)
        endif
     else
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)
#if NDIM==3
        if(sink)call update_sink(ilevel)
#endif
     end if
  else
     call update_time(ilevel)
#if NDIM==3
     if(sink)call update_sink(ilevel)
#endif
  end if

  ! Thermal feedback from stars
#if NDIM==3
  if(hydro.and.star.and.eta_sn>0.and. &
       (trim(feedback_model).eq.'thermal' .or. &
        trim(feedback_model).eq.'agertz'   .or. &
        trim(feedback_model).eq.'kinetic')) then
                               call timer('feedback','start')
               call thermal_feedback(ilevel)
  endif
#endif

  ! Density threshold or Bondi accretion onto sink particle
#if NDIM==3
  if(sink)then
                               call timer('sinks','start')
     call grow_sink(ilevel,.false.)
  end if
#endif
  !-----------
  ! Hydro step
  !-----------
  if((hydro).and.(.not.static_gas))then

     ! Hyperbolic solver
                               call timer('hydro - godunov','start')
     call godunov_fine(ilevel)

     ! Reverse update boundaries
                               call timer('hydro - rev ghostzones','start')
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_reverse_dp(unew(1,ivar),ilevel)
#ifdef SOLVERmhd
     end do
#else
     end do
#endif
     if(momentum_feedback>0)then
        call make_virtual_reverse_dp(pstarnew(1),ilevel)
     endif
     if(pressure_fix)then
        call make_virtual_reverse_dp(enew(1),ilevel)
        call make_virtual_reverse_dp(divu(1),ilevel)
     endif

     ! Set uold equal to unew
                               call timer('hydro - set uold','start')
     call set_uold(ilevel)

     ! Add gravity source term with half time step and old force
     ! in order to complete the time step
                               call timer('poisson','start')
     if(poisson)call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

     ! Restriction operator
                               call timer('hydro upload fine','start')
     call upload_fine(ilevel)

  endif

  !---------------------
  ! Do RT/Chemistry step
  !---------------------
#ifdef RT

#ifdef NTRACEGROUPS
  if(rt) then     
     !if (myid.eq.1) write(*,*) 'a',minind,  halopropz(minind),(1.d0/aexp)-1
     if (rt_tracemassbins(1).gt.-1.d0) then 
        !Find which halo catalog we are looking at
        !If the halo catalog has changed, we need to update ptracegroup for all star particles
        !we also have to load in the new halo catalog
        adiff = 1.d6 !arbirtraily large number
        minind_temp = 1
        do jj=1,size(halopropz)
           adiffnew = abs(((1.d0/aexp)-1.d0)-halopropz(jj))
           if (adiffnew.lt.adiff) then
              adiff = adiffnew
              minind_temp = jj
           endif
        enddo
        if (minind_temp.ne.minind) then
           minind = minind_temp
           !Load in new catalog
           write(hhprops,'(a,a,I0.3,a)') trim(initfile(levelmin)),"/halo_files/ramses_halo_properties_",minind_temp-1,".dat"
           if(myid==1) write(*,*) 'Loading in ', hhprops
           inquire(file=hhprops, exist=okhprop)
           if(.not. okhprop ) then
              if(myid.eq.1) then 
                 write(*,*) 'Cannot read halo properties files...'
                 write(*,*) 'Check if SED-directory contains properties files '
              endif
              call clean_stop
           end if
           open(unit=10,file=hhprops,status='old',action='read')
           read(10,*) hnhs
           deallocate(haloprops)
           allocate(haloprops(hnhs,5))
           do ih = 1,hnhs
              read(10,*) haloprops(ih,:)
           enddo
           close(10)

           file_trigger=.true.
        else
           file_trigger=.false.
        endif
     endif
  endif 
#endif

  if(rt .and. rt_advect) then  
                               call timer('radiative transfer','start')
     call rt_step(ilevel)
  else
     ! Still need a chemistry call if RT is defined but not
     ! actually doing radiative transfer (i.e. rt==false):
                               call timer('cooling','start')
     if(hydro .and. (neq_chem.or.cooling.or.T2_star>0.0))call cooling_fine(ilevel)
  endif
  ! Regular updates and book-keeping:
  if(ilevel==levelmin) then
                               call timer('rt bookkeeping','start')
     if(cosmo) call update_rt_c
     if(cosmo .and. haardt_madau) call update_UVrates(dble(aexp))
     if(cosmo .and. rt_isDiffuseUVsrc) call update_UVsrc
                               call timer('cooling tables','start')
     if(cosmo) call update_coolrates_tables(dble(aexp))
                               call timer('rt stats','start')
     if(ilevel==levelmin) call output_rt_stats
  endif
#else
                               call timer('cooling','start')
  if((hydro).and.(.not.static_gas)) then
    if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
  endif
#endif

  !---------------
  ! Move particles
  !---------------
  if(pic)then
                               call timer('particles','start')
     if(static_dm.or.static_stars)then
        call move_fine_static(ilevel) ! Only remaining particles
     else
        call move_fine(ilevel) ! Only remaining particles
     end if
  end if

  !----------------------------------
  ! Star formation in leaf cells only
  !----------------------------------
#if NDIM==3
                               call timer('star formation','start')
  if(hydro.and.star.and.(.not.static_gas))call star_formation(ilevel)
#endif
  !---------------------------------------
  ! Update physical and virtual boundaries
  !---------------------------------------
  if((hydro).and.(.not.static_gas))then
                               call timer('hydro - ghostzones','start')
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
     if(momentum_feedback>0)call make_virtual_fine_dp(pstarold(1),ilevel)
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif

#ifdef SOLVERmhd
  ! Magnetic diffusion step
  if((hydro).and.(.not.static_gas))then
     if(eta_mag>0d0.and.ilevel==levelmin)then
                               call timer('hydro - diffusion','start')
        call diffusion
     endif

     ! Heat transport step
     if(conduction)then
                               call timer('Electronic conduction','start')
        itemp=1;call conduction_cg(ilevel,icount,itemp)
        if(conduction_ion)then
                               call timer('Ion        conduction','start')
           itemp=2;call conduction_cg(ilevel,itemp)
        endif
     end if
     ! Two temperatures coupling (physical and virtual boundaries are updated within coupling_fine)
     if(twotemp.and.coupling_out_of_conduction)call coupling_fine(ilevel)
#if NCR>0
     ! Cosmic ray diffusion step
     if(cr_diffusion)then
        if(fecr.gt.0. .and. nstar_tot.gt.0) then
                               call timer('Cosmic rays diffusion','start')
           call crdiffusion_cg(ilevel)
        else
           if (myid.eq.1 .and. ilevel.eq.levelmin) &
                write(*,*) 'No stars so no call to crdiffusion_cg'
        end if
     end if
#endif
     
     ! Update boundaries
     ! JOKI: this seems like a repeat of the just-above calls to make_virtual_fine_dp
     !       We could perhaps remove the calls above to speed things up?
     do ivar=1,nvar+3
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  end if
#endif

  !-----------------------
  ! Compute refinement map
  !-----------------------
                               call timer('flag','start')
  if(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)) call flag_fine(ilevel,icount)

  !----------------------------
  ! Merge finer level particles
  !----------------------------
                               call timer('particles','start')
  if(pic)call merge_tree_fine(ilevel)

  !---------------
  ! Radiation step
  !---------------
#ifdef ATON
  if(aton.and.ilevel==levelmin)then
                               call timer('aton','start')
     call rad_step(dtnew(ilevel))
  endif
#endif

  if(sink)then
                               call timer('sinks','start')
     !-------------------------------
     ! Update coarser level sink velocity
     !-------------------------------
     if(ilevel>levelmin)then
        vsold(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel-1)
        if(nsubcycle(ilevel-1)==1)vsnew(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel)
        if(icount==2)vsnew(1:nsink,1:ndim,ilevel-1)= &
             (vsold(1:nsink,1:ndim,ilevel)*dtold(ilevel)+vsnew(1:nsink,1:ndim,ilevel)*dtnew(ilevel))/ &
             (dtold(ilevel)+dtnew(ilevel))
     end if
     !---------------
     ! Sink production
     !---------------
#if NDIM==3
     if(ilevel==levelmin)call create_sink
#endif
  end if

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step

!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################

#ifdef RT
subroutine rt_step(ilevel)
  use amr_parameters, only: dp
  use amr_commons,    only: levelmin, t, dtnew, myid
  use rt_hydro_commons
  use UV_module
  use SED_module,     only: star_RT_feedback
  use mpi_mod
  implicit none
  integer, intent(in) :: ilevel

!--------------------------------------------------------------------------
!  Radiative transfer and chemistry step. Either do one step on ilevel,
!  with radiation field updates in coarser level neighbours, or, if
!  rt_nsubsteps>1, do many substeps in ilevel only, using Dirichlet
!  boundary conditions for the level boundaries.
!--------------------------------------------------------------------------

  real(dp) :: dt_hydro, t_left, dt_rt, t_save
  integer  :: i_substep, ivar

  dt_hydro = dtnew(ilevel)                   ! Store hydro timestep length
  t_left = dt_hydro
  ! We shift the time backwards one hydro-dt, to get evolution of stellar
  ! ages within the hydro timestep, in the case of rt subcycling:
  t_save=t ; t=t-t_left

  i_substep = 0
  do while (t_left > 0)                      !                RT sub-cycle
     i_substep = i_substep + 1
     call get_rt_courant_coarse(dt_rt,ilevel)
     ! Temporarily change timestep length to rt step:
     dtnew(ilevel) = MIN(t_left, dt_rt/2**(ilevel-levelmin))
     t = t + dtnew(ilevel) ! Shift the time forwards one dt_rt

     if (i_substep > 1) call rt_set_unew(ilevel)

                               call timer('stellar radiation','start')
     if(rt_star) call star_RT_feedback(ilevel,dble(dtnew(ilevel)))
                               call timer('radiative transfer','start')

#ifdef NTRACEGROUPS
     !Set file trigger to false so we are not looping over all star particles
     !on every substep.  We only need to do this on the first substep if we 
     !need to move to a new halo catalog
     file_trigger=.false.
#endif

     ! Hyperbolic solver
     if(rt_advect) call rt_godunov_fine(ilevel,dtnew(ilevel))

     call add_rt_sources(ilevel,dtnew(ilevel))

     ! Reverse update boundaries
#ifdef NTRACEGROUPS
     call multi_d_rtunew1 (ilevel)
#endif
     do ivar=1,nrtvartrace
        call make_virtual_reverse_rtdp(rtunew(1,ivar),ilevel)
     end do
#ifdef NTRACEGROUPS
     call div_d_rtunew1 (ilevel)
#endif

     ! Set rtuold equal to rtunew
     call rt_set_uold(ilevel)

                               call timer('cooling','start')
     if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
                               call timer('radiative transfer','start')
     
     do ivar=1,nrtvartrace
        call make_virtual_fine_rtdp(rtuold(1,ivar),ilevel)
     end do
     if(simple_boundary)call rt_make_boundary_hydro(ilevel)

     t_left = t_left - dtnew(ilevel)
  end do                                   !          End RT subcycle loop
  dtnew(ilevel) = dt_hydro                 ! Restore hydro timestep length
  t = t_save       ! Restore original time (otherwise tiny roundoff error)

  ! Restriction operator to update coarser level split cells
  call rt_upload_fine(ilevel)

  if (myid==1 .and. rt_nsubcycle .gt. 1) write(*,901) ilevel, i_substep

901 format (' Performed level', I3, ' RT-step with ', I5, ' subcycles')

end subroutine rt_step

#ifdef NTRACEGROUPS
! Prep tracer groups in rtunew array for reverse update
subroutine multi_d_rtunew (ilevel)
  use amr_commons
  use hydro_commons,only:nvar
  use rt_parameters
  use rt_hydro_commons,only:rtunew
  implicit none
  integer,intent(in)::ilevel
  integer,save::ncache,igrid,ngrid,i,ind_grid,ind,iskip,ind_cell,pgroup,ivar,sizar
  
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid = active(ilevel)%igrid(igrid+i-1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           ind_cell = iskip + ind_grid
           if(son(ind_cell)==0)then
              pgroup = 0
              do ivar=nrtvar+1,nrtvartrace
                 if (ivar-nrtvar.gt.pgroup*ntraceGroups) then
                    pgroup = pgroup + 1
                 end if
                 sizar = iGroups(pgroup)
                 rtunew(ind_cell,ivar) = rtunew(ind_cell,ivar)*max(rtunew(ind_cell,sizar),smallNp)
             end do
           endif
        end do ! ind=1,8
     end do ! i=1,ngrid
  end do
end subroutine multi_d_rtunew

!Fix tracer groups in rtunew array after reverse update
subroutine div_d_rtunew (ilevel)
  use amr_commons
  use hydro_commons,only:nvar
  use rt_parameters
  use rt_hydro_commons,only:rtunew
  implicit none
  integer,intent(in)::ilevel
  integer,save::ncache,igrid,ngrid,i,ind_grid,ind,iskip,ind_cell,ivar,pgroup,sizar,ih,ihh
  real(dp)::totfrac

  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid = active(ilevel)%igrid(igrid+i-1)
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           ind_cell = iskip + ind_grid
           if(son(ind_cell)==0)then
              pgroup = 0
              do ivar=nrtvar+1,nrtvartrace
                 if (ivar-nrtvar.gt.pgroup*ntraceGroups) then
                    pgroup = pgroup + 1
                 end if
                 sizar = iGroups(pgroup)
                 if (rtunew(ind_cell,sizar).gt.0.0d0) then
                    rtunew(ind_cell,ivar) = rtunew(ind_cell,ivar)/rtunew(ind_cell,sizar)
                 else 
                    rtunew(ind_cell,ivar) = 1.d0/ntracegroups
                 endif              
             end do
           endif
        end do ! ind=1,8
     end do ! i=1,ngrid
  end do
end subroutine div_d_rtunew

SUBROUTINE multi_d_rtunew1(ilevel)
!------------------------------------------------------------------------
  use amr_commons
  use rt_hydro_commons
  implicit none
  integer::ilevel,pgroup,sizar
  integer::i,ivar,ind,icpu,iskip
!------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Convert photon tracer fractions to densities in real and virtual cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     pgroup = 0
     do ivar=nrtvar+1,nrtvartrace
        if (ivar-nrtvar.gt.pgroup*ntraceGroups) then
           pgroup = pgroup + 1
        end if
        do i=1,active(ilevel)%ngrid
           rtunew(active(ilevel)%igrid(i)+iskip,ivar) = rtunew(active(ilevel)%igrid(i)+iskip,ivar)*rtunew(active(ilevel)%igrid(i)+iskip,iGroups(pgroup))
        end do
     end do
  end do

  ! Set rtunew to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     pgroup = 0
     do ivar=nrtvar+1,nrtvartrace
        if (ivar-nrtvar.gt.pgroup*ntraceGroups) then
           pgroup = pgroup + 1
        end if
        do i=1,reception(icpu,ilevel)%ngrid
           rtunew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=rtunew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)*rtunew(reception(icpu,ilevel)%igrid(i)+iskip,iGroups(pgroup))
        end do
     end do
  end do
  end do

111 format('   Entering rt_set_unew for level ',i2)

END SUBROUTINE multi_d_rtunew1

SUBROUTINE div_d_rtunew1(ilevel)
!------------------------------------------------------------------------
  use amr_commons
  use rt_hydro_commons
  implicit none
  integer::ilevel,pgroup,sizar
  integer::i,ivar,ind,icpu,iskip
!------------------------------------------------------------------------
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Convert photon tracer fractions from densities to fractions in real cells only
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     pgroup = 0
     do ivar=nrtvar+1,nrtvartrace
        if (ivar-nrtvar.gt.pgroup*ntraceGroups) then
           pgroup = pgroup + 1
        end if
        do i=1,active(ilevel)%ngrid
           if (rtunew(active(ilevel)%igrid(i)+iskip,iGroups(pgroup)).gt. smallNp) then
              !if (rtunew(active(ilevel)%igrid(i)+iskip,iGroups(pgroup)).gt.smallNp) then
              rtunew(active(ilevel)%igrid(i)+iskip,ivar) = rtunew(active(ilevel)%igrid(i)+iskip,ivar)/rtunew(active(ilevel)%igrid(i)+iskip,iGroups(pgroup))
              !else
              !   rtunew(active(ilevel)%igrid(i)+iskip,ivar) = 1.d0/ntracegroups
              !endif
           else
              rtunew(active(ilevel)%igrid(i)+iskip,ivar) = 1.d0/ntracegroups
           endif
        end do
     end do
  end do

111 format('   Entering rt_set_unew for level ',i2)

END SUBROUTINE div_d_rtunew1

#endif

#endif
