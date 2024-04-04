! ---------------------------------------------------------------
!  COND_SPLIT  This routine solves the heat flux following the  
!              anistropic heat conduction.
!
!  inputs/outputs
!  uin         => (const)  input state
!  gravin      => (const)  input gravitational acceleration
!  iu1,iu2     => (const)  first and last index of input array,
!  ju1,ju2     => (const)  cell centered,    
!  ku1,ku2     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx          => (const)  cell width
!  dt          => (const)  time step
!  ngrid       => (const)  number of sub-grids
!  ndim        => (const)  number of dimensions
!
!  uin = (\rho, \rho u, \rho v, \rho w, Etot, A, B, C)
!  the hydro variable are cell-centered
!  whereas the magnetic field B=(A,B,C) are face-centered.
!  Note that here we have 3 components for v and B whatever ndim.
!
!  This routine was written by Yohan Dubois & Benoit Commerçon
! ----------------------------------------------------------------
subroutine cond_split(uin,flux,dx,dt,ngrid,compute,fdx,itemp)
  use amr_parameters
  use const             
  use hydro_parameters
  use constants,ONLY:mH
  implicit none 

  integer,intent(IN) ::ngrid,compute,itemp
  real(dp),intent(IN)::dx,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),intent(IN)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(IN)::fdx

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),intent(OUT)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::Xflux,Yflux,Zflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::kspitzer
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Temp,Temp2,dens
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ffdx

  ! Local scalar variables
  integer::i,j,k,l
  integer::ilo,ihi,jlo,jhi,klo,khi

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_kspitzer
  real(dp)::Tp_loc,scale_kappa
  real(dp)::scale_n

  real(dp)::scale_jsat
  integer ::ind_Teold

  ind_Teold = 3 ! Stored in a different index than in other routines (nvar+3)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kspitzer=scale_d*scale_l**4/scale_t**3/scale_T2
  scale_jsat=scale_d*scale_v**2/scale_t
  scale_kappa=scale_d*scale_l*scale_v**3
  fudge_mfp=262229.17/coulomb_log/scale_l ! 3^1.5*k_B^2/(4*pi^0.5*e^4*coulomb_log)
  if(itemp==1)then
     scale_n=scale_d/mH/mu_electron
  else
     scale_n=scale_d/mH/mu_ion
  endif
!!$  if(myid==1)then
!!$     write(*,'(A,e10.3)')'fudge_mfp',fudge_mfp
!!$     write(*,'(A,e10.3)')'scale_n  ',scale_n
!!$  endif
  
  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  do k = ku1, ku2
  do j = ju1, ju2
  do i = iu1, iu2
     do l = 1, ngrid

        Tp_loc = uin(l,i,j,k,ind_Teold)

        dens(l,i,j,k)=uin(l,i,j,k,1)*scale_n
        
        kspitzer(l,i,j,k)=uin(l,i,j,k,4)
        bf(l,i,j,k,1)=uin(l,i,j,k,6)
        bf(l,i,j,k,2)=uin(l,i,j,k,7)
        bf(l,i,j,k,3)=uin(l,i,j,k,8)
        
        ffdx(l,i,j,k)=fdx(l,i,j,k)

        if(compute==1)Temp(l,i,j,k)=Tp_loc
        if(compute==2)Temp(l,i,j,k)=uin(l,i,j,k,2)

        Temp2(l,i,j,k)=Tp_loc
     end do
  enddo
  enddo
  enddo

  ! Compute the heat flux in X direction  
  call cmpXheatflx(Temp,dens,bf,kspitzer,Xflux,dx,dt,ngrid,compute,ffdx,Temp2)
  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        flux(l,i,j,k,5,1)=Xflux(l,i,j,k)
    enddo
  enddo
  enddo
  enddo

#if NDIM>1
  ! Compute the heat flux in Y direction
  call cmpYheatflx(Temp,dens,bf,kspitzer,Yflux,dx,dt,ngrid,compute,Temp2)
  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
        flux(l,i,j,k,5,2)=Yflux(l,i,j,k)
     enddo
  enddo
  enddo
  enddo
#endif
#if NDIM>2
  ! Compute the heat flux in Z direction
  call cmpZheatflx(Temp,dens,bf,kspitzer,Zflux,dx,dt,ngrid,compute,Temp2)
  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do l = 1, ngrid
        flux(l,i,j,k,5,3)=Zflux(l,i,j,k)
     enddo
  enddo
  enddo
  enddo
#endif
end subroutine cond_split
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmpXheatflx(Temp,dens,bf,kspitzer,myflux,dx,dt,ngrid,compute,ffdx,Temp2)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer,intent(IN) ::ngrid,compute
  real(dp),intent(IN)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2),intent(OUT)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3),intent(IN)::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(IN)::kspitzer
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(IN)::Temp,Temp2,dens
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(IN)::ffdx

  real(dp)::fx,oneovertwodx,oneoverfourdx 
  real(dp)::Bnorm2_1,Bnorm2_2,Bnorm2_3,Bnorm2_4
  real(dp)::dTdx1,dTdx2,dTdx3,dTdx4
  real(dp)::dTdy1,dTdy2,dTdy3,dTdy4
  real(dp)::dTdz1,dTdz2,dTdz3,dTdz4
  real(dp)::bx1,bx2,bx3,bx4
  real(dp)::by1,by2,by3,by4
  real(dp)::bz1,bz2,bz3,bz4
  real(dp)::fx1,fx2,fx3,fx4
  real(dp)::kparx1,kparx2,kparx3,kparx4,oneminuskperp
  real(dp)::dmean1,dmean2,dmean3,dmean4
  real(dp)::Tmean1,Tmean2,Tmean3,Tmean4
  real(dp)::dTolddx1,dTolddx2,dTolddx3,dTolddx4
  real(dp)::dTolddy1,dTolddy2,dTolddy3,dTolddy4
  real(dp)::dTolddz1,dTolddz2,dTolddz3,dTolddz4
  real(dp)::sat_coef1,sat_coef2,sat_coef3,sat_coef4
  real(dp)::ratio1,ratio2,ratio3,ratio4
  real(dp)::famr
#if NDIM==2
  real(dp)::Bnorm
#endif
  ! Local scalar variables
  integer::i,j,k,l
  integer::jlo,jhi,klo,khi

  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  oneovertwodx =0.50d0/dx
  oneoverfourdx=0.25d0/dx
  oneminuskperp=1.0d0-k_perp

if(isotrope_cond)then
  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        dTdx1  =(Temp(l,i,j,k)-Temp(l,i-1,j,k))/dx
        kparx1=0.5d0*(kspitzer(l,i,j,k)+kspitzer(l,i-1,j,k))
        if(saturation)then
           dmean1=0.5d0*(dens (l,i,j,k)+dens (l,i-1,j,k))
           Tmean1=0.5d0*(Temp2(l,i,j,k)+Temp2(l,i-1,j,k))
           dTolddx1=abs(Temp2(l,i,j,k)-Temp2(l,i-1,j,k))/dx
           ratio1=fudge_mfp*Tmean1*dTolddx1/dmean1
           sat_coef1=1d0/(1d0+4.2d0*ratio1)
           kparx1=sat_coef1*kparx1
        endif
        if(compute .ne. 3)then
           fx=kparx1*dTdx1
        else
           fx=kparx1
        endif
!!$        famr=max(ffdx(l,i,j,k),ffdx(l,i-1,j,k))
!!$        fx=fx/famr
        myflux(l,i,j,k)=fx*dt/dx
!!$        if(compute==3)myflux(l,i,j,k)=kparx1*dt/dx**2
     enddo
  enddo
  enddo
  enddo

else

  if (.not.slopelim_cond)then
  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        famr=max(ffdx(l,i,j,k),ffdx(l,i-1,j,k))
#if NDIM==1
!!$        1-------

        if(compute .ne. 3) dTdx1=(Temp(l,i,j,k)-Temp(l,i-1,j,k))/dx
        kparx1=0.5d0*(kspitzer(l,i,j,k)+kspitzer(l,i-1,j,k)) ! arithmetic mean
        if(saturation)then
           dmean1=0.5d0*(dens (l,i,j,k)+dens (l,i-1,j,k))
           Tmean1=0.5d0*(Temp2(l,i,j,k)+Temp2(l,i-1,j,k))
           dTolddx1=abs(Temp2(l,i,j,k)-Temp2(l,i-1,j,k))/dx
           ratio1=fudge_mfp*Tmean1*dTolddx1/dmean1
           sat_coef1=1d0/(1d0+4.2d0*ratio1)
           kparx1=sat_coef1*kparx1
        endif
        if(compute .ne. 3)then
           fx=kparx1*dTdx1
        else 
           fx=kparx1
        endif
        fx=fx/famr

#endif
#if NDIM==2
!!$        2-------
!!$        |      |
!!$        |      |
!!$        |      |
!!$        1-------
        
        if(compute .ne. 3)then   
           dTdx1=(Temp(l,i  ,j  ,k)+Temp(l,i  ,j-1,k) &
                - Temp(l,i-1,j  ,k)-Temp(l,i-1,j-1,k))*oneovertwodx
           dTdx2=(Temp(l,i  ,j+1,k)+Temp(l,i  ,j  ,k) &
                - Temp(l,i-1,j+1,k)-Temp(l,i-1,j  ,k))*oneovertwodx        
           
           dTdy1=(Temp(l,i  ,j  ,k)+Temp(l,i-1,j  ,k) &
                - Temp(l,i  ,j-1,k)-Temp(l,i-1,j-1,k))*oneovertwodx
           dTdy2=(Temp(l,i  ,j+1,k)+Temp(l,i-1,j+1,k) &
                - Temp(l,i  ,j  ,k)-Temp(l,i-1,j  ,k))*oneovertwodx
        endif
        
        bx1=0.5d0*(bf(l,i  ,j-1,k,1)+bf(l,i  ,j  ,k,1))
        bx2=0.5d0*(bf(l,i  ,j  ,k,1)+bf(l,i  ,j+1,k,1))
        
        by1=0.5d0*(bf(l,i-1,j  ,k,2)+bf(l,i  ,j  ,k,2))
        by2=0.5d0*(bf(l,i-1,j+1,k,2)+bf(l,i  ,j+1,k,2))
        
        Bnorm=sqrt(bx1*bx1+by1*by1)
        if(Bnorm.gt.0.0)then
           bx1=bx1/Bnorm
           by1=by1/Bnorm
        endif
        Bnorm=sqrt(bx2*bx2+by2*by2)
        if(Bnorm.gt.0.0)then
           bx2=bx2/Bnorm
           by2=by2/Bnorm
        endif

        kparx1=0.25d0*(kspitzer(l,i  ,j  ,k)+kspitzer(l,i  ,j-1,k) &
             +      kspitzer(l,i-1,j  ,k)+kspitzer(l,i-1,j-1,k))
        kparx2=0.25d0*(kspitzer(l,i  ,j+1,k)+kspitzer(l,i  ,j  ,k) &
             +      kspitzer(l,i-1,j+1,k)+kspitzer(l,i-1,j  ,k))

        if(saturation)then
           dmean1=0.25d0*(dens (l,i  ,j  ,k)+dens (l,i  ,j-1,k) &
                +         dens (l,i-1,j  ,k)+dens (l,i-1,j-1,k))
           dmean2=0.25d0*(dens (l,i  ,j+1,k)+dens (l,i  ,j  ,k) &
                +         dens (l,i-1,j+1,k)+dens (l,i-1,j  ,k))
           Tmean1=0.25d0*(Temp2(l,i  ,j  ,k)+Temp2(l,i  ,j-1,k) &
                +         Temp2(l,i-1,j  ,k)+Temp2(l,i-1,j-1,k))
           Tmean2=0.25d0*(Temp2(l,i  ,j+1,k)+Temp2(l,i  ,j  ,k) &
                +         Temp2(l,i-1,j+1,k)+Temp2(l,i-1,j  ,k))
           dTolddx1=(Temp2(l,i  ,j  ,k)+Temp2(l,i  ,j-1,k) &
                       - Temp2(l,i-1,j  ,k)-Temp2(l,i-1,j-1,k))*oneovertwodx
           dTolddx2=(Temp2(l,i  ,j+1,k)+Temp2(l,i  ,j  ,k) &
                       - Temp2(l,i-1,j+1,k)-Temp2(l,i-1,j  ,k))*oneovertwodx                   
           dTolddy1=(Temp2(l,i  ,j  ,k)+Temp2(l,i-1,j  ,k) &
                       - Temp2(l,i  ,j-1,k)-Temp2(l,i-1,j-1,k))*oneovertwodx
           dTolddy2=(Temp2(l,i  ,j+1,k)+Temp2(l,i-1,j+1,k) &
                       - Temp2(l,i  ,j  ,k)-Temp2(l,i-1,j  ,k))*oneovertwodx
           ratio1=fudge_mfp*Tmean1*dsqrt(MAX(dTolddx1**2+dTolddy1**2,0d0))/dmean1
           ratio2=fudge_mfp*Tmean2*dsqrt(MAX(dTolddx2**2+dTolddy2**2,0d0))/dmean2
           sat_coef1=1d0/(1d0+4.2d0*ratio1)
           sat_coef2=1d0/(1d0+4.2d0*ratio2)
           kparx1=sat_coef1*kparx1
           kparx2=sat_coef2*kparx2
        endif

        if(compute .ne. 3)then   
           fx1=kparx1*(bx1*oneminuskperp*(bx1*dTdx1+by1*dTdy1)+k_perp*dTdx1)
           fx2=kparx2*(bx2*oneminuskperp*(bx2*dTdx2+by2*dTdy2)+k_perp*dTdx2)
        else ! Preconditionner
           fx1=kparx1*(bx1*oneminuskperp*(bx1+by1)+k_perp)
           fx2=kparx2*(bx2*oneminuskperp*(bx2+by2)+k_perp)
        end if
        fx=0.5d0*(fx1+fx2)
#endif
#if NDIM==3
!!$          4--------
!!$         / |      /|
!!$        /  |     / |
!!$        3-------   |
!!$        |  |   |   |
!!$        | /2   |  /
!!$        |/     | /
!!$        1-------
        
        ! Centered symmetric scheme
        if(compute .ne. 3)then   
           dTdx1=(Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  )+Temp(l,i  ,j  ,k-1)+Temp(l,i  ,j-1,k-1) &
                - Temp(l,i-1,j  ,k  )-Temp(l,i-1,j-1,k  )-Temp(l,i-1,j  ,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdx2=(Temp(l,i  ,j+1,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i  ,j+1,k-1)+Temp(l,i  ,j  ,k-1) &
                - Temp(l,i-1,j+1,k  )-Temp(l,i-1,j  ,k  )-Temp(l,i-1,j+1,k-1)-Temp(l,i-1,j  ,k-1))*oneoverfourdx
           dTdx3=(Temp(l,i  ,j  ,k+1)+Temp(l,i  ,j-1,k+1)+Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  ) &
                - Temp(l,i-1,j  ,k+1)-Temp(l,i-1,j-1,k+1)-Temp(l,i-1,j  ,k  )-Temp(l,i-1,j-1,k  ))*oneoverfourdx
           dTdx4=(Temp(l,i  ,j+1,k+1)+Temp(l,i  ,j  ,k+1)+Temp(l,i  ,j+1,k  )+Temp(l,i  ,j  ,k  ) &
                - Temp(l,i-1,j+1,k+1)-Temp(l,i-1,j  ,k+1)-Temp(l,i-1,j+1,k  )-Temp(l,i-1,j  ,k  ))*oneoverfourdx

           dTdy1=(Temp(l,i  ,j  ,k  )+Temp(l,i-1,j  ,k  )+Temp(l,i  ,j  ,k-1)+Temp(l,i-1,j  ,k-1) &
                - Temp(l,i  ,j-1,k  )-Temp(l,i-1,j-1,k  )-Temp(l,i  ,j-1,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdy2=(Temp(l,i  ,j+1,k  )+Temp(l,i-1,j+1,k  )+Temp(l,i  ,j+1,k-1)+Temp(l,i-1,j+1,k-1) &
                - Temp(l,i  ,j  ,k  )-Temp(l,i-1,j  ,k  )-Temp(l,i  ,j  ,k-1)-Temp(l,i-1,j  ,k-1))*oneoverfourdx
           dTdy3=(Temp(l,i  ,j  ,k+1)+Temp(l,i-1,j  ,k+1)+Temp(l,i  ,j  ,k  )+Temp(l,i-1,j  ,k  ) &
                - Temp(l,i  ,j-1,k+1)-Temp(l,i-1,j-1,k+1)-Temp(l,i  ,j-1,k  )-Temp(l,i-1,j-1,k  ))*oneoverfourdx
           dTdy4=(Temp(l,i  ,j+1,k+1)+Temp(l,i-1,j+1,k+1)+Temp(l,i  ,j+1,k  )+Temp(l,i-1,j+1,k  ) &
                - Temp(l,i  ,j  ,k+1)-Temp(l,i-1,j  ,k+1)-Temp(l,i  ,j  ,k  )-Temp(l,i-1,j  ,k  ))*oneoverfourdx
           
           dTdz1=(Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  )+Temp(l,i-1,j  ,k  )+Temp(l,i-1,j-1,k  ) &
                - Temp(l,i  ,j  ,k-1)-Temp(l,i  ,j-1,k-1)-Temp(l,i-1,j  ,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdz2=(Temp(l,i  ,j+1,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i-1,j+1,k  )+Temp(l,i-1,j  ,k  ) &
                - Temp(l,i  ,j+1,k-1)-Temp(l,i  ,j  ,k-1)-Temp(l,i-1,j+1,k-1)-Temp(l,i-1,j  ,k-1))*oneoverfourdx
           dTdz3=(Temp(l,i  ,j  ,k+1)+Temp(l,i-1,j  ,k+1)+Temp(l,i  ,j-1,k+1)+Temp(l,i-1,j-1,k+1) &
                - Temp(l,i  ,j  ,k  )-Temp(l,i-1,j  ,k  )-Temp(l,i  ,j-1,k  )-Temp(l,i-1,j-1,k  ))*oneoverfourdx
           dTdz4=(Temp(l,i  ,j+1,k+1)+Temp(l,i-1,j+1,k+1)+Temp(l,i  ,j  ,k+1)+Temp(l,i-1,j  ,k+1) &
                - Temp(l,i  ,j+1,k  )-Temp(l,i-1,j+1,k  )-Temp(l,i  ,j  ,k  )-Temp(l,i-1,j  ,k  ))*oneoverfourdx
        end if

        bx1=(bf(l,i  ,j-1,k-1,1)+bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j-1,k  ,1)+bf(l,i  ,j  ,k  ,1))
        bx2=(bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j+1,k-1,1)+bf(l,i  ,j  ,k  ,1)+bf(l,i  ,j+1,k  ,1))
        bx3=(bf(l,i  ,j-1,k  ,1)+bf(l,i  ,j  ,k  ,1)+bf(l,i  ,j-1,k+1,1)+bf(l,i  ,j  ,k+1,1))
        bx4=(bf(l,i  ,j  ,k  ,1)+bf(l,i  ,j+1,k  ,1)+bf(l,i  ,j  ,k+1,1)+bf(l,i  ,j+1,k+1,1))

        by1=(bf(l,i-1,j  ,k-1,2)+bf(l,i  ,j  ,k-1,2)+bf(l,i-1,j  ,k  ,2)+bf(l,i  ,j  ,k  ,2))
        by2=(bf(l,i  ,j+1,k-1,2)+bf(l,i-1,j+1,k-1,2)+bf(l,i  ,j+1,k  ,2)+bf(l,i-1,j+1,k  ,2))
        by3=(bf(l,i-1,j  ,k  ,2)+bf(l,i  ,j  ,k  ,2)+bf(l,i-1,j  ,k+1,2)+bf(l,i  ,j  ,k+1,2))
        by4=(bf(l,i  ,j+1,k  ,2)+bf(l,i-1,j+1,k  ,2)+bf(l,i  ,j+1,k+1,2)+bf(l,i-1,j+1,k+1,2))
        
        bz1=(bf(l,i-1,j-1,k  ,3)+bf(l,i-1,j  ,k  ,3)+bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3))
        bz2=(bf(l,i-1,j  ,k  ,3)+bf(l,i-1,j+1,k  ,3)+bf(l,i  ,j  ,k  ,3)+bf(l,i  ,j+1,k  ,3))
        bz3=(bf(l,i  ,j-1,k+1,3)+bf(l,i  ,j  ,k+1,3)+bf(l,i-1,j-1,k+1,3)+bf(l,i-1,j  ,k+1,3))
        bz4=(bf(l,i  ,j  ,k+1,3)+bf(l,i  ,j+1,k+1,3)+bf(l,i-1,j  ,k+1,3)+bf(l,i-1,j+1,k+1,3))

        Bnorm2_1=(bx1*bx1+by1*by1+bz1*bz1)
        Bnorm2_2=(bx2*bx2+by2*by2+bz2*bz2)
        Bnorm2_3=(bx3*bx3+by3*by3+bz3*bz3)
        Bnorm2_4=(bx4*bx4+by4*by4+bz4*bz4)

        kparx1=0.125d0*(kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j-1,k  )+kspitzer(l,i  ,j  ,k-1) &
             +      kspitzer(l,i  ,j-1,k-1)+kspitzer(l,i-1,j  ,k  )+kspitzer(l,i-1,j-1,k  ) &
             +      kspitzer(l,i-1,j  ,k-1)+kspitzer(l,i-1,j-1,k-1))
        kparx2=0.125d0*(kspitzer(l,i  ,j+1,k  )+kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j+1,k-1) &
             +      kspitzer(l,i  ,j  ,k-1)+kspitzer(l,i-1,j+1,k  )+kspitzer(l,i-1,j  ,k  ) &
             +      kspitzer(l,i-1,j+1,k-1)+kspitzer(l,i-1,j  ,k-1))
        kparx3=0.125d0*(kspitzer(l,i  ,j  ,k+1)+kspitzer(l,i  ,j-1,k+1)+kspitzer(l,i  ,j  ,k  ) &
             +      kspitzer(l,i  ,j-1,k  )+kspitzer(l,i-1,j  ,k+1)+kspitzer(l,i-1,j-1,k+1) &
             +      kspitzer(l,i-1,j  ,k  )+kspitzer(l,i-1,j-1,k  ))
        kparx4=0.125d0*(kspitzer(l,i  ,j+1,k+1)+kspitzer(l,i  ,j  ,k+1)+kspitzer(l,i  ,j+1,k  ) &
             +      kspitzer(l,i  ,j  ,k  )+kspitzer(l,i-1,j+1,k+1)+kspitzer(l,i-1,j  ,k+1) &
             +      kspitzer(l,i-1,j+1,k  )+kspitzer(l,i-1,j  ,k  ))

        if(saturation)then
           dmean1=0.125d0*(dens (l,i  ,j  ,k  )+dens (l,i  ,j-1,k  )+dens (l,i  ,j  ,k-1) &
                +          dens (l,i  ,j-1,k-1)+dens (l,i-1,j  ,k  )+dens (l,i-1,j-1,k  ) &
                +          dens (l,i-1,j  ,k-1)+dens (l,i-1,j-1,k-1))
           dmean2=0.125d0*(dens (l,i  ,j+1,k  )+dens (l,i  ,j  ,k  )+dens (l,i  ,j+1,k-1) &
                +          dens (l,i  ,j  ,k-1)+dens (l,i-1,j+1,k  )+dens (l,i-1,j  ,k  ) &
                +          dens (l,i-1,j+1,k-1)+dens (l,i-1,j  ,k-1))
           dmean3=0.125d0*(dens (l,i  ,j  ,k+1)+dens (l,i  ,j-1,k+1)+dens (l,i  ,j  ,k  ) &
                +          dens (l,i  ,j-1,k  )+dens (l,i-1,j  ,k+1)+dens (l,i-1,j-1,k+1) &
                +          dens (l,i-1,j  ,k  )+dens (l,i-1,j-1,k  ))
           dmean4=0.125d0*(dens (l,i  ,j+1,k+1)+dens (l,i  ,j  ,k+1)+dens (l,i  ,j+1,k  ) &
                +          dens (l,i  ,j  ,k  )+dens (l,i-1,j+1,k+1)+dens (l,i-1,j  ,k+1) &
                +          dens (l,i-1,j+1,k  )+dens (l,i-1,j  ,k  ))
           Tmean1=0.125d0*(Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  )+Temp2(l,i  ,j  ,k-1) &
                +          Temp2(l,i  ,j-1,k-1)+Temp2(l,i-1,j  ,k  )+Temp2(l,i-1,j-1,k  ) &
                +          Temp2(l,i-1,j  ,k-1)+Temp2(l,i-1,j-1,k-1))
           Tmean2=0.125d0*(Temp2(l,i  ,j+1,k  )+Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j+1,k-1) &
                +          Temp2(l,i  ,j  ,k-1)+Temp2(l,i-1,j+1,k  )+Temp2(l,i-1,j  ,k  ) &
                +          Temp2(l,i-1,j+1,k-1)+Temp2(l,i-1,j  ,k-1))
           Tmean3=0.125d0*(Temp2(l,i  ,j  ,k+1)+Temp2(l,i  ,j-1,k+1)+Temp2(l,i  ,j  ,k  ) &
                +          Temp2(l,i  ,j-1,k  )+Temp2(l,i-1,j  ,k+1)+Temp2(l,i-1,j-1,k+1) &
                +          Temp2(l,i-1,j  ,k  )+Temp2(l,i-1,j-1,k  ))
           Tmean4=0.125d0*(Temp2(l,i  ,j+1,k+1)+Temp2(l,i  ,j  ,k+1)+Temp2(l,i  ,j+1,k  ) &
                +          Temp2(l,i  ,j  ,k  )+Temp2(l,i-1,j+1,k+1)+Temp2(l,i-1,j  ,k+1) &
                +          Temp2(l,i-1,j+1,k  )+Temp2(l,i-1,j  ,k  ))
           dTolddx1=(Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  )+Temp2(l,i  ,j  ,k-1)+Temp2(l,i  ,j-1,k-1) &
                - Temp2(l,i-1,j  ,k  )-Temp2(l,i-1,j-1,k  )-Temp2(l,i-1,j  ,k-1)-Temp2(l,i-1,j-1,k-1))*oneoverfourdx
           dTolddx2=(Temp2(l,i  ,j+1,k  )+Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j+1,k-1)+Temp2(l,i  ,j  ,k-1) &
                - Temp2(l,i-1,j+1,k  )-Temp2(l,i-1,j  ,k  )-Temp2(l,i-1,j+1,k-1)-Temp2(l,i-1,j  ,k-1))*oneoverfourdx
           dTolddx3=(Temp2(l,i  ,j  ,k+1)+Temp2(l,i  ,j-1,k+1)+Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  ) &
                - Temp2(l,i-1,j  ,k+1)-Temp2(l,i-1,j-1,k+1)-Temp2(l,i-1,j  ,k  )-Temp2(l,i-1,j-1,k  ))*oneoverfourdx
           dTolddx4=(Temp2(l,i  ,j+1,k+1)+Temp2(l,i  ,j  ,k+1)+Temp2(l,i  ,j+1,k  )+Temp2(l,i  ,j  ,k  ) &
                - Temp2(l,i-1,j+1,k+1)-Temp2(l,i-1,j  ,k+1)-Temp2(l,i-1,j+1,k  )-Temp2(l,i-1,j  ,k  ))*oneoverfourdx
           dTolddy1=(Temp2(l,i  ,j  ,k  )+Temp2(l,i-1,j  ,k  )+Temp2(l,i  ,j  ,k-1)+Temp2(l,i-1,j  ,k-1) &
                - Temp2(l,i  ,j-1,k  )-Temp2(l,i-1,j-1,k  )-Temp2(l,i  ,j-1,k-1)-Temp2(l,i-1,j-1,k-1))*oneoverfourdx
           dTolddy2=(Temp2(l,i  ,j+1,k  )+Temp2(l,i-1,j+1,k  )+Temp2(l,i  ,j+1,k-1)+Temp2(l,i-1,j+1,k-1) &
                - Temp2(l,i  ,j  ,k  )-Temp2(l,i-1,j  ,k  )-Temp2(l,i  ,j  ,k-1)-Temp2(l,i-1,j  ,k-1))*oneoverfourdx
           dTolddy3=(Temp2(l,i  ,j  ,k+1)+Temp2(l,i-1,j  ,k+1)+Temp2(l,i  ,j  ,k  )+Temp2(l,i-1,j  ,k  ) &
                - Temp2(l,i  ,j-1,k+1)-Temp2(l,i-1,j-1,k+1)-Temp2(l,i  ,j-1,k  )-Temp2(l,i-1,j-1,k  ))*oneoverfourdx
           dTolddy4=(Temp2(l,i  ,j+1,k+1)+Temp2(l,i-1,j+1,k+1)+Temp2(l,i  ,j+1,k  )+Temp2(l,i-1,j+1,k  ) &
                - Temp2(l,i  ,j  ,k+1)-Temp2(l,i-1,j  ,k+1)-Temp2(l,i  ,j  ,k  )-Temp2(l,i-1,j  ,k  ))*oneoverfourdx
           dTolddz1=(Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  )+Temp2(l,i-1,j  ,k  )+Temp2(l,i-1,j-1,k  ) &
                - Temp2(l,i  ,j  ,k-1)-Temp2(l,i  ,j-1,k-1)-Temp2(l,i-1,j  ,k-1)-Temp2(l,i-1,j-1,k-1))*oneoverfourdx
           dTolddz2=(Temp2(l,i  ,j+1,k  )+Temp2(l,i  ,j  ,k  )+Temp2(l,i-1,j+1,k  )+Temp2(l,i-1,j  ,k  ) &
                - Temp2(l,i  ,j+1,k-1)-Temp2(l,i  ,j  ,k-1)-Temp2(l,i-1,j+1,k-1)-Temp2(l,i-1,j  ,k-1))*oneoverfourdx
           dTolddz3=(Temp2(l,i  ,j  ,k+1)+Temp2(l,i-1,j  ,k+1)+Temp2(l,i  ,j-1,k+1)+Temp2(l,i-1,j-1,k+1) &
                - Temp2(l,i  ,j  ,k  )-Temp2(l,i-1,j  ,k  )-Temp2(l,i  ,j-1,k  )-Temp2(l,i-1,j-1,k  ))*oneoverfourdx
           dTolddz4=(Temp2(l,i  ,j+1,k+1)+Temp2(l,i-1,j+1,k+1)+Temp2(l,i  ,j  ,k+1)+Temp2(l,i-1,j  ,k+1) &
                - Temp2(l,i  ,j+1,k  )-Temp2(l,i-1,j+1,k  )-Temp2(l,i  ,j  ,k  )-Temp2(l,i-1,j  ,k  ))*oneoverfourdx
           ratio1=fudge_mfp*Tmean1*dsqrt(MAX(dTolddx1**2+dTolddy1**2+dTolddz1**2,0.0d0))/dmean1
           ratio2=fudge_mfp*Tmean2*dsqrt(MAX(dTolddx2**2+dTolddy2**2+dTolddz2**2,0.0d0))/dmean2
           ratio3=fudge_mfp*Tmean3*dsqrt(MAX(dTolddx3**2+dTolddy3**2+dTolddz3**2,0.0d0))/dmean3
           ratio4=fudge_mfp*Tmean4*dsqrt(MAX(dTolddx4**2+dTolddy4**2+dTolddz4**2,0.0d0))/dmean4
           sat_coef1=1d0/(1d0+4.2d0*ratio1)
           sat_coef2=1d0/(1d0+4.2d0*ratio2)
           sat_coef3=1d0/(1d0+4.2d0*ratio3)
           sat_coef4=1d0/(1d0+4.2d0*ratio4)
           kparx1=sat_coef1*kparx1
           kparx2=sat_coef2*kparx2
           kparx3=sat_coef3*kparx3
           kparx4=sat_coef4*kparx4
        endif

        if(compute .ne. 3)then   
           fx1=kparx1*(bx1*oneminuskperp*(bx1*dTdx1+by1*dTdy1+bz1*dTdz1)/Bnorm2_1+k_perp*dTdx1)
           fx2=kparx2*(bx2*oneminuskperp*(bx2*dTdx2+by2*dTdy2+bz2*dTdz2)/Bnorm2_2+k_perp*dTdx2)
           fx3=kparx3*(bx3*oneminuskperp*(bx3*dTdx3+by3*dTdy3+bz3*dTdz3)/Bnorm2_3+k_perp*dTdx3)
           fx4=kparx4*(bx4*oneminuskperp*(bx4*dTdx4+by4*dTdy4+bz4*dTdz4)/Bnorm2_4+k_perp*dTdx4)
        else ! Preconditionner
           fx1=kparx1*(bx1*oneminuskperp*(bx1+by1+bz1)/Bnorm2_1+k_perp)
           fx2=kparx2*(bx2*oneminuskperp*(bx2+by2+bz2)/Bnorm2_2+k_perp)
           fx3=kparx3*(bx3*oneminuskperp*(bx3+by3+bz3)/Bnorm2_3+k_perp)
           fx4=kparx4*(bx4*oneminuskperp*(bx4+by4+bz4)/Bnorm2_4+k_perp)
        end if
        fx=0.25d0*(fx1+fx2+fx3+fx4)
#endif
        myflux(l,i,j,k)=fx*dt/dx
     enddo
  enddo
  enddo
  enddo
endif

  if (slopelim_cond)then
     ! TODO
  endif
endif

end subroutine cmpXheatflx
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmpYheatflx(Temp,dens,bf,kspitzer,myflux,dx,dt,ngrid,compute,Temp2)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer,intent(IN) ::ngrid,compute
  real(dp),intent(IN)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2),intent(OUT)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3),intent(IN)::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(IN)::kspitzer
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(IN)::Temp,Temp2,dens

  real(dp)::fy,oneovertwodx,oneoverfourdx
  real(dp)::Bnorm2_1,Bnorm2_2,Bnorm2_3,Bnorm2_4
  real(dp)::dTdx1,dTdx2,dTdx3,dTdx4
  real(dp)::dTdy1,dTdy2,dTdy3,dTdy4
  real(dp)::dTdz1,dTdz2,dTdz3,dTdz4
  real(dp)::bx1,bx2,bx3,bx4
  real(dp)::by1,by2,by3,by4
  real(dp)::bz1,bz2,bz3,bz4
  real(dp)::fy1,fy2,fy3,fy4
  real(dp)::kpary1,kpary2,kpary3,kpary4,oneminuskperp
  real(dp)::dmean1,dmean2,dmean3,dmean4
  real(dp)::Tmean1,Tmean2,Tmean3,Tmean4
  real(dp)::dTolddx1,dTolddx2,dTolddx3,dTolddx4
  real(dp)::dTolddy1,dTolddy2,dTolddy3,dTolddy4
  real(dp)::dTolddz1,dTolddz2,dTolddz3,dTolddz4
  real(dp)::sat_coef1,sat_coef2,sat_coef3,sat_coef4
  real(dp)::ratio1,ratio2,ratio3,ratio4
#if NDIM==2
  real(dp)::Bnorm
#endif
  ! Local scalar variables
  integer::i,j,k,l
  integer::ilo,ihi,klo,khi

  oneovertwodx =0.50d0/dx
  oneoverfourdx=0.25d0/dx
  oneminuskperp=1.0d0-k_perp

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

if(isotrope_cond)then
  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
        dTdy1  =(Temp(l,i,j,k)-Temp(l,i,j-1,k))/dx
        kpary1=0.5d0*(kspitzer(l,i,j,k)+kspitzer(l,i,j-1,k))
        if(saturation)then
           dmean1=0.5d0*(dens (l,i,j,k)+dens (l,i,j-1,k))
           Tmean1=0.5d0*(Temp2(l,i,j,k)+Temp2(l,i,j-1,k))
           dTolddy1=abs(Temp2(l,i,j,k)-Temp2(l,i,j-1,k))/dx
           ratio1=fudge_mfp*Tmean1*dTolddy1/dmean1
           sat_coef1=1d0/(1d0+4.2d0*ratio1)
           kpary1=sat_coef1*kpary1
        endif
        if(compute .ne. 3)then
           fy=kpary1*dTdy1
        else
           fy=kpary1
        endif
!!$        famr=max(ffdx(l,i,j,k),ffdx(l,i,j-1,k))
!!$        fy=fy/famr
        myflux(l,i,j,k)=fy*dt/dx
!!$        if(compute==3)myflux(l,i,j,k)=kpary1*dt/dx**2
     enddo
  enddo
  enddo
  enddo

else

  if (.not.slopelim_cond)then
  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
#if NDIM==2
!!$        --------
!!$        |      |
!!$        |      |
!!$        |      |
!!$        1------2
        
        if(compute .ne. 3)then
           dTdx1=(Temp(l,i  ,j  ,k)+Temp(l,i  ,j-1,k) &
                - Temp(l,i-1,j  ,k)-Temp(l,i-1,j-1,k))*oneovertwodx
           dTdx2=(Temp(l,i+1,j  ,k)+Temp(l,i+1,j-1,k) &
             - Temp(l,i  ,j  ,k)-Temp(l,i  ,j-1,k))*oneovertwodx        
           
           dTdy1=(Temp(l,i  ,j  ,k)+Temp(l,i-1,j  ,k) &
                - Temp(l,i  ,j-1,k)-Temp(l,i-1,j-1,k))*oneovertwodx
           dTdy2=(Temp(l,i  ,j  ,k)+Temp(l,i+1,j  ,k) &
                - Temp(l,i  ,j-1,k)-Temp(l,i+1,j-1,k))*oneovertwodx
        end if

        bx1=0.5d0*(bf(l,i  ,j-1,k,1)+bf(l,i  ,j  ,k,1))
        bx2=0.5d0*(bf(l,i+1,j  ,k,1)+bf(l,i+1,j-1,k,1))
        
        by1=0.5d0*(bf(l,i-1,j  ,k,2)+bf(l,i  ,j  ,k,2))
        by2=0.5d0*(bf(l,i  ,j  ,k,2)+bf(l,i+1,j  ,k,2))
        
        Bnorm=sqrt(bx1*bx1+by1*by1)
        if(Bnorm.gt.0.0)then
           bx1=bx1/Bnorm
           by1=by1/Bnorm
        endif
        Bnorm=sqrt(bx2*bx2+by2*by2)
        if(Bnorm.gt.0.0)then
           bx2=bx2/Bnorm
           by2=by2/Bnorm
        endif

        kpary1=0.25d0*(kspitzer(l,i  ,j  ,k)+kspitzer(l,i  ,j-1,k) &
             +         kspitzer(l,i-1,j  ,k)+kspitzer(l,i-1,j-1,k))
        kpary2=0.25d0*(kspitzer(l,i  ,j  ,k)+kspitzer(l,i+1,j  ,k) &
             +         kspitzer(l,i  ,j-1,k)+kspitzer(l,i+1,j-1,k))

        if(saturation)then
           dmean1=0.25d0*(dens (l,i  ,j  ,k)+dens (l,i  ,j-1,k) &
                +         dens (l,i-1,j  ,k)+dens (l,i-1,j-1,k))
           dmean2=0.25d0*(dens (l,i  ,j  ,k)+dens (l,i+1,j  ,k) &
                +         dens (l,i  ,j-1,k)+dens (l,i+1,j-1,k))
           Tmean1=0.25d0*(Temp2(l,i  ,j  ,k)+Temp2(l,i  ,j-1,k) &
                +         Temp2(l,i-1,j  ,k)+Temp2(l,i-1,j-1,k))
           Tmean2=0.25d0*(Temp2(l,i  ,j+1,k)+Temp2(l,i  ,j  ,k) &
                +         Temp2(l,i-1,j+1,k)+Temp2(l,i-1,j  ,k))
           dTolddx1=(Temp2(l,i  ,j  ,k)+Temp2(l,i  ,j-1,k) &
                       - Temp2(l,i-1,j  ,k)-Temp2(l,i-1,j-1,k))*oneovertwodx
           dTolddx2=(Temp2(l,i+1,j  ,k)+Temp2(l,i+1,j-1,k) &
                       - Temp2(l,i  ,j  ,k)-Temp2(l,i  ,j-1,k))*oneovertwodx        
           dTolddy1=(Temp2(l,i  ,j  ,k)+Temp2(l,i-1,j  ,k) &
                       - Temp2(l,i  ,j-1,k)-Temp2(l,i-1,j-1,k))*oneovertwodx
           dTolddy2=(Temp2(l,i  ,j  ,k)+Temp2(l,i+1,j  ,k) &
                       - Temp2(l,i  ,j-1,k)-Temp2(l,i+1,j-1,k))*oneovertwodx
           ratio1=fudge_mfp*Tmean1*dsqrt(MAX(dTolddx1**2+dTolddy1**2,0d0))/dmean1
           ratio2=fudge_mfp*Tmean2*dsqrt(MAX(dTolddx2**2+dTolddy2**2,0d0))/dmean2
           sat_coef1=1d0/(1d0+4.2d0*ratio1)
           sat_coef2=1d0/(1d0+4.2d0*ratio2)
           kpary1=sat_coef1*kpary1
           kpary2=sat_coef2*kpary2
        endif

        if(compute .ne. 3)then   
           fy1=kpary1*(by1*oneminuskperp*(bx1*dTdx1+by1*dTdy1)+k_perp*dTdy1)
           fy2=kpary2*(by2*oneminuskperp*(bx2*dTdx2+by2*dTdy2)+k_perp*dTdy2)
        else ! Preconditionner
           fy1=kpary1*(by1*oneminuskperp*(bx1+by1)+k_perp)
           fy2=kpary2*(by2*oneminuskperp*(bx2+by2)+k_perp)
        end if
        fy=0.5d0*(fy1+fy2)
#endif
#if NDIM==3
!!$          ---------
!!$         / |      /|
!!$        /  |     / |
!!$        3------4   |
!!$        |  |   |   |
!!$        | /    |  /
!!$        |/     | /
!!$        1------2
        
        ! Centered symmetric scheme
        if(compute .ne. 3)then
           dTdx1=(Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  )+Temp(l,i  ,j  ,k-1)+Temp(l,i  ,j-1,k-1) &
                - Temp(l,i-1,j  ,k  )-Temp(l,i-1,j-1,k  )-Temp(l,i-1,j  ,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdx2=(Temp(l,i+1,j  ,k  )+Temp(l,i+1,j-1,k  )+Temp(l,i+1,j  ,k-1)+Temp(l,i+1,j-1,k-1) &
                - Temp(l,i  ,j  ,k  )-Temp(l,i  ,j-1,k  )-Temp(l,i  ,j  ,k-1)-Temp(l,i  ,j-1,k-1))*oneoverfourdx
           dTdx3=(Temp(l,i  ,j  ,k+1)+Temp(l,i  ,j-1,k+1)+Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  ) &
                - Temp(l,i-1,j  ,k+1)-Temp(l,i-1,j-1,k+1)-Temp(l,i-1,j  ,k  )-Temp(l,i-1,j-1,k  ))*oneoverfourdx
           dTdx4=(Temp(l,i+1,j  ,k+1)+Temp(l,i+1,j-1,k+1)+Temp(l,i+1,j  ,k  )+Temp(l,i+1,j-1,k  ) &
                - Temp(l,i  ,j  ,k+1)-Temp(l,i  ,j-1,k+1)-Temp(l,i  ,j  ,k  )-Temp(l,i  ,j-1,k  ))*oneoverfourdx
           
           dTdy1=(Temp(l,i  ,j  ,k  )+Temp(l,i-1,j  ,k  )+Temp(l,i  ,j  ,k-1)+Temp(l,i-1,j  ,k-1) &
                - Temp(l,i  ,j-1,k  )-Temp(l,i-1,j-1,k  )-Temp(l,i  ,j-1,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdy2=(Temp(l,i+1,j  ,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i+1,j  ,k-1)+Temp(l,i  ,j  ,k-1) &
                - Temp(l,i+1,j-1,k  )-Temp(l,i  ,j-1,k  )-Temp(l,i+1,j-1,k-1)-Temp(l,i  ,j-1,k-1))*oneoverfourdx
           dTdy3=(Temp(l,i  ,j  ,k+1)+Temp(l,i-1,j  ,k+1)+Temp(l,i  ,j  ,k  )+Temp(l,i-1,j  ,k  ) &
                - Temp(l,i  ,j-1,k+1)-Temp(l,i-1,j-1,k+1)-Temp(l,i  ,j-1,k  )-Temp(l,i-1,j-1,k  ))*oneoverfourdx
           dTdy4=(Temp(l,i+1,j  ,k+1)+Temp(l,i  ,j  ,k+1)+Temp(l,i+1,j  ,k  )+Temp(l,i  ,j  ,k  ) &
                - Temp(l,i+1,j-1,k+1)-Temp(l,i  ,j-1,k+1)-Temp(l,i+1,j-1,k  )-Temp(l,i  ,j-1,k  ))*oneoverfourdx
           
           dTdz1=(Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  )+Temp(l,i-1,j  ,k  )+Temp(l,i-1,j-1,k  ) &
                - Temp(l,i  ,j  ,k-1)-Temp(l,i  ,j-1,k-1)-Temp(l,i-1,j  ,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdz2=(Temp(l,i+1,j  ,k  )+Temp(l,i+1,j-1,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  ) &
                - Temp(l,i+1,j  ,k-1)-Temp(l,i+1,j-1,k-1)-Temp(l,i  ,j  ,k-1)-Temp(l,i  ,j-1,k-1))*oneoverfourdx
           dTdz3=(Temp(l,i  ,j  ,k+1)+Temp(l,i-1,j  ,k+1)+Temp(l,i  ,j-1,k+1)+Temp(l,i-1,j-1,k+1) &
                - Temp(l,i  ,j  ,k  )-Temp(l,i-1,j  ,k  )-Temp(l,i  ,j-1,k  )-Temp(l,i-1,j-1,k  ))*oneoverfourdx
           dTdz4=(Temp(l,i+1,j  ,k+1)+Temp(l,i  ,j  ,k+1)+Temp(l,i+1,j-1,k+1)+Temp(l,i  ,j-1,k+1) &
                - Temp(l,i+1,j  ,k  )-Temp(l,i  ,j  ,k  )-Temp(l,i+1,j-1,k  )-Temp(l,i  ,j-1,k  ))*oneoverfourdx
        end if

        bx1=(bf(l,i  ,j-1,k-1,1)+bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j-1,k  ,1)+bf(l,i  ,j  ,k  ,1))
        bx2=(bf(l,i+1,j-1,k-1,1)+bf(l,i+1,j  ,k-1,1)+bf(l,i+1,j-1,k  ,1)+bf(l,i+1,j  ,k  ,1))
        bx3=(bf(l,i  ,j-1,k  ,1)+bf(l,i  ,j  ,k  ,1)+bf(l,i  ,j-1,k+1,1)+bf(l,i  ,j  ,k+1,1))
        bx4=(bf(l,i+1,j-1,k  ,1)+bf(l,i+1,j  ,k  ,1)+bf(l,i+1,j-1,k+1,1)+bf(l,i+1,j  ,k+1,1))

        by1=(bf(l,i-1,j  ,k-1,2)+bf(l,i  ,j  ,k-1,2)+bf(l,i-1,j  ,k  ,2)+bf(l,i  ,j  ,k  ,2))
        by2=(bf(l,i  ,j  ,k-1,2)+bf(l,i+1,j  ,k-1,2)+bf(l,i  ,j  ,k  ,2)+bf(l,i+1,j  ,k  ,2))
        by3=(bf(l,i-1,j  ,k  ,2)+bf(l,i  ,j  ,k  ,2)+bf(l,i-1,j  ,k+1,2)+bf(l,i  ,j  ,k+1,2))
        by4=(bf(l,i  ,j  ,k  ,2)+bf(l,i+1,j  ,k  ,2)+bf(l,i  ,j  ,k+1,2)+bf(l,i+1,j  ,k+1,2))
        
        bz1=(bf(l,i-1,j-1,k  ,3)+bf(l,i-1,j  ,k  ,3)+bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3))
        bz2=(bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3)+bf(l,i+1,j-1,k  ,3)+bf(l,i+1,j  ,k  ,3))
        bz3=(bf(l,i  ,j-1,k+1,3)+bf(l,i  ,j  ,k+1,3)+bf(l,i-1,j-1,k+1,3)+bf(l,i-1,j  ,k+1,3))
        bz4=(bf(l,i+1,j-1,k+1,3)+bf(l,i+1,j  ,k+1,3)+bf(l,i  ,j-1,k+1,3)+bf(l,i  ,j  ,k+1,3))

        Bnorm2_1=(bx1*bx1+by1*by1+bz1*bz1)
        Bnorm2_2=(bx2*bx2+by2*by2+bz2*bz2)
        Bnorm2_3=(bx3*bx3+by3*by3+bz3*bz3)
        Bnorm2_4=(bx4*bx4+by4*by4+bz4*bz4)

        kpary1=0.125d0*(kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j-1,k  )+kspitzer(l,i  ,j  ,k-1) &
             +      kspitzer(l,i  ,j-1,k-1)+kspitzer(l,i-1,j  ,k  )+kspitzer(l,i-1,j-1,k  ) &
             +      kspitzer(l,i-1,j  ,k-1)+kspitzer(l,i-1,j-1,k-1))
        kpary2=0.125d0*(kspitzer(l,i+1,j  ,k  )+kspitzer(l,i+1,j-1,k  )+kspitzer(l,i+1,j  ,k-1) &
             +      kspitzer(l,i+1,j-1,k-1)+kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j-1,k  ) &
             +      kspitzer(l,i  ,j  ,k-1)+kspitzer(l,i  ,j-1,k-1))
        kpary3=0.125d0*(kspitzer(l,i  ,j  ,k+1)+kspitzer(l,i  ,j-1,k+1)+kspitzer(l,i  ,j  ,k  ) &
             +      kspitzer(l,i  ,j-1,k  )+kspitzer(l,i-1,j  ,k+1)+kspitzer(l,i-1,j-1,k+1) &
             +      kspitzer(l,i-1,j  ,k  )+kspitzer(l,i-1,j-1,k  ))
        kpary4=0.125d0*(kspitzer(l,i+1,j  ,k+1)+kspitzer(l,i+1,j-1,k+1)+kspitzer(l,i+1,j  ,k  ) &
             +      kspitzer(l,i+1,j-1,k  )+kspitzer(l,i  ,j  ,k+1)+kspitzer(l,i  ,j-1,k+1) &
             +      kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j-1,k  ))

        if(saturation)then
           dmean1=0.125d0*(dens (l,i  ,j  ,k  )+dens (l,i  ,j-1,k  )+dens (l,i  ,j  ,k-1) &
                +          dens (l,i  ,j-1,k-1)+dens (l,i-1,j  ,k  )+dens (l,i-1,j-1,k  ) &
                +          dens (l,i-1,j  ,k-1)+dens (l,i-1,j-1,k-1))
           dmean2=0.125d0*(dens (l,i+1,j  ,k  )+dens (l,i+1,j-1,k  )+dens (l,i+1,j  ,k-1) &
                +          dens (l,i+1,j-1,k-1)+dens (l,i  ,j  ,k  )+dens (l,i  ,j-1,k  ) &
                +          dens (l,i  ,j  ,k-1)+dens (l,i  ,j-1,k-1))
           dmean3=0.125d0*(dens (l,i  ,j  ,k+1)+dens (l,i  ,j-1,k+1)+dens (l,i  ,j  ,k  ) &
                +          dens (l,i  ,j-1,k  )+dens (l,i-1,j  ,k+1)+dens (l,i-1,j-1,k+1) &
                +          dens (l,i-1,j  ,k  )+dens (l,i-1,j-1,k  ))
           dmean4=0.125d0*(dens (l,i+1,j  ,k+1)+dens (l,i+1,j-1,k+1)+dens (l,i+1,j  ,k  ) &
                +          dens (l,i+1,j-1,k  )+dens (l,i  ,j  ,k+1)+dens (l,i  ,j-1,k+1) &
                +          dens (l,i  ,j  ,k  )+dens (l,i  ,j-1,k  ))
           Tmean1=0.125d0*(Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  )+Temp2(l,i  ,j  ,k-1) &
                +          Temp2(l,i  ,j-1,k-1)+Temp2(l,i-1,j  ,k  )+Temp2(l,i-1,j-1,k  ) &
                +          Temp2(l,i-1,j  ,k-1)+Temp2(l,i-1,j-1,k-1))
           Tmean2=0.125d0*(Temp2(l,i+1,j  ,k  )+Temp2(l,i+1,j-1,k  )+Temp2(l,i+1,j  ,k-1) &
                +          Temp2(l,i+1,j-1,k-1)+Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  ) &
                +          Temp2(l,i  ,j  ,k-1)+Temp2(l,i  ,j-1,k-1))
           Tmean3=0.125d0*(Temp2(l,i  ,j  ,k+1)+Temp2(l,i  ,j-1,k+1)+Temp2(l,i  ,j  ,k  ) &
                +          Temp2(l,i  ,j-1,k  )+Temp2(l,i-1,j  ,k+1)+Temp2(l,i-1,j-1,k+1) &
                +          Temp2(l,i-1,j  ,k  )+Temp2(l,i-1,j-1,k  ))
           Tmean4=0.125d0*(Temp2(l,i+1,j  ,k+1)+Temp2(l,i+1,j-1,k+1)+Temp2(l,i+1,j  ,k  ) &
                +          Temp2(l,i+1,j-1,k  )+Temp2(l,i  ,j  ,k+1)+Temp2(l,i  ,j-1,k+1) &
                +          Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  ))
           dTolddx1=(Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  )+Temp2(l,i  ,j  ,k-1)+Temp2(l,i  ,j-1,k-1) &
                - Temp2(l,i-1,j  ,k  )-Temp2(l,i-1,j-1,k  )-Temp2(l,i-1,j  ,k-1)-Temp2(l,i-1,j-1,k-1))*oneoverfourdx
           dTolddx2=(Temp2(l,i+1,j  ,k  )+Temp2(l,i+1,j-1,k  )+Temp2(l,i+1,j  ,k-1)+Temp2(l,i+1,j-1,k-1) &
                - Temp2(l,i  ,j  ,k  )-Temp2(l,i  ,j-1,k  )-Temp2(l,i  ,j  ,k-1)-Temp2(l,i  ,j-1,k-1))*oneoverfourdx
           dTolddx3=(Temp2(l,i  ,j  ,k+1)+Temp2(l,i  ,j-1,k+1)+Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  ) &
                - Temp2(l,i-1,j  ,k+1)-Temp2(l,i-1,j-1,k+1)-Temp2(l,i-1,j  ,k  )-Temp2(l,i-1,j-1,k  ))*oneoverfourdx
           dTolddx4=(Temp2(l,i+1,j  ,k+1)+Temp2(l,i+1,j-1,k+1)+Temp2(l,i+1,j  ,k  )+Temp2(l,i+1,j-1,k  ) &
                - Temp2(l,i  ,j  ,k+1)-Temp2(l,i  ,j-1,k+1)-Temp2(l,i  ,j  ,k  )-Temp2(l,i  ,j-1,k  ))*oneoverfourdx
           dTolddy1=(Temp2(l,i  ,j  ,k  )+Temp2(l,i-1,j  ,k  )+Temp2(l,i  ,j  ,k-1)+Temp2(l,i-1,j  ,k-1) &
                - Temp2(l,i  ,j-1,k  )-Temp2(l,i-1,j-1,k  )-Temp2(l,i  ,j-1,k-1)-Temp2(l,i-1,j-1,k-1))*oneoverfourdx
           dTolddy2=(Temp2(l,i+1,j  ,k  )+Temp2(l,i  ,j  ,k  )+Temp2(l,i+1,j  ,k-1)+Temp2(l,i  ,j  ,k-1) &
                - Temp2(l,i+1,j-1,k  )-Temp2(l,i  ,j-1,k  )-Temp2(l,i+1,j-1,k-1)-Temp2(l,i  ,j-1,k-1))*oneoverfourdx
           dTolddy3=(Temp2(l,i  ,j  ,k+1)+Temp2(l,i-1,j  ,k+1)+Temp2(l,i  ,j  ,k  )+Temp2(l,i-1,j  ,k  ) &
                - Temp2(l,i  ,j-1,k+1)-Temp2(l,i-1,j-1,k+1)-Temp2(l,i  ,j-1,k  )-Temp2(l,i-1,j-1,k  ))*oneoverfourdx
           dTolddy4=(Temp2(l,i+1,j  ,k+1)+Temp2(l,i  ,j  ,k+1)+Temp2(l,i+1,j  ,k  )+Temp2(l,i  ,j  ,k  ) &
                - Temp2(l,i+1,j-1,k+1)-Temp2(l,i  ,j-1,k+1)-Temp2(l,i+1,j-1,k  )-Temp2(l,i  ,j-1,k  ))*oneoverfourdx
           dTolddz1=(Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  )+Temp2(l,i-1,j  ,k  )+Temp2(l,i-1,j-1,k  ) &
                - Temp2(l,i  ,j  ,k-1)-Temp2(l,i  ,j-1,k-1)-Temp2(l,i-1,j  ,k-1)-Temp2(l,i-1,j-1,k-1))*oneoverfourdx
           dTolddz2=(Temp2(l,i+1,j  ,k  )+Temp2(l,i+1,j-1,k  )+Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  ) &
                - Temp2(l,i+1,j  ,k-1)-Temp2(l,i+1,j-1,k-1)-Temp2(l,i  ,j  ,k-1)-Temp2(l,i  ,j-1,k-1))*oneoverfourdx
           dTolddz3=(Temp2(l,i  ,j  ,k+1)+Temp2(l,i-1,j  ,k+1)+Temp2(l,i  ,j-1,k+1)+Temp2(l,i-1,j-1,k+1) &
                - Temp2(l,i  ,j  ,k  )-Temp2(l,i-1,j  ,k  )-Temp2(l,i  ,j-1,k  )-Temp2(l,i-1,j-1,k  ))*oneoverfourdx
           dTolddz4=(Temp2(l,i+1,j  ,k+1)+Temp2(l,i  ,j  ,k+1)+Temp2(l,i+1,j-1,k+1)+Temp2(l,i  ,j-1,k+1) &
                - Temp2(l,i+1,j  ,k  )-Temp2(l,i  ,j  ,k  )-Temp2(l,i+1,j-1,k  )-Temp2(l,i  ,j-1,k  ))*oneoverfourdx
           ratio1=fudge_mfp*Tmean1*dsqrt(MAX(dTolddx1**2+dTolddy1**2+dTolddz1**2,0.0d0))/dmean1
           ratio2=fudge_mfp*Tmean2*dsqrt(MAX(dTolddx2**2+dTolddy2**2+dTolddz2**2,0.0d0))/dmean2
           ratio3=fudge_mfp*Tmean3*dsqrt(MAX(dTolddx3**2+dTolddy3**2+dTolddz3**2,0.0d0))/dmean3
           ratio4=fudge_mfp*Tmean4*dsqrt(MAX(dTolddx4**2+dTolddy4**2+dTolddz4**2,0.0d0))/dmean4
           sat_coef1=1d0/(1d0+4.2d0*ratio1)
           sat_coef2=1d0/(1d0+4.2d0*ratio2)
           sat_coef3=1d0/(1d0+4.2d0*ratio3)
           sat_coef4=1d0/(1d0+4.2d0*ratio4)
           kpary1=sat_coef1*kpary1
           kpary2=sat_coef2*kpary2
           kpary3=sat_coef3*kpary3
           kpary4=sat_coef4*kpary4
        endif

        if(compute .ne. 3)then        
           fy1=kpary1*(by1*oneminuskperp*(bx1*dTdx1+by1*dTdy1+bz1*dTdz1)/Bnorm2_1+k_perp*dTdy1)
           fy2=kpary2*(by2*oneminuskperp*(bx2*dTdx2+by2*dTdy2+bz2*dTdz2)/Bnorm2_2+k_perp*dTdy2)
           fy3=kpary3*(by3*oneminuskperp*(bx3*dTdx3+by3*dTdy3+bz3*dTdz3)/Bnorm2_3+k_perp*dTdy3)
           fy4=kpary4*(by4*oneminuskperp*(bx4*dTdx4+by4*dTdy4+bz4*dTdz4)/Bnorm2_4+k_perp*dTdy4)
        else
           fy1=kpary1*(by1*oneminuskperp*(bx1+by1+bz1)/Bnorm2_1+k_perp)
           fy2=kpary2*(by2*oneminuskperp*(bx2+by2+bz2)/Bnorm2_2+k_perp)
           fy3=kpary3*(by3*oneminuskperp*(bx3+by3+bz3)/Bnorm2_3+k_perp)
           fy4=kpary4*(by4*oneminuskperp*(bx4+by4+bz4)/Bnorm2_4+k_perp)
        end if
        fy=0.25d0*(fy1+fy2+fy3+fy4)
#endif

        myflux(l,i,j,k)=fy*dt/dx
     enddo
  enddo
  enddo
  enddo
  endif

  if (slopelim_cond)then
     ! TODO
  endif
endif

end subroutine cmpYheatflx
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmpZheatflx(Temp,dens,bf,kspitzer,myflux,dx,dt,ngrid,compute,Temp2)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer,intent(IN) ::ngrid,compute
  real(dp),intent(IN)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2),intent(OUT)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3),intent(IN)::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(IN)::kspitzer
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(IN)::Temp,Temp2,dens

  real(dp)::fz,oneovertwodx,oneoverfourdx
  real(dp)::Bnorm2_1,Bnorm2_2,Bnorm2_3,Bnorm2_4
  real(dp)::dTdx1,dTdx2,dTdx3,dTdx4
  real(dp)::dTdy1,dTdy2,dTdy3,dTdy4
  real(dp)::dTdz1,dTdz2,dTdz3,dTdz4
  real(dp)::bx1,bx2,bx3,bx4
  real(dp)::by1,by2,by3,by4
  real(dp)::bz1,bz2,bz3,bz4
  real(dp)::fz1,fz2,fz3,fz4
  real(dp)::kparz1,kparz2,kparz3,kparz4,oneminuskperp
  real(dp)::dmean1,dmean2,dmean3,dmean4
  real(dp)::Tmean1,Tmean2,Tmean3,Tmean4
  real(dp)::dTolddx1,dTolddx2,dTolddx3,dTolddx4
  real(dp)::dTolddy1,dTolddy2,dTolddy3,dTolddy4
  real(dp)::dTolddz1,dTolddz2,dTolddz3,dTolddz4
  real(dp)::sat_coef1,sat_coef2,sat_coef3,sat_coef4
  real(dp)::ratio1,ratio2,ratio3,ratio4

  ! Local scalar variables
  integer::i,j,k,l
  integer::ilo,ihi,jlo,jhi

  oneovertwodx =0.50d0/dx
  oneoverfourdx=0.25d0/dx
  oneminuskperp=1.0d0-k_perp

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)

if(isotrope_cond)then
  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do l = 1, ngrid
        dTdz1  =(Temp(l,i,j,k)-Temp(l,i,j,k-1))/dx
        kparz1=0.5d0*(kspitzer(l,i,j,k)+kspitzer(l,i,j,k-1))
        if(saturation)then
           dmean1=0.5d0*(dens (l,i,j,k-1)+dens (l,i,j,k))
           Tmean1=0.5d0*(Temp2(l,i,j,k-1)+Temp2(l,i,j,k))
           dTolddz1=abs(Temp2(l,i,j,k)-Temp2(l,i,j,k-1))/dx
           ratio1=fudge_mfp*Tmean1*dTolddz1/dmean1
           sat_coef1=1d0/(1d0+4.2d0*ratio1)
           kparz1=sat_coef1*kparz1
        endif
        if(compute .ne. 3)then
           fz=kparz1*dTdz1
        else
           fz=kparz1
        endif
!!$        famr=max(ffdx(l,i,j,k),ffdx(l,i,j,k-1))
!!$        fz=fz/famr
        myflux(l,i,j,k)=fz*dt/dx
!!$        if(compute==3)myflux(l,i,j,k)=kparz1*dt/dx**2
     enddo
  enddo
  enddo
  enddo

else

  if (.not.slopelim_cond)then
  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do l = 1, ngrid
#if NDIM==3
!!$          ---------
!!$         / |      /|
!!$        /  |     / |
!!$        --------   |
!!$        |  |   |   |
!!$        | /3   |  /4
!!$        |/     | /
!!$        1------2
        
        ! Centered symmetric scheme
        if(compute .ne. 3)then        
           dTdx1=(Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  )+Temp(l,i  ,j  ,k-1)+Temp(l,i  ,j-1,k-1) &
                - Temp(l,i-1,j  ,k  )-Temp(l,i-1,j-1,k  )-Temp(l,i-1,j  ,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdx2=(Temp(l,i+1,j  ,k  )+Temp(l,i+1,j-1,k  )+Temp(l,i+1,j  ,k-1)+Temp(l,i+1,j-1,k-1) &
                - Temp(l,i  ,j  ,k  )-Temp(l,i  ,j-1,k  )-Temp(l,i  ,j  ,k-1)-Temp(l,i  ,j-1,k-1))*oneoverfourdx
           dTdx3=(Temp(l,i  ,j+1,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i  ,j+1,k-1)+Temp(l,i  ,j  ,k-1) &
                - Temp(l,i-1,j+1,k  )-Temp(l,i-1,j  ,k  )-Temp(l,i-1,j+1,k-1)-Temp(l,i-1,j  ,k-1))*oneoverfourdx
           dTdx4=(Temp(l,i+1,j+1,k  )+Temp(l,i+1,j  ,k  )+Temp(l,i+1,j+1,k-1)+Temp(l,i+1,j  ,k-1) &
                - Temp(l,i  ,j+1,k  )-Temp(l,i  ,j  ,k  )-Temp(l,i  ,j+1,k-1)-Temp(l,i  ,j  ,k-1))*oneoverfourdx
           
           dTdy1=(Temp(l,i  ,j  ,k  )+Temp(l,i-1,j  ,k  )+Temp(l,i  ,j  ,k-1)+Temp(l,i-1,j  ,k-1) &
                - Temp(l,i  ,j-1,k  )-Temp(l,i-1,j-1,k  )-Temp(l,i  ,j-1,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdy2=(Temp(l,i+1,j  ,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i+1,j  ,k-1)+Temp(l,i  ,j  ,k-1) &
                - Temp(l,i+1,j-1,k  )-Temp(l,i  ,j-1,k  )-Temp(l,i+1,j-1,k-1)-Temp(l,i  ,j-1,k-1))*oneoverfourdx
           dTdy3=(Temp(l,i  ,j+1,k  )+Temp(l,i-1,j+1,k  )+Temp(l,i  ,j+1,k-1)+Temp(l,i-1,j+1,k-1) &
                - Temp(l,i  ,j  ,k  )-Temp(l,i-1,j  ,k  )-Temp(l,i  ,j  ,k-1)-Temp(l,i-1,j  ,k-1))*oneoverfourdx
           dTdy4=(Temp(l,i+1,j+1,k  )+Temp(l,i  ,j+1,k  )+Temp(l,i+1,j+1,k-1)+Temp(l,i  ,j+1,k-1) &
                - Temp(l,i+1,j  ,k  )-Temp(l,i  ,j  ,k  )-Temp(l,i+1,j  ,k-1)-Temp(l,i  ,j  ,k-1))*oneoverfourdx
           
           dTdz1=(Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  )+Temp(l,i-1,j  ,k  )+Temp(l,i-1,j-1,k  ) &
                - Temp(l,i  ,j  ,k-1)-Temp(l,i  ,j-1,k-1)-Temp(l,i-1,j  ,k-1)-Temp(l,i-1,j-1,k-1))*oneoverfourdx
           dTdz2=(Temp(l,i+1,j  ,k  )+Temp(l,i+1,j-1,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i  ,j-1,k  ) &
                - Temp(l,i+1,j  ,k-1)-Temp(l,i+1,j-1,k-1)-Temp(l,i  ,j  ,k-1)-Temp(l,i  ,j-1,k-1))*oneoverfourdx
           dTdz3=(Temp(l,i  ,j+1,k  )+Temp(l,i  ,j  ,k  )+Temp(l,i-1,j+1,k  )+Temp(l,i-1,j  ,k  ) &
                - Temp(l,i  ,j+1,k-1)-Temp(l,i  ,j  ,k-1)-Temp(l,i-1,j+1,k-1)-Temp(l,i-1,j  ,k-1))*oneoverfourdx
           dTdz4=(Temp(l,i+1,j+1,k  )+Temp(l,i+1,j  ,k  )+Temp(l,i  ,j+1,k  )+Temp(l,i  ,j  ,k  ) &
                - Temp(l,i+1,j+1,k-1)-Temp(l,i+1,j  ,k-1)-Temp(l,i  ,j+1,k-1)-Temp(l,i  ,j  ,k-1))*oneoverfourdx
        end if
        
        bx1=(bf(l,i  ,j-1,k-1,1)+bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j-1,k  ,1)+bf(l,i  ,j  ,k  ,1))
        bx2=(bf(l,i+1,j-1,k-1,1)+bf(l,i+1,j  ,k-1,1)+bf(l,i+1,j-1,k  ,1)+bf(l,i+1,j  ,k  ,1))
        bx3=(bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j+1,k-1,1)+bf(l,i  ,j  ,k  ,1)+bf(l,i  ,j+1,k  ,1))
        bx4=(bf(l,i+1,j  ,k-1,1)+bf(l,i+1,j+1,k-1,1)+bf(l,i+1,j  ,k  ,1)+bf(l,i+1,j+1,k  ,1))

        by1=(bf(l,i-1,j  ,k-1,2)+bf(l,i  ,j  ,k-1,2)+bf(l,i-1,j  ,k  ,2)+bf(l,i  ,j  ,k  ,2))
        by2=(bf(l,i  ,j  ,k-1,2)+bf(l,i+1,j  ,k-1,2)+bf(l,i  ,j  ,k  ,2)+bf(l,i+1,j  ,k  ,2))
        by3=(bf(l,i-1,j+1,k-1,2)+bf(l,i  ,j+1,k-1,2)+bf(l,i-1,j+1,k  ,2)+bf(l,i  ,j+1,k  ,2))
        by4=(bf(l,i  ,j+1,k-1,2)+bf(l,i+1,j+1,k-1,2)+bf(l,i  ,j+1,k  ,2)+bf(l,i+1,j+1,k  ,2))
        
        bz1=(bf(l,i-1,j-1,k  ,3)+bf(l,i-1,j  ,k  ,3)+bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3))
        bz2=(bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3)+bf(l,i+1,j-1,k  ,3)+bf(l,i+1,j  ,k  ,3))
        bz3=(bf(l,i-1,j  ,k  ,3)+bf(l,i-1,j+1,k  ,3)+bf(l,i  ,j  ,k  ,3)+bf(l,i  ,j+1,k  ,3))
        bz4=(bf(l,i  ,j  ,k  ,3)+bf(l,i  ,j+1,k  ,3)+bf(l,i+1,j  ,k  ,3)+bf(l,i+1,j+1,k  ,3))

        Bnorm2_1=(bx1*bx1+by1*by1+bz1*bz1)
        Bnorm2_2=(bx2*bx2+by2*by2+bz2*bz2)
        Bnorm2_3=(bx3*bx3+by3*by3+bz3*bz3)
        Bnorm2_4=(bx4*bx4+by4*by4+bz4*bz4)

        kparz1=0.125d0*(kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j-1,k  )+kspitzer(l,i  ,j  ,k-1) &
             +          kspitzer(l,i  ,j-1,k-1)+kspitzer(l,i-1,j  ,k  )+kspitzer(l,i-1,j-1,k  ) &
             +          kspitzer(l,i-1,j  ,k-1)+kspitzer(l,i-1,j-1,k-1))
        kparz2=0.125d0*(kspitzer(l,i+1,j  ,k  )+kspitzer(l,i+1,j-1,k  )+kspitzer(l,i+1,j  ,k-1) &
             +          kspitzer(l,i+1,j-1,k-1)+kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j-1,k  ) &
             +          kspitzer(l,i  ,j  ,k-1)+kspitzer(l,i  ,j-1,k-1))
        kparz3=0.125d0*(kspitzer(l,i  ,j+1,k  )+kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j+1,k-1) &
             +          kspitzer(l,i  ,j  ,k-1)+kspitzer(l,i-1,j+1,k  )+kspitzer(l,i-1,j  ,k  ) &
             +          kspitzer(l,i-1,j+1,k-1)+kspitzer(l,i-1,j  ,k-1))
        kparz4=0.125d0*(kspitzer(l,i+1,j+1,k  )+kspitzer(l,i+1,j  ,k  )+kspitzer(l,i+1,j+1,k-1) &
             +          kspitzer(l,i+1,j  ,k-1)+kspitzer(l,i  ,j+1,k  )+kspitzer(l,i  ,j  ,k  ) &
             +          kspitzer(l,i  ,j+1,k-1)+kspitzer(l,i  ,j  ,k-1))

        if(saturation)then
           dmean1=0.125d0*(dens (l,i  ,j  ,k  )+dens (l,i  ,j-1,k  )+dens (l,i  ,j  ,k-1) &
                +          dens (l,i  ,j-1,k-1)+dens (l,i-1,j  ,k  )+dens (l,i-1,j-1,k  ) &
                +          dens (l,i-1,j  ,k-1)+dens (l,i-1,j-1,k-1))
           dmean2=0.125d0*(dens (l,i+1,j  ,k  )+dens (l,i+1,j-1,k  )+dens (l,i+1,j  ,k-1) &
                +          dens (l,i+1,j-1,k-1)+dens (l,i  ,j  ,k  )+dens (l,i  ,j-1,k  ) &
                +          dens (l,i  ,j  ,k-1)+dens (l,i  ,j-1,k-1))
           dmean3=0.125d0*(dens (l,i  ,j+1,k  )+dens (l,i  ,j  ,k  )+dens (l,i  ,j+1,k-1) &
                +          dens (l,i  ,j  ,k-1)+dens (l,i-1,j+1,k  )+dens (l,i-1,j  ,k  ) &
                +          dens (l,i-1,j+1,k-1)+dens (l,i-1,j  ,k-1))
           dmean4=0.125d0*(dens (l,i+1,j+1,k  )+dens (l,i+1,j  ,k  )+dens (l,i+1,j+1,k-1) &
                +          dens (l,i+1,j  ,k-1)+dens (l,i  ,j+1,k  )+dens (l,i  ,j  ,k  ) &
                +          dens (l,i  ,j+1,k-1)+dens (l,i  ,j  ,k-1))
           Tmean1=0.125d0*(Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  )+Temp2(l,i  ,j  ,k-1) &
                +          Temp2(l,i  ,j-1,k-1)+Temp2(l,i-1,j  ,k  )+Temp2(l,i-1,j-1,k  ) &
                +          Temp2(l,i-1,j  ,k-1)+Temp2(l,i-1,j-1,k-1))
           Tmean2=0.125d0*(Temp2(l,i+1,j  ,k  )+Temp2(l,i+1,j-1,k  )+Temp2(l,i+1,j  ,k-1) &
                +          Temp2(l,i+1,j-1,k-1)+Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  ) &
                +          Temp2(l,i  ,j  ,k-1)+Temp2(l,i  ,j-1,k-1))
           Tmean3=0.125d0*(Temp2(l,i  ,j+1,k  )+Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j+1,k-1) &
                +          Temp2(l,i  ,j  ,k-1)+Temp2(l,i-1,j+1,k  )+Temp2(l,i-1,j  ,k  ) &
                +          Temp2(l,i-1,j+1,k-1)+Temp2(l,i-1,j  ,k-1))
           Tmean4=0.125d0*(Temp2(l,i+1,j+1,k  )+Temp2(l,i+1,j  ,k  )+Temp2(l,i+1,j+1,k-1) &
                +          Temp2(l,i+1,j  ,k-1)+Temp2(l,i  ,j+1,k  )+Temp2(l,i  ,j  ,k  ) &
                +          Temp2(l,i  ,j+1,k-1)+Temp2(l,i  ,j  ,k-1))
           dTolddx1=(Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  )+Temp2(l,i  ,j  ,k-1)+Temp2(l,i  ,j-1,k-1) &
                   - Temp2(l,i-1,j  ,k  )-Temp2(l,i-1,j-1,k  )-Temp2(l,i-1,j  ,k-1)-Temp2(l,i-1,j-1,k-1))*oneoverfourdx
           dTolddx2=(Temp2(l,i+1,j  ,k  )+Temp2(l,i+1,j-1,k  )+Temp2(l,i+1,j  ,k-1)+Temp2(l,i+1,j-1,k-1) &
                   - Temp2(l,i  ,j  ,k  )-Temp2(l,i  ,j-1,k  )-Temp2(l,i  ,j  ,k-1)-Temp2(l,i  ,j-1,k-1))*oneoverfourdx
           dTolddx3=(Temp2(l,i  ,j+1,k  )+Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j+1,k-1)+Temp2(l,i  ,j  ,k-1) &
                   - Temp2(l,i-1,j+1,k  )-Temp2(l,i-1,j  ,k  )-Temp2(l,i-1,j+1,k-1)-Temp2(l,i-1,j  ,k-1))*oneoverfourdx
           dTolddx4=(Temp2(l,i+1,j+1,k  )+Temp2(l,i+1,j  ,k  )+Temp2(l,i+1,j+1,k-1)+Temp2(l,i+1,j  ,k-1) &
                   - Temp2(l,i  ,j+1,k  )-Temp2(l,i  ,j  ,k  )-Temp2(l,i  ,j+1,k-1)-Temp2(l,i  ,j  ,k-1))*oneoverfourdx
           dTolddy1=(Temp2(l,i  ,j  ,k  )+Temp2(l,i-1,j  ,k  )+Temp2(l,i  ,j  ,k-1)+Temp2(l,i-1,j  ,k-1) &
                   - Temp2(l,i  ,j-1,k  )-Temp2(l,i-1,j-1,k  )-Temp2(l,i  ,j-1,k-1)-Temp2(l,i-1,j-1,k-1))*oneoverfourdx
           dTolddy2=(Temp2(l,i+1,j  ,k  )+Temp2(l,i  ,j  ,k  )+Temp2(l,i+1,j  ,k-1)+Temp2(l,i  ,j  ,k-1) &
                   - Temp2(l,i+1,j-1,k  )-Temp2(l,i  ,j-1,k  )-Temp2(l,i+1,j-1,k-1)-Temp2(l,i  ,j-1,k-1))*oneoverfourdx
           dTolddy3=(Temp2(l,i  ,j+1,k  )+Temp2(l,i-1,j+1,k  )+Temp2(l,i  ,j+1,k-1)+Temp2(l,i-1,j+1,k-1) &
                   - Temp2(l,i  ,j  ,k  )-Temp2(l,i-1,j  ,k  )-Temp2(l,i  ,j  ,k-1)-Temp2(l,i-1,j  ,k-1))*oneoverfourdx
           dTolddy4=(Temp2(l,i+1,j+1,k  )+Temp2(l,i  ,j+1,k  )+Temp2(l,i+1,j+1,k-1)+Temp2(l,i  ,j+1,k-1) &
                   - Temp2(l,i+1,j  ,k  )-Temp2(l,i  ,j  ,k  )-Temp2(l,i+1,j  ,k-1)-Temp2(l,i  ,j  ,k-1))*oneoverfourdx
           dTolddz1=(Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  )+Temp2(l,i-1,j  ,k  )+Temp2(l,i-1,j-1,k  ) &
                   - Temp2(l,i  ,j  ,k-1)-Temp2(l,i  ,j-1,k-1)-Temp2(l,i-1,j  ,k-1)-Temp2(l,i-1,j-1,k-1))*oneoverfourdx
           dTolddz2=(Temp2(l,i+1,j  ,k  )+Temp2(l,i+1,j-1,k  )+Temp2(l,i  ,j  ,k  )+Temp2(l,i  ,j-1,k  ) &
                   - Temp2(l,i+1,j  ,k-1)-Temp2(l,i+1,j-1,k-1)-Temp2(l,i  ,j  ,k-1)-Temp2(l,i  ,j-1,k-1))*oneoverfourdx
           dTolddz3=(Temp2(l,i  ,j+1,k  )+Temp2(l,i  ,j  ,k  )+Temp2(l,i-1,j+1,k  )+Temp2(l,i-1,j  ,k  ) &
                   - Temp2(l,i  ,j+1,k-1)-Temp2(l,i  ,j  ,k-1)-Temp2(l,i-1,j+1,k-1)-Temp2(l,i-1,j  ,k-1))*oneoverfourdx
           dTolddz4=(Temp2(l,i+1,j+1,k  )+Temp2(l,i+1,j  ,k  )+Temp2(l,i  ,j+1,k  )+Temp2(l,i  ,j  ,k  ) &
                   - Temp2(l,i+1,j+1,k-1)-Temp2(l,i+1,j  ,k-1)-Temp2(l,i  ,j+1,k-1)-Temp2(l,i  ,j  ,k-1))*oneoverfourdx
           ratio1=fudge_mfp*Tmean1*dsqrt(MAX(dTolddx1**2+dTolddy1**2+dTolddz1**2,0.0d0))/dmean1
           ratio2=fudge_mfp*Tmean2*dsqrt(MAX(dTolddx2**2+dTolddy2**2+dTolddz2**2,0.0d0))/dmean2
           ratio3=fudge_mfp*Tmean3*dsqrt(MAX(dTolddx3**2+dTolddy3**2+dTolddz3**2,0.0d0))/dmean3
           ratio4=fudge_mfp*Tmean4*dsqrt(MAX(dTolddx4**2+dTolddy4**2+dTolddz4**2,0.0d0))/dmean4
           sat_coef1=1d0/(1d0+4.2d0*ratio1)
           sat_coef2=1d0/(1d0+4.2d0*ratio2)
           sat_coef3=1d0/(1d0+4.2d0*ratio3)
           sat_coef4=1d0/(1d0+4.2d0*ratio4)
           kparz1=sat_coef1*kparz1
           kparz2=sat_coef2*kparz2
           kparz3=sat_coef3*kparz3
           kparz4=sat_coef4*kparz4
        endif

        if(compute .ne. 3)then
           fz1=kparz1*(bz1*oneminuskperp*(bx1*dTdx1+by1*dTdy1+bz1*dTdz1)/Bnorm2_1+k_perp*dTdz1)
           fz2=kparz2*(bz2*oneminuskperp*(bx2*dTdx2+by2*dTdy2+bz2*dTdz2)/Bnorm2_2+k_perp*dTdz2)
           fz3=kparz3*(bz3*oneminuskperp*(bx3*dTdx3+by3*dTdy3+bz3*dTdz3)/Bnorm2_3+k_perp*dTdz3)
           fz4=kparz4*(bz4*oneminuskperp*(bx4*dTdx4+by4*dTdy4+bz4*dTdz4)/Bnorm2_4+k_perp*dTdz4)
        else
           fz1=kparz1*(bz1*oneminuskperp*(bx1+by1+bz1)/Bnorm2_1+k_perp)
           fz2=kparz2*(bz2*oneminuskperp*(bx2+by2+bz2)/Bnorm2_2+k_perp)
           fz3=kparz3*(bz3*oneminuskperp*(bx3+by3+bz3)/Bnorm2_3+k_perp)
           fz4=kparz4*(bz4*oneminuskperp*(bx4+by4+bz4)/Bnorm2_4+k_perp)
        end if
        fz=0.25d0*(fz1+fz2+fz3+fz4)
#endif
        myflux(l,i,j,k)=fz*dt/dx
     enddo
  enddo
  enddo
  enddo
  endif

  if (slopelim_cond)then
     ! TODO
  endif
endif

end subroutine cmpZheatflx
!###########################################################
!###########################################################
!###########################################################
!###########################################################
