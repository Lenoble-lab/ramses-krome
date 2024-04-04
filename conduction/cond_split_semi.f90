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
subroutine cond_split_semi(uin,flux,dx,dt,ngrid,compute,fdx)
  use amr_parameters
  use const             
  use hydro_parameters
!  use radiation_parameters
  implicit none 

  integer ::ngrid,compute
  real(dp)::dx,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::fdx

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::Xflux,Yflux,Zflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3),save::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::kspitzer
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::Temp,Temp2,dens
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ffdx

  ! Local scalar variables
  integer::i,j,k,l
  integer::ilo,ihi,jlo,jhi,klo,khi

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::Tp_loc

  real(dp)::scale_jsat
  integer ::ind_Teold

  ind_Teold = 3 ! Stored in a different index than in other routines (nvar+3)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_jsat=scale_d*scale_v**2d0/scale_t
  
  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  do k = ku1, ku2
  do j = ju1, ju2
  do i = iu1, iu2
     do l = 1, ngrid

        Tp_loc = uin(l,i,j,k,ind_Teold)

        dens(l,i,j,k)=uin(l,i,j,k,1)*scale_d
        
        kspitzer(l,i,j,k)=uin(l,i,j,k,4)
        bf(l,i,j,k,1)=uin(l,i,j,k,6)
        bf(l,i,j,k,2)=uin(l,i,j,k,7)
        bf(l,i,j,k,3)=uin(l,i,j,k,8)
        
        ffdx(l,i,j,k)=fdx(l,i,j,k)

        if(compute==0)Temp(l,i,j,k)=Tp_loc
        if(compute==1)Temp(l,i,j,k)=Tp_loc
        if(compute==2)Temp(l,i,j,k)=uin(l,i,j,k,2)

        Temp2(l,i,j,k)=Tp_loc
     end do
  enddo
  enddo
  enddo
  
  if(compute==0)then

     !================================
     ! Do the explicit transverse flux
     !================================

     ! Compute the heat flux in X direction  
     call cmpXheatflx_trans(Temp,bf,kspitzer,Xflux,dx,dt,ngrid)
     do k=klo,khi
        do j=jlo,jhi
           do i=if1,if2
              do l = 1, ngrid
                 flux(l,i,j,k,5,1)=-Xflux(l,i,j,k)
              enddo
           enddo
        enddo
     enddo
#if NDIM>1
     ! Compute the heat flux in Y direction
     call cmpYheatflx_trans(Temp,bf,kspitzer,Yflux,dx,dt,ngrid)
     do k=klo,khi
        do j=jf1,jf2
           do i=ilo,ihi
              do l = 1, ngrid
                 flux(l,i,j,k,5,2)=-Yflux(l,i,j,k)
              enddo
           enddo
        enddo
     enddo
#endif
#if NDIM>2
     ! Compute the heat flux in Z direction
     call cmpZheatflx_trans(Temp,bf,kspitzer,Zflux,dx,dt,ngrid)
     do k=kf1,kf2
        do j=jlo,jhi
           do i=ilo,ihi
              do l = 1, ngrid
                 flux(l,i,j,k,5,3)=-Zflux(l,i,j,k)
              enddo
           enddo
        enddo
     enddo
#endif
     
  else

     !============================
     ! Do the implicit normal flux
     !============================

     ! Compute the heat flux in X direction  
     call cmpXheatflx_norm(Temp,bf,kspitzer,Xflux,dx,dt,ngrid,compute,ffdx)
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
     call cmpYheatflx_norm(Temp,bf,kspitzer,Yflux,dx,dt,ngrid,compute)
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
     call cmpZheatflx_norm(Temp,bf,kspitzer,Zflux,dx,dt,ngrid,compute)
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

  endif

end subroutine cond_split_semi
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmpXheatflx_trans(Temp,bf,kspitzer,myflux,dx,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::myflux

  ! Primitive variables
  !real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::xloc
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::kspitzer
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Temp

  real(dp)::fx,oneovertwodx,oneoverfourdx 
  real(dp)::Bnorm2_1,Bnorm2_2,Bnorm2_3,Bnorm2_4
  real(dp)::dTdy1,dTdy2,dTdy3,dTdy4
  real(dp)::dTdz1,dTdz2,dTdz3,dTdz4
  real(dp)::bx1,bx2,bx3,bx4
  real(dp)::by1,by2,by3,by4
  real(dp)::bz1,bz2,bz3,bz4
  real(dp)::fx1,fx2,fx3,fx4
  real(dp)::kparx1,kparx2,kparx3,kparx4,oneminuskperp
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

  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
#if NDIM==1
!!$        1-------

        fx=0.0d0

#endif
#if NDIM==2
!!$        2-------
!!$        |      |
!!$        |      |
!!$        |      |
!!$        1-------
        
        dTdy1=(Temp(l,i  ,j  ,k)+Temp(l,i-1,j  ,k) &
             - Temp(l,i  ,j-1,k)-Temp(l,i-1,j-1,k))*oneovertwodx
        dTdy2=(Temp(l,i  ,j+1,k)+Temp(l,i-1,j+1,k) &
             - Temp(l,i  ,j  ,k)-Temp(l,i-1,j  ,k))*oneovertwodx
        
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

        fx1=kparx1*bx1*oneminuskperp*by1*dTdy1
        fx2=kparx2*bx2*oneminuskperp*by2*dTdy2
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

        fx1=kparx1*bx1*oneminuskperp*(by1*dTdy1+bz1*dTdz1)/Bnorm2_1
        fx2=kparx2*bx2*oneminuskperp*(by2*dTdy2+bz2*dTdz2)/Bnorm2_2
        fx3=kparx3*bx3*oneminuskperp*(by3*dTdy3+bz3*dTdz3)/Bnorm2_3
        fx4=kparx4*bx4*oneminuskperp*(by4*dTdy4+bz4*dTdz4)/Bnorm2_4
        fx=0.25d0*(fx1+fx2+fx3+fx4)
#endif
        myflux(l,i,j,k)=fx*dt/dx
     enddo
  enddo
  enddo
  enddo

end subroutine cmpXheatflx_trans
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmpYheatflx_trans(Temp,bf,kspitzer,myflux,dx,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::kspitzer
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Temp

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
        
        dTdx1=(Temp(l,i  ,j  ,k)+Temp(l,i  ,j-1,k) &
             - Temp(l,i-1,j  ,k)-Temp(l,i-1,j-1,k))*oneovertwodx
        dTdx2=(Temp(l,i+1,j  ,k)+Temp(l,i+1,j-1,k) &
             - Temp(l,i  ,j  ,k)-Temp(l,i  ,j-1,k))*oneovertwodx        

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
             +      kspitzer(l,i-1,j  ,k)+kspitzer(l,i-1,j-1,k))
        kpary2=0.25d0*(kspitzer(l,i  ,j  ,k)+kspitzer(l,i+1,j  ,k) &
             +      kspitzer(l,i  ,j-1,k)+kspitzer(l,i+1,j-1,k))

        fy1=kpary1*by1*oneminuskperp*bx1*dTdx1
        fy2=kpary2*by2*oneminuskperp*bx2*dTdx2
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

        fy1=kpary1*by1*oneminuskperp*(bx1*dTdx1+bz1*dTdz1)/Bnorm2_1
        fy2=kpary2*by2*oneminuskperp*(bx2*dTdx2+bz2*dTdz2)/Bnorm2_2
        fy3=kpary3*by3*oneminuskperp*(bx3*dTdx3+bz3*dTdz3)/Bnorm2_3
        fy4=kpary4*by4*oneminuskperp*(bx4*dTdx4+bz4*dTdz4)/Bnorm2_4
        fy=0.25d0*(fy1+fy2+fy3+fy4)
#endif
        myflux(l,i,j,k)=fy*dt/dx
     enddo
  enddo
  enddo
  enddo

end subroutine cmpYheatflx_trans
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmpZheatflx_trans(Temp,bf,kspitzer,myflux,dx,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::kspitzer
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Temp

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

  ! Local scalar variables
  integer::i,j,k,l
  integer::ilo,ihi,jlo,jhi

  oneovertwodx =0.50d0/dx
  oneoverfourdx=0.25d0/dx
  oneminuskperp=1.0d0-k_perp

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)

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
             +      kspitzer(l,i  ,j-1,k-1)+kspitzer(l,i-1,j  ,k  )+kspitzer(l,i-1,j-1,k  ) &
             +      kspitzer(l,i-1,j  ,k-1)+kspitzer(l,i-1,j-1,k-1))
        kparz2=0.125d0*(kspitzer(l,i+1,j  ,k  )+kspitzer(l,i+1,j-1,k  )+kspitzer(l,i+1,j  ,k-1) &
             +      kspitzer(l,i+1,j-1,k-1)+kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j-1,k  ) &
             +      kspitzer(l,i  ,j  ,k-1)+kspitzer(l,i  ,j-1,k-1))
        kparz3=0.125d0*(kspitzer(l,i  ,j+1,k  )+kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j+1,k-1) &
             +      kspitzer(l,i  ,j  ,k-1)+kspitzer(l,i-1,j+1,k  )+kspitzer(l,i-1,j  ,k  ) &
             +      kspitzer(l,i-1,j+1,k-1)+kspitzer(l,i-1,j  ,k-1))
        kparz4=0.125d0*(kspitzer(l,i+1,j+1,k  )+kspitzer(l,i+1,j  ,k  )+kspitzer(l,i+1,j+1,k-1) &
             +      kspitzer(l,i+1,j  ,k-1)+kspitzer(l,i  ,j+1,k  )+kspitzer(l,i  ,j  ,k  ) &
             +      kspitzer(l,i  ,j+1,k-1)+kspitzer(l,i  ,j  ,k-1))

        fz1=kparz1*bz1*oneminuskperp*(bx1*dTdx1+by1*dTdy1)/Bnorm2_1
        fz2=kparz2*bz2*oneminuskperp*(bx2*dTdx2+by2*dTdy2)/Bnorm2_2
        fz3=kparz3*bz3*oneminuskperp*(bx3*dTdx3+by3*dTdy3)/Bnorm2_3
        fz4=kparz4*bz4*oneminuskperp*(bx4*dTdx4+by4*dTdy4)/Bnorm2_4
        fz=0.25d0*(fz1+fz2+fz3+fz4)
#endif
        myflux(l,i,j,k)=fz*dt/dx
     enddo
  enddo
  enddo
  enddo

end subroutine cmpZheatflx_trans
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpXheatflx_norm(Temp,bf,kspitzer,myflux,dx,dt,ngrid,compute,ffdx)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid,compute
  real(dp)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::myflux

  ! Primitive variables
  !real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::xloc
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::kspitzer
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Temp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ffdx

  real(dp)::famr,fx,oneoverdx,oneovertwodx,oneoverfourdx 
  real(dp)::Bnorm2_1,Bnorm2_2,Bnorm2_3,Bnorm2_4
  real(dp)::dTdx1
  real(dp)::bx1,bx2,bx3,bx4,bx
  real(dp)::by1,by2,by3,by4
  real(dp)::bz1,bz2,bz3,bz4
  real(dp)::kparx1,oneminuskperp
#if NDIM==2
  real(dp)::Bnorm
#endif
  ! Local scalar variables
  integer::i,j,k,l
  integer::jlo,jhi,klo,khi

  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  oneoverdx    =1.00d0/dx
  oneovertwodx =0.50d0/dx
  oneoverfourdx=0.25d0/dx
  oneminuskperp=1.0d0-k_perp

  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        famr=max(ffdx(l,i,j,k),ffdx(l,i-1,j,k))
#if NDIM==1
!!$        1-------

        if(compute .ne. 3) dTdx1=(Temp(l,i,j,k)-Temp(l,i-1,j,k))/dx
        kparx1=0.5d0*(kspitzer(l,i,j,k)+kspitzer(l,i-1,j,k)) ! arithmetic mean
        if(compute .ne. 3)then
           fx=kparx1*dTdx1
        else 
           fx=kparx1/dx
        endif
        fx=fx/famr

#endif
#if NDIM==2
!!$        --------
!!$        |      |
!!$        1      |
!!$        |      |
!!$        --------
        
        if(compute .ne. 3)dTdx1=(Temp(l,i  ,j  ,k)-Temp(l,i-1,j  ,k))*oneoverdx
!!$        if(compute.ne.3)then
!!$           dTdx1=(Temp(l,i  ,j  ,k)+Temp(l,i  ,j-1,k) &
!!$                - Temp(l,i-1,j  ,k)-Temp(l,i-1,j-1,k))*oneovertwodx
!!$           dTdx2=(Temp(l,i  ,j+1,k)+Temp(l,i  ,j  ,k) &
!!$                - Temp(l,i-1,j+1,k)-Temp(l,i-1,j  ,k))*oneovertwodx
!!$        endif
        
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

        kparx1=0.5d0*(kspitzer(l,i  ,j  ,k)+kspitzer(l,i-1,j  ,k))
!!$        kparx1=0.25d0*(kspitzer(l,i  ,j  ,k)+kspitzer(l,i  ,j-1,k) &
!!$             +      kspitzer(l,i-1,j  ,k)+kspitzer(l,i-1,j-1,k))
!!$        kparx2=0.25d0*(kspitzer(l,i  ,j+1,k)+kspitzer(l,i  ,j  ,k) &
!!$             +      kspitzer(l,i-1,j+1,k)+kspitzer(l,i-1,j  ,k))
!!$
        if(compute .ne. 3)then   
           fx=kparx1*(bx1*oneminuskperp*(bx1*dTdx1)+k_perp*dTdx1)
        else ! Preconditionner
           fx=kparx1*(bx1*oneminuskperp*(bx1)+k_perp)
        end if
!!$        if(compute .ne. 3)then   
!!$           fx1=kparx1*(bx1*oneminuskperp*(bx1*dTdx1)+k_perp*dTdx1)
!!$           fx2=kparx2*(bx2*oneminuskperp*(bx2*dTdx2)+k_perp*dTdx2)
!!$        else ! Preconditionner
!!$           fx1=kparx1*(bx1*oneminuskperp*(bx1)+k_perp)
!!$           fx2=kparx2*(bx2*oneminuskperp*(bx2)+k_perp)
!!$        end if
!!$        fx=0.5d0*(fx1+fx2)
#endif
#if NDIM==3
!!$          ---------
!!$         / |      /|
!!$        /  |     / |
!!$        -1------   |
!!$        |  |   |   |
!!$        | /    |  /
!!$        |/     | /
!!$        --------
        
        if(compute .ne. 3)dTdx1=(Temp(l,i  ,j  ,k  )-Temp(l,i-1,j  ,k  ))*oneoverdx

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

        bx1=bx1/sqrt(Bnorm2_1)
        bx2=bx2/sqrt(Bnorm2_2)
        bx3=bx3/sqrt(Bnorm2_3)
        bx4=bx4/sqrt(Bnorm2_4)
        bx=0.25d0*(bx1+bx2+bx3+bx4)

        kparx1=0.5d0*(kspitzer(l,i  ,j  ,k  )+kspitzer(l,i-1,j  ,k  ))

        if(compute .ne. 3)then   
           fx=kparx1*(bx*oneminuskperp*(bx*dTdx1)+k_perp*dTdx1)
        else ! Preconditionner
           fx=kparx1*(bx*oneminuskperp*(bx)+k_perp)
        end if
#endif
        myflux(l,i,j,k)=fx*dt/dx
     enddo
  enddo
  enddo
  enddo

end subroutine cmpXheatflx_norm
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmpYheatflx_norm(Temp,bf,kspitzer,myflux,dx,dt,ngrid,compute)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid,compute
  real(dp)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::kspitzer
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Temp

  real(dp)::fy,oneoverdx,oneovertwodx,oneoverfourdx
  real(dp)::Bnorm2_1,Bnorm2_2,Bnorm2_3,Bnorm2_4
  real(dp)::dTdy1
  real(dp)::bx1,bx2,bx3,bx4
  real(dp)::by1,by2,by3,by4,by
  real(dp)::bz1,bz2,bz3,bz4
  real(dp)::kpary1,oneminuskperp
#if NDIM==2
  real(dp)::Bnorm
#endif
  ! Local scalar variables
  integer::i,j,k,l
  integer::ilo,ihi,klo,khi

  oneoverdx    =1.00d0/dx
  oneovertwodx =0.50d0/dx
  oneoverfourdx=0.25d0/dx
  oneminuskperp=1.0d0-k_perp

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
#if NDIM==2
!!$        --------
!!$        |      |
!!$        |      |
!!$        |      |
!!$        ----1---
        
        if(compute .ne. 3)dTdy1=(Temp(l,i  ,j  ,k)-Temp(l,i  ,j-1,k))*oneoverdx
!!$        if(compute .ne. 3)then
!!$           dTdy1=(Temp(l,i  ,j  ,k)+Temp(l,i-1,j  ,k) &
!!$                - Temp(l,i  ,j-1,k)-Temp(l,i-1,j-1,k))*oneovertwodx
!!$           dTdy2=(Temp(l,i  ,j  ,k)+Temp(l,i+1,j  ,k) &
!!$                - Temp(l,i  ,j-1,k)-Temp(l,i+1,j-1,k))*oneovertwodx
!!$        end if

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

        kpary1=0.5d0*(kspitzer(l,i  ,j  ,k)+kspitzer(l,i  ,j-1,k))
!!$        kpary1=0.25d0*(kspitzer(l,i  ,j  ,k)+kspitzer(l,i  ,j-1,k) &
!!$             +      kspitzer(l,i-1,j  ,k)+kspitzer(l,i-1,j-1,k))
!!$        kpary2=0.25d0*(kspitzer(l,i  ,j  ,k)+kspitzer(l,i+1,j  ,k) &
!!$             +      kspitzer(l,i  ,j-1,k)+kspitzer(l,i+1,j-1,k))
!!$
        if(compute .ne. 3)then   
           fy=kpary1*(by1*oneminuskperp*(by1*dTdy1)+k_perp*dTdy1)
        else ! Preconditionner
           fy=kpary1*(by1*oneminuskperp*(by1)+k_perp)
        end if
!!$        if(compute .ne. 3)then   
!!$           fy1=kpary1*(by1*oneminuskperp*(by1*dTdy1)+k_perp*dTdy1)
!!$           fy2=kpary2*(by2*oneminuskperp*(by2*dTdy2)+k_perp*dTdy2)
!!$        else ! Preconditionner
!!$           fy1=kpary1*(by1*oneminuskperp*(by1)+k_perp)
!!$           fy2=kpary2*(by2*oneminuskperp*(by2)+k_perp)
!!$        end if
!!$        fy=0.5d0*(fy1+fy2)
#endif
#if NDIM==3
!!$          ---------
!!$         / |      /|
!!$        /  |     / |
!!$        --------   |
!!$        |  |   |   |
!!$        | / 1  |  /
!!$        |/     | /
!!$        --------
        
        ! Centered symmetric scheme
        if(compute .ne. 3)dTdy1=(Temp(l,i  ,j  ,k  )-Temp(l,i  ,j-1,k  ))*oneoverdx

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
        by1=by1/sqrt(Bnorm2_1)
        by2=by2/sqrt(Bnorm2_2)
        by3=by3/sqrt(Bnorm2_3)
        by4=by4/sqrt(Bnorm2_4)
        by=0.25d0*(by1+by2+by3+by4)

        kpary1=0.5d0*(kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j-1,k  ))

        if(compute .ne. 3)then
           fy=kpary1*(by*oneminuskperp*by*dTdy1+k_perp*dTdy1)
        else
           fy=kpary1*(by*oneminuskperp*by+k_perp)
        end if
#endif
        myflux(l,i,j,k)=fy*dt/dx
     enddo
  enddo
  enddo
  enddo

end subroutine cmpYheatflx_norm
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmpZheatflx_norm(Temp,bf,kspitzer,myflux,dx,dt,ngrid,compute)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid,compute
  real(dp)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::kspitzer
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Temp

  real(dp)::fz,oneoverdx,oneovertwodx,oneoverfourdx
  real(dp)::Bnorm2_1,Bnorm2_2,Bnorm2_3,Bnorm2_4
  real(dp)::dTdz1
  real(dp)::bx1,bx2,bx3,bx4
  real(dp)::by1,by2,by3,by4
  real(dp)::bz1,bz2,bz3,bz4,bz
  real(dp)::kparz1,oneminuskperp

  ! Local scalar variables
  integer::i,j,k,l
  integer::ilo,ihi,jlo,jhi

  oneoverdx    =1.00d0/dx
  oneovertwodx =0.50d0/dx
  oneoverfourdx=0.25d0/dx
  oneminuskperp=1.0d0-k_perp

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)

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
!!$        | /    |  /
!!$        |/   1 | /
!!$        --------
        
        ! Centered symmetric scheme
        if(compute .ne. 3)dTdz1=(Temp(l,i  ,j  ,k  )-Temp(l,i  ,j  ,k-1))*oneoverdx
        
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
        bz1=bz1/sqrt(Bnorm2_1)
        bz2=bz2/sqrt(Bnorm2_2)
        bz3=bz3/sqrt(Bnorm2_3)
        bz4=bz4/sqrt(Bnorm2_4)
        bz=0.25d0*(bz1+bz2+bz3+bz4)

        kparz1=0.5d0*(kspitzer(l,i  ,j  ,k  )+kspitzer(l,i  ,j  ,k-1))

        if(compute .ne. 3)then        
           fz=kparz1*(bz*oneminuskperp*bz*dTdz1+k_perp*dTdz1)
        else
           fz=kparz1*(bz*oneminuskperp*bz+k_perp)
        end if
#endif
        myflux(l,i,j,k)=fz*dt/dx
     enddo
  enddo
  enddo
  enddo

end subroutine cmpZheatflx_norm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
