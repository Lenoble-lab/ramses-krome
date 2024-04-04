!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft,
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft,
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
#if NENER>0
  integer::irad
#endif
#if NVAR>8+NENER
  integer::ivar
#endif
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables

  if(napot.gt.0) then
     call condinit_napot(x,u,dx,nn)
     return
  endif
  
  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do irad=1,nener
     u(1:nn,8+irad)=q(1:nn,8+irad)/(gamma_rad(irad)-1.0d0)
     u(1:nn,5)=u(1:nn,5)+u(1:nn,8+irad)
  enddo
#endif
#if NVAR>8+NENER
  ! passive scalars
  do ivar=9+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy=0.,zz=0.,vx,vy,vz,aa,twopi
!!$  real(dp)::rr,tt,omega
  
  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0+0.*t
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)+0.*dx
#if NDIM > 1
     yy=x(i,2)
#endif
#if NDIM > 2
     zz=x(i,3)
#endif
     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega

     v(i,1)=vx
#if NDIM > 1
     v(i,2)=vy
#endif
#if NDIM > 2
     v(i,3)=vz
#endif
  end do


end subroutine velana

!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit_napot(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use hydro_commons,only:apotx,apoty,apotz,napot
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x,xx ! Cell center position.
  integer::i,ii,jj,kk,inode,jnode,knode
  real(dp)::c00,c01,c10,c11,c0,c1
  real(dp)::xnode,ynode,znode,dx1,dxrand,dxx,dyy,dzz
  real(dp)::divB,coef,smallbox_l,dx_loc
  integer::nap1,iloop,jloop,kloop,i1,i2,j1,j2,k1,k2,napot_loc
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp),dimension(2,2)::axnode,aynode,aznode

  dx1=dx/boxlen
  xx=x/boxlen
  if(napot>napotsmall)then
     napot_loc=napotsmall
  else
     napot_loc=napot
  endif
  smallbox_l=dble(napot_loc)/dble(napot)
  nap1=napot_loc+1
  dxrand=1d0/dble(napot_loc)
  dx_loc=dx1/smallbox_l

  do i=1,nn
     xx(i,1)=(xx(i,1)-int(xx(i,1)/smallbox_l)*smallbox_l)/smallbox_l
     xx(i,2)=(xx(i,2)-int(xx(i,2)/smallbox_l)*smallbox_l)/smallbox_l
     xx(i,3)=(xx(i,3)-int(xx(i,3)/smallbox_l)*smallbox_l)/smallbox_l
     ii=int(xx(i,1)/dxrand)+1
     xnode=(ii-1d0)*dxrand
     i1=int((xx(i,1)-0.5d0*dx_loc)/dxrand)+1
     i2=int((xx(i,1)+0.5d0*dx_loc)/dxrand)+1
     jj=int(xx(i,2)/dxrand)+1
     ynode=(jj-1d0)*dxrand
     j1=int((xx(i,2)-0.5d0*dx_loc)/dxrand)+1
     j2=int((xx(i,2)+0.5d0*dx_loc)/dxrand)+1
#if NDIM>2
     kk=int(xx(i,3)/dxrand)+1
     znode=(kk-1d0)*dxrand
     k1=int((xx(i,3)-0.5d0*dx_loc)/dxrand)+1
     k2=int((xx(i,3)+0.5d0*dx_loc)/dxrand)+1
#else
     kk=1
     k1=1
     k2=1
#endif

     if((i2-i1).gt.1)then
        axnode=0.0d0
        aynode=0.0d0
        aznode=0.0d0
        do iloop=i1,i2
           if(iloop==i1.or.iloop.eq.i2)then
              coef=0.5d0
           else
              coef=1.0d0
           endif
           axnode(1,1)=axnode(1,1)+apotx(icell_apot(iloop,j1,k1,nap1))*coef
           axnode(1,2)=axnode(1,2)+apotx(icell_apot(iloop,j1,k2,nap1))*coef
           axnode(2,1)=axnode(2,1)+apotx(icell_apot(iloop,j2,k1,nap1))*coef
           axnode(2,2)=axnode(2,2)+apotx(icell_apot(iloop,j2,k2,nap1))*coef
        enddo
        axnode=axnode/dble(i2-i1)
        do jloop=j1,j2
           if(jloop==j1.or.jloop.eq.j2)then
              coef=0.5d0
           else
              coef=1.0d0
           endif
           aynode(1,1)=aynode(1,1)+apoty(icell_apot(i1,jloop,k1,nap1))*coef
           aynode(1,2)=aynode(1,2)+apoty(icell_apot(i1,jloop,k2,nap1))*coef
           aynode(2,1)=aynode(2,1)+apoty(icell_apot(i2,jloop,k1,nap1))*coef
           aynode(2,2)=aynode(2,2)+apoty(icell_apot(i2,jloop,k2,nap1))*coef
        enddo
        aynode=aynode/dble(j2-j1)
#if NDIM > 2
        do kloop=k1,k2
           if(kloop==k1.or.kloop.eq.k2)then
              coef=0.5d0
           else
              coef=1.0d0
           endif
           aznode(1,1)=aznode(1,1)+apotz(icell_apot(i1,j1,kloop,nap1))*coef
           aznode(1,2)=aznode(1,2)+apotz(icell_apot(i1,j2,kloop,nap1))*coef
           aznode(2,1)=aznode(2,1)+apotz(icell_apot(i2,j1,kloop,nap1))*coef
           aznode(2,2)=aznode(2,2)+apotz(icell_apot(i2,j2,kloop,nap1))*coef
        enddo
        aznode=aznode/dble(k2-k1)
#else
        aznode(1,1)=apotz(icell_apot(i1,j1,1,nap1))
        aznode(1,2)=apotz(icell_apot(i1,j2,1,nap1))
        aznode(2,1)=apotz(icell_apot(i2,j1,1,nap1))
        aznode(2,2)=apotz(icell_apot(i2,j2,1,nap1))
#endif
     else
        ! Interpolate Ax at (i,j-1/2,k-1/2), (i,j+1/2,k-1/2),
        ! (i,j-1/2,k+1/2), (i,j+1/2,k+1/2) in that order
        dxx=(xx(i,1)-xnode)/dxrand
        c00=apotx(icell_apot(ii  ,jj  ,kk  ,nap1))*(1d0-dxx)+apotx(icell_apot(ii+1,jj  ,kk  ,nap1))*dxx
        c01=apotx(icell_apot(ii  ,jj+1,kk  ,nap1))*(1d0-dxx)+apotx(icell_apot(ii+1,jj+1,kk  ,nap1))*dxx
        c10=apotx(icell_apot(ii  ,jj  ,kk+1,nap1))*(1d0-dxx)+apotx(icell_apot(ii+1,jj  ,kk+1,nap1))*dxx
        c11=apotx(icell_apot(ii  ,jj+1,kk+1,nap1))*(1d0-dxx)+apotx(icell_apot(ii+1,jj+1,kk+1,nap1))*dxx
        do knode=1,2
#if NDIM>2
           dzz=(xx(i,3)+(dble(knode)-1.5d0)*dx_loc-znode)/dxrand
#else
           dzz=0.0d0
#endif
           c0=c00*(1d0-dzz)+c10*dzz
           c1=c01*(1d0-dzz)+c11*dzz
           do jnode=1,2
              dyy=(xx(i,2)+(dble(jnode)-1.5d0)*dx_loc-ynode)/dxrand
              axnode(jnode,knode) = c0*(1d0-dyy)+c1*dyy
           enddo
        enddo
        ! Interpolate Ay at (i-1/2,j,k-1/2), (i+1/2,j,k-1/2),
        ! (i-1/2,j,k+1/2), (i+1/2,j,k+1/2) in that order
        dyy=(xx(i,2)-ynode)/dxrand
        c00=apoty(icell_apot(ii  ,jj  ,kk  ,nap1))*(1d0-dyy)+apoty(icell_apot(ii  ,jj+1,kk  ,nap1))*dyy
        c01=apoty(icell_apot(ii+1,jj  ,kk  ,nap1))*(1d0-dyy)+apoty(icell_apot(ii+1,jj+1,kk  ,nap1))*dyy
        c10=apoty(icell_apot(ii  ,jj  ,kk+1,nap1))*(1d0-dyy)+apoty(icell_apot(ii  ,jj+1,kk+1,nap1))*dyy
        c11=apoty(icell_apot(ii+1,jj  ,kk+1,nap1))*(1d0-dyy)+apoty(icell_apot(ii+1,jj+1,kk+1,nap1))*dyy
        do knode=1,2
#if NDIM>2
           dzz=(xx(i,3)+(dble(knode)-1.5d0)*dx_loc-znode)/dxrand
#else
           dzz=0.0d0
#endif
           c0=c00*(1d0-dzz)+c10*dzz
           c1=c01*(1d0-dzz)+c11*dzz
           do inode=1,2
              dxx=(xx(i,1)+(dble(inode)-1.5d0)*dx_loc-xnode)/dxrand
              aynode(inode,knode) = c0*(1d0-dxx)+c1*dxx
           enddo
        enddo
        ! Interpolate Az at (i-1/2,j-1/2,k), (i+1/2,j-1/2,k), 
        ! (i-1/2,j+1/2,k), (i+1/2,j+1/2,k) in that order
#if NDIM>2
        dzz=(xx(i,3)-znode)/dxrand
#else
        dzz=0.0d0
#endif
        c00=apotz(icell_apot(ii  ,jj  ,kk  ,nap1))*(1d0-dzz)+apotz(icell_apot(ii  ,jj  ,kk+1,nap1))*dzz
        c01=apotz(icell_apot(ii+1,jj  ,kk  ,nap1))*(1d0-dzz)+apotz(icell_apot(ii+1,jj  ,kk+1,nap1))*dzz
        c10=apotz(icell_apot(ii  ,jj+1,kk  ,nap1))*(1d0-dzz)+apotz(icell_apot(ii  ,jj+1,kk+1,nap1))*dzz
        c11=apotz(icell_apot(ii+1,jj+1,kk  ,nap1))*(1d0-dzz)+apotz(icell_apot(ii+1,jj+1,kk+1,nap1))*dzz
        do jnode=1,2
           dyy=(xx(i,2)+(dble(jnode)-1.5d0)*dx_loc-ynode)/dxrand
           c0=c00*(1d0-dyy)+c10*dyy
           c1=c01*(1d0-dyy)+c11*dyy
           do inode=1,2
              dxx=(xx(i,1)+(dble(inode)-1.5d0)*dx_loc-xnode)/dxrand
              aznode(inode,jnode) = c0*(1d0-dxx)+c1*dxx
           enddo
        enddo
     endif

     ! B = rot A
     q(i,6     )=( (aznode(1,2)-aznode(1,1))-(aynode(1,2)-aynode(1,1)) )/dx_loc
     q(i,nvar+1)=( (aznode(2,2)-aznode(2,1))-(aynode(2,2)-aynode(2,1)) )/dx_loc
     q(i,7     )=( (axnode(1,2)-axnode(1,1))-(aznode(2,1)-aznode(1,1)) )/dx_loc
     q(i,nvar+2)=( (axnode(2,2)-axnode(2,1))-(aznode(2,2)-aznode(1,2)) )/dx_loc
     q(i,8     )=( (aynode(2,1)-aynode(1,1))-(axnode(2,1)-axnode(1,1)) )/dx_loc
     q(i,nvar+3)=( (aynode(2,2)-aynode(1,2))-(axnode(2,2)-axnode(1,2)) )/dx_loc

     divB=q(i,nvar+1)-q(i,6)+q(i,nvar+2)-q(i,7)+q(i,nvar+3)-q(i,8)
     if(divB.gt. 1d-10)write(*,*)'divB=',divB
  enddo

  ! Rescale B-field by the user-defined magnitude
  u(1:nn,6:8)=q(1:nn,6:8)*B_ave
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)*B_ave

contains

  function icell_apot(i,j,k,nx)
    implicit none
    integer::i,j,k,nx
    integer::icell_apot

    icell_apot=i+(j-1)*nx+(k-1)*nx*nx
  end function icell_apot


end subroutine condinit_napot
