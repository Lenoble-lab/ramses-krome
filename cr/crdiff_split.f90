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
!  dx          => (const)  (dx,dy,dz)
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
!  Updated for slope limiter (06/2019 Y. Dubois)
! ----------------------------------------------------------------
subroutine crdiff_split(uin,flux,dx,dt,ngrid,compute,fdx,igroup)
  use amr_parameters
  use const          
  use hydro_parameters
  implicit none 

  integer ::ngrid,compute,igroup
  real(dp)::dx,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::fdx

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::Xflux,Yflux,Zflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3),save::bf  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::Ecr,Ecr0,Dpara,kperp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ffdx

  ! Local scalar variables
  integer::i,j,k,l
  integer::ilo,ihi,jlo,jhi,klo,khi
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_kappa,kpar

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kappa = scale_l**2/scale_t
  kpar=Dcr/scale_kappa

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

 ! kappa_coef=1.84d-5/coulomb_log*scale_TK**2.5d0/scale_kspitzer * f_spitzer

  do k = ku1, ku2
  do j = ju1, ju2
  do i = iu1, iu2
     do l = 1, ngrid

        bf(l,i,j,k,1)=uin(l,i,j,k,6)
        bf(l,i,j,k,2)=uin(l,i,j,k,7)
        bf(l,i,j,k,3)=uin(l,i,j,k,8)
        
        ffdx(l,i,j,k)=fdx(l,i,j,k)

        if(compute==0)Ecr(l,i,j,k)=uin(l,i,j,k,igroup)
        if(compute==1)Ecr(l,i,j,k)=uin(l,i,j,k,igroup)
        if(compute==2)Ecr(l,i,j,k)=uin(l,i,j,k,1)
        if(compute==4)Ecr(l,i,j,k)=uin(l,i,j,k,2)
        Dpara(l,i,j,k)=uin(l,i,j,k,3)
        kperp(l,i,j,k)=max(uin(l,i,j,k,4),k_perp)
        Ecr0 (l,i,j,k)=uin(l,i,j,k,igroup) ! store the initial value of the energy
                                           ! must be available at all iterations for slope limiter

     end do
  enddo
  enddo
  enddo

 ! Compute the heat flux in X direction
  call cmpXcrflx(Ecr,Ecr0,bf,Xflux,dx,dt,ngrid,compute,ffdx,kpar,Dpara,kperp)
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
  call cmpYcrflx(Ecr,Ecr0,bf,Yflux,dx,dt,ngrid,compute,ffdx,kpar,Dpara,kperp)
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
  call cmpZcrflx(Ecr,Ecr0,bf,Zflux,dx,dt,ngrid,compute,ffdx,kpar,Dpara,kperp)
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
end subroutine crdiff_split
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmpXcrflx(Temp,Temp0,bf,myflux,dx,dt,ngrid,compute,ffdx,kpar,Dpara,kperp)
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
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Temp,Temp0,Dpara,kperp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ffdx

  real(dp)::Bnorm,fx,oneoverdx,oneovertwodx,oneoverfourdx
  real(dp)::dTdx1,dTdx2,dTdx3,dTdx4,dx_loc
  real(dp)::dTdy1,dTdy2,dTdy3,dTdy4
  real(dp)::dTdz1,dTdz2,dTdz3,dTdz4
  real(dp)::bx1,bx2,bx3,bx4
  real(dp)::by1,by2,by3,by4
  real(dp)::bz1,bz2,bz3,bz4
  real(dp)::fx1,fx2,fx3,fx4
  real(dp)::kpar,kpar2,oneminuskperp,kparax1,kparax2,kparax3,kparax4
  real(dp)::kperpx1,kperpx2,kperpx3,kperpx4
  real(dp)::oneminuskperpx1,oneminuskperpx2,oneminuskperpx3,oneminuskperpx4

  ! Local scalar variables
  integer::i,j,k,l
  integer::jlo,jhi,klo,khi
  real(dp)::a0,b0,c0,d0,a,b,c,d
  real(dp)::dTdy0,dTdyMinus0,dTdyPlus0,dTdz0,dTdzMinus0,dTdzPlus0
  real(dp)::dTdy,dTdyMinus,dTdyPlus,dTdz,dTdzMinus,dTdzPlus

  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  oneoverdx    =1.00d0/dx
  oneovertwodx =0.50d0/dx
  oneoverfourdx=0.25d0/dx
  oneminuskperp=1.0d0-k_perp

if(isotrope_cond)then
  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        dTdx1  =(Temp(l,i,j,k)-Temp(l,i-1,j,k))/dx
        kpar2=0.0d0
        if(streaming_diffusion)kpar2 = 0.5d0*(Dpara(l,i,j,k)+Dpara(l,i-1,j,k))
        fx    =(kpar+kpar2)*dTdx1
        dx_loc=max(ffdx(l,i,j,k),ffdx(l,i-1,j,k))
        fx=fx/dx_loc
        myflux(l,i,j,k)=fx*dt/dx
        if(compute==3)myflux(l,i,j,k)=(kpar+kpar2)*dt/dx**2
     enddo
  enddo
  enddo
  enddo

else

  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
#if NDIM==1
!!$        1-------
        dTdx1=(Temp(l,i,j,k)-Temp(l,i-1,j,k))/dx
        if(alfven_diff_coeff  )kpar  = 0.5d0*(Dpara(l,i,j,k)+Dpara(l,i-1,j,k))
        kpar2=0.0d0
        if(streaming_diffusion)kpar2 = 0.5d0*(Dpara(l,i,j,k)+Dpara(l,i-1,j,k))
        fx=(kpar+kpar2)*dTdx1
!!$        write(*,'(4i9,1e15.7)')l,i,j,k,kpar+kpar2
        if(compute==3)fx=(kpar+kpar2)/dx 
        dx_loc=max(ffdx(l,i,j,k),ffdx(l,i-1,j,k))
        fx=fx/dx_loc
#endif
#if NDIM==2
!!$        2-------
!!$        |      |
!!$        |      |
!!$        |      |
!!$        1-------
        
        if(compute .ne. 3)then   
           if(do_limiter_aniso)then
              ! Limit the normal x-component
              a0=Temp0(l,i,j  ,k)-Temp0(l,i-1,j  ,k)
              b0=Temp0(l,i,j-1,k)-Temp0(l,i-1,j-1,k)
              a =Temp (l,i,j  ,k)-Temp (l,i-1,j  ,k)
              b =Temp (l,i,j-1,k)-Temp (l,i-1,j-1,k)
              call L2_lim(a0,b0,a,b,alpha_limiter,dTdx1)
              b0=Temp0(l,i,j+1,k)-Temp0(l,i-1,j+1,k)
              b =Temp (l,i,j+1,k)-Temp (l,i-1,j+1,k)
              call L2_lim(a0,b0,a,b,alpha_limiter,dTdx2)
              dTdx1=dTdx1*oneoverdx
              dTdx2=dTdx2*oneoverdx

!!$              ! WARNING : this is a test to force asymetric
!!$              dTdx1=a*oneoverdx
!!$              dTdx2=a*oneoverdx

              if(compute==0)then
                 a0=Temp (l,i  ,j  ,k)-Temp (l,i  ,j-1,k)
                 b0=Temp (l,i  ,j+1,k)-Temp (l,i  ,j  ,k)
                 c0=Temp (l,i-1,j  ,k)-Temp (l,i-1,j-1,k)
                 d0=Temp (l,i-1,j+1,k)-Temp (l,i-1,j  ,k)
              else
                 a0=Temp0(l,i  ,j  ,k)-Temp0(l,i  ,j-1,k)
                 b0=Temp0(l,i  ,j+1,k)-Temp0(l,i  ,j  ,k)
                 c0=Temp0(l,i-1,j  ,k)-Temp0(l,i-1,j-1,k)
                 d0=Temp0(l,i-1,j+1,k)-Temp0(l,i-1,j  ,k)
              endif

              a =Temp (l,i  ,j  ,k)-Temp (l,i  ,j-1,k)
              b =Temp (l,i  ,j+1,k)-Temp (l,i  ,j  ,k)
              c =Temp (l,i-1,j  ,k)-Temp (l,i-1,j-1,k)
              d =Temp (l,i-1,j+1,k)-Temp (l,i-1,j  ,k)
              if    (slope_limiter_aniso == 1)then
                 call MinMod_lim(a0,b0,a,b,dTdyPlus0 ,dTdyPlus )
                 call MinMod_lim(c0,d0,c,d,dTdyMinus0,dTdyMinus)
                 call MinMod_lim(dTdyPlus0,dTdyMinus0,dTdyPlus,dTdyMinus,dTdy0,dTdy)
              !elseif(slope_limiter_aniso == 2)then
              !   call MonCen_lim(a0,b0,a,b,dTdyPlus0 ,dTdyPlus )
              !   call MonCen_lim(c0,d0,c,d,dTdyMinus0,dTdyMinus)
              !   call MonCen_lim(dTdyPlus0,dTdyMinus0,dTdyPlus,dTdyMinus,dTdy0,dTdy)
              !   if(abs(dTdy0).le.0.25d0*(a+b+c+d))write(*,*)'Limiting the slope',indpass(l,i,j,k)
              endif
              dTdy1=dTdy*oneoverdx
              dTdy2=dTdy*oneoverdx
           else
              dTdx1=(Temp(l,i  ,j  ,k)+Temp(l,i  ,j-1,k) &
                   - Temp(l,i-1,j  ,k)-Temp(l,i-1,j-1,k))*oneovertwodx
              dTdx2=(Temp(l,i  ,j+1,k)+Temp(l,i  ,j  ,k) &
                   - Temp(l,i-1,j+1,k)-Temp(l,i-1,j  ,k))*oneovertwodx        

              dTdy1=(Temp(l,i  ,j  ,k)+Temp(l,i-1,j  ,k) &
                   - Temp(l,i  ,j-1,k)-Temp(l,i-1,j-1,k))*oneovertwodx
              dTdy2=(Temp(l,i  ,j+1,k)+Temp(l,i-1,j+1,k) &
                   - Temp(l,i  ,j  ,k)-Temp(l,i-1,j  ,k))*oneovertwodx
           endif
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

        if(alfven_diff_coeff)then
           ! arithmetic mean
           kparax1=0.25d0*(Dpara(l,i  ,j  ,k)+Dpara(l,i  ,j-1,k) &
                +      Dpara(l,i-1,j  ,k)+Dpara(l,i-1,j-1,k))
           kparax2=0.25d0*(Dpara(l,i  ,j+1,k)+Dpara(l,i  ,j  ,k) &
                +      Dpara(l,i-1,j+1,k)+Dpara(l,i-1,j  ,k))
!!$        kparx1=4d0/(Dpara(l,i  ,j  ,k)+Dpara(l,i  ,j-1,k) &
!!$             +      Dpara(l,i-1,j  ,k)+Dpara(l,i-1,j-1,k))
!!$        kparx2=4d0/(Dpara(l,i  ,j+1,k)+Dpara(l,i  ,j  ,k) &
!!$             +      Dpara(l,i-1,j+1,k)+Dpara(l,i-1,j  ,k))
           
           kperpx1=0.25d0*(kperp(l,i  ,j  ,k)+kperp(l,i  ,j-1,k) &
                +      kperp(l,i-1,j  ,k)+kperp(l,i-1,j-1,k))
           kperpx2=0.25d0/(kperp(l,i  ,j+1,k)+kperp(l,i  ,j  ,k) &
                +      kperp(l,i-1,j+1,k)+kperp(l,i-1,j  ,k))
           
           oneminuskperpx1 = 1.0d0-kperpx1
           oneminuskperpx2 = 1.0d0-kperpx2
        else if(streaming_diffusion)then
           kparax1=0.25d0*(Dpara(l,i  ,j  ,k)+Dpara(l,i  ,j-1,k) &
                +      Dpara(l,i-1,j  ,k)+Dpara(l,i-1,j-1,k)) + kpar
           kparax2=0.25d0*(Dpara(l,i  ,j+1,k)+Dpara(l,i  ,j  ,k) &
                +      Dpara(l,i-1,j+1,k)+Dpara(l,i-1,j  ,k)) + kpar
           kperpx1 = k_perp
           kperpx2 = k_perp
           oneminuskperpx1 = oneminuskperp
           oneminuskperpx2 = oneminuskperp
        else
           kparax1 = kpar
           kparax2 = kpar
           kperpx1 = k_perp
           kperpx2 = k_perp
           oneminuskperpx1 = oneminuskperp
           oneminuskperpx2 = oneminuskperp
        endif

        if(compute .ne. 3)then
           fx1=kparax1*(bx1*oneminuskperpx1*(bx1*dTdx1+by1*dTdy1)+kperpx1*dTdx1)
           fx2=kparax2*(bx2*oneminuskperpx2*(bx2*dTdx2+by2*dTdy2)+kperpx2*dTdx2)
        else ! Preconditionner
           fx1=kparax1*(bx1*oneminuskperpx1*(bx1+by1)+kperpx1)
           fx2=kparax2*(bx2*oneminuskperpx2*(bx2+by2)+kperpx2)
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
           if(do_limiter_aniso)then
              ! Limit the normal x-component
              ! T:Top , M:Middle, B:Bottom in k
              ! L:Left, C:Center, R:Right  in j
!!$              ii=i-1
!!$              jj=j-1;kk=k+1;TL=Temp(l,i,jj,kk)-Temp(l,ii,jj,kk)
!!$              jj=j         ;TC=Temp(l,i,jj,kk)-Temp(l,ii,jj,kk)
!!$              jj=j+1       ;TR=Temp(l,i,jj,kk)-Temp(l,ii,jj,kk)
!!$              jj=j-1;kk=k  ;ML=Temp(l,i,jj,kk)-Temp(l,ii,jj,kk)
!!$              jj=j         ;MC=Temp(l,i,jj,kk)-Temp(l,ii,jj,kk)
!!$              jj=j+1       ;MR=Temp(l,i,jj,kk)-Temp(l,ii,jj,kk)
!!$              jj=j-1;kk=k-1;BL=Temp(l,i,jj,kk)-Temp(l,ii,jj,kk)
!!$              jj=j         ;BC=Temp(l,i,jj,kk)-Temp(l,ii,jj,kk)
!!$              jj=j+1       ;BR=Temp(l,i,jj,kk)-Temp(l,ii,jj,kk)
              
              dTdx1=(Temp(l,i,j,k)-Temp(l,i-1,j,k))*oneoverdx
              dTdx2=dTdx1
              dTdx3=dTdx1
              dTdx4=dTdx1

              ! Limit the transverse y-component
              if(compute==0)then
                 a0=Temp (l,i  ,j  ,k)-Temp (l,i  ,j-1,k)
                 b0=Temp (l,i  ,j+1,k)-Temp (l,i  ,j  ,k)
                 c0=Temp (l,i-1,j  ,k)-Temp (l,i-1,j-1,k)
                 d0=Temp (l,i-1,j+1,k)-Temp (l,i-1,j  ,k)
              else
                 a0=Temp0(l,i  ,j  ,k)-Temp0(l,i  ,j-1,k)
                 b0=Temp0(l,i  ,j+1,k)-Temp0(l,i  ,j  ,k)
                 c0=Temp0(l,i-1,j  ,k)-Temp0(l,i-1,j-1,k)
                 d0=Temp0(l,i-1,j+1,k)-Temp0(l,i-1,j  ,k)
              endif
              a =Temp (l,i  ,j  ,k)-Temp (l,i  ,j-1,k)
              b =Temp (l,i  ,j+1,k)-Temp (l,i  ,j  ,k)
              c =Temp (l,i-1,j  ,k)-Temp (l,i-1,j-1,k)
              d =Temp (l,i-1,j+1,k)-Temp (l,i-1,j  ,k)
              if(slope_limiter_aniso == 1)then
                 call MinMod_lim(a0,b0,a,b,dTdyPlus0 ,dTdyPlus )
                 call MinMod_lim(c0,d0,c,d,dTdyMinus0,dTdyMinus)
                 call MinMod_lim(dTdyPlus0,dTdyMinus0,dTdyPlus,dTdyMinus,dTdy0,dTdy)
              !elseif(slope_limiter_aniso == 2)then
              !   call MonCen_lim(as,b,dTdyMinus)
              !   call MonCen_lim(c,d,dTdyPlus )
              !   call MonCen_lim(dTdyMinus,dTdyPlus,dTdy)
              endif
              dTdy1=dTdy*oneoverdx
              dTdy2=dTdy1
              dTdy3=dTdy1
              dTdy4=dTdy1

              ! Limit the transverse z-component
              if(compute==0)then
                 a0=Temp (l,i  ,j,k  )-Temp (l,i  ,j,k-1)
                 b0=Temp (l,i  ,j,k+1)-Temp (l,i  ,j,k  )
                 c0=Temp (l,i-1,j,k  )-Temp (l,i-1,j,k-1)
                 d0=Temp (l,i-1,j,k+1)-Temp (l,i-1,j,k  )
              else
                 a0=Temp0(l,i  ,j,k  )-Temp0(l,i  ,j,k-1)
                 b0=Temp0(l,i  ,j,k+1)-Temp0(l,i  ,j,k  )
                 c0=Temp0(l,i-1,j,k  )-Temp0(l,i-1,j,k-1)
                 d0=Temp0(l,i-1,j,k+1)-Temp0(l,i-1,j,k  )
              endif
              a =Temp (l,i  ,j,k  )-Temp (l,i  ,j,k-1)
              b =Temp (l,i  ,j,k+1)-Temp (l,i  ,j,k  )
              c =Temp (l,i-1,j,k  )-Temp (l,i-1,j,k-1)
              d =Temp (l,i-1,j,k+1)-Temp (l,i-1,j,k  )
              if(slope_limiter_aniso == 1)then
                 call MinMod_lim(a0,b0,a,b,dTdzPlus0 ,dTdzPlus )
                 call MinMod_lim(c0,d0,c,d,dTdzMinus0,dTdzMinus)
                 call MinMod_lim(dTdzPlus0,dTdzMinus0,dTdzPlus,dTdzMinus,dTdz0,dTdz)
              !elseif(slope_limiter_aniso == 2)then
              !   call MonCen_lim(a,b,dTdzMinus)
              !   call MonCen_lim(c,d,dTdzPlus )
              !   call MonCen_lim(dTdzMinus,dTdzPlus,dTdz)
              endif
              dTdz1=dTdz*oneoverdx
              dTdz2=dTdz1
              dTdz3=dTdz1
              dTdz4=dTdz1

           else
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
           endif
        end if

        bx1=0.25d0*(bf(l,i  ,j-1,k-1,1)+bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j-1,k  ,1)+bf(l,i  ,j  ,k  ,1))
        bx2=0.25d0*(bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j+1,k-1,1)+bf(l,i  ,j  ,k  ,1)+bf(l,i  ,j+1,k  ,1))
        bx3=0.25d0*(bf(l,i  ,j-1,k  ,1)+bf(l,i  ,j  ,k  ,1)+bf(l,i  ,j-1,k+1,1)+bf(l,i  ,j  ,k+1,1))
        bx4=0.25d0*(bf(l,i  ,j  ,k  ,1)+bf(l,i  ,j+1,k  ,1)+bf(l,i  ,j  ,k+1,1)+bf(l,i  ,j+1,k+1,1))

        by1=0.25d0*(bf(l,i-1,j  ,k-1,2)+bf(l,i  ,j  ,k-1,2)+bf(l,i-1,j  ,k  ,2)+bf(l,i  ,j  ,k  ,2))
        by2=0.25d0*(bf(l,i  ,j+1,k-1,2)+bf(l,i-1,j+1,k-1,2)+bf(l,i  ,j+1,k  ,2)+bf(l,i-1,j+1,k  ,2))
        by3=0.25d0*(bf(l,i-1,j  ,k  ,2)+bf(l,i  ,j  ,k  ,2)+bf(l,i-1,j  ,k+1,2)+bf(l,i  ,j  ,k+1,2))
        by4=0.25d0*(bf(l,i  ,j+1,k  ,2)+bf(l,i-1,j+1,k  ,2)+bf(l,i  ,j+1,k+1,2)+bf(l,i-1,j+1,k+1,2))
        
        bz1=0.25d0*(bf(l,i-1,j-1,k  ,3)+bf(l,i-1,j  ,k  ,3)+bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3))
        bz2=0.25d0*(bf(l,i-1,j  ,k  ,3)+bf(l,i-1,j+1,k  ,3)+bf(l,i  ,j  ,k  ,3)+bf(l,i  ,j+1,k  ,3))
        bz3=0.25d0*(bf(l,i  ,j-1,k+1,3)+bf(l,i  ,j  ,k+1,3)+bf(l,i-1,j-1,k+1,3)+bf(l,i-1,j  ,k+1,3))
        bz4=0.25d0*(bf(l,i  ,j  ,k+1,3)+bf(l,i  ,j+1,k+1,3)+bf(l,i-1,j  ,k+1,3)+bf(l,i-1,j+1,k+1,3))

        Bnorm=sqrt(bx1*bx1+by1*by1+bz1*bz1)
        if(Bnorm.gt.0.0)then
           bx1=bx1/Bnorm
           by1=by1/Bnorm
           bz1=bz1/Bnorm
        endif
        Bnorm=sqrt(bx2*bx2+by2*by2+bz2*bz2)
        if(Bnorm.gt.0.0)then
           bx2=bx2/Bnorm
           by2=by2/Bnorm
           bz2=bz2/Bnorm
        endif
        Bnorm=sqrt(bx3*bx3+by3*by3+bz3*bz3)
        if(Bnorm.gt.0.0)then
           bx3=bx3/Bnorm
           by3=by3/Bnorm
           bz3=bz3/Bnorm
        endif
        Bnorm=sqrt(bx4*bx4+by4*by4+bz4*bz4)
        if(Bnorm.gt.0.0)then
           bx4=bx4/Bnorm
           by4=by4/Bnorm
           bz4=bz4/Bnorm
        endif

        if(alfven_diff_coeff)then
           kparax1=0.125d0*(Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  )+Dpara(l,i  ,j  ,k-1) &
                +      Dpara(l,i  ,j-1,k-1)+Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ) &
                +      Dpara(l,i-1,j  ,k-1)+Dpara(l,i-1,j-1,k-1))
           kparax2=0.125d0*(Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j+1,k-1) &
                +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i-1,j+1,k  )+Dpara(l,i-1,j  ,k  ) &
                +      Dpara(l,i-1,j+1,k-1)+Dpara(l,i-1,j  ,k-1))
           kparax3=0.125d0*(Dpara(l,i  ,j  ,k+1)+Dpara(l,i  ,j-1,k+1)+Dpara(l,i  ,j  ,k  ) &
                +      Dpara(l,i  ,j-1,k  )+Dpara(l,i-1,j  ,k+1)+Dpara(l,i-1,j-1,k+1) &
                +      Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ))
           kparax4=0.125d0*(Dpara(l,i  ,j+1,k+1)+Dpara(l,i  ,j  ,k+1)+Dpara(l,i  ,j+1,k  ) &
                +      Dpara(l,i  ,j  ,k  )+Dpara(l,i-1,j+1,k+1)+Dpara(l,i-1,j  ,k+1) &
                +      Dpara(l,i-1,j+1,k  )+Dpara(l,i-1,j  ,k  ))
           
           kperpx1=0.125d0*(kperp(l,i  ,j  ,k  )+kperp(l,i  ,j-1,k  )+kperp(l,i  ,j  ,k-1) &
                +      kperp(l,i  ,j-1,k-1)+kperp(l,i-1,j  ,k  )+kperp(l,i-1,j-1,k  ) &
                +      kperp(l,i-1,j  ,k-1)+kperp(l,i-1,j-1,k-1))
           kperpx2=0.125d0*(kperp(l,i  ,j+1,k  )+kperp(l,i  ,j  ,k  )+kperp(l,i  ,j+1,k-1) &
                +      kperp(l,i  ,j  ,k-1)+kperp(l,i-1,j+1,k  )+kperp(l,i-1,j  ,k  ) &
                +      kperp(l,i-1,j+1,k-1)+kperp(l,i-1,j  ,k-1))
           kperpx3=0.125d0*(kperp(l,i  ,j  ,k+1)+kperp(l,i  ,j-1,k+1)+kperp(l,i  ,j  ,k  ) &
                +      kperp(l,i  ,j-1,k  )+kperp(l,i-1,j  ,k+1)+kperp(l,i-1,j-1,k+1) &
                +      kperp(l,i-1,j  ,k  )+kperp(l,i-1,j-1,k  ))
           kperpx4=0.125d0*(kperp(l,i  ,j+1,k+1)+kperp(l,i  ,j  ,k+1)+kperp(l,i  ,j+1,k  ) &
                +      kperp(l,i  ,j  ,k  )+kperp(l,i-1,j+1,k+1)+kperp(l,i-1,j  ,k+1) &
                +      kperp(l,i-1,j+1,k  )+kperp(l,i-1,j  ,k  ))
           
           oneminuskperpx1 = 1.0d0-kperpx1
           oneminuskperpx2 = 1.0d0-kperpx2
           oneminuskperpx3 = 1.0d0-kperpx3
           oneminuskperpx4 = 1.0d0-kperpx4
        else if(streaming_diffusion)then
           kparax1=0.125d0*(Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  )+Dpara(l,i  ,j  ,k-1) &
                +      Dpara(l,i  ,j-1,k-1)+Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ) &
                +      Dpara(l,i-1,j  ,k-1)+Dpara(l,i-1,j-1,k-1)) + kpar
           kparax2=0.125d0*(Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j+1,k-1) &
                +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i-1,j+1,k  )+Dpara(l,i-1,j  ,k  ) &
                +      Dpara(l,i-1,j+1,k-1)+Dpara(l,i-1,j  ,k-1)) + kpar
           kparax3=0.125d0*(Dpara(l,i  ,j  ,k+1)+Dpara(l,i  ,j-1,k+1)+Dpara(l,i  ,j  ,k  ) &
                +      Dpara(l,i  ,j-1,k  )+Dpara(l,i-1,j  ,k+1)+Dpara(l,i-1,j-1,k+1) &
                +      Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  )) + kpar
           kparax4=0.125d0*(Dpara(l,i  ,j+1,k+1)+Dpara(l,i  ,j  ,k+1)+Dpara(l,i  ,j+1,k  ) &
                +      Dpara(l,i  ,j  ,k  )+Dpara(l,i-1,j+1,k+1)+Dpara(l,i-1,j  ,k+1) &
                +      Dpara(l,i-1,j+1,k  )+Dpara(l,i-1,j  ,k  )) + kpar
           kperpx1 = k_perp
           kperpx2 = k_perp
           kperpx3 = k_perp
           kperpx4 = k_perp
           oneminuskperpx1 = oneminuskperp
           oneminuskperpx2 = oneminuskperp
           oneminuskperpx3 = oneminuskperp
           oneminuskperpx4 = oneminuskperp           
        else
           kparax1 = kpar
           kparax2 = kpar
           kparax3 = kpar
           kparax4 = kpar
           kperpx1 = k_perp
           kperpx2 = k_perp
           kperpx3 = k_perp
           kperpx4 = k_perp
           oneminuskperpx1 = oneminuskperp
           oneminuskperpx2 = oneminuskperp
           oneminuskperpx3 = oneminuskperp
           oneminuskperpx4 = oneminuskperp
        endif

        if(compute .ne. 3)then   
           fx1=kparax1*(bx1*oneminuskperpx1*(bx1*dTdx1+by1*dTdy1+bz1*dTdz1)+kperpx1*dTdx1)
           fx2=kparax2*(bx2*oneminuskperpx2*(bx2*dTdx2+by2*dTdy2+bz2*dTdz2)+kperpx2*dTdx2)
           fx3=kparax3*(bx3*oneminuskperpx3*(bx3*dTdx3+by3*dTdy3+bz3*dTdz3)+kperpx3*dTdx3)
           fx4=kparax4*(bx4*oneminuskperpx4*(bx4*dTdx4+by4*dTdy4+bz4*dTdz4)+kperpx4*dTdx4)
        else ! Preconditionner
           fx1=kparax1*(bx1*oneminuskperpx1*(bx1+by1+bz1)+kperpx1)
           fx2=kparax2*(bx2*oneminuskperpx2*(bx2+by2+bz2)+kperpx2)
           fx3=kparax3*(bx3*oneminuskperpx3*(bx3+by3+bz3)+kperpx3)
           fx4=kparax4*(bx4*oneminuskperpx4*(bx4+by4+bz4)+kperpx4)
        end if
        fx=0.25d0*(fx1+fx2+fx3+fx4)
#endif
        myflux(l,i,j,k)=fx*dt/dx
     enddo
  enddo
  enddo
  enddo

endif

end subroutine cmpXcrflx
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmpYcrflx(Temp,Temp0,bf,myflux,dx,dt,ngrid,compute,ffdx,kpar,Dpara,kperp)
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
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Temp,Temp0,Dpara,kperp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ffdx

  real(dp)::Bnorm,fy,oneoverdx,oneovertwodx,oneoverfourdx,dx_loc
  real(dp)::dTdx1,dTdx2,dTdx3,dTdx4
  real(dp)::dTdy1,dTdy2,dTdy3,dTdy4
  real(dp)::dTdz1,dTdz2,dTdz3,dTdz4
  real(dp)::bx1,bx2,bx3,bx4
  real(dp)::by1,by2,by3,by4
  real(dp)::bz1,bz2,bz3,bz4
  real(dp)::fy1,fy2,fy3,fy4
  real(dp)::kpar,kpar2,oneminuskperp,kparay1,kparay2,kparay3,kparay4
  real(dp)::kperpy1,kperpy2,kperpy3,kperpy4
  real(dp)::oneminuskperpy1,oneminuskperpy2,oneminuskperpy3,oneminuskperpy4

  ! Local scalar variables
  integer::i,j,k,l
  integer::ilo,ihi,klo,khi

  real(dp)::a0,b0,c0,d0,a,b,c,d
  real(dp)::dTdx0,dTdxMinus0,dTdxPlus0,dTdz0,dTdzMinus0,dTdzPlus0
  real(dp)::dTdx,dTdxMinus,dTdxPlus,dTdz,dTdzMinus,dTdzPlus

  oneoverdx    =1.00d0/dx
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
        kpar2=0.0d0
        if(streaming_diffusion)kpar2 = 0.5d0*(Dpara(l,i,j,k)+Dpara(l,i,j-1,k))
        fy    =(kpar+kpar2)*dTdy1
        dx_loc=max(ffdx(l,i,j,k),ffdx(l,i,j-1,k))
        fy=fy/dx_loc
        myflux(l,i,j,k)=fy*dt/dx
        if(compute==3)myflux(l,i,j,k)=(kpar+kpar2)*dt/dx**2
     enddo
  enddo
  enddo
  enddo

else

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
           if(do_limiter_aniso)then
              ! Limit the normal component
              a0=Temp0(l,i  ,j,k)-Temp0(l,i  ,j-1,k)
              b0=Temp0(l,i-1,j,k)-Temp0(l,i-1,j-1,k)
              a =Temp (l,i  ,j,k)-Temp (l,i  ,j-1,k)
              b =Temp (l,i-1,j,k)-Temp (l,i-1,j-1,k)
              call L2_lim(a0,b0,a,b,alpha_limiter,dTdy1)
              b0=Temp0(l,i+1,j,k)-Temp0(l,i+1,j-1,k)
              b =Temp (l,i+1,j,k)-Temp (l,i+1,j-1,k)
              call L2_lim(a0,b0,a,b,alpha_limiter,dTdy2)
              dTdy1=dTdy1*oneoverdx
              dTdy2=dTdy2*oneoverdx

!!$              ! WARNING : this is a test to force asymetric
!!$              dTdy1=a*oneoverdx
!!$              dTdy2=a*oneoverdx

              
              ! Limit the transverse x-component
!!$              a0=Temp0(l,i  ,j-1,k)-Temp0(l,i-1,j-1,k)
!!$              b0=Temp0(l,i+1,j-1,k)-Temp0(l,i  ,j-1,k)
!!$              c0=Temp0(l,i  ,j,k  )-Temp0(l,i-1,j,k  )
!!$              d0=Temp0(l,i+1,j,k  )-Temp0(l,i  ,j,k  )
!!$              a =Temp (l,i  ,j-1,k)-Temp (l,i-1,j-1,k)
!!$              b =Temp (l,i+1,j-1,k)-Temp (l,i  ,j-1,k)
!!$              c =Temp (l,i  ,j,k  )-Temp (l,i-1,j,k  )
!!$              d =Temp (l,i+1,j,k  )-Temp (l,i  ,j,k  )
              if(compute==0)then
                 a0=Temp (l,i  ,j  ,k)-Temp (l,i-1,j  ,k)
                 b0=Temp (l,i+1,j  ,k)-Temp (l,i  ,j  ,k)
                 c0=Temp (l,i  ,j-1,k)-Temp (l,i-1,j-1,k)
                 d0=Temp (l,i+1,j-1,k)-Temp (l,i  ,j-1,k)
              else
                 a0=Temp0(l,i  ,j  ,k)-Temp0(l,i-1,j  ,k)
                 b0=Temp0(l,i+1,j  ,k)-Temp0(l,i  ,j  ,k)
                 c0=Temp0(l,i  ,j-1,k)-Temp0(l,i-1,j-1,k)
                 d0=Temp0(l,i+1,j-1,k)-Temp0(l,i  ,j-1,k)
              endif
              a =Temp (l,i  ,j  ,k)-Temp (l,i-1,j  ,k)
              b =Temp (l,i+1,j  ,k)-Temp (l,i  ,j  ,k)
              c =Temp (l,i  ,j-1,k)-Temp (l,i-1,j-1,k)
              d =Temp (l,i+1,j-1,k)-Temp (l,i  ,j-1,k)
              if(slope_limiter_aniso == 1)then
                 call MinMod_lim(a0,b0,a,b,dTdxPlus0 ,dTdxPlus )
                 call MinMod_lim(c0,d0,c,d,dTdxMinus0,dTdxMinus)
                 call MinMod_lim(dTdxPlus0,dTdxMinus0,dTdxPlus,dTdxMinus,dTdx0,dTdx)
              !elseif(slope_limiter_aniso == 2)then
              !   call MonCen_lim(a0,b0,a,b,dTdxPlus0 ,dTdxPlus )
              !   call MonCen_lim(c0,d0,c,d,dTdxMinus0,dTdxMinus)
              !   call MonCen_lim(dTdxPlus0,dTdxMinus0,dTdxPlus,dTdxMinus,dTdx0,dTdx)
              endif
              dTdx1=dTdx*oneoverdx
              dTdx2=dTdx*oneoverdx
           else
              dTdx1=(Temp(l,i  ,j  ,k)+Temp(l,i  ,j-1,k) &
                   - Temp(l,i-1,j  ,k)-Temp(l,i-1,j-1,k))*oneovertwodx
              dTdx2=(Temp(l,i+1,j  ,k)+Temp(l,i+1,j-1,k) &
                   - Temp(l,i  ,j  ,k)-Temp(l,i  ,j-1,k))*oneovertwodx        
           
              dTdy1=(Temp(l,i  ,j  ,k)+Temp(l,i-1,j  ,k) &
                   - Temp(l,i  ,j-1,k)-Temp(l,i-1,j-1,k))*oneovertwodx
              dTdy2=(Temp(l,i  ,j  ,k)+Temp(l,i+1,j  ,k) &
                   - Temp(l,i  ,j-1,k)-Temp(l,i+1,j-1,k))*oneovertwodx
           end if
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

!!$        kpar=4d0/(Dpara(l,i  ,j  ,k)+Dpara(l,i  ,j-1,k) &
!!$             +      Dpara(l,i-1,j  ,k)+Dpara(l,i-1,j-1,k))
!!$        kpar=4d0/(Dpara(l,i  ,j  ,k)+Dpara(l,i+1,j  ,k) &
!!$             +      Dpara(l,i  ,j-1,k)+Dpara(l,i+1,j-1,k))
        if(alfven_diff_coeff)then
           kparay1=0.25d0*(Dpara(l,i  ,j  ,k)+Dpara(l,i  ,j-1,k) &
                +      Dpara(l,i-1,j  ,k)+Dpara(l,i-1,j-1,k))
           kparay2=0.25d0*(Dpara(l,i  ,j  ,k)+Dpara(l,i+1,j  ,k) &
                +      Dpara(l,i  ,j-1,k)+Dpara(l,i+1,j-1,k))
           
           kperpy1=0.25d0*(kperp(l,i  ,j  ,k)+kperp(l,i  ,j-1,k) &
                +      kperp(l,i-1,j  ,k)+kperp(l,i-1,j-1,k))
           kperpy2=0.25d0*(kperp(l,i  ,j  ,k)+kperp(l,i+1,j  ,k) &
                +      kperp(l,i  ,j-1,k)+kperp(l,i+1,j-1,k))
           
           oneminuskperpy1 = 1.0d0-kperpy1
           oneminuskperpy2 = 1.0d0-kperpy2
        else if(streaming_diffusion)then
           kparay1=0.25d0*(Dpara(l,i  ,j  ,k)+Dpara(l,i  ,j-1,k) &
                +      Dpara(l,i-1,j  ,k)+Dpara(l,i-1,j-1,k)) + kpar
           kparay2=0.25d0*(Dpara(l,i  ,j  ,k)+Dpara(l,i+1,j  ,k) &
                +      Dpara(l,i  ,j-1,k)+Dpara(l,i+1,j-1,k)) + kpar
           kperpy1 = k_perp
           kperpy2 = k_perp
           oneminuskperpy1 = oneminuskperp
           oneminuskperpy2 = oneminuskperp
        else
           kparay1 = kpar
           kparay2 = kpar
           kperpy1 = k_perp
           kperpy2 = k_perp
           oneminuskperpy1 = oneminuskperp
           oneminuskperpy2 = oneminuskperp
        endif

        if(compute .ne. 3)then   
           fy1=kparay1*(by1*oneminuskperpy1*(bx1*dTdx1+by1*dTdy1)+kperpy1*dTdy1)
           fy2=kparay2*(by2*oneminuskperpy2*(bx2*dTdx2+by2*dTdy2)+kperpy2*dTdy2)
        else ! Preconditionner
           fy1=kparay1*(by1*oneminuskperpy1*(bx1+by1)+kperpy1)
           fy2=kparay2*(by2*oneminuskperpy2*(bx2+by2)+kperpy2)
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
           if(do_limiter_aniso)then
              ! Limit the normal y-component
              ! T:Top , M:Middle, B:Bottom in k
              ! L:Left, C:Center, R:Right  in i
!!$              jj=j-1
!!$              ii=i-1;kk=k+1;TL=Temp(l,ii,j,kk)-Temp(l,ii,jj,kk)
!!$              ii=i         ;TC=Temp(l,ii,j,kk)-Temp(l,ii,jj,kk)
!!$              ii=i+1       ;TR=Temp(l,ii,j,kk)-Temp(l,ii,jj,kk)
!!$              ii=i-1;kk=k  ;ML=Temp(l,ii,j,kk)-Temp(l,ii,jj,kk)
!!$              ii=i         ;MC=Temp(l,ii,j,kk)-Temp(l,ii,jj,kk)
!!$              ii=i+1       ;MR=Temp(l,ii,j,kk)-Temp(l,ii,jj,kk)
!!$              ii=i-1;kk=k-1;BL=Temp(l,ii,j,kk)-Temp(l,ii,jj,kk)
!!$              ii=i         ;BC=Temp(l,ii,j,kk)-Temp(l,ii,jj,kk)
!!$              ii=i+1       ;BR=Temp(l,ii,j,kk)-Temp(l,ii,jj,kk)
              dTdy1=(Temp(l,i,j,k)-Temp(l,i,j-1,k))*oneoverdx
              dTdy2=dTdy1
              dTdy3=dTdy1
              dTdy4=dTdy1

              ! Limit the transverse x-component
              if(compute==0)then
                 a0=Temp (l,i  ,j  ,k)-Temp (l,i-1,j  ,k)
                 b0=Temp (l,i+1,j  ,k)-Temp (l,i  ,j  ,k)
                 c0=Temp (l,i  ,j-1,k)-Temp (l,i-1,j-1,k)
                 d0=Temp (l,i+1,j-1,k)-Temp (l,i  ,j-1,k)
              else
                 a0=Temp0(l,i  ,j  ,k)-Temp0(l,i-1,j  ,k)
                 b0=Temp0(l,i+1,j  ,k)-Temp0(l,i  ,j  ,k)
                 c0=Temp0(l,i  ,j-1,k)-Temp0(l,i-1,j-1,k)
                 d0=Temp0(l,i+1,j-1,k)-Temp0(l,i  ,j-1,k)
              endif
              a =Temp (l,i  ,j  ,k)-Temp (l,i-1,j  ,k)
              b =Temp (l,i+1,j  ,k)-Temp (l,i  ,j  ,k)
              c =Temp (l,i  ,j-1,k)-Temp (l,i-1,j-1,k)
              d =Temp (l,i+1,j-1,k)-Temp (l,i  ,j-1,k)
              if(slope_limiter_aniso == 1)then
                 call MinMod_lim(a0,b0,a,b,dTdxPlus0 ,dTdxPlus )
                 call MinMod_lim(c0,d0,c,d,dTdxMinus0,dTdxMinus)
                 call MinMod_lim(dTdxPlus0,dTdxMinus0,dTdxPlus,dTdxMinus,dTdx0,dTdx)
              !elseif(slope_limiter_aniso == 2)then
              !   call MonCen_lim(a,b,dTdxMinus)
              !   call MonCen_lim(c,d,dTdxPlus )
              !   call MonCen_lim(dTdxMinus,dTdxPlus,dTdx)
              endif
              dTdx1=dTdx*oneoverdx
              dTdx2=dTdx1
              dTdx3=dTdx1
              dTdx4=dTdx1

              ! Limit the transverse z-component
              if(compute==0)then
                 a0=Temp (l,i,j  ,k  )-Temp (l,i,j  ,k-1)
                 b0=Temp (l,i,j  ,k+1)-Temp (l,i,j  ,k  )
                 c0=Temp (l,i,j-1,k  )-Temp (l,i,j-1,k-1)
                 d0=Temp (l,i,j-1,k+1)-Temp (l,i,j-1,k  )
              else
                 a0=Temp0(l,i,j  ,k  )-Temp0(l,i,j  ,k-1)
                 b0=Temp0(l,i,j  ,k+1)-Temp0(l,i,j  ,k  )
                 c0=Temp0(l,i,j-1,k  )-Temp0(l,i,j-1,k-1)
                 d0=Temp0(l,i,j-1,k+1)-Temp0(l,i,j-1,k  )
              endif
              a =Temp (l,i,j  ,k  )-Temp (l,i,j  ,k-1)
              b =Temp (l,i,j  ,k+1)-Temp (l,i,j  ,k  )
              c =Temp (l,i,j-1,k  )-Temp (l,i,j-1,k-1)
              d =Temp (l,i,j-1,k+1)-Temp (l,i,j-1,k  )
              if(slope_limiter_aniso == 1)then
                 call MinMod_lim(a0,b0,a,b,dTdzPlus0 ,dTdzPlus )
                 call MinMod_lim(c0,d0,c,d,dTdzMinus0,dTdzMinus)
                 call MinMod_lim(dTdzPlus0,dTdzMinus0,dTdzPlus,dTdzMinus,dTdz0,dTdz)
              !elseif(slope_limiter_aniso == 2)then
              !   call MonCen_lim(a,b,dTdzMinus)
              !   call MonCen_lim(c,d,dTdzPlus )
              !   call MonCen_lim(dTdzMinus,dTdzPlus,dTdz)
              endif
              dTdz1=dTdz*oneoverdx
              dTdz2=dTdz1
              dTdz3=dTdz1
              dTdz4=dTdz1

           else
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
           endif
        end if

        bx1=0.25d0*(bf(l,i  ,j-1,k-1,1)+bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j-1,k  ,1)+bf(l,i  ,j  ,k  ,1))
        bx2=0.25d0*(bf(l,i+1,j-1,k-1,1)+bf(l,i+1,j  ,k-1,1)+bf(l,i+1,j-1,k  ,1)+bf(l,i+1,j  ,k  ,1))
        bx3=0.25d0*(bf(l,i  ,j-1,k  ,1)+bf(l,i  ,j  ,k  ,1)+bf(l,i  ,j-1,k+1,1)+bf(l,i  ,j  ,k+1,1))
        bx4=0.25d0*(bf(l,i+1,j-1,k  ,1)+bf(l,i+1,j  ,k  ,1)+bf(l,i+1,j-1,k+1,1)+bf(l,i+1,j  ,k+1,1))

        by1=0.25d0*(bf(l,i-1,j  ,k-1,2)+bf(l,i  ,j  ,k-1,2)+bf(l,i-1,j  ,k  ,2)+bf(l,i  ,j  ,k  ,2))
        by2=0.25d0*(bf(l,i  ,j  ,k-1,2)+bf(l,i+1,j  ,k-1,2)+bf(l,i  ,j  ,k  ,2)+bf(l,i+1,j  ,k  ,2))
        by3=0.25d0*(bf(l,i-1,j  ,k  ,2)+bf(l,i  ,j  ,k  ,2)+bf(l,i-1,j  ,k+1,2)+bf(l,i  ,j  ,k+1,2))
        by4=0.25d0*(bf(l,i  ,j  ,k  ,2)+bf(l,i+1,j  ,k  ,2)+bf(l,i  ,j  ,k+1,2)+bf(l,i+1,j  ,k+1,2))
        
        bz1=0.25d0*(bf(l,i-1,j-1,k  ,3)+bf(l,i-1,j  ,k  ,3)+bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3))
        bz2=0.25d0*(bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3)+bf(l,i+1,j-1,k  ,3)+bf(l,i+1,j  ,k  ,3))
        bz3=0.25d0*(bf(l,i  ,j-1,k+1,3)+bf(l,i  ,j  ,k+1,3)+bf(l,i-1,j-1,k+1,3)+bf(l,i-1,j  ,k+1,3))
        bz4=0.25d0*(bf(l,i+1,j-1,k+1,3)+bf(l,i+1,j  ,k+1,3)+bf(l,i  ,j-1,k+1,3)+bf(l,i  ,j  ,k+1,3))

        Bnorm=sqrt(bx1*bx1+by1*by1+bz1*bz1)
        if(Bnorm.gt.0.0)then
           bx1=bx1/Bnorm
           by1=by1/Bnorm
           bz1=bz1/Bnorm
        endif
        Bnorm=sqrt(bx2*bx2+by2*by2+bz2*bz2)
        if(Bnorm.gt.0.0)then
           bx2=bx2/Bnorm
           by2=by2/Bnorm
           bz2=bz2/Bnorm
        endif
        Bnorm=sqrt(bx3*bx3+by3*by3+bz3*bz3)
        if(Bnorm.gt.0.0)then
           bx3=bx3/Bnorm
           by3=by3/Bnorm
           bz3=bz3/Bnorm
        endif
        Bnorm=sqrt(bx4*bx4+by4*by4+bz4*bz4)
        if(Bnorm.gt.0.0)then
           bx4=bx4/Bnorm
           by4=by4/Bnorm
           bz4=bz4/Bnorm
        endif

!!$        kpar=8d0/(Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  )+Dpara(l,i  ,j  ,k-1) &
!!$             +      Dpara(l,i  ,j-1,k-1)+Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ) &
!!$             +      Dpara(l,i-1,j  ,k-1)+Dpara(l,i-1,j-1,k-1))
!!$        kpar=8d0/(Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j-1,k  )+Dpara(l,i+1,j  ,k-1) &
!!$             +      Dpara(l,i+1,j-1,k-1)+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  ) &
!!$             +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i  ,j-1,k-1))
!!$        kpar=8d0/(Dpara(l,i  ,j  ,k+1)+Dpara(l,i  ,j-1,k+1)+Dpara(l,i  ,j  ,k  ) &
!!$             +      Dpara(l,i  ,j-1,k  )+Dpara(l,i-1,j  ,k+1)+Dpara(l,i-1,j-1,k+1) &
!!$             +      Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ))
!!$        kpar=8d0/(Dpara(l,i+1,j  ,k+1)+Dpara(l,i+1,j-1,k+1)+Dpara(l,i+1,j  ,k  ) &
!!$             +      Dpara(l,i+1,j-1,k  )+Dpara(l,i  ,j  ,k+1)+Dpara(l,i  ,j-1,k+1) &
!!$             +      Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  ))
        if(alfven_diff_coeff)then
           kparay1=0.125d0*(Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  )+Dpara(l,i  ,j  ,k-1) &
                +      Dpara(l,i  ,j-1,k-1)+Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ) &
                +      Dpara(l,i-1,j  ,k-1)+Dpara(l,i-1,j-1,k-1))
           kparay2=0.125d0*(Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j-1,k  )+Dpara(l,i+1,j  ,k-1) &
                +      Dpara(l,i+1,j-1,k-1)+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  ) &
                +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i  ,j-1,k-1))
           kparay3=0.125d0*(Dpara(l,i  ,j  ,k+1)+Dpara(l,i  ,j-1,k+1)+Dpara(l,i  ,j  ,k  ) &
                +      Dpara(l,i  ,j-1,k  )+Dpara(l,i-1,j  ,k+1)+Dpara(l,i-1,j-1,k+1) &
                +      Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ))
           kparay4=0.125d0*(Dpara(l,i+1,j  ,k+1)+Dpara(l,i+1,j-1,k+1)+Dpara(l,i+1,j  ,k  ) &
                +      Dpara(l,i+1,j-1,k  )+Dpara(l,i  ,j  ,k+1)+Dpara(l,i  ,j-1,k+1) &
                +      Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  ))

           kperpy1=0.125d0*(kperp(l,i  ,j  ,k  )+kperp(l,i  ,j-1,k  )+kperp(l,i  ,j  ,k-1) &
             +      kperp(l,i  ,j-1,k-1)+kperp(l,i-1,j  ,k  )+kperp(l,i-1,j-1,k  ) &
             +      kperp(l,i-1,j  ,k-1)+kperp(l,i-1,j-1,k-1))
           kperpy2=0.125d0*(kperp(l,i+1,j  ,k  )+kperp(l,i+1,j-1,k  )+kperp(l,i+1,j  ,k-1) &
                +      kperp(l,i+1,j-1,k-1)+kperp(l,i  ,j  ,k  )+kperp(l,i  ,j-1,k  ) &
                +      kperp(l,i  ,j  ,k-1)+kperp(l,i  ,j-1,k-1))
           kperpy3=0.125d0*(kperp(l,i  ,j  ,k+1)+kperp(l,i  ,j-1,k+1)+kperp(l,i  ,j  ,k  ) &
                +      kperp(l,i  ,j-1,k  )+kperp(l,i-1,j  ,k+1)+kperp(l,i-1,j-1,k+1) &
                +      kperp(l,i-1,j  ,k  )+kperp(l,i-1,j-1,k  ))
           kperpy4=0.125d0*(kperp(l,i+1,j  ,k+1)+kperp(l,i+1,j-1,k+1)+kperp(l,i+1,j  ,k  ) &
                +      kperp(l,i+1,j-1,k  )+kperp(l,i  ,j  ,k+1)+kperp(l,i  ,j-1,k+1) &
                +      kperp(l,i  ,j  ,k  )+kperp(l,i  ,j-1,k  ))

           oneminuskperpy1 = 1.0d0-kperpy1
           oneminuskperpy2 = 1.0d0-kperpy2
           oneminuskperpy3 = 1.0d0-kperpy3
           oneminuskperpy4 = 1.0d0-kperpy4

        else if(streaming_diffusion)then
           kparay1=0.125d0*(Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  )+Dpara(l,i  ,j  ,k-1) &
                +      Dpara(l,i  ,j-1,k-1)+Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ) &
                +      Dpara(l,i-1,j  ,k-1)+Dpara(l,i-1,j-1,k-1)) + kpar
           kparay2=0.125d0*(Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j-1,k  )+Dpara(l,i+1,j  ,k-1) &
                +      Dpara(l,i+1,j-1,k-1)+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  ) &
                +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i  ,j-1,k-1)) + kpar
           kparay3=0.125d0*(Dpara(l,i  ,j  ,k+1)+Dpara(l,i  ,j-1,k+1)+Dpara(l,i  ,j  ,k  ) &
                +      Dpara(l,i  ,j-1,k  )+Dpara(l,i-1,j  ,k+1)+Dpara(l,i-1,j-1,k+1) &
                +      Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  )) + kpar
           kparay4=0.125d0*(Dpara(l,i+1,j  ,k+1)+Dpara(l,i+1,j-1,k+1)+Dpara(l,i+1,j  ,k  ) &
                +      Dpara(l,i+1,j-1,k  )+Dpara(l,i  ,j  ,k+1)+Dpara(l,i  ,j-1,k+1) &
                +      Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  )) + kpar
           kperpy1 = k_perp
           kperpy2 = k_perp
           kperpy3 = k_perp
           kperpy4 = k_perp
           oneminuskperpy1 = oneminuskperp
           oneminuskperpy2 = oneminuskperp
           oneminuskperpy3 = oneminuskperp
           oneminuskperpy4 = oneminuskperp
        else
           kparay1 = kpar
           kparay2 = kpar
           kparay3 = kpar
           kparay4 = kpar
           kperpy1 = k_perp
           kperpy2 = k_perp
           kperpy3 = k_perp
           kperpy4 = k_perp
           oneminuskperpy1 = oneminuskperp
           oneminuskperpy2 = oneminuskperp
           oneminuskperpy3 = oneminuskperp
           oneminuskperpy4 = oneminuskperp
        end if
        
        if(compute .ne. 3)then        
           fy1=kparay1*(by1*oneminuskperpy1*(bx1*dTdx1+by1*dTdy1+bz1*dTdz1)+kperpy1*dTdy1)
           fy2=kparay2*(by2*oneminuskperpy2*(bx2*dTdx2+by2*dTdy2+bz2*dTdz2)+kperpy2*dTdy2)
           fy3=kparay3*(by3*oneminuskperpy3*(bx3*dTdx3+by3*dTdy3+bz3*dTdz3)+kperpy3*dTdy3)
           fy4=kparay4*(by4*oneminuskperpy4*(bx4*dTdx4+by4*dTdy4+bz4*dTdz4)+kperpy4*dTdy4)
        else
           fy1=kparay1*(by1*oneminuskperpy1*(bx1+by1+bz1)+kperpy1)
           fy2=kparay2*(by2*oneminuskperpy2*(bx2+by2+bz2)+kperpy2)
           fy3=kparay3*(by3*oneminuskperpy3*(bx3+by3+bz3)+kperpy3)
           fy4=kparay4*(by4*oneminuskperpy4*(bx4+by4+bz4)+kperpy4)
        end if
        fy=0.25d0*(fy1+fy2+fy3+fy4)
#endif

        myflux(l,i,j,k)=fy*dt/dx
     enddo
  enddo
  enddo
  enddo
endif

end subroutine cmpYcrflx
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cmpZcrflx(Temp,Temp0,bf,myflux,dx,dt,ngrid,compute,ffdx,kpar,Dpara,kperp)
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
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Temp,Temp0,Dpara,kperp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ffdx

  real(dp)::Bnorm,fz,oneoverdx,oneovertwodx,oneoverfourdx,dx_loc
  real(dp)::dTdx1,dTdx2,dTdx3,dTdx4
  real(dp)::dTdy1,dTdy2,dTdy3,dTdy4
  real(dp)::dTdz1,dTdz2,dTdz3,dTdz4
  real(dp)::bx1,bx2,bx3,bx4
  real(dp)::by1,by2,by3,by4
  real(dp)::bz1,bz2,bz3,bz4
  real(dp)::fz1,fz2,fz3,fz4
  real(dp)::kpar,kpar2,oneminuskperp,kparaz1,kparaz2,kparaz3,kparaz4
  real(dp)::kperpz1,kperpz2,kperpz3,kperpz4
  real(dp)::oneminuskperpz1,oneminuskperpz2,oneminuskperpz3,oneminuskperpz4

  real(dp)::a0,b0,c0,d0,a,b,c,d
  real(dp)::dTdx0,dTdxMinus0,dTdxPlus0,dTdy0,dTdyMinus0,dTdyPlus0
  real(dp)::dTdx,dTdxMinus,dTdxPlus,dTdy,dTdyMinus,dTdyPlus

  ! Local scalar variables
  integer::i,j,k,l
  integer::ilo,ihi,jlo,jhi

  oneoverdx    =1.00d0/dx
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
        kpar2=0.0d0
        if(streaming_diffusion)kpar2 = 0.5d0*(Dpara(l,i,j,k)+Dpara(l,i,j,k-1))
        fz    =(kpar+kpar2)*dTdz1
        dx_loc=max(ffdx(l,i,j,k),ffdx(l,i,j,k-1))
        fz=fz/dx_loc
        myflux(l,i,j,k)=fz*dt/dx
        if(compute==3)myflux(l,i,j,k)=(kpar+kpar2)*dt/dx**2
     enddo
  enddo
  enddo
  enddo

else

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
           if(do_limiter_aniso)then
              ! Limit the normal z-component
              ! T:Top , M:Middle, B:Bottom in j
              ! L:Left, C:Center, R:Right  in i
!!$              kk=k-1
!!$              ii=i-1;jj=j+1;TL=Temp(l,ii,jj,k)-Temp(l,ii,jj,kk)
!!$              ii=i         ;TC=Temp(l,ii,jj,k)-Temp(l,ii,jj,kk)
!!$              ii=i+1       ;TR=Temp(l,ii,jj,k)-Temp(l,ii,jj,kk)
!!$              ii=i-1;jj=j  ;ML=Temp(l,ii,jj,k)-Temp(l,ii,jj,kk)
!!$              ii=i         ;MC=Temp(l,ii,jj,k)-Temp(l,ii,jj,kk)
!!$              ii=i+1       ;MR=Temp(l,ii,jj,k)-Temp(l,ii,jj,kk)
!!$              ii=i-1;jj=j-1;BL=Temp(l,ii,jj,k)-Temp(l,ii,jj,kk)
!!$              ii=i         ;BC=Temp(l,ii,jj,k)-Temp(l,ii,jj,kk)
!!$              ii=i+1       ;BR=Temp(l,ii,jj,k)-Temp(l,ii,jj,kk)
              dTdz1=(Temp(l,i,j,k)-Temp(l,i,j,k-1))*oneoverdx
              dTdz2=dTdz1
              dTdz3=dTdz1
              dTdz4=dTdz1

              ! Limit the transverse x-component
              if(compute==0)then
                 a0=Temp (l,i  ,j,k  )-Temp (l,i-1,j,k  )
                 b0=Temp (l,i+1,j,k  )-Temp (l,i  ,j,k  )
                 c0=Temp (l,i  ,j,k-1)-Temp (l,i-1,j,k-1)
                 d0=Temp (l,i+1,j,k-1)-Temp (l,i  ,j,k-1)
              else
                 a0=Temp0(l,i  ,j,k  )-Temp0(l,i-1,j,k  )
                 b0=Temp0(l,i+1,j,k  )-Temp0(l,i  ,j,k  )
                 c0=Temp0(l,i  ,j,k-1)-Temp0(l,i-1,j,k-1)
                 d0=Temp0(l,i+1,j,k-1)-Temp0(l,i  ,j,k-1)
              endif
              a =Temp (l,i  ,j,k  )-Temp (l,i-1,j,k  )
              b =Temp (l,i+1,j,k  )-Temp (l,i  ,j,k  )
              c =Temp (l,i  ,j,k-1)-Temp (l,i-1,j,k-1)
              d =Temp (l,i+1,j,k-1)-Temp (l,i  ,j,k-1)
              if(slope_limiter_aniso == 1)then
                 call MinMod_lim(a0,b0,a,b,dTdxPlus0 ,dTdxPlus )
                 call MinMod_lim(c0,d0,c,d,dTdxMinus0,dTdxMinus)
                 call MinMod_lim(dTdxPlus0,dTdxMinus0,dTdxPlus,dTdxMinus,dTdx0,dTdx)
              !elseif(slope_limiter_aniso == 2)then
              !   call MonCen_lim(a,b,dTdxMinus)
              !   call MonCen_lim(c,d,dTdxPlus )
              !   call MonCen_lim(dTdxMinus,dTdxPlus,dTdx)
              endif
              dTdx1=dTdx*oneoverdx
              dTdx2=dTdx1
              dTdx3=dTdx1
              dTdx4=dTdx1

              ! Limit the transverse y-component
              if(compute==0)then
                 a0=Temp (l,i,j  ,k  )-Temp (l,i,j-1,k  )
                 b0=Temp (l,i,j+1,k  )-Temp (l,i,j  ,k  )
                 c0=Temp (l,i,j  ,k-1)-Temp (l,i,j-1,k-1)
                 d0=Temp (l,i,j+1,k-1)-Temp (l,i,j  ,k-1)
              else
                 a0=Temp0(l,i,j  ,k  )-Temp0(l,i,j-1,k  )
                 b0=Temp0(l,i,j+1,k  )-Temp0(l,i,j  ,k  )
                 c0=Temp0(l,i,j  ,k-1)-Temp0(l,i,j-1,k-1)
                 d0=Temp0(l,i,j+1,k-1)-Temp0(l,i,j  ,k-1)
              endif
              a =Temp (l,i,j  ,k  )-Temp (l,i,j-1,k  )
              b =Temp (l,i,j+1,k  )-Temp (l,i,j  ,k  )
              c =Temp (l,i,j  ,k-1)-Temp (l,i,j-1,k-1)
              d =Temp (l,i,j+1,k-1)-Temp (l,i,j  ,k-1)
              if(slope_limiter_aniso == 1)then
                 call MinMod_lim(a0,b0,a,b,dTdyPlus0 ,dTdyPlus )
                 call MinMod_lim(c0,d0,c,d,dTdyMinus0,dTdyMinus)
                 call MinMod_lim(dTdyPlus0,dTdyMinus0,dTdyPlus,dTdyMinus,dTdy0,dTdy)
              !elseif(slope_limiter_aniso == 2)then
              !   call MonCen_lim(a,b,dTdyMinus)
              !   call MonCen_lim(c,d,dTdyPlus )
              !   call MonCen_lim(dTdyMinus,dTdyPlus,dTdy)
              endif
              dTdy1=dTdy*oneoverdx
              dTdy2=dTdy1
              dTdy3=dTdy1
              dTdy4=dTdy1

           else
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
        endif
        
        bx1=0.25d0*(bf(l,i  ,j-1,k-1,1)+bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j-1,k  ,1)+bf(l,i  ,j  ,k  ,1))
        bx2=0.25d0*(bf(l,i+1,j-1,k-1,1)+bf(l,i+1,j  ,k-1,1)+bf(l,i+1,j-1,k  ,1)+bf(l,i+1,j  ,k  ,1))
        bx3=0.25d0*(bf(l,i  ,j  ,k-1,1)+bf(l,i  ,j+1,k-1,1)+bf(l,i  ,j  ,k  ,1)+bf(l,i  ,j+1,k  ,1))
        bx4=0.25d0*(bf(l,i+1,j  ,k-1,1)+bf(l,i+1,j+1,k-1,1)+bf(l,i+1,j  ,k  ,1)+bf(l,i+1,j+1,k  ,1))

        by1=0.25d0*(bf(l,i-1,j  ,k-1,2)+bf(l,i  ,j  ,k-1,2)+bf(l,i-1,j  ,k  ,2)+bf(l,i  ,j  ,k  ,2))
        by2=0.25d0*(bf(l,i  ,j  ,k-1,2)+bf(l,i+1,j  ,k-1,2)+bf(l,i  ,j  ,k  ,2)+bf(l,i+1,j  ,k  ,2))
        by3=0.25d0*(bf(l,i-1,j+1,k-1,2)+bf(l,i  ,j+1,k-1,2)+bf(l,i-1,j+1,k  ,2)+bf(l,i  ,j+1,k  ,2))
        by4=0.25d0*(bf(l,i  ,j+1,k-1,2)+bf(l,i+1,j+1,k-1,2)+bf(l,i  ,j+1,k  ,2)+bf(l,i+1,j+1,k  ,2))
        
        bz1=0.25d0*(bf(l,i-1,j-1,k  ,3)+bf(l,i-1,j  ,k  ,3)+bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3))
        bz2=0.25d0*(bf(l,i  ,j-1,k  ,3)+bf(l,i  ,j  ,k  ,3)+bf(l,i+1,j-1,k  ,3)+bf(l,i+1,j  ,k  ,3))
        bz3=0.25d0*(bf(l,i-1,j  ,k  ,3)+bf(l,i-1,j+1,k  ,3)+bf(l,i  ,j  ,k  ,3)+bf(l,i  ,j+1,k  ,3))
        bz4=0.25d0*(bf(l,i  ,j  ,k  ,3)+bf(l,i  ,j+1,k  ,3)+bf(l,i+1,j  ,k  ,3)+bf(l,i+1,j+1,k  ,3))

        Bnorm=sqrt(bx1*bx1+by1*by1+bz1*bz1)
        if(Bnorm.gt.0.0)then
           bx1=bx1/Bnorm
           by1=by1/Bnorm
           bz1=bz1/Bnorm
        endif
        Bnorm=sqrt(bx2*bx2+by2*by2+bz2*bz2)
        if(Bnorm.gt.0.0)then
           bx2=bx2/Bnorm
           by2=by2/Bnorm
           bz2=bz2/Bnorm
        endif
        Bnorm=sqrt(bx3*bx3+by3*by3+bz3*bz3)
        if(Bnorm.gt.0.0)then
           bx3=bx3/Bnorm
           by3=by3/Bnorm
           bz3=bz3/Bnorm
        endif
        Bnorm=sqrt(bx4*bx4+by4*by4+bz4*bz4)
        if(Bnorm.gt.0.0)then
           bx4=bx4/Bnorm
           by4=by4/Bnorm
           bz4=bz4/Bnorm
        endif

!!$        kpar=8d0/(Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  )+Dpara(l,i  ,j  ,k-1) &
!!$             +      Dpara(l,i  ,j-1,k-1)+Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ) &
!!$             +      Dpara(l,i-1,j  ,k-1)+Dpara(l,i-1,j-1,k-1))
!!$        kpar=8d0/(Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j-1,k  )+Dpara(l,i+1,j  ,k-1) &
!!$             +      Dpara(l,i+1,j-1,k-1)+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  ) &
!!$             +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i  ,j-1,k-1))
!!$        kpar=8d0/(Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j+1,k-1) &
!!$             +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i-1,j+1,k  )+Dpara(l,i-1,j  ,k  ) &
!!$             +      Dpara(l,i-1,j+1,k-1)+Dpara(l,i-1,j  ,k-1))
!!$        kpar=8d0/(Dpara(l,i+1,j+1,k  )+Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j+1,k-1) &
!!$             +      Dpara(l,i+1,j  ,k-1)+Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  ) &
!!$             +      Dpara(l,i  ,j+1,k-1)+Dpara(l,i  ,j  ,k-1))
!!$
        if(alfven_diff_coeff)then
           kparaz1=0.125d0*(Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  )+Dpara(l,i  ,j  ,k-1) &
                +      Dpara(l,i  ,j-1,k-1)+Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ) &
                +      Dpara(l,i-1,j  ,k-1)+Dpara(l,i-1,j-1,k-1))
           kparaz2=0.125d0*(Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j-1,k  )+Dpara(l,i+1,j  ,k-1) &
                +      Dpara(l,i+1,j-1,k-1)+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  ) &
                +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i  ,j-1,k-1))
           kparaz3=0.125d0*(Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j+1,k-1) &
                +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i-1,j+1,k  )+Dpara(l,i-1,j  ,k  ) &
                +      Dpara(l,i-1,j+1,k-1)+Dpara(l,i-1,j  ,k-1))
           kparaz4=0.125d0*(Dpara(l,i+1,j+1,k  )+Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j+1,k-1) &
                +      Dpara(l,i+1,j  ,k-1)+Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  ) &
                +      Dpara(l,i  ,j+1,k-1)+Dpara(l,i  ,j  ,k-1))

           kperpz1=0.125d0*(kperp(l,i  ,j  ,k  )+kperp(l,i  ,j-1,k  )+kperp(l,i  ,j  ,k-1) &
                +      kperp(l,i  ,j-1,k-1)+kperp(l,i-1,j  ,k  )+kperp(l,i-1,j-1,k  ) &
                +      kperp(l,i-1,j  ,k-1)+kperp(l,i-1,j-1,k-1))
           kperpz2=0.125d0*(kperp(l,i+1,j  ,k  )+kperp(l,i+1,j-1,k  )+kperp(l,i+1,j  ,k-1) &
                +      kperp(l,i+1,j-1,k-1)+kperp(l,i  ,j  ,k  )+kperp(l,i  ,j-1,k  ) &
                +      kperp(l,i  ,j  ,k-1)+kperp(l,i  ,j-1,k-1))
           kperpz3=0.125d0*(kperp(l,i  ,j+1,k  )+kperp(l,i  ,j  ,k  )+kperp(l,i  ,j+1,k-1) &
                +      kperp(l,i  ,j  ,k-1)+kperp(l,i-1,j+1,k  )+kperp(l,i-1,j  ,k  ) &
                +      kperp(l,i-1,j+1,k-1)+kperp(l,i-1,j  ,k-1))
           kperpz4=0.125d0*(kperp(l,i+1,j+1,k  )+kperp(l,i+1,j  ,k  )+kperp(l,i+1,j+1,k-1) &
                +      kperp(l,i+1,j  ,k-1)+kperp(l,i  ,j+1,k  )+kperp(l,i  ,j  ,k  ) &
                +      kperp(l,i  ,j+1,k-1)+kperp(l,i  ,j  ,k-1))

           oneminuskperpz1 = 1.0d0-kperpz1
           oneminuskperpz2 = 1.0d0-kperpz2
           oneminuskperpz3 = 1.0d0-kperpz3
           oneminuskperpz4 = 1.0d0-kperpz4

        else if(streaming_diffusion)then
           kparaz1=0.125d0*(Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  )+Dpara(l,i  ,j  ,k-1) &
                +      Dpara(l,i  ,j-1,k-1)+Dpara(l,i-1,j  ,k  )+Dpara(l,i-1,j-1,k  ) &
                +      Dpara(l,i-1,j  ,k-1)+Dpara(l,i-1,j-1,k-1)) + kpar
           kparaz2=0.125d0*(Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j-1,k  )+Dpara(l,i+1,j  ,k-1) &
                +      Dpara(l,i+1,j-1,k-1)+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j-1,k  ) &
                +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i  ,j-1,k-1)) + kpar
           kparaz3=0.125d0*(Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  )+Dpara(l,i  ,j+1,k-1) &
                +      Dpara(l,i  ,j  ,k-1)+Dpara(l,i-1,j+1,k  )+Dpara(l,i-1,j  ,k  ) &
                +      Dpara(l,i-1,j+1,k-1)+Dpara(l,i-1,j  ,k-1)) + kpar
           kparaz4=0.125d0*(Dpara(l,i+1,j+1,k  )+Dpara(l,i+1,j  ,k  )+Dpara(l,i+1,j+1,k-1) &
                +      Dpara(l,i+1,j  ,k-1)+Dpara(l,i  ,j+1,k  )+Dpara(l,i  ,j  ,k  ) &
                +      Dpara(l,i  ,j+1,k-1)+Dpara(l,i  ,j  ,k-1)) + kpar
           kperpz1 = k_perp
           kperpz2 = k_perp
           kperpz3 = k_perp
           kperpz4 = k_perp
           oneminuskperpz1 = oneminuskperp
           oneminuskperpz2 = oneminuskperp
           oneminuskperpz3 = oneminuskperp
           oneminuskperpz4 = oneminuskperp
        else
           kparaz1 = kpar
           kparaz2 = kpar
           kparaz3 = kpar
           kparaz4 = kpar
           kperpz1 = k_perp
           kperpz2 = k_perp
           kperpz3 = k_perp
           kperpz4 = k_perp
           oneminuskperpz1 = oneminuskperp
           oneminuskperpz2 = oneminuskperp
           oneminuskperpz3 = oneminuskperp
           oneminuskperpz4 = oneminuskperp
        end if

        if(compute .ne. 3)then        
           fz1=kparaz1*(bz1*oneminuskperpz1*(bx1*dTdx1+by1*dTdy1+bz1*dTdz1)+kperpz1*dTdz1)
           fz2=kparaz2*(bz2*oneminuskperpz2*(bx2*dTdx2+by2*dTdy2+bz2*dTdz2)+kperpz2*dTdz2)
           fz3=kparaz3*(bz3*oneminuskperpz3*(bx3*dTdx3+by3*dTdy3+bz3*dTdz3)+kperpz3*dTdz3)
           fz4=kparaz4*(bz4*oneminuskperpz4*(bx4*dTdx4+by4*dTdy4+bz4*dTdz4)+kperpz4*dTdz4)
        else
           fz1=kparaz1*(bz1*oneminuskperpz1*(bx1+by1+bz1)+kperpz1)
           fz2=kparaz2*(bz2*oneminuskperpz2*(bx2+by2+bz2)+kperpz2)
           fz3=kparaz3*(bz3*oneminuskperpz3*(bx3+by3+bz3)+kperpz3)
           fz4=kparaz4*(bz4*oneminuskperpz4*(bx4+by4+bz4)+kperpz4)
        end if
        fz=0.25d0*(fz1+fz2+fz3+fz4)
#endif

        myflux(l,i,j,k)=fz*dt/dx
     enddo
  enddo
  enddo
  enddo
endif

end subroutine cmpZcrflx
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine L2_lim(a,b,a2,b2,alpha,l2result)
  use amr_parameters,only :dp
  implicit none
  real(dp)::a,b,a2,b2,alpha,l2result
  real(dp)::mean,aoveralpha,aalpha,mymin,mymax
  
  mean=0.5d0*(a+b)
  aoveralpha=a/alpha
  aalpha=a*alpha

  mymin=MIN(aalpha,aoveralpha)
  mymax=MAX(aalpha,aoveralpha)

  if    (mean.le.mymin)then
!!$     l2result=MIN(a2/alpha,a2*alpha)
     if(aoveralpha.lt.aalpha)then
        l2result=a2/alpha
     else
        l2result=a2*alpha
     endif
  elseif(mean.ge.mymax)then
!!$     l2result=MAX(a2/alpha,a2*alpha)
     if(aoveralpha.gt.aalpha)then
        l2result=a2/alpha
     else
        l2result=a2*alpha
     endif
  else
     l2result=0.5d0*(a2+b2)
  endif

end subroutine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine L4_lim(a,b,c,d,alpha,l4result)
  use amr_parameters,only :dp
  implicit none
  real(dp)::a,b,c,d,alpha,l4result
  real(dp)::mean,aoveralpha,aalpha,mymin,mymax
  
  mean=0.125d0*(a+b+c+d)
  aoveralpha=a/alpha
  aalpha=a*alpha

  mymin=MIN(aalpha,aoveralpha)
  mymax=MAX(aalpha,aoveralpha)

  l4result= MAX(mymin,MIN(mean,mymax))

end subroutine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine MinMod_lim(a,b,a2,b2,minmodresult,minmodresult2)
  use amr_parameters,only :dp
  implicit none
  real(dp)::a,b,a2,b2,minmodresult,minmodresult2
  
  if(a*b.le.0d0)then 
     minmodresult =0d0
     minmodresult2=0d0
  else
     if(a.gt.0.0d0)then
        if(a.lt.b)then
           minmodresult =a
           minmodresult2=a2
        else
           minmodresult =b
           minmodresult2=b2
        endif
     else
        if(a.lt.b)then
           minmodresult =b
           minmodresult2=b2
        else
           minmodresult =a
           minmodresult2=a2
        endif
     endif
  endif

end subroutine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine MonCen_lim(a,b,a2,b2,rr1,moncenresult)
  use amr_parameters,only :dp
  implicit none
  real(dp)::a,b,a2,b2,aa,bb,r1,moncenresult,rr1,aa2,bb2,r2
  
  call MinMod_lim(a,b,a2,b2,r1,r2)

  aa =2d0*r1
  aa2=2d0*r2
  bb =0.5d0*(a +b )
  bb2=0.5d0*(a2+b2)
  call MinMod_lim(aa,bb,aa2,bb2,rr1,moncenresult)

end subroutine
