
!############### MODULE ##############
module krome_commons
  implicit none

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2022-06-08 12:22:52
  !  Changeset d0f9dcf
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************
  integer,parameter::idx_E=1
  integer,parameter::idx_Hk=2
  integer,parameter::idx_H=3
  integer,parameter::idx_H2=4
  integer,parameter::idx_D=5
  integer,parameter::idx_HD=6
  integer,parameter::idx_HE=7
  integer,parameter::idx_Hj=8
  integer,parameter::idx_H2j=9
  integer,parameter::idx_Dj=10
  integer,parameter::idx_HEHj=11
  integer,parameter::idx_CR=12
  integer,parameter::idx_g=13
  integer,parameter::idx_Tgas=14
  integer,parameter::idx_dummy=15
  integer,parameter::nrea=25
  integer,parameter::nmols=11
  integer,parameter::nspec=15
  integer,parameter::natoms=4
  integer,parameter::ndust=0
  integer,parameter::ndustTypes=0
  integer,parameter::nPhotoBins=0
  integer,parameter::nPhotoRea=0

  !cooling index
  integer,parameter::idx_cool_h2 = 1
  integer,parameter::idx_cool_h2gp = 2
  integer,parameter::idx_cool_atomic = 3
  integer,parameter::idx_cool_cen = 3
  integer,parameter::idx_cool_hd = 4
  integer,parameter::idx_cool_z = 5
  integer,parameter::idx_cool_metal = 5
  integer,parameter::idx_cool_dh = 6
  integer,parameter::idx_cool_enthalpic = 6
  integer,parameter::idx_cool_dust = 7
  integer,parameter::idx_cool_compton = 8
  integer,parameter::idx_cool_cie = 9
  integer,parameter::idx_cool_continuum = 10
  integer,parameter::idx_cool_cont = 10
  integer,parameter::idx_cool_exp = 11
  integer,parameter::idx_cool_expansion = 11
  integer,parameter::idx_cool_ff = 12
  integer,parameter::idx_cool_bss = 12
  integer,parameter::idx_cool_custom = 13
  integer,parameter::idx_cool_co = 14
  integer,parameter::idx_cool_zcie = 15
  integer,parameter::idx_cool_zcienouv = 16
  integer,parameter::idx_cool_zextend = 17
  integer,parameter::idx_cool_gh = 18
  integer,parameter::idx_cool_oh = 19
  integer,parameter::idx_cool_h2o = 20
  integer,parameter::idx_cool_hcn = 21
  integer,parameter::ncools = 21

  !heating index
  integer,parameter::idx_heat_chem = 1
  integer,parameter::idx_heat_compress = 2
  integer,parameter::idx_heat_compr = 2
  integer,parameter::idx_heat_photo = 3
  integer,parameter::idx_heat_dh = 4
  integer,parameter::idx_heat_enthalpic = 4
  integer,parameter::idx_heat_photoav = 5
  integer,parameter::idx_heat_av = 5
  integer,parameter::idx_heat_cr = 6
  integer,parameter::idx_heat_dust = 7
  integer,parameter::idx_heat_xray = 8
  integer,parameter::idx_heat_visc = 9
  integer,parameter::idx_heat_viscous = 9
  integer,parameter::idx_heat_custom = 10
  integer,parameter::idx_heat_zcie = 11
  integer,parameter::nheats = 11

  real*8::arr_k(nrea)

  !commons for rate tables
  !modify ktab_n according to the required precision
  integer,parameter::ktab_n=int(1e3)
  real*8::ktab(nrea,ktab_n),ktab_logTlow, ktab_logTup, ktab_T(ktab_n)
  real*8::inv_ktab_T(ktab_n-1), inv_ktab_idx

  !thermo toggle (when >0 do cooling/heating)
  integer::krome_thermo_toggle
  !$omp threadprivate(krome_thermo_toggle)

  !debug bit flag, print and array with fallback values for extreme environments
  integer:: red_flag
  real*8::n_global(nspec)
  integer, save :: nprint_negative=10
  !$omp threadprivate(n_global,nprint_negative,red_flag)

  !commons for implicit RHS
  integer::arr_r1(nrea)
  integer::arr_r2(nrea)
  integer::arr_r3(nrea)
  integer::arr_p1(nrea)
  integer::arr_p2(nrea)
  integer::arr_p3(nrea)

  !commons for reduction
  integer::arr_u(nrea)
  real*8::arr_flux(nrea)

  !commons for frequency bins

  ! Draine dust absorption data loaded from file, via load_kabs
  ! in krome_photo module
  real*8::find_Av_draine_kabs(nPhotoBins)

  !commons for H2 photodissociation (Solomon)
  ! note: paramters here are set depending on the data
  ! but if you have a different file you should modify them
  integer,parameter::H2pdData_nvibX=15
  integer,parameter::H2pdData_nvibB=37
  real*8::H2pdData_dE(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_pre(H2pdData_nvibX,H2pdData_nvibB)
  real*8::H2pdData_EX(H2pdData_nvibX)
  integer::H2pdData_binMap(H2pdData_nvibX,H2pdData_nvibB)

  !commons for dust optical properties

  !square of turbulence velocity for broadening
  real*8::broadeningVturb2

  !mpi rank of process. If 0, ignored
  integer::krome_mpi_rank=0, krome_omp_thread
  !$omp threadprivate(krome_omp_thread)

  !user-defined commons variables from the reaction file
  real*8::user_Trad
  !$omp threadprivate(user_Trad)

  !commons for anytab

  !physical commons
  real*8::phys_Tcmb
  real*8::phys_zredshift
  real*8::phys_orthoParaRatio
  real*8::phys_metallicity
  real*8::phys_Tfloor
  !$omp threadprivate(phys_Tcmb)
  !$omp threadprivate(phys_zredshift)
  !$omp threadprivate(phys_orthoParaRatio)
  !$omp threadprivate(phys_metallicity)
  !$omp threadprivate(phys_Tfloor)

  !machine precision
  real*8::krome_epsilon

  !xrayJ21 for tabulated heating and rate
  real*8::J21xray

  !total metallicity relative to solar Z/Z_solar
  real*8::total_Z
  real*8::dust2gas_ratio

  !commons for dust tabs (cool,H2,Tdust)
  integer,parameter::dust_tab_imax=50, dust_tab_jmax=50
  real*8::dust_tab_ngas(dust_tab_imax)
  real*8::dust_tab_Tgas(dust_tab_jmax)
  real*8::dust_mult_Tgas,dust_mult_ngas
  real*8::dust_table_AvVariable_log

  real*8::dust_tab_cool(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_heat(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_Tdust(dust_tab_imax, dust_tab_jmax)
  real*8::dust_tab_H2(dust_tab_imax, dust_tab_jmax)

  !commons for exp(-a) table
  integer,parameter::exp_table_na=int(1d5)
  real*8,parameter::exp_table_aMax=1d4,exp_table_aMin=0d0
  real*8,parameter::exp_table_multa=(exp_table_na-1) &
      / (exp_table_aMax-exp_table_aMin)
  real*8,parameter::exp_table_da=1d0/exp_table_multa
  real*8::exp_table(exp_table_na)

  !stores the last evaluation of the rates in the fex
  real*8::last_coe(nrea)
  !$omp threadprivate(last_coe)

  !xsecs from file variables

  ! Gibbs free energy data from file variables

  !partition function from file
  integer,parameter::zpart_nCO=641
  integer,parameter::zpart_nH2even=2000
  integer,parameter::zpart_nH2odd=2000
  real*8::zpart_CO(zpart_nCO),minpart_CO,partdT_CO
  real*8::zpart_H2even(zpart_nH2even),minpart_H2even,partdT_H2even
  real*8::zpart_H2odd(zpart_nH2odd),minpart_H2odd,partdT_H2odd

  !Habing flux for the photoelectric heating by dust
  ! and clumping factor for H2 formation
  ! on dust by Jura/Gnedin
  real*8::GHabing,Ghabing_thin,clump_factor
  !$omp threadprivate(GHabing,GHabing_thin)

  !partition functions common vars

  !verbatim reactions
  character*50::reactionNames(nrea)

end module krome_commons

!############### MODULE ##############
module krome_constants
  implicit none

  !constants
  real*8,parameter::boltzmann_eV = 8.617332478d-5 !eV / K
  real*8,parameter::boltzmann_J = 1.380648d-23 !J / K
  real*8,parameter::boltzmann_erg = 1.380648d-16 !erg / K
  real*8,parameter::iboltzmann_eV = 1d0/boltzmann_eV !K / eV
  real*8,parameter::iboltzmann_erg = 1d0/boltzmann_erg !K / erg
  real*8,parameter::planck_eV = 4.135667516d-15 !eV s
  real*8,parameter::planck_J = 6.62606957d-34 !J s
  real*8,parameter::planck_erg = 6.62606957d-27 !erg s
  real*8,parameter::iplanck_eV = 1d0/planck_eV !1 / eV / s
  real*8,parameter::iplanck_J = 1d0/planck_J !1 / J / s
  real*8,parameter::iplanck_erg = 1d0/planck_erg !1 / erg / s
  real*8,parameter::gravity = 6.674d-8 !cm3 / g / s2
  real*8,parameter::e_mass = 9.10938188d-28 !g
  real*8,parameter::p_mass = 1.67262158d-24 !g
  real*8,parameter::n_mass = 1.674920d-24 !g
  real*8,parameter::ip_mass = 1d0/p_mass !1/g
  real*8,parameter::clight = 2.99792458e10 !cm/s
  real*8,parameter::pi = 3.14159265359d0 !#
  real*8,parameter::eV_to_erg = 1.60217646d-12 !eV -> erg
  real*8,parameter::ry_to_eV = 13.60569d0 !rydberg -> eV
  real*8,parameter::ry_to_erg = 2.179872d-11 !rydberg -> erg
  real*8,parameter::seconds_per_year = 365d0*24d0*3600d0 !yr -> s
  real*8,parameter::km_to_cm = 1d5 !km -> cm
  real*8,parameter::cm_to_Mpc = 1.d0/3.08d24 !cm -> Mpc
  real*8,parameter::kvgas_erg = 8.d0*boltzmann_erg/pi/p_mass !
  real*8,parameter::pre_kvgas_sqrt = sqrt(8.d0*boltzmann_erg/pi) !
  real*8,parameter::pre_planck = 2.d0*planck_erg/clight**2 !erg/cm2*s3
  real*8,parameter::exp_planck = planck_erg / boltzmann_erg !s*K
  real*8,parameter::stefboltz_erg = 5.670373d-5 !erg/s/cm2/K4
  real*8,parameter::N_avogadro = 6.0221d23 !#
  real*8,parameter::Rgas_J = 8.3144621d0 !J/K/mol
  real*8,parameter::Rgas_kJ = 8.3144621d-3 !kJ/K/mol
  real*8,parameter::hubble = 0.704d0 !dimensionless
  real*8,parameter::Omega0 = 1.0d0 !dimensionless
  real*8,parameter::Omegab = 0.0456d0 !dimensionless
  real*8,parameter::Hubble0 = 1.d2*hubble*km_to_cm*cm_to_Mpc !1/s

end module krome_constants

!############### MODULE ##############
module krome_fit
contains

  !*****************************
  subroutine init_anytab3D(filename,x,y,z,f,xmul,ymul,zmul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8::rout(4)
    integer::i,j,k,ios,unit

    !check the size of the X input array
    if(size(x).ne.size(f,1)) then
      print *,"ERROR: in init_anytab3D x size differs from f(x,y,z)"
      stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(f,2)) then
      print *,"ERROR: in init_anytab3D y size differs from f(x,y,z)"
      stop
    end if

    !check the size of the Z input array
    if(size(z).ne.size(f,3)) then
      print *,"ERROR: in init_anytab3D z size differs from f(x,y,z)"
      stop
    end if

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: in init_anytab3D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(unit,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
      print *,"ERROR: file "//filename//" should"
      print *," contain the number of grid points"
      print *," per dimension in the format"
      print *,"  XX, YY, ZZ"
      print *,row_string
      stop
    end if

    !loop to read file (3rd dimension of f() is
    ! first in the tables. i.e. tables are z,x,y,
    ! while f() is x,y,z
    do i=1,size(z)
      do j=1,size(x)
        do k=1,size(y)
          read(unit,*,iostat=ios) rout(:)
          y(k) = rout(3)
          f(j,k,i) = rout(4)
        end do
        x(j) = rout(2)
        read(unit,*,iostat=ios) !skip blanks
      end do
      z(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios.ne.0) exit
    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))
    zmul = 1d0/(z(2)-z(1))

  end subroutine init_anytab3D

  !********************************************
  !load 2d tables from filename
  subroutine init_anytab2D(filename,x,y,z,xmul,ymul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),z(:,:),xmul,ymul
    real*8::rout(3)
    integer::i,j,ios,unit

    !check the size of the X input array
    if(size(x).ne.size(z,1)) then
      print *,"ERROR: in init_anytab2D x size differs from z"
      stop
    end if

    !check the size of the Y input array
    if(size(y).ne.size(z,2)) then
      print *,"ERROR: in init_anytab2D y size differs from z"
      stop
    end if

    if (krome_mpi_rank<=1) print *,"Reading tables from "//trim(filename)

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: in init_anytab2D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(unit,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    !check if first line is OK
    if(scan(row_string,",")==0) then
      print *,"ERROR: file "//filename//" should"
      print *," contain the number of rows and "
      print *," columns in the format"
      print *,"  RR, CC"
      print *,row_string
      stop
    end if

    !loop to read file
    do i=1,size(x)
      do j=1,size(y)
        read(unit,*,iostat=ios) rout(:)
        y(j) = rout(2)
        z(i,j) = rout(3)
      end do
      x(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios.ne.0) exit
    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))
    ymul = 1d0/(y(2)-y(1))

  end subroutine init_anytab2D

  !********************************************
  !load 1d tables from filename
  subroutine init_anytab1D(filename,x,y,xmul)
    use krome_commons
    implicit none
    character(len=*),intent(in)::filename
    character(len=60)::row_string
    real*8,intent(out)::x(:),y(:),xmul
    real*8::rout(2)
    integer::i,ios,unit

    !check the size of the X input array
    if(size(x) /= size(y)) then
      print *,"ERROR: in init_anytab1D x size differs from y"
      stop
    end if

    if (krome_mpi_rank <= 1) print *,"Reading tables from "//trim(filename)

    !open file and check if it exists
    open(newunit=unit,file=trim(filename),status="old",iostat=ios)
    if(ios /= 0) then
      print *,"ERROR: in init_anytab1D file ",trim(filename)," not found!"
      stop
    end if

    !skip the comments and the first line and the sizes of the data
    ! which are already known from the pre-processing
    do
      read(unit,'(a)') row_string
      if(row_string(1:1)/="#") exit
    end do

    ! !check if first line is OK
    ! if(scan(row_string,",")==0) then
    !    print *,"ERROR: file "//filename//" should"
    !    print *," contain the number of rows and "
    !    print *," columns in the format"
    !    print *,"  RR, CC"
    !    print *,row_string
    !    stop
    ! end if

    !loop to read file
    do i=1,size(x)
      read(unit,*,iostat=ios) rout(:)
      y(i) = rout(2)
      x(i) = rout(1)
      read(unit,*,iostat=ios) !skip blanks
      if(ios /= 0) exit

    end do
    close(unit)

    xmul = 1d0/(x(2)-x(1))

  end subroutine init_anytab1D

  !******************************
  !test 2d fit and save to file
  subroutine test_anytab2D(fname,x,y,z,xmul,ymul)
    implicit none
    integer::i,j,unit1,unit2
    real*8,intent(in)::x(:),y(:),z(:,:),xmul,ymul
    real*8::xx,yy,zz
    character(len=*),intent(in)::fname

    open(newunit=unit1,file=fname//".fit",status="replace")
    open(newunit=unit2,file=fname//".org",status="replace")
    do i=1,size(x)
      do j=1,size(y)
        xx = x(i)
        yy = y(j)
        zz = fit_anytab2D(x(:),y(:),z(:,:),xmul,ymul,xx,yy)
        write(unit1,*) xx,yy,zz
        write(unit2,*) x(i),y(j),z(i,j)
      end do
      write(unit1,*)
      write(unit2,*)
    end do
    close(unit1)
    close(unit2)
    print *,"original file wrote in ",fname//".org"
    print *,"fit test file wrote in ",fname//".fit"

  end subroutine test_anytab2D

  !*****************************
  function fit_anytab3D(x,y,z,f,xmul,ymul,zmul,xx,yy,zz)
    implicit none
    real*8,intent(in)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8,intent(in)::xx,yy,zz
    real*8::fleft(size(x),size(y)), fright(size(x),size(y))
    real*8::fit_anytab3D,fl,fr
    integer::ipos,i1,i2

    ipos = (zz-z(1)) * zmul + 1
    i1 = min(max(ipos,1), size(z)-1)
    i2 = i1 + 1
    fleft(:,:) = f(:,:,i1)
    fright(:,:) = f(:,:,i2)

    fl = fit_anytab2D(x(:), y(:), fleft(:,:), xmul, ymul, xx, yy)
    fr = fit_anytab2D(x(:), y(:), fright(:,:), xmul, ymul, xx, yy)

    fit_anytab3D = (zz-z(i1))*zmul*(fr-fl)+fl

  end function fit_anytab3D

  !******************************
  !return 2d fit at xx,yy
  function fit_anytab2D(x,y,z,xmul,ymul,xx,yy)
    implicit none
    real*8::fit_anytab2D
    real*8,intent(in)::x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D(x(:),zright(:),xmul,xx)

    fit_anytab2D = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D

  !*********************
  !return 1d fit at xx
  function fit_anytab1D(x,z,xmul,xx)
    real*8,intent(in)::x(:),z(:),xmul,xx
    real*8::fit_anytab1D,p
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    fit_anytab1D = p * (z(i2) - z(i1)) + z(i1)

  end function fit_anytab1D

  !*****************************
  function fit_anytab3D_linlinlog(x,y,z,f,xmul,ymul,zmul,xx,yy,zz)
    implicit none
    real*8,intent(in)::x(:),y(:),z(:),f(:,:,:),xmul,ymul,zmul
    real*8,intent(in)::xx,yy,zz
    real*8::fleft(size(x),size(y)), fright(size(x),size(y))
    real*8::fit_anytab3D_linlinlog,fl,fr
    integer::ipos,i1,i2

    ipos = (zz-z(1)) * zmul + 1
    i1 = min(max(ipos,1), size(z)-1)
    i2 = i1 + 1
    fleft(:,:) = f(:,:,i1)
    fright(:,:) = f(:,:,i2)

    fl = fit_anytab2D_linlog(x(:), y(:), fleft(:,:), xmul, ymul, xx, yy)
    fr = fit_anytab2D_linlog(x(:), y(:), fright(:,:), xmul, ymul, xx, yy)

    fit_anytab3D_linlinlog = (zz-z(i1))*zmul*(fr-fl)+fl

  end function fit_anytab3D_linlinlog

  !***************************
  function fit_anytab2D_linlog(x,y,z,xmul,ymul,xx,yy)
    real*8::fit_anytab2D_linlog,x(:),y(:),z(:,:),xmul,ymul,xx,yy
    real*8::zleft(size(x)),zright(size(x)),zl,zr
    integer::ipos,i1,i2

    ipos = (yy-y(1)) * ymul + 1
    i1 = min(max(ipos,1),size(y)-1)
    i2 = i1 + 1
    zleft(:) = z(:,i1)
    zright(:) = z(:,i2)

    zl = fit_anytab1D_linlog(x(:),zleft(:),xmul,xx)
    zr = fit_anytab1D_linlog(x(:),zright(:),xmul,xx)

    fit_anytab2D_linlog = (yy-y(i1))*ymul*(zr-zl)+zl

  end function fit_anytab2D_linlog

  !*********************
  function fit_anytab1D_linlog(x,z,xmul,xx)
    real*8::fit_anytab1D_linlog,x(:),z(:),xmul,xx,p,z2,z1
    integer::ipos,i1,i2

    ipos = (xx-x(1)) * xmul + 1
    i1 = min(max(ipos,1),size(x)-1)
    i2 = i1 + 1

    p = (xx-x(i1)) * xmul

    z2 = z(i2)
    z1 = z(i1)
    if(z1<0d0 .and. z2<0d0) then
      z1 = log10(-z1)
      z2 = log10(-z2)
      fit_anytab1D_linlog = -1d1**(p * (z2 - z1) + z1)
      return
    end if

    if(z1>0d0 .and. z2>0d0) then
      z1 = log10(z1)
      z2 = log10(z2)
      fit_anytab1D_linlog = 1d1**(p * (z2 - z1) + z1)
      return
    end if

    fit_anytab1D_linlog = (p * (z2 - z1) + z1)

  end function fit_anytab1D_linlog

  !*****************************
  !spline interpolation at t using array  x,y (size n) as data
  function fspline(x,y,t)
    implicit none
    real*8::fspline,x(:),y(:),b(size(x)),c(size(x)),d(size(x)),t
    integer::n

    n = size(x)
    call spline(x(:),y(:),b(:),c(:),d(:),n)
    fspline = ispline(t,x(:),y(:),b(:),c(:),d(:),n)

  end function fspline

  !*******************************+
  subroutine spline(x, y, b, c, d, n)
    !======================================================================
    !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
    !  for cubic spline interpolation
    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    !  for  x(i) <= x <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    !  input..
    !  x = the arrays of data abscissas (in strictly increasing order)
    !  y = the arrays of data ordinates
    !  n = size of the arrays xi() and yi() (n>=2)
    !  output..
    !  b, c, d  = arrays of spline coefficients
    !  comments ...
    !  spline.f90 program is based on fortran version of program spline.f
    !  the accompanying function fspline can be used for interpolation
    !======================================================================
    implicit none
    integer::n
    real*8::x(n), y(n), b(n), c(n), d(n)
    integer::i, j, gap
    real*8::h

    gap = n-1

    !check input
    if(n<2) return
    if(n<3) then
      b(1) = (y(2)-y(1))/(x(2)-x(1)) !linear interpolation
      c(1) = 0d0
      d(1) = 0d0
      b(2) = b(1)
      c(2) = 0d0
      d(2) = 0d0
      return
    end if

    !step 1: preparation
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2d0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
    end do

    ! step 2: end conditions
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0d0
    c(n) = 0d0
    if(n.ne.3) then
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if

    ! step 3: forward elimination
    do i = 2, n
      h = d(i-1)/b(i-1)
      b(i) = b(i) - h*d(i-1)
      c(i) = c(i) - h*c(i-1)
    end do

    ! step 4: back substitution
    c(n) = c(n)/b(n)
    do j = 1, gap
      i = n-j
      c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do

    ! step 5: compute spline coefficients
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2d0*c(n))
    do i = 1, gap
      b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2d0*c(i))
      d(i) = (c(i+1) - c(i))/d(i)
      c(i) = 3d0*c(i)
    end do
    c(n) = 3d0*c(n)
    d(n) = d(n-1)
  end subroutine spline

  !*******************************
  function ispline(u, x, y, b, c, d, n)
    !======================================================================
    ! function ispline evaluates the cubic spline interpolation at point z
    ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
    ! where  x(i) <= u <= x(i+1)
    !  Alexadner L Godunov (ODU): January 2010
    !
    !  http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90
    !----------------------------------------------------------------------
    ! input..
    ! u       = the abscissa at which the spline is to be evaluated
    ! x, y    = the arrays of given data points
    ! b, c, d = arrays of spline coefficients computed by spline
    ! n       = the number of data points
    ! output:
    ! ispline = interpolated value at point u
    !=======================================================================
    implicit none
    real*8::ispline
    integer::n
    real*8::u, x(n), y(n), b(n), c(n), d(n)
    integer::i, j, k
    real*8::dx

    ! if u is ouside the x() interval take a boundary value (left or right)
    if(u<=x(1)) then
      ispline = y(1)
      return
    end if

    if(u>=x(n)) then
      ispline = y(n)
      return
    end if

    ! binary search for for i, such that x(i) <= u <= x(i+1)
    i = 1
    j = n+1
    do while (j>i+1)
      k = (i+j)/2
      if(u<x(k)) then
        j=k
      else
        i=k
      end if
    end do

    ! evaluate spline interpolation
    dx = u - x(i)
    ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

  end function ispline

end module krome_fit
!This module contains useful routines to get physical
! quantities, like mean molecular weight, mass density,
! mass, jeans length, etc. etc.

!############### MODULE ##############
module krome_getphys
contains

  !*****************************
  !get the mean molecular weight
  function get_mu(n)
    use krome_commons
    use krome_constants
    implicit none
    real*8::n(:),get_mu,m(nspec)
    m(:) = get_mass()

    !ip_mass is 1/proton_mass_in_g
    get_mu = max(sum(n(1:nmols)*m(1:nmols)),1d-40) &
        / max(sum(n(1:nmols)),1d-40) * ip_mass

  end function get_mu

  !***************************
  !get mean molecular weight
  function get_mu_rho(n,rhogas)
    use krome_commons
    use krome_constants
    implicit none
    real*8::get_mu_rho,rhogas,n(:)

    !ip_mass is 1/proton_mass_in_g
    get_mu_rho = rhogas / max(sum(n(1:nmols)),1d-40) * ip_mass

  end function get_mu_rho

  !************************
  !get species masses (g)
  function get_mass()
    use krome_commons
    implicit none
    real*8::get_mass(nspec)

    get_mass(1) = 9.10938188d-28	!E
    get_mass(2) = 1.6744434563759998d-24	!H-
    get_mass(3) = 1.673532518188d-24	!H
    get_mass(4) = 3.347065036376d-24	!H2
    get_mass(5) = 3.346032518188d-24	!D
    get_mass(6) = 5.0195650363760004d-24	!HD
    get_mass(7) = 6.692065036376d-24	!HE
    get_mass(8) = 1.67262158d-24	!H+
    get_mass(9) = 3.346154098188d-24	!H2+
    get_mass(10) = 3.34512158d-24	!D+
    get_mass(11) = 8.364686616375999d-24	!HEH+
    get_mass(12) = 0.d0	!CR
    get_mass(13) = 0.d0	!g
    get_mass(14) = 0.d0	!Tgas
    get_mass(15) = 0.d0	!dummy

  end function get_mass

  !************************
  !get sqrt of the inverse of the masses (1/sqrt(g))
  function get_imass_sqrt()
    use krome_commons
    implicit none
    real*8::get_imass_sqrt(nspec)

    get_imass_sqrt(1) = 33132602150543.92	!E
    get_imass_sqrt(2) = 772795806394.0071	!H-
    get_imass_sqrt(3) = 773006102110.9268	!H
    get_imass_sqrt(4) = 546597856701.2171	!H2
    get_imass_sqrt(5) = 546682184736.54565	!D
    get_imass_sqrt(6) = 446341179966.5731	!HD
    get_imass_sqrt(7) = 386562679981.0883	!HE
    get_imass_sqrt(8) = 773216569600.4055	!H+
    get_imass_sqrt(9) = 546672253002.7066	!H2+
    get_imass_sqrt(10) = 546756615481.1794	!D+
    get_imass_sqrt(11) = 345760328884.13617	!HEH+
    get_imass_sqrt(12) = 0.d0	!CR
    get_imass_sqrt(13) = 0.d0	!g
    get_imass_sqrt(14) = 0.d0	!Tgas
    get_imass_sqrt(15) = 0.d0	!dummy

  end function get_imass_sqrt

  !************************
  !get inverse of the species masses (1/g)
  function get_imass()
    use krome_commons
    implicit none
    real*8::get_imass(nspec)

    get_imass(1) = 1.0977693252662275d+27	!E
    get_imass(2) = 5.9721335838016375d+23	!H-
    get_imass(3) = 5.9753843390072855d+23	!H
    get_imass(4) = 2.9876921695036427d+23	!H2
    get_imass(5) = 2.9886141110832266d+23	!D
    get_imass(6) = 1.9922044893395282d+23	!HD
    get_imass(7) = 1.4943070555416133d+23	!HE
    get_imass(8) = 5.978638635046188d+23	!H+
    get_imass(9) = 2.988505522030552d+23	!H2+
    get_imass(10) = 2.9894279657243426d+23	!D+
    get_imass(11) = 1.1955020503006603d+23	!HEH+
    get_imass(12) = 0.d0	!CR
    get_imass(13) = 0.d0	!g
    get_imass(14) = 0.d0	!Tgas
    get_imass(15) = 0.d0	!dummy

  end function get_imass

  !************************
  !species binding energies (surface=BARE), K
  function get_EbindBare()
    use krome_commons
    implicit none
    real*8::get_EbindBare(nspec)

    get_EbindBare(:) = 1d99

    get_EbindBare(idx_H) = 500.0d0
    get_EbindBare(idx_H2) = 300.0d0

  end function get_EbindBare

  !************************
  !species binding energies (surface=ICE), K
  function get_EbindIce()
    use krome_commons
    implicit none
    real*8::get_EbindIce(nspec)

    get_EbindIce(:) = 1d99

    get_EbindIce(idx_H) = 650.0d0
    get_EbindIce(idx_H2) = 300.0d0

  end function get_EbindIce

  !************************
  function get_kevap70()
    use krome_commons
    implicit none
    real*8::get_kevap70(nspec)

    get_kevap70(idx_E) = 0d0
    get_kevap70(idx_Hk) = 0d0
    get_kevap70(idx_H) = 790490323.1199661
    get_kevap70(idx_H2) = 13763786733.050402
    get_kevap70(idx_D) = 0d0
    get_kevap70(idx_HD) = 0d0
    get_kevap70(idx_HE) = 0d0
    get_kevap70(idx_Hj) = 0d0
    get_kevap70(idx_H2j) = 0d0
    get_kevap70(idx_Dj) = 0d0
    get_kevap70(idx_HEHj) = 0d0
    get_kevap70(idx_CR) = 0d0
    get_kevap70(idx_g) = 0d0
    get_kevap70(idx_Tgas) = 0d0
    get_kevap70(idx_dummy) = 0d0

  end function get_kevap70

  !************************
  !get verbatim reaction names
  function get_rnames()
    use krome_commons
    implicit none
    character*50::get_rnames(nrea)

    !reaction names are loaded from file
    get_rnames(:) = reactionNames(:)

  end function get_rnames

  !************************
  !get species names
  function get_names()
    use krome_commons
    implicit none
    character*16::get_names(nspec)

    get_names(1) = "E"
    get_names(2) = "H-"
    get_names(3) = "H"
    get_names(4) = "H2"
    get_names(5) = "D"
    get_names(6) = "HD"
    get_names(7) = "HE"
    get_names(8) = "H+"
    get_names(9) = "H2+"
    get_names(10) = "D+"
    get_names(11) = "HEH+"
    get_names(12) = "CR"
    get_names(13) = "g"
    get_names(14) = "Tgas"
    get_names(15) = "dummy"

  end function get_names

  !************************
  !get cooling names list (empty element if cooling not present)
  function get_cooling_names()
    use krome_commons
    implicit none
    character*16::get_cooling_names(ncools)

    get_cooling_names(:) = ""

    get_cooling_names(idx_cool_h2) = "H2"
    get_cooling_names(idx_cool_h2gp) = "H2GP"
    get_cooling_names(idx_cool_atomic) = "ATOMIC"
    get_cooling_names(idx_cool_cen) = "CEN"
    get_cooling_names(idx_cool_hd) = "HD"
    get_cooling_names(idx_cool_z) = "Z"
    get_cooling_names(idx_cool_metal) = "METAL"
    get_cooling_names(idx_cool_dh) = "DH"
    get_cooling_names(idx_cool_enthalpic) = "ENTHALPIC"
    get_cooling_names(idx_cool_dust) = "DUST"
    get_cooling_names(idx_cool_compton) = "COMPTON"
    get_cooling_names(idx_cool_cie) = "CIE"
    get_cooling_names(idx_cool_continuum) = "CONTINUUM"
    get_cooling_names(idx_cool_cont) = "CONT"
    get_cooling_names(idx_cool_exp) = "EXP"
    get_cooling_names(idx_cool_expansion) = "EXPANSION"
    get_cooling_names(idx_cool_ff) = "FF"
    get_cooling_names(idx_cool_bss) = "BSS"
    get_cooling_names(idx_cool_custom) = "CUSTOM"
    get_cooling_names(idx_cool_co) = "CO"
    get_cooling_names(idx_cool_zcie) = "ZCIE"
    get_cooling_names(idx_cool_zcienouv) = "ZCIENOUV"
    get_cooling_names(idx_cool_zextend) = "ZEXTEND"
    get_cooling_names(idx_cool_gh) = "GH"
    get_cooling_names(idx_cool_oh) = "OH"
    get_cooling_names(idx_cool_h2o) = "H2O"
    get_cooling_names(idx_cool_hcn) = "HCN"

  end function get_cooling_names

  !************************
  !get heating names list (empty element if heating not present)
  function get_heating_names()
    use krome_commons
    implicit none
    character*16::get_heating_names(nheats)

    get_heating_names(:) = ""

    get_heating_names(idx_heat_chem) = "CHEM"
    get_heating_names(idx_heat_compress) = "COMPRESS"
    get_heating_names(idx_heat_compr) = "COMPR"
    get_heating_names(idx_heat_photo) = "PHOTO"
    get_heating_names(idx_heat_dh) = "DH"
    get_heating_names(idx_heat_enthalpic) = "ENTHALPIC"
    get_heating_names(idx_heat_photoav) = "PHOTOAV"
    get_heating_names(idx_heat_av) = "AV"
    get_heating_names(idx_heat_cr) = "CR"
    get_heating_names(idx_heat_dust) = "DUST"
    get_heating_names(idx_heat_xray) = "XRAY"
    get_heating_names(idx_heat_visc) = "VISC"
    get_heating_names(idx_heat_viscous) = "VISCOUS"
    get_heating_names(idx_heat_custom) = "CUSTOM"
    get_heating_names(idx_heat_zcie) = "ZCIE"

  end function get_heating_names

  !******************************
  !get the total number of H nuclei
  function get_Hnuclei(n)
    use krome_commons
    real*8::n(:),get_Hnuclei,nH

    nH = n(idx_Hk) + &
        n(idx_H) + &
        n(idx_H2)*2d0 + &
        n(idx_HD) + &
        n(idx_Hj) + &
        n(idx_H2j)*2d0 + &
        n(idx_HEHj)
    get_Hnuclei = nH

  end function get_Hnuclei

  !***************************
  function get_zatoms()
    use krome_commons
    implicit none
    integer::get_zatoms(nspec)

    get_zatoms(1) = 0	!E
    get_zatoms(2) = 1	!H-
    get_zatoms(3) = 1	!H
    get_zatoms(4) = 2	!H2
    get_zatoms(5) = 1	!D
    get_zatoms(6) = 2	!HD
    get_zatoms(7) = 2	!HE
    get_zatoms(8) = 1	!H+
    get_zatoms(9) = 2	!H2+
    get_zatoms(10) = 1	!D+
    get_zatoms(11) = 3	!HEH+
    get_zatoms(12) = 0	!CR
    get_zatoms(13) = 0	!g
    get_zatoms(14) = 0	!Tgas
    get_zatoms(15) = 0	!dummy

  end function get_zatoms

  !******************************
  function get_qeff()
    use krome_commons
    implicit none
    real*8::get_qeff(nrea)

    get_qeff(:) = 0e0

  end function get_qeff

  !**************************
  function get_free_fall_time(n)
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),m(nspec)
    real*8::rhogas,get_free_fall_time

    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols))
    get_free_fall_time = sqrt(3d0*pi/32d0/gravity/rhogas)

  end function get_free_fall_time

  !**************************
  function get_free_fall_time_rho(rhogas)
    use krome_constants
    implicit none
    real*8::rhogas,get_free_fall_time_rho

    get_free_fall_time_rho = sqrt(3d0*pi/32d0/gravity/rhogas)

  end function get_free_fall_time_rho

  !********************************
  function get_jeans_length(n,Tgas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::m(nspec),get_jeans_length
    m(:) = get_mass()
    rhogas = max(sum(n(1:nmols)*m(1:nmols)),1d-40)
    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length

  !********************************
  function get_jeans_length_rho(n,Tgas,rhogas)
    !get jeans length in cm
    use krome_constants
    use krome_commons
    implicit none
    real*8::n(:),Tgas,mu,rhogas
    real*8::get_jeans_length_rho

    mu = get_mu_rho(n(:),rhogas)
    get_jeans_length_rho = sqrt(pi*boltzmann_erg*Tgas/rhogas&
        /p_mass/gravity/mu)

  end function get_jeans_length_rho

  !***************************
  !number density to column density conversion
  function num2col(ncalc,n)
    use krome_commons
    implicit none
    real*8::num2col,ncalc,n(:),Tgas
    Tgas = max(n(idx_Tgas),phys_Tcmb)

    num2col = 1.87d21*(max(ncalc,1d-40)*1d-3)**(2./3.)

  end function num2col

  !***********************
  !column density to number density conversion
  function col2num(ncalc,n)
    use krome_commons
    implicit none
    real*8::col2num,ncalc,n(:),Tgas
    Tgas = max(n(idx_Tgas),phys_Tcmb)

    col2num = 1d3 * (max(ncalc,1d-40)/1.87d21)**1.5

  end function col2num

  !************************
  !get electrons by balancing charges
  function get_electrons(n)
    use krome_commons
    implicit none
    real*8::get_electrons,n(nspec)

    get_electrons =  - n(idx_Hk) &
        + n(idx_Hj) &
        + n(idx_H2j) &
        + n(idx_Dj) &
        + n(idx_HEHj)
    get_electrons = max(get_electrons,0d0)

  end function get_electrons

  !************************
  !get species charges
  function get_charges()
    use krome_commons
    implicit none
    integer::get_charges(nspec)

    get_charges(1) = -1.d0 	!E
    get_charges(2) = -1.d0 	!H-
    get_charges(3) = 0.d0 	!H
    get_charges(4) = 0.d0 	!H2
    get_charges(5) = 0.d0 	!D
    get_charges(6) = 0.d0 	!HD
    get_charges(7) = 0.d0 	!HE
    get_charges(8) = 1.d0 	!H+
    get_charges(9) = 1.d0 	!H2+
    get_charges(10) = 1.d0 	!D+
    get_charges(11) = 1.d0 	!HEH+
    get_charges(12) = 0.d0 	!CR
    get_charges(13) = 0.d0 	!g
    get_charges(14) = 0.d0 	!Tgas
    get_charges(15) = 0.d0 	!dummy

  end function get_charges

end module krome_getphys
!This module contains the functions and subroutines
! needed to evaluate the adiabatic index.

!############### MODULE ##############
module krome_gadiab
contains

  !#KROME_header

  !**************************
  !compute 1/(gamma-1) at Tgasin using the partition function
  ! provided in the array_part with a temperature step dT_part
  ! and a minimum Tgas value min_part
  function gamma_pop(array_part,dT_part,min_part,Tgasin)
    implicit none
    real*8::array_part(:),dT_part
    real*8::min_part,Tgas,gamma_pop,Tgas2,Tgasin
    real*8::logz,logz1,logz2,emed1,emed2,Cv,inTgas,T2,T1,Cv1,Cv2
    integer::idx

    !temperature above minimum data point
    inTgas = max(Tgasin,min_part)

    !data index
    idx = (inTgas-min_part)/dT_part+1
    !corresponding Tgas
    Tgas = (idx-1)*dT_part+min_part
    !store Tgas
    T1 = Tgas

    !ln of partition functions (3 points forward)
    logz = log(array_part(idx))
    logz1 = log(array_part(idx+1))
    logz2 = log(array_part(idx+2))

    !derivative for mean energy (2 points forward)
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !derivative for 1/(gamma-1)
    Cv1 = (emed2-emed1)/dT_part

    !next point temperature
    Tgas = (idx)*dT_part+min_part
    !store Tgas
    T2 = Tgas
    !ln of partition functions
    logz = logz1
    logz1 = logz2
    logz2 = log(array_part(idx+3))

    !derivative for mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !derivative for 1/(gamma-1)
    Cv2 = (emed2-emed1)/dT_part

    !interpolation for 1/(gamma-1)
    Cv = (Cv2-Cv1)*(inTgas-T1)/(T2-T1)+Cv1

    !returns result
    gamma_pop = Cv

  end function gamma_pop

  !*****************************
  !compute 1/(gamma-1) at Tgasin using the partition function
  ! provided in the array_part with a temperature step dT_part
  ! and a minimum Tgas value min_part, for H2 with a ortho/para
  ! ratio of opratio. Needs even and odd partition functions.
  function gamma_pop_H2(array_part_even,array_part_odd,dT_part,&
        min_part,Tgasin,opratio)
    implicit none
    real*8::array_part_even(:),array_part_odd(:),dT_part,zcut(4)
    real*8::min_part,Tgas,opratio,gamma_pop_H2,Tgas2,a,b,Tgasin
    real*8::logz,logz1,logz2,emed1,emed2,Cv,inTgas,T2,T1,Cv1,Cv2
    integer::idx

    !Tgas above the data limit
    inTgas = max(Tgasin,min_part)

    !exponents for ortho/para ratio
    a = opratio/(opratio+1d0) !exponent zo
    b = 1d0-a !exponent zp

    !index in the data for the given Tgas
    idx = (inTgas-min_part)/dT_part+1
    !get the corresponding Tgas
    Tgas = (idx-1)*dT_part+min_part
    !store Tgas
    T1 = Tgas

    !needed for ortho partition function (see Boley+2007)
    zcut(1) = exp(2d0*85.4/Tgas)
    zcut(2) = exp(2d0*85.4/(Tgas+dT_part))
    zcut(3) = exp(2d0*85.4/(Tgas+2d0*dT_part))
    zcut(4) = exp(2d0*85.4/(Tgas+3d0*dT_part))

    !ln of the composite partition function
    logz = log(array_part_even(idx)**b*(3d0*array_part_odd(idx)*zcut(1))**a)
    logz1 = log(array_part_even(idx+1)**b*(3d0*array_part_odd(idx+1)*zcut(2))**a)
    logz2 = log(array_part_even(idx+2)**b*(3d0*array_part_odd(idx+2)*zcut(3))**a)
    !derivative for mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !get 1/(gamma-1) for the left point
    Cv1 = (emed2-emed1)/dT_part

    !Tgas of the right point
    Tgas = (idx)*dT_part+min_part
    !store Tgas
    T2 = Tgas
    !ln of the composite function
    logz = logz1
    logz1 = logz2
    logz2 = log(array_part_even(idx+3)**b*(3d0*array_part_odd(idx+3)*zcut(4))**a)
    !derivative for the mean energy
    emed1 = Tgas**2*(logz1-logz)/dT_part
    emed2 = (Tgas+dT_part)**2*(logz2-logz1)/dT_part

    !get 1/(gamma-1) for the right point
    Cv2 = (emed2-emed1)/dT_part

    !interpolation of 1/(gamma-1)
    Cv = (Cv2-Cv1)*(inTgas-T1)/(T2-T1)+Cv1

    !returns the result
    gamma_pop_H2 = Cv
  end function gamma_pop_H2

  !**************************
  !function to get the partition function
  ! of H2 at Tgas with a orto-para ratio
  ! equal to opratio
  function zfop(Tgas,opratio)
    implicit none
    real*8::Tgas,zfop,brot,ibTgas
    real*8::a,b,zo,zp,opratio
    integer::j,jmax,j1
    brot = 85.4d0 !H2 rotational constant in K
    zo = 0d0 !sum for ortho partition function
    zp = 0d0 !sum for para partition function
    jmax = 10 !number of terms in sum

    ibTgas = brot/Tgas !pre-calc

    !loop over levels
    do j=0,jmax,2 !step 2
      j1 = j + 1
      zp = zp + (2d0*j+1d0) * exp(-j*(j+1d0)*ibTgas)
      zo = zo + 3d0 * (2d0*j1+1d0) * exp(-j1*(j1+1d0)*ibTgas)
    end do

    a = opratio/(opratio+1d0) !exponent zo
    b = 1d0-a !exponent zp

    zfop = (zp**b * (zo*exp(2d0*ibTgas))**a) !final partition f

  end function zfop

  !*********************
  !get the partition function at Tgas
  ! of a diatom with rotational constant
  ! brot in K
  function zf(Tgas,brot)
    real*8::Tgas,zf,brot,z,ibTgas
    integer::j,jmax
    jmax = 10 !number of levels

    ibTgas = brot/Tgas !store
    z = 0d0
    !loop on levels
    do j=0,jmax
      z = z + (2d0*j+1d0)*exp(-j*(j+1d0)*ibTgas)
    end do

    zf = z

  end function zf

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of H2 with
  ! an ortho-para ratio of opratio
  function gamma_rotop(Tgas_in,opratio)
    implicit none
    real*8::gamma_rotop,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,opratio

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zfop(Tgas+dT,opratio)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zfop(Tgas,opratio)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zfop(Tgas+dT+dT,opratio))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rotop = (prot2-prot1)*idT

  end function gamma_rotop

  !***********************
  !get the degrees of freedom at Tgas for
  ! the rotational component of a diatom
  ! with rotational constant brot in K
  function gamma_rot(Tgas_in,brot)
    implicit none
    real*8::gamma_rot,Tgas,dT,Tgas_in
    real*8::idT,dlog1,prot1,dlog2,prot2
    real*8::logp1,brot

    Tgas = max(Tgas_in,1d1)

    dT = Tgas*1d-5 !dT for derivative
    idT =  1d0/dT !stored for numeric derivative
    logp1 = log(zf(Tgas+dT,brot)) !store since used twice

    !derivative dlog(T)/dT = f(T)
    dlog1 = (logp1-log(zf(Tgas,brot)))*idT
    prot1 = dlog1*Tgas**2

    !derivative dlog(T+dT)/dT = f(T+dT)
    dlog2 = (log(zf(Tgas+dT+dT,brot))-logp1)*idT
    prot2 = dlog2*(Tgas+dT)**2

    !derivative df(T)/dT
    gamma_rot = (prot2-prot1)*idT

  end function gamma_rot

  !*********************
  !get gamma
  function gamma_index(n)
    use krome_commons
    implicit none
    real*8::n(:),gamma_index,krome_gamma

    krome_gamma = (5.d0*(n(idx_E) + n(idx_H) + n(idx_HE)) + 7.d0*(n(idx_H2))) / &
        (3.d0*(n(idx_E) + n(idx_H) + n(idx_HE)) + 5.d0*(n(idx_H2)))

    gamma_index = krome_gamma
  end function gamma_index

end module krome_gadiab
!This module contains functions and subroutines
! for the surface chemistry, including adsorption, desorption, chemisorption
! and icy grains.

!############### MODULE ##############
module krome_grfuncs
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2022-06-08 12:22:52
  !  Changeset d0f9dcf
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !**********************
  !get Tdust from tables, K
  function get_table_Tdust(n) result(Tdust)
    use krome_commons
    use krome_fit
    implicit none
    real*8,intent(in)::n(nspec)
    real*8::ntot,Tdust,Tgas

    Tgas = n(idx_Tgas)

    !default, K
    Tdust = 1d0

    !total densitym, cm-3
    ntot = sum(n(1:nmols))

    !zero density returns default
    if(ntot==0d0) return

    !get dust temperature from table, K
    Tdust = 1d1**fit_anytab2D(dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_Tdust(:,:), dust_mult_ngas, &
        dust_mult_Tgas, &
        log10(ntot), log10(Tgas))

  end function get_table_Tdust

  !**********************
  !adsorpion rate Hollenbach+McKee 1979, Cazaux+2010, Hocuk+2014
  function dust_adsorption_rate(nndust,ims,stick,adust2,sqrTgas)
    use krome_constants
    implicit none
    real*8::dust_adsorption_rate,nndust,ims,stick,adust2,sqrTgas

    dust_adsorption_rate = nndust * pi * adust2 &
        * pre_kvgas_sqrt * ims * sqrTgas &
        * stick

  end function dust_adsorption_rate

  !*****************************
  !desorption rate Cazaux+2010, Hocuk+2014
  function dust_desorption_rate(fice,expEice,expEbare)
    implicit none
    real*8::dust_desorption_rate
    real*8::fice,expEice,expEbare,nu0,fbare

    nu0 = 1d12 !1/s
    fbare = 1d0 - fice
    dust_desorption_rate = nu0 * (fbare * expEbare &
        + fice * expEice)

  end function dust_desorption_rate

  !**************************
  function dust_2body_rate(p,invphi,fice,expEice1,expEice2,&
        expEbare1,expEbare2,pesc_ice,pesc_bare)
    use krome_constants
    implicit none
    real*8::fice,expEice1,expEice2,expEbare1,expEbare2,invphi
    real*8::nu0,p,dust_2body_rate,fbare,pesc_ice,pesc_bare

    !no need to calculate this if the dust is not present
    dust_2body_rate = 0d0

    fbare = 1d0-fice
    nu0 = 1d12 ! 1/s
    dust_2body_rate = fbare * (expEbare1 + expEbare2) * pesc_bare &
        + fice * (expEice1 + expEice2) * pesc_ice
    dust_2body_rate = dust_2body_rate * p * nu0 * invphi

  end function dust_2body_rate

  !******************
  function krate_2bodySi(alpha,Ea,idx1,idx2,n,Tdust) result(krate)
    use krome_commons
    implicit none
    real*8,intent(in)::n(nspec),Ea,Tdust,alpha
    integer,intent(in)::idx1,idx2
    real*8::krate,amin,amax,pexp,d2g,rho0

    !some default values OK for silicates
    amin = 5d-7 !cm
    amax = 2.5d-5 !cm
    pexp = -3.5
    rho0 = 3d0 !g/cm3
    d2g = 1d-2

    krate = krate_2body(n(:),idx1,idx2,alpha,amin,amax,pexp,d2g,rho0,Ea,Tdust)

  end function krate_2bodySi

  !********************
  !This routine has been modified to accomodate
  !Semenov framework  for surface chemistry.
  function krate_2body(n,idx1,idx2,alpha,amin,amax,pexp,d2g,rho0, &
        Ea,Tdust) result(krate)
    use krome_commons
    use krome_constants
    use krome_getphys
    implicit none
    integer,intent(in)::idx1,idx2
    real*8,intent(in)::n(nspec),amin,amax,pexp,d2g,rho0,Ea,Tdust,alpha
    real*8::rhog,p3,p4,ndns,krate,mred,fice,fbare,Preac,p
    real*8::iTd23,Ebare(nspec),Eice(nspec),mass(nspec)
    real*8,parameter::app2=(3d-8)**2 !cm^2 (Hocuk+2015)
    real*8,parameter::nu0=1d12 !1/s
    real*8,parameter::hbar=planck_erg/2d0/pi !erg*s
    real*8,parameter::ar=1d-8 !cm

    mass(:) = get_mass()

    !gas density, g/cm3
    rhog = sum(mass(1:nmols)*n(1:nmols))

    !exponentes
    p3 = pexp + 3d0
    p4 = pexp + 4d0

    !number of sites cm-3/mly
    ndns = rhog/(4d0/3d0*rho0*app2)*(amax**p3-amin**p3) &
        / (amax**p4-amin**p4) * p4 / p3

    !reduced mass
    mred = mass(idx1)*mass(idx2)/(mass(idx1)+mass(idx2))

    !tunneling probability
    Preac = exp(-2d0*ar/hbar*sqrt(2d0*mred*Ea*boltzmann_erg))
    !exponent
    iTd23 = 2d0/3d0/Tdust

    !get Ebind, K
    Ebare(:) = get_Ebind_bare()

    !ice/bare fraction
    fbare = 1d0

    !compute rate
    krate = fbare*(exp(-Ebare(idx1)*iTd23)+exp(-Ebare(idx2)*iTd23))

    !rate in cm3/s
    krate = nu0*Preac/ndns*krate

  end function krate_2body

  !*************************
  function dust_get_inv_phi(asize2,nndust)
    use krome_commons
    use krome_constants
    implicit none
    real*8::iapp2,dust_get_inv_phi(ndust),asize2(ndust)
    real*8::nndust(ndust),dephi
    integer::i

    iapp2 = (3d-8)**2 !1/cm2
    do i=1,ndust
      dust_get_inv_phi(i) = 0d0
      dephi = (4d0 * nndust(i) * pi * asize2(i))
      if(dephi.le.0d0) cycle
      dust_get_inv_phi(i) = iapp2 / dephi
    end do

  end function dust_get_inv_phi

  !****************************
  !returns an array with the sticking coefficient for each bin
  ! following Hollenbach+McKee 1979
  function dust_stick_array(Tgas,Tdust)
    use krome_commons
    implicit none
    real*8::dust_stick_array(ndust),Tgas,Tdust(ndust)
    real*8::Tg100,Td100
    integer::i

    Tg100 = Tgas * 1d-2
    do i=1,ndust
      Td100 = Tdust(i) * 1d-2
      dust_stick_array(i) = 1d0/(1d0+.4d0*sqrt(Tg100+Td100) &
          + .2d0*Tg100 + 0.08d0*Tg100**2)
    end do

  end function dust_stick_array

  !*************************
  function dust_stick(Tgas,Tdust)
    implicit none
    real*8,intent(in)::Tgas,Tdust
    real*8::dust_stick
    real*8::Tg100,Td100

    Tg100 = Tgas * 1d-2
    Td100 = Tdust * 1d-2
    dust_stick = 1d0/(1d0 + 0.4d0*sqrt(Tg100+Td100) &
        + 0.2d0*Tg100 + 0.08d0*Tg100**2)

  end function dust_stick

  !****************************
  !sticking rate (1/s), assuming power-law dust distribution
  ! example rate is
  !  @format:idx,R,P,rate
  !  1,CO,CO_ice,krate_stick(n(:),idx_CO,1d-7,1d-5,-3.5,3d0,1d-2)
  ! n(:): internal status array (number densities, temeperature, etc...)
  ! idx : index of the sticking species, e.g. idx_CO
  ! Tdust: dust temperature (assume same for all bins), K
  ! amin: min grain size, cm
  ! amax: max grain size, cm
  ! pexp: power-law exponent, usually -3.5
  ! rho0: bulk material density, g/cm3, e.g. 3 g/cm3 for silicates
  ! d2g: dust to gass mass ratio, usually 0.01
  function krate_stick(n,idx,Tdust,amin,amax,pexp,rho0,d2g) result(k)
    use krome_constants
    use krome_commons
    use krome_getphys
    implicit none
    real*8,intent(in)::n(nspec),Tdust,amin,amax,pexp,rho0,d2g
    real*8::k,imass(nspec),p4,p3,mass(nspec),rhod
    integer,intent(in)::idx

    !get inverse mass squared
    imass(:) = get_imass_sqrt()
    !get masses
    mass(:) = get_mass()
    !derived exponents
    p3 = pexp + 3.
    p4 = pexp + 4.

    !total dust density, g/cm3
    rhod = sum(n(1:nmols)*mass(1:nmols))*d2g

    !compute rate (1/s) coefficient assuming normalization
    k = pre_kvgas_sqrt*sqrt(n(idx_Tgas)) * imass(idx) &
        * rhod / (4./3.*rho0) * p4 / p3 &
        * (amax**p3-amin**p3) / (amax**p4-amin**p4) &
        * dust_stick(n(idx_Tgas),Tdust)

  end function krate_stick

  !********************************
  !compact version of krate_stick
  function krate_stickSi(n,idx,Tdust) result(k)
    use krome_commons
    implicit none
    integer,intent(in)::idx
    real*8,intent(in)::n(nspec),Tdust
    real*8::k,amin,amax,d2g,rho0,pexp

    !some default values OK for silicates
    amin = 5d-7 !cm
    amax = 2.5d-5 !cm
    pexp = -3.5
    rho0 = 3d0 !g/cm3
    d2g = 1d-2

    k = krate_stick(n(:),idx,Tdust,amin,amax,pexp,rho0,d2g)

  end function krate_stickSi

  !***************************
  !evaporation rate, 1/s
  function krate_evaporation(n,idx,Tdust) result(k)
    use krome_commons
    use krome_getphys
    implicit none
    integer,intent(in)::idx
    real*8,intent(in)::n(nspec),Tdust
    real*8::k,Ebind(nspec),nu0

    nu0 = 1d12 !1/s
    Ebind(:) = get_EbindBare()

    k = nu0 * exp(-Ebind(idx)/Tdust)

  end function krate_evaporation

  !***************************
  !non-thermal evaporation rate (1/s) following Hollenbach 2009,
  ! http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:0809.1642
  !Gnot is the habing flux (1.78 is Draine)
  !Av is the visual extinction
  !crflux the ionization flux of cosmic rays, 1/s
  !yield is the efficiency of the photons to desorb the given molecule
  function krate_nonthermal_evaporation(idx, Gnot, Av, crflux, yield) result(k)
    use krome_commons
    use krome_getphys
    implicit none
    integer,intent(in)::idx
    real*8,parameter::crnot=1.3d-17
    real*8,parameter::Fnot=1d8 !desorbing photons flux, 1/s
    real*8,parameter::ap2=(3d-8)**2 !sites separation squared, cm2
    real*8,intent(in)::Gnot, Av, crflux, yield
    real*8::k,f70,kevap70(nspec)

    f70 = 3.16d-19*crflux/crnot
    kevap70(:) = get_kevap70()

    k = Gnot*Fnot*ap2*yield*exp(-1.8*Av)
    k = k + f70*kevap70(idx)

  end function krate_nonthermal_evaporation

  !***************************
  function dust_ice_fraction_array(invphi,nH2O)
    use krome_constants
    use krome_commons
    implicit none
    integer::i
    real*8::dust_ice_fraction_array(ndust)
    real*8::invphi(ndust),nH2O(ndust)

    do i=1,ndust
      dust_ice_fraction_array(i) = min(nH2O(i) * invphi(i), 1d0)
    end do

  end function dust_ice_fraction_array

  !*****************************
  function get_Ebareice_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice_exp_array(:) = 0d0

  end function get_Ebareice_exp_array

  !*****************************
  function get_Ebareice23_exp_array(invTdust)
    use krome_commons
    implicit none
    real*8::get_Ebareice23_exp_array(2*nspec),invTdust(ndust)

    get_Ebareice23_exp_array(:) = 0d0

  end function get_Ebareice23_exp_array

  !************************
  !returns the binding energy for ice coated grain (K)
  function get_Ebind_ice()
    use krome_commons
    implicit none
    real*8::get_Ebind_ice(nspec)

    get_Ebind_ice(:) = 0d0

  end function get_Ebind_ice

  !************************
  !returns the binding energy for bare grain (K)
  function get_Ebind_bare()
    use krome_commons
    implicit none
    real*8::get_Ebind_bare(nspec)

    get_Ebind_bare(:) = 0d0

  end function get_Ebind_bare

  !************************
  !returns the index of the parent dust bin (0 if none)
  function get_parent_dust_bin()
    use krome_commons
    implicit none
    integer::get_parent_dust_bin(nspec)

    get_parent_dust_bin(:) = 0

  end function get_parent_dust_bin

  !*****************************
  function get_exp_table(ain,invT)
    use krome_commons
    implicit none
    integer::ia
    real*8::get_exp_table,a,invT,ain
    real*8::x1a,f1,f2

    a = ain*invT
    a = min(a, exp_table_aMax - exp_table_da)

    ia = (a-exp_table_aMin) * exp_table_multa + 1
    ia = max(ia,1)

    x1a = (ia-1)*exp_table_da

    f1 = exp_table(ia)
    f2 = exp_table(ia+1)

    get_exp_table = (a-x1a) * exp_table_multa * (f2-f1) + f1

  end function get_exp_table

end module krome_grfuncs
!This module mainly contains shielding routine and
! function to initialize radiation background (e.g. Planck).

!############### MODULE ##############
module krome_phfuncs
contains

  !****************************
  !dust shielding factor
  function shield_dust(n,Tgas,gam)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::shield_dust,n(:),Tgas,gam,eff_d2g
    real*8::sigma_d,NHtot

    eff_d2g = dust2gas_ratio
    sigma_d = 2d-21*eff_d2g*gam !Richings et al. 2014
    !sigma_d = 2d-21 !Glover+2007
    !sigma_d = 4d-22 !Richings+ 2014
    !sigma_d = 4d-21 !Gnedin 2009

    NHtot = 0d0
    NHtot  = NHtot + num2col(n(idx_H),n(:))
    NHtot  = NHtot + num2col(n(idx_Hj),n(:))
    NHtot  = NHtot + 2d0 * num2col(n(idx_H2),n(:))

    shield_dust = exp(-sigma_d*NHtot)

  end function shield_dust

  !**********************
  !planck function in eV/s/cm2/Hz/sr
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB(x,Tbb)
    use krome_constants
    implicit none
    real*8::Tbb,x,xexp,planckBB

    !exponent
    xexp = x/boltzmann_eV/Tbb

    !default value
    planckBB = 0d0

    !limit exp overflow
    if(xexp<3d2.and.x>1d-10) then
      planckBB = 2d0*x**3/planck_eV**2/clight**2 &
          / (exp(xexp)-1d0)
    end if

  end function planckBB

  !********************
  !planck function dTdust differential
  ! in eV/s/cm2/Hz/sr/K, where
  ! x is the energy in eV, Tbb the black body
  ! temperature in K
  function planckBB_dT(x,Tbb)
    use krome_constants
    real*8::a,b,x,Tbb,xexp,planckBB_dT

    b = 1d0/boltzmann_eV
    xexp = b*x/Tbb

    planckBB_dT = 0d0

    if(xexp<3d2) then
      a = 2d0/planck_eV**2/clight**2
      planckBB_dT = a*b*x**4/Tbb/Tbb * exp(xexp)/(exp(xexp)-1d0)**2
    end if

  end function planckBB_dT

  !***********************
  !shielding function selected with -shield option
  function krome_fshield(n,Tgas)
    implicit none
    real*8::krome_fshield,n(:),Tgas

    krome_fshield = 1d0 !default shielding value

  end function krome_fshield

  !**************************
  !shielding function for H2O+ and H3O+
  ! following Glover+2010 MNRAS sect 2.2 eqn.4
  function fHnOj(Av)
    implicit none
    real*8::fHnOj,Av
    if(Av.le.15d0) then
      fHnOj = exp(-2.55*Av+0.0165*Av**2)
    else
      fHnOj = exp(-2.8*Av)
    end if
  end function fHnOj

  !******************************
  !self-shielding for H2
  ! following Glover+2010 MNRAS sect2.2 eqn.6
  ! N: column density (cm-2)
  ! b: doppler broadening (cm/s)
  function fselfH2(N, b)
    implicit none
    real*8::fselfH2,N,b,x,b5

    x = N*2d-15 !normalized column density (#)
    b5 = b*1d-5 !normalized doppler broadening (#)

    fselfH2 = 0.965d0/(1+x/b5)**2 + &
        0.035d0/sqrt(1d0+x) * &
        exp(max(-8.5d-4*sqrt(1+x),-250.))

  end function fselfH2

end module krome_phfuncs

!############### MODULE ##############
module krome_subs
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2022-06-08 12:22:52
  !  Changeset d0f9dcf
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !************************
  !compute reaction rates cm^3(n-1)/s
  function coe(n)
    use krome_commons
    use krome_constants
    use krome_user_commons
    use krome_getphys
    use krome_grfuncs
    use krome_phfuncs
    use krome_fit
    implicit none
    real*8::coe(nrea),k(nrea),Tgas,n(nspec),kmax
    real*8::small,nmax
    integer::i
    real*8::T !preproc from coevar
    !Tgas is in K
    Tgas = max(n(idx_Tgas), phys_Tcmb)
    Tgas = min(Tgas,1d9)

    !maxn initialization can be removed and small can be
    ! replaced with a proper value according to the environment
    nmax = max(maxval(n(1:nmols)),1d0)
    small = 1d-40/(nmax*nmax*nmax)

    T = Tgas

    k(:) = small !inizialize coefficients

    !H+ + E -> H
    k(1) = small + (3.3d-12&
        *((Tgas&
        /300.0)**(-0.71))*exp(-2.8/Tgas))

    !H -> H+ + E
    k(2) = small + (5.6d7&
        *((user_Trad&
        /300.0)**(0.67))*exp(-1.57803d5/user_Trad))

    !H- -> H + E
    k(3) = small + (1.5d4&
        *((user_Trad&
        /300.0)**(2.29))*exp(-8763/user_Trad))

    !H + E -> H-
    k(4) = small + (2.9d-16*((Tgas&
        /300.0)**(0.78))&
        *exp(-12.0/Tgas))

    !H- + H -> H2 + E
    k(5) = small + (3.6d-9&
        *((Tgas&
        /300.0)**(-0.311))*exp(-14.7/Tgas))

    !H- + H+ -> H + H
    k(6) = small + (1.9d-7&
        *((Tgas&
        /300.0)**(-0.442)))

    !H + H+ -> H2+
    k(7) = small + (5.5d-19&
        *((Tgas&
        /300.0)**(1.80))*exp(-24.0/Tgas))

    !H2+ -> H + H+
    k(8) = small + (4d5&
        *((user_Trad&
        /300.0)**(1.21))*exp(-3.134d4/user_Trad))

    !H2+ + H -> H2 + H+
    k(9) = small + (8.9d-10)

    !H2 + H+ -> H2+ + H
    k(10) = small + (2.44d-10&
        *exp(-21356.0&
        /Tgas))

    !D+ + E -> D
    k(11) = small + (3.3d-12&
        *((Tgas&
        /300.0)**(-0.71))*exp(-2.8/Tgas))

    !D -> D+ + E
    k(12) = small + (5.6d7&
        *((user_Trad&
        /300.0)**(0.67))*exp(-1.57803d5/user_Trad))

    !D + H+ -> D+ + H
    k(13) = small + (1.96d-9&
        *((Tgas&
        /300.0)**(0.410))*exp(-36.4/Tgas))

    !D+ + H -> D + H+
    k(14) = small + (1.95d-9&
        *((Tgas&
        /300.0)**(0.410))*exp(6.5/Tgas))

    !D+ + H2 -> HD + H+
    k(15) = small + (1.66d-9&
        *exp(-5.6&
        /Tgas))

    !HD + H+ -> D+ + H2
    k(16) = small + (1.34d-9&
        *exp(-504.0&
        /Tgas))

    !D + H2 -> HD + H
    k(17) = small + (2d-12&
        *((Tgas&
        /300.0)**(2.16))*exp(-2.60d3/Tgas))

    !HD + H -> D + H2
    k(18) = small + (3.4d-12&
        *((Tgas&
        /300.0)**(1.52))*exp(-3.29d3/Tgas))

    !H- + D -> HD + E
    k(19) = small + (3.6d-9&
        *((Tgas&
        /300.0)**(-0.311))*exp(-14.7/Tgas))

    !HE + H+ -> HEH+
    k(20) = small + (7d-20&
        *((Tgas&
        /300.0)**(-0.55))*exp(-9.0/Tgas))

    !HEH+ + H -> HE + H2+
    k(21) = small + (1.06d-9&
        *((Tgas&
        /300.0)**(0.252))*exp(-24.2/Tgas))

    !HEH+ + E -> HE + H
    k(22) = small + (1.8d-9&
        *((Tgas&
        /300.0)**(-0.18)))

    !HEH+ -> HE + H+
    k(23) = small + (2d6&
        *((user_Trad&
        /300.0)**(-1.9))*exp(-2.56d4/user_Trad))

    !H + H + H -> H2 + H
    k(24) = small + (2.6d-32&
        *((Tgas&
        /300.)**(-0.35))*exp(1.6/Tgas))

    !H + H2 -> H + H + H
    k(25) = small + (5.8d-8&
        *((Tgas&
        /300.)**(-0.38))*exp(-5.2d4/Tgas))

    coe(:) = k(:) !set coefficients to return variable

    !!uncomment below to check coefficient values
    !kmax = 1d0
    !if(maxval(k)>kmax.or.minval(k)<0d0) then
    !   print *,"***************"
    !   do i=1,size(k)
    !      if(k(i)<0d0.or.k(i)>kmax) print *,i,k(i)
    !   end do
    !end if
  end function coe

  !*************************
  subroutine loadReactionsVerbatim()
    use krome_commons
    implicit none
    character*255::fname,line
    integer::ios,i,nunit

    ! Verbatim reactions filename defaults to `reactions_verbatim.dat`
    fname = "reactions_verbatim.dat"

    reactionNames(1) = "H+ + E -> H"
    reactionNames(2) = "H -> H+ + E"
    reactionNames(3) = "H- -> H + E"
    reactionNames(4) = "H + E -> H-"
    reactionNames(5) = "H- + H -> H2 + E"
    reactionNames(6) = "H- + H+ -> H + H"
    reactionNames(7) = "H + H+ -> H2+"
    reactionNames(8) = "H2+ -> H + H+"
    reactionNames(9) = "H2+ + H -> H2 + H+"
    reactionNames(10) = "H2 + H+ -> H2+ + H"
    reactionNames(11) = "D+ + E -> D"
    reactionNames(12) = "D -> D+ + E"
    reactionNames(13) = "D + H+ -> D+ + H"
    reactionNames(14) = "D+ + H -> D + H+"
    reactionNames(15) = "D+ + H2 -> HD + H+"
    reactionNames(16) = "HD + H+ -> D+ + H2"
    reactionNames(17) = "D + H2 -> HD + H"
    reactionNames(18) = "HD + H -> D + H2"
    reactionNames(19) = "H- + D -> HD + E"
    reactionNames(20) = "HE + H+ -> HEH+"
    reactionNames(21) = "HEH+ + H -> HE + H2+"
    reactionNames(22) = "HEH+ + E -> HE + H"
    reactionNames(23) = "HEH+ -> HE + H+"
    reactionNames(24) = "H + H + H -> H2 + H"
    reactionNames(25) = "H + H2 -> H + H + H"
    return

    !verbatim reactions are loaded from file
    ! to increase compilation speed
    open(newunit=nunit,file=trim(fname),status="old",iostat=ios)
    if(ios/=0) then
      print *,"ERROR: "//trim(fname)//" file not present!"
      stop
    end if

    !load reactions from file
    do i=1,nrea
      read(nunit,'(a)',iostat=ios) line
      if(ios/=0) then
        print *,"ERROR: problem reading "//trim(fname)
        stop
      end if
      reactionNames(i) = trim(line)
    end do
    close(nunit)

  end subroutine loadReactionsVerbatim

  !*******************
  !The following functions compute the recombination rate
  ! on dust for H+, He+, C+, Si+, and O+. See Weingartner&Draine 2001
  ! dust2gas_ratio, D/D_sol, default is assumed equal to Z/Z_sol
  function H_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::H_recombination_on_dust

    H_recombination_on_dust = 0d0

    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    H_recombination_on_dust =  1.225d-13*dust2gas_ratio &
        /(1.d0+8.074d-6*psi**(1.378)*(1.d0+5.087d2 &
        *Tgas**(0.01586)*psi**(-0.4723-1.102d-5*log(Tgas))))

  end function H_recombination_on_dust

  !******************
  function He_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::He_recombination_on_dust

    He_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    He_recombination_on_dust = 5.572d-14*dust2gas_ratio&
        /(1.d0+3.185d-7*psi**(1.512)*(1.d0+5.115d3&
        *Tgas**(3.903d-7)*psi**(-0.4956-5.494d-7*log(Tgas))))

  end function He_recombination_on_dust

  !*******************
  function C_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::C_recombination_on_dust

    C_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    C_recombination_on_dust = 4.558d-13*dust2gas_ratio&
        /(1.d0+6.089d-3*psi**(1.128)*(1.d0+4.331d2&
        *Tgas**(0.04845)*psi**(-0.8120-1.333d-4*log(Tgas))))

  end function C_recombination_on_dust

  !******************
  function Si_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,psi
    real*8::Si_recombination_on_dust

    Si_recombination_on_dust = 0d0
    if(n(idx_E)<1d-20.or.GHabing<=0.d0) return

    psi = GHabing*sqrt(Tgas)/n(idx_E)

    if(psi<=0) return

    Si_recombination_on_dust = 2.166d-14*dust2gas_ratio&
        /(1.d0+5.678d-8*psi**(1.874)*(1.d0+4.375d4&
        *Tgas**(1.635d-6)*psi**(-0.8964-7.538d-5*log(Tgas))))

  end function Si_recombination_on_dust

  !********************
  function O_recombination_on_dust(n,Tgas)
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,k_H
    real*8::O_recombination_on_dust

    k_H = H_recombination_on_dust(n(:),Tgas)
    O_recombination_on_dust = 0.25d0*k_H

  end function O_recombination_on_dust

  !*********************
  !This function returns the
  ! photorate of H2 occurring in the
  ! Lyman-Werner bands following the approximation
  ! provided by Glover&Jappsen 2007. Rate in 1/s.
  !Approximation valid at low-density, it assumes H2(nu = 0).
  !It also stores the rate as a common, needed for the photoheating
  function H2_solomonLW(myflux)
    use krome_commons
    use krome_constants
    implicit none
    real*8::H2_solomonLW,myflux

    !myflux is the radiation background at E = 12.87 eV
    !should be converted to erg
    H2_solomonLW = 1.38d9*myflux*eV_to_erg

  end function H2_solomonLW

  !****************************
  !tanh smoothing function that
  ! increses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_increase(xarg,xpos,slope)
    implicit none
    real*8::smooth_increase,xarg,xpos,slope

    smooth_increase = .5d0 * (tanh(slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_increase

  !****************************
  !tanh smoothing function that
  ! decreses when xarg increases.
  ! xpos is the position of the transition point.
  ! slope is the steepness of the curve.
  function smooth_decrease(xarg,xpos,slope)
    implicit none
    real*8::smooth_decrease,xarg,xpos,slope

    smooth_decrease = .5d0 * (tanh(-slope * (xarg - xpos)) &
        + 1d0)

  end function smooth_decrease

  !*********************
  !sign: return 1d0 if x>=0d0,
  ! else return -1d0
  function get_sgn(x)
    implicit none
    real*8::x,get_sgn

    get_sgn = 1d0
    if(x==0d0) return
    get_sgn = x/abs(x)

  end function get_sgn

  !*********************
  function conserve(n,ni)
    use krome_commons
    implicit none
    real*8::conserve(nspec),n(nspec),ni(nspec),no(nspec)
    real*8::ntot,nitot,factor

    no(:) = n(:)

    conserve(:) = 0d0
    conserve(:) = no(:)

  end function conserve

  !*************************
  !this subroutine changes the x(:) mass fractions of the species
  ! to force conservation according to the reference ref(:)
  subroutine conserveLin_x(x,ref)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::x(nmols),ref(natoms)
    real*8::A(natoms,natoms),B(natoms),m(nspec)

    m(:) = get_mass()
    A(:,:) = 0d0
    B(:) = ref(:)

    !charge conservation
    x(idx_E) = m(idx_E)*(- 1d0*x(idx_Hk) / m(idx_Hk) &
        + 1d0*x(idx_Hj) / m(idx_Hj) &
        + 1d0*x(idx_H2j) / m(idx_H2j) &
        + 1d0*x(idx_Dj) / m(idx_Dj) &
        + 1d0*x(idx_HEHj) / m(idx_HEHj))
    !check if charge conservation goes wrong
    if(x(idx_E)<0d0) then
      print *,"ERROR in conserveLin, electrons < 0"
      stop
    end if

  end subroutine conserveLin_x

  !***************************
  !compute the total reference mass atom type by atom type
  function conserveLinGetRef_x(x)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::conserveLinGetRef_x(natoms),x(nmols)
    real*8::m(nspec)

    m(:) = get_mass()
    conserveLinGetRef_x(:) = 0d0

  end function conserveLinGetRef_x

  !***************************
  !Ref: Sasaki & Takahara (1993)
  !This function evaluate the recombination rate
  ! for H+ + e --> H + gamma and the same
  ! for D+ + e --> D + gamma
  function elec_recomb_ST93(nabund,nelec,ntot,nucleiH,Trad)
    use krome_commons
    use krome_constants
    implicit none
    real*8::nabund,nelec,Trad
    real*8::nucleiH,elec_recomb_ST93
    real*8::al,ak,rc2,r2c
    real*8::a0,b0,c0,d0,e0
    real*8::a1,b1,c1,d1,e1,f1,g1,h1
    real*8::ntot,ratio

    al = 8.227d0
    ak = 22.06d0 / (hubble  *(1d0 + phys_zredshift) &
        * sqrt(1d0 + Omega0 * phys_zredshift))
    !Rc2 evaluation
    rc2 = 8.76d-11 * (1d0 + phys_zredshift)**(-0.58)
    !R2c evaluation
    r2c = (1.80d10 * Trad)**(1.5) &
        * exp(-3.9472d4 / Trad) * rc2

    !coefficients
    a0 = nucleiH * rc2
    b0 = ak * al * nucleiH
    c0 = ak * rc2 * nucleiH * nucleiH
    d0 = r2c * exp(-1.18416d5/Trad)
    e0 = ak * r2c * nucleiH

    !polynomial terms
    a1 = -d0 * (1d0 + b0)
    b1 = d0 * (1d0 + 2d0 * b0)
    c1 = a0 + b0 * (a0 - d0)
    d1 = -a0 * b0
    e1 = a0 * c0
    f1 = 1d0 + b0 + e0
    g1 = -(b0 + e0)
    h1 = c0

    ratio = nabund / ntot

    elec_recomb_ST93 = ntot*(a1 + b1*ratio + c1*ratio**2 + d1*ratio**3 &
        + e1*ratio**4) / (f1 + g1*ratio + h1*ratio**2)

    elec_recomb_ST93 = elec_recomb_ST93 / (nabund * nelec)

  end function elec_recomb_ST93

  !********************
  subroutine load_parts()
    use krome_commons
    implicit none

  end subroutine load_parts

  !*************************
  subroutine load_part(fname,array_part,min_part,dT_part)
    character(len=*)::fname
    integer::ios,icount,i,cv
    real*8,allocatable::array_part(:),emed(:)
    real*8::min_part,dT_part,Told,array_tmp(int(1e5)),rout(2)

    open(33,file=trim(fname),status="old",iostat=ios)
    if(ios.ne.0) then
      print *,"ERROR: partition function not found"
      print *," in file "//fname
      stop
    end if

    print *,"loading partition function from "//fname
    icount = 0
    min_part = 1d99
    Told = 0d0
    do
      read(33,*,iostat=ios) rout(:)
      if(ios<0) exit
      if(ios.ne.0) cycle
      icount = icount + 1
      min_part = min(min_part,rout(1))
      array_tmp(icount) = rout(2)
      dT_part = rout(1) - Told
      Told = rout(1)
    end do
    close(33)

    allocate(array_part(icount),emed(icount))
    array_part(:) = array_tmp(1:icount)

  end subroutine load_part

  !**********************
  function troe_falloff(k0,kinf,Fc,m)
    implicit none
    real*8::troe_falloff,k0,kinf,Fc,m,rm,xexp
    rm = k0*m/kinf
    xexp = 1d0/(1d0+log10(rm)**2)
    troe_falloff = k0*m/(1d0+rm)*Fc**xexp
  end function troe_falloff

  !*************************
  function k3body(k0,kinf,Fc,nM)
    implicit none
    real*8::k3body,k0,kinf,Fc,nM
    real*8::c,n,d,Pr,xexp,F

    c = -0.4d0-0.67d0*log10(Fc)
    n = 0.75d0-1.27d0*log10(Fc)
    d = 0.14d0
    Pr = k0*nM/kinf
    xexp = (log10(Pr)+c)/(n-d*(log10(Pr)+c))
    F = 1d1**(log10(Fc)/(1d0+xexp**2))
    k3body = kinf*(Pr/(1d0+Pr)) * F

  end function k3body

  !***********************
  !see http://kida.obs.u-bordeaux1.fr/help
  function KIDA3body(ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,&
        kcFc,kdFc,npart,Tgas,pmin,pmax)
    implicit none
    real*8::ka0,kb0,kc0,kaInf,kbInf,kcInf,kaFc,kbFc,kcFc,kdFc
    real*8::KIDA3body,kinf,p,f,npart,Tgas,fc,fexp,invT
    real*8::k0,cc,dd,nn,pmin,pmax

    KIDA3body = 0d0

    invT = 1d0/Tgas
    k0 = ka0*(Tgas/3d2)**kb0*exp(-kc0*invT)
    kinf = kainf*(Tgas/3d2)**kbinf*exp(-kcinf*invT)

    p = k0*npart/kinf
    if(p<pmin) return
    if(p>pmax) return

    fc = (1d0-kaFc)*exp(-Tgas/kbFc) + kaFc*exp(-Tgas/kbFc) &
        + exp(-kdFc*invT)

    cc = -0.4d0 - 0.67d0 *log10(fc)
    dd = 0.14d0
    nn = 0.75d0 - 1.27d0*log10(fc)
    fexp = 1d0 + ((log10(p)+cc)/(nn-dd*(log10(p)+cc)))**2

    f = fc**(1d0/fexp)

    KIDA3body = kinf*(p/(1d0+p))*f

  end function KIDA3body

  !******************************
  !collisional ionization rate from Verner+96
  ! unit: cm3/s
  function colion_v96(Tgas,dE,P,A,X,K)
    implicit none
    real*8::colion_v96,Tgas,dE,A,X,K,U,Te,P

    Te = Tgas * 8.621738d-5 !K to eV
    U = dE / Te
    colion_v96 = A * (1d0 + P*sqrt(U)) * U**K * exp(-U) / (X+U)

  end function colion_v96

  !****************************
  !radiative recombination rates from
  ! Verner routine, standard fit, cm3/s
  function recV96(Tgas,a,b)
    implicit none
    real*8::recV96,Tgas,a,b

    recV96 = a*(1d4/Tgas)**b

  end function recV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, new fit, cm3/s
  function recNewV96(Tgas,r1,r2,r3,r4)
    implicit none
    real*8::recNewV96,Tgas,r1,r2,r3,r4,tt

    tt = sqrt(Tgas/r3)
    recNewV96 = r1/(tt*(tt + 1d0)**(1.-r2) &
        * (1d0 + sqrt(Tgas/r4))**(1.+r2))

  end function recNewV96

  !****************************
  !radiative recombination rates from
  ! Verner routine, iron only, cm3/s
  function recFeV96(Tgas,r1,r2,r3)
    implicit none
    real*8::recFeV96,Tgas,r1,r2,r3,tt

    tt = sqrt(Tgas*1d-4)
    recFeV96 = r1/tt**(r2 + r3 + log10(tt))

  end function recFeV96

  !******************************
  !radiative recombination rates from Verner+96
  ! unit: cm3/s
  function radrec_v96(Tgas,a,b,T0,T1)
    implicit none
    real*8::Tgas,a,b,T0,T1,radrec_v96,iT0

    iT0 = 1d0/T0
    radrec_v96 = a/(sqrt(Tgas*iT0) + (1d0*sqrt(Tgas*iT0))**(1.-b) &
        * (1d0+sqrt(Tgas/T1))**(1+b))

  end function radrec_v96

  !*******************************
  !radiative recombination rates low-temp fit, Verner+96
  ! unit: cm3/s
  function radrec_low_v96(Tgas,a,b,c,d,f)
    implicit none
    real*8::Tgas,a,b,c,d,f,radrec_low_v96,t,invt

    t = Tgas*1d-4
    invt = 1d0/t

    radrec_low_v96 = 1d-12 * (a*invt + b + c*t + d*t**2) &
        * t**(-1.5) * exp(-f*invt)

    radrec_low_v96 = max(0d0,radrec_low_v96)

  end function radrec_low_v96

  !***************************
  !Collisional dissociation rate (cm-3/s) by Martin et al. 1996
  ! H2+H->H+H+H
  !NOTE: the use of this rate is suggested
  ! for high-density regime and in the presence of UV backgrounds.
  ! if necessary it must be included in the reaction file as
  ! H2,H,,H,H,H,,NONE,NONE,dissH2_Martin96(n,Tgas)
  function dissH2_Martin96(n,Tgas)
    use krome_commons
    use krome_getphys
    integer::i
    real*8::n(nspec),Tgas,dissH2_Martin96
    real*8::CDrates,logTv(4),k_CIDm(21,2),k_CID,invT,logT,n_c1,n_c2,n_H
    real*8::logk_h1,logk_h2,logk_l1,logk_l2,logn_c1,logn_c2,p,logk_CID
    real*8::logT2,logT3

    !k_CID = collision-induced dissociation + dissociative tunneling

    !Collisional dissociation of H2
    k_CIDm(:,1) = (/-178.4239d0, -68.42243d0, 43.20243d0, -4.633167d0, &
        69.70086d0, 40870.38d0, -23705.70d0, 128.8953d0, -53.91334d0, &
        5.315517d0, -19.73427d0, 16780.95d0, -25786.11d0, 14.82123d0, &
        -4.890915d0, 0.4749030d0, -133.8283d0, -1.164408d0, 0.8227443d0,&
        0.5864073d0, -2.056313d0/)

    !Dissociative tunneling of H2
    k_CIDm(:,2) = (/-142.7664d0, 42.70741d0, -2.027365d0, -0.2582097d0, &
        21.36094d0, 27535.31d0, -21467.79d0, 60.34928d0, -27.43096d0, &
        2.676150d0, -11.28215d0, 14254.55d0, -23125.20d0, 9.305564d0, &
        -2.464009d0, 0.1985955d0, 743.0600d0, -1.174242d0, 0.7502286d0, &
        0.2358848d0, 2.937507d0/)

    n_H  = get_Hnuclei(n(:))
    logT = log10(Tgas)
    invT = 1.0d0/Tgas
    logT2 = logT*logT
    logT3 = logT2*logT
    logTv = (/1.d0, logT, logT2, logT3/)
    k_CID = 0.d0
    do i=1,2
      logk_h1 = k_CIDm(1,i)*logTv(1) + k_CIDm(2,i)*logTv(2) + &
          k_CIDm(3,i)*logTv(3) + k_CIDm(4,i)*logTv(4) + &
          k_CIDm(5,i)*log10(1.d0+k_CIDm(6,i)*invT)
      logk_h2 = k_CIDm(7,i)*invT
      logk_l1 = k_CIDm(8,i)*logTv(1) + k_CIDm(9,i)*logTv(2) + &
          k_CIDm(10,i)*logTv(3) + k_CIDm(11,i)*log10(1.d0+k_CIDm(12,i)*invT)
      logk_l2 = k_CIDm(13,i)*invT
      logn_c1 = k_CIDm(14,i)*logTv(1) + k_CIDm(15,i)*logTv(2) &
          + k_CIDm(16,i)*logTv(3) + k_CIDm(17,i)*invT
      logn_c2 = k_CIDm(18,i) + logn_c1
      p = k_CIDm(19,i) + k_CIDm(20,i)*exp(-Tgas/1.850d3) &
          + k_CIDm(21,i)*exp(-Tgas/4.40d2)
      n_c1 = 1d1**(logn_c1)
      n_c2 = 1d1**(logn_c2)
      logk_CID = logk_h1 - (logk_h1 - logk_l1) / (1.d0 + (n_H/n_c1)**p) &
          + logk_h2 - (logk_h2 - logk_l2) / (1.d0 + (n_H/n_c2)**p)
      k_CID = k_CID + 1.d1**logk_CID
    enddo

    dissH2_Martin96 = k_CID

  end function dissH2_Martin96

  !***********************************
  subroutine init_exp_table()
    use krome_commons
    implicit none
    integer::i
    real*8::a

    do i=1,exp_table_na
      a = (i-1)*(exp_table_aMax-exp_table_aMin)/(exp_table_na-1) + exp_table_aMin
      exp_table(i) = exp(-a)
    end do

  end subroutine init_exp_table

  !***************************
  !get the index of the specie name
  function get_index(name)
    use krome_commons
    use krome_getphys
    integer::get_index,i
    character*16::names(nspec)
    character*(*)::name
    names(:) = get_names()
    get_index = -1 !default index
    !loop on species to found the specie named name
    do i=1,nspec
      !when found store and break loop
      if(trim(names(i))== trim(name)) then
        get_index = i !store index
        exit
      end if
    end do

    !error if species not found
    if(get_index<0) then
      print *,"ERROR: can't find the index of ",name
      stop
    end if

  end function get_index

  !*****************************
  !computes revers kinetics from reaction and
  ! product indexes
  ! k_rev = k_for * revKc
  ! Note that reaction constant revKc is calculated with
  ! reactants and products from reverse reaction
  function revKc(Tgas,ridx,pidx)
    use krome_constants
    use krome_commons
    implicit none
    real*8::revKc,Tgas,dgibss,stoichiometricChange
    integer::ridx(:),pidx(:),i

    ! when considering forward reaction:
    ! Kc = (P°)**(p+p-r-r) * exp(-dGibss_forward°)
    ! where ° means at standard conditions of
    ! P° = 1 bar = (kb*T/1e6) dyn/cm^2 (cgs)
    ! when considering reverse:
    ! 1/Kc = revKc = (kb*T/1e6)**(p+p-r-r) * exp(-dGibss_reverse°)
    ! kb*T/1e6 is to go from 1 atm pressure to number density cm^-3
    ! When not at standard pressure this does not change:
    ! revKc = P**(p+p-r-r) *exp(-dGibss_reverse° - (p+p-r-r)*ln(P/P°))
    !       = (P°)**(p+p-r-r) * exp(-dGibss_reverse°)

    dgibss = 0.d0 ! Gibbs free energy/(R*T)
    stoichiometricChange = 0d0

    do i=1,size(pidx)
      dgibss = dgibss + revHS(Tgas,pidx(i))
      stoichiometricChange = stoichiometricChange + 1
    end do

    do i=1,size(ridx)
      dgibss = dgibss - revHS(Tgas,ridx(i))
      stoichiometricChange = stoichiometricChange - 1
    end do

    revKc = (boltzmann_erg * Tgas * 1e-6)**(-stoichiometricChange)&
        * exp(-dgibss)

  end function revKc

  !*****************************
  !compute H-S for species with index idx
  ! when temperature is Tgas
  function revHS(Tgas,idx)
    use krome_commons
    use krome_constants
    use krome_fit
    real*8::revHS,Tgas,Tgas2,Tgas3,Tgas4,invT,lnT,H,S
    real*8::Tnist,Tnist2,Tnist3,Tnist4,invTnist,invTnist2,lnTnist
    real*8::p1_nasa(15,7), p2_nasa(15,7), Tlim_nasa(15,3), p(7)
    real*8::p1_nist(15,7), p2_nist(15,7), Tlim_nist(15,3)
    integer::idx

    p(:) = 0.d0
    p1_nasa(:,:) = 0.d0
    p2_nasa(:,:) = 0.d0
    Tlim_nasa(:,:) = 0.d0
    p1_nist(:,:) = 0.d0
    p2_nist(:,:) = 0.d0
    Tlim_nist(:,:) = 0.d0
    Tgas2 = Tgas * Tgas
    Tgas3 = Tgas2 * Tgas
    Tgas4 = Tgas3 * Tgas
    invT = 1d0/Tgas
    lnT = log(Tgas)
    ! NIST polynomials are quite differernt
    ! it doesn't like easy stuff...
    Tnist = Tgas * 1.d-3
    Tnist2 = Tnist * Tnist
    Tnist3 = Tnist2 * Tnist
    Tnist4 = Tnist3 * Tnist2
    invTnist = 1d0/Tnist
    invTnist2 = invTnist * invTnist
    lnTnist = log(Tnist)

    p1_nasa(idx_E,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        -11.7208122d0/)
    p1_nasa(idx_Hk,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        15976.167d0,&
        -1.1390139d0/)
    p1_nasa(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p1_nasa(idx_H2,:)  = (/2.34433112d0,&
        0.00798052075d0,&
        -1.9478151d-05,&
        2.01572094d-08,&
        -7.37611761d-12,&
        -917.935173d0,&
        0.683010238d0/)
    p1_nasa(idx_D,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25921.2596d0,&
        0.591714338d0/)
    p1_nasa(idx_HD,:)  = (/3.43752369d0,&
        0.000617471555d0,&
        -1.85267846d-06,&
        2.32581486d-09,&
        -8.35140695d-13,&
        -17.7564616d0,&
        -2.41112115d0/)
    p1_nasa(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p1_nasa(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p1_nasa(idx_H2j,:)  = (/3.77256072d0,&
        -0.0019574659d0,&
        4.54812047d-06,&
        -2.82152141d-09,&
        5.33969209d-13,&
        178694.654d0,&
        -3.96609192d0/)
    p1_nasa(idx_Dj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184512.004d0,&
        -0.101841452d0/)
    p1_nasa(idx_HEHj,:)  = (/3.46125597d0,&
        0.000458215398d0,&
        -1.84588136d-06,&
        2.85619143d-09,&
        -1.19055292d-12,&
        162267.791d0,&
        -1.4734683d0/)
    p2_nasa(idx_E,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        -11.7208122d0/)
    p2_nasa(idx_Hk,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        15976.167d0,&
        -1.1390139d0/)
    p2_nasa(idx_H,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25473.66d0,&
        -0.44668285d0/)
    p2_nasa(idx_H2,:)  = (/2.93286575d0,&
        0.000826608026d0,&
        -1.46402364d-07,&
        1.54100414d-11,&
        -6.888048d-16,&
        -813.065581d0,&
        -1.02432865d0/)
    p2_nasa(idx_D,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        25921.2596d0,&
        0.591714338d0/)
    p2_nasa(idx_HD,:)  = (/2.80029834d0,&
        0.0011562336d0,&
        -3.06064442d-07,&
        4.51518392d-11,&
        -2.62838877d-15,&
        238.213151d0,&
        1.23069947d0/)
    p2_nasa(idx_HE,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        -745.375d0,&
        0.928723974d0/)
    p2_nasa(idx_Hj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184021.488d0,&
        -1.14064664d0/)
    p2_nasa(idx_H2j,:)  = (/3.44204765d0,&
        0.000599083239d0,&
        6.69133685d-08,&
        -3.43574373d-11,&
        1.97626599d-15,&
        178650.236d0,&
        -2.79499055d0/)
    p2_nasa(idx_Dj,:)  = (/2.5d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        0.0d0,&
        184512.004d0,&
        -0.101841452d0/)
    p2_nasa(idx_HEHj,:)  = (/2.7886533d0,&
        0.00131211424d0,&
        -4.20210058d-07,&
        6.21496385d-11,&
        -3.47860467d-15,&
        162499.316d0,&
        2.24055051d0/)
    Tlim_nasa(idx_E,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Hk,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_D,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HD,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HE,:)  = (/200.0d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Hj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_H2j,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_Dj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)
    Tlim_nasa(idx_HEHj,:)  = (/298.15d0,&
        1000.0d0,&
        6000.0d0/)

    ! pick NASA data if present for species
    if (Tlim_nasa(idx,2) /= 0.d0) then
      !select set of NASA polynomials using temperature
      if(Tlim_nasa(idx,1).le.Tgas .and. Tgas.le.Tlim_nasa(idx,2)) then
        p(:) = p1_nasa(idx,:)

      else if(Tlim_nasa(idx,2)<Tgas .and. Tgas.le.Tlim_nasa(idx,3)) then
        p(:) = p2_nasa(idx,:)

        ! currently no option when Tgas not in Tlim range p(:) = 0
      end if

      !compute NASA polynomials for enthalpy and enthropy (unitless)
      H = p(1) + p(2)*0.5d0*Tgas + p(3)*Tgas2/3.d0 + p(4)*Tgas3*0.25d0 + &
          p(5)*Tgas4*0.2d0 + p(6)*invT
      S = p(1)*lnT + p(2)*Tgas + p(3)*Tgas2*0.5d0 + p(4)*Tgas3/3.d0 + &
          p(5)*Tgas4*0.25d0 + p(7)

      revHS = H - S

      ! else pick NIST data (if present)
    else if (Tlim_nist(idx,2) /= 0.d0) then
      if (Tlim_nist(idx,1) < Tgas .and. Tgas < Tlim_nist(idx,2)) then
        p(:) = p1_nist(idx,:)

      else if (Tlim_nist(idx,2) < Tgas .and. Tgas < Tlim_nist(idx,3)) then
        p(:) = p2_nist(idx,:)

        ! currently no option when Tgas not in Tlim range p(:) = 0
      end if

      !compute NIST polynomials for enthalpy and enthropy
      ! H in (kJ/mol)
      H = p(1)*Tnist + p(2)*0.5d0*Tnist2 + p(3)*Tnist3/3.d0 + p(4)*Tnist4*0.25d0&
          - p(5)*invTnist + p(6)
      !  Unitsless
      H = H / (Rgas_kJ * Tgas)

      ! S in (J/mol*K)
      S = p(1)*lnTnist + p(2)*Tnist + p(3)*Tnist2*0.5d0 + p(4)*Tnist3/3.d0&
          - p(5)*invTnist2*0.5d0 + p(7)
      !  Unitless. Note: do not use Tnist
      S = S / Rgas_J

      revHS = H - S

      ! return zero is no data exists
    else
      print *, "No thermochemical data of species index", idx
      revHS = 0.d0

    end if

  end function revHS

  !******************************
  subroutine print_best_flux(n,Tgas,nbestin)
    !print the first nbestin fluxes
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea)
    integer::nbest,idx(nrea),i,nbestin
    character*50::name(nrea)

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux

  !******************************
  subroutine print_best_flux_frac(n,Tgas,frac)
    !print the first nbestin fluxes
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),frac
    integer::idx(nrea),i
    character*50::name(nrea)

    if(frac>1d0) then
      print *,"ERROR: fraction in krome_print_best_flux should be <=1!"
      stop
    end if

    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nrea
      if(flux(idx(i))<flux(idx(1))*frac) exit
      print '(I8,a1,a50,E17.8)',idx(i)," ",name(idx(i)),flux(idx(i))
    end do

  end subroutine print_best_flux_frac

  !******************************
  subroutine print_best_flux_spec(n,Tgas,nbestin,idx_found)
    !print the first nbestin fluxes for the reactions
    ! that contains the species with index idx_found
    use krome_commons
    use krome_getphys
    implicit none
    real*8::n(nspec),Tgas,flux(nrea),maxflux
    integer::nbest,idx(nrea),i,nbestin,idx_found
    character*50::name(nrea)
    logical::found

    nbest = min(nbestin,nrea) !cannot exceed the number of reactions
    maxflux = 0d0
    flux(:) = get_flux(n(:),Tgas) !get fluxes
    name(:) = get_rnames() !get reaction names
    do i=1,nrea
      found = .false.
      if(arr_r1(i) == idx_found) found = .true.
      if(arr_r2(i) == idx_found) found = .true.
      if(arr_r3(i) == idx_found) found = .true.
      if(arr_p1(i) == idx_found) found = .true.
      if(arr_p2(i) == idx_found) found = .true.
      if(arr_p3(i) == idx_found) found = .true.
      maxflux = max(maxflux,flux(i))
      if(.not.found) flux(i) = 0d0
    end do

    !call the sorting algorithm (bubblesort)
    idx(:) = idx_sort(flux(:))

    !print to screen
    print *,"***************"
    do i=1,nbest
      print '(I8,a1,a50,2E17.8)',idx(i)," ",name(idx(i)),flux(idx(i)),&
          flux(idx(i))/maxflux
    end do

  end subroutine print_best_flux_spec

  !*****************************
  function idx_sort(fin)
    !sorting algorithm: requires an array of real values fin
    ! and returns the sorted index list. descending.
    ! bubblesort: not very efficient, replace with what you prefer
    implicit none
    real*8::fin(:),f(size(fin)),ftmp
    integer::idx_sort(size(fin)),n,itmp,i
    logical::found

    f(:) = fin(:) !copy to local

    n = size(f)
    !init indexes
    do i=1,n
      idx_sort(i) = i
    end do

    !loop to sort
    do
      found = .false. !swapped something flag
      do i=2,n
        !> for descending, < for ascending
        if(f(i)>f(i-1)) then
          found = .true.
          !swap real value
          ftmp = f(i)
          f(i) = f(i-1)
          f(i-1) = ftmp
          !swap index
          itmp = idx_sort(i)
          idx_sort(i) = idx_sort(i-1)
          idx_sort(i-1) = itmp
        end if
      end do
      !if nothing swapped exit
      if(.not.found) exit
    end do

  end function idx_sort

  !******************************
  function get_flux(n,Tgas)
    !get the flux k*n*n*... of the rates
    use krome_commons
    implicit none
    integer::i
    integer::r1,r2,r3
    real*8::get_flux(nrea),n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      r3 = arr_r3(i)
      arr_flux(i) = k(i)*n(r1)*n(r2)*n(r3)
    end do
    get_flux(:) = arr_flux(:)

  end function get_flux

  !*****************************
  subroutine load_arrays()
    !load the array containing reactants
    ! and product index
    use krome_commons

    arr_r1(1:25) = (/8,3,2,3,2,2,3,9,9,4,10,5,5,10,10,6,5,6,2,7&
        ,11,11,11,3,3/)
    arr_r2(1:25) = (/1,15,15,1,3,8,8,15,3,8,1&
        ,15,8,3,4,8,4,3,5,8,3,1,15,3,4/)
    arr_r3(1:25) = (/15,15,15&
        ,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15&
        ,3,15/)
    arr_p1(1:25) = (/3,8,3,2,4,3,9,3,4,9,5,10,10,5,6&
        ,10,6,5,6,11,7,7,7,4,3/)
    arr_p2(1:25) = (/15,1,1,15,1,3,15&
        ,8,8,3,15,1,3,8,8,4,3,4,1,15,9,3,8,3&
        ,3/)
    arr_p3(1:25) = (/15,15,15,15,15,15,15,15,15,15,15,15&
        ,15,15,15,15,15,15,15,15,15,15,15,15,3/)

  end subroutine load_arrays

  ! ************************************
  ! solves linear least squares
  subroutine llsq(n, x, y, a, b)

    !****************************************************
    !
    !! LLSQ solves a linear least squares problem matching a line to data.
    !
    !  Discussion:
    !
    !    A formula for a line of the form Y = A * X + B is sought, which
    !    will minimize the root-mean-square error to N data points
    !    ( X(I), Y(I) );
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 March 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    In: N, the number of data values.
    !
    !    In: X(N), Y(N), the coordinates of the data points.
    !
    !    Out: A, B, the slope and Y-intercept of the
    !    least-squares approximant to the data.
    !
    implicit none
    integer,intent(in)::n
    real*8,intent(out)::a, b
    real*8,intent(in)::x(n), y(n)
    real*8::bot, top, xbar, ybar

    ! special case
    if(n == 1) then
      a = 0d0
      b = y(1)
      return
    end if

    ! average X and Y
    xbar = sum(x) / n
    ybar = sum(y) / n

    ! compute beta
    top = dot_product(x(:) - xbar, y(:) - ybar)
    bot = dot_product(x(:) - xbar, x(:) - xbar)

    ! if top is zero a is zero
    if(top==0d0) then
      a = 0d0
    else
      a = top / bot
    end if

    b = ybar - a * xbar

  end subroutine llsq

end module krome_subs

!############### MODULE ##############
module krome_stars

end module krome_stars

!############### MODULE ##############
module krome_dust
contains

  !***********************
  subroutine init_dust_tabs()
    use krome_commons
    use krome_fit
    implicit none

    call init_anytab2D("dust_table_cool.dat",dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_cool(:,:), dust_mult_ngas, &
        dust_mult_Tgas)
    call init_anytab2D("dust_table_Tdust.dat",dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_Tdust(:,:), dust_mult_ngas, &
        dust_mult_Tgas)
    call init_anytab2D("dust_table_H2.dat",dust_tab_ngas(:), &
        dust_tab_Tgas(:), dust_tab_H2(:,:), dust_mult_ngas, &
        dust_mult_Tgas)

  end subroutine init_dust_tabs

end module krome_dust

!############### MODULE ##############
module krome_photo
contains

end module krome_photo

!############### MODULE ##############
module krome_tabs
contains

  !***********************+
  function coe_tab(n)
    !interface to tabs
    use krome_subs
    use krome_getphys
    use krome_phfuncs
    use krome_grfuncs
    use krome_constants
    use krome_commons
    use krome_user_commons
    implicit none
    integer::idx,j
    real*8::Tgas, coe_tab(nrea),n(nspec),small

    Tgas = max(n(idx_Tgas),phys_Tcmb)
    small = 0d0

    coe_tab(:) = coe(n(:))

  end function coe_tab

end module krome_tabs

!############### MODULE ##############
module KROME_coolingGH
end module KROME_coolingGH


!############### MODULE ##############
module KROME_cooling
  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2022-06-08 12:22:52
  !  Changeset d0f9dcf
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************
  integer,parameter::coolTab_n=int(1e2)
  integer,parameter::nZrate=0
  real*8::coolTab(nZrate,coolTab_n),coolTab_logTlow, coolTab_logTup
  real*8::coolTab_T(coolTab_n),inv_coolTab_T(coolTab_n-1),inv_coolTab_idx
contains

  !*******************
  function cooling(n,inTgas)
    use krome_commons
    implicit none
    real*8::n(:),inTgas,cooling,Tgas

    Tgas = inTgas
    cooling = sum(get_cooling_array(n(:),Tgas))

  end function cooling

  !*******************************
  function get_cooling_array(n, Tgas)
    use krome_commons
    implicit none
    real*8::n(:), Tgas
    real*8::get_cooling_array(ncools),cools(ncools)
    real*8::f1,f2,smooth

    f1 = 1d0
    f2 = 1d0

    !returns cooling in erg/cm3/s
    cools(:) = 0d0

    cools(idx_cool_H2) = cooling_H2(n(:), Tgas)

    cools(idx_cool_HD) = cooling_HD(n(:), Tgas) - cooling_HD(n(:), phys_Tcmb)

    cools(idx_cool_atomic) = cooling_Atomic(n(:), Tgas)

    cools(idx_cool_compton) = cooling_compton(n(:), Tgas)

    cools(idx_cool_cont) = cooling_continuum(n(:), Tgas)

    cools(idx_cool_custom) = cooling_custom(n(:),Tgas)

    get_cooling_array(:) = cools(:)

  end function get_cooling_array

  !*****************************
  function cooling_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,cooling_custom

    cooling_custom = 0d0

  end function cooling_custom

  !**********************************
  function kpla(n,Tgas)
    !Planck opacity mean fit (Lenzuni+1996)
    !only denisity dependent (note that the
    ! fit provided by Lenzuni is wrong)
    ! valid for T<3e3 K
    !use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::kpla,rhogas,Tgas,n(:),y
    real*8::a0,a1,m(nspec)

    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3

    kpla = 0.d0
    !opacity is zero under 1e-12 g/cm3
    if(rhogas<1d-12) return

    !fit coefficients
    a0 = 1.000042d0
    a1 = 2.14989d0

    !log density cannot exceed 0.5 g/cm3
    y = log10(min(rhogas,0.5d0))

    kpla = 1d1**(a0*y + a1) !fit density only

  end function kpla

  !*****************************
  function coolingChem(n,Tgas)
    implicit none
    real*8::coolingChem,n(:),Tgas

    !note that this function is a dummy.
    ! For chemical cooling you should see
    ! heatingChem function in krome_heating.f90

    coolingChem = 0.d0

  end function coolingChem

  !**********************************
  function cooling_Continuum(n,Tgas)
    !cooling from continuum for a thin gas (no opacity)
    !see Omukai+2000 for details
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_getphys
    implicit none
    real*8::n(:),Tgas,cooling_Continuum,kgas,rhogas
    real*8::lj,tau,beta,m(nspec)

    m(:) = get_mass()
    rhogas = sum(n(1:nmols)*m(1:nmols)) !g/cm3
    kgas = kpla(n(:),Tgas) !planck opacity cm2/g (Omukai+2000)
    lj = get_jeans_length(n(:), Tgas) !cm
    tau = lj * kgas * rhogas + 1d-40 !opacity
    beta = min(1.d0,tau**(-2)) !beta escape (always <1.)
    cooling_Continuum = 4.d0 * stefboltz_erg * Tgas**4 &
        * kgas * rhogas * beta !erg/s/cm3

  end function cooling_Continuum

  !*******************************
  function cooling_compton(n, Tgas)
    !compton cooling erg/cm3/s from Cen1992
    use krome_user_commons
    use krome_commons
    real*8::cooling_compton,n(:),Tgas

    !note that redhsift is a common variable and
    ! should be provided by the user, otherwise the default is zero
    cooling_compton = 5.65d-36 * (1.d0 + phys_zredshift)**4 &
        * (Tgas - 2.73d0 * (1.d0 + phys_zredshift)) * n(idx_e) !erg/s/cm3

  end function cooling_compton

  !*****************
  !sigmoid function with x0 shift and s steepness
  function sigmoid(x,x0,s)
    implicit none
    real*8::sigmoid,x,x0,s

    sigmoid = 1d1/(1d1+exp(-s*(x-x0)))

  end function sigmoid

  !*******************
  !window function for H2 cooling to smooth limits
  function wCool(logTgas,logTmin,logTmax)
    implicit none
    real*8::wCool,logTgas,logTmin,logTmax,x

    x = (logTgas-logTmin)/(logTmax-logTmin)
    wCool = 1d1**(2d2*(sigmoid(x,-2d-1,5d1)*sigmoid(-x,-1.2d0,5d1)-1d0))
    if(wCool<1d-199) wCool = 0d0
    if(wCool>1d0) then
      print *,"ERROR: wCool>1"
      stop
    end if

  end function wCool

  !ALL THE COOLING FUNCTIONS ARE FROM GLOVER & ABEL, MNRAS 388, 1627, 2008
  !FOR LOW DENSITY REGIME: CONSIDER AN ORTHO-PARA RATIO OF 3:1
  !UPDATED TO THE DATA REPORTED BY GLOVER 2015, MNRAS
  !EACH SINGLE FUNCTION IS IN erg/s
  !FINAL UNITS = erg/cm3/s
  !*******************************
  function cooling_H2(n, Tgas)
    use krome_commons
    use krome_subs
    use krome_getphys
    real*8::n(:),Tgas
    real*8::temp,logt3,logt,cool,cooling_H2,T3
    real*8::LDL,HDLR,HDLV,HDL
    real*8::logt32,logt33,logt34,logt35,logt36,logt37,logt38
    real*8::dump14,fH2H,fH2e,fH2H2,fH2Hp,fH2He,w14,w24
    integer::i
    character*16::names(nspec)

    temp = Tgas
    cooling_H2 = 0d0
    !if(temp<2d0) return

    T3 = temp * 1.d-3
    logt3 = log10(T3)
    logt = log10(temp)
    cool = 0d0

    logt32 = logt3 * logt3
    logt33 = logt32 * logt3
    logt34 = logt33 * logt3
    logt35 = logt34 * logt3
    logt36 = logt35 * logt3
    logt37 = logt36 * logt3
    logt38 = logt37 * logt3

    w14 = wCool(logt, 1d0, 4d0)
    w24 = wCool(logt, 2d0, 4d0)

    !//H2-H
    if(temp<=1d2) then
      fH2H = 1.d1**(-16.818342D0 +3.7383713D1*logt3 &
          +5.8145166D1*logt32 +4.8656103D1*logt33 &
          +2.0159831D1*logt34 +3.8479610D0*logt35 )*n(idx_H)
    elseif(temp>1d2 .and. temp<=1d3) then
      fH2H = 1.d1**(-2.4311209D1 +3.5692468D0*logt3 &
          -1.1332860D1*logt32 -2.7850082D1*logt33 &
          -2.1328264D1*logt34 -4.2519023D0*logt35)*n(idx_H)
    elseif(temp>1.d3.and.temp<=6d3) then
      fH2H = 1d1**(-2.4311209D1 +4.6450521D0*logt3 &
          -3.7209846D0*logt32 +5.9369081D0*logt33 &
          -5.5108049D0*logt34 +1.5538288D0*logt35)*n(idx_H)
    else
      fH2H = 1.862314467912518E-022*wCool(logt,1d0,log10(6d3))*n(idx_H)
    end if
    cool = cool + fH2H

    !//H2-Hp
    if(temp>1d1.and.temp<=1d4) then
      fH2Hp = 1d1**(-2.2089523d1 +1.5714711d0*logt3 &
          +0.015391166d0*logt32 -0.23619985d0*logt33 &
          -0.51002221d0*logt34 +0.32168730d0*logt35)*n(idx_Hj)
    else
      fH2Hp = 1.182509139382060E-021*n(idx_Hj)*w14
    endif
    cool = cool + fH2Hp

    !//H2-H2
    fH2H2 = w24*1d1**(-2.3962112D1 +2.09433740D0*logt3 &
        -.77151436D0*logt32 +.43693353D0*logt33 &
        -.14913216D0*logt34 -.033638326D0*logt35)*n(idx_H2) !&
        cool = cool + fH2H2

    !//H2-e
    fH2e = 0d0
    if(temp<=5d2) then
      fH2e = 1d1**(min(-2.1928796d1 + 1.6815730d1*logt3 &
          +9.6743155d1*logt32 +3.4319180d2*logt33 &
          +7.3471651d2*logt34 +9.8367576d2*logt35 &
          +8.0181247d2*logt36 +3.6414446d2*logt37 &
          +7.0609154d1*logt38,3d1))*n(idx_e)
    elseif(temp>5d2)  then
      fH2e = 1d1**(-2.2921189D1 +1.6802758D0*logt3 &
          +.93310622D0*logt32 +4.0406627d0*logt33 &
          -4.7274036d0*logt34 -8.8077017d0*logt35 &
          +8.9167183*logt36 + 6.4380698*logt37 &
          -6.3701156*logt38)*n(idx_e)
    end if
    cool = cool + fH2e*w24

    !//H2-He
    if(temp>1d1.and.temp<=1d4)then
      fH2He = 1d1**(-2.3689237d1 +2.1892372d0*logt3&
          -.81520438d0*logt32 +.29036281d0*logt33 -.16596184d0*logt34 &
          +.19191375d0*logt35)*n(idx_He)
    else
      fH2He = 1.002560385050777E-022*n(idx_He)*w14
    endif
    cool = cool + fH2He

    !check error
    if(cool>1.d30) then
      print *,"ERROR: cooling >1.d30 erg/s/cm3"
      print *,"cool (erg/s/cm3): ",cool
      names(:) = get_names()
      do i=1,size(n)
        print '(I3,a18,E11.3)',i,names(i),n(i)
      end do
      stop
    end if

    !this to avoid negative, overflow and useless calculations below
    if(cool<=0d0) then
      cooling_H2 = 0d0
      return
    end if

    !high density limit from HM79, GP98 below Tgas = 2d3
    !UPDATED USING GLOVER 2015 for high temperature corrections, MNRAS
    !IN THE HIGH DENSITY REGIME LAMBDA_H2 = LAMBDA_H2(LTE) = HDL
    !the following mix of functions ensures the right behaviour
    ! at low (T<10 K) and high temperatures (T>2000 K) by
    ! using both the original Hollenbach and the new Glover data
    ! merged in a smooth way.
    if(temp.lt.2d3)then
      HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*exp(-(0.13/t3)**3)+&
          3.e-24*exp(-0.51/t3)) !erg/s
      HDLV = (6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3)) !erg/s
      HDL  = HDLR + HDLV !erg/s
    elseif(temp>=2d3 .and. temp<=1d4)then
      HDL = 1d1**(-2.0584225d1 + 5.0194035*logt3 &
          -1.5738805*logt32 -4.7155769*logt33 &
          +2.4714161*logt34 +5.4710750*logt35 &
          -3.9467356*logt36 -2.2148338*logt37 &
          +1.8161874*logt38)
    else
      dump14 = 1d0 / (1d0 + exp(min((temp-3d4)*2d-4,3d2)))
      HDL = 5.531333679406485E-019*dump14
    endif

    LDL = cool !erg/s
    if (HDL==0.) then
      cooling_H2 = 0.d0
    else
      cooling_H2 = n(idx_H2)/(1.d0/HDL+1.d0/LDL) &
          * min(1.d0, max(1.25d-10 * sum(n(1:nmols)),1d-40)**(-.45)) !erg/cm3/s
    endif

  end function cooling_H2

  !Atomic COOLING  Cen ApJS, 78, 341, 1992
  !UNITS = erg/s/cm3
  !*******************************
  function cooling_Atomic(n, Tgas)
    use krome_commons
    use krome_subs
    real*8::Tgas,cooling_atomic,n(:)
    real*8::temp,T5,cool

    temp = max(Tgas,10d0) !K
    T5 = temp/1d5 !K
    cool = 0d0 !erg/cm3/s

    !COLLISIONAL IONIZATION: H, He, He+, He(2S)
    cool = cool+ 1.27d-21*sqrt(temp)/(1.d0+sqrt(T5))&
        *exp(-1.578091d5/temp)*n(idx_e)*n(idx_H)

    cool = cool+ 9.38d-22*sqrt(temp)/(1.d0+sqrt(T5))&
        *exp(-2.853354d5/temp)*n(idx_e)*n(idx_He)

    !RECOMBINATION: H+, He+,He2+
    cool = cool+ 8.7d-27*sqrt(temp)*(temp/1.d3)**(-0.2)&
        /(1.d0+(temp/1.d6)**0.7)*n(idx_e)*n(idx_Hj)

    ! !DIELECTRONIC RECOMBINATION: He
    ! cool = cool+ 1.24d-13*temp**(-1.5)*exp(-4.7d5/temp)&
    !     *(1.d0+0.3d0*exp(-9.4d4/temp))*n(idx_e)*n(idx_Hej)

    !COLLISIONAL EXCITATION:
    !H(all n), He(n=2,3,4 triplets), He+(n=2)
    cool = cool+ 7.5d-19/(1.d0+sqrt(T5))*exp(-1.18348d5/temp)*n(idx_e)*n(idx_H)

    cooling_atomic = max(cool, 0d0)  !erg/cm3/s

  end function cooling_Atomic

  !HD COOLING LIPOVKA ET AL. MNRAS, 361, 850, (2005)
  !UNITS=erg/cm3/s
  !*******************************
  function cooling_HD(n, inTgas)
    use krome_commons
    use krome_subs
    implicit none
    integer::i,j
    integer, parameter::ns=4
    real*8::cooling_HD
    real*8::n(:),Tgas,logTgas,lognH,inTgas
    real*8::c(0:ns,0:ns),logW,W,dd,lhj

    !default HD cooling value
    cooling_HD = 0.0d0 !erg/cm3/s

    !this function does not have limits on density
    ! and temperature, even if the original paper do.
    ! However, we extrapolate the limits.

    !exit on low temperature
    if(inTgas<phys_Tcmb) return
    !extrapolate higher temperature limit
    Tgas = min(inTgas,1d4)

    !calculate density
    dd = n(idx_H) !sum(n(1:nmols))
    !exit if density is out of Lipovka bounds (uncomment if needed)
    !if(dd<1d0 .or. dd>1d8) return

    !extrapolate density limits
    dd = min(max(dd,1d-2),1d10)

    !POLYNOMIAL COEFFICIENT: TABLE 1 LIPOVKA
    c(0,:) = (/-42.56788d0, 0.92433d0, 0.54962d0, -0.07676d0, 0.00275d0/)
    c(1,:) = (/21.93385d0, 0.77952d0, -1.06447d0, 0.11864d0, -0.00366d0/)
    c(2,:) = (/-10.19097d0, -0.54263d0, 0.62343d0, -0.07366d0, 0.002514d0/)
    c(3,:) = (/2.19906d0, 0.11711d0, -0.13768d0, 0.01759d0, -0.00066631d0/)
    c(4,:) = (/-0.17334d0, -0.00835d0, 0.0106d0, -0.001482d0, 0.00006192d0/)

    logTgas = log10(Tgas)
    lognH   = log10(dd)

    !loop to compute coefficients
    logW = 0.d0
    do j = 0, ns
      lHj = lognH**j
      do i = 0, ns
        logW = logW + c(i,j)*logTgas**i*lHj !erg/s
      enddo
    enddo

    W = 10.d0**(logW)
    cooling_HD = W * n(idx_HD) !erg/cm3/s

  end function cooling_HD

  !***********************
  subroutine mylin2(a,b)
    !solve Ax=B analytically for a 2-levels system
    implicit none
    integer,parameter::n=2
    real*8::a(n,n),b(n),c(n),iab

    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,1)) then
    !   print *,"ERROR: a22=a21 in mylin2"
    !   stop
    !end if
    iab = b(1)/(a(2,2)-a(2,1))
    c(1) = a(2,2) * iab
    c(2) = -a(2,1) * iab
    b(:) = c(:)

  end subroutine mylin2

  !************************
  subroutine mylin3(a,b)
    !solve Ax=B analytically for a 3-levels system
    implicit none
    integer,parameter::n=3
    real*8::iab,a(n,n),b(n),c(n)

    !uncomment this: safer but slower function
    !if(a(2,2)==a(2,3)) then
    !   print *,"ERROR: a22=a23 in mylin3"
    !   stop
    !end if

    !uncomment this: safer but slower
    !if(a(2,1)*a(3,2)+a(2,2)*a(3,3)+a(2,3)*a(3,1) == &
        !     a(2,1)*a(3,3)+a(2,2)*a(3,1)+a(2,3)*a(3,2)) then
    !   print *,"ERROR: division by zero in mylin3"
    !   stop
    !end if

    iab = b(1) / (a(2,1)*(a(3,3)-a(3,2)) + a(2,2)*(a(3,1)-a(3,3)) &
        + a(2,3)*(a(3,2)-a(3,1)))
    c(1) = (a(2,3)*a(3,2)-a(2,2)*a(3,3)) * iab
    c(2) = -(a(2,3)*a(3,1)-a(2,1)*a(3,3)) * iab
    c(3) = (a(3,1)*a(2,2)-a(2,1)*a(3,2)) * iab
    b(:) = c(:)

  end subroutine mylin3

  !************************************
  subroutine plot_cool(n)
    !routine to plot cooling at runtime
    real*8::n(:),Tgas,Tmin,Tmax
    real*8::cool_atomic,cool_H2,cool_HD,cool_tot, cool_totGP,cool_H2GP
    real*8::cool_dH,cool_Z
    integer::i,imax
    imax = 1000
    Tmin = log10(1d1)
    Tmax = log10(1d8)
    print *,"plotting cooling..."
    open(33,file="KROME_cooling_plot.dat",status="replace")
    do i=1,imax
      Tgas = 1d1**(i*(Tmax-Tmin)/imax+Tmin)
      cool_H2 = 0.d0
      cool_H2GP = 0.d0
      cool_HD = 0.d0
      cool_atomic = 0.d0
      cool_Z = 0.d0
      cool_dH = 0.d0
      cool_H2 = cooling_H2(n(:),Tgas)
      cool_atomic = cooling_atomic(n(:),Tgas)
      cool_HD = cooling_HD(n(:),Tgas)
      cool_tot = cool_H2 + cool_atomic + cool_HD + cool_Z + cool_dH
      cool_totGP = cool_H2GP + cool_atomic + cool_HD + cool_Z + cool_dH
      write(33,'(99E12.3e3)') Tgas, cool_tot, cool_totGP, cool_H2, &
          cool_atomic, cool_HD, cool_H2GP, cool_Z, cool_dH
    end do
    close(33)
    print *,"done!"

  end subroutine plot_cool

  !***********************************
  !routine to dump cooling in unit nfile
  subroutine dump_cool(n,Tgas,nfile)
    use krome_commons
    implicit none
    real*8::Tgas,n(:),cools(ncools)
    integer::nfile

    cools(:) = get_cooling_array(n(:),Tgas)
    write(nfile,'(99E14.5e3)') Tgas, sum(cools), cools(:)

  end subroutine dump_cool

end module KROME_cooling

!############### MODULE ##############
module KROME_heating
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2022-06-08 12:22:52
  !  Changeset d0f9dcf
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !************************
  function heating(n,inTgas,k,nH2dust)
    implicit none
    real*8::n(:), Tgas, inTgas, k(:), nH2dust
    real*8::heating

    Tgas = inTgas
    heating = sum(get_heating_array(n(:),Tgas,k(:), nH2dust))

  end function heating

  !*******************************
  function get_heating_array(n, Tgas, k, nH2dust)
    use krome_commons
    implicit none
    real*8::n(:), Tgas, k(:), nH2dust
    real*8::get_heating_array(nheats),heats(nheats)
    real*8::smooth,f1,f2
    !returns heating in erg/cm3/s

    heats(:) = 0.d0

    f2 = 1.

    heats(idx_heat_chem) = heatingChem(n(:), Tgas, k(:), nH2dust)

    heats(idx_heat_custom) = heat_custom(n(:),Tgas)

    get_heating_array(:) = heats(:)

  end function get_heating_array

  !*************************
  function heat_custom(n,Tgas)
    use krome_commons
    use krome_subs
    use krome_constants
    implicit none
    real*8::n(:),Tgas,heat_custom

    heat_custom = 0d0

  end function heat_custom

  !H2 FORMATION HEATING and other exo/endothermic
  ! processes (including H2 on dust) in erg/cm3/s
  !krome builds the heating/cooling term according
  ! to the chemical network employed
  !*******************************
  function heatingChem(n, Tgas, k, nH2dust)
    use krome_constants
    use krome_commons
    use krome_dust
    use krome_subs
    use krome_getphys
    implicit none
    real*8::heatingChem, n(:), Tgas,k(:),nH2dust
    real*8::h2heatfac,HChem,yH,yH2
    real*8::ncr,ncrn,ncrd1,ncrd2,dd,n2H,small,nmax
    dd = get_Hnuclei(n(:))

    !replace small according to the desired enviroment
    ! and remove nmax if needed
    nmax = maxval(n(1:nmols))
    small = 1d-40/(nmax*nmax*nmax)

    heatingChem = 0.d0

    ncrn  = 1.0d6*(Tgas**(-0.5d0))
    ncrd1 = 1.6d0*exp(-(4.0d2/Tgas)**2)
    ncrd2 = 1.4d0*exp(-1.2d4/(Tgas+1.2d3))

    yH = n(idx_H)/dd   !dimensionless
    yH2= n(idx_H2)/dd  !dimensionless

    ncr = ncrn/(ncrd1*yH+ncrd2*yH2)      !1/cm3
    h2heatfac = 1.0d0/(1.0d0+ncr/dd)     !dimensionless

    HChem = 0.d0 !inits chemical heating
    n2H = n(idx_H) * n(idx_H)

    !H- + H -> H2 + E (heating)
    HChem = HChem + k(5) * (3.53d0*h2heatfac*n(idx_Hk)*n(idx_H))
    !H2+ + H -> H2 + H+ (heating)
    HChem = HChem + k(9) * (1.83d0*h2heatfac*n(idx_H2j)*n(idx_H))
    !H + H + H -> H2 + H (heating)
    HChem = HChem + k(24) * (4.48d0*h2heatfac*n(idx_H)*n(idx_H)*n(idx_H))
    !H + H2 -> H + H + H (cooling)
    HChem = HChem + k(25) * (-4.48d0*n(idx_H)*n(idx_H2))

    heatingChem = HChem * eV_to_erg  !erg/cm3/s

  end function heatingChem

end module KROME_heating

!############### MODULE ##############
module krome_ode
contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2022-06-08 12:22:52
  !  Changeset d0f9dcf
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  subroutine fex(neq,tt,nin,dn)
    use krome_commons
    use krome_constants
    use krome_subs
    use krome_cooling
    use krome_heating
    use krome_tabs
    use krome_photo
    use krome_gadiab
    use krome_getphys
    use krome_phfuncs
    use krome_fit
    implicit none
    integer::neq,idust
    real*8::tt,dn(neq),n(neq),k(nrea),krome_gamma
    real*8::gamma,Tgas,vgas,ntot,nH2dust,nd,nin(neq)
    real*8::rr
    integer::i,r1,r2,r3,p1,p2,p3

    n(:) = nin(:)

    nH2dust = 0.d0
    n(idx_CR) = 1.d0
    n(idx_g)  = 1.d0
    n(idx_dummy) = 1.d0

    dn(:) = 0.d0 !initialize differentials
    n(idx_Tgas) = max(n(idx_tgas),2.73d0)
    n(idx_Tgas) = min(n(idx_tgas),1d9)
    Tgas = n(idx_Tgas) !get temperature

    k(:) = coe_tab(n(:)) !compute coefficients

    !E
    !E
    dn(idx_E) = &
        -k(1)*n(idx_Hj)*n(idx_E) &
        +k(2)*n(idx_H) &
        +k(3)*n(idx_Hk) &
        -k(4)*n(idx_H)*n(idx_E) &
        +k(5)*n(idx_Hk)*n(idx_H) &
        -k(11)*n(idx_Dj)*n(idx_E) &
        +k(12)*n(idx_D) &
        +k(19)*n(idx_Hk)*n(idx_D) &
        -k(22)*n(idx_HEHj)*n(idx_E)

    !H-
    !H-
    dn(idx_Hk) = &
        -k(3)*n(idx_Hk) &
        +k(4)*n(idx_H)*n(idx_E) &
        -k(5)*n(idx_Hk)*n(idx_H) &
        -k(6)*n(idx_Hk)*n(idx_Hj) &
        -k(19)*n(idx_Hk)*n(idx_D)

    !H
    !H
    dn(idx_H) = &
        +k(1)*n(idx_Hj)*n(idx_E) &
        -k(2)*n(idx_H) &
        +k(3)*n(idx_Hk) &
        -k(4)*n(idx_H)*n(idx_E) &
        -k(5)*n(idx_Hk)*n(idx_H) &
        +2.d0*k(6)*n(idx_Hk)*n(idx_Hj) &
        -k(7)*n(idx_H)*n(idx_Hj) &
        +k(8)*n(idx_H2j) &
        -k(9)*n(idx_H2j)*n(idx_H) &
        +k(10)*n(idx_H2)*n(idx_Hj) &
        +k(13)*n(idx_D)*n(idx_Hj) &
        -k(14)*n(idx_Dj)*n(idx_H) &
        +k(17)*n(idx_D)*n(idx_H2) &
        -k(18)*n(idx_HD)*n(idx_H) &
        -k(21)*n(idx_HEHj)*n(idx_H) &
        +k(22)*n(idx_HEHj)*n(idx_E) &
        -3.d0*k(24)*n(idx_H)*n(idx_H)*n(idx_H) &
        +k(24)*n(idx_H)*n(idx_H)*n(idx_H) &
        -k(25)*n(idx_H)*n(idx_H2) &
        +3.d0*k(25)*n(idx_H)*n(idx_H2)

    !H2
    !H2
    dn(idx_H2) = &
        +k(5)*n(idx_Hk)*n(idx_H) &
        +k(9)*n(idx_H2j)*n(idx_H) &
        -k(10)*n(idx_H2)*n(idx_Hj) &
        -k(15)*n(idx_Dj)*n(idx_H2) &
        +k(16)*n(idx_HD)*n(idx_Hj) &
        -k(17)*n(idx_D)*n(idx_H2) &
        +k(18)*n(idx_HD)*n(idx_H) &
        +k(24)*n(idx_H)*n(idx_H)*n(idx_H) &
        -k(25)*n(idx_H)*n(idx_H2)

    !D
    !D
    dn(idx_D) = &
        +k(11)*n(idx_Dj)*n(idx_E) &
        -k(12)*n(idx_D) &
        -k(13)*n(idx_D)*n(idx_Hj) &
        +k(14)*n(idx_Dj)*n(idx_H) &
        -k(17)*n(idx_D)*n(idx_H2) &
        +k(18)*n(idx_HD)*n(idx_H) &
        -k(19)*n(idx_Hk)*n(idx_D)

    !HD
    !HD
    dn(idx_HD) = &
        +k(15)*n(idx_Dj)*n(idx_H2) &
        -k(16)*n(idx_HD)*n(idx_Hj) &
        +k(17)*n(idx_D)*n(idx_H2) &
        -k(18)*n(idx_HD)*n(idx_H) &
        +k(19)*n(idx_Hk)*n(idx_D)

    !HE
    !HE
    dn(idx_HE) = &
        -k(20)*n(idx_HE)*n(idx_Hj) &
        +k(21)*n(idx_HEHj)*n(idx_H) &
        +k(22)*n(idx_HEHj)*n(idx_E) &
        +k(23)*n(idx_HEHj)

    !H+
    !H+
    dn(idx_Hj) = &
        -k(1)*n(idx_Hj)*n(idx_E) &
        +k(2)*n(idx_H) &
        -k(6)*n(idx_Hk)*n(idx_Hj) &
        -k(7)*n(idx_H)*n(idx_Hj) &
        +k(8)*n(idx_H2j) &
        +k(9)*n(idx_H2j)*n(idx_H) &
        -k(10)*n(idx_H2)*n(idx_Hj) &
        -k(13)*n(idx_D)*n(idx_Hj) &
        +k(14)*n(idx_Dj)*n(idx_H) &
        +k(15)*n(idx_Dj)*n(idx_H2) &
        -k(16)*n(idx_HD)*n(idx_Hj) &
        -k(20)*n(idx_HE)*n(idx_Hj) &
        +k(23)*n(idx_HEHj)

    !H2+
    !H2+
    dn(idx_H2j) = &
        +k(7)*n(idx_H)*n(idx_Hj) &
        -k(8)*n(idx_H2j) &
        -k(9)*n(idx_H2j)*n(idx_H) &
        +k(10)*n(idx_H2)*n(idx_Hj) &
        +k(21)*n(idx_HEHj)*n(idx_H)

    !D+
    !D+
    dn(idx_Dj) = &
        -k(11)*n(idx_Dj)*n(idx_E) &
        +k(12)*n(idx_D) &
        +k(13)*n(idx_D)*n(idx_Hj) &
        -k(14)*n(idx_Dj)*n(idx_H) &
        -k(15)*n(idx_Dj)*n(idx_H2) &
        +k(16)*n(idx_HD)*n(idx_Hj)

    !HEH+
    !HEH+
    dn(idx_HEHj) = &
        +k(20)*n(idx_HE)*n(idx_Hj) &
        -k(21)*n(idx_HEHj)*n(idx_H) &
        -k(22)*n(idx_HEHj)*n(idx_E) &
        -k(23)*n(idx_HEHj)

    !CR

    !CR
    dn(idx_CR) = 0.d0

    !g

    !g
    dn(idx_g) = 0.d0

    !Tgas

    !Tgas
    dn(idx_Tgas) = 0.d0

    !dummy

    !dummy
    dn(idx_dummy) = 0.d0

    krome_gamma = gamma_index(n(:))

    dn(idx_Tgas) = (heating(n(:), Tgas, k(:), nH2dust) &
        - cooling(n(:), Tgas)  ) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

    last_coe(:) = k(:)

  end subroutine fex

  !***************************
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use krome_commons
    use krome_subs
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_gadiab
    use krome_getphys
    implicit none
    integer::neq, j, ian, jan, r1, r2, p1, p2, p3, i
    real*8::tt, n(neq), pdj(neq), dr1, dr2, kk,k(nrea),Tgas
    real*8::nn(neq),dn0,dn1,dnn,nH2dust,dn(neq),krome_gamma

    nH2dust = 0.d0
    Tgas = n(idx_Tgas)

    krome_gamma = gamma_index(n(:))

    k(:) = last_coe(:) !get rate coefficients

    if(j==1) then
    elseif(j==1) then
      pdj(1) =  &
          -k(1)*n(idx_Hj)  &
          -k(4)*n(idx_H)  &
          -k(11)*n(idx_Dj)  &
          -k(22)*n(idx_HEHj)
      pdj(2) =  &
          +k(4)*n(idx_H)
      pdj(3) =  &
          +k(1)*n(idx_Hj)  &
          -k(4)*n(idx_H)  &
          +k(22)*n(idx_HEHj)
      pdj(5) =  &
          +k(11)*n(idx_Dj)
      pdj(7) =  &
          +k(22)*n(idx_HEHj)
      pdj(8) =  &
          -k(1)*n(idx_Hj)
      pdj(10) =  &
          -k(11)*n(idx_Dj)
      pdj(11) =  &
          -k(22)*n(idx_HEHj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(1)*1d-3
      if(dnn>0.d0) then
        nn(1) = n(1) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==2) then
      pdj(1) =  &
          +k(3)  &
          +k(5)*n(idx_H)  &
          +k(19)*n(idx_D)
      pdj(2) =  &
          -k(3)  &
          -k(5)*n(idx_H)  &
          -k(6)*n(idx_Hj)  &
          -k(19)*n(idx_D)
      pdj(3) =  &
          +k(3)  &
          -k(5)*n(idx_H)  &
          +2.d0*k(6)*n(idx_Hj)
      pdj(4) =  &
          +k(5)*n(idx_H)
      pdj(5) =  &
          -k(19)*n(idx_D)
      pdj(6) =  &
          +k(19)*n(idx_D)
      pdj(8) =  &
          -k(6)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(2)*1d-3
      if(dnn>0.d0) then
        nn(2) = n(2) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==3) then
      pdj(1) =  &
          +k(2)  &
          -k(4)*n(idx_E)  &
          +k(5)*n(idx_Hk)
      pdj(2) =  &
          +k(4)*n(idx_E)  &
          -k(5)*n(idx_Hk)
      pdj(3) =  &
          -k(2)  &
          -k(4)*n(idx_E)  &
          -k(5)*n(idx_Hk)  &
          -k(7)*n(idx_Hj)  &
          -k(9)*n(idx_H2j)  &
          -k(14)*n(idx_Dj)  &
          -k(18)*n(idx_HD)  &
          -k(21)*n(idx_HEHj)  &
          -9.d0*k(24)*n(idx_H)*n(idx_H)  &
          +3.d0*k(24)*n(idx_H)*n(idx_H)  &
          -k(25)*n(idx_H2)  &
          +3.d0*k(25)*n(idx_H2)
      pdj(4) =  &
          +k(5)*n(idx_Hk)  &
          +k(9)*n(idx_H2j)  &
          +k(18)*n(idx_HD)  &
          +3.d0*k(24)*n(idx_H)*n(idx_H)  &
          -k(25)*n(idx_H2)
      pdj(5) =  &
          +k(14)*n(idx_Dj)  &
          +k(18)*n(idx_HD)
      pdj(6) =  &
          -k(18)*n(idx_HD)
      pdj(7) =  &
          +k(21)*n(idx_HEHj)
      pdj(8) =  &
          +k(2)  &
          -k(7)*n(idx_Hj)  &
          +k(9)*n(idx_H2j)  &
          +k(14)*n(idx_Dj)
      pdj(9) =  &
          +k(7)*n(idx_Hj)  &
          -k(9)*n(idx_H2j)  &
          +k(21)*n(idx_HEHj)
      pdj(10) =  &
          -k(14)*n(idx_Dj)
      pdj(11) =  &
          -k(21)*n(idx_HEHj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(3)*1d-3
      if(dnn>0.d0) then
        nn(3) = n(3) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==4) then
      pdj(3) =  &
          +k(10)*n(idx_Hj)  &
          +k(17)*n(idx_D)  &
          -k(25)*n(idx_H)  &
          +3.d0*k(25)*n(idx_H)
      pdj(4) =  &
          -k(10)*n(idx_Hj)  &
          -k(15)*n(idx_Dj)  &
          -k(17)*n(idx_D)  &
          -k(25)*n(idx_H)
      pdj(5) =  &
          -k(17)*n(idx_D)
      pdj(6) =  &
          +k(15)*n(idx_Dj)  &
          +k(17)*n(idx_D)
      pdj(8) =  &
          -k(10)*n(idx_Hj)  &
          +k(15)*n(idx_Dj)
      pdj(9) =  &
          +k(10)*n(idx_Hj)
      pdj(10) =  &
          -k(15)*n(idx_Dj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(4)*1d-3
      if(dnn>0.d0) then
        nn(4) = n(4) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==5) then
      pdj(1) =  &
          +k(12)  &
          +k(19)*n(idx_Hk)
      pdj(2) =  &
          -k(19)*n(idx_Hk)
      pdj(3) =  &
          +k(13)*n(idx_Hj)  &
          +k(17)*n(idx_H2)
      pdj(4) =  &
          -k(17)*n(idx_H2)
      pdj(5) =  &
          -k(12)  &
          -k(13)*n(idx_Hj)  &
          -k(17)*n(idx_H2)  &
          -k(19)*n(idx_Hk)
      pdj(6) =  &
          +k(17)*n(idx_H2)  &
          +k(19)*n(idx_Hk)
      pdj(8) =  &
          -k(13)*n(idx_Hj)
      pdj(10) =  &
          +k(12)  &
          +k(13)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(5)*1d-3
      if(dnn>0.d0) then
        nn(5) = n(5) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==6) then
      pdj(3) =  &
          -k(18)*n(idx_H)
      pdj(4) =  &
          +k(16)*n(idx_Hj)  &
          +k(18)*n(idx_H)
      pdj(5) =  &
          +k(18)*n(idx_H)
      pdj(6) =  &
          -k(16)*n(idx_Hj)  &
          -k(18)*n(idx_H)
      pdj(8) =  &
          -k(16)*n(idx_Hj)
      pdj(10) =  &
          +k(16)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(6)*1d-3
      if(dnn>0.d0) then
        nn(6) = n(6) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==7) then
      pdj(7) =  &
          -k(20)*n(idx_Hj)
      pdj(8) =  &
          -k(20)*n(idx_Hj)
      pdj(11) =  &
          +k(20)*n(idx_Hj)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(7)*1d-3
      if(dnn>0.d0) then
        nn(7) = n(7) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==8) then
      pdj(1) =  &
          -k(1)*n(idx_E)
      pdj(2) =  &
          -k(6)*n(idx_Hk)
      pdj(3) =  &
          +k(1)*n(idx_E)  &
          +2.d0*k(6)*n(idx_Hk)  &
          -k(7)*n(idx_H)  &
          +k(10)*n(idx_H2)  &
          +k(13)*n(idx_D)
      pdj(4) =  &
          -k(10)*n(idx_H2)  &
          +k(16)*n(idx_HD)
      pdj(5) =  &
          -k(13)*n(idx_D)
      pdj(6) =  &
          -k(16)*n(idx_HD)
      pdj(7) =  &
          -k(20)*n(idx_HE)
      pdj(8) =  &
          -k(1)*n(idx_E)  &
          -k(6)*n(idx_Hk)  &
          -k(7)*n(idx_H)  &
          -k(10)*n(idx_H2)  &
          -k(13)*n(idx_D)  &
          -k(16)*n(idx_HD)  &
          -k(20)*n(idx_HE)
      pdj(9) =  &
          +k(7)*n(idx_H)  &
          +k(10)*n(idx_H2)
      pdj(10) =  &
          +k(13)*n(idx_D)  &
          +k(16)*n(idx_HD)
      pdj(11) =  &
          +k(20)*n(idx_HE)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(8)*1d-3
      if(dnn>0.d0) then
        nn(8) = n(8) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==9) then
      pdj(3) =  &
          +k(8)  &
          -k(9)*n(idx_H)
      pdj(4) =  &
          +k(9)*n(idx_H)
      pdj(8) =  &
          +k(8)  &
          +k(9)*n(idx_H)
      pdj(9) =  &
          -k(8)  &
          -k(9)*n(idx_H)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(9)*1d-3
      if(dnn>0.d0) then
        nn(9) = n(9) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==10) then
      pdj(1) =  &
          -k(11)*n(idx_E)
      pdj(3) =  &
          -k(14)*n(idx_H)
      pdj(4) =  &
          -k(15)*n(idx_H2)
      pdj(5) =  &
          +k(11)*n(idx_E)  &
          +k(14)*n(idx_H)
      pdj(6) =  &
          +k(15)*n(idx_H2)
      pdj(8) =  &
          +k(14)*n(idx_H)  &
          +k(15)*n(idx_H2)
      pdj(10) =  &
          -k(11)*n(idx_E)  &
          -k(14)*n(idx_H)  &
          -k(15)*n(idx_H2)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(10)*1d-3
      if(dnn>0.d0) then
        nn(10) = n(10) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==11) then
      pdj(1) =  &
          -k(22)*n(idx_E)
      pdj(3) =  &
          -k(21)*n(idx_H)  &
          +k(22)*n(idx_E)
      pdj(7) =  &
          +k(21)*n(idx_H)  &
          +k(22)*n(idx_E)  &
          +k(23)
      pdj(8) =  &
          +k(23)
      pdj(9) =  &
          +k(21)*n(idx_H)
      pdj(11) =  &
          -k(21)*n(idx_H)  &
          -k(22)*n(idx_E)  &
          -k(23)
      dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      nn(:) = n(:)
      dnn = n(11)*1d-3
      if(dnn>0.d0) then
        nn(11) = n(11) + dnn
        dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
            * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
        pdj(idx_Tgas) = (dn1-dn0)/dnn
      end if

    elseif(j==12) then
      pdj(14) = 0.d0
    elseif(j==13) then
      pdj(14) = 0.d0
    elseif(j==14) then
      !use fex to compute temperature-dependent Jacobian
      dnn = n(idx_Tgas)*1d-3
      nn(:) = n(:)
      nn(idx_Tgas) = n(idx_Tgas) + dnn
      call fex(neq,tt,nn(:),dn(:))
      do i=1,neq-1
        pdj(i) = dn(i) / dnn
      end do
    elseif(j==15) then
      pdj(14) = 0.d0
    end if

    return
  end subroutine jes

  !*************************
  subroutine jex(neq,t,n,ml,mu,pd,npd)
    use krome_commons
    use krome_tabs
    use krome_cooling
    use krome_heating
    use krome_constants
    use krome_subs
    use krome_gadiab
    implicit none
    real*8::n(neq),pd(neq,neq),t,k(nrea),dn0,dn1,dnn,Tgas
    real*8::krome_gamma,nn(neq),nH2dust
    integer::neq,ml,mu,npd

    Tgas = n(idx_Tgas)
    npd = neq
    k(:) = coe_tab(n(:))
    pd(:,:) = 0d0
    krome_gamma = gamma_index(n(:))

    !d[E_dot]/d[E]
    pd(1,1) =  &
        -k(1)*n(idx_Hj)  &
        -k(4)*n(idx_H)  &
        -k(11)*n(idx_Dj)  &
        -k(22)*n(idx_HEHj)

    !d[H-_dot]/d[E]
    pd(2,1) =  &
        +k(4)*n(idx_H)

    !d[H_dot]/d[E]
    pd(3,1) =  &
        +k(1)*n(idx_Hj)  &
        -k(4)*n(idx_H)  &
        +k(22)*n(idx_HEHj)

    !d[D_dot]/d[E]
    pd(5,1) =  &
        +k(11)*n(idx_Dj)

    !d[HE_dot]/d[E]
    pd(7,1) =  &
        +k(22)*n(idx_HEHj)

    !d[H+_dot]/d[E]
    pd(8,1) =  &
        -k(1)*n(idx_Hj)

    !d[D+_dot]/d[E]
    pd(10,1) =  &
        -k(11)*n(idx_Dj)

    !d[HEH+_dot]/d[E]
    pd(11,1) =  &
        -k(22)*n(idx_HEHj)

    !d[Tgas_dot]/d[E]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(1)*1d-3
    if(dnn>0.d0) then
      nn(1) = n(1) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,1) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H-]
    pd(1,2) =  &
        +k(3)  &
        +k(5)*n(idx_H)  &
        +k(19)*n(idx_D)

    !d[H-_dot]/d[H-]
    pd(2,2) =  &
        -k(3)  &
        -k(5)*n(idx_H)  &
        -k(6)*n(idx_Hj)  &
        -k(19)*n(idx_D)

    !d[H_dot]/d[H-]
    pd(3,2) =  &
        +k(3)  &
        -k(5)*n(idx_H)  &
        +2.d0*k(6)*n(idx_Hj)

    !d[H2_dot]/d[H-]
    pd(4,2) =  &
        +k(5)*n(idx_H)

    !d[D_dot]/d[H-]
    pd(5,2) =  &
        -k(19)*n(idx_D)

    !d[HD_dot]/d[H-]
    pd(6,2) =  &
        +k(19)*n(idx_D)

    !d[H+_dot]/d[H-]
    pd(8,2) =  &
        -k(6)*n(idx_Hj)

    !d[Tgas_dot]/d[H-]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(2)*1d-3
    if(dnn>0.d0) then
      nn(2) = n(2) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,2) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H]
    pd(1,3) =  &
        +k(2)  &
        -k(4)*n(idx_E)  &
        +k(5)*n(idx_Hk)

    !d[H-_dot]/d[H]
    pd(2,3) =  &
        +k(4)*n(idx_E)  &
        -k(5)*n(idx_Hk)

    !d[H_dot]/d[H]
    pd(3,3) =  &
        -k(2)  &
        -k(4)*n(idx_E)  &
        -k(5)*n(idx_Hk)  &
        -k(7)*n(idx_Hj)  &
        -k(9)*n(idx_H2j)  &
        -k(14)*n(idx_Dj)  &
        -k(18)*n(idx_HD)  &
        -k(21)*n(idx_HEHj)  &
        -9.d0*k(24)*n(idx_H)*n(idx_H)  &
        +3.d0*k(24)*n(idx_H)*n(idx_H)  &
        -k(25)*n(idx_H2)  &
        +3.d0*k(25)*n(idx_H2)

    !d[H2_dot]/d[H]
    pd(4,3) =  &
        +k(5)*n(idx_Hk)  &
        +k(9)*n(idx_H2j)  &
        +k(18)*n(idx_HD)  &
        +3.d0*k(24)*n(idx_H)*n(idx_H)  &
        -k(25)*n(idx_H2)

    !d[D_dot]/d[H]
    pd(5,3) =  &
        +k(14)*n(idx_Dj)  &
        +k(18)*n(idx_HD)

    !d[HD_dot]/d[H]
    pd(6,3) =  &
        -k(18)*n(idx_HD)

    !d[HE_dot]/d[H]
    pd(7,3) =  &
        +k(21)*n(idx_HEHj)

    !d[H+_dot]/d[H]
    pd(8,3) =  &
        +k(2)  &
        -k(7)*n(idx_Hj)  &
        +k(9)*n(idx_H2j)  &
        +k(14)*n(idx_Dj)

    !d[H2+_dot]/d[H]
    pd(9,3) =  &
        +k(7)*n(idx_Hj)  &
        -k(9)*n(idx_H2j)  &
        +k(21)*n(idx_HEHj)

    !d[D+_dot]/d[H]
    pd(10,3) =  &
        -k(14)*n(idx_Dj)

    !d[HEH+_dot]/d[H]
    pd(11,3) =  &
        -k(21)*n(idx_HEHj)

    !d[Tgas_dot]/d[H]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(3)*1d-3
    if(dnn>0.d0) then
      nn(3) = n(3) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,3) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[H2]
    pd(3,4) =  &
        +k(10)*n(idx_Hj)  &
        +k(17)*n(idx_D)  &
        -k(25)*n(idx_H)  &
        +3.d0*k(25)*n(idx_H)

    !d[H2_dot]/d[H2]
    pd(4,4) =  &
        -k(10)*n(idx_Hj)  &
        -k(15)*n(idx_Dj)  &
        -k(17)*n(idx_D)  &
        -k(25)*n(idx_H)

    !d[D_dot]/d[H2]
    pd(5,4) =  &
        -k(17)*n(idx_D)

    !d[HD_dot]/d[H2]
    pd(6,4) =  &
        +k(15)*n(idx_Dj)  &
        +k(17)*n(idx_D)

    !d[H+_dot]/d[H2]
    pd(8,4) =  &
        -k(10)*n(idx_Hj)  &
        +k(15)*n(idx_Dj)

    !d[H2+_dot]/d[H2]
    pd(9,4) =  &
        +k(10)*n(idx_Hj)

    !d[D+_dot]/d[H2]
    pd(10,4) =  &
        -k(15)*n(idx_Dj)

    !d[Tgas_dot]/d[H2]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(4)*1d-3
    if(dnn>0.d0) then
      nn(4) = n(4) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,4) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[D]
    pd(1,5) =  &
        +k(12)  &
        +k(19)*n(idx_Hk)

    !d[H-_dot]/d[D]
    pd(2,5) =  &
        -k(19)*n(idx_Hk)

    !d[H_dot]/d[D]
    pd(3,5) =  &
        +k(13)*n(idx_Hj)  &
        +k(17)*n(idx_H2)

    !d[H2_dot]/d[D]
    pd(4,5) =  &
        -k(17)*n(idx_H2)

    !d[D_dot]/d[D]
    pd(5,5) =  &
        -k(12)  &
        -k(13)*n(idx_Hj)  &
        -k(17)*n(idx_H2)  &
        -k(19)*n(idx_Hk)

    !d[HD_dot]/d[D]
    pd(6,5) =  &
        +k(17)*n(idx_H2)  &
        +k(19)*n(idx_Hk)

    !d[H+_dot]/d[D]
    pd(8,5) =  &
        -k(13)*n(idx_Hj)

    !d[D+_dot]/d[D]
    pd(10,5) =  &
        +k(12)  &
        +k(13)*n(idx_Hj)

    !d[Tgas_dot]/d[D]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(5)*1d-3
    if(dnn>0.d0) then
      nn(5) = n(5) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,5) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[HD]
    pd(3,6) =  &
        -k(18)*n(idx_H)

    !d[H2_dot]/d[HD]
    pd(4,6) =  &
        +k(16)*n(idx_Hj)  &
        +k(18)*n(idx_H)

    !d[D_dot]/d[HD]
    pd(5,6) =  &
        +k(18)*n(idx_H)

    !d[HD_dot]/d[HD]
    pd(6,6) =  &
        -k(16)*n(idx_Hj)  &
        -k(18)*n(idx_H)

    !d[H+_dot]/d[HD]
    pd(8,6) =  &
        -k(16)*n(idx_Hj)

    !d[D+_dot]/d[HD]
    pd(10,6) =  &
        +k(16)*n(idx_Hj)

    !d[Tgas_dot]/d[HD]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(6)*1d-3
    if(dnn>0.d0) then
      nn(6) = n(6) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,6) = (dn1-dn0)/dnn
    end if

    !d[HE_dot]/d[HE]
    pd(7,7) =  &
        -k(20)*n(idx_Hj)

    !d[H+_dot]/d[HE]
    pd(8,7) =  &
        -k(20)*n(idx_Hj)

    !d[HEH+_dot]/d[HE]
    pd(11,7) =  &
        +k(20)*n(idx_Hj)

    !d[Tgas_dot]/d[HE]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(7)*1d-3
    if(dnn>0.d0) then
      nn(7) = n(7) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,7) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[H+]
    pd(1,8) =  &
        -k(1)*n(idx_E)

    !d[H-_dot]/d[H+]
    pd(2,8) =  &
        -k(6)*n(idx_Hk)

    !d[H_dot]/d[H+]
    pd(3,8) =  &
        +k(1)*n(idx_E)  &
        +2.d0*k(6)*n(idx_Hk)  &
        -k(7)*n(idx_H)  &
        +k(10)*n(idx_H2)  &
        +k(13)*n(idx_D)

    !d[H2_dot]/d[H+]
    pd(4,8) =  &
        -k(10)*n(idx_H2)  &
        +k(16)*n(idx_HD)

    !d[D_dot]/d[H+]
    pd(5,8) =  &
        -k(13)*n(idx_D)

    !d[HD_dot]/d[H+]
    pd(6,8) =  &
        -k(16)*n(idx_HD)

    !d[HE_dot]/d[H+]
    pd(7,8) =  &
        -k(20)*n(idx_HE)

    !d[H+_dot]/d[H+]
    pd(8,8) =  &
        -k(1)*n(idx_E)  &
        -k(6)*n(idx_Hk)  &
        -k(7)*n(idx_H)  &
        -k(10)*n(idx_H2)  &
        -k(13)*n(idx_D)  &
        -k(16)*n(idx_HD)  &
        -k(20)*n(idx_HE)

    !d[H2+_dot]/d[H+]
    pd(9,8) =  &
        +k(7)*n(idx_H)  &
        +k(10)*n(idx_H2)

    !d[D+_dot]/d[H+]
    pd(10,8) =  &
        +k(13)*n(idx_D)  &
        +k(16)*n(idx_HD)

    !d[HEH+_dot]/d[H+]
    pd(11,8) =  &
        +k(20)*n(idx_HE)

    !d[Tgas_dot]/d[H+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(8)*1d-3
    if(dnn>0.d0) then
      nn(8) = n(8) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,8) = (dn1-dn0)/dnn
    end if

    !d[H_dot]/d[H2+]
    pd(3,9) =  &
        +k(8)  &
        -k(9)*n(idx_H)

    !d[H2_dot]/d[H2+]
    pd(4,9) =  &
        +k(9)*n(idx_H)

    !d[H+_dot]/d[H2+]
    pd(8,9) =  &
        +k(8)  &
        +k(9)*n(idx_H)

    !d[H2+_dot]/d[H2+]
    pd(9,9) =  &
        -k(8)  &
        -k(9)*n(idx_H)

    !d[Tgas_dot]/d[H2+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(9)*1d-3
    if(dnn>0.d0) then
      nn(9) = n(9) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,9) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[D+]
    pd(1,10) =  &
        -k(11)*n(idx_E)

    !d[H_dot]/d[D+]
    pd(3,10) =  &
        -k(14)*n(idx_H)

    !d[H2_dot]/d[D+]
    pd(4,10) =  &
        -k(15)*n(idx_H2)

    !d[D_dot]/d[D+]
    pd(5,10) =  &
        +k(11)*n(idx_E)  &
        +k(14)*n(idx_H)

    !d[HD_dot]/d[D+]
    pd(6,10) =  &
        +k(15)*n(idx_H2)

    !d[H+_dot]/d[D+]
    pd(8,10) =  &
        +k(14)*n(idx_H)  &
        +k(15)*n(idx_H2)

    !d[D+_dot]/d[D+]
    pd(10,10) =  &
        -k(11)*n(idx_E)  &
        -k(14)*n(idx_H)  &
        -k(15)*n(idx_H2)

    !d[Tgas_dot]/d[D+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(10)*1d-3
    if(dnn>0.d0) then
      nn(10) = n(10) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,10) = (dn1-dn0)/dnn
    end if

    !d[E_dot]/d[HEH+]
    pd(1,11) =  &
        -k(22)*n(idx_E)

    !d[H_dot]/d[HEH+]
    pd(3,11) =  &
        -k(21)*n(idx_H)  &
        +k(22)*n(idx_E)

    !d[HE_dot]/d[HEH+]
    pd(7,11) =  &
        +k(21)*n(idx_H)  &
        +k(22)*n(idx_E)  &
        +k(23)

    !d[H+_dot]/d[HEH+]
    pd(8,11) =  &
        +k(23)

    !d[H2+_dot]/d[HEH+]
    pd(9,11) =  &
        +k(21)*n(idx_H)

    !d[HEH+_dot]/d[HEH+]
    pd(11,11) =  &
        -k(21)*n(idx_H)  &
        -k(22)*n(idx_E)  &
        -k(23)

    !d[Tgas_dot]/d[HEH+]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(11)*1d-3
    if(dnn>0.d0) then
      nn(11) = n(11) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,11) = (dn1-dn0)/dnn
    end if

    !d[Tgas_dot]/d[CR]
    pd(14,12) = 0.d0

    !d[Tgas_dot]/d[g]
    pd(14,13) = 0.d0

    !d[E_dot]/d[Tgas]
    pd(1,14) = 0.d0

    !d[H-_dot]/d[Tgas]
    pd(2,14) = 0.d0

    !d[H_dot]/d[Tgas]
    pd(3,14) = 0.d0

    !d[H2_dot]/d[Tgas]
    pd(4,14) = 0.d0

    !d[D_dot]/d[Tgas]
    pd(5,14) = 0.d0

    !d[HD_dot]/d[Tgas]
    pd(6,14) = 0.d0

    !d[HE_dot]/d[Tgas]
    pd(7,14) = 0.d0

    !d[H+_dot]/d[Tgas]
    pd(8,14) = 0.d0

    !d[H2+_dot]/d[Tgas]
    pd(9,14) = 0.d0

    !d[D+_dot]/d[Tgas]
    pd(10,14) = 0.d0

    !d[HEH+_dot]/d[Tgas]
    pd(11,14) = 0.d0

    !d[CR_dot]/d[Tgas]
    pd(12,14) = 0.d0

    !d[g_dot]/d[Tgas]
    pd(13,14) = 0.d0

    !d[Tgas_dot]/d[Tgas]
    dn0 = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
    nn(:) = n(:)
    dnn = n(14)*1d-3
    if(dnn>0.d0) then
      nn(14) = n(14) + dnn
      dn1 = (heating(nn(:), Tgas, k(:), nH2dust) - cooling(nn(:), Tgas)) &
          * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))
      pd(idx_Tgas,14) = (dn1-dn0)/dnn
    end if

    !d[dummy_dot]/d[Tgas]
    pd(15,14) = 0.d0

    !d[Tgas_dot]/d[dummy]
    pd(14,15) = 0.d0

  end subroutine jex

end module krome_ode

!############### MODULE ##############
module krome_user
  implicit none

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2022-06-08 12:22:52
  !  Changeset d0f9dcf
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  integer,parameter::KROME_idx_E = 1	!E
  integer,parameter::KROME_idx_Hk = 2	!H-
  integer,parameter::KROME_idx_H = 3	!H
  integer,parameter::KROME_idx_H2 = 4	!H2
  integer,parameter::KROME_idx_D = 5	!D
  integer,parameter::KROME_idx_HD = 6	!HD
  integer,parameter::KROME_idx_HE = 7	!HE
  integer,parameter::KROME_idx_Hj = 8	!H+
  integer,parameter::KROME_idx_H2j = 9	!H2+
  integer,parameter::KROME_idx_Dj = 10	!D+
  integer,parameter::KROME_idx_HEHj = 11	!HEH+
  integer,parameter::KROME_idx_CR = 12	!CR
  integer,parameter::KROME_idx_g = 13	!g
  integer,parameter::KROME_idx_Tgas = 14	!Tgas
  integer,parameter::KROME_idx_dummy = 15	!dummy

  integer,parameter::krome_idx_cool_h2 = 1
  integer,parameter::krome_idx_cool_h2gp = 2
  integer,parameter::krome_idx_cool_atomic = 3
  integer,parameter::krome_idx_cool_cen = 3
  integer,parameter::krome_idx_cool_hd = 4
  integer,parameter::krome_idx_cool_z = 5
  integer,parameter::krome_idx_cool_metal = 5
  integer,parameter::krome_idx_cool_dh = 6
  integer,parameter::krome_idx_cool_enthalpic = 6
  integer,parameter::krome_idx_cool_dust = 7
  integer,parameter::krome_idx_cool_compton = 8
  integer,parameter::krome_idx_cool_cie = 9
  integer,parameter::krome_idx_cool_continuum = 10
  integer,parameter::krome_idx_cool_cont = 10
  integer,parameter::krome_idx_cool_exp = 11
  integer,parameter::krome_idx_cool_expansion = 11
  integer,parameter::krome_idx_cool_ff = 12
  integer,parameter::krome_idx_cool_bss = 12
  integer,parameter::krome_idx_cool_custom = 13
  integer,parameter::krome_idx_cool_co = 14
  integer,parameter::krome_idx_cool_zcie = 15
  integer,parameter::krome_idx_cool_zcienouv = 16
  integer,parameter::krome_idx_cool_zextend = 17
  integer,parameter::krome_idx_cool_gh = 18
  integer,parameter::krome_idx_cool_oh = 19
  integer,parameter::krome_idx_cool_h2o = 20
  integer,parameter::krome_idx_cool_hcn = 21
  integer,parameter::krome_ncools = 21

  integer,parameter::krome_idx_heat_chem = 1
  integer,parameter::krome_idx_heat_compress = 2
  integer,parameter::krome_idx_heat_compr = 2
  integer,parameter::krome_idx_heat_photo = 3
  integer,parameter::krome_idx_heat_dh = 4
  integer,parameter::krome_idx_heat_enthalpic = 4
  integer,parameter::krome_idx_heat_photoav = 5
  integer,parameter::krome_idx_heat_av = 5
  integer,parameter::krome_idx_heat_cr = 6
  integer,parameter::krome_idx_heat_dust = 7
  integer,parameter::krome_idx_heat_xray = 8
  integer,parameter::krome_idx_heat_visc = 9
  integer,parameter::krome_idx_heat_viscous = 9
  integer,parameter::krome_idx_heat_custom = 10
  integer,parameter::krome_idx_heat_zcie = 11
  integer,parameter::krome_nheats = 11

  integer,parameter::krome_nrea=25
  integer,parameter::krome_nmols=11
  integer,parameter::krome_nspec=15
  integer,parameter::krome_natoms=4
  integer,parameter::krome_ndust=0
  integer,parameter::krome_ndustTypes=0
  integer,parameter::krome_nPhotoBins=0
  integer,parameter::krome_nPhotoRates=0

  real*8,parameter::krome_boltzmann_eV = 8.617332478d-5 !eV / K
  real*8,parameter::krome_boltzmann_J = 1.380648d-23 !J / K
  real*8,parameter::krome_boltzmann_erg = 1.380648d-16 !erg / K
  real*8,parameter::krome_iboltzmann_eV = 1d0/krome_boltzmann_eV !K / eV
  real*8,parameter::krome_iboltzmann_erg = 1d0/krome_boltzmann_erg !K / erg
  real*8,parameter::krome_planck_eV = 4.135667516d-15 !eV s
  real*8,parameter::krome_planck_J = 6.62606957d-34 !J s
  real*8,parameter::krome_planck_erg = 6.62606957d-27 !erg s
  real*8,parameter::krome_iplanck_eV = 1d0/krome_planck_eV !1 / eV / s
  real*8,parameter::krome_iplanck_J = 1d0/krome_planck_J !1 / J / s
  real*8,parameter::krome_iplanck_erg = 1d0/krome_planck_erg !1 / erg / s
  real*8,parameter::krome_gravity = 6.674d-8 !cm3 / g / s2
  real*8,parameter::krome_e_mass = 9.10938188d-28 !g
  real*8,parameter::krome_p_mass = 1.67262158d-24 !g
  real*8,parameter::krome_n_mass = 1.674920d-24 !g
  real*8,parameter::krome_ip_mass = 1d0/krome_p_mass !1/g
  real*8,parameter::krome_clight = 2.99792458e10 !cm/s
  real*8,parameter::krome_pi = 3.14159265359d0 !#
  real*8,parameter::krome_eV_to_erg = 1.60217646d-12 !eV -> erg
  real*8,parameter::krome_ry_to_eV = 13.60569d0 !rydberg -> eV
  real*8,parameter::krome_ry_to_erg = 2.179872d-11 !rydberg -> erg
  real*8,parameter::krome_seconds_per_year = 365d0*24d0*3600d0 !yr -> s
  real*8,parameter::krome_km_to_cm = 1d5 !km -> cm
  real*8,parameter::krome_cm_to_Mpc = 1.d0/3.08d24 !cm -> Mpc
  real*8,parameter::krome_kvgas_erg = 8.d0*krome_boltzmann_erg/krome_pi/krome_p_mass !
  real*8,parameter::krome_pre_kvgas_sqrt = sqrt(8.d0*krome_boltzmann_erg/krome_pi) !
  real*8,parameter::krome_pre_planck = 2.d0*krome_planck_erg/krome_clight**2 !erg/cm2*s3
  real*8,parameter::krome_exp_planck = krome_planck_erg / krome_boltzmann_erg !s*K
  real*8,parameter::krome_stefboltz_erg = 5.670373d-5 !erg/s/cm2/K4
  real*8,parameter::krome_N_avogadro = 6.0221d23 !#
  real*8,parameter::krome_Rgas_J = 8.3144621d0 !J/K/mol
  real*8,parameter::krome_Rgas_kJ = 8.3144621d-3 !kJ/K/mol
  real*8,parameter::krome_hubble = 0.704d0 !dimensionless
  real*8,parameter::krome_Omega0 = 1.0d0 !dimensionless
  real*8,parameter::krome_Omegab = 0.0456d0 !dimensionless
  real*8,parameter::krome_Hubble0 = 1.d2*krome_hubble*krome_km_to_cm*krome_cm_to_Mpc !1/s

contains

  !*******************
  subroutine krome_set_user_Trad(argset)
    use krome_commons
    implicit none
    real*8 :: argset
    user_Trad = argset
  end subroutine krome_set_user_Trad

  !*******************
  function krome_get_user_Trad()
    use krome_commons
    implicit none
    real*8 :: krome_get_user_Trad
    krome_get_user_Trad = user_Trad
  end function krome_get_user_Trad

  !************************
  !returns the Tdust averaged over the number density
  ! as computed in the tables
  function krome_get_table_Tdust(x,Tgas)
    use krome_commons
    use krome_grfuncs
    implicit none
    real*8 :: Tgas
    real*8 :: x(nmols), krome_get_table_Tdust
    real*8::n(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_get_table_Tdust = get_table_Tdust(n(:))

  end function krome_get_table_Tdust

  !**********************
  !convert from MOCASSIN abundances to KROME
  ! xmoc(i,j): MOCASSIN matrix (note: cm-3, real*4)
  !  i=species, j=ionization level
  ! imap: matrix position index map, integer
  ! returns KROME abundances (cm-3, real*8)
  function krome_convert_xmoc(xmoc,imap) result(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*4,intent(in):: xmoc(:,:)
    real*8::x(nmols),n(nspec)
    integer,intent(in)::imap(:)

    x(:) = 0d0

    x(idx_H) = xmoc(imap(1), 1)
    x(idx_D) = xmoc(imap(1), 1)
    x(idx_HE) = xmoc(imap(2), 1)
    x(idx_Hj) = xmoc(imap(1), 2)
    x(idx_Dj) = xmoc(imap(1), 2)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    x(idx_e) = get_electrons(n(:))

  end function krome_convert_xmoc

  !*************************
  !convert from KROME abundances to MOCASSIN
  ! x: KROME abuances (cm-3, real*8)
  ! imap: matrix position index map, integer
  ! xmoc(i,j): MOCASSIN matrix (note: cm-3, real*4)
  !  i=species, j=ionization level
  subroutine krome_return_xmoc(x,imap,xmoc)
    use krome_commons
    implicit none
    real*8,intent(in)::x(nmols)
    real*4,intent(out)::xmoc(:,:)
    integer,intent(in)::imap(:)

    xmoc(:,:) = 0d0

    xmoc(imap(1), 1) = x(idx_H)
    xmoc(imap(1), 1) = x(idx_D)
    xmoc(imap(2), 1) = x(idx_HE)
    xmoc(imap(1), 2) = x(idx_Hj)
    xmoc(imap(1), 2) = x(idx_Dj)

  end subroutine krome_return_xmoc

  !**********************
  !convert number density (cm-3) into column
  ! density (cm-2) using the specific density
  ! column method (see help for option
  ! -columnDensityMethod)
  ! num is the number density, x(:) is the species
  ! array, Tgas is the gas temperature
  ! If the method is not JEANS, x(:) and Tgas
  ! are dummy variables
  function krome_num2col(num,x,Tgas)
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8 :: x(nmols),krome_num2col
    real*8 :: Tgas,num
    real*8::n(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_num2col = num2col(num,n(:))

  end function krome_num2col

  !***********************
  !print on screen the current values of all phys variables
  subroutine krome_print_phys_variables()
    use krome_commons
    implicit none

    print *, "Tcmb:", phys_Tcmb
    print *, "zredshift:", phys_zredshift
    print *, "orthoParaRatio:", phys_orthoParaRatio
    print *, "metallicity:", phys_metallicity
    print *, "Tfloor:", phys_Tfloor

  end subroutine krome_print_phys_variables

  !*******************
  subroutine krome_set_Tcmb(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_Tcmb = arg
  end subroutine krome_set_Tcmb

  !*******************
  function krome_get_Tcmb()
    use krome_commons
    implicit none
    real*8 :: krome_get_Tcmb
    krome_get_Tcmb = phys_Tcmb
  end function krome_get_Tcmb

  !*******************
  subroutine krome_set_zredshift(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_zredshift = arg
  end subroutine krome_set_zredshift

  !*******************
  function krome_get_zredshift()
    use krome_commons
    implicit none
    real*8 :: krome_get_zredshift
    krome_get_zredshift = phys_zredshift
  end function krome_get_zredshift

  !*******************
  subroutine krome_set_orthoParaRatio(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_orthoParaRatio = arg
  end subroutine krome_set_orthoParaRatio

  !*******************
  function krome_get_orthoParaRatio()
    use krome_commons
    implicit none
    real*8 :: krome_get_orthoParaRatio
    krome_get_orthoParaRatio = phys_orthoParaRatio
  end function krome_get_orthoParaRatio

  !*******************
  subroutine krome_set_metallicity(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_metallicity = arg
  end subroutine krome_set_metallicity

  !*******************
  function krome_get_metallicity()
    use krome_commons
    implicit none
    real*8 :: krome_get_metallicity
    krome_get_metallicity = phys_metallicity
  end function krome_get_metallicity

  !*******************
  subroutine krome_set_Tfloor(arg)
    use krome_commons
    implicit none
    real*8 :: arg
    phys_Tfloor = arg
  end subroutine krome_set_Tfloor

  !*******************
  function krome_get_Tfloor()
    use krome_commons
    implicit none
    real*8 :: krome_get_Tfloor
    krome_get_Tfloor = phys_Tfloor
  end function krome_get_Tfloor

  !*****************************
  !dump the data for restart (UNDER DEVELOPEMENT!)
  !arguments: the species array and the gas temperature
  subroutine krome_store(x,Tgas,dt)
    use krome_commons
    implicit none
    integer::nfile,i
    real*8 :: x(nmols)
    real*8 :: Tgas,dt

    nfile = 92

    open(nfile,file="krome_dump.dat",status="replace")
    !dump temperature
    write(nfile,*) Tgas
    write(nfile,*) dt
    !dump species
    do i=1,nmols
      write(nfile,*) x(i)
    end do
    close(nfile)

  end subroutine krome_store

  !*****************************
  !restore the data from a dump (UNDER DEVELOPEMENT!)
  !arguments: the species array and the gas temperature
  subroutine krome_restore(x,Tgas,dt)
    use krome_commons
    implicit none
    integer::nfile,i
    real*8 :: x(nmols)
    real*8 :: Tgas,dt

    nfile = 92

    open(nfile,file="krome_dump.dat",status="old")
    !restore temperature
    read(nfile,*) Tgas
    read(nfile,*) dt
    !restore species
    do i=1,nmols
      read(nfile,*) x(i)
    end do
    close(nfile)

  end subroutine krome_restore

  !****************************
  !switch on the thermal calculation
  subroutine krome_thermo_on()
    use krome_commons
    krome_thermo_toggle = 1
  end subroutine krome_thermo_on

  !****************************
  !switch off the thermal calculation
  subroutine krome_thermo_off()
    use krome_commons
    krome_thermo_toggle = 0
  end subroutine krome_thermo_off

  !***************************
  !alias for coe in krome_subs
  ! returns the coefficient array of size krome_nrea
  ! for a given Tgas
  function krome_get_coef(Tgas,x)
    use krome_commons
    use krome_subs
    use krome_tabs
    real*8 :: krome_get_coef(nrea),x(nmols)
    real*8,value:: Tgas
    real*8::n(nspec)
    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas

    krome_get_coef(:) = coe(n(:))

  end function krome_get_coef

  !****************************
  !get the mean molecular weight from
  ! mass fractions
  function krome_get_mu_x(xin)
    use krome_commons
    implicit none
    real*8 :: xin(nmols), krome_get_mu_x
    real*8::n(nmols)
    n(:) = krome_x2n(xin(:),1d0)
    krome_get_mu_x = krome_get_mu(n(:))
  end function krome_get_mu_x

  !****************************
  !return the adiabatic index from mass fractions
  ! and temperature in K
  function krome_get_gamma_x(xin,inTgas)
    use krome_commons
    implicit none
    real*8 :: inTgas
    real*8 :: xin(nmols), krome_get_gamma_x
    real*8::x(nmols),Tgas,rhogas

    Tgas = inTgas
    x(:) = krome_x2n(xin(:),1d0)
    krome_get_gamma_x = krome_get_gamma(x(:),Tgas)

  end function krome_get_gamma_x

  !***************************
  !normalize mass fractions and
  ! set charge to zero
  subroutine krome_consistent_x(x)
    use krome_commons
    use krome_constants
    implicit none
    real*8 :: x(nmols)
    real*8::isumx,sumx,xerr,imass(nmols),ee

    !1. charge consistency
    imass(:) = krome_get_imass()

    x(idx_e) = 0.d0

    ee = sum(krome_get_charges()*x(:)*imass(:))
    ee = max(ee*e_mass,0d0)
    x(idx_e) = ee

    !2. mass fraction consistency
    sumx = sum(x)

    !NOTE: uncomment here if you want some additional control
    !conservation error threshold: rise an error if above xerr
    !xerr = 1d-2
    !if(abs(sum-1d0)>xerr) then
    !   print *,"ERROR: some problem with conservation!"
    !   print *,"|sum(x)-1|=",abs(sum-1d0)
    !   stop
    !end if

    isumx = 1d0/sumx
    x(:) = x(:) * isumx

  end subroutine krome_consistent_x

  !*********************
  !return an array sized krome_nmols containing
  ! the mass fractions (#), computed from the number
  ! densities (1/cm3) and the total density in g/cm3
  function krome_n2x(n,rhogas)
    use krome_commons
    implicit none
    real*8 :: n(nmols),krome_n2x(nmols)
    real*8,value :: rhogas

    krome_n2x(:) = n(:) * krome_get_mass() / rhogas

  end function krome_n2x

  !********************
  !return an array sized krome_nmols containing
  ! the number densities (1/cm3), computed from the mass
  ! fractions and the total density in g/cm3
  function krome_x2n(x,rhogas)
    use krome_commons
    implicit none
    real*8 :: x(nmols),krome_x2n(nmols)
    real*8,value :: rhogas

    !compute densities from fractions
    krome_x2n(:) = rhogas * x(:) * krome_get_imass()

  end function krome_x2n

  !******************
  !returns free-fall time using the number density
  ! abundances of array x(:)
  function krome_get_free_fall_time(x)
    use krome_commons
    use krome_getphys
    implicit none
    real*8::krome_get_free_fall_time
    real*8::x(:),n(nspec)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    krome_get_free_fall_time = get_free_fall_time(n(:))

  end function krome_get_free_fall_time

  !******************
  !returns free-fall time using the total mass density
  !  of gas, rhogas (g/cm3)
  function krome_get_free_fall_time_rho(rhogas)
    use krome_getphys
    implicit none
    real*8::krome_get_free_fall_time_rho
    real*8::rhogas

    krome_get_free_fall_time_rho = get_free_fall_time_rho(rhogas)

  end function krome_get_free_fall_time_rho

  !*******************
  !do only cooling and heating
  subroutine krome_thermo(x,Tgas,dt)
    use krome_commons
    use krome_cooling
    use krome_heating
    use krome_subs
    use krome_tabs
    use krome_constants
    use krome_gadiab
    implicit none
    real*8 :: x(nmols)
    real*8 :: Tgas,dt
    real*8::n(nspec),nH2dust,dTgas,k(nrea),krome_gamma

    nH2dust = 0d0
    n(:) = 0d0
    n(idx_Tgas) = Tgas
    n(1:nmols) = x(:)
    k(:) = coe_tab(n(:)) !compute coefficients
    krome_gamma = gamma_index(n(:))

    dTgas = (heating(n(:), Tgas, k(:), nH2dust) - cooling(n(:), Tgas)) &
        * (krome_gamma - 1.d0) / boltzmann_erg / sum(n(1:nmols))

    Tgas = Tgas + dTgas*dt !update gas

  end subroutine krome_thermo

  !*************************
  !get heating (erg/cm3/s) for a given species
  ! array x(:) and Tgas
  function krome_get_heating(x,inTgas)
    use krome_heating
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: inTgas
    real*8 :: x(nmols), krome_get_heating
    real*8::Tgas,k(nrea),nH2dust,n(nspec)
    n(1:nmols) = x(:)
    Tgas = inTgas
    n(idx_Tgas) = Tgas
    k(:) = coe(n(:))
    nH2dust = 0d0
    krome_get_heating = heating(n(:),Tgas,k(:),nH2dust)
  end function krome_get_heating

  !*****************************
  ! get an array containing individual heatings (erg/cm3/s)
  ! the array has size krome_nheats. see heatcool.gps
  ! for index list
  function krome_get_heating_array(x,inTgas)
    use krome_heating
    use krome_subs
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas,k(nrea),nH2dust
    real*8 :: x(nmols),krome_get_heating_array(nheats)
    real*8,value :: inTgas

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = inTgas
    !#KROME_Tdust_copy
    k(:) = coe(n(:))
    Tgas = inTgas
    nH2dust = 0d0
    krome_get_heating_array(:) = get_heating_array(n(:),Tgas,k(:),nH2dust)

  end function krome_get_heating_array

  !*************************
  !get cooling (erg/cm3/s) for x(:) species array
  ! and Tgas
  function krome_get_cooling(x,inTgas)
    use krome_cooling
    use krome_commons
    implicit none
    real*8 :: inTgas
    real*8 :: x(nmols), krome_get_cooling
    real*8::Tgas,n(nspec)
    n(1:nmols) = x(:)
    Tgas = inTgas
    n(idx_Tgas) = Tgas
    krome_get_cooling = cooling(n,Tgas)
  end function krome_get_cooling

  !*****************************
  ! get an array containing individual coolings (erg/cm3/s)
  ! the array has size krome_ncools. see heatcool.gps
  ! for index list
  function krome_get_cooling_array(x,inTgas)
    use krome_cooling
    use krome_commons
    implicit none
    real*8::n(nspec),Tgas
    real*8 :: x(nmols),krome_get_cooling_array(ncools)
    real*8,value :: inTgas

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = inTgas
    !#KROME_Tdust_copy
    Tgas = inTgas
    krome_get_cooling_array(:) = get_cooling_array(n(:),Tgas)

  end function krome_get_cooling_array

  !******************
  !alias of plot_cool
  subroutine krome_plot_cooling(n)
    use krome_cooling
    implicit none
    real*8 :: n(krome_nmols)

    call plot_cool(n(:))

  end subroutine krome_plot_cooling

  !****************
  !alias for dumping cooling in the unit nfile_in
  subroutine krome_dump_cooling(n,Tgas,nfile_in)
    use krome_cooling
    use krome_commons
    implicit none
    real*8 :: n(nmols)
    real*8 :: Tgas
    real*8::x(nspec)
    integer, optional :: nfile_in
    integer::nfile
    nfile = 31
    x(:) = 0.d0
    x(1:nmols) = n(:)
    if(present(nfile_in)) nfile = nfile_in
    call dump_cool(x(:),Tgas,nfile)

  end subroutine krome_dump_cooling

  !************************
  !conserve the total amount of nucleii,
  ! alias for conserveLin_x in subs
  subroutine krome_conserveLin_x(x,ref)
    use krome_commons
    use krome_subs
    implicit none
    real*8 :: x(nmols),ref(natoms)

    call conserveLin_x(x(:),ref(:))

  end subroutine krome_conserveLin_x

  !************************
  !conserve the total amount of nucleii,
  ! alias for conserveLin_x in subs
  function krome_conserveLinGetRef_x(x)
    use krome_commons
    use krome_subs
    implicit none
    real*8 :: x(nmols),krome_conserveLinGetRef_x(natoms)

    krome_conserveLinGetRef_x(:) = &
        conserveLinGetRef_x(x(:))

  end function krome_conserveLinGetRef_x

  !*************************
  !force conservation to array x(:)
  !using xi(:) as initial abundances.
  !alias for conserve in krome_subs
  function krome_conserve(x,xi)
    use krome_subs
    implicit none
    real*8 :: x(krome_nmols),xi(krome_nmols),krome_conserve(krome_nmols)
    real*8::n(krome_nspec),ni(krome_nspec)

    n(:) = 0d0
    ni(:) = 0d0
    n(1:krome_nmols) = x(1:krome_nmols)
    ni(1:krome_nmols) = xi(1:krome_nmols)
    n(:) = conserve(n(:), ni(:))
    krome_conserve(:) = n(1:krome_nmols)

  end function krome_conserve

  !***************************
  !get the adiabatic index for x(:) species abundances
  ! and Tgas.
  ! alias for gamma_index in krome_subs
  function krome_get_gamma(x,Tgas)
    use krome_subs
    use krome_commons
    use krome_gadiab
    real*8 :: Tgas
    real*8 :: x(nmols), krome_get_gamma
    real*8::n(nspec)
    n(:) = 0.d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    krome_get_gamma = gamma_index(n(:))
  end function krome_get_gamma

  !***************************
  !get an integer array containing the atomic numbers Z
  ! of the spcecies.
  ! alias for get_zatoms
  function krome_get_zatoms()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    integer :: krome_get_zatoms(nmols)
    integer::zatoms(nspec)

    zatoms(:) = get_zatoms()
    krome_get_zatoms(:) = zatoms(1:nmols)

  end function krome_get_zatoms

  !****************************
  !get the mean molecular weight from
  ! number density and mass density.
  ! alias for get_mu in krome_subs module
  function krome_get_mu(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8 :: x(nmols), krome_get_mu
    real*8::n(1:nspec)
    n(:) = 0d0
    n(1:nmols) = x(:)
    krome_get_mu = get_mu(n(:))
  end function krome_get_mu

  !***************************
  !get the names of the reactions as a
  ! character*50 array of krome_nrea
  ! elements
  !! !! cannot yet be called from C
  function krome_get_rnames()
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    character*50 :: krome_get_rnames(nrea)

    krome_get_rnames(:) = get_rnames()

  end function krome_get_rnames

  !*****************
  !get an array of double containing the masses in g
  ! of the species.
  ! alias for get_mass in krome_subs
  function krome_get_mass()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_mass(nmols)
    tmp(:) = get_mass()
    krome_get_mass = tmp(1:nmols)
  end function krome_get_mass

  !*****************
  !get an array of double containing the inverse
  ! of the mass (1/g) of the species
  !alias for get_imass in krome_subs
  function krome_get_imass()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_imass(nmols)
    tmp(:) = get_imass()
    krome_get_imass = tmp(1:nmols)
  end function krome_get_imass

  !***********************
  !get the total number of H nuclei
  function krome_get_Hnuclei(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    real*8::n(nspec)
    real*8 :: krome_get_Hnuclei, x(nmols)
    n(:) = 0d0
    n(1:nmols) = x(:)

    krome_get_Hnuclei = get_Hnuclei(n(:))

  end function krome_get_Hnuclei

  !*****************
  !get an array of size krome_nmols containing the
  ! charges of the species.
  ! alias for get_charges
  function krome_get_charges()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    real*8::tmp(nspec)
    real*8 :: krome_get_charges(nmols)
    tmp(:) = get_charges()
    krome_get_charges = tmp(1:nmols)
  end function krome_get_charges

  !*****************
  !get an array of character*16 and size krome_nmols
  ! containing the names of all the species.
  ! alias for get_names
  !!  !! cannot yet be called from C
  function krome_get_names()
    use krome_subs
    use krome_commons
    use krome_getphys
    implicit none
    character*16 :: krome_get_names(nmols)
    character*16::tmp(nspec)
    tmp(:) = get_names()
    krome_get_names = tmp(1:nmols)
  end function krome_get_names

  !********************
  !get space-separated header of chemical species
  function krome_get_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    character*33::krome_get_names_header
    character*16::tmp(nspec)
    integer::i

    tmp(:) = get_names()

    krome_get_names_header = ""
    do i=1,nmols
      krome_get_names_header = trim(krome_get_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_names_header

  !********************
  !get space-separated header of coolings
  function krome_get_cooling_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    character*141::krome_get_cooling_names_header
    character*16::tmp(ncools)
    integer::i

    tmp(:) = get_cooling_names()

    krome_get_cooling_names_header = ""
    do i=1,ncools
      if(trim(tmp(i))=="") cycle
      krome_get_cooling_names_header = trim(krome_get_cooling_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_cooling_names_header

  !********************
  !get space-separated header of heatings
  function krome_get_heating_names_header()
    use krome_commons
    use krome_getphys
    implicit none
    character*87::krome_get_heating_names_header
    character*16::tmp(nheats)
    integer::i

    tmp(:) = get_heating_names()

    krome_get_heating_names_header = ""
    do i=1,nheats
      if(trim(tmp(i))=="") cycle
      krome_get_heating_names_header = trim(krome_get_heating_names_header)//" "//trim(tmp(i))
    end do

  end function krome_get_heating_names_header

  !*****************
  !get the index of the species with name name.
  ! alias for get_index
  !! !! cannot yet be called from C
  function krome_get_index(name)
    use krome_subs
    implicit none
    integer :: krome_get_index
    character*(*) :: name
    krome_get_index = get_index(name)
  end function krome_get_index

  !*******************
  !get the total density of the gas in g/cm3
  ! giving all the number densities n(:)
  function krome_get_rho(n)
    use krome_commons
    real*8 :: krome_get_rho, n(nmols)
    real*8::m(nmols)
    m(:) = krome_get_mass()
    krome_get_rho = sum(m(:)*n(:))
  end function krome_get_rho

  !*************************
  !scale the abundances of the metals contained in n(:)
  ! to Z according to Asplund+2009.
  ! note that this applies only to neutral atoms.
  subroutine krome_scale_Z(x,Z)
    use krome_commons
    use krome_getphys
    real*8 :: x(nmols)
    real*8 :: Z
    real*8::Htot,n(nspec)

    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0

    Htot = get_Hnuclei(n(:))

  end subroutine krome_scale_Z

  !*************************
  !set the total metallicity
  ! in terms of Z/Z_solar
  subroutine krome_set_Z(xarg)
    use krome_commons
    real*8 :: xarg

    total_Z = xarg

  end subroutine krome_set_Z

  !*************************
  !set D is in terms of D_solar (D/D_sol).
  subroutine krome_set_dust_to_gas(xarg)
    use krome_commons
    real*8 :: xarg

    dust2gas_ratio = xarg

  end subroutine

  !*************************
  !set the clumping factor
  subroutine krome_set_clump(xarg)
    use krome_commons
    real*8 :: xarg

    clump_factor = xarg

  end subroutine krome_set_clump

  !***********************
  !get the number of electrons assuming
  ! total neutral charge (cations-anions)
  function krome_get_electrons(x)
    use krome_commons
    use krome_subs
    use krome_getphys
    real*8 :: x(nmols), krome_get_electrons
    real*8::n(nspec)
    n(1:nmols) = x(:)
    n(nmols+1:nspec) = 0d0
    krome_get_electrons = get_electrons(n(:))
  end function krome_get_electrons

  !**********************
  !print on screen the first nbest highest reaction fluxes
  subroutine krome_print_best_flux(xin,Tgas,nbest)
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: xin(nmols)
    real*8 :: Tgas
    real*8::x(nmols),n(nspec)
    integer :: nbest
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux(n,Tgas,nbest)

  end subroutine krome_print_best_flux

  !*********************
  !print only the highest fluxes greater than a fraction frac
  ! of the maximum flux
  subroutine krome_print_best_flux_frac(xin,Tgas,frac)
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: xin(nmols)
    real*8 :: Tgas,frac
    real*8::n(nspec)
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux_frac(n,Tgas,frac)

  end subroutine krome_print_best_flux_frac

  !**********************
  !print the highest nbest fluxes for reactions involving
  !a given species using the index idx_find (e.g. krome_idx_H2)
  subroutine krome_print_best_flux_spec(xin,Tgas,nbest,idx_find)
    use krome_subs
    use krome_commons
    implicit none
    real*8 :: xin(nmols)
    real*8 :: Tgas
    real*8::n(nspec)
    integer :: nbest,idx_find
    n(1:nmols) = xin(:)
    n(idx_Tgas) = Tgas
    call print_best_flux_spec(n,Tgas,nbest,idx_find)
  end subroutine krome_print_best_flux_spec

  !*******************************
  !get an array of size krome_nrea with
  ! the fluxes of all the reactions in cm-3/s
  function krome_get_flux(n,Tgas)
    use krome_commons
    use krome_subs
    real*8 :: krome_get_flux(nrea),n(nmols)
    real*8,value :: Tgas
    real*8::x(nspec)
    x(:) = 0.d0
    x(1:nmols) = n(:)
    x(idx_Tgas) = Tgas
    krome_get_flux(:) = get_flux(x(:), Tgas)
  end function krome_get_flux

  !*****************************
  !store the fluxes to the file unit ifile
  ! using the chemical composition x(:), and the
  ! gas temperature Tgas. xvar is th value of an
  ! user-defined independent variable that
  ! can be employed for plots.
  ! the file columns are as follow
  ! rate number, xvar, absolute flux,
  !  flux/maxflux, flux fraction wrt total,
  !  reaction name (*50 string)
  subroutine krome_explore_flux(x,Tgas,ifile,xvar)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8 :: x(nmols)
    real*8 :: Tgas,xvar
    real*8::flux(nrea),fluxmax,sumflux,n(nspec)
    integer :: ifile
    integer::i
    character*50::rname(nrea)

    !get reaction names
    rname(:) = get_rnames()
    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = Tgas
    !get fluxes
    flux(:) = get_flux(n(:), Tgas)
    fluxmax = maxval(flux) !maximum flux
    sumflux = sum(flux) !sum of all the fluxes
    !loop on reactions
    do i=1,nrea
      write(ifile,'(I8,5E17.8e3,a3,a50)') i,xvar,Tgas,flux(i),&
          flux(i)/fluxmax, flux(i)/sumflux," ",rname(i)
    end do
    write(ifile,*)

  end subroutine krome_explore_flux

  !*********************
  !get nulcear qeff for the reactions
  function krome_get_qeff()
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    real*8 :: krome_get_qeff(nrea)

    krome_get_qeff(:) = get_qeff()

  end function krome_get_qeff

  !************************
  !dump the fluxes to the file unit nfile
  subroutine krome_dump_flux(n,Tgas,nfile)
    use krome_commons
    real*8 :: n(nmols)
    real*8 :: Tgas
    real*8::flux(nrea)
    integer :: nfile
    integer::i

    flux(:) = krome_get_flux(n(:),Tgas)
    do i=1,nrea
      write(nfile,'(I8,E17.8e3)') i,flux(i)
    end do
    write(nfile,*)

  end subroutine krome_dump_flux

  !************************
  !dump all the evaluation of the coefficient rates in
  ! the file funit, in the range inTmin, inTmax, using
  ! imax points
  subroutine krome_dump_rates(inTmin,inTmax,imax,funit)
    use krome_commons
    use krome_subs
    implicit none
    integer::i,j
    integer :: funit,imax
    real*8 :: inTmin,inTmax
    real*8::Tmin,Tmax,Tgas,k(nrea),n(nspec)

    Tmin = log10(inTmin)
    Tmax = log10(inTmax)

    n(:) = 1d-40
    do i=1,imax
      Tgas = 1d1**((i-1)*(Tmax-Tmin)/(imax-1)+Tmin)
      n(idx_Tgas) = Tgas
      k(:) = coe(n(:))
      do j=1,nrea
        write(funit,'(E17.8e3,I8,E17.8e3)') Tgas,j,k(j)
      end do
      write(funit,*)
    end do

  end subroutine krome_dump_rates

  !************************
  !print species informations on screen
  subroutine krome_get_info(x, Tgas)
    use krome_commons
    use krome_subs
    use krome_getphys
    implicit none
    integer::i,charges(nspec)
    real*8 :: x(nmols)
    real*8 :: Tgas
    real*8::masses(nspec)
    character*16::names(nspec)

    names(:) = get_names()
    charges(:) = get_charges()
    masses(:) = get_mass()

    print '(a4,a10,a11,a5,a11)',"#","Name","m (g)","Chrg","x"
    do i=1,size(x)
      print '(I4,a10,E11.3,I5,E11.3)',i," "//names(i),masses(i),charges(i),x(i)
    end do
    print '(a30,E11.3)'," sum",sum(x)

    print '(a14,E11.3)',"Tgas",Tgas
  end subroutine krome_get_info

  !*****************************
  subroutine krome_set_mpi_rank(xarg)
    use krome_commons
    implicit none
    integer :: xarg
    krome_mpi_rank=xarg
  end subroutine krome_set_mpi_rank

  !**************************
  function krome_get_jacobian(j,x,Tgas)
    use krome_ode
    use krome_commons
    implicit none
    integer, value :: j
    real*8,value :: Tgas
    real*8 :: x(nmols),krome_get_jacobian(nspec)
    integer::ian, jan, i
    real*8::tt, n(nspec)
    real*8::pdj(nspec)

    n(:) = 0d0
    n(1:nmols) = x(:)
    n(idx_Tgas) = tgas

    tt = 0d0
    ian = 0
    jan = 0

    call jes(nspec, tt, n, j, ian, jan, pdj)
    krome_get_jacobian(:) = pdj(:)

  end function krome_get_jacobian

end module krome_user

!############### MODULE ##############
module krome_reduction
contains

  !**************************
  function fex_check(n,Tgas)
    use krome_commons
    use krome_tabs
    implicit none
    integer::i
    integer::r1,r2,r3
    real*8::fex_check,n(nspec),k(nrea),rrmax,Tgas

    k(:) = coe_tab(n(:))
    rrmax = 0.d0
    n(idx_dummy) = 1.d0
    n(idx_g) = 1.d0
    n(idx_CR) = 1.d0
    do i=1,nrea
      r1 = arr_r1(i)
      r2 = arr_r2(i)
      r3 = arr_r3(i)
      arr_flux(i) = k(i)*n(r1)*n(r2)*n(r3)
      rrmax = max(rrmax, arr_flux(i))
    end do
    fex_check = rrmax

  end function fex_check

end module krome_reduction

!############### MODULE ##############
module krome_main

  integer::krome_call_to_fex
  !$omp threadprivate(krome_call_to_fex)

contains

  ! *************************************************************
  !  This file has been generated with:
  !  KROME 14.08.dev on 2022-06-08 12:22:52
  !  Changeset d0f9dcf
  !  see http://kromepackage.org
  !
  !  Written and developed by Tommaso Grassi and Stefano Bovino
  !
  !  Contributors:
  !  J.Boulangier, T.Frostholm, D.Galli, F.A.Gianturco, T.Haugboelle,
  !  A.Lupi, J.Prieto, J.Ramsey, D.R.G.Schleicher, D.Seifried, E.Simoncini,
  !  E.Tognelli
  !  KROME is provided "as it is", without any warranty.
  ! *************************************************************

  !*******************************
  !KROME main (interface to the solver library)

  subroutine krome(x,Tgas,dt  )
    use krome_commons
    use krome_subs
    use krome_ode
    use krome_reduction
    use krome_dust
    use krome_getphys
    use krome_tabs
    implicit none
    real*8 :: Tgas,dt
    real*8 :: x(nmols)
    real*8 :: rhogas

    real*8::mass(nspec),n(nspec),tloc,xin
    real*8::rrmax,totmass,n_old(nspec),ni(nspec),invTdust(ndust)
    integer::icount,i,icount_max
    integer:: ierr

    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer::neq(1),itol,itask,istate,iopt,lrw,liw,mf
    integer::iwork(153)
    real*8::atol(nspec),rtol(nspec)
    real*8::rwork(760)
    logical::got_error,equil

    !****************************
    !init DLSODES (see DLSODES manual)
    call XSETF(0)!toggle solver verbosity
    got_error = .false.
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are scalar
    rtol(:) = 1.000000d-04 !relative tolerance
    atol(:) = 1.000000d-10 !absolute tolerance
    icount_max = 10000 !maximum number of iterations

    itask = 1
    iopt = 0

    !MF=
    !  = 222 internal-generated JAC and sparsity
    !  = 121 user-provided JAC and internal generated sparsity
    !  =  22 internal-generated JAC but sparsity user-provided
    !  =  21 user-provided JAC and sparsity
    MF = 222
    !end init DLSODES
    !****************************

    ierr = 0 !error flag, zero==OK!
    n(:) = 0d0 !initialize densities

    n(1:nmols) = x(:)

    n(idx_Tgas) = Tgas !put temperature in the input array

    icount = 0 !count solver iterations
    istate = 1 !init solver state
    tloc = 0.d0 !set starting time

    !store initial values
    ni(:) = n(:)
    n_global(:) = n(:)

    n_old(:) = -1d99
    krome_call_to_fex = 0
    do
      icount = icount + 1
      !solve ODE
      CALL DLSODES(fex, NEQ(:), n(:), tloc, dt, &
          ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, &
          LIW, JES, MF)

      krome_call_to_fex = krome_call_to_fex + IWORK(12)
      !check DLSODES exit status
      if(istate==2) then
        exit !sucsessful integration
      elseif(istate==-1) then
        istate = 1 !exceeded internal max iterations
      elseif(istate==-5 .or. istate==-4) then
        istate = 3 !wrong sparsity recompute
      elseif(istate==-3) then
        n(:) = ni(:)
        istate = 1
      else
        got_error = .true.
      end if

      if(got_error.or.icount>icount_max) then
        if (krome_mpi_rank>0) then
          print *,krome_mpi_rank,"ERROR: wrong solver exit status!"
          print *,krome_mpi_rank,"istate:",istate
          print *,krome_mpi_rank,"iter count:",icount
          print *,krome_mpi_rank,"max iter count:",icount_max
          print *,krome_mpi_rank,"SEE KROME_ERROR_REPORT file"
        else
          print *,"ERROR: wrong solver exit status!"
          print *,"istate:",istate
          print *,"iter count:",icount
          print *,"max iter count:",icount_max
          print *,"SEE KROME_ERROR_REPORT file"
        end if
        call krome_dump(n(:), rwork(:), iwork(:), ni(:))
        stop
      end if

    end do

    !avoid negative species
    do i=1,nspec
      n(i) = max(n(i),0d0)
    end do

    !returns to user array
    x(:) = n(1:nmols)

    Tgas = n(idx_Tgas) !get new temperature

  end subroutine krome

  !*********************************
  !integrates to equilibrium using constant temperature
  subroutine krome_equilibrium(x,Tgas,verbosity)
    use krome_ode
    use krome_subs
    use krome_commons
    use krome_constants
    use krome_getphys
    use krome_tabs
    implicit none
    integer::mf,liw,lrw,itol,meth,iopt,itask,istate,neq(1)
    integer::i,imax
    integer,optional::verbosity
    integer::verbose
    real*8 :: Tgas
    real*8 :: x(nmols)
    real*8 :: rhogas
    real*8::tloc,n(nspec),mass(nspec),ni(nspec)
    real*8::dt,xin
    integer::iwork(153)
    real*8::atol(nspec),rtol(nspec)
    real*8::rwork(760)
    real*8::ertol,eatol,max_time,t_tot,ntot_tol,err_species
    logical::converged

    integer, save :: ncall=0
    integer, parameter :: ncall_print_frequency=20000
    integer :: ncallp
    integer::charges(nspec)
    real*8::masses(nspec)
    character*16::names(nspec)

    !set verbosity from argument
    verbose = 1 !default is verbose
    if(present(verbosity)) verbose = verbosity

    call XSETF(0)!toggle solver verbosity
    meth = 2
    neq = nspec !number of eqns
    liw = size(iwork)
    lrw = size(rwork)
    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are scalar
    rtol(:) = 1d-6 !relative tolerance
    atol(:) = 1d-20 !absolute tolerance

    ! Switches to decide when equilibrium has been reached
    ertol = 1d-5  ! relative min change in a species
    eatol = 1d-12 ! absolute min change in a species
    max_time=seconds_per_year*5d8 ! max time we will be integrating for

    !for DLSODES options see its manual
    iopt = 0
    itask = 1
    istate = 1

    mf = 222 !internally evaluated sparsity and jacobian
    tloc = 0d0 !initial time

    n(:) = 0d0 !initialize densities
    !copy into array
    n(nmols+1:) = 0d0
    n(1:nmols) = x(:)

    n(idx_Tgas) = Tgas

    !store previous values
    ni(:) = n(:)
    n_global(:) = ni(:)

    imax = 1000

    dt = seconds_per_year * 1d2
    t_tot = dt
    converged = .false.
    do while (.not. converged)
      do i=1,imax
        !solve ODE
        CALL DLSODES(fcn_tconst, NEQ(:), n(:), tloc, dt, ITOL, RTOL, ATOL,&
            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, jcn_dummy, MF)
        if(istate==2) then
          exit
        else
          istate=1
        end if
      end do
      !check errors
      if(istate.ne.2) then
        print *,"ERROR: no equilibrium found!"
        stop
      end if

      !avoid negative species
      do i=1,nspec
        n(i) = max(n(i),0d0)
      end do

      ! check if we have converged by comparing the error in any species with an relative abundance above eatol
      converged = maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) .lt. ertol &
          .or. t_tot .gt. max_time

      ! Increase integration time by a reasonable factor
      if(.not. converged) then
        dt = dt * 3.
        t_tot = t_tot + dt
        ni = n
        n_global = n
      endif
    enddo
    !returns to user array
    x(:) = n(1:nmols)

    if(t_tot > max_time .and. &
        maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols)))) > 0.2 .and. verbose>0) then
    print *, 'krome_equilibrium: Did not converge in ', max_time / seconds_per_year, ' years.'
    print *, 'Tgas :', Tgas
    names(:) = get_names()
    charges(:) = get_charges()
    masses(:) = get_mass()

    print '(a4,a10,a11,a5,a16)',"#","Name","m (g)","Chrg","  Current / Last"
    do i=1,nmols
      print '(I4,a10,E11.3,I5,2E14.6,E11.3)',i," "//names(i),masses(i),charges(i),n(i),ni(i),abs(n(i) - ni(i)) / max(n(i),eatol*sum(n(1:nmols)))
    end do
    print '(a30,2E14.6)'," sum",sum(n(1:nmols)),sum(ni(1:nmols))
    print *, 'Fractional error :', maxval(abs(n(1:nmols) - ni(1:nmols)) / max(n(1:nmols),eatol*sum(n(1:nmols))))
    print *, 'Absolute and relative floors:', eatol, ertol
  end if

  ! Print info ever so often
  !$omp critical
  ncall=ncall+1
  ncallp = ncall
  !$omp end critical

  if(modulo(ncallp,ncall_print_frequency)==0 .and. verbose>0) then
    print *, 'Found equilibrium for ', ncallp, ' cells.'
  end if

end subroutine krome_equilibrium

!********************
!dummy jacobian
subroutine jcn_dummy()
  implicit none
end subroutine jcn_dummy

!*******************
!dn/dt where dT/dt=0
subroutine fcn_tconst(n,tt,x,f)
  use krome_commons
  use krome_ode
  implicit none
  integer::n,ierr
  real*8::x(n),f(n),tt
  call fex(n,tt,x(:),f(:))
  f(idx_Tgas) = 0d0
end subroutine fcn_tconst

!*******************************
subroutine krome_dump(n,rwork,iwork,ni)
  use krome_commons
  use krome_subs
  use krome_tabs
  use krome_reduction
  use krome_ode
  use krome_getphys
  integer::fnum,i,iwork(:),idx(nrea),j
  real*8::n(:),rwork(:),rrmax,k(nrea),kmax,rperc,kperc,dn(nspec),tt,ni(:)
  real*8::m(nspec)
  character*16::names(nspec),FMTi,FMTr
  character*50::rnames(nrea),fname,prex
  integer,save::mx_dump=1000 ! max nr of reports before terminating
  fnum = 99
  if (krome_mpi_rank>0) then
    write(fname,'(a,i5.5)') "KROME_ERROR_REPORT_",krome_mpi_rank
  else
    fname = "KROME_ERROR_REPORT"
  endif
  open(fnum,FILE=trim(fname),status="replace")
  tt = 0d0
  names(:) = get_names()
  rnames(:) = get_rnames()
  call fex(nspec,tt,n(:),dn(:))

  write(fnum,*) "KROME ERROR REPORT"
  write(fnum,*)
  !SPECIES
  write(fnum,*) "Species abundances"
  write(fnum,*) "**********************"
  write(fnum,'(a5,a20,3a12)') "#","name","qty","dn/dt","ninit"
  write(fnum,*) "**********************"
  do i=1,nspec
    write(fnum,'(I5,a20,3E12.3e3)') i,names(i),n(i),dn(i),ni(i)
  end do
  write(fnum,*) "**********************"

  !F90 FRIENDLY RESTART
  write(fnum,*)
  write(fnum,*) "**********************"
  write(fnum,*) "F90-friendly species"
  write(fnum,*) "**********************"
  do i=1,nspec
    write(prex,'(a,i3,a)') "x(",i,") = "
    write(fnum,*) trim(prex),ni(i),"!"//names(i)
  end do

  write(fnum,*) "**********************"

  !RATE COEFFIECIENTS
  k(:) = coe_tab(n(:))
  idx(:) = idx_sort(k(:))
  kmax = maxval(k)
  m(:) = get_mass()
  write(fnum,*)
  write(fnum,*) "Rate coefficients (sorted) at Tgas",n(idx_Tgas)
  write(fnum,*) "Initial density (g/cm**-3)", sum(m(1:nmols)*ni(1:nmols))
  write(fnum,*) "Density at icount = icount_max (g/cm**-3)", sum(m(1:nmols)*n(1:nmols))
  write(fnum,*) "**********************"
  write(fnum,'(a5,2a12,a10)') "#","k","k %","  name"
  write(fnum,*) "**********************"
  do j=1,nrea
    i = idx(j)
    kperc = 0.d0
    if(kmax>0.d0) kperc = k(i)*1d2/kmax
    write(fnum,'(I5,2E12.3e3,a2,a50)') i,k(i),kperc,"  ", rnames(i)
  end do
  write(fnum,*) "**********************"
  write(fnum,*)

  !FLUXES
  call load_arrays
  rrmax = fex_check(n(:), n(idx_Tgas))
  idx(:) = idx_sort(arr_flux(:))
  write(fnum,*)
  write(fnum,*) "Reaction magnitude (sorted) [k*n1*n2*n3*...]"
  write(fnum,*) "**********************"
  write(fnum,'(a5,2a12,a10)') "#","flux","flux %","  name"
  write(fnum,*) "**********************"
  do j=1,nrea
    i = idx(j)
    rperc = 0.d0
    if(rrmax>0.d0) rperc = arr_flux(i)*1d2/rrmax
    write(fnum,'(I5,2E12.3e3,a2,a50)') i,arr_flux(i),rperc,"  ",rnames(i)
  end do
  write(fnum,*) "**********************"
  write(fnum,*)

  !SOLVER
  FMTr = "(a30,E16.7e3)"
  FMTi = "(a30,I10)"
  write(fnum,*) "Solver-related information:"
  write(fnum,FMTr) "step size last",rwork(11)
  write(fnum,FMTr) "step size attempt",rwork(12)
  write(fnum,FMTr) "time current",rwork(13)
  write(fnum,FMTr) "tol scale factor",rwork(14)
  write(fnum,FMTi) "numeber of steps",iwork(11)
  write(fnum,FMTi) "call to fex",iwork(12)
  write(fnum,FMTi) "call to jex",iwork(13)
  write(fnum,FMTi) "last order used",iwork(14)
  write(fnum,FMTi) "order attempt",iwork(15)
  write(fnum,FMTi) "idx largest error",iwork(16)
  write(fnum,FMTi) "RWORK size required",iwork(17)
  write(fnum,FMTi) "IWORK size required",iwork(18)
  write(fnum,FMTi) "NNZ in Jac",iwork(19)
  write(fnum,FMTi) "extra fex to compute jac",iwork(20)
  write(fnum,FMTi) "number of LU decomp",iwork(21)
  write(fnum,FMTi) "base address in RWORK",iwork(22)
  write(fnum,FMTi) "base address of IAN",iwork(23)
  write(fnum,FMTi) "base address of JAN",iwork(24)
  write(fnum,FMTi) "NNZ in lower LU",iwork(25)
  write(fnum,FMTi) "NNZ in upper LU",iwork(21)
  write(fnum,*) "See DLSODES manual for further details on Optional Outputs"
  write(fnum,*)
  write(fnum,*) "END KROME ERROR REPORT"
  write(fnum,*)
  close(fnum)

  mx_dump = mx_dump - 1
  if (mx_dump==0) stop

end subroutine krome_dump

!********************************
subroutine krome_init()
  use krome_commons
  use krome_tabs
  use krome_subs
  use krome_reduction
  use krome_dust
  use krome_cooling
  use krome_photo
  use krome_fit
  use krome_getphys

  !init phys common variables
  !$omp parallel
  phys_Tcmb = 2.73d0
  phys_zredshift = 0d0
  phys_orthoParaRatio = 3d0
  phys_metallicity = 0d0
  phys_Tfloor = 2.73d0
  !$omp end parallel

  !init metallicity default
  !assuming solar
  total_Z = 1d0

  !default D/D_sol = Z/Z_sol
  !assuming linear scaling
  dust2gas_ratio = total_Z

  !default broadening turubulence velocity
  broadeningVturb2 = 0d0

  !default clumping factor for
  ! H2 formation on dust by Jura/Gnedin
  clump_factor = 1d0

  !default for thermo toggle is ON
  !$omp parallel
  krome_thermo_toggle = 1
  !$omp end parallel

  !load arrays with ractants/products indexes
  call load_arrays()

  !initialize the table for exp(-a/T) function
  call init_exp_table()

  call load_parts()

  !init photo reactants indexes

  !get machine precision
  krome_epsilon = epsilon(0d0)

  !load verbatim reactions
  call loadReactionsVerbatim()

end subroutine krome_init

!****************************
function krome_get_coe(x,Tgas)
  !krome_get_coe: public interface to obtain rate coefficients
  use krome_commons
  use krome_subs
  use krome_tabs
  implicit none
  real*8 :: krome_get_coe(nrea), x(nmols), Tgas
  real*8::n(nspec)

  n(:) = 0d0
  n(1:nmols) = x(:)
  n(idx_Tgas) = Tgas
  krome_get_coe(:) = coe_tab(n(:))

end function krome_get_coe

!****************************
function krome_get_coeT(Tgas)
  !krome_get_coeT: public interface to obtain rate coefficients
  ! with argument Tgas only
  use krome_commons
  use krome_subs
  use krome_tabs
  implicit none
  real*8 :: krome_get_coeT(nrea),Tgas
  real*8::n(nspec)
  n(idx_Tgas) = Tgas
  krome_get_coeT(:) = coe_tab(n(:))
end function krome_get_coeT

end module krome_main
