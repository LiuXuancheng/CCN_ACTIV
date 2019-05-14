program Parcel
        !aerosol activation for parcel of 1cm^3 in volumeume
        !///////////////////////////////////////////////////////////////////////////////!!
        ! list of variables 
        ! qvpp, qvs are the water vapor mixing ratio & its repective saturated value
        ! esat water vapor pressure
        ! temperature, temp_old, thetapp are new & old temperature & potential temperature
        ! delt: time step 
        ! time_prep: spin up time
        ! rhoa, rhow, rho_ccn, density of air, water & solute
        ! sp, seq are supersaturation wrt water & equilibrium value over solute particle
        ! h & volume are width & volume of parcel box
        ! lwc: liquid water content
        ! dsd: droplet size distribution with nbinsout width
        ! ndrop: total droplet number
        ! rad, rad_ccn: wet & dry radius of the particle
        ! kappa : hygroscopicity, see petters & Kredenweis 2007
        ! solu : solute term calculated either from kappa (isolu=1) or from classical format,
        ! see Jensen & Nugent 2017, eqn (1)
        ! disp: DSD dispersion flags
        ! flag_activ: to judge CCN activation
        ! activ_ccn: activated CCN concentration
        !///////////////////////////////////////////////////////////////////////////////!!
      implicit none
      real*8            :: qvpp,qvs,esat,ks,temperature,temp_old,time,delt,time_prep
      real*8            :: rhoa,sp,volume,h,thetapp,seq,seq1,seq2,sp_wet
      real*8            :: exner,racp,p1,PP,sumrp
      real*8            :: rm,rm0,curv,solu
      real*8            :: deltaqp,lwc,cql,cdp
      character*6       :: name
      real(8), allocatable, dimension(:) :: rad,rad_ccn,dsd,dsd_temp,dr3,rad_wet,kappa
      real(8), allocatable, dimension(:,:,:,:) ::frac_activ,time_table
      integer           :: nbins,nbinsout,iinit,ifinal, isolu,crit_bin
      integer           :: disp,GCCN
      integer           :: iter,ntmic,ntot,i,i_ta,j_ta,k_ta,s_ta,t_ta,ii,iii,i_ta_max,j_ta_max,k_ta_max,s_ta_max,t_ta_max
      real*8            :: ndrop
      real*8            :: diffvnd1,diffvnd2,ka1,ka2
      real, parameter :: p0=1.0e5 !reference pressure
      real, parameter :: diffvnd = 2.55d-5              ! Coefficient of diffusion of water vapour in air [m**2/s]
      real, parameter :: ka = 2.48d-2                   ! Thermal conductivity of air [J/msK]
      real,parameter :: pi=3.1415
      real,parameter ::  KK = 8.54d-11
      real,parameter :: grav=9.8
      real,parameter :: visc = 1.78e05!0.16d-4!1.78e-5
      real,parameter :: lat = 2.5e-6!2.477d6!2.5e-6
      real,parameter :: ra=287.0 
      real,parameter :: cp=1004.0!1005.0 
      real,parameter :: rv= 467!461.5!467
      real,parameter :: m_w=18.d-3 !molecular weight of water
      real,parameter :: rhow=1000.0 !density of water
      real,parameter :: eps=ra/rv
      real,parameter :: latovercp=lat/cp
      real*8,parameter :: sigma_sa=7.61d-2
      real,parameter :: upp=0.0 !spin-up updraft
      real :: up,vh !updraft velo and van Hoff factor 
      real :: m_s !molecular weight of solute; ammonium sulfate=132.14d-3; NaCl = 58.44d-3
      real :: rho_ccn!kg/m**3 for ammonium sulfate =1726.d0 for NaCl=2160.d0
      real :: flag_activ!to judge CCN activation
      real*8 ::activ_ccn!activated CCN concentration
      real*8, dimension(7):: temp_i,na_i,ns_i,meanD_i
      real*8, dimension(9):: up_i
      real*8, dimension(4):: dt_i
      data na_i/10.0, 31.6, 100.0, 316.0, 1000.0, 3160.0, 10000.0/
      data ns_i/1.0, 3.16, 10.00, 31.60, 100.00, 316.00, 1000.00/
      data up_i/0.01, 0.0316, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0/
      data temp_i/243.15, 253.15, 263.15, 273.15, 283.15, 293.15, 303.15/
      
196   format(1x,6(f16.8,2x))
145   format(1x,100(e16.8,2x))
146   format(1x,7(e16.8,2x))
      
      i_ta_max = size(up_i)
      j_ta_max = size(temp_i)
      k_ta_max = size(na_i)
      s_ta_max = size(ns_i)
      allocate(frac_activ(s_ta_max,k_ta_max,i_ta_max,j_ta_max))
      allocate(time_table(s_ta_max,k_ta_max,i_ta_max,j_ta_max))
      time_table=0.0
      time_prep=.0
      delt= 2.d-04 
      ntot=10./delt
      name = 'simple'
      !------------------------------setup output files-------------------------------!!
      OPEN(UNIT=16,FILE=name//'.out')!parcel mean variables
      OPEN(UNIT=50,FILE=name//'.dsd')!number concentration of each bin
      OPEN(UNIT=51,FILE=name//'.rad')!droplet size of each bin
      OPEN(UNIT=60,FILE=name//'.test')!output the tested variables
      OPEN(UNIT=52,FILE=name//'.test_relax_time')!output the tested variables
      OPEN(UNIT=61,FILE='CCNACTIV_10S_PRISTINE.BIN',FORM='UNFORMATTED')!output the tested variables
      OPEN(UNIT=62,FILE='CCNACTIV_10S_PRISTINE')!output the tested variables
      !--------------------------------------initialize aerosols----------------------!!
      !------------------------------initial size distribution------------------------!!
      !////GCCN=1 insert Jensen & Nugent 2017 Giant CCN           ////////////////////!!
      !//////   =2 monotonic seeding r=1micron;                   ////////////////////!!
      !//////=3 add 3 mode lognormal seeding distribution by Cooper et al. 1997///////!!
      !///disp=35 Xue10 urban, 30 Xue10 marine, 31 JN17 polluted, 32 NJ17 pristine////!!
      !////disp=1 monotonic initial sizes ////////////////////////////////////////////!!
      !-------------------------------------------------------------------------------!!
      disp = 32
      GCCN = 1 
      isolu = 1            !=1 kappa form solute term; =2 classical solute term
      nbins = 100
      allocate (dsd(nbins),dsd_temp(nbins),rad(nbins),rad_ccn(nbins),dr3(nbins),rad_wet(nbins),kappa(nbins))
      
      if (disp .le. 2) then !bi-disperse
	      rho_ccn=1726.d0!ammonium sulfate
	      m_s=132.14d-3
	      vh = 3.!2.
      elseif (disp .eq. 30 .or. disp .eq. 35 ) then !Xue10  case
	      rho_ccn=1726.d0!2160.d0!1726.d0 !ammonium sulfate
	      m_s=132.14d-3!58.44d-3!132.14d-3 
	      vh = 3.!2.
      elseif (disp .eq. 31 .or. disp .eq. 32) then !NJ17
	      rho_ccn=2160.d0 !NaCl
	      m_s=58.44d-3
	      vh = 2.
	      GCCN=1
      endif
      
      call iaerosol(disp,rad_ccn,dsd,nbins,ndrop,rho_ccn,rm,nbinsout,crit_bin,GCCN)
      dsd_temp=dsd
      
      do s_ta = 1,s_ta_max
      do k_ta = 1,k_ta_max
          
      print*,'dispersion type ',disp
      dsd=dsd_temp
      dsd(1:crit_bin)=dsd(1:crit_bin)*na_i(k_ta)/sum(dsd(1:crit_bin))
      dsd(crit_bin+1:nbins)=dsd(crit_bin+1:nbins)*ns_i(s_ta)/sum(dsd(crit_bin+1:nbins))
      ndrop = na_i(k_ta)+ns_i(s_ta)
      
      do i_ta = 1,i_ta_max
      do j_ta = 1,j_ta_max
      ii=0
      iii=0
      up = up_i(i_ta)!2.0             !  updraft velocity
      sp = 0.0!-14.39d-2       !supersaturation
      sp_wet = -0.01!-14.39d-2
      temperature=temp_i(j_ta)!284.3d0  !initial temperature
      p1=93850.0d0         !initial pressure
      activ_ccn = 0.0
      !-------------------------------------------------------!
      rm =0.d0
      rad=rm
      dr3=0.d0
      ! to get dry radius rad_ccn
      print*,'maximum binsize is',nbinsout
      if (GCCN .ne. 0) print*,'GCCN is on, value = ',GCCN
      write(50,145) 0.,(dsd(i), i=1,nbinsout)
      write(51,145) 0.,(rad(i), i=1,nbinsout)
      print*,'number of drops ',sum(dsd)!na_i(k_ta)!ndrop
      
!-------------initialize variables------------
      if(isolu .eq. 1) then
         kappa(1:nbinsout)=vh*m_w/m_s*rho_ccn/rhow
         print*, 'kappa = ', kappa(1)
      endif
      h = 1.d-2 !.01m=1cm
      volume = h**3
      pp=p1
      sumrp=0.d0
      racp= ra/cp
      exner=(pp/p0)**racp
      thetapp=temperature/exner
      cql=4.0d0*pi*rhow/volume
      esat = 2.53d11*exp(-5.42d3/temperature)
      ks = 1.d0/(rhow*Rv*temperature/(esat*diffvnd)+rhow*Lat/(Ka*temperature)*(Lat/(Rv*temperature)-1))
      qvs = eps*esat/(PP-esat)
      qvpp= (sp+1.d0)*qvs !kg/m^3
      rhoa=pp/(Ra*(1+18.d0/29.d0*qvpp)*temperature)
!--------------first guess the radius of wet aerosol--------!
!--------------start with r_wet=1.8*r_d-----------------------!
      if(s_ta .eq. 1 .and. k_ta .eq. 1 .and. i_ta .eq. 1)then
      rad_wet=rad_ccn*1.8d0
      diffvnd1=1.d-5*(0.015*temperature-1.9)
      ka1=1.5d-11*temperature**3-4.8d-8*temperature**2+1.d-4*temperature-3.9d-4
      do i=1,nbinsout
         seq=sp_wet+.01d0
	      seq1=seq+.01d0
	      seq2=seq+.01d0
      do while(abs(seq-sp_wet) .gt. 1.d-7 .and. seq2 .ne. seq) 
	      seq2=seq1
	      seq1=seq
         !diffvnd2=diffvnd1*1.d0/(rad_wet(i)/(rad_wet(i)+0.104d-6)+diffvnd1/(rad_wet(i)*0.036)*sqrt(2.d0*pi/(Ra*temperature)))
         !ka2=ka1*1.d0/(rad_wet(i)/(rad_wet(i)+.216d-6)+ka1/(rad_wet(i)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temperature)))
         !ks = 1.d0/(rhow*Rv*temperature/(esat*diffvnd2)+rhow*Lat/(Ka2*temperature)*(Lat/(Rv*temperature)-1))
         curv=2.d0*sigma_sa/(Rv*rhow*temperature*rad_wet(i))  

         if(isolu .eq. 1) then !kappa
            solu=(rad_wet(i)**3-rad_ccn(i)**3)/(rad_wet(i)**3-(1-kappa(i))*rad_ccn(i)**3)
         elseif (isolu .eq. 2) then !classical solute term
            solu=exp(-vh*m_w/m_s*rho_ccn/rhow*rad_ccn(i)**3/(rad_wet(i)**3-rad_ccn(i)**3)) !solute effect coefficient ms=132.14 for ammonium sulfate !m^3
         endif !isolu

         seq = solu * exp(curv)-1.d0           
         if(seq .gt. sp_wet .and. rad_wet(i) .gt. rad_ccn(i)) then
	         rad_wet(i)=rad_wet(i)-rad_ccn(i)*1.d-3
	      elseif(seq .lt. sp_wet) then 
	         rad_wet(i)=rad_wet(i)+rad_ccn(i)*1.d-3
         endif
         ii =ii+1
      enddo
      rad(i)=rad_wet(i) !droplet radius
      enddo
      print*, 'wet equlibrium: ',ii
      endif
!--------------spin-up------------------
if(time_prep .ne. 0) then
      do 200 ntmic=1,int(time_prep/delt*2.d0)
         time = ntmic*delt/2.d0-time_prep
         pp=rhoa*Ra*(1.d0+18.d0/29.d0*qvpp)*temperature
         exner = (PP/P0)**RACP
         sumrp=sum(dr3(1:nbinsout)*dsd(1:nbinsout))/3.d0
         deltaqp=cql*sumrp!new cdp condensation
         temp_old=temperature
         temperature=temp_old-grav/cp*delt/2.0d0*upp+latovercp*deltaqp!mark new
         esat=2.53d11*exp(-5.42d3/temperature)
         qvs = eps*esat/(PP-esat)
         rhoa= rhoa*(-grav*upp/(Ra*temperature)*delt/2.d0-(temperature-temp_old)/temperature)+rhoa
          thetapp=thetapp+latovercp*deltaqp/exner
          qvpp    = qvpp - deltaqp
          sp = qvpp/qvs-1.d0 !mark new sp
          rm=0.d0
          diffvnd1=1.d-5*(0.015*temperature-1.9)
          ka1=1.5d-11*temperature**3-4.8d-8*temperature**2+1.d-4*temperature-3.9d-4
          do i = 1,nbinsout
             curv=2.d0*sigma_sa/(Rv*rhow*temperature*rad(i)) !curvature effect coefficient !unit in m
            if(isolu .eq. 1) then !kappa
               solu=(rad(i)**3-rad_ccn(i)**3)/(rad(i)**3-(1-kappa(i))*rad_ccn(i)**3)
            elseif (isolu .eq. 2) then !classical solute term
               solu=exp(-vh*m_w/m_s*rho_ccn*rad_ccn(i)**3/rhow/(rad(i)**3-rad_ccn(i)**3))
            endif !isolu
            
            diffvnd2=diffvnd1*1.d0/(rad(i)/(rad(i)+0.104d-6)+diffvnd1/(rad(i)*0.036)*sqrt(2.d0*pi/(Ra*temperature)))     
            ka2=ka1*1.d0/(rad(i)/(rad(i)+.216d-6)+ka1/(rad(i)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temperature)))
            ks =1.d0/(rhow*Rv*temperature/(esat*diffvnd2)+rhow*Lat/(Ka2*temperature)*(Lat/(Rv*temperature)-1))
            !!caculate equillibrium supersat.
             seq = solu * exp(curv)-1.d0
            print*,'solu=',solu, 'seq=',seq, 'supersat=',sp-seq!mark
             rm0=rad(i)*3.0d0*delt*ks*(sp-seq)+rad(i)**3 !!r^3 scheme
             if(rm0 .gt. rad_ccn(i)**3 ) then
                dr3(i)=rm0-rad(i)**3
                rad(i)=rm0**(1.d0/3.d0)
             else
                dr3(i)=0.d0
                rad(i)=rad_ccn(i)
             endif
          enddo
          rm=(sum(rad(1:nbinsout)**3*dsd(1:nbinsout))/ndrop)**(1.d0/3.0d0)
          if(mod(ntmic,int(1.d0/delt)) .eq. 0 .or. int(time_prep/delt*2.d0) .le. 1000) then
             write(16,*) time,0,Sp,0.,rad(15),rad(40),pp,temperature,281.5-grav/cp*delt*up*ntmic,&
               thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
             write(50,145) time,(dsd(i), i=1,nbinsout)
             write(51,145) time,(rad(i), i=1,nbinsout)
	  endif
  200 enddo
endif ! spin_up
lwc=sum(4.d0/3.d0*pi*rad(1:nbinsout)**3*rhow*dsd(1:nbinsout))
write(16,*) 0.d0-time_prep,up*delt*ntmic-287.6,Sp,lwc,rad(15),rad(40),pp,temperature,281.5-grav/cp*delt*up,& !8
            thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
!--------------evolumeution------------------
      do 100 ntmic=1,ntot
         time = ntmic*delt
         pp=rhoa*Ra*(1+m_w/29.d-3*qvpp)*temperature
         exner = (PP/P0)**RACP         
         sumrp=sum(dr3(1:nbinsout)*dsd(1:nbinsout))/3.d0
	      deltaqp=cql*sumrp!cdp condensation
         temp_old=temperature
         temperature=temp_old-grav/cp*delt*up+latovercp*deltaqp!mark new
         esat=2.53d11*exp(-5.42d3/temperature)
         qvs = eps*esat/(PP-esat)
         rhoa= rhoa*(-grav*up/(Ra*temperature)*delt-(temperature-temp_old)/temperature)+rhoa
         thetapp=thetapp+latovercp*deltaqp/exner
         qvpp    = qvpp - deltaqp
	      sp = qvpp/qvs-1.d0 !mark new sp
         rm=0.d0
         diffvnd1=1.d-5*(0.015*temperature-1.9)
         ka1=1.5d-11*temperature**3-4.8d-8*temperature**2+1.d-4*temperature-3.9d-4
          do i = 1,nbinsout
             curv=2.d0*sigma_sa/(Rv*rhow*temperature*rad(i)) !curvature effect coefficient !in m
            if(isolu .eq. 1) then !kappa
               solu=(rad(i)**3-rad_ccn(i)**3)/(rad(i)**3-(1-kappa(i))*rad_ccn(i)**3)
            elseif (isolu .eq. 2) then !classical solute term
               solu=exp(-vh*m_w/m_s*rho_ccn*rad_ccn(i)**3/rhow/(rad(i)**3-rad_ccn(i)**3)) !solute effect coefficient ms=132.14 for ammonium sulfate !m^3
            endif !isolu
            diffvnd2=diffvnd1*1.d0/(rad(i)/(rad(i)+0.104d-6)+diffvnd1/(rad(i)*0.036)*sqrt(2.d0*pi/(Ra*temperature)))
            ka2=ka1*1.d0/(rad(i)/(rad(i)+.216d-6)+ka1/(rad(i)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temperature)))
            ks =1.d0/(rhow*Rv*temperature/(esat*diffvnd2)+rhow*Lat/(Ka2*temperature)*(Lat/(Rv*temperature)-1))
!!--------------caculate equillibrium supersat.
            seq=solu*exp(curv)-1.d0
            rm0=rad(i)*3.0d0*delt*ks*(sp-seq)+rad(i)**3 !!r^3 
	         if(rm0 .gt. rad_ccn(i)**3) then
	            dr3(i)=rm0-rad(i)**3
	            rad(i)=rm0**(1.d0/3.d0)
	         else
		         dr3(i)=0.d0
               rad(i)=rad_ccn(i)
	         endif
             
             if(iii .eq. 0) call ACTIVE_OR_NOT(rad_ccn(i),rad(i),temperature,kappa(i),flag_activ)
             if(flag_activ .eq. 1. .and. iii .eq. 0) THEN
             PRINT*, 'ACTIVATION TIME IS: ',TIME
             time_table(s_ta,k_ta,i_ta,j_ta) = dble(time)
             iii=1
             endif
             
             if(ntmic .eq. ntot) then!judege whether CCN is activated
                 call ACTIVE_OR_NOT(rad_ccn(i),rad(i),temperature,kappa(i),flag_activ)
                 if(flag_activ .eq. 1.) activ_ccn =activ_ccn + dsd(i)
                 if(i==nbinsout)print*,'ACTIV_CCN :',activ_ccn
             endif
             
          enddo !size loop
            rm=(sum(rad(1:nbinsout)**3*dsd(1:nbinsout))/ndrop)**(1.d0/3.0d0)
            if(mod(ntmic,int(1./delt)) .eq. 0 .or. ntot .le. 1000) then
               lwc=sum(4.d0/3.d0*pi*rad(1:nbinsout)**3*rhow*dsd(1:nbinsout))
              write(16,*) time,up*delt*ntmic-287.6,Sp,lwc,rad(15),rad(40),pp,temperature,281.5-grav/cp*delt*up*ntmic,&
                thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
              write(50,145) time,(dsd(i), i=1,nbinsout)
              write(51,145) time,(rad(i), i=1,nbinsout)
            endif
      100 enddo!time loop
      frac_activ(s_ta,k_ta,i_ta,j_ta) = activ_ccn/ndrop
      enddo!j_ta
      write(52,146) (time_table(s_ta,k_ta,i_ta,j_ta),j_ta = 1,j_ta_max)
      write(62,146) (frac_activ(s_ta,k_ta,i_ta,j_ta),j_ta = 1,j_ta_max)
      enddo!i_ta
      enddo!k_ta
      enddo!s_ta
      
      !do k_ta = 1,k_ta_max
      !write(61,'(A,F8.2)')'Number = ', na_i(k_ta)
      write(61,err=522)frac_activ
522   continue
      !enddo

      close(unit=16)
      deallocate (dsd, rad, rad_ccn,dr3,rad_wet,frac_activ,time_table)
      close(unit=50)
      close(unit=51)
      close(unit=60)
      close(unit=61)
      close(unit=62)
      close(unit=52)


      end program parcel

    SUBROUTINE IAEROSOL(disp,rad,nrad,nbins,ndrop,rho_ccn,rm,nbinsout,crit_bin,GCCN)
!---- This subroutine determines the initial position and size of all droplets
  implicit none

  ! --- argument ---
  integer :: nbins,nbinsout,GCCN,crit_bin
  real(8), dimension(nbins) :: rad,nrad
  integer :: disp
  real(8) :: rm,ndrop

  ! --- local --

  integer :: i,iinit,ifinal
  real(8), allocatable, dimension(:) :: wid,dNdlogr,dNdr
  real(8) :: r1,n1,logsig,rmin,rmax
  real :: rho_ccn
  real(8) :: logrmin,logrmax,rad_power,bin_factor
  real,parameter :: pi=3.14159265
 145   format(1x,100(e16.8,2x))
 111 format(a20,i3)

! Set everything to zero.
rad = 0.0
ndrop=0.0 
rm=0.d0
nrad=0.
allocate(wid(nbins),dNdlogr(nbins),dNdr(nbins))
dNdlogr = 0.0
dNdr = 0.0d0
wid =0.d0
!size dispersion
  if (disp .eq. 1) then !mono disperse
      nbinsout=1
      rad(nbinsout)=1.d-7
      ndrop=100
      nrad=ndrop
  elseif (disp .eq. 30) then !lulin 2010 maritime case
     nbinsout=39
     rmin = 6.d-9
     rad(1)=rmin
     bin_factor=2.d0
     wid(1)=rad(1)*(bin_factor**(1.d0/3.d0)-1.d0)
     do i=2,nbinsout
        rad_power=real(i)/3.d0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
     enddo
     do i=1,nbinsout
        n1=133.d0
        r1=0.0039d-6
        logsig=.657d0
        logsig=log(10.d0)*logsig
        dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=66.6d0
        r1=.133d-6
        logsig=.21d0
        logsig=log(10.d0)*logsig
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=3.06d0
        r1=.29d-6
        logsig=.396d0
        logsig=log(10.d0)*logsig
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
	      if (GCCN .eq. 3) then !add 3 mode lognormal seeding distribution by Cooper et al. 1997
	      !mode 1
	         n1=100.d0
	         r1=.15d-6
	         logsig=.2d0
	         logsig=log(10.d0)*logsig
            dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
	      !mode 2
            n1=100.d0*1.7d-4
            r1=.5d-6
            logsig=.4d0
            logsig=log(10.d0)*logsig
            dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
	      !mode 3
            n1=100.d0*3.d-7
            r1=5.d-6
            logsig=.6d0
            logsig=log(10.d0)*logsig
            dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
	      endif!GCCN=3
     enddo
     do i=1,nbinsout 
        dNdr(i)=dNdlogr(i)/rad(i)
        nrad(i)=dNdr(i)*wid(i)
     enddo
  elseif (disp .eq. 31) then !Jensen&Nugent 2017 modified polluted case
     nbinsout=30
     rmin = 1.d-8
     rad(1)=rmin
     bin_factor=50.**0.1!2.0d0 !mass increment
     wid(1)=rad(1)*(bin_factor**(1.d0/3.d0)-1.d0)
     do i=2,nbinsout
        rad_power=real(i)/3.d0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
     enddo
     do i=1,nbinsout
        n1=48.d0!160!48.d0
        r1=.029d-6
        logsig=1.36d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) *exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=125.d0!380!125.d0
        r1=.071d-6
        logsig=1.57d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig)*exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
     enddo
     dNdr(1:nbinsout)=dNdlogr(1:nbinsout)/rad(1:nbinsout)
     nrad(1:nbinsout)=dNdr(1:nbinsout)*wid(1:nbinsout)


  elseif (disp .eq. 32) then !Jensen&Nugent 2017 pristine case
     nbinsout=30
     rmin = 1.d-8
     rad(1)=rmin
     bin_factor=50.**0.1!2.0d0 !mass increment
     wid(1)=rad(1)*(bin_factor**(1.d0/3.d0)-1.d0)
     do i=2,nbinsout
        rad_power=real(i)/3.d0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
     enddo
     do i=1,nbinsout
        n1=125.d0
        r1=.011d-6
        logsig=1.2d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig)*exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=65.d0
        r1=.06d-6
        logsig=1.7d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig)*exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
     enddo
        dNdr(1:nbinsout) = dNdlogr(1:nbinsout)/rad(1:nbinsout)
        nrad(1:nbinsout) = dNdr(1:nbinsout)*wid(1:nbinsout)
  elseif (disp .eq. 35) then !Lulin 2010 rural
     rmin = 6.d-9
     rad(1)=rmin
     bin_factor=2.d0 !mass increment
     wid(1)=rad(1)*(bin_factor**(1.d0/3.0d0)-1)
     nbinsout=35
     do i=2,nbinsout
        rad_power=real(i)/3.0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
        n1=6650.d0
        r1= 0.00739d-6
        logsig= .225d0
        logsig=log(10.d0)*logsig
        dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1= 147.d0    
        r1= .0269d-6
        logsig= .557
        logsig=log(10.d0)*logsig
        dNdlogr(i) = dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1= 1990.d0
        r1=.0419d-6
        logsig= .266
        logsig=log(10.d0)*logsig
        dNdlogr(i) = dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
     enddo
     do i=1,35
        dNdr(i)=dNdlogr(i)/rad(i)
        nrad(i)=dNdr(i)*wid(i)
     enddo
 
  elseif(disp .eq. 39 ) then !Jaenicke1988
     rmin = 6.d-9
     rmax = 5.d-6
     logrmin=10.d0**floor(log10(rmin))
     logrmax=10.d0**floor(log10(rmax))
     iinit=Nint(rmin/logrmin)
     ifinal=(floor(log10(rmax))-floor(log10(rmin)))*9+floor(rmax/logrmax)
     r1=0.
     n1=0.
     logsig=0.
     nbinsout=ifinal
     do i=iinit,nbinsout
          rad(i) = real(mod(i,9))*10.d0**(-7+i/9)
          if (mod(i,9) .eq. 0) then
             rad(i) = 9.0d0*10.d0**(-7+i/9-1)
          endif
          n1=133.d0
          r1=0.0039d-6
          logsig=.657d0
          dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          n1=66.6d0
          r1=.133d-6
          logsig=.21d0
          dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          n1=3.06d0
          r1=.29d-6
          logsig=.396d0
          dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
     enddo
  endif ! disp
  crit_bin=nbinsout
  if(GCCN .eq. 1) then !JN17 GCCN
     do i=1,42
        rad(nbinsout+i)=.8d-6+0.2d-6*real(i-1)
     enddo
        nrad(nbinsout+1)=0.1118d0
        nrad(nbinsout+2)=.06849d0!1micron
        nrad(nbinsout+3)=.0384d0
        nrad(nbinsout+4)=.02182d0
        nrad(nbinsout+5)=.0133d0
        nrad(nbinsout+6)=.8496d-2
        nrad(nbinsout+7)=.5486d-2!2micron
        nrad(nbinsout+8)=.3805d-2
        nrad(nbinsout+9)=.2593d-2
        nrad(nbinsout+10)=.1919d-2
        nrad(nbinsout+11)=.1278d-2
        nrad(nbinsout+12)=.9884d-3!3micron
        nrad(nbinsout+13)=.7779d-3
        nrad(nbinsout+14)=.5195d-3
        nrad(nbinsout+15)=.4005d-3
        nrad(nbinsout+16)=.3769d-3
        nrad(nbinsout+17)=.2653d-3!4micron
        nrad(nbinsout+18)=.2124d-3
        nrad(nbinsout+19)=.1378d-3
        nrad(nbinsout+20)=.1214d-3
        nrad(nbinsout+21)=.1009d-3
        nrad(nbinsout+22)=.1222d-3!5micron
        nrad(nbinsout+23)=.5064d-4
        nrad(nbinsout+24)=.383d-4
        nrad(nbinsout+25)=.5547d-4
        nrad(nbinsout+26)=.2145d-4
        nrad(nbinsout+27)=.1295d-4!6micron
        nrad(nbinsout+28)=.4323d-4
        nrad(nbinsout+29)=.2626d-4
        nrad(nbinsout+30)=.305d-4
        nrad(nbinsout+31)=.4385d-5
        nrad(nbinsout+32)=.4372d-5!7micron
        nrad(nbinsout+33)=.4465d-5
        nrad(nbinsout+34)=.4395d-5
        nrad(nbinsout+35)=.4427d-5
        nrad(nbinsout+36)=.4411d-5
        nrad(nbinsout+37)=.0d0    !8micron
        nrad(nbinsout+38)=.0d0
        nrad(nbinsout+39)=.0d0
        nrad(nbinsout+40)=.4522d-5
        nrad(nbinsout+41)=.0d0
        nrad(nbinsout+42)=.4542d-5!9micron
        nbinsout=nbinsout+42
  elseif (GCCN==2) then !some simple one size GCCN
        rad(nbinsout+1)=1.d-6
        nrad(nbinsout+1)=10
        nbinsout=nbinsout+1
  endif!GCCN
  ndrop=sum(nrad(1:nbinsout)) !total number
  write(50,145) 0.,(dNdlogr(i), i=1,nbinsout) !initial dry size distribution
  rm = (sum(rad(1:nbinsout)**3*nrad(1:nbinsout))/ndrop)**(1.d0/3.d0)
  print*,"rm", rm, 'nbins',nbins,'nbinsout',nbinsout,'ndrop',ndrop
deallocate(wid,dNdlogr,dNdr)

 299   format(1x,(e12.6),3(f8.6))
 199   format(1x,3(f8.6))


  end SUBROUTINE IAEROSOL

real function log2(x)
  implicit none
  real(8), intent(in) :: x

  log2 = log(x) / log(2.d0)
end function

SUBROUTINE ACTIVE_OR_NOT(Rd,R,T,K,FLAG_ACTIV)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: Rd, R, T,K
    REAL*8 :: R_CRI,A,B,C,D
    REAL*8 :: rv= 467.!461.5!467
    REAL*8 :: rhow=1000.0 !density of water
    REAL*8 :: sigma_sa=7.61d-2
    REAL,INTENT(INOUT) :: FLAG_ACTIV
    FLAG_ACTIV = 0.0
    R_CRI = 0.0
    A = K*(1-K)
    B = -(R**3.0+3*K*rv*rhow*T/2.0/sigma_sa*R**4.0)
    C = R**6.0
    D = B**2.-4.*A*C
    IF(D .GE. 0.0)THEN
        R_CRI = ((-B-(B**2.-4.*A*C)**0.5)/2./A)**0.333333
        !print*, 'R_CRI: ',R_CRI
        IF (Rd .LE. R_CRI)THEN
            FLAG_ACTIV = 1.
        ENDIF
    ENDIF
    END SUBROUTINE
