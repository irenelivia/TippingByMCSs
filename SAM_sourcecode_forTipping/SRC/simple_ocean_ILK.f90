module simple_ocean

!------------------------------------------------------------
! Purpose:
!
! A collection of routines used to specify fixed 
! or compute interactive SSTs, like slab-ocean model, etc.
!
! Author: Marat Khairoutdinov
! Based on dynamic ocean impelemntation from the UW version of SAM.
!------------------------------------------------------------

use grid
implicit none

public set_sst     ! set SST 
public sst_evolve ! evolve SST according to a model set by the ocean_type

CONTAINS


SUBROUTINE set_sst()

 use vars, only: sstxy, t00, sst_cor
 use params, only: tabs_s, delta_sst, ocean_type

! parameters of the sinusoidal SST destribution 
! along the X for Walker-type simulatons( ocean-type = 1):

 real(8) tmpx(nx), pii, lx, yy, ly
 integer i,j, it,jt

 ! Initialize the slab layer integral corrector to zero
 sst_cor = 0.

 select case (ocean_type)

   case(0) ! fixed constant SST

      sstxy = tabs_s - t00  ! NOTE: sstxy is the perturbation from t00

   case(1) ! Sinusoidal distribution along the x-direction:

     lx = float(nx_gl)*dx
     do i = 1,nx
        tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
     end do
     pii = atan2(0.d0,-1.d0)
     do j=1,ny
       do i=1,nx
         sstxy(i,j) = tabs_s-delta_sst*cos(2.*pii*tmpx(i)/lx) - t00
       end do
     end do
   
   case(2) ! Sinusoidal distribution along the y-direction, peak at center:
     
     call task_rank_to_index(rank,it,jt)
     
     pii = atan2(0.d0,-1.d0)
     ly = float(ny_gl)*dy
     do j=1,ny
        yy = dy*(j+jt-(ny_gl+YES3D-1)/2-1)
       do i=1,nx
         sstxy(i,j) = tabs_s+delta_sst*(2.*cos(pii*yy/ly)-1.) - t00 ! Half cos
         !sstxy(i,j) = tabs_s+delta_sst*cos(2.*pii*yy/ly) - t00 ! Full cos
       end do
     end do

   case(3) !  Aquaplanet Experiment (APE) "QOBS" (Neale and Hoskins 2001)

     call task_rank_to_index(rank,it,jt)

     pii = atan2(0.d0,-1.d0)
     lx = float(nx_gl)*dx
     ly = 0.5*float(ny_gl)*dy
     do j=1,ny
        yy = dy*(j+jt-(ny_gl+YES3D-1)/2-1)
       do i=1,nx
         sstxy(i,j) = tabs_s+delta_sst*(1-0.5*(sin(1.5*pii/180.*yy/40000000.*360.)**2 &
                                              +sin(1.5*pii/180.*yy/40000000.*360.)**4))-t00
       end do
     end do

   case(4) ! fixed SST that will evolve with time at timestep == timesimpleocean

      sstxy = tabs_s - t00  ! NOTE: sstxy is the perturbation from t00




   case default

     if(masterproc) then
         print*, 'unknown ocean type in set_sst. Exitting...'
         call task_abort
     end if

 end select

end subroutine set_sst



SUBROUTINE sst_evolve
 use vars, only: sstxy, t00, fluxbt, fluxbq, latitude, rhow, qocean_xy, sst_cor
 use params, only: cp, pi, lcond, tabs_s, ocean_type, dossthomo, &
                   depth_slab_ocean, Szero, deltaS, delta_sst, timesimpleland, timesimpleocean, ocean_type
 use rad, only: swnsxy, lwnsxy

 real, parameter :: rhor = 1000. ! density of water (kg/m3)
 real, parameter :: cw = 4187.   ! Liquid Water heat capacity = 4187 J/kg/K
 real, parameter :: gain_proportional = 0.00001389 ! gain of the proportional correction: 0.1 x 2h relaxation timescale
 real, parameter :: gain_integral = 0.000000002315 ! gain of the integral correction: 0.0001 x 12h relaxation timescale
 real factor_cp, factor_lc, qoceanxy
 real yy, ly, pii
 real(8) sss(1),ssss(1)
 integer i,j,it,jt

 if(time.lt.timesimpleland) return   

    call task_rank_to_index(rank,it,jt)
 
    pii = atan2(0.d0,-1.d0)
    ly = float(ny_gl)*dy
      ! If the ocean_type equals 4 (starting constant everywhere), then we enforce a sinusoidal diurnal cycle
      ! with amplitude deltaS and a zero phase, that is max SST at midday
    if (ocean_type.eq.4) then

      call set_sst()
    
      sstxy = sstxy - delta_sst*cos(2*pii*time/86400)

      if (time.gt.timesimpleocean) then
        sstxy = tabs_s - t00 
      end if

      ! If the ocean_type equals 2 (starting with sin profile along y-dir), then we enforce a sinusoidal diurnal cycle
      ! with amplitude deltaS and a zero phase, that is max SST at midday
    !else if (ocean_type.eq.2) then

     ! call set_sst()
    
      !do j=1,ny
       ! do i=1,nx
        !  sstxy(i,j) = sstxy(i,j) - delta_sst*cos(2*pii*time/86400)
        !enddo
      !enddo




  else ! Then we use the slab layer model

    ! Define weight factors for the mixed layer heating due to
    ! the model's sensible and latent heat flux.
    factor_cp = rhow(1)*cp
    factor_lc = rhow(1)*lcond

    ! Use forward Euler to integrate the differential equation
    ! for the ocean mixed layer temperature: dT/dt = S - E.

    ! NEW: Vary the slab layer depth and heat sink along the y-coordinate.
    ! Both layer and heat sink increase from the center at specific rates
    ! which are given in the input file. Note that by default they equal 0, i.e.
    ! the surface layer is homogene. As ny needs to be even, the n/2 and n/2+1 element
    ! share the same values anyway.

    ! Romain: We now apply a proportional-integral correction towards tabs_s - hardcoded rates
    ! Note that the integral corrector is applied on the domain-avg SST, so we evaluate it first
    sss = 0.
    do j=1,ny
      do i=1,nx
        sss(1) = sss(1) + sstxy(i,j)
      end do
    end do
    sss(1) = sss(1) / dble(nx*ny)
    if(dompi) then
        call task_sum_real8(sss,ssss,1)
        sss = ssss /float(nsubdomains)
    end if

    sst_cor = sst_cor - dtn*gain_integral*sss(1)

    do j=1,ny

       yy = dy*(j+jt-(ny_gl+YES3D-1)/2-1)
       qoceanxy = Szero + deltaS*cos(2*pii*yy/ly)

       do i=1,nx
          qocean_xy(i,j) = qocean_xy(i,j) + qoceanxy*dtfactor

          sstxy(i,j) = sstxy(i,j) &
               + dtn*(swnsxy(i,j)          & ! SW Radiative Heating
               - lwnsxy(i,j)               & ! LW Radiative Heating
               - factor_cp*fluxbt(i,j)     & ! Sensible Heat Flux
               - factor_lc*fluxbq(i,j)     & ! Latent Heat Flux
               + qoceanxy)                 & ! Ocean Heating
               /(rhor*cw*depth_slab_ocean)   ! Convert W/m^2 Heating to K/s

          sstxy(i,j) = sstxy(i,j) - dtn*gain_proportional*sstxy(i,j) + sst_cor

       end do
    end do

    if(dossthomo) then
      sss = 0.
      do j=1,ny
       do i=1,nx
         sss(1) = sss(1) + sstxy(i,j)
       end do
      end do
      sss(1) = sss(1) / dble(nx*ny)
      if(dompi) then
          call task_sum_real8(sss,ssss,1)
          sss = ssss /float(nsubdomains)
      end if ! dompi
      if(ocean_type.eq.2) then
          tabs_s = sss(1) + t00
          call set_sst()
      else
         sstxy(:,:) = sss(1)
      end if
    end if
  endif

end subroutine sst_evolve


end module simple_ocean
