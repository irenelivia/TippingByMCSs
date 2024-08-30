 subroutine write_all()

  use vars
  implicit none
  character *4 rankchar
  character *256 filename
  character *256 filename_backup
  character *10 nstepchar
  integer irank
  integer, external :: lenstr

  call t_startf ('restart_out')

  if(masterproc) then
   write(nstepchar,'(i10)') nstep !new line for writing iteration number in restart filename
   print*,'Writing restart file ...'
   filename =  './RESTART/'//trim(case)//'_'//trim(caseid)//'_misc_restart.bin'
   !filename_backup = './RESTART/'//trim(case)//'_'//trim(caseid)//'_misc_restart_'//nstepchar(11-lenstr(nstepchar):10)//'.bin' 
                      !new name with iteration number so restart files aren't overwritten
   
   open(66,file=trim(filename), status='unknown',form='unformatted')
   !open(76,file=trim(filename_backup), status='unknown',form='unformatted') !file number 76+ used for backup restart files
  end if


	if(restart_sep) then
          write(nstepchar,'(i10)') nstep !new line for writing iteration number in restart filename
          write(rankchar,'(i4)') rank

          filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'
          !filename_backup = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
          !      rankchar(5-lenstr(rankchar):4)//'_restart_'//nstepchar(11-lenstr(nstepchar):10)//'.bin'


          open(65,file=trim(filename), status='unknown',form='unformatted')
          !open(75,file=trim(filename_backup), status='unknown',form='unformatted')
          write(65) nsubdomains, nsubdomains_x, nsubdomains_y
          !write(75) nsubdomains, nsubdomains_x, nsubdomains_y


	  call write_statement(.false.)


	else
    write(nstepchar,'(i10)') nstep !new line for writing iteration number in restart filename
	  write(rankchar,'(i4)') nsubdomains
	  filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart.bin'
    filename_backup = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart_'//nstepchar(11-lenstr(nstepchar):10)//'.bin'

	  do irank=0,nsubdomains-1
	
	     call task_barrier()

	     if(irank.eq.rank) then

	       if(masterproc) then
	      
	        open(65,file=trim(filename), status='unknown',form='unformatted')
          !open(75,file=trim(filename_backup), status='unknown',form='unformatted')
	        write(65) nsubdomains, nsubdomains_x, nsubdomains_y
          !write(75) nsubdomains, nsubdomains_x, nsubdomains_y

	       else

                open(65,file=trim(filename), status='unknown',form='unformatted',&
                   position='append')
                !open(75,file=trim(filename_backup), status='unknown',form='unformatted',&
                   !position='append')


	       end if

               call write_statement(.false.)

             end if
	  end do

	end if ! restart_sep

	call task_barrier()

        call t_stopf ('restart_out')

        return
        end


subroutine write_backups()

  use vars
  implicit none
  character *4 rankchar
  character *256 filename
  character *256 filename_backup
  character *10 nstepchar
  integer irank
  integer, external :: lenstr

  call t_startf ('restart_out')

  if(masterproc) then
   write(nstepchar,'(i10)') nstep !new line for writing iteration number in restart filename
   print*,'Writing backup restart file ...'
   !filename =  './RESTART/'//trim(case)//'_'//trim(caseid)//'_misc_restart.bin'
   filename_backup = './RESTART/'//trim(case)//'_'//trim(caseid)//'_misc_restart_'//nstepchar(11-lenstr(nstepchar):10)//'.bin' 
                      !new name with iteration number so restart files aren't overwritten
   
   open(66,file=trim(filename_backup), status='unknown',form='unformatted')
  end if


  if(restart_sep) then
          write(nstepchar,'(i10)') nstep !new line for writing iteration number in restart filename
          write(rankchar,'(i4)') rank

          !filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
          !      rankchar(5-lenstr(rankchar):4)//'_restart.bin'
          filename_backup = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart_'//nstepchar(11-lenstr(nstepchar):10)//'.bin'


          open(65,file=trim(filename_backup), status='unknown',form='unformatted')
          
          write(65) nsubdomains, nsubdomains_x, nsubdomains_y
          !write(75) nsubdomains, nsubdomains_x, nsubdomains_y


    call write_statement(.false.)


  else
    write(nstepchar,'(i10)') nstep !new line for writing iteration number in restart filename
    write(rankchar,'(i4)') nsubdomains
    !filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
    !            rankchar(5-lenstr(rankchar):4)//'_restart.bin'
    filename_backup = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                rankchar(5-lenstr(rankchar):4)//'_restart_'//nstepchar(11-lenstr(nstepchar):10)//'.bin'

    do irank=0,nsubdomains-1
  
       call task_barrier()

       if(irank.eq.rank) then

         if(masterproc) then
        
          open(65,file=trim(filename_backup), status='unknown',form='unformatted')
          write(65) nsubdomains, nsubdomains_x, nsubdomains_y

         else

                open(65,file=trim(filename_backup), status='unknown',form='unformatted',&
                   position='append')
            

         end if

               call write_statement(.false.)

             end if
    end do

  end if ! restart_sep

  call task_barrier()

        call t_stopf ('restart_out')

        return
        end
 
 
 
 
     
  subroutine read_all()

  use vars
  implicit none
  character *4 rankchar
  character *256 filename
  integer irank, ii
  integer, external :: lenstr

  if(masterproc) print*,'Reading restart file ...'

  if(nrestart.ne.2) then
    filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_misc_restart.bin'
  else
    filename = './RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//'_misc_restart.bin'
  end if
  open(66,file=trim(filename), status='unknown',form='unformatted')


	if(restart_sep) then

    write(rankchar,'(i4)') rank

    if(nrestart.ne.2) then
     filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
          rankchar(5-lenstr(rankchar):4)//'_restart.bin'
    else
     filename = './RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//'_'//&
          rankchar(5-lenstr(rankchar):4)//'_restart.bin'
    end if


    open(65,file=trim(filename), status='unknown',form='unformatted')
    read(65)

   call read_statement(.false.)



	else

    write(rankchar,'(i4)') nsubdomains

    if(nrestart.ne.2) then
      filename='./RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart.bin'
    else
      filename='./RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//'_'//&
                  rankchar(5-lenstr(rankchar):4)//'_restart.bin'
    end if
    
    open(65,file=trim(filename), status='unknown',form='unformatted')

    do irank=0,nsubdomains-1

       call task_barrier()

	     if(irank.eq.rank) then

         read (65)
 
         do ii=0,irank-1 ! skip records
                 read(65)
	       end do

         call read_statement(.false.)

        end if

	  end do

	end if ! restart_sep

  call task_barrier()

  dtfactor = -1.

! update the boundaries 
! (just in case when some parameterization initializes and needs boundary points)

  call boundaries(1)
  call boundaries(2)
  call boundaries(3)
  call boundaries(4)

  return
  end



subroutine write_allg()

  use vars
  implicit none
  character *10 nstepchar
  character *256 filename
  integer, external :: lenstr

  call t_startf ('restart_out')

  if(masterproc) then
   print*,'Writing general restart file ...'
   write(nstepchar,'(i10)') nstep
   filename = './RESTART/'//trim(case)//'_'//trim(caseid)//'_'//&
          nstepchar(11-lenstr(nstepchar):10)//'_grestart.bin'
   open(65,file=trim(filename), status='unknown',form='unformatted')
  end if

  call write_statement(.true.)

  call task_barrier()
  call t_stopf ('restart_out')

  return

end


subroutine read_allg()

  use vars
  implicit none
  character *256 filename

  if(masterproc) print*,'Reading general restart file ...'
  filename = './RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//&
                        '_grestart.bin'
  open(65,file=trim(filename), status='old',form='unformatted')
  if(masterproc) print*,'restart file:',trim(filename)

  call read_statement(.true.)

  call task_barrier()

  dtfactor = -1.

  ! update the boundaries 
  ! (just in case when some parameterization initializes and needs boundary points)

  call boundaries(1)
  call boundaries(2)
  call boundaries(3)
  call boundaries(4)

  return

end


  subroutine write_statement(dogenrestart)

  use vars
  use microphysics, only: micro_field, nmicro_fields
  use sgs, only: sgs_field, nsgs_fields, sgs_field_diag, nsgs_fields_diag
  use tracers
  use params
  use movies, only: irecc

  implicit none

  logical, intent(in) ::  dogenrestart
  integer ntape

  if(.not.dogenrestart) then
    write(65)  rank, &
      u, v, w, t, p, qv, qcl, qci, qpl, qpi, dudt, dvdt, dwdt, &
      tracer, micro_field, sgs_field, sgs_field_diag, sstxy, precinst
    close(65)
    if(masterproc) write(66) dt,dx,dy
    ntape = 66
  else  
    if(masterproc) write(65) nx_gl, ny_gl, nzm
    call restart_write_field3D (65,u(1:nx,1:ny,1:nzm),nzm)
    call restart_write_field3D (65,v(1:nx,1:ny,1:nzm),nzm)
    call restart_write_field3D (65,w(1:nx,1:ny,1:nzm),nzm)
    call restart_write_field3D (65,t(1:nx,1:ny,1:nzm),nzm)
    call restart_write_field3D (65,p(1:nx,1:ny,1:nzm),nzm)
    call restart_write_field3D (65,qv(1:nx,1:ny,1:nzm),nzm)
    call restart_write_field3D (65,qcl(1:nx,1:ny,1:nzm),nzm)
    call restart_write_field3D (65,qci(1:nx,1:ny,1:nzm),nzm)
    call restart_write_field3D (65,qpl(1:nx,1:ny,1:nzm),nzm)
    call restart_write_field3D (65,qpi(1:nx,1:ny,1:nzm),nzm)
    call restart_write_field3D (65,dudt(1:nx,1:ny,1:nzm,1:3),nzm*3)
    call restart_write_field3D (65,dvdt(1:nx,1:ny,1:nzm,1:3),nzm*3)
    call restart_write_field3D (65,dwdt(1:nx,1:ny,1:nzm,1:3),nzm*3)
    call restart_write_field3D (65,tracer(1:nx,1:ny,1:nzm,1:ntracers),nzm*ntracers)
    call restart_write_field3D (65,micro_field(1:nx,1:ny,1:nzm,1:nmicro_fields), &
                               nzm*nmicro_fields)
    call restart_write_field3D (65,sgs_field(1:nx,1:ny,1:nzm,1:nsgs_fields), &
                               nzm*nsgs_fields)
    call restart_write_field3D (65,sgs_field_diag(1:nx,1:ny,1:nzm,1:nsgs_fields_diag), &
                               nzm*nsgs_fields_diag)
    call restart_write_field3D (65,sstxy(1:nx,1:ny),1)
    call restart_write_field3D (65,precinst(1:nx,1:ny),1)
    ntape = 65
  end if

  if(masterproc) then
     write(ntape) version, &
      nx, ny, nz, irecc, z, pres, prespot, presi, prespoti, rho, rhow, bet, &
      at, bt, ct, dtn, dt3, time, dz, doconstdz, &
      day, day0, nstep, na, nb, nc, caseid, case, &
      dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo, &
      dolongwave, doshortwave, dosgs, dosubsidence, dotracers,  dosmoke, &
      docoriolis, dosurface, dolargescale,doradforcing, dossthomo, &
      dosfcforcing, doradsimple, donudging_uv, donudging_tq, &
      dowallx, dowally, doperpetual, doseasons, &
      docup, docolumn, soil_wetness, dodynamicocean, ocean_type, &
      delta_sst, depth_slab_ocean, Szero, deltaS, timesimpleocean, &
      pres0, ug, vg, sst_cor, fcor, fcorz, tabs_s, z0, fluxt0, fluxq0, tau0, &
      tauls, tautqls, timelargescale, epsv, nudging_uv_z1, nudging_uv_z2, &
      donudging_t, donudging_q, doisccp, domodis, domisr, dosimfilesout, & 
      dosolarconstant, solar_constant, zenith_angle, notracegases, &
      doSAMconditionals, dosatupdnconditionals, LES_S, &
      nudging_t_z1, nudging_t_z2, nudging_q_z1, nudging_q_z2, &
      ocean, land, sfc_flx_fxd, sfc_tau_fxd, &
      nrad, nxco2, latitude0, longitude0, dofplane, SLM, &
      docoriolisz, doradlon, doradlat, doseawater, salt_factor, &
      ntracers, nmicro_fields, nsgs_fields, nsgs_fields_diag
      close(ntape)
  end if
  if(rank.eq.nsubdomains-1) then
      print *,'Restart file was saved. nstep=',nstep
  endif

  return
  end




  subroutine read_statement(dogenrestart)

  use vars
  use microphysics, only: micro_field, nmicro_fields
  use sgs, only: sgs_field, nsgs_fields, sgs_field_diag, nsgs_fields_diag
  use tracers
  use params
  use movies, only: irecc

  implicit none
  
  real     dx1, dy1
  integer  ntape, nx1, ny1, nz1, rank1, ntr, nmic, nsgs, nsgsd
  character(100) case1,caseid1
  character(10) version1
  logical, intent(in) ::  dogenrestart

  if(.not.dogenrestart) then
    read(65) rank1, &
      u, v, w, t, p, qv, qcl, qci, qpl, qpi, dudt, dvdt, dwdt, &
      tracer, micro_field, sgs_field, sgs_field_diag, sstxy, precinst
    close(65)
    read(66) dt,dx,dy
    ntape = 66
  else
    read(65) nx1, ny1, nz1
    call restart_read_field3D (65,u(1:nx,1:ny,1:nzm),nx1,ny1,nz1)
    call restart_read_field3D (65,v(1:nx,1:ny,1:nzm),nx1,ny1,nz1)
    call restart_read_field3D (65,w(1:nx,1:ny,1:nzm),nx1,ny1,nz1)
    call restart_read_field3D (65,t(1:nx,1:ny,1:nzm),nx1,ny1,nz1)
    call restart_read_field3D (65,p(1:nx,1:ny,1:nzm),nx1,ny1,nz1)
    call restart_read_field3D (65,qv(1:nx,1:ny,1:nzm),nx1,ny1,nz1)
    call restart_read_field3D (65,qcl(1:nx,1:ny,1:nzm),nx1,ny1,nz1)
    call restart_read_field3D (65,qci(1:nx,1:ny,1:nzm),nx1,ny1,nz1)
    call restart_read_field3D (65,qpl(1:nx,1:ny,1:nzm),nx1,ny1,nz1)
    call restart_read_field3D (65,qpi(1:nx,1:ny,1:nzm),nx1,ny1,nz1)
    call restart_read_field3D (65,dudt(1:nx,1:ny,1:nzm,1:3),nx1,ny1,nz1*3)
    call restart_read_field3D (65,dvdt(1:nx,1:ny,1:nzm,1:3),nx1,ny1,nz1*3)
    call restart_read_field3D (65,dwdt(1:nx,1:ny,1:nzm,1:3),nx1,ny1,nz1*3)
    call restart_read_field3D (65,tracer(1:nx,1:ny,1:nzm,1:ntracers),nx1,ny1,nz1*ntracers)
    call restart_read_field3D (65,micro_field(1:nx,1:ny,1:nzm,1:nmicro_fields), &
                               nx1,ny1,nz1*nmicro_fields)
    call restart_read_field3D (65,sgs_field(1:nx,1:ny,1:nzm,1:nsgs_fields), &
                               nx1,ny1,nz1*nsgs_fields)
    call restart_read_field3D (65,sgs_field_diag(1:nx,1:ny,1:nzm,1:nsgs_fields_diag), &
                               nx1,ny1,nz1*nsgs_fields_diag)
    call restart_read_field3D (65,sstxy(1:nx,1:ny),nx1,ny1,1)
    call restart_read_field3D (65,precinst(1:nx,1:ny),nx1,ny1,1)
    ntape = 65
  end if

  read(ntape) version1, &
      nx1, ny1, nz1, irecc, z, pres, prespot, presi, prespoti, rho, rhow, bet, &
      at, bt, ct, dtn, dt3, time, dz, doconstdz, &
      day, day0, nstep, na, nb, nc, caseid1(1:sizeof(caseid)), case1(1:sizeof(case)), &
      dodamping, doupperbound, docloud, doprecip, doradhomo, dosfchomo, &
      dolongwave, doshortwave, dosgs, dosubsidence, dotracers,  dosmoke, &
      docoriolis, dosurface, dolargescale,doradforcing, dossthomo, &
      dosfcforcing, doradsimple, donudging_uv, donudging_tq, &
      dowallx, dowally, doperpetual, doseasons, &
      docup, docolumn, soil_wetness, dodynamicocean, ocean_type, &
      delta_sst, depth_slab_ocean, Szero, deltaS, timesimpleocean, &
      pres0, ug, vg, sst_cor, fcor, fcorz, tabs_s, z0, fluxt0, fluxq0, tau0, &
      tauls, tautqls, timelargescale, epsv, nudging_uv_z1, nudging_uv_z2, &
      donudging_t, donudging_q, doisccp, domodis, domisr, dosimfilesout, &
      dosolarconstant, solar_constant, zenith_angle, notracegases, &
      doSAMconditionals, dosatupdnconditionals, LES_S, &
      nudging_t_z1, nudging_t_z2, nudging_q_z1, nudging_q_z2, &
      ocean, land, sfc_flx_fxd, sfc_tau_fxd, &
      nrad, nxco2, latitude0, longitude0, dofplane, SLM, &
      docoriolisz, doradlon, doradlat, doseawater, salt_factor, &
      ntr, nmic, nsgs, nsgsd
  close(ntape)

  if(version1.ne.version) then
    !print *,'Wrong restart file!'
    print *,'CAUTION: Version of SAM that wrote the restart files:',version1
    print *,'Current version of SAM',version
    !call task_abort()
  end if
  if(.not.dogenrestart) then
    if(rank.ne.rank1) then
       print *,'Error: rank of restart data is not the same as rank of the process'
       print *,'rank1=',rank1,'   rank=',rank
       call task_abort()
    endif
    if(nx.ne.nx1.or.ny.ne.ny1.or.nz.ne.nz1) then
       print *,'Error: domain dims (nx,ny,nz) set by grid.f'
       print *,' not correspond to ones in the restart file.'
       print *,'in executable:   nx, ny, nz:',nx,ny,nz
       print *,'in restart file: nx, ny, nz:',nx1,ny1,nz1
       print *,'Exiting...'
       call task_abort()
    endif
  end if
  if(nmic.ne.nmicro_fields) then
     print*,'Error: number of micro_field in restart file is not the same as nmicro_fields'
     print*,'nmicro_fields=',nmicro_fields,'   in file=',nmic
     print*,'Exiting...'
     call task_abort()
  end if
  if(nsgs.ne.nsgs_fields.or.nsgsd.ne.nsgs_fields_diag) then
     print*,'Error: number of sgs_field in restart file is not the same as nsgs_fields'
     print*,'nsgs_fields=',nsgs_fields,'   in file=',nsgs
     print*,'nsgs_fields_diag=',nsgs_fields_diag,'   in file=',nsgsd
     print*,'Exiting...'
     call task_abort()
  end if
  if(ntr.ne.ntracers) then
     print*,'Error: number of tracers in restart file is not the same as ntracers.'
     print*,'ntracers=',ntracers,'   ntracers(in file)=',ntr
     print*,'Exiting...'
     call task_abort()
  end if
  close(65)
  if(rank.eq.nsubdomains-1) then
     print *,'Case:',caseid
     print *,'Restarting at step:',nstep
     print *,'Time(s):',nstep*dt
  endif

  return
  end
 

subroutine restart_write_field3D(ntape,f,nlevs)

! write full global arrays into a general
! (not dependent on number of PEs) restart file

use grid, only: masterproc, nx, ny, nx_gl, ny_gl, rank, nsubdomains
implicit none
! Input:
integer ntape
real f(nx,ny,nlevs)
integer nlevs,req
! Local:
real fld(nx_gl,ny_gl), buf(nx,ny)
integer irank,ranksource,tagsource,ii,jj,k

if(masterproc) write(ntape) nlevs
do k=1,nlevs
   if(masterproc) fld(1:nx,1:ny) = f(:,:,k)
   do irank = 1, nsubdomains-1
       call task_barrier()
       if(irank.eq.rank) then
         call task_bsend_float(0,f(:,:,k),nx*ny,irank)
       end if
       if(masterproc) then
         call task_receive_float(buf(:,:),nx*ny,req)
         call task_wait(req,ranksource,tagsource)
         call task_rank_to_index(ranksource,ii,jj)
         fld(ii+1:ii+nx,jj+1:jj+ny) = buf(:,:)
       end if
   end do
   if(masterproc) write(ntape) fld
end do

end



subroutine restart_read_field3D(ntape,f,nxr,nyr,nlevs)

! write full global arrays into a general
! (not dependent on number of PEs) restart file

use params, only: nstepgrestart
use grid, only: nx, ny, nx_gl, ny_gl, rank, dx, dy
implicit none
! Input:
integer ntape
real f(nx,ny,nlevs)
integer nxr,nyr ! size of arrays in restart file
integer nlevs
! Local:
real,allocatable ::  fld(:,:)
integer nlevs1,ii,jj,iii,jjj,k,i,j,i1,i2,j1,j2
real f11,f12,f21,f22,ax,by

allocate(fld(nxr,nyr))

call task_rank_to_index(rank,iii,jjj)

read(ntape) nlevs1
if(nlevs.ne.nlevs1) then
  print*,'Restart_read_field3D: error: nlevs.ne.nlevs1. nlevs=',nlevs,' nlevs1=',nlevs1,' quitting...'
  call task_abort
end if

if (nstepgrestart.eq.-1) then !Special case which replicates the PERIODIC solution to increase domain size

  if (nxr.le.nx_gl.or.nyr.le.ny_gl) then

    if ((mod(nx_gl,nxr).ne.0).or.(mod(ny_gl,nyr).ne.0)) then
      print*,'Restart_read_field3D: error: the periodic replication requires proportional scaling, quitting...'
      call task_abort
    end if
    

    do k=1,nlevs
      read(ntape) fld

      do j=1,ny
       do i=1,nx

         i1 = mod(iii+i,nxr)
         if (i1.eq.0) i1 = nxr
         j1 = mod(jjj+j,nyr)
         if (j1.eq.0) j1 = nyr

         f(i,j,k)=fld(i1,j1)
       end do
      end do

    end do

  else

    print*,'Restart_read_field3D: error: doing a special periodic general restart with fewer elements, quitting...'
    call task_abort

  end if

else !usual grid refinement or coarsening

  if (nxr.le.nx_gl.or.nyr.le.ny_gl) then !grid refinement

    ii=nx_gl/nxr
    jj=ny_gl/nyr

    do k=1,nlevs
      read(ntape) fld

      do j=1,ny
       do i=1,nx
         i1 = iii/ii+(i-1)/ii+1
         i2 = iii/ii+(i+1)/ii+1
         if (i2.gt.nxr) i2 = i2-nxr
         j1 = jjj/jj+(j-1)/jj+1
         j2 = jjj/jj+(j+1)/jj+1
         if (j2.gt.nyr) j2 = j2-nyr

         ax = (1.0-mod(iii/ii+i-1,ii)/real(ii,8))
         by = (1.0-mod(jjj/jj+j-1,jj)/real(jj,8))

         f11 = fld(i1,j1)
         f21 = fld(i2,j1)
         f12 = fld(i1,j2)
         f22 = fld(i2,j2)

         f(i,j,k)=by*(ax*f11+(1.0-ax)*f21)+(1.0-by)*(ax*f12+(1.0-ax)*f22)

       end do
      end do

    end do

  else !grid coarsening

    ii=nxr/nx_gl
    jj=nyr/ny_gl

    do k=1,nlevs
      read(ntape) fld

      f(1:nx,1:ny,k) = fld(iii*ii+1:(iii+nx)*ii:ii,jjj*jj+1:(jjj+ny)*jj:jj)

    end do

  end if

end if


deallocate(fld)

end

