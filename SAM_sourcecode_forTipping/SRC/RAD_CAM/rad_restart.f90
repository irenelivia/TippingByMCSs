	subroutine write_rad()
	
	use rad
        use radae, only: abstot_3d, absnxt_3d, emstot_3d
	implicit none
	character *4 rankchar
  character *256 filename
  character *256 filename_backup
  character *10 nstepchar
	integer irank
        integer lenstr
        external lenstr

        if(masterproc) print*,'Writting radiation restart file CAM...'

        if(restart_sep) then
          write(rankchar,'(i4)') rank
          open(56,file='./RESTART/'//case(1:len_trim(case))//'_'// &
               caseid(1:len_trim(caseid))//'_'// &
               rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
               status='unknown',form='unformatted')
          !open(66,file='./RESTART/'//case(1:len_trim(case))//'_'// &
          !     caseid(1:len_trim(caseid))//'_'// &
          !     rankchar(5-lenstr(rankchar):4)//'_restart_rad'// &
          !     nstepchar(11-lenstr(nstepchar):10)//'.bin', &
          !     status='unknown',form='unformatted')
          
          write(56) nsubdomains
          !write(66) nsubdomains
          
	        write(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qv_rad,qs_rad, &
                    cld_rad, rel_rad, rei_rad, res_rad, &
                    swdsvisxy,swdsnirxy,swdsvisdxy,swdsnirdxy,coszrsxy,lwdsxy,swdsxy, &
	                  qrad,absnxt_3d,abstot_3d,emstot_3d
          !write(66) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qv_rad,qs_rad, &
                    cld_rad, rel_rad, rei_rad, res_rad, &
                    swdsvisxy,swdsnirxy,swdsvisdxy,swdsnirdxy,coszrsxy,lwdsxy,swdsxy, &
                    qrad,absnxt_3d,abstot_3d,emstot_3d 
          close(56)
          !close(66)

        else !

          write(rankchar,'(i4)') nsubdomains
          write(nstepchar,'(i10)') nstep 
          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               if(masterproc) then

                filename='./RESTART/'//case(1:len_trim(case))//'_'// &
                      caseid(1:len_trim(caseid))//'_'// &
                      rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin'
                !filename_backup = './RESTART/'//case(1:len_trim(case))//'_'// &
                !                   caseid(1:len_trim(caseid))//'_'// &
                !                   rankchar(5-lenstr(rankchar):4)//'_restart_rad_'// &
                !                   nstepchar(11-lenstr(nstepchar):10)//'.bin'

                  open(56,file= trim(filename),status='unknown',form='unformatted')
                  !open(66,file= trim(filename_backup),status='unknown',form='unformatted')
                  write(56) nsubdomains
                  !write(66) nsubdomains

               else
                  filename = './RESTART/'//case(1:len_trim(case))//'_'// & 
                      caseid(1:len_trim(caseid))//'_'// &
                      rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin'
                  !filename_backup = './RESTART/'//case(1:len_trim(case))//'_'// & 
                  !                   caseid(1:len_trim(caseid))//'_'// &
                  !                   rankchar(5-lenstr(rankchar):4)//'_restart_rad_'// &
                  !                   nstepchar(11-lenstr(nstepchar):10)//'.bin'
                  open(56,file= file= trim(filename),status='unknown',form='unformatted', position='append')
                  !open(66,file= file= trim(filename_backup),status='unknown',form='unformatted', position='append')

               end if

	       write(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qv_rad,qs_rad,  &
                         cld_rad, rel_rad, rei_rad, res_rad, &
                    swdsvisxy,swdsnirxy,swdsvisdxy,swdsnirdxy,coszrsxy,lwdsxy,swdsxy, &
	      	   qrad,absnxt_3d,abstot_3d,emstot_3d 
               close(56)
        !write(66) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qv_rad,qs_rad,  &
        !                 cld_rad, rel_rad, rei_rad, res_rad, &
        !            swdsvisxy,swdsnirxy,swdsvisdxy,swdsnirdxy,coszrsxy,lwdsxy,swdsxy, &
        !     qrad,absnxt_3d,abstot_3d,emstot_3d 
        !       close(66)
           end if
        end do

        end if ! restart_sep

	if(masterproc) then
           print *,'Saved radiation restart file. nstep=',nstep
	endif

        call task_barrier()

        return
        end
 
 
 
 
     
	subroutine read_rad()
	
	use rad
        use radae, only: abstot_3d, absnxt_3d, emstot_3d
	implicit none
	character *4 rankchar
	integer irank,ii
        integer lenstr
        external lenstr
	
        if(masterproc) print*,'Reading radiation restart file...'

        if(restart_sep) then

          write(rankchar,'(i4)') rank

          if(nrestart.ne.2) then
            open(56,file='./RESTART/'//trim(case)//'_'//trim(caseid)//'_'// &
              rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
              status='unknown',form='unformatted')
          else
            open(56,file='./RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//'_'// &
              rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
              status='unknown',form='unformatted')
          end if
          read (56)
	  read(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qv_rad,qs_rad,  &
                         cld_rad, rel_rad, rei_rad, res_rad, &
                    swdsvisxy,swdsnirxy,swdsvisdxy,swdsnirdxy,coszrsxy,lwdsxy,swdsxy, &
	     qrad,absnxt_3d,abstot_3d,emstot_3d 
          close(56)

        else

          write(rankchar,'(i4)') nsubdomains

          if(nrestart.ne.2) then
            open(56,file='./RESTART/'//trim(case)//'_'//trim(caseid)//'_'// &
              rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
              status='unknown',form='unformatted')
          else
            open(56,file='./RESTART/'//trim(case_restart)//'_'//trim(caseid_restart)//'_'// &
              rankchar(5-lenstr(rankchar):4)//'_restart_rad.bin', &
              status='unknown',form='unformatted')
          end if

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               read (56)

               do ii=0,irank-1 ! skip records
                 read(56)
               end do

	       read(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qv_rad,qs_rad, &
                         cld_rad, rel_rad, rei_rad, res_rad, &
                    swdsvisxy,swdsnirxy,swdsvisdxy,swdsnirdxy,coszrsxy,lwdsxy,swdsxy, &
	  	 qrad,absnxt_3d,abstot_3d,emstot_3d 
               close(56)
             end if

          end do

        end if ! restart_sep

        if(rank.eq.nsubdomains-1) then
             print *,'Case:',caseid
             print *,'Restart radiation at step:',nstep
             print *,'Time:',nstep*dt
        endif

        call task_barrier()


        return
        end
