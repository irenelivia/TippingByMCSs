! Initialize 2D output

subroutine stat_2Dinit(ResetStorage)

use vars
implicit none

integer :: ResetStorage

     prec_xy(:,:) = 0.
     shf_xy(:,:) = 0.
     lhf_xy(:,:) = 0.
     lwnt_xy(:,:) = 0.
     swnt_xy(:,:) = 0.
     lwntc_xy(:,:) = 0.
     swntc_xy(:,:) = 0.
     pw_xy(:,:) = 0.
     cw_xy(:,:) = 0.
     iw_xy(:,:) = 0.
     cld_xy(:,:) = 0.
     u200_xy(:,:) = 0.
     v200_xy(:,:) = 0.
     usfc_xy(:,:) = 0.
     vsfc_xy(:,:) = 0.
     w500_xy = 0.
     lwns_xy(:,:) = 0.
     swns_xy(:,:) = 0.
     solin_xy(:,:) = 0.
     lwnsc_xy(:,:) = 0.
     swnsc_xy(:,:) = 0.
     qocean_xy(:,:) = 0.

!===================================
! RF ADDITIONS: SAMPLING AT Z-LEVELS

     qt1300_xy = 0.
     t1300_xy = 0.
     p1300_xy = 0.
     u1300_xy = 0.
     v1300_xy = 0.
     w1300_xy = 0.
     qt3500_xy = 0.
     t3500_xy = 0.
     p3500_xy = 0.
     u3500_xy = 0.
     v3500_xy = 0.
     w3500_xy = 0.
     wsrfc_xy = 0.
     qvsrfc_xy = 0.
     tsrfc_xy = 0.


! END RF ADDITIONS

!===================================
! UW ADDITIONS: MOSTLY 2D STATISTICS

     if(ResetStorage.eq.1) then
       !bloss: store initial profiles for computation of storage terms in budgets
       ustor(:) = u0(1:nzm)
       vstor(:) = v0(1:nzm)
       tstor(:) = t0(1:nzm)
       qstor(:) = q0(1:nzm)

       utendcor(:) = 0.
       vtendcor(:) = 0.
     end if

     psfc_xy(:,:) = 0.

     u850_xy(:,:) = 0.
     v850_xy(:,:) = 0.

     swvp_xy(:,:) = 0.

     cloudtopheight(:,:) = 0.
     echotopheight(:,:) = 0.
     cloudtoptemp(:,:) = 0.

! END UW ADDITIONS

end

