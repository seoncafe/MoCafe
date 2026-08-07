program shg_bench
   !--- Camps et al. (2015, A&A 580, A87) stochastic-heating dust-emission
   !--- benchmark, run through the SEDust library with the same
   !--- Zubko et al. (2004) BARE-GR-S dust model the benchmark uses.  For each
   !--- Mathis-ISRF scaling U we compute the emission spectrum and write it so
   !--- it can be overlaid on the published multi-code reference (CRT, DIRTY,
   !--- SKIRT, DustEM, ...).  Validates the emission engine, not the RT.
   use constants, only: wp
   use radfield,  only: J_Mathis
   use dust_lib,  only: dust_model_t, build_dust, dust_emission, dust_nlam, dust_lambda
   implicit none
   !--- One data directory is all this names, the same par%sed_datadir a MoCafe
   !--- run names: build_dust resolves the Zubko model from SEDust/data/zubko/
   !--- under it.  The path is relative to this directory, so a fresh checkout
   !--- runs the benchmark as it stands.
   character(len=*), parameter :: DATA = '../../SEDust/data'
   type(dust_model_t)    :: m
   real(wp), allocatable :: J(:), total(:), lam(:)
   real(wp) :: Ulist(9), U
   integer  :: iu, il, unit
   character(len=64) :: fname
   character(len=8)  :: utag

   Ulist = [1.0e-2_wp, 1.0e-1_wp, 1.0e0_wp, 1.0e1_wp, 1.0e2_wp, &
            1.0e3_wp, 1.0e4_wp, 1.0e5_wp, 1.0e6_wp]

   !--- Whatever SEDust serves for 'zubko' is what is measured.  Which optics
   !--- that model is made of -- the Zubko et al. tables as distributed, or
   !--- SEDust's own recomputation of them -- is SEDust's decision, not this
   !--- code's, and reproducing the published model is SEDust's business to get
   !--- right.  SEDust now makes it a named argument, zubko_optics, whose default
   !--- 'zda' is the distributed tables the seven codes read; this program leaves
   !--- it unset, so what it measures is the library's own default, which is also
   !--- what a MoCafe run with par%dust_model = 'zubko' gets.
   !---
   !--- include_euv = .true. keeps the whole 1201-point axis of the Zubko optics
   !--- tables, from 1.0e-3 um, which is the range the benchmark's own reference
   !--- spectra cover.  The non-ionizing view would cut it at the Lyman limit --
   !--- right for a transport run, whose field stops there, but not for this.
   call build_dust(m, 'zubko', DATA, 300, 2.7_wp, 5.0e3_wp, include_euv=.true.)
   lam = dust_lambda(m)
   allocate(J(dust_nlam(m)), total(dust_nlam(m)))
   print '(a,i0)', 'SEDust zubko model built, NLAM = ', dust_nlam(m)

   do iu = 1, size(Ulist)
      U = Ulist(iu)
      call J_Mathis(U, lam, J)
      call dust_emission(m, J, total)
      write(utag,'(es8.0e2)') U
      write(fname,'(a,es7.0e2,a)') 'sedust_zubko_Mathis_U_', U, '.dat'
      open(newunit=unit, file=trim(adjustl(fname)), status='replace')
      write(unit,'(a)') '# lambda[um]   lambda*emission (SEDust zubko BARE-GR-S, per H)'
      do il = 1, dust_nlam(m)
         write(unit,'(es14.6,es16.7)') lam(il), total(il)
      end do
      close(unit)
      print '(a,es9.1,a,es12.4,a,f7.1,a)', '  U=', U, ' : peak lamI=', &
         maxval(total), ' at lam=', lam(maxloc(total,1)), ' um'
   end do
end program shg_bench
