! Add PBHs to a Gadget format 0 IC file
#define DEBUG
!#define OUTPBH
#define PERTURB
!#define STREAMING ! not relevant for this
!#define FLIPID
#define PBH
#define DAMPING
#define VELOCITY
#define TRUNCATE
!#define CENTER_PBH
#define ALLPBH
#define PERIOD
! change SELFS
! change NONLINEAR
! undefine both -> strong perturbation
!#define SELFS ! recommended for the matter-dominated regime!
!#define NONLINEAR ! initial nonlinear structures around PBHs
!#define PERTURB_PBH ! perturbations on PBHs by themselves
!#define MASSFUNC_LOGNORMAL
! Extended (log-normal) PBH mass function instead of a monochromatic one.
! When undefined the code is bit-for-bit identical to the monochromatic version.
! Parameterisation (see Mc_pbh / sigma_ln below):
!   ln M ~ Normal(ln Mc_pbh, sigma_ln**2),  <M> = Mc_pbh*exp(sigma_ln**2/2)
! <M> replaces mpbh in the PBH number normalisation; fpbh keeps its meaning as
! the mass fraction of dark matter in PBHs.

program pbh_ic

   use random_object
   use kdtree

   implicit none

   character*200, parameter :: path='./'
   !character*200, parameter :: icbase='ics_N128L014_nu2_zoom_lowres'
   !character*200, parameter :: icbase='ics_N128L014_zoom_lowres'
   !character*200, parameter :: icbase='ics_N128L1'
   !character*200, parameter :: icbase='ics_N256L10'
   character*200, parameter :: icbase='ics_N32L025_nu3'
   !character*200, parameter :: icbase='ics_N32L025_nu3_dmo_z1e6'
   !character*200, parameter :: icbase='ics_N32L5_dmo'
   !character*200, parameter :: icbase='ics_N32L5_dmo'

   character*200 filename, outname

   integer*4 npart(0:5), nall(0:5)
   real*8    massarr(0:5)
   real*8    a,redshift,boxsize
   real*8    Om0,Ol,h
   integer*4 unused(24)
   integer*4 fn,i,j,k,l,ipbh,nstart,flag_sfr,flag_feedback,flag_cool,nfs
   integer*4 Ngas,Ndm,N,Nm,Nhres,i_dummy,rand_seed,n_refine,nonlinear


   real*4,allocatable    :: pos(:,:),vel(:,:),lmass(:),iegas(:)
   integer*4,allocatable :: id(:), idpbh(:)
   real*4,allocatable    :: pos1(:,:),vel1(:,:),lmass1(:)
   integer*4,allocatable :: id1(:)

   real*4,allocatable    :: phaseDM(:,:), lpbh(:,:), ldpbh(:,:) !, lpbh0(:,:)
#ifdef MASSFUNC_LOGNORMAL
   ! 2-3: per-PBH mass in PHYSICAL Msun, allocated/filled in lockstep with lpbh.
   ! Slots beyond Ngen (the "parked at 10*bs" padding) are left at 0.0, which
   ! is a neutral placeholder: it contributes nothing to any force sum.
   real*4,allocatable    :: mpbh_list(:)
   real*8 :: mm_avg, mm_std, lnm_avg, lnm_std
   integer*4 :: Nm_hres
#endif
   real*4 :: p_dummy(6), dx(3), fac, dx_pert, rpert, fac_pert, disp(3), bs, esft, rpert0
#ifdef VELOCITY
   real*4 :: dv_pert
#endif
   type(kd_root) :: kdroot, pbhroot
   integer*4 :: pbh_type, nx, nxdm, nnn, npert, Ncell, Ngen, Nthis, xind(3) !, Ngen0
   integer*4, allocatable :: nnindex(:) !, n_dummy
   real*4, allocatable :: nndis(:), dmgrid(:,:,:)
   real*4 :: dismax, vol, re(3), le(3), rvec(3), xdis(3), ftrunc, H0, mref, dpertmax
   real*4 :: dummy, facavg0, facavg1, fshrink, d_dummy
   real*4 :: delta, fscale, fscalev!, rnearest, dnearest(3)
   type(rng) :: rng_obj
   
   real*4 :: msun, kpc, gra, PI
   parameter(msun=1.989e33, kpc=3.085678e21, gra=6.672e-8, PI=3.1415926536)
   
   real*8 :: mpbh, sigma_rms, z_rec, fstream!, a_target
   real*8 :: mpbh_mean   ! <M>; equals mpbh in the monochromatic case
   ! Realised (not target) PBH mass budget used to deplete massarr(1); see the
   ! "mass-conservation fix" block below. All real*8, no real*4 intermediates.
   real*8 :: m_dm_msun, M_tot_msun, f_real, mass_removed_msun
#ifdef MASSFUNC_LOGNORMAL
   ! 2-2: log-normal mass function parameters.
   !   Mc_pbh   median mass [Msun]
   !   sigma_ln standard deviation of ln M
   ! Mc_pbh=1e3, sigma_ln=1.0 puts 95% of the PBHs (+-2 sigma) in
   !   1e3*exp(-2) .. 1e3*exp(+2) = 135 .. 7389 Msun,
   ! and gives <M> = 1e3*exp(1.0**2/2) = 1648.7 Msun.
   real*8 :: Mc_pbh, sigma_ln
   parameter(Mc_pbh=1.0d3, sigma_ln=1.0d0)
#endif
   real*4 :: fpbh, Npbh, Or0, aeq, vrel(3), fc, Ob0, a_, gfac, gfac_gas, vfac!, facv
   parameter(sigma_rms=30.0, z_rec=1100.0, fstream=0.8)!, a_target=1.0/31.0)
   !parameter(nnn=1, pbh_type=5, mpbh=1e10, fpbh=1e-4)
   parameter(nnn=1, pbh_type=3, mpbh=1e6, fpbh=1e-3)
   !parameter(nnn=1, pbh_type=3, mpbh=1e6, fpbh=1e-2)
   !parameter(nnn=1, pbh_type=3, mpbh=1e9, fpbh=1e-2)
   

   logical :: zoomin
   parameter(zoomin=.false., rand_seed=2333, n_refine=1) !use rand_seed=2333 for 10 PBH with fpbh=1e-3, mpbh=1e10 in a 10 Mpc/h box
   parameter(esft=1.0/30.0, Or0=9.117e-5, Ob0=0.04864)
   parameter(ftrunc=1.0, fac_pert=100.0, rpert0=0)!300)
   parameter(delta=18.0*PI**2, fscale=1.0, fscalev=1.0)
   !parameter(delta=18.0*PI**2, fscale=0.16260587661983628, fscalev=1.0)
   !parameter(delta=18.0*PI**2, fscale=1, fscalev=0.16260587661983628)!0.21138763960578716)
   
#ifdef TEST_MASSFUNC
! Mass-sampler self-test: draws Ntest masses with the exact same ran3/Box-Muller
! routine and (Mc_pbh, sigma_ln) used by the real PBH generation loop below,
! and reports the realised statistics against the analytic targets. Runs
! standalone (own rng_obj%init) before any IC file is touched, so it works
! without a large Ngen and without needing an actual snapshot on disk.
#ifndef MASSFUNC_LOGNORMAL
   write(0,*) 'TEST_MASSFUNC requires MASSFUNC_LOGNORMAL to also be defined', &
     & ' (it tests that sampler''s parameters, Mc_pbh/sigma_ln).'
   stop 1
#else
   block
     integer*4, parameter :: Ntest = 1000000
     real*8 :: lnm, s_lnm, s_lnm2, s_m
     integer*4 :: itest
     call rng_obj%init(rand_seed)
     s_lnm  = 0.0d0
     s_lnm2 = 0.0d0
     s_m    = 0.0d0
     do itest=1,Ntest
       lnm = log(Mc_pbh) + sigma_ln*dble(gauss_bm(rng_obj))
       s_lnm  = s_lnm  + lnm
       s_lnm2 = s_lnm2 + lnm**2
       s_m    = s_m    + exp(lnm)
     enddo
     write(0,*) 'TEST_MASSFUNC: N=',Ntest,' Mc_pbh=',Mc_pbh,' sigma_ln=',sigma_ln
     write(0,*) '  mean ln M =', s_lnm/dble(Ntest), &
       & '  (target ln Mc_pbh =', log(Mc_pbh), ')'
     write(0,*) '  std  ln M =', sqrt(s_lnm2/dble(Ntest) - (s_lnm/dble(Ntest))**2), &
       & '  (target sigma_ln  =', sigma_ln, ')'
     write(0,*) '  mean     M =', s_m/dble(Ntest), &
       & '  (target Mc_pbh*exp(sigma_ln**2/2) =', Mc_pbh*exp(0.5d0*sigma_ln**2), ')'
   end block
   stop
#endif
#endif

   npert=1

   allocate(nnindex(max(nnn,npert)))
   allocate(nndis(max(nnn,npert)))

   filename= path(1:len_trim(path)) // icbase(1:len_trim(icbase)) // '_cdm.dat'
   outname= path(1:len_trim(path)) // icbase(1:len_trim(icbase))
   print *,'opening...  '//filename(1:len_trim(filename))
   ! read in the header
   Nm = 0
   open (1, file=filename, form='unformatted')
   read (1) npart, massarr ,a, redshift, flag_sfr,flag_feedback, nall, flag_cool, nfs, boxsize, Om0, Ol, h, unused
     print *,'Redshift=',redshift,'h=',h,'Om=',Om0
     do i=0,5
        if (massarr(i) .eq. 0) then
          Nm = Nm + nall(i)
        endif
        print *,'Type=',i,' Particles=', nall(i), 'mass=', massarr(i)
     end do
#ifdef MASSFUNC_LOGNORMAL
     ! 2-5: number of entries the ORIGINAL lmass block contributes for types
     ! 0/1.  The PBHs are inserted into pos1/vel1/id1 right after the type 0/1
     ! particles, so in the variable-mass block their masses must go right
     ! after these Nm_hres entries and before those of types >= 2.
     Nm_hres = 0
     do i=0,1
        if (massarr(i) .eq. 0) Nm_hres = Nm_hres + nall(i)
     end do
#endif
     Ngas = nall(0)
     Ndm = nall(1)
     Nhres = Ngas+Ndm
     N=sum(npart)
     ! read in data
     !write(0,*) Ngas, N
      allocate(pos(1:3,1:N))
      allocate(vel(1:3,1:N))
      allocate(id(1:N))
      allocate(idpbh(1:N))
      if (Ngas .gt. 0) then
        allocate(iegas(1:Ngas))
      endif
      read (1) pos
      read (1) vel
      read (1) id
      if (Nm .gt. 0) then
        allocate(lmass(1:Nm))
        read (1) lmass
      endif
      if (Ngas .gt. 0) then
        read (1) iegas
      endif
   close (1)
   print *,'Done with reading.'

  allocate(phaseDM(1:6,1:npart(1)))
  phaseDM(1:3,:) = pos(1:3,Ngas+1:Ngas+Ndm)
  phaseDM(4:6,:) = vel(1:3,Ngas+1:Ngas+Ndm)
  
  if(zoomin) then
    vol = 1.0
    do i=1,3
      le(i) = minval(phaseDM(i,:))
      re(i) = maxval(phaseDM(i,:))
      vol = vol*(re(i)-le(i))
    enddo
  else
    le = 0.0
    re = boxsize
    vol = boxsize**3
  endif
  
#ifdef STREAMING
  if (Ngas .gt. 0) then
	call rng_obj%init(rand_seed*2) 
    rvec(1) = rng_obj%get_ran()
    rvec(2) = rng_obj%get_ran()
    rvec(3) = rng_obj%get_ran()
    vrel = rvec/sum(rvec**2)**0.5*sigma_rms*(1.0+redshift)/(1.0+z_rec) *fstream/a**0.5
    write(0,*) 'Add DM-gas relative velocity: ',vrel
    do i=1,Ngas
      vel(1:3,i) = vel(1:3,i) + vrel
    enddo
  endif
#endif
  
  call rng_obj%init(rand_seed) 
  Ngen = 0
  bs=real(boxsize)
  
#ifdef PBH
  ! generate a grid in which each cell contains ~ 1 PBH on average
  nxdm = int(float(Ndm)**(1.0/3.0)+0.5)
  dismax = (vol/float(nall(1)))**(1.0/3.0)
  if(zoomin) then
    le = le-dismax*0.5 *(1.0+1.0/(float(nxdm)-1.0))
    re = re+dismax*0.5 *(1.0+1.0/(float(nxdm)-1.0))
    vol = 1.0
    do i=1,3
      re(i) = min(bs, re(i))
      le(i) = max(0.0, le(i))
      vol = vol*(re(i)-le(i))
    enddo
    dismax = (vol/float(nall(1)))**(1.0/3.0)
  endif
  write(0,*) 'Grid spacing = ', dismax
  write(0,*) 'Target region extent: le=', le,'re=', re
  
  ! 2-2: the PBH number is set by the total PBH mass budget fpbh*M_dm divided
  ! by the MEAN PBH mass.  For a monochromatic function that is mpbh itself.
#ifdef MASSFUNC_LOGNORMAL
  mpbh_mean = Mc_pbh*exp(0.5d0*sigma_ln**2)
#else
  mpbh_mean = mpbh
#endif
  if (Ngas .gt. 0) then
    Npbh = float(Ndm) * massarr(1)*1e10*fpbh/h/mpbh_mean
  else
    Npbh = float(Ndm) * massarr(1)*1e10*fpbh/h/mpbh_mean * (Om0-Ob0)/Om0
  endif
  if (rpert0 .gt. 0) then
    rpert = rpert0
  else
    rpert = fac_pert * (vol/Npbh)**(1.0/3.0)
  endif

  aeq = Or0/Om0
  H0 = 100*h*1e5/1e3/kpc
  fc = (Om0-Ob0)/Om0
  a_ = 0.25*(sqrt(1+24*fc)-1)
  fshrink = 1-(delta)**(-1.0/3.0)
  
  gfac = Dgrow(a/aeq, fc, a_) * fscale
  gfac_gas = Dgrow(a*(1+z_rec), fc, a_) * fscale
  vfac = Dgrowv(a, aeq, fc, a_) * fscalev
  
!#ifdef IHALO
  ! diagnostic only; uses the mean PBH mass so that it stays meaningful for an
  ! extended mass function (mpbh_mean == mpbh in the monochromatic case)
  dummy = (8*PI/3.0*gra*mpbh_mean*msun * gfac/(Om0*H0**2))**(1.0/3.0) * h/kpc
!  dpertmax = dummy !*(4*PI/3.0)*(Dgrow(a/aeq, fc, a_)/Dgrow(a_target/aeq, fc, a_))**(2.0/3.0)
  write(0,*) 'Nonlinear scale = ', dummy, dummy/dismax, gfac!, dummy
!  dpertmax = max(dpertmax, dismax)
!#else
  dpertmax = dismax
!#endif

  ! make a kdtree of DM particles
  call kd_init(kdroot, phaseDM, bs)
  write(0,*) 'Done KD tree construction'

#ifdef CENTER_PBH
  Ngen = 1
  Ncell = 1
  write(0,*) 'Add one PBH at the center'
  allocate(lpbh(1:6,2*Ncell))
#ifdef MASSFUNC_LOGNORMAL
  allocate(mpbh_list(2*Ncell))
  mpbh_list = 0.0                       ! neutral placeholder for unused slots
  mpbh_list(1) = Mc_pbh*exp(sigma_ln*gauss_bm(rng_obj))
#endif
  do i=1,3
    lpbh(i,1) = 0.5*(re(i)+le(i))
  enddo
  call kd_nnearest(kdroot, lpbh(1:3,1), nnn, nnindex, nndis, i_dummy)
  p_dummy = 0.0
  do l=1,nnn
    p_dummy(4:6) = p_dummy(4:6) + phaseDM(4:6,nnindex(l))
  enddo
  p_dummy(4:6) = p_dummy(4:6)/float(nnn)
  lpbh(4:6,1) = p_dummy(4:6) !phaseDM(4:6,nnindex(1))
#else  
  nx = int(Npbh**(1.0/3.0)+1.0)*n_refine
  Ncell = float(nx)**3
  write(0,*) '1D grid num.:',nx
  write(0,*) 'Num. of cells for PBH sampling:',nx**3
  write(0,*) 'PBH num. per cell:', Npbh/Ncell
  write(0,*) 'Num. of PDM particles per cell:', float(Ndm)/float(Ncell)
  
  ! mapping DM particles to a grid
  allocate(dmgrid(nx, nx, nx))
  dmgrid = 0.0
  dx = (re-le)/float(nx)
  do i=1,Ndm
    do j=1,3
      xind(j) = int((phaseDM(j,i)-le(j))/dx(j))+1
    enddo
    dmgrid(xind(1), xind(2), xind(3)) = dmgrid(xind(1), xind(2), xind(3)) + 1.0
  enddo
  
  ! sample PBHs
  write(0,*) 'Adding PBHs...'
  allocate(lpbh(1:6,2*Ncell))
#ifdef MASSFUNC_LOGNORMAL
  allocate(mpbh_list(2*Ncell))
  mpbh_list = 0.0                       ! neutral placeholder for unused slots
#endif
  do i=1,nx
    do j=1,nx
      do k=1,nx
        if (dmgrid(i,j,k) .gt. 0) then
        Nthis = poisson(dmgrid(i,j,k)/float(Ndm)*Npbh, rng_obj)
        do ipbh=1,Nthis
          rvec(1) = rng_obj%get_ran()
          rvec(2) = rng_obj%get_ran()
          rvec(3) = rng_obj%get_ran()
          lpbh(1:3,Ngen+ipbh) = rvec*dx + dx*(/float(i-1), float(j-1), float(k-1)/) + le
          call kd_nnearest(kdroot, lpbh(1:3,Ngen+ipbh), nnn, nnindex, nndis, i_dummy)
          p_dummy = 0.0
          do l=1,nnn
            p_dummy(4:6) = p_dummy(4:6) + phaseDM(4:6,nnindex(l))
          enddo
          p_dummy(4:6) = p_dummy(4:6)/float(nnn)
          lpbh(4:6,Ngen+ipbh) = p_dummy(4:6) !phaseDM(4:6,nnindex(1))
#ifdef MASSFUNC_LOGNORMAL
          ! 2-3: independent mass draw for this PBH, M = Mc*exp(sigma*z)
          ! with z a standard normal from the Box-Muller transform.
          ! NOTE: the spatial sampling above is untouched; mass and position
          ! are drawn independently of each other.
          mpbh_list(Ngen+ipbh) = Mc_pbh*exp(sigma_ln*gauss_bm(rng_obj))
#endif
        enddo
        Ngen = Ngen + Nthis
        !else
        !  write(0,*) 'empty cell:',i,j,k
        endif
      enddo
    enddo
  enddo
#endif
  write(0,*) 'Num. of generated PBH(s): ',Ngen, ', ratio: ',float(Ngen)/Npbh
#ifdef MASSFUNC_LOGNORMAL
  ! 2-6: statistics of the sampled mass function
  if (Ngen .gt. 0) then
    mm_avg  = sum(real(mpbh_list(1:Ngen),8))/dble(Ngen)
    mm_std  = sqrt(sum((real(mpbh_list(1:Ngen),8)-mm_avg)**2)/dble(Ngen))
    lnm_avg = sum(log(real(mpbh_list(1:Ngen),8)))/dble(Ngen)
    lnm_std = sqrt(sum((log(real(mpbh_list(1:Ngen),8))-lnm_avg)**2)/dble(Ngen))
    write(0,*) 'PBH mass function: LOG-NORMAL, Mc=',Mc_pbh,' Msun, sigma=',sigma_ln
    write(0,*) '  mean   M  =',mm_avg,' Msun   (target <M> =',mpbh_mean,')'
    write(0,*) '  std    M  =',mm_std,' Msun   (target      =', &
      & mpbh_mean*sqrt(exp(sigma_ln**2)-1.0d0),')'
    write(0,*) '  min/max M =',minval(mpbh_list(1:Ngen)),maxval(mpbh_list(1:Ngen)),' Msun'
    write(0,*) '  mean ln M =',lnm_avg,'        (target ln Mc =',log(Mc_pbh),')'
    write(0,*) '  std  ln M =',lnm_std,'        (target sigma =',sigma_ln,')'
    write(0,*) '  total PBH mass =',sum(real(mpbh_list(1:Ngen),8)),' Msun'
  endif
#endif
    
  allocate(ldpbh(1:6,Ngen))
    
  if (Ngen .gt. 0) then
  
#ifndef ALLPBH
  if (Ngen .le. npert) then
    write(0,*) 'Only a few PBHs! You can use ALLPBH if rpert is large!'
    !stop  
  endif
#endif
#ifdef DEBUG
  write(0,*) 'PBH position range:', minval(lpbh(1:3,1:Ngen)), maxval(lpbh(1:3,1:Ngen))
  write(0,*) 'PDM velocity range:', minval(phaseDM(4:6,1:Ndm)), maxval(phaseDM(4:6,1:Ndm))
  write(0,*) 'PBH velocity range:', minval(lpbh(4:6,1:Ngen)), maxval(lpbh(4:6,1:Ngen))
#endif

#ifdef OUTPBH
  if (Ngen>0) then
  open(3, file=trim(path) // trim(icbase) //'_pbhlist.txt', status='unknown')
    do i=1,Ngen
      write(3,'(E14.6,E14.6,E14.6,E14.6,E14.6,E14.6)') lpbh(:,i)
    enddo
  close(3)
  endif
#endif

#ifdef PERTURB
! set the positions of the remaining elements of the PBH list to far away from
! the simulation box, so that they will not be included in perturbations
do i=Ngen+1,2*Ncell 
  lpbh(1:3,i) = (/10*bs,10*bs,10*bs/)
#ifdef MASSFUNC_LOGNORMAL
  mpbh_list(i) = 0.0   ! 2-3: placeholder mass, contributes nothing to any sum
#endif
enddo
! BUG FIX: idpbh is allocated with N (one entry per particle, only 1..Nhres is
! ever read or written), but this loop used to run to 2*Ncell, which is the
! size of the PBH list lpbh - clearly copy-pasted from the loop just above.
! 2*Ncell only depends on Npbh (i.e. on the total DM mass, fpbh and the mean
! PBH mass), NOT on the particle count, so whenever the PBH sampling produces
! more than N/2 candidate cells the write runs past the end of idpbh and
! -fbounds-check aborts.  With an extended mass function at ~19 Msun this is
! the normal case (2*Ncell = 182250 vs N = 65536).
! Conversely, when 2*Ncell < N the tail idpbh(2*Ncell+1:N) was left
! UNINITIALISED and then read below in "if (idpbh(i) .ne. j)"; it only
! happened to behave because fresh pages read back as zero and 0 is not a
! valid PBH index.
do i=1,N
  idpbh(i) = -1
enddo
write(0,*) 'Perturbing particles within', rpert
  dx_pert = 0.0
#ifdef VELOCITY
  dv_pert = 0.0
#endif
  npert=min(npert, Ngen)
  fac = 1.0
  !facv = 1.0
  facavg0 = 0.0
  facavg1 = 0.0
  ! perturb existing particles
  call kd_init(pbhroot, lpbh, bs)

#ifdef NONLINEAR  
  nonlinear = 1
  do i=1,Nhres
    if (Ngen .gt. 1) then
      call kd_nnearest(pbhroot, pos(1:3,i), 1, nnindex, nndis, i_dummy)
      if (nnindex(1) .gt. Ngen) then
        write(0,*) 'Fake particle from kdtree:', nnindex(1), pos(1,i), pos(2,i), pos(3,i)
        stop
      endif
      idpbh(i) = nnindex(1)
    else
      idpbh(i) = 1
    endif
    call dx_periodic(pos(1:3,i), lpbh(1:3,idpbh(i)), disp, bs)
    dummy = dist0(disp,(/0.,0.,0./),bs)
! 2-4: with an extended mass function the mass factor must sit INSIDE the sum
! over PBHs (each PBH carries its own mass); only the constant
! 4*PI*gra*msun/kpc/kpc stays outside.
#ifdef MASSFUNC_LOGNORMAL
    xdis = mpbh_list(idpbh(i)) * disp/dummy/(dummy**2 + (esft*dismax)**2)
    xdis = xdis * 4*PI*gra*msun/kpc/kpc !* fscale
#else
    xdis = disp/dummy/(dummy**2 + (esft*dismax)**2)
    xdis = xdis * 4*PI*gra*mpbh*msun/kpc/kpc !* fscale
#endif
    !xdis = xdis * 4/(Om0*H0**2)*log(0.5*(1.0+sqrt(1.0+a/aeq)))/kpc
    xdis = xdis * gfac/1.5/(Om0*H0**2)/kpc *h**3
    d_dummy = dist0(xdis,(/0.,0.,0./),bs)
    fac = min(1.0, fshrink*dummy)
    dx_pert = dx_pert + d_dummy*fac
    pos(1:3,i) = pos(1:3,i) - xdis*fac
    dv_pert = dv_pert + d_dummy * vfac *H0*kpc/h/1e5/sqrt(a) *fac
    vel(1:3,i) = vel(1:3,i) - xdis *vfac *H0*kpc/h/1e5/sqrt(a) *fac
  enddo
#else:
  nonlinear = 0
#endif
  
  do i=1,Nhres
    xdis = 0.0
    !rnearest = bs
    !dnearest = (/10*bs,10*bs,10*bs/)
  if (Ngen .gt. nonlinear) then
#ifdef ALLPBH
    do j=1,Ngen
      !call dx_periodic(pos(1:3,i), lpbh(1:3,j), disp, bs)
      !dummy = dist0(disp,(/0.,0.,0./),bs)
      !if (dummy .lt. rpert) then
      if (idpbh(i) .ne. j) then
        call dx_periodic(pos(1:3,i), lpbh(1:3,j), disp, bs)
        dummy = dist0(disp,(/0.,0.,0./),bs)
#ifdef MASSFUNC_LOGNORMAL
        xdis = xdis + mpbh_list(j) * disp/dummy/(dummy**2 + (esft*dismax)**2)
#else
        xdis = xdis + disp/dummy/(dummy**2 + (esft*dismax)**2)
#endif
        !if (dummy .lt. rnearest) then
        !  rnearest = dummy
        !  dnearest = disp
        !endif
      endif
    enddo
#else
    !if (i .eq. 1) then
    !  write(0,*) 'Perturbing particles within', rpert
    !endif
    call kd_nnearest(pbhroot, pos(1:3,i), npert, nnindex, nndis, i_dummy)
    do j=1,npert
      !if (nndis(j) .lt. rpert) then
      if (idpbh(i) .ne. nnindex(j)) then
        call dx_periodic(pos(1:3,i), lpbh(1:3,nnindex(j)), disp, bs)
        if (nnindex(j) .gt. Ngen) then
          write(0,*) 'Fake particle from kdtree:', nnindex(j), pos(1,i), pos(2,i), pos(3,i)
          stop
        endif
        dummy = dist0(disp,(/0.,0.,0./),bs)
#ifdef MASSFUNC_LOGNORMAL
        xdis = xdis + mpbh_list(nnindex(j)) * disp/dummy/(dummy**2 + (esft*dismax)**2)
#else
        xdis = xdis + disp/dummy/(dummy**2 + (esft*dismax)**2)  !dist0(pos(1:3,i), lpbh0(1:3,j))**3
#endif
        !if (dummy .lt. rnearest) then
        !  rnearest = dummy
        !  dnearest = disp
        !endif
      endif
    enddo
#endif
! 2-4: with an extended mass function the mass factor must sit INSIDE the sum
! over PBHs (each PBH carries its own mass); only the constant
! 4*PI*gra*msun/kpc/kpc stays outside.
#ifdef MASSFUNC_LOGNORMAL
    xdis = xdis * 4*PI*gra*msun/kpc/kpc !* fscale
#else
    xdis = xdis * 4*PI*gra*mpbh*msun/kpc/kpc !* fscale
#endif
    !xdis = xdis * 4/(Om0*H0**2)*log(0.5*(1.0+sqrt(1.0+a/aeq)))/kpc
    xdis = xdis * gfac/1.5/(Om0*H0**2)/kpc *h**3
    
    if (i .le. Ngas) then
      mref = massarr(0)
#ifdef DAMPING
      xdis = xdis * gfac_gas/gfac
#endif
    elseif (i .le. Nhres) then
      mref = massarr(1)
    else
      mref = lmass(i-Nhres)
    endif
    !if (mpbh .lt. mref*1e10/h) then
    !  fac = fac * mpbh*h/1e10/mref
    !endif
    d_dummy = dist0(xdis,(/0.,0.,0./),bs)
#ifdef TRUNCATE
    !dummy = p1dotp2(dnearest, xdis)
    !if (dummy .gt. 0) then
    !  dummy = min(ftrunc*dpertmax/d_dummy, fshrink*rnearest**2/dummy)
      !dummy = min(ftrunc*dpertmax/d_dummy, fshrink*rnearest/d_dummy)
    !else
      dummy = ftrunc*dpertmax/d_dummy
    !endif
    fac = min(dummy, 1.0)
    if (i .gt. Ngas) then
      facavg1 = facavg1+fac
   	else  
      facavg0 = facavg0+fac
    endif
#endif
    dx_pert = dx_pert + d_dummy*fac
    pos(1:3,i) = pos(1:3,i) - xdis*fac
#ifdef VELOCITY
    !dv_pert = dv_pert + d_dummy * H0*sqrt(Om0) * kpc/h/1e5/a *fac
    !vel(1:3,i) = vel(1:3,i) - xdis *H0*sqrt(Om0) * kpc/h/1e5/a *fac
    dv_pert = dv_pert + d_dummy * vfac *H0*kpc/h/1e5/sqrt(a) *fac
    vel(1:3,i) = vel(1:3,i) - xdis *vfac *H0*kpc/h/1e5/sqrt(a) *fac
    ! The internal velocity in GIZMO/GADGET is actually u=v*sqrt(a) 
    ! given the physical velocity v and scale factor a
#endif
  endif !if (Ngen .gt. nonlinear)
    do k=1,3
      if (pos(k,i) .gt. bs) then
        pos(k,i) = pos(k,i)-bs
      elseif (pos(k,i) .lt. 0) then
        pos(k,i) = pos(k,i)+bs
      endif
    enddo
  enddo
! NaN guard: each of these averages over a particle subset (DM, gas, or
! DM+gas) that can be empty (Ngas=0 is the common DM-only case), in which case
! the division is 0/0 = NaN. Skip the division and say so instead; no physics
! changes here, diagnostics only.
  if (Nhres .gt. 0) then
    write(0,*) 'Average displacement by PBHs:',dx_pert/float(Nhres)
  else
    write(0,*) 'Average displacement by PBHs: n/a (Nhres=0)'
  endif
  if (Ndm .gt. 0) then
    write(0,*) 'Average truncation factor (DM)     :',facavg1/Ndm
  else
    write(0,*) 'Average truncation factor (DM)     : n/a (Ndm=0)'
  endif
  if (Ngas .gt. 0) then
    write(0,*) 'Average truncation factor (gas)    :',facavg0/Ngas
  else
    write(0,*) 'Average truncation factor (gas)    : n/a (Ngas=0)'
  endif
  if (Nhres .gt. 0) then
    write(0,*) 'Average truncation factor (DM+gas) :',(facavg0+facavg1)/Nhres
  else
    write(0,*) 'Average truncation factor (DM+gas) : n/a (Nhres=0)'
  endif
#ifdef VELOCITY
  if (Nhres .gt. 0) then
    write(0,*) 'Average velocity induced by PBHs:',dv_pert/float(Nhres)
  else
    write(0,*) 'Average velocity induced by PBHs: n/a (Nhres=0)'
  endif
#endif

#ifdef PERTURB_PBH
  !if (Ngen .gt. 1) then
  do i=1,Ngen
    xdis = 0.0
    ldpbh(:,i) = 0.0
    do j=1,Ngen
      if (i .ne. j) then
        call dx_periodic(lpbh(1:3,i), lpbh(1:3,j), disp, bs)
        dummy = dist0(disp,(/0.,0.,0./),bs)
#ifdef MASSFUNC_LOGNORMAL
        xdis = xdis + mpbh_list(j) * disp/dummy/(dummy**2 + (esft*dismax)**2)
#else
        xdis = xdis + disp/dummy/(dummy**2 + (esft*dismax)**2)
#endif
      endif
    enddo
! 2-4: with an extended mass function the mass factor must sit INSIDE the sum
! over PBHs (each PBH carries its own mass); only the constant
! 4*PI*gra*msun/kpc/kpc stays outside.
#ifdef MASSFUNC_LOGNORMAL
    xdis = xdis * 4*PI*gra*msun/kpc/kpc !* fscale
#else
    xdis = xdis * 4*PI*gra*mpbh*msun/kpc/kpc !* fscale
#endif
    !xdis = xdis * 4/(Om0*H0**2)*log(0.5*(1.0+sqrt(1.0+a/aeq)))/kpc
    xdis = xdis * gfac/1.5/(Om0*H0**2)/kpc *h**3
    d_dummy = dist0(xdis,(/0.,0.,0./),bs)
#ifdef TRUNCATE
    !dummy = p1dotp2(dnearest, xdis)
    !if (dummy .gt. 0) then
    !  dummy = min(ftrunc*dpertmax/d_dummy, fshrink*rnearest**2/dummy)
      !dummy = min(ftrunc*dpertmax/d_dummy, fshrink*rnearest/d_dummy)
    !else
      dummy = ftrunc*dpertmax/d_dummy
    !endif
    fac = min(dummy, 1.0)
#endif
    ldpbh(1:3,i) = - xdis*fac
#ifdef VELOCITY
    !dv_pert = dv_pert + d_dummy * H0*sqrt(Om0) * kpc/h/1e5/a *fac
    !vel(1:3,i) = vel(1:3,i) - xdis *H0*sqrt(Om0) * kpc/h/1e5/a *fac
    ldpbh(4:6,i) = - xdis *vfac *H0*kpc/h/1e5/sqrt(a) *fac
    ! The internal velocity in GIZMO/GADGET is actually u=v*sqrt(a) 
    ! given the physical velocity v and scale factor a
#endif
  enddo
  !endif
  do i=1,Ngen
    lpbh(:,i) = lpbh(:,i) + ldpbh(:,i)
    do k=1,3
      if (lpbh(k,i) .gt. bs) then
        lpbh(k,i) = lpbh(k,i)-bs
      elseif (lpbh(k,i) .lt. 0) then
        lpbh(k,i) = lpbh(k,i)+bs
      endif
    enddo
  enddo
#endif

#endif ! #ifdef PERTURB
  endif ! if (Ngen .gt. 0)
#endif ! #ifdef PBH

  ! output data
if (Ngen .gt. 0) then
! --- mass-conservation fix ----------------------------------------------
! Previously: massarr(1) = massarr(1)*(1.0-fpbh). That used the TARGET input
! fraction fpbh instead of the REALISED PBH mass budget - Ngen is a Poisson
! draw around Npbh, and under MASSFUNC_LOGNORMAL each PBH's mass is itself a
! random draw around mpbh_mean, so sum(M_j) != fpbh*Ndm*m_dm in general. It
! also mixed real*4 fpbh into real*8 massarr(1) arithmetic ((1.0-fpbh) was
! evaluated at real*4 precision before being promoted).
!
! f_real mirrors the DM-mass convention used in the Npbh normalisation above:
! m_dm_msun gets the (Om0-Ob0)/Om0 factor only for a DM-only IC (Ngas=0),
! where a type-1 particle represents the combined baryon+DM budget and only
! the (Om0-Ob0)/Om0 slice of it is "DM" for the purpose of the fpbh target.
! f_real is therefore comparable to the input fpbh (diagnostic only).
!
! The mass actually DEBITED from massarr(1) must not use that factor, though:
! there is only one mass pool per type-1 particle (GIZMO cannot track a
! "DM-only" sub-mass within a combined particle), and what has to balance is
! the total mass handed to the PBHs, sum(M_j) = M_tot_msun, against the total
! mass removed from type-1, Ndm*(massarr_old-massarr_new)*1e10/h. Dividing
! M_tot_msun by (Om0-Ob0)/Om0 before subtracting it (i.e. applying f_real
! directly to massarr(1)) would remove MORE than M_tot_msun whenever Ngas=0
! (by a factor 1/((Om0-Ob0)/Om0) ~ 1.19 for Planck-like parameters) and silently
! break total mass conservation - caught by tools/check_snapshot.py's
! DM+PBH-vs-input-DM check. So the actual subtraction uses M_tot_msun directly.
  if (Ngas .gt. 0) then
    m_dm_msun = massarr(1)*1.0d10/h
  else
    m_dm_msun = massarr(1)*1.0d10/h * (dble(Om0)-dble(Ob0))/dble(Om0)
  endif
#ifdef MASSFUNC_LOGNORMAL
  M_tot_msun = sum(real(mpbh_list(1:Ngen),8))
#else
  M_tot_msun = dble(Ngen)*mpbh
#endif
  f_real = M_tot_msun / (dble(Ndm)*m_dm_msun)
  mass_removed_msun = M_tot_msun
  massarr(1) = massarr(1) - mass_removed_msun*h/(dble(Ndm)*1.0d10)
  write(0,*) 'PBH mass budget: f_pbh(input)=',fpbh,' f_real(realised, vs true DM mass)=',f_real
  write(0,*) '  sum M_j (PBHs)   [Msun] =', M_tot_msun
  write(0,*) '  DM mass removed  [Msun] =', mass_removed_msun
  write(0,*) '  relative mismatch (removed-added)/added =', &
    & (mass_removed_msun-M_tot_msun)/M_tot_msun
! -------------------------------------------------------------------------
#ifdef MASSFUNC_LOGNORMAL
  ! 2-5: with an extended mass function every PBH has its own mass, so the
  ! header slot must be 0 and the masses move into the variable-mass block.
  ! Ordering: the variable-mass block follows particle-type order, and the
  ! PBHs were inserted into pos1/vel1/id1 between the type 0/1 block and the
  ! type >= 2 block, so their masses go at offset Nm_hres.
  ! Units: mpbh_list is in physical Msun, the snapshot uses 1e10 Msun/h,
  ! hence the same h/1e10 conversion that massarr(pbh_type) used to apply.
  massarr(pbh_type) = 0.0
  allocate(lmass1(1:Nm+Ngen))
  if (Nm_hres .gt. 0) lmass1(1:Nm_hres) = lmass(1:Nm_hres)
  lmass1(Nm_hres+1:Nm_hres+Ngen) = mpbh_list(1:Ngen) * real(h/1e10)
  if (Nm .gt. Nm_hres) lmass1(Nm_hres+Ngen+1:Nm+Ngen) = lmass(Nm_hres+1:Nm)
  Nm = Nm + Ngen
  write(0,*) 'Variable-mass block: Nm=',Nm,' (',Nm_hres,' from types 0/1 +', &
    & Ngen,' PBHs + ',Nm-Nm_hres-Ngen,' from types >=2), massarr(pbh_type)=0'
#else
  massarr(pbh_type) = mpbh*h/1e10
#endif
  nall(pbh_type) = Ngen
  npart(pbh_type) = Ngen
  allocate(pos1(1:3,1:N+Ngen))
  pos1(1:3,1:Nhres) = pos(1:3,1:Nhres)
  pos1(1:3,Nhres+1:Nhres+Ngen) = lpbh(1:3,1:Ngen)
  pos1(1:3,Nhres+Ngen+1:N+Ngen) = pos(1:3,Nhres+1:N)
  deallocate(pos)
  allocate(vel1(1:3,1:N+Ngen))
  vel1(1:3,1:Nhres) = vel(1:3,1:Nhres)
  vel1(1:3,Nhres+1:Nhres+Ngen) = lpbh(4:6,1:Ngen)
  vel1(1:3,Nhres+Ngen+1:N+Ngen) = vel(1:3,Nhres+1:N)
  deallocate(vel)
  allocate(id1(1:N+Ngen))
  do i=1,Ngen
    id1(Nhres+i) = N+i-1
  enddo
  id1(Nhres+Ngen+1:N+Ngen) = id(Nhres+1:N)
else
  allocate(pos1(1:3,1:N))
  pos1=pos
  deallocate(pos)
  allocate(vel1(1:3,1:N))
  vel1=vel
  deallocate(vel)
  allocate(id1(1:N))
  id1(Nhres+1:N) = id(Nhres+1:N)
  !deallocate(id)
#ifdef PBH
  !write(0,*) 'Error: no PBH!'
  write(0,*) 'No PBH! Please modify PBH parameters or the random seed!'
#endif
endif

#ifdef PBH
! --- PBH catalogue -------------------------------------------------------
! Written unconditionally whenever PBH is defined (tiny file). Records
! lpbh/mpbh_list AFTER PERTURB_PBH (if enabled) and the id1 value actually
! assigned above - i.e. exactly what pos1/vel1/id1 carry into the snapshot
! for the type-pbh_type particles, not the pre-perturbation sample.
! This is independent of OUTPBH, which (if enabled) instead dumps the
! UNPERTURBED lpbh (no ID, no mass) to <icbase>_pbhlist.txt right after PBH
! generation, before PERTURB_PBH; the two are not redundant, but this file
! is the complete, canonical record of what was written to disk.
  open(4, file=trim(path)//trim(icbase)//'_pbh_list.txt', status='unknown')
  write(4,'(A)') '# j ID x y z[kpc/h] vx vy vz[u=v*sqrt(a)] M[Msun]'
  if (Ngen .gt. 0) then
    do i=1,Ngen
#ifdef MASSFUNC_LOGNORMAL
      write(4,'(I10,1X,I12,1X,7ES16.8)') i, id1(Nhres+i), real(lpbh(1:6,i),8), &
        & real(mpbh_list(i),8)
#else
      write(4,'(I10,1X,I12,1X,7ES16.8)') i, id1(Nhres+i), real(lpbh(1:6,i),8), mpbh
#endif
    enddo
#ifdef MASSFUNC_LOGNORMAL
    write(4,'(A,ES16.8)') '# total M [Msun] = ', sum(real(mpbh_list(1:Ngen),8))
#else
    write(4,'(A,ES16.8)') '# total M [Msun] = ', dble(Ngen)*mpbh
#endif
  else
    write(4,'(A,ES16.8)') '# total M [Msun] = ', 0.0d0
  endif
  close(4)
! -------------------------------------------------------------------------
#endif

#ifdef FLIPID
  ! flip the IDs of gas and dark matter particles (may be needed for arepo)
  id1(1:Ngas) = id(1:Ngas)+Ndm
  id1(Ngas+1:Nhres) = id(Ngas+1:Nhres)-Ndm
#else
  id1(1:Ngas) = id(1:Ngas)
  id1(Ngas+1:Nhres) = id(Ngas+1:Nhres)
#endif
  deallocate(id)
  !write(0,*) id1(1), id1(Ngas+1)
  
  do i=0,5
     print *,'Type=',i,' Particles=', nall(i), 'mass=', massarr(i)
  end do
  !write(0,*) 'check PBH mass:',massarr(pbh_type)*1e10/h
#ifdef PBH
  outname = trim(outname) // '_pbh'
#ifndef PERTURB
  outname = trim(outname) // '_uptb'
#endif
#else
  outname = trim(outname) // '_cdm'
#endif
#ifdef STREAMING
  outname = trim(outname) // '_streaming'
#endif
  open(2, file=trim(outname) // '.dat', form='unformatted')
    write(2) npart, massarr ,a, redshift, flag_sfr,flag_feedback, nall, flag_cool, nfs, boxsize, Om0, Ol, h, unused
    write(2) pos1
    write(2) vel1
    write(2) id1
    if (Nm .gt. 0) then
#ifdef MASSFUNC_LOGNORMAL
      if (Ngen .gt. 0) then
        write(2) lmass1     ! 2-5: type 0/1 masses + PBH masses + type >=2 masses
      else
        write(2) lmass
      endif
#else
      write(2) lmass
#endif
    endif
    if (Ngas .gt. 0) then
      write(2) iegas
    endif
  close(2)

contains

#ifdef MASSFUNC_LOGNORMAL
! 2-3: one standard normal deviate via the Box-Muller transform, built from
! two uniforms of rng_obj%get_ran().  get_ran() returns [0,1), so u1 can come
! out as exactly 0; redraw in that case to keep log(u1) finite.
real function gauss_bm(rng_obj)
    use random_object
    implicit none
    type(rng) :: rng_obj
    real :: u1, u2
    u1 = rng_obj%get_ran()
    do while (u1 .le. 0.0)
      u1 = rng_obj%get_ran()
    end do
    u2 = rng_obj%get_ran()
    gauss_bm = sqrt(-2.0*log(u1)) * cos(2.0*PI*u2)
end function gauss_bm
#endif

integer function poisson(mu, rng_obj)
    use random_object
    implicit none
    real, intent(in) :: mu
    type(rng) :: rng_obj
    real L, p

    L = exp(-mu)
    p = 1.0
    poisson = 0
    do 
      p = p*rng_obj%get_ran()
      if (p .le. L) exit
      poisson = poisson + 1
    end do
end function poisson

  pure function dist0(p1, p2, bs)
    real, intent(in) :: bs, p1(:), p2(:)
    real :: dist0
    integer :: n
    real :: r1, r2
    r1 = 0.0
    do n=1, 3 !size(p1)
      r2 = p1(n)-p2(n)
!#ifdef PERIOD
!      if (r2 .gt. bs*0.5) then
!        r2 = r2-bs
!      endif
!      if (r2 .lt. -0.5*bs) then
!        r2 = r2+bs
!      endif
!#endif
      r1 = r1 + r2*r2
    enddo
    dist0 = sqrt(r1)
  end function

  subroutine dx_periodic(p1, p2, dp, bs)
    real, intent(in) :: bs, p1(:), p2(:)
    real, intent(out) :: dp(:)
    integer :: n
    real :: r2
    do n=1, 3
#ifdef PERIOD
      r2 = p1(n) - p2(n)
      if (r2 .gt. bs*0.5) then
        dp(n) = r2-bs
      elseif (r2 .lt. -0.5*bs) then
        dp(n) = r2+bs
      else
        dp(n) = r2
      endif
#else
      dp(n) = p1(n) - p2(n)
#endif
    enddo
  end subroutine

  pure function Dgrow(s, fc, a_)
    real, intent(in) :: fc, a_
    real*8, intent(in) :: s
    real :: Dgrow
#ifdef SELFS !try to comment out this part, reduce to the same initial condition as we did 2yrs ago
    if (s .gt. 1) then
      !Dgrow = s*6.0*log(0.5+0.5**0.5)
      Dgrow = s*3.0/(4*PI)/1.3
    else
      Dgrow = 6*log(0.5+0.5*(1+s)**0.5) !(1.0+1.5*fc/a_*s)**a_ - 1.0
      !Dgrow = s*3.0/(4*PI)
    endif
#else
    Dgrow = (1.0+1.5*fc/a_*s)**a_ - 1.0
#endif
  end function

  pure function Dgrowv(a, aeq, fc, a_)
    real, intent(in) :: fc, a_, aeq
    real*8, intent(in) :: a
    real :: Dgrowv, hubble
!#ifdef SELFS
    hubble = sqrt(Om0*(1.0/a**3+aeq/a**4)+Ol)
    !Dgrowv = a**2*hubble*1.5/aeq/(1.5*a/aeq+1.0)
    Dgrowv = a*hubble
!#else
!    Dgrowv = sqrt(Om0/a)
!#endif
  end function

  pure function p1dotp2(p1, p2)
    real, intent(in) :: p1(:), p2(:)
    real :: p1dotp2
    integer :: n
    p1dotp2 = 0.0
    do n=1,3
      p1dotp2 = p1dotp2 + p1(n)*p2(n)
    enddo
  end function

end program




