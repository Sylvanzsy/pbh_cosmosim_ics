! Add PBHs to a Gadget format 0 IC file
#define DEBUG
!#define OUTPBH ! write an output file for the initial position/velocity of PBHs
!#define PERTURB
#define STREAMING ! add a baryon-DM relative streaming velocity
!#define FLIPID
#define PBH
!#define DAMPING
!#define VELOCITY
!#define TRUNCATE
#define CENTER_PBH ! place one single PBH in the center of the box
!#define ALLPBH
#define PERIOD
!#define NONLINEAR ! initial nonlinear structures around PBHs
!#define PERTURB_PBH ! perturbations on PBHs by themselves
#define NO_GAS_PERTURB ! no perturbations to gas particles
!#define noPBH ! handle to not include PBH generated but only include perturbation induced by PBHs in IC

program pbh_ic

   use random_object
   use kdtree

   implicit none


   character*200, parameter :: path='../ics_SIDM_cdm/' 
   !character*200, parameter :: icbase='ics_N512L1'
   !character*200, parameter :: icbase='ics_N256L1'
   !character*200, parameter :: icbase='ics_N128L014_nu2_zoom_lowres'
   !character*200, parameter :: icbase='ics_N128L014_zoom_lowres'
   !character*200, parameter :: icbase='ics_N512L060'
   !character*200, parameter :: icbase='ics_N512L060_z1000'
   !character*200, parameter :: icbase='ics_N1024L060_z1000'
   !character*200, parameter :: icbase='ics_N256L025'
   character*200, parameter :: icbase='ics_N256L025_z1100'
   !character*200, parameter :: icbase='ics_N256L1_z1100'
   !character*200, parameter :: icbase='ics_N256L1_z1100_sstr'
   !character*200, parameter :: icbase='ics_N256L1_z1100_str'
   !character*200, parameter :: icbase='ics_N256L1_z3400'
   !character*200, parameter :: icbase='ics_N1024L060_dmo'
   !character*200, parameter :: icbase='ics_N128L2'
   !character*200, parameter :: icbase='ics_N256L10'
   !character*200, parameter :: icbase='ics_N32L025_nu3_dmo_z1e6'
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
   real*4 :: p_dummy(6), dx(3), fac, dx_pert, rpert, fac_pert, disp(3), bs, esft, rpert0
#ifdef VELOCITY
   real*4 :: dv_pert
#endif
   type(kd_root) :: kdroot, pbhroot
   integer*4 :: pbh_type, nx, nxdm, nnn, npert, Ncell, Ngen, Nthis, xind(3) !, Ngen0
   integer*4, allocatable :: nnindex(:) !, n_dummy
   real*4, allocatable :: nndis(:), dmgrid(:,:,:)
   real*4 :: dismax, vol, re(3), le(3), rvec(3), xdis(3), ftrunc, H0, mref, dpertmax
   real*4 :: dummy, facavg0, facavg1, fshrink, d_dummy, vdis(3)
   real*4 :: delta, fscale, fscalev!, rnearest, dnearest(3)
   type(rng) :: rng_obj
   
   real*4 :: msun, kpc, gra, PI
   parameter(msun=1.989e33, kpc=3.085678e21, gra=6.672e-8, PI=3.1415926536)
   
   real*8 :: mpbh, sigma_rms, z_rec, fstream, fpbh_, d_nl!, a_target
   real*4 :: fpbh, Npbh, Or0, aeq, vrel(3), fc, Ob0, a_, gfac, gfac_gas, vfac!, facv
   parameter(sigma_rms=30.0, z_rec=1100.0, fstream=0.4)!, a_target=1.0/31.0)
   !parameter(nnn=1, pbh_type=5, mpbh=5e10, fpbh=5e-4)
   !parameter(nnn=1, pbh_type=3, mpbh=1e6, fpbh=1e-3)
   parameter(nnn=3000, pbh_type=3, mpbh=1e6, fpbh=1.0)
   !parameter(nnn=1, pbh_type=3, mpbh=1e10, fpbh=1e-3)
   logical :: zoomin
   parameter(zoomin=.false., rand_seed=2333, n_refine=1) !use rand_seed=2333 for 10 PBH with fpbh=1e-3, mpbh=1e10 in a 10 Mpc/h box
   parameter(esft=1.0/30.0, Or0=9.117e-5, Ob0=0.04864)
   parameter(ftrunc=1.0, fac_pert=100.0, rpert0=0)!300)
   parameter(delta=18.0*PI**2, fscalev=1.0)
   !parameter(delta=18.0*PI**2, fscale=0.16260587661983628, fscalev=1.0)
   !parameter(delta=18.0*PI**2, fscale=1, fscalev=0.16260587661983628)!0.21138763960578716)
   
   npert=64 !1
   
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
  
  if (Ngas .gt. 0) then
    Npbh = float(Ndm) * massarr(1)*1e10*fpbh/h/mpbh
  else
    Npbh = float(Ndm) * massarr(1)*1e10*fpbh/h/mpbh * (Om0-Ob0)/Om0
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
  !fscale = 1.0 !2.0/Om0
  fscale = 4.4
  
  gfac = Dgrow(a/aeq, fc, a_) * fscale
#ifdef NO_GAS_PERTURB
  gfac_gas = 0.0
#else  
  gfac_gas = Dgrow(a*(1+z_rec), fc, a_) * fscale
#endif
  vfac = 1.0 !Dgrowv(a, aeq, fc, a_) * fscalev
  
!#ifdef IHALO
  dummy = (8*PI/3.0*gra*mpbh*msun * gfac/(Om0*H0**2))**(1.0/3.0) * h/kpc
  d_nl = q_nl(mpbh, a, fscale)
!  dpertmax = dummy !*(4*PI/3.0)*(Dgrow(a/aeq, fc, a_)/Dgrow(a_target/aeq, fc, a_))**(2.0/3.0)
  write(0,*) 'Nonlinear scale = ', dummy, dummy/dismax, ' (old)', d_nl, d_nl/dismax, ' (new)'!, dummy
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
  do i=1,3
    lpbh(i,1) = 0.5*(re(i)+le(i))
  enddo


  ! Add random perturbation to the PBH position
  !call rng_obj%init(rand_seed)
  !do i=1,3
    ! The perturbation magnitude should be smaller than the average particle separation (dismax)
  !  lpbh(i,1) = lpbh(i,1) + (2.0*rng_obj%get_ran() - 1.0) * 0.2 * dismax
  !enddo
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
  fpbh_ = fpbh * float(Ngen)/Npbh
    
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

! set the positions of the remaining elements of the PBH list to far away from
! the simulation box, so that they will not be included in perturbations
#ifdef PERTURB
  do i=Ngen+1,2*Ncell 
    lpbh(1:3,i) = (/10*bs,10*bs,10*bs/)
  enddo
  do i=1,2*Ncell
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
    if (dummy .gt. d_nl) then
      xdis = disp/dummy *psi_l(dummy**2 + (esft*dismax)**2, a, mpbh, fscale)
      vdis = disp/dummy *psi_dot_l(dummy**2 + (esft*dismax)**2, a, mpbh, fscale)/sqrt(a)
    else
      xdis = disp/dummy *psi_nl(dummy**2 + (esft*dismax)**2, a, mpbh, fscale)
      vdis = disp/dummy *psi_dot_nl(dummy**2 + (esft*dismax)**2, a, mpbh, fscale)/sqrt(a)
    endif
    
    if (i .le. Ngas) then
#if defined(DAMPING) || defined(NO_GAS_PERTURB)
      xdis = xdis * gfac_gas/gfac
#endif
    endif    

    d_dummy = dist0(xdis,(/0.,0.,0./),bs)
    fac = d_dummy/psi_l(dummy**2 + (esft*dismax)**2, a, mpbh, fscale)
    
    if (i .gt. Ngas) then
      facavg1 = facavg1+fac
   	else  
      facavg0 = facavg0+fac
    endif
    
    dx_pert = dx_pert + d_dummy !*fac
    pos(1:3,i) = pos(1:3,i) - xdis !*fac
    dv_pert = dv_pert + dist0(vdis,(/0.,0.,0./),1e20)
    vel(1:3,i) = vel(1:3,i) - vdis
  enddo
#else:
  nonlinear = 0
#endif
  
  do i=1,Nhres
    xdis = 0.0
    vdis = 0.0
!#ifdef NO_GAS_PERTURB    
!    if (i .le. Ngas) then
!      cycle
!    endif
!#endif
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
        xdis = xdis + disp/dummy *psi_l(dummy**2 + (esft*dismax)**2, a, mpbh, fscale)
        vdis = vdis + disp/dummy *psi_dot_l(dummy**2 + (esft*dismax)**2, a, mpbh, fscale)/sqrt(a)
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
      if (nndis(j) .lt. rpert) then
      if (idpbh(i) .ne. nnindex(j)) then
        call dx_periodic(pos(1:3,i), lpbh(1:3,nnindex(j)), disp, bs)
        if (nnindex(j) .gt. Ngen) then
          write(0,*) 'Fake particle from kdtree:', nnindex(j), pos(1,i), pos(2,i), pos(3,i)
          stop
        endif
        dummy = dist0(disp,(/0.,0.,0./),bs)
        xdis = xdis + disp/dummy *psi_l(dummy**2 + (esft*dismax)**2, a, mpbh, fscale)  !dist0(pos(1:3,i), lpbh0(1:3,j))**3
        vdis = vdis + disp/dummy *psi_dot_l(dummy**2 + (esft*dismax)**2, a, mpbh, fscale)/sqrt(a)
        !if (dummy .lt. rnearest) then
        !  rnearest = dummy
        !  dnearest = disp
        !endif
      endif
      endif
    enddo
#endif
    !xdis = xdis * 4*PI*gra*mpbh*msun/kpc/kpc !* fscale
    !xdis = xdis * 4/(Om0*H0**2)*log(0.5*(1.0+sqrt(1.0+a/aeq)))/kpc
    !xdis = xdis * gfac/1.5/(Om0*H0**2)/kpc *h**3
    
    if (i .le. Ngas) then
      mref = massarr(0)
#if defined(DAMPING) || defined(NO_GAS_PERTURB)
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
    dv_pert = dv_pert + dist0(vdis,(/0.,0.,0./),1e20)*fac
    vel(1:3,i) = vel(1:3,i) - vdis*fac
    ! The internal velocity in GIZMO/GADGET is actually u=v/sqrt(a) 
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
  write(0,*) 'Average displacement by PBHs:',dx_pert/float(Nhres)
  write(0,*) 'Average truncation factor:',facavg1/Ndm,facavg0/Ngas,(facavg0+facavg1)/Nhres
#ifdef VELOCITY
  write(0,*) 'Average velocity induced by PBHs:',dv_pert/float(Nhres)
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
        xdis = xdis + disp/dummy *psi_l(dummy**2 + (esft*dismax)**2, a, mpbh, fscale)
        vdis = vdis + disp/dummy *psi_dot_l(dummy**2 + (esft*dismax)**2, a, mpbh, fscale)/sqrt(a)
      endif
    enddo
    !xdis = xdis * 4*PI*gra*mpbh*msun/kpc/kpc !* fscale
    !xdis = xdis * 4/(Om0*H0**2)*log(0.5*(1.0+sqrt(1.0+a/aeq)))/kpc
    !xdis = xdis * gfac/1.5/(Om0*H0**2)/kpc *h**3
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
    ldpbh(4:6,i) = - vdis*fac
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

#ifdef noPBH
  if (Ngen .gt. 0) then
    Ngen = 0
    write(0,*) 'Remove generated PBH but keep perturbation'
  endif
#endif

  ! output data
if (Ngen .gt. 0) then
  massarr(1) = massarr(1)*(1.0-fpbh_)
  massarr(pbh_type) = mpbh*h/1e10
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
endif ! if (Ngen .gt. 0)

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
#ifdef CENTER_PBH
  outname = trim(outname) // '_pbhcen'
#endif
#ifndef PERTURB
  outname = trim(outname) // '_uptb'
#endif
#ifdef NO_GAS_PERTURB
  outname = trim(outname) // '_ngasptb'
#endif
! #ifndef SELFS
! outname = trim(outname) // '_nself'
#endif
#ifndef NONLINEAR
  outname = trim(outname) // '_nlinptb'
#endif
#ifdef noPBH
  outname = trim(outname) // '_pertonly'
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
      write(2) lmass
    endif
    if (Ngas .gt. 0) then
      write(2) iegas
    endif
  close(2)

contains

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

! Jiao Hao's model for point-source perturbations

! linear comoving displacement [kpc/h]
! dis2: comving distance squared [kpc/h]
! a: scale factor
! mpbh: point-source mass [Msun]
! fscale: normalization factor
  pure function psi_l(dis2, a, mpbh, fscale)
    real*8, intent(in) :: a, mpbh 
    real*4, intent(in) :: dis2, fscale
    real*4 :: psi_l, t0
    t0 = 1/H0
    psi_l = 1.5 * gra*mpbh*msun * t0**2 / dis2 /(kpc/h)**2
    psi_l = psi_l * (0.4*(aeq/a)**1.5 + 0.6*a/aeq - 1) /(kpc/h) * fscale
  end function

! linear physical velocity change [km/s]
  pure function psi_dot_l(dis2, a, mpbh, fscale)
    real*8, intent(in) :: a, mpbh 
    real*4, intent(in) :: dis2, fscale
    real*4 :: psi_dot_l, t0
    t0 = 1/H0
    psi_dot_l = 0.6 * gra*mpbh*msun *t0/aeq**1.5 / dis2 /(kpc/h)**2
    psi_dot_l = psi_dot_l * ((aeq/a)**(0.5) -(aeq/a)**3) *a/1e5 * fscale!/ sqrt(a)
  end function

! non-linear scale: critical comoving distance [kpc/h] for the turnaround shell
  pure function q_nl(mpbh, a, fscale)
    real*8, intent(in) :: mpbh, a
    real*4, intent(in) :: fscale
    real*4 :: q_nl, t0
    t0 = 1/H0
    q_nl = (1.8 * gra*mpbh*msun *t0**2 *a/aeq)**(1.0/3.0) /(kpc/h) !* fscale
  end function

! comoving displacement [kpc/h] in the nonlinear (virialized) region
  pure function psi_nl(dis2, a, mpbh, fscale)
    real*8, intent(in) :: a, mpbh 
    real*4, intent(in) :: dis2, fscale
    real*4 :: psi_nl, t0, ata, q
    q = sqrt(dis2)/(kpc/h)**2
    t0 = 1/H0
    ata = aeq * 5 * q**3/(9.0 * gra*mpbh*msun * t0**2)
    psi_nl = q * (1.0 - ata/a/4.0) /(kpc/h) !* fscale
  end function

! physical velocity change [km/s] in the nonlinear (virialized) region
  pure function psi_dot_nl(dis2, a, mpbh, fscale)
    real*8, intent(in) :: a, mpbh 
    real*4, intent(in) :: dis2, fscale
    real*4 :: psi_dot_nl, psi, t, q
    q = sqrt(dis2)/(kpc/h)**2
    t = a**1.5/H0
    psi = q - psi_nl(dis2, a, mpbh, fscale)
    psi_dot_nl = 2.0/3.0 * psi / t
  end function

  pure function Dgrow(s, fc, a_)
    real, intent(in) :: fc, a_
    real*8, intent(in) :: s
    real :: Dgrow
!#ifdef SELFS ! No longer in use since we have used Jiao's formulism
!    if (s .gt. 1) then
!      !Dgrow = s*6.0*log(0.5+0.5**0.5)
!      Dgrow = s*3.0/(4*PI)/1.3
!    else
!      Dgrow = 6*log(0.5+0.5*(1+s)**0.5) !(1.0+1.5*fc/a_*s)**a_ - 1.0
      !Dgrow = s*3.0/(4*PI)
!    endif
!#else
    Dgrow = (1.0+1.5*fc/a_*s)**a_ - 1.0
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




