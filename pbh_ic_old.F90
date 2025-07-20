! Add PBHs to a Gadget format 0 IC file
#define DEBUG
!#define OUTPBH
#define PERTURB
!#define STREAMING
!#define FLIPID
#define PBH
#define DAMPING
#define VELOCITY
#define TRUNCATE
#define ALLPBH
#define PERIOD

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

   character*200 filename, outname

   integer*4 npart(0:5), nall(0:5)
   real*8    massarr(0:5)
   real*8    a,redshift,boxsize
   real*8    Om0,Ol,h
   integer*4 unused(24)
   integer*4 fn,i,j,k,l,ipbh,nstart,flag_sfr,flag_feedback,flag_cool,nfs
   integer*4 Ngas,Ndm,N,Nm,Nhres,i_dummy,rand_seed
! a, redshift, boxsize, Om0, Ol, h: These store the cosmological parameters.
! unused: unused columns from input files


   real*4,allocatable    :: pos(:,:),vel(:,:),lmass(:),iegas(:)
   integer*4,allocatable :: id(:)
   real*4,allocatable    :: pos1(:,:),vel1(:,:),lmass1(:)
   integer*4,allocatable :: id1(:)
! pos and vel: 2D arrays holding the positions and velocities of all particles in the simulation, respectively.

   real*4,allocatable    :: phaseDM(:,:), lpbh(:,:) !, lpbh0(:,:)
! lpbh: 2D array that will hold the positions and velocities of the PBHs.
! phaseDM: 2D array that holds the phase space information (positions and velocities) of the dark matter particles.


   real*4 :: p_dummy(6), dx(3), fac, dx_pert, rpert, fac_pert, disp(3), bs, esft
#ifdef VELOCITY
   real*4 :: dv_pert
#endif
   type(kd_root) :: kdroot, pbhroot
! pbhroot and kdroot: Root nodes of two KD-trees for efficient nearest-neighbor searches, one for the dark matter particles and the other for the PBHs.

   integer*4 :: pbh_type, nx, nxdm, nnn, npert, Ncell, Ngen, Nthis, xind(3) !, Ngen0
! npert: sets the number of nearest neighbors used in a KDTree nearest neighbor search.
! Ngen: Number of PBHs actually generated.
! Ncell: Total number of cells in the simulation box.
! pbh_type: type index for PBHs.0: gas; 1: Halo (dark matter); 2 to 5: different types of matter.

   integer*4, allocatable :: nnindex(:) !, n_dummy
   real*4, allocatable :: nndis(:), dmgrid(:,:,:)
   real*4 :: dismax, vol, re(3), le(3), rvec(3), xdis(3), ftrunc, H0, mref
   real*4 :: dummy, facavg0, facavg1
! bs, dismax, rvec: Variables related to the size and shape of the simulation volume and the individual cells within it.
! le(3), re(3), vol: These are the left edge, right edge, and volume of the region of interest.

   type(rng) :: rng_obj
! rng_obj: An object of a random number generator class, used to generate random numbers
   
   real*4 :: msun, kpc, gra, PI
   parameter(msun=1.989e33, kpc=3.085678e21, gra=6.672e-8, PI=3.1415926536)
   
   real*8 :: mpbh, sigma_rms, z_rec, fstream
   real*4 :: fpbh, Npbh, Or0, aeq, vrel(3), fc, Ob0, a_
!Expected number of PBHs in the simulation
!mpbh and fpbh: Mass of PBHs and fraction of dark matter in PBHs
!vrel: Relative velocity of dark matter and gas


   parameter(sigma_rms=30.0, z_rec=1100.0, fstream=0.8)
   !parameter(nnn=1, pbh_type=5, mpbh=1e10, fpbh=1e-4)
   !parameter(nnn=1, pbh_type=3, mpbh=1e6, fpbh=1e-3)
   !parameter(nnn=1, pbh_type=3, mpbh=1e6, fpbh=1e-2)
   parameter(nnn=1, pbh_type=3, mpbh=33, fpbh=1e-3)
!parameter: indicate that these are constants and cannot be changed during the execution of the program.
   logical :: zoomin
! boolean type variable, as a switch for zoom-in mode

   parameter(zoomin=.false., rand_seed=23333) !use rand_seed=2333 for 10 PBH with fpbh=1e-3, mpbh=1e10 in a 10 Mpc/h box
   parameter(ftrunc=1.0, esft=1.0/30.0, Or0=9.117e-5, fac_pert=2.0, Ob0=0.04864)
   
   npert=64
   
   allocate(nnindex(max(nnn,npert)))
   allocate(nndis(max(nnn,npert)))

   filename= path(1:len_trim(path)) // icbase(1:len_trim(icbase)) // '_cdm.dat'
   outname= path(1:len_trim(path)) // icbase(1:len_trim(icbase))
   print *,'opening...  '//filename(1:len_trim(filename))
   ! read in the header
   Nm = 0
   open (1, file=filename, form='unformatted')
   read (1) npart, massarr ,a, redshift, flag_sfr,flag_feedback, nall, flag_cool, nfs, boxsize, Om0, Ol, h, unused
!massarr: An array that holds the masses of the particles in the simulation, split by type.
!    massarr(0): gas particles, massarr(1): dark matter particle
!Om0, Ol, and h: Cosmological parameters - matter density (Omega_m), dark energy density (Omega_lambda), and Hubble constant (h),
!npart and nall: Arrays that contain the number of particles of each type in the file and in the total simulation, respectively.


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
! Ngas, Ndm, and N: Number of gas particles, dark matter particles, and total particles, respectively.
     ! read in data
     !write(0,*) Ngas, N


      allocate(pos(1:3,1:N))
      allocate(vel(1:3,1:N))
      allocate(id(1:N)) ! allocate statement dynamically set the array sizes at runtime
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
  
!zoomin:  .true.: it calculates the volume (vol) of a parallelepiped (likely in 3D space) using the data from the array phaseDM. Otherwise, it assigns simple values based on the variable boxsize.
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
  
! STREAMING: add a velocity to the gas particles.
#ifdef STREAMING
  if (Ngas .gt. 0) then
	call rng_obj%init(rand_seed*2) !random number generator object rng_obj is initialized with the seed rand_seed * 2.
    rvec(1) = rng_obj%get_ran()!A random 3D vector rvec is then generated using the rng_obj%get_ran() function
    rvec(2) = rng_obj%get_ran()
    rvec(3) = rng_obj%get_ran()
    vrel = rvec/sum(rvec**2)**0.5*sigma_rms*(1.0+redshift)/(1.0+z_rec) *fstream/a**0.5
! randomized stream velocity vector with dispersion of 'sigma_rms', and scale as (1+z)
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
! calculates the expected number of PBHs based on the mass fraction fpbh and generates a spatial distribution of these PBHs.
! constructs a kd-tree of DM particles.
! generates the PBHs and maps their velocities to the local DM velocities.
  ! generate a grid in which each cell contains ~ 1 PBH on average
  nxdm = int(float(Ndm)**(1.0/3.0)+0.5) ! computes an approximate number of grid cells along one dimension for DM particles
  dismax = (vol/float(nall(1)))**(1.0/3.0)! computes the maximum distance between DM particles

! if in "zoom-in" region, le and re is adjusted by half the maximum DM distance
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
  write(0,*) 'Maximum distance = ', dismax
  write(0,*) 'Target region extent: le=', le,'re=', re
  
! calculate # of PBHs
! change this part to apply for extended mass case
  if (Ngas .gt. 0) then
    Npbh = float(Ndm) * massarr(1)*1e10*fpbh/h/mpbh ! massarr: DM particle mass, units?
  else
    Npbh = float(Ndm) * massarr(1)*1e10*fpbh/h/mpbh * (Om0-Ob0)/Om0 ! DM only case
  endif


!nx:  an approximate number of grid cells along one dimension for PBH particles based on Npbh.
!Ncell: the total number of cells for PBHs in 3D.
!rpert:  some perturbation radius for PBHs based on the volume and number of PBHs.
  nx = int(Npbh**(1.0/3.0)+1.0)
  Ncell = float(nx)**3
  rpert = fac_pert * (vol/Npbh)**(1.0/3.0)
  write(0,*) '1D grid num.:',nx
  write(0,*) 'Num. of cells for PBH sampling:',nx**3
  write(0,*) 'PBH num. per cell:', Npbh/Ncell
  write(0,*) 'Num. of PDM particles per cell:', float(Ndm)/float(Ncell)
  
  ! mapping DM particles to a grid
  allocate(dmgrid(nx, nx, nx))
  dmgrid = 0.0 ! allocate and initialize grid
  dx = (re-le)/float(nx) !spacing
  do i=1,Ndm
    do j=1,3
      xind(j) = int((phaseDM(j,i)-le(j))/dx(j))+1 ! compute its position in the phaseDM array and determine which grid cell (xind) it falls into, index run from 1
    enddo
    dmgrid(xind(1), xind(2), xind(3)) = dmgrid(xind(1), xind(2), xind(3)) + 1.0 ! count total # of DM particle in grid
  enddo
  
  ! make a kdtree of DM particles
  call kd_init(kdroot, phaseDM, bs)
  write(0,*) 'Done KD tree construction'
  
  ! sample PBHs, could be altered to work on different PBH senario
  write(0,*) 'Adding PBHs...'
  allocate(lpbh(1:6,2*Ncell))
  do i=1,nx
    do j=1,nx
      do k=1,nx ! traverse 3-d grid cells
        if (dmgrid(i,j,k) .gt. 0) then
        Nthis = poisson(dmgrid(i,j,k)/float(Ndm)*Npbh, rng_obj)! # of PBHs generated in that cell
        do ipbh=1,Nthis ! random positions
          rvec(1) = rng_obj%get_ran()
          rvec(2) = rng_obj%get_ran()
          rvec(3) = rng_obj%get_ran()
          lpbh(1:3,Ngen+ipbh) = rvec*dx + dx*(/float(i-1), float(j-1), float(k-1)/) + le
        ! Nearest neighboring points to the PBH in a KD tree (kd_nnearest function) are found, and the velocities of those neighbors are averaged and assigned to the PBH.
          call kd_nnearest(kdroot, lpbh(1:3,Ngen+ipbh), nnn, nnindex, nndis, i_dummy)
          p_dummy = 0.0
          do l=1,nnn ! run over # of neighbors and add velocities
            p_dummy(4:6) = p_dummy(4:6) + phaseDM(4:6,nnindex(l))
          enddo
          p_dummy(4:6) = p_dummy(4:6)/float(nnn)
          lpbh(4:6,Ngen+ipbh) = p_dummy(4:6) !phaseDM(4:6,nnindex(1))
        enddo
        Ngen = Ngen + Nthis ! record total #s of PBHs generated
        !else
        !  write(0,*) 'empty cell:',i,j,k
        endif
      enddo
    enddo
  enddo
  write(0,*) 'Num. of generated PBHs: ',Ngen, ', ratio: ',float(Ngen)/Npbh
  
  if (Ngen .gt. 0) then
  
#ifndef ALLPBH
  if (Ngen .le. npert) then
    write(0,*) 'Only a few PBHs! You can use ALLPBH if fac_pert is large!'
    !stop  
  endif
#endif
!the program outputs the position and velocity ranges of the PBHs and dark matter particles
#ifdef DEBUG
  write(0,*) 'PBH position range:', minval(lpbh(1:3,1:Ngen)), maxval(lpbh(1:3,1:Ngen))
  write(0,*) 'PDM velocity range:', minval(phaseDM(4:6,1:Ndm)), maxval(phaseDM(4:6,1:Ndm))
  write(0,*) 'PBH velocity range:', minval(lpbh(4:6,1:Ngen)), maxval(lpbh(4:6,1:Ngen))
#endif

! write the positions and velocities of PBHs into a file.
#ifdef OUTPBH
  if (Ngen>0) then
  open(3, file=trim(path) // trim(icbase) //'_pbhlist.txt', status='unknown')
    do i=1,Ngen
      write(3,'(E14.6,E14.6,E14.6,E14.6,E14.6,E14.6)') lpbh(:,i) ! data written in exponential notation,  in a field of 14 characters, including 6 digits after the decimal point
    enddo
  close(3)
  endif
#endif

! perturb the velocities of particles around the PBHs.
#ifdef PERTURB
! set the positions of the remaining elements of the PBH list to far away from
! the simulation box, so that they will not be included in perturbations
do i=Ngen+1,2*Ncell 
  lpbh(1:3,i) = (/10*bs,10*bs,10*bs/) ! place PBHs away
enddo

  dx_pert = 0.0
#ifdef VELOCITY
  dv_pert = 0.0
#endif
  fc = (Om0-Ob0)/Om0
  a_ = 0.25*(sqrt(1+24*fc)-1)
  npert=min(npert, Ngen)
  fac = 1.0
  facavg0 = 0.0
  facavg1 = 0.0
  ! perturb existing particles
  aeq = Or0/Om0
  H0 = 100*h*1e5/1e3/kpc
  call kd_init(pbhroot, lpbh, bs) ! initialize kd-tree with list of particles including PBHs
  do i=1,Nhres
    xdis = 0.0 ! calculate displacement of particles from PBHs
#ifdef ALLPBH !  calculates the displacement from all PBHs. Otherwise, it limits the number of PBHs that are considered based on a radius
    do j=1,Ngen
      call dx_periodic(pos(1:3,i), lpbh(1:3,j), disp, bs)
      dummy = dist0(disp,(/0.,0.,0./),bs)
      xdis = xdis + disp/dummy/(dummy**2 + (esft*dismax)**2)
    enddo
#else
    if (i .eq. 1) then
      write(0,*) 'Perturbing particles within', rpert
    endif
    call kd_nnearest(pbhroot, pos(1:3,i), npert, nnindex, nndis, i_dummy) ! find nearest # of particles as npert to perturb the velocity
    do j=1,npert
      if (nndis(j) .lt. rpert) then
        call dx_periodic(pos(1:3,i), lpbh(1:3,nnindex(j)), disp, bs)
        if (nnindex(j) .gt. Ngen) then ! halts the program if any of those neighbors is outside the intended range
          write(0,*) 'Fake particle from kdtree:', nnindex(j), pos(1,i), pos(2,i), pos(3,i)
          stop
        endif
        dummy = dist0(disp,(/0.,0.,0./),bs)
        xdis = xdis + disp/dummy/(dummy**2 + (esft*dismax)**2)  !dist0(pos(1:3,i), lpbh0(1:3,j))**3
      endif
    enddo
#endif
    xdis = xdis * 4*PI*gra*mpbh*msun/kpc/kpc
    !xdis = xdis * 4/(Om0*H0**2)*log(0.5*(1.0+sqrt(1.0+a/aeq)))/kpc
    xdis = xdis * Dgrow(a/aeq, fc, a_)/1.5/(Om0*H0**2)/kpc *h**3
    
    if (i .le. Ngas) then
      mref = massarr(0) ! categorize to gas particle
#ifdef DAMPING ! x displacement scaled by growth factor
      xdis = xdis * Dgrow(a*(1+z_rec), fc, a_)/Dgrow(a/aeq, fc, a_)
#endif
    elseif (i .le. Nhres) then
      mref = massarr(1) ! case of  DM particle
    else
      mref = lmass(i-Nhres)
    endif
    !if (mpbh .lt. mref*1e10/h) then
    !  fac = fac * mpbh*h/1e10/mref
    !endif
#ifdef TRUNCATE
    fac = min(ftrunc*dismax/dist0(xdis,(/0.,0.,0./),bs), 1.0)
    if (i .gt. Ngas) then
      facavg1 = facavg1+fac
   	else  
      facavg0 = facavg0+fac
    endif
#endif
    dx_pert = dx_pert + dist0(xdis,(/0.,0.,0./),bs)*fac
    pos(1:3,i) = pos(1:3,i) - xdis*fac
#ifdef VELOCITY
    dv_pert = dv_pert + dist0(xdis,(/0.,0.,0./),bs) * H0*sqrt(Om0) * kpc/h/1e5/a *fac
    vel(1:3,i) = vel(1:3,i) - xdis *H0*sqrt(Om0) * kpc/h/1e5/a *fac
    ! The internal velocity in GIZMO/GADGET is actually u=v*sqrt(a) 
    ! given the physical velocity v and scale factor a
#endif
    do k=1,3 ! make sure the position is within the box after perturbation
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

#endif 
  endif ! if (Ngen .gt. 0)
#endif

  ! output data
if (Ngen .gt. 0) then
  massarr(1) = massarr(1)*(1.0-fpbh)
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
endif

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
      write(2) lmass
    endif
    if (Ngas .gt. 0) then
      write(2) iegas
    endif
  close(2)

contains
! functions and subroutines used.
! poisson(): generate a Poisson-distributed random number
! dist(): calculate the distance between two points considering periodic boundary conditions
! dx_periodic(): calculate the displacement vector between two points considering periodic boundary conditions
! Dgrow(): calculate the growth factor for a given scale factor

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
#ifdef PERIOD
      if (r2 .gt. bs*0.5) then
        r2 = r2-bs
      endif
      if (r2 .lt. -0.5*bs) then
        r2 = r2+bs
      endif
#endif
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

! contain periodic boundary case
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
    Dgrow = (1.0+1.5*fc/a_*s)**a_ - 1.0
  end function

end program




