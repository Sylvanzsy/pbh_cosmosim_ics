! This is an object oriented version of the ran3 algorithm from numerical recipes
! Before using a rng type object it needs to be initialized by calling rng%init(seed) where seed is an integer seed
! Afterwards rng%get_ran() will return a random real in the range [0,1)
#define RNG_SEED 161803398
module random_object

  implicit none
  type rng

    integer IFF,INEXT,INEXTP,MA(55),MJ
  contains
    procedure, pass(this) :: get_ran
    procedure, pass(this) :: init
  end type
contains
  subroutine init(this, IDUM)
      class(rng) :: this
      integer, intent(in) :: IDUM
      integer MBIG,MZ
      integer :: jran3 ,MK ,II ,K
      parameter (MBIG=1000000000,MZ=0)
      this%IFF = 0
      this%INEXT = 0
      this%INEXTP = 0
      this%MA = 0
      this%MJ = 0
        
      this%IFF=1
      this%MJ=RNG_SEED-iabs(IDUM)
      this%MJ=mod(this%MJ,MBIG)
      if (this%MJ.lt.MZ) this%MJ=this%MJ+MBIG
      this%MA(55)=this%MJ
      MK=1
      do jran3=1,54
        II=mod(21*jran3,55)
        this%MA(II)=MK
        MK=this%MJ-MK
        if (MK.lt.MZ) MK=MK+MBIG
        this%MJ=this%MA(II)
      end do
      do K=1,4
        do jran3=1,55
          this%MA(jran3)=this%MA(jran3)-this%MA(1+mod(jran3+30,55))
          if (this%MA(jran3).lt.MZ) this%MA(jran3)=this%MA(jran3)+MBIG
        end do
      end do
      this%INEXT=0
      this%INEXTP=31
  end subroutine
  real function get_ran(this)
      implicit none
      !
      ! Integers
      class(rng) :: this
      integer MBIG,MZ
      ! Floats
      real FAC
      !
      ! Parameters
      parameter (MBIG=1000000000,FAC=1.0/MBIG,MZ=0)
      
      this%INEXT=this%INEXT+1
      if (this%INEXT.eq.56) this%INEXT=1
      this%INEXTP=this%INEXTP+1
      if (this%INEXTP.eq.56) this%INEXTP=1
      this%MJ=this%MA(this%INEXT)-this%MA(this%INEXTP)
      if (this%MJ.lt.MZ) this%MJ=this%MJ+MBIG
      this%MA(this%INEXT)=this%MJ
      get_ran=this%MJ*FAC
      return
  end function get_ran

end module random_object 
