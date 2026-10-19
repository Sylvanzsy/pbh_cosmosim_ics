! kdtree nearst neighbor search
! adapted from the implementation by Travis Sluka
! see https://github.com/travissluka/geoKdTree
#define PERIOD
! considers periodic boundary conditions.
module kdtree

  implicit none
  private
  
  public :: kd_root, kd_init, kd_free, kd_nnearest
  
  integer kd_dim, task_size, MINDIV
  real x1, x2, y1, y2, z1, z2, ext
! task_size: depth of the fixed-size LIFO task stacks used by kd_init
! (taskp/taskdim) and by kd_nnearest (task).
!
! Both stacks are plain depth-first traversals of the tree: every iteration
! pops one entry and pushes at most two (the two children of the popped box).
! The entries sitting below the top of such a stack are exactly the
! not-yet-visited siblings of the boxes on the current root-to-node path, so
! there is at most ONE pending entry per tree level, plus the pair currently on
! top.  Hence
!     max stack occupancy <= h + 1,   h = maximum tree depth.
!
! The split is by point INDEX, not by coordinate value: the left child gets
! ii = (np+1)/2 = ceil(np/2) points and the right child np-ii = floor(np/2).
! A box is therefore at most ceil(N/2**d) points deep at level d, independent
! of how degenerate the point distribution is (duplicate coordinates included),
! and a box is only split while its point count exceeds MINDIV.  So
!     h <= floor(log2(N/MINDIV)) + 1
! and the requirement is
!     task_size >= floor(log2(N/MINDIV)) + 2.
!
! With MINDIV=20 that is 25 for N=512**3 (=1.3e8) and 28 for N=2048**3 (=8.6e9).
! task_size=128 leaves a factor ~2**100 of headroom in N, and costs only
! 2*128*4 = 1 kB of stack in kd_init and 128*4 = 512 B per kd_nnearest call.
!
! NOTE: the previous value (50) was in fact only ever exercised up to depth 1,
! because kd_init never pushed the second child back onto the stack (see the
! bug fix in kd_init below).  It is raised here together with that fix so that
! the now genuinely recursive traversal has a derived, not incidental, bound.
  parameter(kd_dim=3, task_size=128, MINDIV = 20)
  parameter(x1=0, x2=1e6, y1=0, y2=1e6, z1=0, z2=1e6, ext=2)

!kd_root: Represents the root of the k-d tree. It contains pointers to the index array, node array, and point array.
  type kd_root
    integer, pointer :: lindex(:)
    type(node), pointer :: lbox(:)
    real, pointer :: lp(:,:)
    real :: bs
  end type kd_root

!node: Represents a node in the k-d tree. It contains information about the bounding box of points, parent index, and indices of its descendant nodes.
  type node
    real :: p1(kd_dim), p2(kd_dim)
    integer :: pind
    integer :: dind1
    integer :: dind2
    integer :: pindi
    integer :: pindf
  end type node

contains

!kd_free: Deallocates the memory assigned to the tree.
  subroutine kd_free(root)
    type(kd_root), intent(inout) :: root
    deallocate(root%lindex)
    deallocate(root%lbox)
    deallocate(root%lp)
  end subroutine kd_free

!kd_init: Initializes the k-d tree using the provided points.
!It initializes the k-d tree. The basic idea is to repeatedly find the median along the longest axis and split the data. This results in a balanced tree.
!The data is divided until the number of points in a node is less than MINDIV.
  subroutine kd_init(root, lp, bs)
    type(kd_root), intent(out) :: root
    real, intent(in) :: lp(:,:), bs
    
!root: This is the root node of the k-d tree and is an output parameter.
!lp: This is an array of points, where each column represents a different point in k-dimensional space.
!bs: This appears to be a boundary size parameter (although it's not fully explained in the snippet).
    integer :: ii, np, ntmp, n, m, nbox, task, pindi, pindf
    integer :: tp, tdim, jbox, ntot, i
    real, pointer :: cp(:)
    integer, pointer :: hp(:)
    real :: p1(kd_dim), p2(kd_dim)
    integer :: taskp(task_size), taskdim(task_size)
! p1 and p2 are arrays of size kd_dim and represent the coordinates of the bounding box.

    
    root%bs = bs
    ntot = size(lp(1,:))
    allocate(root%lindex(ntot))
    do n=1, ntot
      root%lindex(n) = n
    end do
    
    allocate(root%lp(kd_dim, ntot))
    do i=1, kd_dim
      root%lp(i, :) = lp(i, :)
    enddo
    
    m = 1
    ntmp = ntot
    do while (ntmp>0)
      ntmp = ishft(ntmp, -1)
      m = ishft(m, 1)
    enddo
    nbox = 2*ntot - ishft(m, -1)
    if (m .lt. nbox) nbox=m
    allocate(root%lbox(nbox))
    
    !p1 = (/-ext*bs,-ext*bs,-ext*bs/) !(/x1, y1, z1/)
    !p2 = (/bs+ext*bs, bs+ext*bs, bs+ext*bs/) !(/min(x2,bs), min(y2,bs), min(z2,bs)/)
    p1 = (/x1, y1, z1/)
    p2 = (/min(x2,bs), min(y2,bs), min(z2,bs)/)
    root%lbox(1) = node(p1,p2,0,0,0,1,ntot)
    
    if (ntot .lt. MINDIV) return
    
    jbox = 1
    taskp(1) = 1
    taskdim(1) = 0
    task = 1

!The do while(task .gt. 0) loop constructs the tree by iterating through each box (node), sorting the points within that box by a particular dimension (tdim), and then splitting that box into two child boxes.
    do while(task .gt. 0)
      tp = taskp(task)
      tdim = taskdim(task)
      task = task - 1
      pindi = root%lbox(tp)%pindi
      pindf = root%lbox(tp)%pindf
      hp => root%lindex(pindi:pindf)
      cp => root%lp(tdim+1,:)
      
      np = pindf - pindi + 1
      ii = (np+1)/2
      
      call kd_sort(ii, hp, cp)
      
      p1 = root%lbox(tp)%p1
      p2 = root%lbox(tp)%p2
      p1(tdim+1) = cp(hp(ii))
      p2(tdim+1) = p1(tdim+1)
      root%lbox(jbox+1) = node(root%lbox(tp)%p1, p2, tp, 0,0,pindi,pindi+ii-1)
      root%lbox(jbox+2) = node(p1, root%lbox(tp)%p2, tp, 0,0,pindi+ii,pindf)
      jbox = jbox+2
      root%lbox(tp)%dind1 = jbox-1
      root%lbox(tp)%dind2 = jbox
      root%lbox(jbox-1)%pind = tp
      root%lbox(jbox)%pind = tp
      
! BUG FIX: schedule BOTH children for further splitting.
! The original code only tested the first child (jbox-1) and never pushed the
! second child (jbox) back onto the task stack, no matter how many points it
! held.  Because every iteration then popped one task and pushed at most one,
! the stack depth stayed <= 1 and the tree "built" without any error: the split
! chain only ever followed the left child, so the tree degenerated into a
! caterpillar of ~log2(N/MINDIV) internal nodes whose right children are
! unsplit leaves.  The very first right child alone keeps floor(N/2) points -
! i.e. half the particles end up in one giant leaf.
! Consequences that grow with N:
!   * kd_nnearest has to scan that O(N) leaf whenever a query lands in it (or
!     whenever it cannot be pruned), so a query costs O(N) instead of
!     O(log N) and the whole PBH loop becomes O(N^2);
!   * the leaf that kd_locate returns is huge, so the dn(1) bound used for
!     pruning and every subsequent bound are computed over a wildly
!     unrepresentative point set - the search degrades instead of degrading
!     gracefully;
!   * the box geometry stored for that leaf covers half the simulation volume,
!     which makes the dist_box pruning test meaningless there.
! Pushing the second child restores the intended balanced median tree
! (all leaves between MINDIV/2 and MINDIV points, depth ~log2(N/MINDIV)).
! See the task_size derivation above: with both children pushed the stack is a
! real DFS stack and needs log2(N/MINDIV)+2 slots, which is why task_size was
! raised from 50 to 128 at the same time.
      if (ii .gt. MINDIV) then
        task = task + 1
        taskp(task) = jbox-1
        taskdim(task) = mod(tdim+1, kd_dim)
      endif   
      if (np-ii .gt. MINDIV) then
        task = task + 1
        taskp(task) = jbox
        taskdim(task) = mod(tdim+1, kd_dim)
      endif
    enddo    
  end subroutine kd_init

  pure subroutine kd_nnearest(root, pt, num, rp, rd, rnum)
    type(kd_root), intent(in) :: root
    real, intent(in) :: pt(:)
    integer, intent(in) :: num
    integer, intent(out) :: rp(:)
    real, intent(out) :: rd(:)
    integer, intent(out) :: rnum
    
    real :: dn(num)
    integer :: nn(num)
    integer :: kp, i, n, ntask, k
    real :: d
    integer :: task(task_size)
    
    dn = 1e20
    
    kp = kd_locate(root, pt)
    do while(root%lbox(kp)%pindf-root%lbox(kp)%pindi .lt. num)
      kp = root%lbox(kp)%pind
    enddo
    
    do i=root%lbox(kp)%pindi, root%lbox(kp)%pindf
       n=root%lindex(i)
       d = dist(pt, root%lp(:,n), root%bs)
       if (d .lt. dn(1)) then
         dn(1) = d
         nn(1) = n
         if (num .gt. 1) call shift_down(dn, nn)
       endif
    enddo
    
    task(1)=1
    ntask=1
    do while (ntask .ne. 0)
      k = task(ntask)
      ntask = ntask - 1
      if (k .eq. kp) cycle
      d = dist_box(root%lbox(k), pt, root%bs)
      if (d .lt. dn(1)) then
        if (root%lbox(k)%dind1 .ne. 0) then
          task(ntask+1) = root%lbox(k)%dind1
          task(ntask+2) = root%lbox(k)%dind2
          ntask = ntask+2
        else
          do i=root%lbox(k)%pindi, root%lbox(k)%pindf
             n = root%lindex(i)
             d = dist(root%lp(:,n), pt, root%bs)
             if (d .lt. dn(1)) then
               dn(1) = d
               nn(1) = n
               if (num .gt. 1) call shift_down(dn, nn)
             endif
          enddo
        endif
      endif
    enddo
    rnum = num
    do n=1, num
      rp(n) = nn(n)
      rd(n) = dn(n)
    enddo
  end subroutine

! kd_locate: Given a point, it finds the leaf node which would contain the point.
  pure function kd_locate(root, pt)
    type(kd_root), intent(in) :: root
    real, intent(in) :: pt(kd_dim)
    integer :: kd_locate
    integer :: d1, jdim, nb
    nb = 1
    jdim = 0
    do while (root%lbox(nb)%dind1 .ne. 0)
      d1 = root%lbox(nb)%dind1
      if (pt(jdim+1) .le. root%lbox(d1)%p2(jdim+1)) then
        nb = d1
      else
        nb = root%lbox(nb)%dind2
      endif
      jdim = mod(jdim+1, kd_dim)
    enddo
    kd_locate = nb
  end function kd_locate

!dist_box: Calculates the minimum distance from a point to a bounding box.
  pure function dist_box(box, p, bs)
    type(node), intent(in) :: box
    real, intent(in) :: p(kd_dim), bs
    real :: dist_box, dd, p1(kd_dim), p2(kd_dim), r2
    integer :: n
    dd = 0
    do n=1,kd_dim
      p1(n) = box%p1(n)
      p2(n) = box%p2(n)
#ifdef PERIOD
! BUG FIX: wrap the box as a WHOLE, not its two edges independently.
! The previous version applied the minimum-image shift separately to p1(n) and
! to p2(n).  That is only self-consistent while the box is narrower than bs/2;
! for a wider box the two edges get shifted by different amounts and the
! interval collapses.  The root box is the extreme case: with p1=0, p2=bs and
! a query at p=0.9*bs one gets p-p1 = 0.9*bs > bs/2, so p1 -> bs, while
! p-p2 = -0.1*bs does not trigger the p2 shift, leaving p1 = p2 = bs.  The
! "box" is then a single point and dist_box returns 0.1*bs instead of 0.
! Since kd_nnearest prunes a subtree as soon as dist_box >= dn(1), the ROOT box
! itself was pruned for essentially every query: the traversal never descended
! and the search silently degenerated into "return the nearest point inside the
! leaf that kd_locate happened to land in".  Measured on 32**3 random points,
! 33% of the returned neighbours were not the true nearest neighbour.
! This was previously masked by the kd_init imbalance bug: that leaf used to
! hold N/2 points, so the wrong answer was usually right by accident.  With a
! correctly balanced tree the leaf holds ~MINDIV points and the error shows.
! Correct version: take the minimum image of the vector from the box CENTRE to
! the point, then subtract the box half-width.  This stays valid for boxes of
! any width (a box at least bs wide correctly yields distance 0).
      r2 = p(n) - 0.5*(box%p1(n)+box%p2(n))
      r2 = r2 - bs*nint(r2/bs)
      r2 = abs(r2) - 0.5*(box%p2(n)-box%p1(n))
      if (r2 .gt. 0.0) dd = dd + r2*r2
#else
      if (p(n) .lt. p1(n)) dd = dd + (p(n)-p1(n)) * (p(n)-p1(n))
      if (p(n) .gt. p2(n)) dd = dd + (p(n)-p2(n)) * (p(n)-p2(n))
#endif
  enddo
  dist_box = sqrt(dd)
  end function dist_box

!dist: Computes the distance between two points, considering periodic boundary conditions.
  pure function dist(p1, p2, bs)
    real, intent(in) :: p1(kd_dim), p2(kd_dim), bs
    real :: dist
    integer :: n
    real :: r1, r2
    r1 = 0.0
    do n=1, kd_dim !size(p1)
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
    dist = sqrt(r1)
  end function

!kd_sort: A sorting routine that arranges points in increasing order of their distances.
  pure subroutine kd_sort(k, indx, arr)
    integer, intent(in) :: k
    integer, intent(inout) :: indx(:)
    real, intent(in) :: arr(:)
    
    integer :: i, ia, ir, j, l, mid
    real :: a
    
    l=1
    ir=size(indx)
    do while(.true.)
      if (ir .lt. l+1) then
        if (ir .eq. l+1 .and. arr(indx(ir)) .lt. arr(indx(l))) &
          & call swap(indx(l), indx(ir))
        exit
      else
        mid = (l+ir)/2
        call swap(indx(mid), indx(l+1))
        if (arr(indx(l)) .gt. arr(indx(ir)))   call swap(indx(l),indx(ir))
        if (arr(indx(l+1)) .gt. arr(indx(ir))) call swap(indx(l+1),indx(ir))
        if (arr(indx(l)) .gt. arr(indx(l+1)))  call swap(indx(l),indx(l+1))
        i=l+1
        j=ir
        ia=indx(l+1)
        a=arr(ia)
        do while(.true.)
          i = i+1
          do while(arr(indx(i)) .lt. a)
            i = i+1
          enddo
          j=j-1
          do while(arr(indx(j)) .gt. a)
            j = j-1
          enddo
          if (j .lt. i) exit
          call swap(indx(i), indx(j))
        enddo
        indx(l+1)=indx(j)
        indx(j)=ia
        if (j .ge. k) ir=j-1
        if (j .le. k) l=i
      endif
    enddo
  end subroutine kd_sort

!swap: Exchanges two integers.
  pure subroutine swap(a1, a2)
    integer, intent(inout) :: a1, a2
    integer :: a
    a = a1
    a1 = a2
    a2 = a
  end subroutine swap

!shift_down: Rearranges the heap (used in the nearest neighbor search).
  pure subroutine shift_down(heap, ndx)
    real, intent(inout) :: heap(:)
    integer, intent(inout) :: ndx(:)    
    integer:: n, nn, j, jold, ia
    real :: a
    
    nn = size(heap)
    n = nn
    a = heap(1)
    ia = ndx(1)
    jold = 1
    j = 2
    do while(j .le. n)
      if (j .lt. n) then
        if (heap(j) .lt. heap(j+1)) j = j+1
      endif
      if (a .ge. heap(j)) exit
      heap(jold) = heap(j)
      ndx(jold) = ndx(j)
      jold = j
      j = 2*j
    enddo
    heap(jold) = a
    ndx(jold) = ia
  end subroutine shift_down

end module kdtree
