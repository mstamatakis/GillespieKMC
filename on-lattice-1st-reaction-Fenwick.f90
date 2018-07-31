module queuing_class
   implicit none
   type, public :: queu
      integer                                    :: nevents
      real*8,  allocatable, dimension(:)         :: propens,bit
      integer, allocatable, dimension(:)         :: label
      integer, allocatable, dimension(:),private :: index   !! rename it
   contains
      procedure, public :: initialise_queu
      procedure, public :: insert
      procedure, public :: update
      procedure, public :: remove
      procedure, public :: get_next_direct_method
      procedure, public :: get_next_direct_method_BIT !! only in BIT
      procedure, public :: init_bit                   !! only in BIT
      procedure, public :: get_next_first_reaction_method
      procedure, public :: prop_sum
      procedure, public :: print_all
   end type

contains
   subroutine initialise_queu(this,n)
      implicit none
      class(queu) :: this
      integer     :: n
      
      this%nevents = 0
      allocate(this%propens(n))
      allocate(this%label  (n))
      allocate(this%index  (n))
      allocate(this%bit    (n)) !! only for BIT
   end

   subroutine insert(this,reaction_rate,reaction_type)
      implicit none
      class(queu) :: this
      real*8  :: reaction_rate
      integer :: reaction_type

      this%nevents               = this%nevents + 1
      this%propens(this%nevents) = reaction_rate
      this%label  (this%nevents) = reaction_type
      this%index(reaction_type)  = this%nevents
   end

   subroutine update(this,new_reaction_rate,reaction_type)
      implicit none
      class(queu) :: this
      real*8  :: new_reaction_rate, add_value
      integer :: reaction_type,   i, length

      !================================================================= ONLY IN BIT
      i = this%index(reaction_type)
      length = size(this%bit)
      add_value = new_reaction_rate - this%propens(this%index(reaction_type))
      do while (i <= length)
         this%bit(i) = this%bit(i) + add_value
         i = i + (iand(i,-i))
      enddo
      !================================================================= ONLY IN BIT
      this%propens(this%index(reaction_type))=new_reaction_rate

   end
   
   subroutine remove(this,event_to_remove)
      implicit none
      class(queu) :: this
      integer :: event_to_remove

      this%propens(this%index(event_to_remove)) = this%propens(this%nevents) ! removed propensity is replaced by the last one
      this%propens(this%nevents) = 0.                                        ! duplicate propensity set to zero, NOT REALLY needed

      this%label(this%index(event_to_remove)) = this%label(this%nevents) ! removed label is replaced by the last one
      this%index(this%label(this%nevents)) = this%index(event_to_remove) ! index of the relocated event type is updated 
      this%label(this%nevents) = 0                                       ! duplicate label set to zero, NOT REALLY needed
      
      this%index(event_to_remove) = 0  ! index of removed event is set to zero, NOT replaced, events are NOT renumbered (should they?)
      this%nevents = this%nevents - 1
   end

   subroutine get_next_first_reaction_method(this,dt,reaction_occured)
      implicit none
      class (queu) :: this
      real*8, dimension(this%nevents) :: rand_nums, event_times
      real*8                          :: dt
      integer                         :: reaction_occured

      call random_number(rand_nums)                         ! generate the random numbers
      event_times(:) = -LOG(rand_nums(:)) / this%propens(:) ! generate the random times τ_i
      call FIND_MIN(event_times, dt, reaction_occured)      ! returns both MIN_VAL and MIN_LOC in one pass
!~       dt = minval(event_times)                      ! finds smallest time
!~       reaction_occured = MINLOC(event_times, DIM=1) ! finds the reaction that corresponds to smallest time
   end

   subroutine get_next_direct_method(this, dt, reaction_occured)
      implicit none
      class(queu) :: this
      real*8      :: a0, r1, r2, dt, r2a0, temp_sum
      integer     :: i, reaction_occured

      a0 = this%prop_sum()
      CALL RANDOM_NUMBER(r1)
      dt = LOG(1/r1) / a0
      CALL RANDOM_NUMBER(r2)
      r2a0 = r2 * a0
      temp_sum=0
      do i=1, this%nevents
         temp_sum = temp_sum + this%propens(i)
         if (temp_sum > r2a0) then
            reaction_occured = this%label(i)
            exit
         endif
      enddo
   end

   subroutine get_next_direct_method_BIT(this, dt, reaction_occured)    !! only in BIT
      implicit none
      class(queu) :: this
      real*8      :: a0, r1, r2, dt, r2a0, temp_sum
      integer     :: i, reaction_occured
      a0 = this%prop_sum()
      CALL RANDOM_NUMBER(r1)
      dt = LOG(1/r1) / a0
      CALL RANDOM_NUMBER(r2)
      r2a0 = r2 * a0 !! this is the sum I need to find
      i = 0       !! could use "reaction_occured" directly
      temp_sum=0
!~       print*, this%bit
!~       print*, r2a0, a0
      do while(temp_sum < r2a0)
         i = i + 1
         call get_sum_up_to_element(this, i, temp_sum)
      enddo
      reaction_occured = i
!~       print*, i, char(10), "============================================"
   end

   subroutine get_sum_up_to_element(this, indx, sum_up_to)       !! only in BIT
      implicit none
      class(queu) :: this
      integer :: i, indx
      real*8  :: sum_up_to
      i = indx
      sum_up_to = 0
      do while(i>0)
         sum_up_to = sum_up_to + this%bit(i)
         i = i - (iand(i,-i))
      enddo
   end

   subroutine init_bit(this) !! only in BIT
      implicit none
      class(queu) :: this
      integer :: length, i, j
      length = size(this%bit)
      do i=1,length
         j=i
         do while(j <= length)
            this%bit(j) = this%bit(j) + this%propens(i)
            j = j + (iand(j,-j))
         enddo
      enddo
!~       print*,"BIT = ", this%bit, CHAR(10)
!~       print*,"PROP= ", this%propens
   end

   function prop_sum(this)
      implicit none
      class(queu) :: this
      real*8      :: prop_sum
      prop_sum = 0
      prop_sum = sum(this%propens(1:this%nevents))
   end

   subroutine print_all(this)
      implicit none
      class(queu) :: this
      print*, "events", this%nevents !, CHAR(10) !to enter a newline character
      print*, "propen", this%propens
      print*, "labels", this%label
      print*, "index ", this%index
      print*, "------------------------------------------------------------------------"
   end
   
   subroutine FIND_MIN(array, min_val, min_loc)
      implicit none
      real*8, dimension(:) :: array
      real*8               :: min_val
      integer              :: min_loc, length, i

      length = size(array)
      min_val = array(1)
      min_loc = 1
      do i=2,length
         if ( array(i) < min_val ) then
            min_val = array(i)
            min_loc = i
         endif
      enddo
   end
   
end module queuing_class

module lattice_KMC
   use queuing_class
   implicit none
   integer                               :: Ns, c, empty, occupied
   integer,  allocatable, dimension(:)   :: reaction_effect
   integer,  allocatable, dimension(:,:) :: lattice, coords, neighbours
   real*8,   allocatable, dimension(:)   :: reaction_const
   real*8,   allocatable, dimension(:,:) :: propensities
   integer*8,allocatable, dimension(:)   :: h

contains
   subroutine init(qu,n,percentage)
      implicit none
      type(queu) :: qu
      integer :: n, i
      real*8  :: percentage

      Ns = n*n ! total number of sites

      allocate(lattice(n,n))
      allocate(coords(Ns,2))
      allocate(neighbours(Ns,4))
      allocate(propensities(Ns,6))
      
      call qu%initialise_queu(6*Ns)
      call find_coords(n)
      call find_neighbours(n)
      call randomize_coverage(percentage)

      ! ALL events are inserted in the queu, even those with zero propensity
      do i=1, Ns  ! signature: insert(propensity, event-type / event-identifier)
         call qu%insert(propensities(i,1),  6*i-5) ! diffusion North
         call qu%insert(propensities(i,2),  6*i-4) ! diffusion South
         call qu%insert(propensities(i,3),  6*i-3) ! diffusion East
         call qu%insert(propensities(i,4),  6*i-2) ! diffusion West
         call qu%insert(propensities(i,5),  6*i-1) ! desorption
         call qu%insert(propensities(i,6),  6*i  ) ! adsorption 
      enddo
      call qu%init_bit()                                                !! only in BIT
   end

   subroutine execute_reaction(qu,reaction_occured)
      implicit none
      type(queu) :: qu
      integer    :: reaction_occured, site_or, site_des, r_type,i,j ! site of origin & destination
!~       real*8     :: des_const, ad_const

      !specify which site is affected and the reaction_type (diff,des,ads)
      site_or = int((reaction_occured-1)/6) + 1 ! ∈ [1, Ns]
      r_type  = mod( reaction_occured-1, 6) + 1 ! ∈ [1,  6]

      ! find the destination site and update the state of the lattice
      ! update matrix of propensites that lattice_KMC "sees" AND do the checks, then pass changes to QUEU
      select case(r_type)
         case(1:4)! diffusions
            site_des = neighbours(site_or,r_type) ! neighbour according to reaction_type
            lattice(coords(site_or, 1), coords(site_or, 2)) = 0 ! change state of ORIGIN site
            lattice(coords(site_des,1), coords(site_des,2)) = 1 ! change state of DESTINATION site
            !----------------ORIGIN SITE-----------------------------------------------------------
            propensities(site_or, 1:4) = 0. _8 ! further diffusions are disabled
            propensities(site_or, 5  ) = 0. _8 ! DEsorption is  disabled
            propensities(site_or, 6  ) = 1. _8 ! ADsorption is   enabled
            call enable_diffusions(site_or)    ! diffusion BACK is permitted by default but needs special treatment
            !----------------DESTINATION SITE------------------------------------------------------
            propensities(site_des, 5 ) = 1. _8 ! DEsorption is  enabled
            propensities(site_des, 6 ) = 0. _8 ! ADsorption is disabled
            call check_permitted_diffusions(site_des)
         case(5)! DEsorption
            site_des = site_or !! needed
            lattice(coords(site_or, 1), coords(site_or, 2)) = 0 ! change state of ORIGIN site
            propensities(site_or, 1:4     ) = 0. _8 ! diffusions are disabled
            propensities(site_or, r_type  ) = 0. _8 ! DEsorption is  disabled
            propensities(site_or, r_type+1) = 1. _8 ! ADsorption is   enabled
            call enable_diffusions(site_or)!diffusion of the neighbours to the empty site HAS TO BE enabled
         case(6)! ADsorption
            site_des = site_or !! needed
            lattice(coords(site_or, 1), coords(site_or, 2)) = 1
            call check_permitted_diffusions(site_or) ! OK
            propensities(site_or, r_type-1) = 1. _8  ! DEsorption is  enabled
            propensities(site_or, r_type  ) = 0. _8  ! ADsorption is disabled
      end select

      ! update QUEU of events according to NEW positions/propensities
      ! signature: update(new_reaction_rate, reaction_type)
      do i=1,4
         call qu%update( propensities(site_or,  i), (site_or -1)*6 + i )
         call qu%update( propensities(site_des, i), (site_des-1)*6 + i )
         do j=1,4
            call qu%update( propensities(neighbours(site_or, j), i), (neighbours(site_or, j)-1)*6 + i )
            call qu%update( propensities(neighbours(site_des,j), i), (neighbours(site_des,j)-1)*6 + i )
         enddo
      enddo
      do i=5,6
         call qu%update( propensities(site_or,  i), (site_or -1)*6 + i )
         call qu%update( propensities(site_des, i), (site_des-1)*6 + i )
      enddo
! have also to update the involved neighbours of the sites:
! 4 if AD or DES occurs, 6 if diffusion occurs.
   end

   subroutine check_permitted_diffusions(site_of_interest)
      implicit none
      integer :: i, site_of_interest, neighbour, row, col
      integer, dimension(4) :: mir=(/2,1,4,3/)

      do i=1,4 !check the 4 neighbours for permitted diffusions
         neighbour = neighbours(site_of_interest, i) ! destination site
         row = coords(neighbour,1)
         col = coords(neighbour,2)
!~          print*,neighbour,lattice(row,col),mir(i),propensities(site_of_interest,i),propensities(neighbour,mir(i))
!~          if(lattice( coords(neighbours(site_of_interest,i),1),coords(neighbours(site_of_interest,i),2) ) == 0) then
         if (lattice(row, col)==0) then
            propensities(site_of_interest,i) = 5.0_8
         else ! if lattice site is already occupied
            propensities(site_of_interest,i) = 0. _8
            ! need to disable mirrored reaction, 1-2, 2-1, 3-4, 4-3
            propensities(neighbour,mir(i))=0. _8
         endif
!~          print*,neighbour,lattice(row,col),mir(i),propensities(site_of_interest,i),propensities(neighbour,mir(i))
      enddo
   end

   subroutine enable_diffusions(site_of_interest)
      implicit none
      integer :: i, site_of_interest, neighbour, row, col
      integer, dimension(4) :: mir=(/2,1,4,3/)
      do i=1,4
         neighbour = neighbours(site_of_interest, i) ! destination site
         row = coords(neighbour,1)
         col = coords(neighbour,2)
         if (lattice(row, col)==1) then
            propensities(neighbour,mir(i))=5. _8
         endif
      enddo
   end

   subroutine find_neighbours(n)
      implicit none
      integer :: n, i
      ! find all NORTH--------------------------------------------------
      do i=1, Ns
         if (mod(i,n)==1) then
            neighbours(i,1) = i + n - 1
         else
            neighbours(i,1) = i - 1
         endif
      enddo
      ! find all SOUTH--------------------------------------------------
      do i=1, Ns
         if (mod(i,n)==0) then
            neighbours(i,2) = i - n + 1
         else
            neighbours(i,2) = i + 1
         endif
      enddo
      ! find all EAST---------------------------------------------------
      do i=1, Ns-n
         neighbours(i,3) = i + n
      enddo
      do i=1, n
         neighbours(Ns-n+i,3) = i
      enddo
      ! find all WEST---------------------------------------------------
      do i=n+1, Ns
         neighbours(i,4) = i - n
      enddo
      do i=1, n
         neighbours(i,4) = Ns - n + i
      enddo
   end

   subroutine find_coords(n)
      implicit none
      integer :: i,n
      do i=1, Ns ! fill coords array
         coords(i,1) = mod( i-1,n)  + 1 ! row of i-lattice site
         coords(i,2) = int((i-1)/n) + 1 ! col of i-lattice site
      enddo      
   end
   
   subroutine randomize_coverage(covrg)
      implicit none
      integer :: i,j, site, c
      real*8  :: randN, covrg
      integer, dimension(Ns) :: rand_positions
      
      ! initialize an empty lattice, only adsorption can occur
      lattice = 0
      propensities(:, 1:5) = 0. _8 ! diffusions + DEsorption
      propensities(:, 6  ) = 1. _8 ! ADsorption
      if ( covrg > 0. ) then
         do i=1, Ns ! fill shuffled array
            call random_number(randN)
            j = 1 + FLOOR(i*randN)
            if(j /= i) rand_positions(i) = rand_positions(j)
            rand_positions(j) = i
         enddo
         c = NINT(covrg*Ns) ! number of randomly selected iccupied 
         print*,"Initially Occupied: ", c, " out of", Ns
         do i=1, c ! insert c sites in lattice
            site = rand_positions(i)
!~             print*,"Site Occupied: ", site
            lattice(coords(site,1), coords(site,2)) = 1
            call check_permitted_diffusions(site)
            propensities(site, 5) = 1. _8 ! DEsorption is  enabled
            propensities(site, 6) = 0. _8 ! ADsorption is disabled
         enddo
      endif
   end

end

program on_lattice
   use queuing_class
   use lattice_KMC
   implicit none
   type(queu) :: qu
   character(len=20) :: read_file
   character(len=20) :: write_file
   real*8    :: t, dt, t1, t2, cov, ts, tf!, r1, r2
   integer   :: i, j, iters, reaction_occured, Ldim
   call cpu_time(ts) !-----------------global start time----------------
   read_file= 'init.txt'
   write_file= 'lattice.txt'
   Ldim = 20
   cov = 0.00 ! fractional lattice coverage: cov ∈ [0, 1]
   call init(qu,Ldim,cov)
   t=0
   iters = 10
   reaction_occured=-1
   print*,"START:", reaction_occured, "a0 = ", qu%prop_sum()
   open(unit=88,file=write_file)
   do j=1, Ldim
      write(88,*) lattice(j,:)
   enddo
   do i=1, iters
      call cpu_time(t1)
!~       call qu%get_next_first_reaction_method(dt,reaction_occured)
      call qu%get_next_direct_method(dt, reaction_occured)
!~       call qu%get_next_direct_method_BIT(dt, reaction_occured)
      call cpu_time(t2)
!~       print*,"time to FIND next reaction: ",t2-t1,i,iters
      t = t + dt

!~       print*,"BEFORE:",reaction_occured, qu%prop_sum()
      call cpu_time(t1)
      call execute_reaction(qu,reaction_occured)
      call cpu_time(t2)
!~       print*,"time to EXEC next reaction: ",t2-t1
!~       print*,"AFTER:",reaction_occured, qu%prop_sum()
!~       if (mod(i,5000)==0) then
         do j=1, Ldim
            write(88,*) lattice(j,:)
         enddo
!~       endif
   enddo
   call cpu_time(tf) !-----------------global finish time---------------
   print*,"Total time: ",tf-ts, " seconds"
end
