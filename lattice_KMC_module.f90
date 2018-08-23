module lattice_KMC

   use execution_queue
   implicit none
   integer                               :: Ns, c, empty, occupied
   integer,  allocatable, dimension(:)   :: reaction_effect
   integer,  allocatable, dimension(:,:) :: lattice, coords, neighbours
   real*8,   allocatable, dimension(:)   :: reaction_const
   real*8,   allocatable, dimension(:,:) :: propensities
   integer*8,allocatable, dimension(:)   :: h
   real*8                                :: ads_const  = 5., des_const= 5.
   real*8                                :: diff_const = 20.

contains
   subroutine init(queue_struct,n,percentage)
      implicit none
!~       type(PropensPartialSums_Type) :: queue_struct
      class(ProcessQueue_Type) ::  queue_struct
      integer :: n!, i
      real*8  :: percentage

      Ns = n*n ! total number of sites

      allocate(lattice(n,n))
      allocate(coords(Ns,2))
      allocate(neighbours(Ns,4))
      allocate(propensities(Ns,6))
      
      call queue_struct%initialize(6*Ns)
      call find_coords(n)
      call find_neighbours(n)
      call randomize_coverage(percentage)

!~       call insert_events_to_sums_tree(queue_struct)!using sums tree
      call insert_events_to_heap(queue_struct)! using binary heap
!      print*, queue_struct%queue_elements, char(10),char(10)
!~       print*, queue_struct%heap_labels 
   end

   subroutine execute_reaction(queue_struct,reaction_occured,t_kmc)
      implicit none
!~       type(PropensPartialSums_Type) :: queue_struct
      class(ProcessQueue_Type) ::  queue_struct
      real*8     :: t_kmc
      integer    :: reaction_occured, site_or, site_des, r_type!,i,j ! site of origin & destination

      !specify which site is affected and the reaction_type (diff,des,ads)
      site_or = int((reaction_occured-1)/6) + 1 ! ∈ [1, Ns]
      r_type  = mod( reaction_occured-1, 6) + 1 ! ∈ [1,  6]

      ! find the destination site and update the state of the lattice
      ! update matrix of propensites that lattice_KMC "sees" AND do the checks
      select case(r_type)
         case(1:4)! diffusions
            site_des = neighbours(site_or,r_type) ! neighbour according to reaction_type
            lattice(coords(site_or, 1), coords(site_or, 2)) = 0 ! change state of ORIGIN site
            lattice(coords(site_des,1), coords(site_des,2)) = 1 ! change state of DESTINATION site
            !----------------ORIGIN SITE-----------------------------------------------------------
            propensities(site_or, 1:4) = 0.0_8     ! further diffusions are disabled
            propensities(site_or, 5  ) = 0.0_8     ! DEsorption is  disabled
            propensities(site_or, 6  ) = ads_const ! ADsorption is   enabled
            call enable_diffusions(site_or)    ! diffusion BACK is permitted by default but needs special treatment
            !----------------DESTINATION SITE------------------------------------------------------
            propensities(site_des, 5 ) = des_const ! DEsorption is  enabled
            propensities(site_des, 6 ) = 0.0_8     ! ADsorption is disabled
            call check_permitted_diffusions(site_des)
         case(5)! DEsorption
            site_des = site_or !! needed
            lattice(coords(site_or, 1), coords(site_or, 2)) = 0 ! change state of ORIGIN site
            propensities(site_or, 1:4     ) = 0.0_8     ! diffusions are disabled
            propensities(site_or, r_type  ) = 0.0_8     ! DEsorption is  disabled
            propensities(site_or, r_type+1) = ads_const ! ADsorption is   enabled
            call enable_diffusions(site_or)!diffusion of the neighbours to the empty site HAS TO BE enabled
         case(6)! ADsorption
            site_des = site_or !! needed
            lattice(coords(site_or, 1), coords(site_or, 2)) = 1
            call check_permitted_diffusions(site_or) ! OK
            propensities(site_or, r_type-1) = des_const ! DEsorption is  enabled
            propensities(site_or, r_type  ) = 0.0_8     ! ADsorption is disabled
      end select

      ! update queuing structure (tree or heap) according to NEW positions/propensities
!~       call update_sums_tree(queue_struct, site_or, site_des)
      call update_heap(queue_struct, site_or, site_des,t_kmc)
!~       print*, queue_struct%heap_elements, char(10)
   end

!----------------------- Internally Used Subroutines ----------------------
   !subroutine insert_events_to_sums_tree(queue_struct)
   !   implicit none
   !   type(PropensPartialSums_Type) :: queue_struct
   !   integer :: i
   !   ! ALL events are inserted in the tree, even those with zero propensity
   !   do i=1, Ns ! insert(inserted_element,inserted_details_in (OPTIONAL) )
   !      call queue_struct%insert(propensities(i,1)) ! diffusion North
   !      call queue_struct%insert(propensities(i,2)) ! diffusion South
   !      call queue_struct%insert(propensities(i,3)) ! diffusion East
   !      call queue_struct%insert(propensities(i,4)) ! diffusion West
   !      call queue_struct%insert(propensities(i,5)) ! desorption
   !      call queue_struct%insert(propensities(i,6)) ! adsorption 
   !   enddo
   !end

   !subroutine update_sums_tree(queue_struct, site_or, site_des)
   !   implicit none
   !   type(PropensPartialSums_Type) :: queue_struct
   !   integer :: i,j, site_or, site_des
   !   ! update(updated_element_label,updated_element_value,updated_element_details (OPTIONAL))
   !   do i=1,4
   !      call queue_struct%update((site_or -1)*6 + i, propensities(site_or, i))
   !      call queue_struct%update((site_des-1)*6 + i, propensities(site_des,i))
   !      do j=1,4
   !         call queue_struct%update((neighbours(site_or, j)-1)*6+i, propensities(neighbours(site_or, j),i) )
   !         call queue_struct%update((neighbours(site_des,j)-1)*6+i, propensities(neighbours(site_des,j),i) )
   !      enddo
   !   enddo
   !   do i=5,6
   !      call queue_struct%update((site_or -1)*6+i, propensities(site_or, i) )
   !      call queue_struct%update((site_des-1)*6+i, propensities(site_des,i) )
   !   enddo
   !
   !end

   subroutine insert_events_to_heap(queue_struct)
      implicit none
      class(ProcessQueue_Type) ::  queue_struct
      real*8, dimension(6) :: rand_nums
      integer :: i,j
      ! binary heap stores the occurence times of the events.
      do i=1, Ns ! insert(inserted_element,inserted_details_in (OPTIONAL) )
         call random_number(rand_nums)
         do j=1,6
            call queue_struct%insert(-LOG(rand_nums(j)) / propensities(i,j))
         enddo
      enddo
   end

   subroutine update_heap(queue_struct,site_or, site_des, t_kmc)
      implicit none
      class(ProcessQueue_Type) ::  queue_struct
      real*8, dimension(-1:9,4) :: rand_nums
      real*8  :: t_kmc
      integer :: i,j, site_or, site_des
      call random_number(rand_nums)
      do i=1,4 ! update origin & destination propensities 
         call queue_struct%update((site_or -1)*6 + i, t_kmc - LOG(rand_nums(-1,i)) / propensities(site_or, i))
         call queue_struct%update((site_des-1)*6 + i, t_kmc - LOG(rand_nums( 0,i)) / propensities(site_des,i))
         do j=1,4 ! update neighbours of origin & destination sites
            call queue_struct%update((neighbours(site_or, j)-1)*6+i, &
            t_kmc - LOG(rand_nums(2*j-1,i)) / propensities(neighbours(site_or, j),i) )
            call queue_struct%update((neighbours(site_des,j)-1)*6+i, &
            t_kmc - LOG(rand_nums(2*j  ,i)) / propensities(neighbours(site_des,j),i) )
         enddo
      enddo
      do i=5,6 ! update AD & DES propensities of origin and destination sites
         j = i-4 ! j=1,2
         call queue_struct%update((site_or -1)*6+i, t_kmc - LOG(rand_nums(9, 2*j-1)) / propensities(site_or, i) ) ! rands = (9,1) (9,3)
         call queue_struct%update((site_des-1)*6+i, t_kmc - LOG(rand_nums(9, 2*j  )) / propensities(site_des,i) ) ! rands = (9,2) (9,4)
      enddo
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
            propensities(site_of_interest,i) = diff_const
         else ! if lattice site is already occupied
            propensities(site_of_interest,i) = 0.0_8
            ! need to disable mirrored reaction, 1-2, 2-1, 3-4, 4-3
            propensities(neighbour,mir(i))=0.0_8
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
            propensities(neighbour,mir(i)) = diff_const
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
      propensities(:, 1:5) = 0.0_8     ! diffusions + DEsorption
      propensities(:, 6  ) = ads_const ! ADsorption
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
            propensities(site, 5) = des_const ! DEsorption is  enabled
            propensities(site, 6) = 0.0_8     ! ADsorption is disabled
         enddo
      endif
   end

end module lattice_KMC
