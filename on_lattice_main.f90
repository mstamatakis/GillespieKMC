program on_lattice_MAIN
   use lattice_KMC
   use execution_queue
   use execution_queue_binary_heap
   implicit none
!~    type(PropensPartialSums_Type) :: queue_struct
   class(ProcessQueue_Type), allocatable ::  queue_struct
   
   character(len=20) :: write_file
   real*8    :: t_kmc, t, dt, t1, t2, cov, ts, tf
   integer   :: i, j, iters, reaction_occured, Ldim
   call cpu_time(ts) !-----------------global start time----------------
   write_file= 'lattice.txt'
   write(*,*) 'Give lattice Dimension and # of iterations:'
   read*, Ldim, iters
!~    Ldim = 100
!~    iters = 20
   cov = 0.00 ! fractional lattice coverage: cov ∈ [0, 1]
   write(*,*) 'Which queuing system you want to use? 1 for linear vector, 2 for binary heap:'
   read*, i
   
   select case (i)
   case (1)
       allocate(queue_struct)
   case (2)
       allocate(ProcessQueueBinaryHeap_Type::queue_struct)
   case default
       write(*,*) 'Invalid selection'
       stop
   end select
   
   call init(queue_struct,Ldim,cov)
   t_kmc = 0
   reaction_occured=-1
!~    print*,"START:", reaction_occured, "a0 = ", queue_struct%tree_elements(queue_struct%head_node_indx)
   print*,"START:", reaction_occured ! a0 does not exist in Binary Heap
   open(unit=88,file=write_file)
   do j=1, Ldim
      write(88,*) lattice(j,:)
   enddo
   do i=1, iters
!~       call find_using_sums_tree(queue_struct,dt, reaction_occured)
!~       t = t + dt
      call find_using_execution_queue(queue_struct,t_kmc, reaction_occured)

      if (mod(i,1000) == 0) then
          print*,i, t_kmc, reaction_occured!, char(10)
      endif
      call cpu_time(t1)
      call execute_reaction(queue_struct,reaction_occured,t_kmc)
      call cpu_time(t2)
!~       print*,"time to EXEC next reaction: ",t2-t1
!~       print*,"a0 AFTER reaction execution:",reaction_occured, queue_struct%tree_elements(sums_tree%head_node_indx)
       if (mod(i,5000)==0) then
         do j=1, Ldim
            write(88,*) lattice(j,:)
         enddo
       endif
   enddo
   call cpu_time(tf) !-----------------global finish time---------------
   print*,"Total time: ",tf-ts, " seconds"
   
contains

   !subroutine find_using_sums_tree(queue_struct, dt, reaction_occured)
   !   implicit none
   !   type(PropensPartialSums_Type) :: queue_struct
   !   integer :: reaction_occured
   !   real*8 :: dt, r1, r2, a0, r2a0, t1, t2 
   !   
   !   CALL RANDOM_NUMBER(r1)
   !   CALL RANDOM_NUMBER(r2)
   !   a0 = queue_struct%tree_elements(queue_struct%head_node_indx)
   !   dt = -LOG(r1) / a0
   !   r2a0 = r2 * a0
   !   
   !   call cpu_time(t1)
   !   reaction_occured = queue_struct%first_greater(r2a0)
   !   call cpu_time(t2)
   !   print*,"time to FIND next reaction: ",t2-t1
   !   print*,"a0 of the current lattice state:",queue_struct%tree_elements(queue_struct%head_node_indx)
   !   print*,"reaction chosen:",reaction_occured
   !end

   subroutine find_using_execution_queue(queue_struct,t_kmc,reaction_occured)
      implicit none
      class(ProcessQueue_Type) ::  queue_struct
      integer :: reaction_occured
      real*8  :: t_kmc
      reaction_occured = queue_struct%highest_priority_label()
      t_kmc            = queue_struct%key_value_of(reaction_occured)
   end

end
