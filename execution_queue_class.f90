!> Contains structures and methods for an execution queue implementation.
module execution_queue

!> The default implementation is a simple collection of arrays in which 
!! the search of the imminent process is dones by a brute-force minval algorithm.

use parser_module, only: int2str, remove_leading_blanks, striccompare

implicit none

! Make everything private by default
private
! Explictly define public types, data, and procedures
public :: ProcessQueue_Type, d_QNaN, d_QpInf

! Class-wide parameters
real(8), parameter :: d_QNaN = transfer(z'7ff8000000000000',1.0_8)
real(8), parameter :: d_QpInf = transfer(z'7FF0000000000000',1.0_8)

type ProcessQueue_Type
   integer, public :: nsize ! number of elements stored
   real(8), dimension(:), allocatable, public :: queue_elements ! key values of these elements (occurrence times of elementary processes in KMC)
                                                                ! To get the time of occurrence for elementary process ilabel, 
                                                                ! we evaluate: queue_elements(ilabel). Indexes of elements that have value NaN correspond to 
                                                                ! invalid labels
    contains
     procedure, public :: initialize
     procedure, public :: populate
     ! elements are referenced by their labels in the following functions ->
     procedure, public :: update
     procedure, public :: remove
     procedure, public :: insert
     procedure, public :: insert_explicit
     procedure, public :: key_value_of
      procedure, public :: highest_priority_label
     ! <- elements are referenced by their labels in preceding functions
     procedure, public :: check_status
     procedure, private :: error
      procedure, public :: save_restart_info
      procedure, public :: load_restart_info
end type

contains
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine initialize(this,n0)

implicit none

class(ProcessQueue_Type), intent(inout) :: this

integer n0 ! the maximum number of elements stored in the heap

this%nsize = 0

allocate(this%queue_elements(n0), source = d_QNaN)

return

end subroutine initialize

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine populate(this,values,labels)

implicit none

class(ProcessQueue_Type), intent(inout) :: this
real(8), intent(in) :: values(:)
integer, intent(in) :: labels(:)

integer i, j

if (size(values) > size(this%queue_elements)) then
   call this%error(881014,'Occurred when trying to populate a queue of size ' // &
   trim(int2str(size(this%queue_elements))) // ' with ' // trim(int2str(size(values))) // &
   ' elements.')
endif

if (size(values) /= size(labels)) then
   call this%error(881001)
end if

if (minval(labels) < 1) then
   call this%error(881003,'A label is lower than one.')
elseif (maxval(labels) > size(labels)) then
   call this%error(881003,'A label is higher than size(values).')
else
   do i = 1,size(labels)
      do j = i+1,size(labels)
         if (labels(i) == labels(j)) then
            call this%error(881003,'Value ' // &
            trim(int2str(labels(i))) // ' repeats.')
         endif
      enddo
   enddo
endif

this%nsize = maxval(values)
do i = 1,this%nsize
    this%queue_elements(labels(i)) = values(i)
enddo

return
end subroutine populate

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine update(this,updated_element_label,updated_element_value)

implicit none

class(ProcessQueue_Type), intent(inout) :: this
integer, intent(in) ::  updated_element_label
real(8), intent(in) ::  updated_element_value

! This subroutine updates item with label updated_element_label
if (isnan(this%queue_elements(updated_element_label))) then
   call this%error(881013,'Element with label ' // &
   trim(int2str(updated_element_label)) // ' cannot be updated.')
endif

this%queue_elements(updated_element_label) = updated_element_value

return

end subroutine update

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine remove(this,removed_element_label)

implicit none

class(ProcessQueue_Type), intent(inout) :: this
integer, intent(in) :: removed_element_label

! This subroutine removes item with label updated_element_label
if (isnan(this%queue_elements(removed_element_label))) then
   call this%error(881013,'Element with label ' // &
   trim(int2str(removed_element_label)) // ' cannot be removed.')
endif

if (this%nsize > removed_element_label) then ! copy the last element into the one to be removed (thereby eliminating the latter and makign the last element to get the label of the removed one)
    this%queue_elements(removed_element_label) = this%queue_elements(this%nsize)
endif

this%queue_elements(this%nsize) = d_QNaN
this%nsize = this%nsize-1

return

end subroutine remove

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine insert(this,inserted_element)

class(ProcessQueue_Type), intent(inout) :: this
real(8), intent(in) :: inserted_element

integer inserted_label

inserted_label = this%nsize + 1

call insert_explicit(this,inserted_element,inserted_label)

return

end subroutine insert

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine insert_explicit(this,inserted_element,inserted_label)

implicit none

class(ProcessQueue_Type), intent(inout) :: this
integer, intent(in) :: inserted_label
real(8), intent(in) :: inserted_element

if (inserted_label > size(this%queue_elements)) then
   call this%error(881014,'Occurred when trying to add one more element ' // &
         'in a heap of size ' // trim(int2str(size(this%queue_elements))) // '.')
endif
if (.not.isnan(this%queue_elements(inserted_label))) then
   call this%error(881016,'Occurred when trying to add one more element ' // &
         'with the already used label ' // trim(int2str(inserted_label)) // '.')
endif

this%queue_elements(inserted_label) = inserted_element

this%nsize = max(inserted_label,this%nsize)

end subroutine insert_explicit

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   

real(8) function key_value_of(this,req_element_label)

implicit none

class(ProcessQueue_Type), intent(inout) :: this
integer, intent(in) ::  req_element_label

if (req_element_label == 0) then
    key_value_of = huge(1.d0)
    return
endif

key_value_of = this%queue_elements(req_element_label)

return

end function key_value_of

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function highest_priority_label(this)

implicit none

class(ProcessQueue_Type), intent(inout) :: this

if (this%nsize == 0) then
    highest_priority_label = 0
    return
endif

highest_priority_label = minloc(this%queue_elements(1:this%nsize),1, &
    mask=.not.isnan(this%queue_elements(1:this%nsize)))

return

end function highest_priority_label

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine check_status(this)

implicit none

class(ProcessQueue_Type), intent(inout) :: this

integer indx

if (.not.all(isnan(this%queue_elements(this%nsize:)))) then
    call this%error(881050,'Problem in heap node with index ' // trim(int2str(indx)) // &
    ': incorrect sorting - right left has priority!')
endif

return

end subroutine check_status     

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine error(this,ierror,moreinfo,myerrorout)

use constants_module, only : iwrite

implicit none

class(ProcessQueue_Type), intent(inout) :: this
character(*), optional, intent(in) :: moreinfo
integer, optional, intent(in) :: myerrorout

integer ierror
integer errorout

if (present(myerrorout)) then
    errorout = myerrorout
else
    errorout = iwrite
endif

write(errorout,'(/,a)') '***************'

select case (ierror)

    case (881001)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
      ': in queue constructor the size of values must be equal to size of labels.'
    case (881002)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
      ': in queue constructor the size of values must be equal to size of details.'
    case (881003)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
      ': in queue constructor labels must be an array of consecutive integers from 1 to size(values).'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)
    case (881013)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
      ': queue provided does not point to an element in the heap!'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)
    case (881014)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
      ': queue maxed out!'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)
    case (881015)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
      ': re-labeling not allowed: new label already in use.'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)
    case (881016)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
      ': insertion not allowed: new label already in use.'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)
    case (881050)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
      ': queue did not pass self-consistency test.'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)
    case (881060)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
      ': attempted to save queue restart information but failed.'
    case (881061)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
      ': attempted to load queue restart information but failed.'

    case default
        write(errorout,'(/,a)') 'Error - Unspecified error code.'
       
end select

write(errorout,'(/,a)') '***************'

write(errorout,'(/,a)') '> ABNORMAL TERMINATION DUE TO FATAL ERROR <'

flush(iwrite)
stop

end subroutine error

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine save_restart_info(this,irestart)

implicit none

class(ProcessQueue_Type), intent(inout) :: this
integer, intent(in) :: irestart
integer i

write(irestart,'(2I10)',err=100) this%nsize, size(this%queue_elements)

do i = 1,this%nsize
    write(irestart,'(ES32.16E3)',err=100) this%queue_elements(i)
enddo

return

100 call this%error(881060)

end subroutine save_restart_info

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine load_restart_info(this,irestart)

implicit none

class(ProcessQueue_Type), intent(inout) :: this
integer, intent(in) :: irestart
integer i, io, nsize, n0
character(256) buf256
character(256) recinput

read(irestart,'(2I10)',err=101) nsize, n0

call this%initialize(n0)
this%nsize = nsize

recinput = ' '
buf256 = ' '

do i = 1,this%nsize
    read(irestart,'(a' // int2str(len(recinput)) // ')', iostat = io) recinput
    buf256 = recinput(1:32)
    call remove_leading_blanks(buf256)
    if (striccompare(trim(buf256),'inf') .or. striccompare(trim(buf256),'infinity')) then
        ! MA -- not portable. Replacing with large number.
        this%queue_elements(i) = huge(1.d0)
    else
        read(recinput,'(ES32.16E3)',err=101) this%queue_elements(i)
    endif
    
enddo

return

101 call this%error(881061)

end subroutine load_restart_info

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module execution_queue
