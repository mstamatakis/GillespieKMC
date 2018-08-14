!> Contains structures and methods for an execution queue implementation.
module execution_queue

!> The default implementation is a simple collection of arrays in which 
!! the search of the imminent process is dones by a brute-force minval algorithm.

use parser_module, only: int2str
use process_details_module

implicit none

! Make everything private by default
private
! Explictly define public types, data, and procedures
public :: ProcessQueue_Type, d_QNaN, d_QpInf, d_QnInf

! Class-wide parameters
real(8), parameter :: d_QNaN = transfer(z'7ff8000000000000',1.0_8)
real(8), parameter :: d_QpInf = transfer(z'7FF0000000000000',1.0_8)
real(8), parameter :: d_QnInf = transfer(z'FFF0000000000000',1.0_8)

type ProcessQueue_Type
	integer, public :: nsize ! number of elements stored
	real(8), dimension(:), allocatable, public :: queue_elements ! key values of these elements (occurrence times of elementary processes in KMC)
                                                                ! To get the time of occurrence for elementary process ilabel, 
                                                                ! we evaluate: queue_elements(ilabel). Indexes of elements that have value NaN correspond to 
                                                                ! invalid labels
    type(QProcessDetails_Type), dimension(:), allocatable, public :: queue_details ! elementary process details (type defined in process_details_module so the user can put whatever information they want there).
                                                                ! To get the details for lattice provess ilabel, 
                                                                ! we evaluate: queue_details(ilabel)
    contains
	  procedure, public :: initialize
	  procedure, public :: populate
	  ! elements are referenced by their labels in the following functions ->
	  procedure, public :: update
	  procedure, public :: remove
	  procedure, public :: insert
	  procedure, public :: insert_explicit
	  procedure, public :: details_of
	  procedure, public :: key_value_of
      procedure, public :: highest_priority_label
	  ! <- elements are referenced by their labels in preceding functions
	  procedure, public :: check_status
	  procedure, private :: error
end type

contains
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine initialize(this,n0)

implicit none

class(ProcessQueue_Type), intent(inout) :: this

integer n0 ! the maximum number of elements stored in the heap

this%nsize = 0

allocate(this%queue_elements(n0), source = d_QNaN)

allocate(this%queue_details(n0))
this%queue_details = nullProcessDetails()

return

end subroutine initialize

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine populate(this,values,labels,details_in)

implicit none

class(ProcessQueue_Type), intent(inout) :: this
real(8), intent(in) :: values(:)
integer, intent(in) :: labels(:)
type(QProcessDetails_Type), optional :: details_in(:)

integer i, j
type(QProcessDetails_Type) :: details(size(values))

if (size(values) > size(this%queue_elements)) then
	call this%error(881014,'Occurred when trying to populate a queue of size ' // &
	trim(int2str(size(this%queue_elements))) // ' with ' // trim(int2str(size(values))) // &
	' elements.')
endif

if (size(values) /= size(labels)) then
	call this%error(881001)
end if

if (present(details_in)) then
	if (size(values) /= size(details_in)) then
		call this%error(881002)
	end if
	details = details_in
else
	details = nullProcessDetails()
endif

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
    this%queue_details(labels(i)) = details(i)
enddo

return
end subroutine populate

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine update(this,updated_element_label,updated_element_value,updated_element_details)

implicit none

class(ProcessQueue_Type), intent(inout) :: this
integer, intent(in) ::  updated_element_label
real(8), intent(in) ::  updated_element_value
type(QProcessDetails_Type), intent(in), optional ::  updated_element_details

! This subroutine updates item with label updated_element_label
if (isnan(this%queue_elements(updated_element_label))) then
	call this%error(881013,'Element with label ' // &
	trim(int2str(updated_element_label)) // ' cannot be updated.')
endif

this%queue_elements(updated_element_label) = updated_element_value

if (present(updated_element_details)) then
    this%queue_details(updated_element_label) = updated_element_details
endif

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
    this%queue_details(removed_element_label) = this%queue_details(this%nsize)
endif

this%queue_elements(this%nsize) = d_QNaN
this%queue_details(this%nsize) = nullProcessDetails()
this%nsize = this%nsize-1

return

end subroutine remove

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine insert(this,inserted_element,inserted_details_in)

class(ProcessQueue_Type), intent(inout) :: this
real(8), intent(in) :: inserted_element
type(QProcessDetails_Type), optional :: inserted_details_in

integer inserted_label
type(QProcessDetails_Type) :: inserted_details

if (present(inserted_details_in)) then
	inserted_details = inserted_details_in
else
	inserted_details = nullProcessDetails()
endif

inserted_label = this%nsize + 1

call insert_explicit(this,inserted_element,inserted_label,inserted_details)

return

end subroutine insert

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine insert_explicit(this,inserted_element,inserted_label,inserted_details)

implicit none

class(ProcessQueue_Type), intent(inout) :: this
integer, intent(in) :: inserted_label
real(8), intent(in) :: inserted_element
type(QProcessDetails_Type), intent(in) :: inserted_details

if (inserted_label > size(this%queue_elements)) then
	call this%error(881014,'Occurred when trying to add one more element ' // &
			'in a heap of size ' // trim(int2str(size(this%queue_elements))) // '.')
endif
if (.not.isnan(this%queue_elements(inserted_label))) then
	call this%error(881016,'Occurred when trying to add one more element ' // &
			'with the already used label ' // trim(int2str(inserted_label)) // '.')
endif

this%queue_elements(inserted_label) = inserted_element
this%queue_details(inserted_label) = inserted_details

this%nsize = max(inserted_label,this%nsize)

end subroutine insert_explicit

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  

function details_of(this,req_element_label)

implicit none

type(QProcessDetails_Type) :: details_of

class(ProcessQueue_Type), intent(inout) :: this
integer, intent(in) ::  req_element_label

details_of = this%queue_details(req_element_label)

return

end function details_of

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  

real(8) function key_value_of(this,req_element_label)

implicit none

class(ProcessQueue_Type), intent(inout) :: this
integer, intent(in) ::  req_element_label

key_value_of = this%queue_elements(req_element_label)

return

end function key_value_of

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function highest_priority_label(this)

implicit none

class(ProcessQueue_Type), intent(inout) :: this

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

    case default
        write(errorout,'(/,a)') 'Error - Unspecified error code.'
       
end select

write(errorout,'(/,a)') '***************'

write(errorout,'(/,a)') '> ABNORMAL TERMINATION DUE TO FATAL ERROR <'

flush(iwrite)
stop

end subroutine error

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module execution_queue
