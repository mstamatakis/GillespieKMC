!> Contains structures and methods for an execution queue implementation.
module execution_queue_binary_heap

!> This is implemented as a binary heap. The field "queue_elements" is
!! compulsory and is the key according to which the partial sorting happens. 
!! Every parent has smaller key value than any child in the tree structure. 
!! In Zacros's implementation this is the time of occurrence for a lattice process.
!! The field "heap_labels" returns integers identifying each element; in Zacros's
!! implementation this is the "serial number" (index) of a lattice process.
!! The remaining fields (a) are optional; in Zacros's implementation they
!! provide more details of lattice processes (e.g. type of elementary event,
!! sites on which they happen etc.).
!! Finally, "array_indexes" is the inverse mapping of heap_labels, showing
!! each element's position in the array that stores the elements of the heap.

use parser_module, only: int2str
use process_details_module
use execution_queue

implicit none

! Make everything private by default
private

public :: ProcessQueueBinaryHeap_Type

type, extends(ProcessQueue_Type) :: ProcessQueueBinaryHeap_Type
	! queue_elements is inherited from ProcessQueue_Type. It stores key values of these elements 
    ! (occurrence times of elementary processes in KMC)
    ! To get the time of occurrence for elementary process ilabel, 
    ! we evaluate: queue_elements(array_indexes(ilabel))
	integer, dimension(:), allocatable, public :: heap_labels ! labels of these elements (elementary process indexes)
    ! queue_details is inherited from ProcessQueue_Type. It stores elementary process details 
    ! (type defined in process_details_module so the user can put whatever information they want there).
    ! To get the details for lattice provess ilabel, we evaluate: queue_details(ilabel)
	integer, dimension(:), allocatable, private :: array_indexes ! inverse mapping of labels, to be used privately in the class
    contains
	  procedure, public :: initialize => heap_initialize
	  procedure, public :: populate => heap_populate
      procedure, private :: fullsort => heap_fullsort
	  ! elements are referenced by their labels in the following functions ->
	  procedure, public :: update => heap_update_element
	  procedure, public :: remove => heap_remove_element
	  procedure, public :: insert => heap_insert_element
	  procedure, public :: insert_explicit => heap_insert_element_with_explicit_label
	  procedure, public :: details_of => heap_get_details_of_element
	  procedure, public :: key_value_of => heap_get_key_value_of_element
	  procedure, public :: highest_priority_label => heap_get_label_of_highest_priority_element
      procedure, private :: relabel => heap_relabel
	  ! <- elements are referenced by their labels in preceding functions
	  procedure, public :: check_status => heap_check_status
	  procedure, private :: error => heap_error
end type

contains
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_initialize(this,n0)

implicit none

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this

integer n0 ! the maximum number of elements stored in the heap

this%nsize = 0

allocate(this%queue_elements(n0))
this%queue_elements(n0) = d_QNaN

allocate(this%heap_labels(n0))
this%heap_labels = 0

allocate(this%queue_details(n0))
this%queue_details = nullProcessDetails()

allocate(this%array_indexes(n0))
this%array_indexes = 0

return

end subroutine heap_initialize

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_populate(this,values,labels,details_in)

implicit none

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this
real(8), intent(in) :: values(:)
integer, intent(in) :: labels(:)
type(QProcessDetails_Type), optional :: details_in(:)

integer i, j, n_process
type(QProcessDetails_Type) :: details(size(values))

if (size(values) > size(this%queue_elements)) then
	call this%error(888014,'Occurred when trying to populate a heap of size ' // &
	trim(int2str(size(this%queue_elements))) // ' with ' // trim(int2str(size(values))) // &
	' elements.')
endif

if (size(values) /= size(labels)) then
	call this%error(888001)
end if

if (present(details_in)) then
	if (size(values) /= size(details_in)) then
		call this%error(888002)
	end if
	details = details_in
else
	details = nullProcessDetails()
endif

if (minval(labels) < 1) then
	call this%error(888003,'A label is lower than one.')
elseif (maxval(labels) > size(labels)) then
	call this%error(888003,'A label is higher than size(values).')
else
	do i = 1,size(labels)
		do j = i+1,size(labels)
			if (labels(i) == labels(j)) then
				call this%error(888003,'Value ' // &
				trim(int2str(labels(i))) // ' repeats.')
			endif
		enddo
	enddo
endif

this%nsize = size(values)
do i = 1,this%nsize
    this%queue_elements(i) = values(i)
    this%heap_labels(i) = labels(i)
    this%array_indexes(i) = i
enddo
! queue_details is addressed by label (the lattice process index)
! This is done to avoid queue_details being involved in too many 
! swap operations when (re)sorting the queue. It will be involved
! however in relabelling (when deleting a process).
do i = 1,this%nsize
    this%queue_details(labels(i)) = details(i)
enddo

call heap_fullsort(this)

return
end subroutine heap_populate

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_fullsort(this)

implicit none

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this

integer itemp, i, j
real(8) dtemp

do i = 1,this%nsize
    do j = i+1,this%nsize
        if ( .not.(has_priority(this%queue_elements(i),this%queue_elements(j))) ) then

            dtemp = this%queue_elements(i)
            this%queue_elements(i) = this%queue_elements(j)
            this%queue_elements(j) = dtemp

            itemp = this%heap_labels(i)
            this%heap_labels(i) = this%heap_labels(j)
            this%heap_labels(j) = itemp

            this%array_indexes(this%heap_labels(i)) = i
            this%array_indexes(this%heap_labels(j)) = j

        endif
    enddo
enddo

return
end subroutine heap_fullsort

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_update_element(this,updated_element_label,updated_element_value,updated_element_details)

implicit none

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this
integer, intent(in) ::  updated_element_label
real(8), intent(in) ::  updated_element_value
type(QProcessDetails_Type), intent(in), optional ::  updated_element_details

integer indx, newindx
logical boolv1, boolv2, boolv3, boolv4

! This subroutine updates item with label updated_element_label
! and resorts the heap (using partial sorting)

if (present(updated_element_details)) then
    this%queue_details(updated_element_label) = updated_element_details
endif

indx = this%array_indexes(updated_element_label)

if ( (indx <= this%nsize) .and. (indx >= 1) ) then

    ! boolv1 evaluates to .true. iff
    ! (we are updating an internal element) .AND. 
    !    (the updated element has priority over 
    !     the parent of the updated element)
    boolv1 = (elem_parent(indx) >= 1)
    if (boolv1) then
        boolv1 = has_priority(updated_element_value,this%queue_elements(elem_parent(indx)))
    endif

    if (boolv1) then ! we need to perform upheap
    
        do while ( elem_parent(indx) >= 1 )

            if ( has_priority(updated_element_value,this%queue_elements(elem_parent(indx))) ) then

                this%queue_elements(indx) = this%queue_elements(elem_parent(indx))
                this%heap_labels(indx) = this%heap_labels(elem_parent(indx))

                this%array_indexes(this%heap_labels(elem_parent(indx))) = indx

                indx = elem_parent(indx)

            else

                exit

            endif

        enddo
    
        this%queue_elements(indx) = updated_element_value
        this%heap_labels(indx) = updated_element_label

        this%array_indexes(updated_element_label) = indx

    else  ! we need to perform downheap

        do
        
            ! boolv2 evaluates to .true. iff
            ! (there exists a left child) .AND. 
            !    (the left child of the updated element has priority 
            !    over the updated element of the heap)
            boolv2 = (elem_left_child(indx) <= this%nsize)
            if (boolv2) then
                boolv2 = has_priority(this%queue_elements(elem_left_child(indx)),updated_element_value)
            endif

            ! boolv3 evaluates to .true. iff
            ! (there exists a right child) .AND. 
            !    (the right child of the updated element has priority 
            !    over the updated element of the heap)
            boolv3 = (elem_right_child(indx) <= this%nsize)
            if (boolv3) then
                boolv3 = has_priority(this%queue_elements(elem_right_child(indx)),updated_element_value)
            endif

            if (.not.(boolv2 .or. boolv3)) then
                exit
            endif
    
            ! boolv4 evaluates to .true. iff
            ! (there does not exist a right child) .OR. 
            !    (the left child of the updated element has priority 
            !    over the right child of the updated element)
            boolv4 = (elem_right_child(indx) > this%nsize)
            if (.not.(boolv4)) then
                boolv4 = has_priority(this%queue_elements(elem_left_child(indx)),this%queue_elements(elem_right_child(indx)))
            endif
                            
            if (boolv4) then
                newindx = elem_left_child(indx)
            else
                newindx = elem_right_child(indx)
            endif
            
            this%queue_elements(indx) = this%queue_elements(newindx)
            this%heap_labels(indx) = this%heap_labels(newindx)

            this%array_indexes(this%heap_labels(newindx)) = indx

            indx = newindx

        enddo
        
        this%queue_elements(indx) = updated_element_value
        this%heap_labels(indx) = updated_element_label

        this%array_indexes(updated_element_label) = indx
    
    endif
    
else 
	call this%error(888013,moreinfo='Element with label ' // &
	int2str(updated_element_label) // ' cannot be updated.')
endif

return
end subroutine heap_update_element

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_remove_element(this,removed_element_label)

implicit none

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this
integer, intent(in) :: removed_element_label

integer indx, newindx
logical boolv1, boolv2, boolv3, boolv4

! This subroutine removes item with label removed_element_label
! from the heap. It also swaps the label of the removed element with the 
! highest label (if the removed is not the highest label), to avoid "holes"
! in the structure.

indx = this%array_indexes(removed_element_label)

if ( (indx <= this%nsize) .and. (indx >= 1) ) then

    ! boolv1 evaluates to .true. iff
    ! (we are removing an internal element) .AND. 
    !    (the last element of the heap has priority over 
    !     the parent of the removed element)
    boolv1 = (elem_parent(indx) >= 1)
    if (boolv1) then
        boolv1 = has_priority(this%queue_elements(this%nsize),this%queue_elements(elem_parent(indx)))
    endif

    if (boolv1) then ! we need to perform upheap
    
        do while ( elem_parent(indx) >= 1 )

            if ( has_priority(this%queue_elements(this%nsize),this%queue_elements(elem_parent(indx))) ) then

                this%queue_elements(indx) = this%queue_elements(elem_parent(indx))
                this%heap_labels(indx) = this%heap_labels(elem_parent(indx))

                this%array_indexes(this%heap_labels(elem_parent(indx))) = indx

                indx = elem_parent(indx)

            else

                exit

            endif

        enddo
    
        this%queue_elements(indx) = this%queue_elements(this%nsize)
        this%heap_labels(indx) = this%heap_labels(this%nsize)

        this%array_indexes(this%heap_labels(this%nsize)) = indx

    else  ! we need to perform downheap

        do
        
            ! boolv2 evaluates to .true. iff
            ! (there exists a left child) .AND. 
            !    (the left child of the removed element has priority 
            !    over the last element of the heap)
            boolv2 = (elem_left_child(indx) <= this%nsize)
            if (boolv2) then
                boolv2 = has_priority(this%queue_elements(elem_left_child(indx)),this%queue_elements(this%nsize))
            endif

            ! boolv3 evaluates to .true. iff
            ! (there exists a right child) .AND. 
            !    (the right child of the removed element has priority 
            !    over the last element of the heap)
            boolv3 = (elem_right_child(indx) <= this%nsize)
            if (boolv3) then
                boolv3 = has_priority(this%queue_elements(elem_right_child(indx)),this%queue_elements(this%nsize))
            endif

            if (.not.(boolv2 .or. boolv3)) then
                exit
            endif
    
            ! boolv4 evaluates to .true. iff
            ! (there does not exist a right child) .OR. 
            !    (the left child of the removed element has priority 
            !    over the right child of the removed element)
            boolv4 = (elem_right_child(indx) > this%nsize)
            if (.not.(boolv4)) then
                boolv4 = has_priority(this%queue_elements(elem_left_child(indx)),this%queue_elements(elem_right_child(indx)))
            endif
                            
            if (boolv4) then
                newindx = elem_left_child(indx)
            else
                newindx = elem_right_child(indx)
            endif
            
            this%queue_elements(indx) = this%queue_elements(newindx)
            this%heap_labels(indx) = this%heap_labels(newindx)

            this%array_indexes(this%heap_labels(newindx)) = indx

            indx = newindx

        enddo
        
        if (indx < this%nsize) then

            this%queue_elements(indx) = this%queue_elements(this%nsize)
            this%heap_labels(indx) = this%heap_labels(this%nsize)

            this%array_indexes(this%heap_labels(this%nsize)) = indx

        endif
    
    endif
    
    this%array_indexes(removed_element_label) = 0
    this%queue_elements(this%nsize) = 0
    this%heap_labels(this%nsize) = 0
	! Don't forget that queue_details is indexed by labels:
    this%queue_details(removed_element_label) = nullProcessDetails()

	! Relabel if necessary: the last label will now take the removed label to avoid "holes"
	! (gaps) in the labelling system...
    if (this%nsize /= removed_element_label) then
	    call this%relabel(this%nsize,removed_element_label)
    endif
    
    this%nsize = this%nsize - 1

else 
	call this%error(888013,'Element with label ' // &
	trim(int2str(removed_element_label)) // ' cannot be removed.')
endif

return
end subroutine heap_remove_element

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_insert_element(this,inserted_element,inserted_details_in)

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this
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

call heap_insert_element_with_explicit_label(this,inserted_element,inserted_label,inserted_details)

return

end subroutine heap_insert_element

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_insert_element_with_explicit_label(this,inserted_element,inserted_label,inserted_details)

implicit none

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this
integer, intent(in) :: inserted_label
real(8), intent(in) :: inserted_element
type(QProcessDetails_Type), intent(in) :: inserted_details

integer indx

indx = this%nsize + 1

if (indx > size(this%queue_elements)) then
	call this%error(888014,'Occurred when trying to add one more element ' // &
			'in a heap of size ' // trim(int2str(size(this%queue_elements))) // '.')
endif
if (this%array_indexes(inserted_label) /= 0) then
	call this%error(888016,'Occurred when trying to add one more element ' // &
			'with the already used label ' // trim(int2str(inserted_label)) // '.')
endif

do while (indx > 1)
    
    if (has_priority(inserted_element,this%queue_elements(elem_parent(indx)))) then

        this%queue_elements(indx) = this%queue_elements(elem_parent(indx))
        this%heap_labels(indx) = this%heap_labels(elem_parent(indx))

        this%array_indexes(this%heap_labels(elem_parent(indx))) = indx

        indx = elem_parent(indx)

    else

        exit

    endif

enddo

this%queue_elements(indx) = inserted_element
this%heap_labels(indx) = inserted_label
this%array_indexes(this%heap_labels(indx)) = indx

this%queue_details(inserted_label) = inserted_details

this%nsize = this%nsize + 1

end subroutine heap_insert_element_with_explicit_label

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  

function heap_get_details_of_element(this,req_element_label)

implicit none

type(QProcessDetails_Type) :: heap_get_details_of_element

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this
integer, intent(in) ::  req_element_label

heap_get_details_of_element = this%queue_details(req_element_label)

return

end function heap_get_details_of_element

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  

real(8) function heap_get_key_value_of_element(this,req_element_label)

implicit none

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this
integer, intent(in) ::  req_element_label

heap_get_key_value_of_element = this%queue_elements(this%array_indexes(req_element_label))

return

end function heap_get_key_value_of_element

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  

integer function heap_get_label_of_highest_priority_element(this)

implicit none

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this

heap_get_label_of_highest_priority_element = this%heap_labels(1)

return

end function heap_get_label_of_highest_priority_element

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	  

subroutine heap_relabel(this,cur_element_label,new_element_label)

implicit none

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this

integer, intent(in) ::  new_element_label, cur_element_label

! This subroutine is intended for filling holes in this%array_indexes after 
! removal of an element from the heap. It will re-label heap item
! this%array_indexes(cur_element_label) to new_element_label
! The subroutine will return an error if new_element_label is in use

if (this%array_indexes(new_element_label) /= 0) then

	call this%error(888015,'Attempted to re-label element with label ' // &
	trim(int2str(cur_element_label)) // ' but the new label ' // &
	trim(int2str(new_element_label)) // ' is already in use.')

else

    this%heap_labels(this%array_indexes(cur_element_label)) = new_element_label
    this%array_indexes(new_element_label) = this%array_indexes(cur_element_label)
    this%array_indexes(cur_element_label) = 0
    
	this%queue_details(new_element_label) = this%queue_details(cur_element_label)
	this%queue_details(cur_element_label) = nullProcessDetails()
	
endif

end subroutine heap_relabel

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_check_status(this)

implicit none

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this

integer indx

do indx = 1,this%nsize
    
    if (elem_left_child(indx) <= this%nsize) then

        if ( has_priority(this%queue_elements(elem_left_child(indx)),this%queue_elements(indx)) ) then
			call this%error(888050,'Problem in heap node with index ' // trim(int2str(indx)) // &
			': incorrect sorting - right left has priority!')
        endif

    endif

    if (elem_right_child(indx) <= this%nsize) then

        if ( has_priority(this%queue_elements(elem_right_child(indx)),this%queue_elements(indx)) ) then
			call this%error(888050,'Problem in heap node with index ' // trim(int2str(indx)) // &
			': incorrect sorting - right child has priority!')
        endif

    endif

    if (this%array_indexes(this%heap_labels(indx)) /= indx) then
			call this%error(888050,'Problem in heap node with label ' // &
			trim(int2str(this%heap_labels(indx))) // ': inverse index mapping returns ' // &
            trim(int2str(this%array_indexes(indx))) // ' which is not the actual index ' // &
			trim(int2str(indx)) // '!')
    endif
    
enddo

return

end subroutine heap_check_status	  

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function elem_parent(indx)

implicit none

integer indx

elem_parent = int(floor(dble(indx)/2.D0))

return

end function elem_parent

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function elem_left_child(indx)

implicit none

integer indx

elem_left_child = 2*indx

return

end function elem_left_child

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

integer function elem_right_child(indx)

implicit none

integer indx

elem_right_child = 2*indx + 1

return

end function elem_right_child

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

logical function has_priority(element1,element2)

implicit none

real(8) element1, element2

if (element1 < element2) then
    has_priority = .true.
else
    has_priority = .false.
endif

return

end function has_priority

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine heap_error(this,ierror,moreinfo,myerrorout)

use constants_module, only : iwrite

implicit none

class(ProcessQueueBinaryHeap_Type), intent(inout) :: this
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

    case (888001)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
		': in heap constructor the size of values must be equal to size of labels.'
    case (888002)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
		': in heap constructor the size of values must be equal to size of details.'
    case (888003)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
		': in heap constructor labels must be an array of consecutive integers from 1 to size(values).'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)
    case (888013)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
		': label provided does not point to an element in the heap!'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)
    case (888014)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
		': heap maxed out!'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)
    case (888015)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
		': re-labeling not allowed: new label already in use.'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)
    case (888016)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
		': insertion not allowed: new label already in use.'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)
    case (888050)
        write(errorout,'(/,a)') 'Internal error code ' // trim(int2str(ierror)) // &
		': heap did not pass self-consistency test.'
        write(errorout,'(/,a)') 'More information: '
        write(errorout,'(a)') trim(moreinfo)

    case default
        write(errorout,'(/,a)') 'Error - Unspecified error code.'
       
end select

write(errorout,'(/,a)') '***************'

write(errorout,'(/,a)') '> ABNORMAL TERMINATION DUE TO FATAL ERROR <'

flush(iwrite)
stop

end subroutine heap_error

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

end module execution_queue_binary_heap
