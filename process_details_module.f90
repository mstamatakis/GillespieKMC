module process_details_module
    
    private
    
    integer :: elemstepnmxsites
    logical :: isInit=.FALSE.
    
    type QProcessDetails_Type
        integer, public :: eventtype ! type of elementary event
        integer, dimension(:), allocatable, public :: sites ! sites on which it will/may happen
        real(8), public :: propenst0 ! propensity at initial time
        real(8), public :: deltaenrg ! delta energy of event (final-initial state)
    end type
    
    type(QProcessDetails_Type) :: QProcsDetails
    
    public QProcsDetails, nullProcessDetails, QProcessDetails_Type
    
    contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    subroutine populate_process_details_module_commons()
    
    implicit none

    elemstepnmxsites = 2
    isInit=.TRUE.
    
    return
    
    end subroutine populate_process_details_module_commons

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    function nullProcessDetails()

    implicit none
        
    type(QProcessDetails_Type) :: nullProcessDetails
    
    real(8), parameter :: d_QNaN = transfer(z'7ff8000000000000',1.0_8)

    ! Make sure common private variables are initialised
    if (.not.isInit) then
        call populate_process_details_module_commons()
    endif

    allocate(nullProcessDetails%sites(elemstepnmxsites))
    nullProcessDetails%eventtype = 0
    nullProcessDetails%sites = 0
    nullProcessDetails%propenst0 = d_QNaN
    nullProcessDetails%deltaenrg = d_QNaN

    return

    end function
    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
end module process_details_module