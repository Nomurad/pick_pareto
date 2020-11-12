module readcsv
    implicit none
    real(8), allocatable :: csvdata(:,:)
    logical, allocatable :: mask(:)

    contains
    subroutine read_csv(fname)
        implicit none 
        character(*), intent(in) :: fname
        character(256) :: line
        integer :: ios
        integer :: i, j, dim_counter, n_rows
        integer :: unitnum = 11
        integer :: skip_rows = 1
        
        open(unitnum, file=trim(fname), status="old")
            read(unitnum, "(a)", iostat=ios) line
            read(unitnum, "(a)", iostat=ios) line
            !print *, line 

            i = 1
            dim_counter = 1
            do
                j = index(line(i:), ",")
                if (j == 0) exit
                dim_counter = dim_counter + 1
                i = i + j 
            end do
            dim_counter = dim_counter + 1
            print *, "num of values  = ", dim_counter
            
            n_rows = 1
            do while(ios == 0)
                read(unitnum, "(a)", iostat=ios) line 
                n_rows = n_rows + 1 
            end do
            !print *, "reach EOF."
            print *, n_rows, " rows"
            
            allocate(csvdata(n_rows, dim_counter))
            rewind(unitnum)
            read(unitnum, *)

            do i = 1, n_rows
                read(unitnum, "(a)", iostat=ios) line 
                read(line, *, iostat=ios) csvdata(i, :)
                ! print *, csvdata(i, :)
            end do

        close(unitnum)
        
    end subroutine read_csv

end module readcsv


module pick_pareto
    use omp_lib
    implicit none 
    private
    type :: individual
        real(8),allocatable :: values(:)
    end type
    integer :: n_pareto
    real(8), allocatable :: pareto_data(:)
    type(individual), allocatable :: indiv(:)
    type(individual), allocatable :: feasible_indiv(:)

    public :: indiv, feasible_indiv, n_pareto
    public :: individual
    public :: pick_feasible_solution, pick_pareto_solution
    public :: qsort_indiv, de_duplication

    contains
    subroutine pick_feasible_solution(data, objdim, retdata)
        implicit none 
        real(8), intent(in) :: data(:,:)
        integer, intent(in) :: objdim
        real(8), allocatable, intent(out) :: retdata(:,:)

        integer :: i, j 
        integer :: datasize(2), retsize(2)
        integer :: n_feasible = 0
        real(8), allocatable :: subdata(:)
        logical, allocatable :: parammask(:), indivmask(:)
        logical, allocatable :: indivmask_multi(:,:)
        type(individual), allocatable :: tmpindiv(:)

        datasize = shape(data)
        ! print *, datasize
        allocate(subdata(datasize(2)))
        allocate(indivmask(datasize(1)))
        allocate(parammask(datasize(2)-(objdim*2)))
        allocate(indivmask_multi(datasize(1),datasize(2)))

        call alloc_indivs(indiv, datasize(1), datasize(2))

        parammask = .false.
        indivmask = .false.

        print *, "objdim*2 + 1 = ", objdim*2 + 1
        !do i = datasize(1)-1, datasize(1)
        do i = 1, datasize(1)
            indiv(i)%values = data(i,:)
            subdata = indiv(i)%values
            parammask = (subdata(objdim*2+1:)<=0)
            indivmask(i) = all(parammask)
            ! if (indivmask(i) .eqv. .true.) then
            !     n_feasible = n_feasible + 1
            ! end if
            indivmask_multi(i, :) = indivmask(i)
            ! print *, indivmask_multi(i,:)
            ! print *, parammask 
        end do

        n_feasible = count(indivmask)
        print *, n_feasible, count(indivmask)
        retsize = shape(reshape(pack(data, indivmask_multi), (/n_feasible, datasize(2)/)))
        print *, "num of feasible indivs", retsize(1)
        allocate(retdata(retsize(1), retsize(2)))
        retdata = reshape(pack(data, indivmask_multi), retsize)
        !print *, retdata(1,:)
        call alloc_indivs(feasible_indiv, n_feasible, datasize(2))
        feasible_indiv = pack(indiv, indivmask)
        

    end subroutine pick_feasible_solution

    subroutine alloc_indivs(idv, d1, d2)
        type(individual), allocatable, intent(inout) :: idv(:)
        integer, intent(in) :: d1, d2
        integer :: i 

        allocate(idv(d1))
        do i = 1, d1 
            allocate(idv(i)%values(d2))
        end do

    end subroutine alloc_indivs

    ! subroutine get_maskedarray(arr, mask, ret)
    !     implicit none
    !     real(8), intent(in) :: arr(:,:)
    !     logical, intent(in) :: mask(:), mask_multi
    !     real(8), allocatable, intent(out) :: ret(:,:)
    !     
    !     integer :: i
    !     integer :: n_true
    !     integer :: arrsize(2), retsize(2)

    !     arrsize = shape(arr)
    !     allocate(ret(arrsize(1), arrsize(2)))
    !     do i = 1, arrsize(1)
    !         if (mask(i) .eqv. .true.) n_true = n_true + 1
    !     end do
    !     retsize = shape(reshape(pack(arr, mask), (/n_true, datasize(2)/)))
    !     allocate(retdata(retsize(1), retsize(2))
    !     retdata = reshape(pack(arr, mask_multi), retsize)

    ! end subroutine

    function dominated(idv1, idv2) result(res)
        implicit none 
        type(individual), intent(in) :: idv1, idv2
        type(individual) :: res

        if (all(idv1%values == idv2%values)) then 
            return
        end if
        
        if (all(idv1%values <= idv2%values)) then
            res = idv1 
        else 
            res = idv2 
        end if
        
    end function

    subroutine pick_pareto_solution(indivs, objdim, valdim,  paretoidv)
        implicit none 
        type(individual), intent(in) :: indivs(:)
        integer, intent(in) :: objdim, valdim
        ! real(8), allocatable, intent(out) :: pareto(:,:)
        type(individual), allocatable, intent(out) :: paretoidv(:)
        
        integer :: i, j, k
        integer :: n_thread, sum_thread
        integer :: datasize(2), retsize(2)
        logical, allocatable :: mask(:,:), submask(:)
        logical, allocatable :: indivmask(:)

        datasize(1) = size(indivs)
        datasize(2) = valdim
        print *, "indiv value size", datasize
        allocate(indivmask(datasize(1)))
        n_pareto = 0

        
        indivmask = .true.
        sum_thread = dble(omp_get_num_threads())
        !$omp parallel
        !$omp do
        do i = 1, datasize(1)
            n_thread = dble(omp_get_thread_num())
            if((mod(i,1000) == 0)) then 
                print "(i7' /'2i7)", i, datasize(1) , int(n_thread)
                ! print "(i7' /'2i7)", &
                ! print "(f10.5' /'2i7)", &
                ! & (i-(datasize(1)/dble(sum_thread))) / (datasize(1)/sum_thread), &
                ! & datasize(1), omp_get_thread_num()
            end if
            do j = 1, datasize(1)
                !if (i==j) then 
                !    continue
                !else if (indivmask(j) .eqv. .false.) then 
                !    continue
                !end if

                !if(any(indivs(i)%values(1:2) <= indivs(j)%values(1:2)) .eqv. .true.) then
                !    indivmask(i) = .true.
                !    
                !    !indivmask(j) = .true.
                if(all(indivs(i)%values(1:objdim) > indivs(j)%values(1:objdim))) then
                    indivmask(i) = .false.
                    exit
                end if
            end do
        end do
        !$omp end do

        ! !$omp do
        ! do i = 1, datasize(1)
        !     do j = 1, datasize(1)
        !         if(all(indivs(i)%values == indivs(j)%values)) then
        !             indivmask(j) = .false.
        !         end if
        !     end do
        ! end do
        ! !$omp end do
        !$omp end parallel


        !allocate(paretoidv(count(indivmask)))
        paretoidv = pack(indivs, indivmask)
        print *, count(indivmask)

    end subroutine pick_pareto_solution

    
    recursive subroutine qsort_indiv(indivlist, first, last, i_dim)
        implicit none
        type(individual), intent(inout):: indivlist(:)
        integer, intent(in):: first, last, i_dim

        type(individual) :: x, t
        integer :: i, j 

        x = indivlist((first+last)/2)
        i = first 
        j = last
        do 
            do while(indivlist(i)%values(i_dim) < x%values(i_dim))
                i = i + 1
            end do
            do while(x%values(i_dim) < indivlist(j)%values(i_dim))
                j = j - 1
            end do
            if (i >= j) exit 
            t = indivlist(i)
            indivlist(i) = indivlist(j)
            indivlist(j) = t 
            i = i + 1
            j = j - 1
        end do
        if (first < i-1) call qsort_indiv(indivlist, first, i-1, i_dim)
        if (j+1 < last) call qsort_indiv(indivlist, j+1, last, i_dim)

        return 

    end subroutine qsort_indiv 

    subroutine de_duplication(indivlist)
        implicit none 
        type(individual), intent(in) :: indivlist(:)
        type(individual), allocatable :: tmplist(:), tmplist2(:)
        integer :: i, j, k

        allocate(tmplist(2))
        tmplist(1) = indivlist(1)
        tmplist(2) = indivlist(2)

        do i = 3, size(indivlist)
            do j = 1, size(tmplist)
                if (.not. (all(indivlist(i)%values == tmplist(j)%values))) then
                    tmplist2 = tmplist 
                    deallocate(tmplist)
                    allocate(tmplist(size(tmplist2)+1))
                    do k = 1, size(tmplist2)
                        tmplist(k) = tmplist2(k)
                    end do
                    tmplist(size(tmplist)) = indivlist(i)
                    exit
                end if
            end do
        end do
    end subroutine de_duplication

end module pick_pareto


program main
    use readcsv, only: csvdata, read_csv
    use pick_pareto
    implicit none 
    real(8), allocatable :: retdata(:,:)
    type(individual), allocatable :: paretoidv(:)
    integer :: i
    integer :: objdim = 2
    character(64) :: formats

    ! call read_csv("const_opt_result.csv")
    call read_csv("duped.csv")

    write(*, "(a)", advance="no") "input n_object ->"
    read(*, *) objdim 
    write(*,"(i2)") objdim
    !print *, csvdata(size(csvdata, 1),:)

    call pick_feasible_solution(csvdata(:,2:), objdim, retdata)

    ! call de_duplication(feasible_indiv)
    call pick_pareto_solution(feasible_indiv, objdim, size(feasible_indiv(1)%values), paretoidv)
    do i = 1, objdim
        call qsort_indiv(paretoidv, 1, size(paretoidv), i)
    end do

    print *, "feasible indivs shape =", size(feasible_indiv), size(feasible_indiv(1)%values)
    print *, "pareto size", size(paretoidv)

    open(12, file="pareto.txt")
    write(12, "(a)") "#"

    write(formats, "(i2)") size(paretoidv(1)%values)-1
    write(formats, "(i2)") size(paretoidv(1)%values(1:objdim))-1
    formats = trim(formats)//"(f12.7' ')"
    formats = "("//trim(formats)
    formats = trim(formats)//",f12.7)"
    !print *, formats

    write(12, formats) (paretoidv(i)%values(1:objdim), i=1,size(paretoidv))
    ! do i = 1, size(paretoidv)
    !     !write(12, "(9(f12.7','),f12.7)") paretoidv(i)%values
    !     write(12, formats) paretoidv(i)%values(1:objdim)
    ! end do
    write(12, "(a)") "#"
    !print *, (feasible_indiv(1)%values(1:2)),(feasible_indiv(110000)%values(1:2))
    !print *, (feasible_indiv(1)%values(1:2) > feasible_indiv(110000)%values(1:2))
    print *, all(paretoidv(1)%values == paretoidv(2)%values)
end program
