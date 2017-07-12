Program FDS2SLCF

implicit none

    interface
        subroutine parse(buffer,sb_toks,se_toks,n_toks)
            implicit none
            character(*), intent(inout) :: buffer
            integer, dimension(*), intent(out) :: sb_toks, se_toks
            integer, intent(out) :: n_toks
        end subroutine parse
    end interface

    character(255), parameter :: f2aversion='2.1.0'
    integer, parameter :: fb = selected_real_kind(6)
    integer, parameter :: file_dim = 500
    integer, parameter :: time_dim = 1010
    integer :: ierr, NMESHES, nm, noc, i, j, k, l, ix, iy, iz
    integer :: idum, ifile, nsam, nv, mv
    integer :: imin,imax,jmin,jmax,kmin,kmax
    integer :: i1, i2, j1, j2, k1, k2, t1, t2, i3, j3, k3  ! dimension
    integer :: i10, i20, j10, j20, k10, k20
    integer :: ncount, ior_input, npatch, ijbar, jkbar
    integer :: ii, nxp, nyp, nzp, n, i_sample,jj
    real(fb) :: xs, xf, ys, yf, zs, zf  ! user defined bounds
    real(fb) :: d1, d2, d3, d4, tbeg, tend
    character(256) :: auto_slice_label
    integer :: auto_slice_flag, n_auto_slices, iauto, iztemp
    integer, dimension(1) :: izmin1,ixmin1,iymin1
    integer :: izmin,ixmin,iymin,ixoff,iyoff,izoff
    integer, dimension(file_dim) :: auto_slice_lists,temp_slice_lists
    real(fb), dimension(file_dim) :: auto_slice_z,auto_slice_y,auto_slice_x
    real(fb) :: ax1, ay1, az1, az2
    real(fb) :: eps, zmin


    type mesh_type
        real(fb), pointer, dimension(:) :: x,y,z
        real(fb) :: d1,d2,d3,d4
        integer :: ibar,jbar,kbar,ierr,nxp,nyp,nzp
    end type mesh_type

    type (mesh_type), dimension(:), allocatable, target :: mesh
    type(mesh_type), pointer :: m

    real(fb), dimension(60) :: sum
    integer, allocatable, dimension(:) :: ix1,ix2,iy1,iy2,iz1,iz2
    integer, allocatable, dimension(:) :: auto_slices
    real(fb), allocatable, dimension(:,:,:,:) :: quantity
    real(fb), allocatable, dimension(:,:,:,:) :: q
    real(fb), allocatable, dimension(:,:,:) :: f
    real(fb), allocatable, dimension(:) :: time
    logical, allocatable, dimension(:,:,:) :: already_used
    logical :: new_plot3d=.true.
    integer ior_slcf, stationID
    character(1) :: ans,seq
    character(4) :: ext1, ext2
    character(4) :: choice
    character(30) :: unitjunk
    character(256) gridfile,qname,chid,qfile,junk,frmt,outfile,outfile1,outfile2,slcf_label_dummy
    character(256), dimension(file_dim) :: pl3d_file,slcf_file,bndf_file,slcf_text,bndf_text,slcf_unit,bndf_unit,slcf_label
    character(256), dimension(file_dim) :: slice_labels
    character(20), dimension(file_dim) :: bndf_type
    character(20) :: bndf_type_chosen
    integer :: nslice_labels, slice_exist, new_slice_seq
    real(fb), dimension(file_dim) :: x1, x2, y1, y2, z1, z2
    integer,  dimension(60) :: ib,is
    integer,  dimension(file_dim) :: slcf_mesh
    logical, dimension(file_dim) :: file_exists
    integer :: rcode
    integer :: nfiles_exist
    real(fb) :: zoffset
    integer :: zoffset_flag
    logical :: exists
    integer :: lu_in, nargs, iarg, lenstring 
    character(256) :: buffer, filein
    integer, dimension(256) :: sb_toks, se_toks
    integer :: n_toks
    integer :: Error_status
    character(256) :: arg

    zoffset=0.0
    zoffset_flag=0
    eps=0.001
    filein='stdin'

    ! parse command line arguments

    if (filein.eq.'stdin') then
        lu_in=5
    else
        lu_in=7
        inquire(file=filein,exist=exists)
        if (.not.exists) then
            write(6,*)"*** fatal Error: the file: ",trim(filein), " does not exist"
            stop
        endif
        open(lu_in,file=filein,status='old',form='formatted')
    endif

    write(6,*) ' Enter job id string (chid):'

    read(lu_in,'(a)') chid

    ! check to see if the .smv file exists

    gridfile = trim(chid)//'.smv'
    inquire(file=gridfile,exist=exists)
    if (.not.exists) then
        write(6,*)"*** fatal Error: the file: ",trim(gridfile)," does not exist"
        stop
    endif

    ! open the .smv file

    open(11,file=gridfile,status='old',form='formatted')
    write(*,*) "Open SMV file: ",trim(gridfile)
    print *


    ! determine the number of meshes

    rewind(11)

    call search('NMESHES',7,11,ierr)
    if (ierr.eq.1) then
        write(6,*) ' warning: assuming 1 mesh'
        NMESHES = 1
    else
        read(11,*) NMESHES
        write(*,*) "Number of mesh is: ",NMESHES
        print *
    endif

    allocate(mesh(NMESHES))

    ! get the coordinates of the meshes

    rewind(11)
    read_smv: do nm=1,NMESHES

    m=>mesh(nm)
    call search('GRID',4,11,ierr)
    read(11,*) m%ibar,m%jbar,m%kbar
    allocate(m%x(0:m%ibar))
    allocate(m%y(0:m%jbar))
    allocate(m%z(0:m%kbar))

    call search('TRNX',4,11,ierr)
    read(11,*) noc
    do i=1,noc
    read(11,*)
    enddo
    do i=0,m%ibar
    read(11,*) idum,m%x(i)
    enddo

    call search('TRNY',4,11,ierr)
    read(11,*) noc
    do i=1,noc
    read(11,*)
    enddo
    do j=0,m%jbar
    read(11,*) idum,m%y(j)
    enddo

    call search('TRNZ',4,11,ierr)
    read(11,*) noc
    do i=1,noc
    read(11,*)
    enddo
    do k=0,m%kbar
    read(11,*) idum,m%z(k)
    enddo

    enddo read_smv
    
    write(6,*) ' domain selection:'
    write(6,*) '   y - domain size is limited'
    write(6,*) '   n - domain size is not limited'
    write(6,*) '   z - domain size is not limited and z levels are offset'
    read(lu_in,'(a)') ans
    print *
    call toupper(ans,ans)
    if (ans(1:1).eq.'y') then
        write(6,*) ' Enter min/max x, y and z'
        read(lu_in,*) xs,xf,ys,yf,zs,zf
    else
        xs = -100000.
        xf =  100000.
        ys = -100000.
        yf =  100000.
        zs = -100000.
        zf =  100000.
    endif

    ! if ans is z or z then subtract zoffset from z data (for multi-level cases)

    if (ans(1:1).eq.'z') then
        zoffset_flag=1
    else
        zoffset_flag=0
    endif

    
    !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    !
    !  Key ieteration search slcf (i), list slcf file and variables
    !		1. Print all exist slice file and their contains.
    !  		2. Generate unique slice variable list.
    !
    !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    
    ! Set some initial values
    slcf_label = 'null'
    slcf_mesh = 1
    rewind(11)
    nfiles_exist=0  ! record variables cover multiple slice files
    nslice_labels=0
    
    Search_slcf: do i=1,file_dim    

        ! General check and read SLCF headers
        call search2('SLCF',4,'SLCC',4,11,ierr,choice)
        if (ierr.eq.1) exit search_slcf
        backspace(11)
        read(11,*) junk,slcf_mesh(i)  ! Slcf dimensions
        read(11,'(a)') slcf_file(i)  ! slcf file name
        read(11,'(a)') slcf_text(i)  ! slcf variable defination
        
        ! Check if slice file exist
        read(11,*)    ! skip dummy line
        read(11,'(a)') slcf_unit(i)
        open(12,file=trim(slcf_file(i)),form='unformatted',status='old', iostat=rcode)
        if (rcode.ne.0) then
            close(12)
            cycle ! jump to search_slcf for another run
        endif
        
        ! record of exist slice files
        nfiles_exist=nfiles_exist+1  
        
        ! TASK 1: Read and Print what individual slice file contains.
   
        read(12) unitjunk
        read(12) unitjunk
        read(12) unitjunk
        read(12) i1,i2,j1,j2,k1,k2
        close(12)
    
        nm=slcf_mesh(i)
        m=>mesh(nm)
        x1(i)=m%x(i1)
        x2(i)=m%x(i2)
        y1(i)=m%y(j1)
        y2(i)=m%y(j2)
        z1(i)=m%z(k1)
        z2(i)=m%z(k2)
    
        write(6,'(i3,1x,a,1x,a)')i,trim(slcf_text(i)),trim(slcf_file(i))
        write(6,'(3x,a,6(1x,f8.2))')'slice bounds:',m%x(i1),m%x(i2),m%y(j1),m%y(j2),m%z(k1),m%z(k2)
        write(6,*) 
        
        ! TASK 2: Create unique list of slice types      
    
        slice_exist=0
        
        do ii=1, nslice_labels ! loop over labels kind
        if(trim(slcf_text(i)).eq.slice_labels(ii))then
            slice_exist=1 ! set flag slice already recorded
            exit  
        endif
        enddo
    
        if (slice_exist.eq.0) then ! new slice
            nslice_labels=nslice_labels+1  ! add record 
            slice_labels(nslice_labels)=trim(slcf_text(i)) ! Record new slice variable
        endif

    Enddo Search_slcf
    
    ! Error code: slice file not found
    write(*,*) "nfiles_exist =",nfiles_exist
    print *
    IF (nfiles_exist.EQ.0)THEN
        WRITE(6,*)"There are no slice files to convert"
        STOP
    ENDIF
    
    close(11)   ! close .svm file
    
    !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    !
    !  1. Print unique slice variable list 
    !  2. Select their stationary dimension 
    !  3. Adjust the mesh sequence to insure dimmension increases in all direction
    !
    !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    
    ! Task 1: Print the unique slice variable list.
    N_AUTO_SLICES=0
    write(6,*)' Enter slice file type index'
    DO II=1, nslice_labels
       write(6,'(2x,I3,1x,A)')II,TRIM(slice_labels(II))
    ENDDO
    READ(LU_IN,'(I3)')II
    AUTO_SLICE_LABEL=TRIM(slice_labels(II))

    ! Task 2: Select stationary dimension.
    print *
    write(*,*) "Choose the 2D slice plate (1/2/3):"
    write(*,'(3X,A)')  "1   Y-Z Plate"
    write(*,'(3X,A)')  "2   X-Z Plate"
    write(*,'(3X,A)')  "3   X-Y Plate"
    read(lu_in,*) StationID
    print *
    write(*,*) "The original Mesh sequence is: "
    
    n_auto_slices = 0
    
    Select case(stationID)
    
        Case(1) ! Y-Z plate
        boundloop1: DO I=1,nfiles_exist ! loop over slice files, equal to (I) in Search_slcf
           IF(TRIM(AUTO_SLICE_LABEL).NE.TRIM(SLCF_TEXT(I))) CYCLE boundloop1
           IF(x2(I)-x1(I).GT.EPS)CYCLE boundloop1
           n_auto_slices = n_auto_slices + 1
           nm=slcf_mesh(I)
           m=>mesh(nm)
           write(6,'(A,I3,5(a,f8.3))')"Mesh:",nm,", ymin=",y1(I),", ymax=",y2(I),", zmin=",z1(I),", zmax=",z2(I),", x=",x1(I)
           write(6,'(3(3X,A,I3))') "The mesh size on each axis is: NX =",m%ibar,", NY=",m%jbar,", NZ=",m%kbar
           auto_slice_lists(n_auto_slices) = I
        ENDDO boundloop1
        
        case(2) ! X-Z plate
        boundloop2: DO I=1,nfiles_exist ! loop over slice files, equal to (I) in Search_slcf
           IF(TRIM(AUTO_SLICE_LABEL).NE.TRIM(SLCF_TEXT(I))) CYCLE boundloop2
           IF(y2(I)-y1(I).GT.EPS)CYCLE boundloop2
           n_auto_slices = n_auto_slices + 1
           nm=slcf_mesh(I)
           m=>mesh(nm)
           write(6,'(A,I3,5(a,f8.3))')"Mesh:",nm,", xmin=",x1(I),", xmax=",x2(I),", zmin=",z1(I),", zmax=",z2(I),", y=",y1(I)
           write(6,'(3(3X,A,I3))') "The mesh size on each axis is: IX =",m%ibar,", IY=",m%jbar,", IZ=",m%kbar
           auto_slice_lists(n_auto_slices) = I
        ENDDO boundloop2
        
        case(3) ! X-Y plate
        boundloop3: DO I=1,nfiles_exist ! loop over slice files, equal to (I) in Search_slcf
           IF(TRIM(AUTO_SLICE_LABEL).NE.TRIM(SLCF_TEXT(I))) CYCLE boundloop3
           IF(z2(I)-z1(I).GT.EPS)CYCLE boundloop3
           n_auto_slices = n_auto_slices + 1
           nm=slcf_mesh(I)
           m=>mesh(nm)
           write(6,'(A,I3,5(a,f8.3))')"Mesh:",nm,", xmin=",x1(I),", xmax=",x2(I),", ymin=",y1(I),", ymax=",y2(I),", z=",z1(I)
           write(6,'(3(3X,A,I3))') "The mesh size on each axis is: IX =",m%ibar,", IY=",m%jbar,", IZ=",m%kbar
           auto_slice_lists(n_auto_slices) = I
        ENDDO boundloop3
        
    end select
    
    ! Task 3:  Adjust mesh list sequence
    
    allocate(ix1(n_auto_slices))
    allocate(ix2(n_auto_slices))
    allocate(iy1(n_auto_slices))
    allocate(iy2(n_auto_slices))
    allocate(iz1(n_auto_slices))
    allocate(iz2(n_auto_slices))
    
	print *	
    write(*,*) "Need to adjust the coordinates offset? (y/n)"
    read(lu_in,*) seq
    
    if (seq .eq. 'y' .or. seq .eq. 'Y') then
		write(*,*) "Provide approprate coordinates offset (xoff,yoff,zoff)"
		read(lu_in,*) ixoff,iyoff,izoff
		print *
		write(*,*) "New mesh coordinates after adjustment"
		
		! Select case
		
		select case(stationID)
		
		case(1)
		do J=1,n_auto_slices
           I=auto_slice_lists(J)
			ix1(J)=int(x1(I))
			ix2(J)=int(x2(I))
			iy1(J)=int(y1(I))-iyoff
			iy2(J)=int(y2(I))-iyoff
			iz1(J)=int(z1(I))-izoff
			iz2(J)=int(z2(I))-izoff
           write(6,'(A,I3,A,I3,5(a,f8.3))')"Seq:",J,", Mesh:",slcf_mesh(I),", ymin=",y1(I),", ymax=",y2(I),", zmin=",z1(I),", zmax=",z2(I),", x=",x1(I)
           write(6,'(3(3X,A,I3,2X,I3))') "The mesh size on each axis is: IX =",ix1(j),ix2(j),", IY=",iy1(j),iy2(j),", IZ=",iz1(j),iz2(J)
           if (J .eq. n_auto_slices) then
           allocate(quantity(0:ix2(J),0:iy2(J),0:iz2(J),time_dim))
           endif
        enddo
        
		case(2)
		do J=1,n_auto_slices
           I=auto_slice_lists(J)
			ix1(J)=int(x1(I))-ixoff
			ix2(J)=int(x2(I))-ixoff
			iy1(J)=int(y1(I))
			iy2(J)=int(y2(I))
			iz1(J)=int(z1(I))-izoff
			iz2(J)=int(z2(I))-izoff
           write(6,'(A,I3,A,I3,5(a,f8.3))')"Seq:",J,", Mesh:",slcf_mesh(I),", xmin=",x1(I),", xmax=",x2(I),", zmin=",z1(I),", zmax=",z2(I),", y=",y1(I)
           write(6,'(3(3X,A,I3,2X,I3))') "The mesh size on each axis is: IX =",ix1(j),ix2(j),", IY=",iy1(j),iy2(j),", IZ=",iz1(j),iz2(J)
           if (J .eq. n_auto_slices) then
           allocate(quantity(0:ix2(J),0:iy2(J),0:iz2(J),time_dim))
           endif
        enddo
        
		case(3)
		do J=1,n_auto_slices
           I=auto_slice_lists(J)
			ix1(J)=int(x1(I))-ixoff
			ix2(J)=int(x2(I))-ixoff
			iy1(J)=int(y1(I))-iyoff
			iy2(J)=int(y2(I))-iyoff
			iz1(J)=int(z1(I))
			iz2(J)=int(z2(I))
           write(6,'(A,I3,A,I3,5(a,f8.3))')"Seq:",J,", Mesh:",slcf_mesh(I),", xmin=",x1(I),", xmax=",x2(I),", ymin=",y1(I),", ymax=",y2(I),", z=",z1(I)
           write(6,'(3(3X,A,I3,2X,I3))') "The mesh size on each axis is: IX =",ix1(j),ix2(j),", IY=",iy1(j),iy2(j),", IZ=",iz1(j),iz2(J)
           if (J .eq. n_auto_slices) then
           allocate(quantity(0:ix2(J),0:iy2(J),0:iz2(J),time_dim))
           endif
        enddo
        
        end select
        else
        ixoff=0
        iyoff=0
        izoff=0
    endif
		
		
    !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    !
	! LOOP:  Read through N_Auto_Slices and write to output file
	!   	1. Calculate the total dimension for output
	!		2. Read and write file to output
	!
    !#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    
    ! Task 1: Input Dimension of IX,IY,IZ
    
    
        allocate(time(time_dim))
    
	Auto_List: DO JJ=1, n_auto_slices
	
    	I=AUTO_SLICE_LISTS(JJ) ! I indicate file ID in slcf_file
    	
    	nm=slcf_mesh(I)
        m=>mesh(nm)
!        allocate(q(0:m%ibar,0:m%jbar,0:m%kbar,nv))
        allocate(f(0:m%ibar,0:m%jbar,0:m%kbar))
        if ( JJ .eq. 1) then
        endif
        f = 0.
!        q = 0.
    	
    	
        IS(MV) = I
        QFILE = SLCF_FILE(I)	
        
        OPEN(12,FILE=QFILE,FORM='UNFORMATTED',STATUS='OLD')
        READ(12)
        READ(12)
        READ(12)
        READ(12) I1,I2,J1,J2,K1,K2                     
        
        ! Read each slcf file with specific variable
        ncount = 1
        t1=ncount
    
        Read_loop: do
            read(12,end=99) time(ncount)
            read(12,end=99) (((f(i,j,k),i=i1,i2),j=j1,j2),k=k1,k2)
            t2=ncount ! indicating T_end
            print *
            write(*,*) " Successfully reading file :", trim(qfile)
            print *
            write(*,*) "quantity index:"
            write(*,'(3(I3,A,I3,2X))') ix1(JJ),":",ix2(JJ),iy1(JJ),":",iy2(JJ),iz1(JJ),":",iz2(JJ)
            print *
            write(*,*) "f index:"
            write(*,'(3(I3,A,I3,2X))') i1,":",i2,j1,":",j2,k1,":",k2
            print *
            quantity(ix1(JJ):ix2(JJ),iy1(JJ):iy2(JJ),iz1(JJ):iz2(JJ),ncount) = f(i1:i2,j1:j2,k1:k2)
            if (ncount .ge. time_dim) then
                write(*,"(A,I4,A)") " *** Fatal error: preset time_dim = ",time_dim," is smaller than needed."
                write(*,"(A,f7.2)") " *** Try enlarge the parameter time_dim. Current time is: ",time(ncount)
                write(*,"(A,I4)") " *** Current iteration reaches: ",ncount
                stop
            endif
            ncount = ncount + 1
            
            deallocate(f)
            
        enddo Read_loop 
        write(*,*) "End read_loop"
        
        99  close(12)
    
    
    enddo Auto_list
    
    
    
    
    
    
    
    
    
End Program FDS2SLCF



!#################################################################
! 
!                     Subroutine Part 
!
!#################################################################

! *********************** search *******************************

subroutine search(string,length,lu,ierr)

    implicit none
    character(*), intent(in) :: string
    integer, intent(out) :: ierr
    integer, intent(in) :: lu, length
    character(20) :: junk

    search_loop: do 
    read(lu,'(a)',end=10) junk
    if (junk(1:length).eq.string(1:length)) exit search_loop
    enddo search_loop

    ierr = 0
    return

    10 ierr = 1
    return

end subroutine search

subroutine search2(string,length,string2,length2,lu,ierr,choice)

    implicit none
    character(*), intent(in) :: string,string2
    character(*), intent(out) :: choice
    integer, intent(out) :: ierr
    integer, intent(in) :: lu, length, length2
    character(20) :: junk

    search_loop: do 
    read(lu,'(a)',end=10) junk
    if (junk(1:length).eq.string(1:length).or.junk(1:length2).eq.string2(1:length2)) then
        if (junk(1:length) .eq.string(1:length))   choice = junk(1:length)
        if (junk(1:length2).eq.string2(1:length2)) choice = junk(1:length2)
        exit search_loop
    endif
    enddo search_loop

    ierr = 0
    return

    10 ierr = 1
    return

end subroutine search2

! *********************** parse *******************************

subroutine parse(buffer,sb_toks,se_toks,n_toks)

    ! parse buffer into a series of tokens each separated by blanks

    implicit none
    character(*), intent(inout) :: buffer
    integer, dimension(*), intent(out) :: sb_toks, se_toks
    integer, intent(out) :: n_toks
    integer :: i, intok, inquote, lenbuf
    character(len=1) :: c

    n_toks=0
    lenbuf = len(trim(buffer))
    if(lenbuf.eq.0)return ! buffer is blank so there are no tokens
    intok=0
    inquote=0
    do i = 1, lenbuf
    if(intok.eq.0)then  
        if(buffer(i:i).ne.' ')then  ! beginning of a new token since previous char
            intok=1                   ! was not in a token and this one is
            n_toks=n_toks + 1
            sb_toks(n_toks)=i
            if(buffer(i:i).eq."'")inquote=1
        endif
    else
        if(inquote.eq.1)then
            if(buffer(i:i).eq."'")then
                se_toks(n_toks)=i
                intok=0
                inquote=0
            endif
        endif
        if(buffer(i:i).eq.' ')then
            se_toks(n_toks)=i-1       ! previous char was in a token, this one is not
            intok=0                   ! so previous char is end of token
        endif
    endif
    end do
    if(buffer(lenbuf:lenbuf).ne.' ')se_toks(n_toks)=lenbuf ! last char in buffer is not blank
    ! so it is end of last token

    ! strip out single or double quotes if present
    do i = 1, n_toks
    c = buffer(sb_toks(i):sb_toks(i))
    if(c.eq."'")sb_toks(i)=sb_toks(i)+1
    c = buffer(se_toks(i):se_toks(i))
    if(c.eq."'")se_toks(i)=se_toks(i)-1
    if(se_toks(i).lt.sb_toks(i))then
        se_toks(i)=sb_toks(i)
        buffer(sb_toks(i):se_toks(i))=' '
    endif
    end do
end subroutine parse

! *********************** toupper *******************************

subroutine toupper(bufferin, bufferout)
    character(len=*), intent(in) :: bufferin
    character(len=*), intent(out) :: bufferout
    character(len=1) :: c

    integer :: lenbuf, i

    lenbuf=min(len(trim(bufferin)),len(bufferout))
    do i = 1, lenbuf
    c = bufferin(i:i)
    if(c.ge.'a'.and.c.le.'z')c=char(ichar(c)+ichar('a')-ichar('a'))
    bufferout(i:i)=c
    end do

end subroutine toupper
