Module CMD_Progress
     Implicit None
  private
	  Logical , parameter , public :: CMD_PROGRESS_ABSOLUTE = .true.
  Type , public :: CLS_CMD_Progress
    Integer , private :: N , lens , i
    Character :: M = "*" , O = "."
    Character(len=64) :: Prefix
  Contains
    Procedure :: Set
    Procedure :: Put
  End Type CLS_CMD_Progress
   
contains
 
  Subroutine Set( this , N , L )
    Class( CLS_CMD_Progress ) :: this
    Integer , Intent( IN ) :: N , L
    this % N    = N
    this % lens = L
    this % i = 0
    this % Prefix = " Progress: " !//
  End Subroutine Set
   
  Subroutine Put( this , K , bAbsol )
    Class( CLS_CMD_Progress ) :: this
    Integer , Intent( IN ) :: K
    Logical , optional :: bAbsol
    Character(len=1) :: br
    integer :: jm
    this % i = this % i + K
    if ( present( bAbsol ) ) then
      if ( bAbsol ) this % i = K
    end if
    if ( this % i > this % n ) this % i = this % n   
    jm = Nint( real( this%i * this%lens ) / real( this%N ) )
    if ( this%i < this%n ) then
      br = char(13)
    else
      br = char(10)
    end if
    !write( * , '(5a,f6.2,2a)',advance="no") trim(this%Prefix) , '[' , &
    write( * , '(5a,f6.2,2a\)') trim(this%Prefix) , '[' , & !// 如您的编译器不支持，请用上方语句代替
      repeat(this%M , jm ) , repeat( this%O , this%lens-jm ) , '] ' , this%i*100.0/this%N , "%" , br
  End Subroutine Put
   
End Module CMD_Progress
 
Program www_fcode_cn
  use CMD_Progress
  Implicit None
  type( CLS_CMD_Progress ) ::Progress
  integer :: i , j
  call Progress % Set( N = 1700 , L = 25 )!// 1700次，显示长度25
  Progress % Prefix = "Test Progress:  "  !// 前方提示文字，不是必须
  Progress % M = "#" !// 已完成部分的字符，不是必须
  Progress % O = "." !// 未完成部分的字符，不是必须
  write(*,*) 'Begin'
  Do i = 0 , 1700 , 50
    call Progress % Put( i , CMD_PROGRESS_ABSOLUTE ) !// 绝对方式
    !call Progress % Put( 50 ) !// 也可用相对方式
    call sleep(1) !// 如您的编译器不支持，请用下方的循环代替
    !Do j = 1 , 10000000
    !End Do
  End Do
  write(*,*) 'End'
End Program www_fcode_cn
