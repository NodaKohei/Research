module parameter
    implicit none

    ! 系のパラメータ
    integer, parameter :: NN = 108 ! 粒子数
    double precision, parameter :: length = 4.762d0 ! 箱のサイズ
    double precision, parameter :: T = 1.0d0 ! 温度

    ! ちっこいパラメータ
    double precision, parameter :: eta = 0.09d0 ! displacement coefficient


    ! 入出力に関するパラメータ

    

end module parameter

program LJ_mc
    use parameter
    use mt19937
    implicit none

    ! 宣言
    integer :: i

    ! 粒子（構造体）
    type particle
        double precision :: x,y,z
    end type particle

    type(particle) :: ptcl(NN)
    type(particle) :: pre_ptcl(NN)




    ! メインの部分


    ! 初期配置
    call sgrnd(1)

    do i = 1, NN
        write(*,*) grnd()
        ptcl(i)%x = grnd()*2
        ptcl(i)%y = grnd()*2
        ptcl(i)%z = grnd()*2
    end do
        
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initial()
    use parameter
    implicit none


end program LJ_mc

