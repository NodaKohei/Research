! parameter を設定するモジュール
module parameter
    implicit none

    ! 系のパラメータ
    integer, parameter :: NN = 108 ! 粒子数
    double precision, parameter :: length = 4.762d0 ! 箱のサイズ
    double precision, parameter :: T = 1.0d0 ! 温度

    ! ちっこいパラメータ
    double precision, parameter :: eta = 0.09d0 ! displacement coefficient
    integer, parameter :: seed = 1 ! seed of random number
    integer, parameter :: MC_sweeps = 100 ! モンテカルロスウィープ数 １スウィープはNNとする


    ! 入出力に関するパラメータ(ファイル名)
    character(64), parameter :: filepath = "./"
    character(64), parameter :: filename = "sample.xyz"
    

end module parameter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program LJ_mc
    use parameter
    use mt19937
    implicit none

    !!! 宣言 !!!
    ! 変数たち
    integer :: i, j
    integer :: N_select
    double precision :: l, r, E_pre, E_new, E_pair

    ! 粒子（構造体）
    type particle
        double precision :: x,y,z
    end type particle

    type(particle) :: p(NN)
    type(particle) :: p_memo(NN)



    !!! プログラム !!!
    ! 初期配置
    call sgrnd(seed)
    call initial(p)

    call output(p, filepath, filename)

    ! 初期エネルギーの計算
    call energy(p, E_new)

    ! mainのループ
    do i = 1, MC_sweeps
        do j = 1, NN ! 1 sweepはNNとする
            ! エネルギーと座標をメモする
            E_pre = E_new
            p_memo = p

            ! 微小変化させる粒子を選択
            N_select = int( NN * grnd() ) + 1

            ! ちょっと動かす
            p(N_select)%x = p(N_select)%x + eta * ( 2.0d0 * grnd() - 1.0d0 ) 
            p(N_select)%y = p(N_select)%y + eta * ( 2.0d0 * grnd() - 1.0d0 ) 
            p(N_select)%z = p(N_select)%z + eta * ( 2.0d0 * grnd() - 1.0d0 ) 

            ! 周期境界条件の修正
            ! x
            l = int( ( 2.0d0 * p(N_select)%x - length ) / length)
            p(N_select)%x = p(N_select)%x - length * l
            ! y
            l = int( ( 2.0d0 * p(N_select)%y - length ) / length)
            p(N_select)%y = p(N_select)%y - length * l
            ! z
            l = int( ( 2.0d0 * p(N_select)%z - length ) / length)
            p(N_select)%z = p(N_select)%z - length * l

            ! エネルギー変化の計算

            write(*,*) (i-1)*NN+j, N_select
        end do
    end do
        
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 初期配置をランダムに作成するサブルーチン
subroutine initial(p_init)
    use parameter
    implicit none
    
    integer :: i
    type(particle), intent(out) :: p_init(NN)

    do i = 1, NN 
      p_init(i)%x = grnd()*length 
      p_init(i)%y = grnd()*length 
      p_init(i)%z = grnd()*length 
    end do
    
end subroutine initial

! 座標を書き出すサブルーチン(vmdの形式)
subroutine output(p, fpath, fname)
    use parameter
    implicit none

    type(particle), intent(in) :: p(NN)
    character(len=64), intent(in) :: fpath, fname
    integer :: i
    character(len=128) :: name

    ! 文字列を結合している TRIMでそれぞれの文字列の空白を消している
    name = TRIM(fpath)//TRIM(fname)
    open(50, file = name, status = "replace")
    write(50,*) NN 
    write(50,*)
    do i = 1, NN
        write(50,*) "Ar", p(i)
    end do
    close(50)
end subroutine output

! 周期境界条件を考慮して粒子間距離を計算するサブルーチン
subroutine distance(p, i, j, r)
    use parameter
    implicit none

    type(particle), intent(in) :: p(NN)
    integer, intent(in) :: i, j
    double precision, intent(out) :: r

    double precision :: x_ij, y_ij, z_ij, l

    ! x方向の計算
    x_ij = p(i)%x - p(j)%x
    l = int(2.0d0 * x_ij / length)
    x_ij = x_ij - l * length
    
    ! y方向の計算
    y_ij = p(i)%y - p(j)%y
    l = int(2.0d0 * y_ij / length)
    y_ij = y_ij - l * length

    ! z方向の計算
    z_ij = p(i)%z - p(j)%z
    l = int(2.0d0 * z_ij / length)
    z_ij = z_ij - l * length

    ! 距離計算
    r = sqrt(x_ij**2.0d0 + y_ij**2.0d0 + z_ij**2.0d0)
end subroutine distance

! 2体間のエネルギーを計算するサブルーチン(Hamiltonianの設定)
subroutine energy_pair(r, E)
    use parameter
    implicit none

    double precision, intent(in) :: r ! diatance
    double precision, intent(out) :: E  ! Energy

    E = 4.0d0 * (r**(-12.0d0) - r**(-6.0d0))

end subroutine energy_pair

! 全系のエネルギーを計算するサブルーチン
subroutine energy(p, E)
    use parameter
    implicit none

    integer :: i, j 
    double precision :: E_pair
    type(particle), intent(in) :: p(NN)
    double precision, intent(out) :: E
    
    ! 念のための初期化
    E = 0.0d0

    ! エネルギー計算
    do i = 1, NN 
        do j = i + 1, NN
            call distance(p, i, j, r)
            call energy_pair(r, E_pair)
            E = E + E_pair
        end do
    end do

    ! エネルギーを二倍にする
    E = 2 * E

end subroutine energy

! 1粒子動かしたことによるエネルギー変化を計算するサブルーチン
subroutine delta_energy(p, p_pre, number, dE)
    use parameter
    implicit none


end program LJ_mc

