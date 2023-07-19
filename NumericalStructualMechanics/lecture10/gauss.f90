! ---------------------------------------------------------------------
!   Gauss elimination
!   LU分解とガウスの消去法を行う．
!
!   Directory Structure
!   - lecture10 (curdir)
!       - gauss.f90
!       - a.out
!       - rst.txt
! ---------------------------------------------------------------------
program Gauss_elimination

    implicit none

    integer, parameter :: N = 4 ! 配列の確保のため

    integer :: i
    
    double precision :: A(N, N), X(N), B(N), AL(N, N), AD(N, N), AU(N, N)

    open(2, file='rst.txt')

    call inp(A, B, N)

    call gauss(A, X, B, AL, AD, AU, N)

    write(2, *) 'AL'
    do i = 1, N
        write(2, '(10E12.4)') AL(i, 1:N)
    end do

    write(2, *) 'AD'
    do i = 1, N
        write(2, '(10E12.4)') AD(i, 1:N)
    end do

    write(2, *) 'AU'
    do i = 1, N
        write(2, '(10E12.4)') AU(i, 1:N)
    end do

    write(2, *) 'A'
    do i = 1, N
        write(2, '(10E12.4)') A(i, 1:N)
    end do

    write(2, *) 'X'
    do i = 1, N
        write(2, '(E12.4)') X(i)
    end do

    close(2)

    stop
end

! ---------------------------------------------------------------------
!   analytical model input
! ---------------------------------------------------------------------
subroutine inp(a, b, N)
    implicit none
    integer ::  N
    double precision :: a(N, *), b(*)

    a(1,1) = 5D0
    a(2,1) = -4D0
    a(3,1) = 1D0
    a(4,1) = 0D0
    a(1,2) = a(2,1)
    a(2,2) = 6D0
    a(3,2) = -4D0
    a(4,2) = 1D0
    a(1,3) = a(3,1)
    a(2,3) = a(3,2)
    a(3,3) = 6D0
    a(4,3) = -4D0
    a(1,4) = a(4,1)
    a(2,4) = a(4,2)
    a(3,4) = a(4,3)
    a(4,4) = 5D0

    b(1) = 0D0
    b(2) = 1D0
    b(3) = 0D0
    b(4) = 0D0

    return 
end

! ---------------------------------------------------------------------
!   LU decomposition and Gaussian elimination
! ---------------------------------------------------------------------
subroutine gauss(a, x, b, al, ad, au, N)
    implicit none
    integer :: N, i, j, k
    double precision :: a(N, *), x(*), b(*), al(N, *), ad(N, *), au(N, *), sum_au, sum_al, sum_ad, y(N), temp, z(N)

    
    ! 下三角行列 al, 対角行列 ad, 上三角行列 au 成分を左上から右下に向かって逐次的に求めていく ---------------
    j = 1
    ad(j, j) = a(j, j)

    j = 2
    au(j-1, j) = a(j-1, j) / ad(j-1, j-1)
    al(j, j-1) = a(j, j-1) / ad(j-1, j-1)
    ad(j, j) = a(j, j) - al(j, j-1) * ad(j-1, j-1) * au(j-1, j)

    do j = 3, N ! n: 行列[A]の次元
        au(1, j) = a(1, j) / ad(1, 1)
        al(j, 1) = a(j, 1) / ad(1, 1)

        do i = 2, j-1
            sum_au = 0.D0 ! 初期化
            do k = 1, i-1
                sum_au = sum_au + al(i, k) * ad(k, k) * au(k, j)
            end do
            
            sum_al = 0.D0 ! 初期化
            do k = 1, i-1
                sum_al = sum_al + al(j, k) * ad(k, k) * au(k, i)
            end do
            
            au(i, j) = (a(i, j) - sum_au) / ad(i, i)
            al(j, i) = (a(j, i) - sum_al) / ad(i, i)
        end do
        sum_ad = 0.D0 ! 初期化
        do k = 1, j-1
            sum_ad = sum_ad + al(j, k) * ad(k, k) * au(k, j)
        end do
        ad(j, j) = a(j, j) - sum_ad
    end do

    ! 未知ベクトルを求めていく ---------------
    y(1:N) = 0D0 ! 初期化
    ! 前進代入
    y(1) = b(1) ! 1行目
    do i = 2, N
        temp = 0.D0 ! 初期化
        do j = 1, i-1
            temp = temp + al(i, j) * y(j)
        end do  
        y(i) = b(i) - temp
    end do
    
    ! 対角項
    do i = 1, N
        z(i) = y(i) / ad(i, i)
    end do

    ! 後退消去
    x(N) = z(N) ! n行目
    do j = N, 2, -1 ! 逆順のインクリメント．3つ目がstep size
        do k = 1, j-1
            z(k) = z(k) - au(k, j) * x(j)
        end do
        x(j-1) = z(j-1)
    end do

    return
end










