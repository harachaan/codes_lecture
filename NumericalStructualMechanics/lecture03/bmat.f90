! ----------------------------------------------------------------------
!   Finite Element Analysis
!   Element type: 3-node triangle elements (2-dimensional)
!
!   Directory Structure
!   - inputFile
!       - inp_concentrated-f_triangle.txt
!       - param.inc
!   - lecture03 (curdir)
!       - bmat.f90
!       
! ----------------------------------------------------------------------

program FEM_triangle_element    
    implicit none

    integer NNODE, NDOF, NELEM, MAXNODE, MAXELEM, MAXBC, lelem
    double precision THICK
    include '../inputFile/param.inc' ! MAXNODE, MAXELEM, MAXBCをinclude
    double precision X(MAXNODE), Y(MAXNODE), B(3, 6)
    integer LNODS(3, MAXELEM) ! 配列の領域確保？

    open(2, file='rst.txt')
    open(3, file='../inputFile/inp_concentrated-f_triangle.txt')
    open(4, file='bmat.txt')

    call inp(NNODE, NDOF, NELEM, X, Y, THICK, LNODS)
    close(3)

    ! do lelem = 1, NELEM
    do lelem = 1, 5
        call bmat(lelem, B, X, Y, LNODS)
    end do

    close(2)
    close(4)

    stop
end

! ----------------------------------------------------------------------
!   analytical model input
! ----------------------------------------------------------------------
subroutine inp(nnode, ndof, nelem, x, y, thick, lnods)
    implicit none
    
    integer nnode, ndof, nelem ! (節点数，総自由度数，要素数)
    integer inode, lelem ! for do文
    ! double precision :: thick = 1.D0 ! なんか初期値いれれんかった．．
    double precision thick
    integer lnods(3, *)
    double precision x(*), y(*) ! (節点座標(x, y), 要素を構成する節点（数字）のリスト)
    
    thick = 1.D0

    read(3, *) nnode, nelem ! inputファイルより接点数，要素数の読み込み

    do inode = 1, nnode
        read(3, *) x(inode), y(inode) ! 節点座標の読み込み
        write(2, *) x(inode), y(inode) ! 節点座標の書き出し
    end do

    ndof = 2 * nnode ! 総自由度数（未知変位成分の数）

    do lelem = 1, nelem
        read(3, *) lnods(1:3, lelem) ! 要素を構成する節点のリスト（必ず反時計回り）
    end do

    return 
end

! ----------------------------------------------------------------------
!   create B matrix (for each element)
! ----------------------------------------------------------------------
subroutine bmat(lelem, b, x, y, lnods)
    implicit none

    integer lelem, lnods(3, *) ! (対象となる要素の番号, 対象となる要素の節点3つ)
    integer inode, m ! for do文
    double precision x(*), y(*), b(3, *), xe(2, 3), Ax2

    ! 配列xeに対象要素の座標を入れている．
    do m = 1, 3
        inode = lnods(m, lelem)
        xe(1, m) = x(inode) ! xe: shape=(2次元平面＊要素の頂点の数)
        xe(2, m) = y(inode)
    end do

    Ax2 = xe(1, 2) * xe(2, 3) + xe(1, 1) * xe(2, 2) + xe(1, 3) * xe(2, 1) - xe(1, 2) * xe(2, 2) - xe(1, 3) * xe(2, 2) - xe(1, 1) * xe(2, 3)! 式(8)の計算

    b(1:3, 1:6) = 0.D0 ! 初期化
    ! 脳筋代入
    b(1, 1) = (xe(2, 2) - xe(2, 3)) / Ax2
    b(1, 3) = (xe(2, 3) - xe(2, 1)) / Ax2
    b(1, 5) = (xe(2, 1) - xe(2, 2)) / Ax2
    b(2, 2) = (xe(1, 3) - xe(1, 2)) / Ax2
    b(2, 4) = (xe(1, 1) - xe(1, 3)) / Ax2
    b(2, 6) = (xe(1, 2) - xe(1, 1)) / Ax2
    b(3, 1) = (xe(1, 3) - xe(1, 2)) / Ax2
    b(3, 2) = (xe(2, 2) - xe(2, 3)) / Ax2
    b(3, 3) = (xe(1, 1) - xe(1, 3)) / Ax2
    b(3, 4) = (xe(2, 3) - xe(2, 1)) / Ax2
    b(3, 5) = (xe(1, 3) - xe(1, 1)) / Ax2
    b(3, 6) = (xe(2, 1) - xe(2, 2)) / Ax2
    ! b = b / Ax2

    do m = 1, 3
        write(4, 100) b(m, 1:6)
    end do 
    100 format(6(E12.4, ','))

    return
end