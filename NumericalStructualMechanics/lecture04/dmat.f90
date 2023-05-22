! ---------------------------------------------------------------------
!   Finite Element Analysis
!   Element type: 3-node triangle elements (2-dimensional)
!
!   Directory Structure
!   - inputFile
!       - inp_concentrated-f_quadrilateral.txt
!       - inp_concentrated-f_triangle.txt
!       - param.inc
!   - lecture04 (curdir)
!       - dmat.f90
!       - a.out
!       - rst.txt
!       - dmat.txt
!       - bmat.txt
! ---------------------------------------------------------------------

program FEM_triangle_element
    implicit none

    integer :: MAXNODE, MAXELEM, MAXBC
    include '../inputFile/param.inc' ! MAXNODE, MAXELEM, MAXBCをinclude

    integer :: NNODE, NDOF, NELEM, LNODS(3, MAXELEM), lelem
        ! NNODE: 節点数，NDOF: 総自由度数，NELEM: 要素数，LNODS: 要素を構成する接点のリスト？, lelem: 対象要素の番号
    double precision :: X(MAXNODE), Y(MAXNODE), B(3, 6), D(3, 3), E, VNU, THICK
        ! (X, Y): 節点座標, B: Bmat, D: Dmat, E: ヤング率, VNU: ポアソン比, THICK: 厚さ？


    open(2, file = 'rst.txt') ! 節点座標が読み込めてるか確認するための書き出し先ファイル
    open(3, file = '../inputFile/inp_concentrated-f_triangle.txt') ! 節点座標が記述してあるファイル
    open(4, file = 'bmat.txt') ! bmatを記述するファイル
    open(5, file = 'dmat.txt') ! dmatを記述するファイル

    call inp(NNODE, NDOF, NELEM, X, Y, THICK, LNODS, E, VNU) ! 座標系，パラメータのinput
    close(3)

    call dmat(D, E, VNU) ! dmat 作成

    do lelem = 1, NELEM
        call bmat(lelem, B, X, Y, LNODS) ! bmat 作成 
    end do 

    close(2)
    close(4)
    close(5)

    stop
end

! ---------------------------------------------------------------------
!   analytical model input
! ---------------------------------------------------------------------
subroutine inp(nnode, ndof, nelem, x, y, thick, lnods, e, vnu)
    implicit none

    integer :: nnode, ndof, nelem, lnods(3, *)
    double precision :: thick, x(*), y(*), e, vnu

    integer :: inode, lelem
        ! for do文
    

    ! 機械的特性
    e = 200D3! ヤング率, [MPa]=[N/mm^2], `D3`は10^3を表す
    vnu = 0.4D0 ! ポアソン比[-]

    thick = 1.D0

    read(3, *) nnode, nelem ! inputファイルから節点数，要素数の読み込み

    do inode = 1, nnode
        read(3, *) x(inode), y(inode) ! 各節点の座標の読み込み
        write(2, *) x(inode), y(inode) ! rst.txtに記述
    end do

    ndof = 2 * nnode ! 総自由度数（未知変位成分の数）

    do lelem = 1, nelem
        read(3, *) lnods(1:3, lelem) ! 要素を構成する節点のリスト（必ず反時計回り？）
    end do

    return
end

! ---------------------------------------------------------------------
!   create B matrix (for each element)
!       - B matrix：shape=(ひずみ成分の数, 要素を構成する節点数x自由度)
! ---------------------------------------------------------------------
subroutine bmat(lelem, b, x, y, lnods)
    implicit none

    integer :: lelem, lnods(3, *)
        ! lelem: 対象要素の番号，lnods: 対象要素の節点のリスト
    double precision :: x(*), y(*), b(3, *)
        ! (x, y): 各節点の座標のリスト, b: bmat

    integer :: inode, m
    double precision :: xe(2, 3), Ax2
        ! xe: 対象要素の節点座標を格納shape=(2次元平面, 要素の頂点の数), Ax2: 要素の面積の2倍


    
    ! 配列xeに対象要素の座標を格納
    do m = 1, 3
        inode = lnods(m, lelem)
        xe(1, m) = x(inode) ! 対象要素のx座標を格納
        xe(2, m) = y(inode) ! 対象要素のy座標を格納
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

    ! なんかあとからbmatをAx2で割ろうとするとエラー出る．なんで．．
    ! b(1, 1) = (xe(2, 2) - xe(2, 3))
    ! b(1, 3) = (xe(2, 3) - xe(2, 1))
    ! b(1, 5) = (xe(2, 1) - xe(2, 2))
    ! b(2, 2) = (xe(1, 3) - xe(1, 2))
    ! b(2, 4) = (xe(1, 1) - xe(1, 3))
    ! b(2, 6) = (xe(1, 2) - xe(1, 1))
    ! b(3, 1) = (xe(1, 3) - xe(1, 2))
    ! b(3, 2) = (xe(2, 2) - xe(2, 3))
    ! b(3, 3) = (xe(1, 1) - xe(1, 3))
    ! b(3, 4) = (xe(2, 3) - xe(2, 1))
    ! b(3, 5) = (xe(1, 3) - xe(1, 1))
    ! b(3, 6) = (xe(2, 1) - xe(2, 2))
    ! b = b / Ax2

    do m = 1, 3
        write(4, 100) b(m, 1:6) ! bmatをbmat.txtに記述
    end do
    100 format(6(E12.4, ','))

    return
end

! ---------------------------------------------------------------------
!   create D matrix (for each element)
!       - D matrix：shape=(3, 3), 応力-ひずみ間形式の行列
! ---------------------------------------------------------------------
subroutine dmat(d, e, vnu)
    implicit none

    integer :: m
    double precision :: d(3,3), e, vnu, dCoef

    ! 資料04の式(12)をそのまま
    dCoef = e / (1 - vnu**2)
    d(1:3, 1:3) = 0.D0 ! 初期化
    d(1, 1) = 1
    d(1, 2) = vnu
    d(2, 1) = vnu
    d(2, 2) = 1
    d(3, 3) = (1 - vnu) / 2
    d = dCoef * d

    do m = 1, 3
        write(5, 100) d(m, 1:3) ! dmatをdmat.txtに記述
    end do 
    100 format(6(E12.4, ','))

    return
end












