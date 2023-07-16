! ---------------------------------------------------------------------
!   Finite Element Analysis
!   Element type: 3-node triangle elements (2-dimensional)
!   全体剛性マトリックスKを求める．
!
!   Directory Structure
!   - inputFile
!       - inp_concentrated-f_quadrilateral.txt
!       - inp_concentrated-f_triangle.txt
!       - param.inc
!   - lecture05 (curdir)
!       - Kmat.f90
!       - a.out
!       - rst.txt
!       - kmat.txt
!       - Kmat_diag.txt
!       - dmat.txt
!       - bmat.txt
! ---------------------------------------------------------------------

program FEM_triangle_element
    implicit none

    integer, parameter :: DIM = 2

    integer :: MAXNODE, MAXELEM, MAXBC
    include '../inputFile/param.inc' ! MAXNODE, MAXELEM, MAXBCをinclude

    integer :: NNODE, NDOF, NELEM, LNODS(3, MAXELEM), lelem, i
        ! NNODE: 節点数，NDOF: 総自由度数，NELEM: 要素数, LNODS: 要素を構成する接点のリスト, lelem: 対象要素の番号
    double precision :: X(MAXNODE), Y(MAXNODE), B(3, 6), D(3, 3), TK(2*MAXNODE, 2*MAXNODE), E, VNU, THICK
        ! (X, Y): 節点座標, B: ひずみ-変位mat, D: 応力-ひずみmat, TK: 全体剛性mat, E: ヤング率, VNU: ポアソン比, THICK: 厚さ
    

    open(2, file = 'rst.txt') ! 節点座標が読み込めてるか確認するための書き出し先ファイル
    open(3, file = '../inputFile/inp_concentrated-f_triangle.txt') ! 節点座標が記述してあるファイル
    open(4, file = 'bmat.txt') ! bmatを記述するファイル
    open(5, file = 'dmat.txt') ! dmatを記述するファイル
    open(6, file = 'kmat.txt') ! 要素剛性マトリックスを記述するファイル
    open(7, file = 'Kmat_diag.txt') ! 全体剛性マトリックスの体格高を記述するファイル

    call inp(NNODE, NDOF, NELEM, X, Y, THICK, LNODS, E, VNU, DIM) ! 座標系，パラメータのinput
    close(3)

    call dmat(D, E, VNU) ! dmat 作成

    TK(1:NDOF, 1:NDOF) = 0.D0 ! 全体剛性マトリックスを初期化

    do lelem = 1, nelem
        call bmat(lelem, B, X, Y, LNODS) ! bmat 作成
        call mergekmat(TK, lelem, LNODS, B, X, Y, D, THICK, MAXNODE)
    end do

    write(7, *) 'TK(i, i)'
    ! write(7, *) TK
    do i = 1, NDOF
        write(7, '(F12.4)') TK(i, i)
    end do

    close(2)
    close(4)
    close(5)
    close(6)
    close(7)

    stop
end

! ---------------------------------------------------------------------
!   analytical model input
! ---------------------------------------------------------------------
subroutine inp(nnode, ndof, nelem, x, y, thick, lnods, e, vnu, dim)
    implicit none

    integer :: nnode, ndof, nelem, lnods(3, *), dim
    double precision :: thick, x(*), y(*), e, vnu

    integer :: inode, lelem
        ! for do文


    ! 機械的特性
    e = 200D3 ! ヤング率, [MPa]=[N/mm^2], `D3`は10^3を表している
    vnu = 0.4D0 ! ポアソン比, [-]

    thick = 1.D0

    read(3, *) nnode, nelem ! inputファイルから節点数，要素数の読み込み

    do inode = 1, nnode
        read(3, *) x(inode), y(inode) ! 各節点の座標の読み込み
        write(2, *) x(inode), y(inode) ! rst.txtに記述
    end do

    ndof = dim * nnode ! 総自由度数（未知変位成分の数）


    do lelem = 1, nelem
        read(3, *) lnods(1:3, lelem) ! 要素を構成する節点のリスト（必ず反時計回り）
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
        ! lelem: 対象要素の番号, lnods: 対象要素の節点のリスト
    double precision :: x(*), y(*), b(3, *)
        ! (x, y): 各節点の座標のリスト, b: bmat

    integer :: inode, matrix, m
    double precision :: xe(2, 3), Ax2
        ! xe: 対象要素の節点座標を格納, shape=(2次元平面，要素の頂点の数), Ax2: 要素の面積の2倍

    
    ! 配列xeに対象要素の座標を格納
    do m = 1, 3
        inode = lnods(m, lelem)
        xe(1, m) = x(inode) ! 対象要素のx座標を格納
        xe(2, m) = y(inode) ! 対象要素のy座標を格納
    end do

    Ax2 = xe(1, 2) * xe(2, 3) + xe(1, 1) * xe(2, 2) + xe(1, 3) * xe(2, 1) - xe(1, 2) * xe(2, 1) - xe(1, 3) * xe(2, 2) - xe(1, 1) * xe(2, 3)! 式(8)の計算
    
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
    b(3, 5) = (xe(1, 2) - xe(1, 1)) / Ax2
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
!   create D matrix 
!       - D matrix：shape=(3, 3), 応力-ひずみ間形式の行列
! ---------------------------------------------------------------------
subroutine dmat(d, e, vnu)
    implicit none

    double precision :: d(3, 3), e, vnu
    
    integer :: m
        ! for do文
    double precision :: dCoef


    ! 資料04の式(12)をそのまま
    ! dCoef = e / (1 - vnu**2)
    ! d(1:3, 1:3) = 0.D0 ! 初期化
    ! d(1, 1) = 1
    ! d(1, 2) = vnu
    ! d(2, 1) = vnu
    ! d(2, 2) = 1
    ! d(3, 3) = (1 - vnu) / 2
    ! d = dCoef * d

    ! 資料04の式(6)をべたうち
    dCoef = e / (1 - 2 * vnu) / (1 + vnu)
    d(1:3, 1:3) = 0.D0 ! 初期化
    d(1, 1) = 1 - vnu
    d(1, 2) = vnu
    d(2, 1) = vnu
    d(2, 2) = 1 - vnu
    d(3, 3) = (1 - 2 * vnu) / 2
    d = dCoef * d

    do m = 1, 3
        write(5, 100) d(m, 1:3) ! dmatをdmat.txtに記述
    end do 
    100 format(6(E12.4, ','))

    return
end

! ---------------------------------------------------------------------
!   create K matrix (for each element)
!       - D matrix：shape=(3, 3), 応力-ひずみ間形式の行列
! ---------------------------------------------------------------------
subroutine mergekmat(tk, lelem, lnods, b, x, y, d, thick, maxnode)
    implicit none

    integer :: lelem, lnods(3, *), node_no, ip(6), maxnode
        ! lelem: 対象要素の要素番号, lnods: 対象要素の節点のリスト, node_no: 対象節点の番号, 
        ! ip(6): 要素剛性matの行・列と，全体剛性matに組み込む場所の対応関係, db(3, 6): [D][B]を一時記憶
    double precision :: thick, Ax2, Ae, x(*), y(*), xe(2, 3), b(3, 6), d(3, 3), db(3, 6), dk(6, 6), tk(2*maxnode, 2*maxnode)
        ! thick: 厚さ, Ae: 要素の面積, Ax2: Aeの2倍

    integer :: inode, idof, i, j, k, m
        ! for do文

    ! 要素剛性[dk]と全体剛性[tk]の対応関係
    ip(1:6) = 0.D0 ! 初期化
    do inode = 1, 3 ! 要素の3個の節点
        node_no = lnods(inode, lelem) ! 対象節点の番号
        do idof = 1, 2 ! 自由度（degree of freedom）
            ip((inode-1)*2 + idof) = (node_no-1)*2 + idof ! 要素の節点変位ベクトルでの順番 = 全体剛性での自由度番号
        end do
    end do

    ! 配列xeに対象要素の座標を格納
    do m = 1, 3
        inode = lnods(m, lelem)
        xe(1, m) = x(inode) ! 対象要素のx座標を格納
        xe(2, m) = y(inode) ! 対象要素のy座標を格納
    end do

    Ax2 = xe(1, 2) * xe(2, 3) + xe(1, 1) * xe(2, 2) + xe(1, 3) * xe(2, 1) - xe(1, 2) * xe(2, 1) - xe(1, 3) * xe(2, 2) - xe(1, 1) * xe(2, 3)  ! 式(8)の計算
    Ae = Ax2 / 2.D0 ! 要素の面積の2倍

    ! まず，[D(3, 3)][B(3, 6)]を計算し，[db(3, 6)]に記憶
    do i = 1, 3 ! dbの行
        do j = 1, 6 ! db の列
            db(i, j) = 0.D0 ! 初期化
            do k = 1, 3 ! the product of d's row and b's column
                db(i, j) = db(i, j) + d(i, k) * b(k, j)
            end do
        end do
    end do

    ! 次に，要素剛性mat, [dk(6, 6)] = [B(3, 6)]^T [db(3, 6)] * Ae * thick
    do i = 1, 6 ! dk's row
        do j = 1, 6 ! dk's column
            dk(i, j) = 0.D0 ! 初期化
            do k = 1, 3 ! the product of b^T's row and db's column
                dk(i, j) = dk(i, j) + b(k, i) * db(k, j) ! 転置に注意
            end do
            dk(i, j) = dk(i, j) * Ae * thick
        end do
    end do

    ! 要素番号1--5の要素剛性matを記述したい
    if (lelem .le. 5) then
        write(6, *) '要素番号', lelem, 'の要素剛性'
        write(6, 100) dk
        100 format(6(E12.4, ','))
    end if



    ! 要素剛性マトリックスを全体剛性マトリックスに組み込む（ミスりがちらしいので注意）
    do i = 1, 6 ! 要素剛性matの行
        do j = 1, 6 ! 要素剛性matの列
            tk(ip(i), ip(j)) = tk(ip(i), ip(j)) + dk(i, j)
        end do 
    end do 

    return
end
