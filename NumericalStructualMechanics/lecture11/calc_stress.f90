! ---------------------------------------------------------------------
!   Finite Element Analysis
!   Element type: 3-node triangle elements (2-dimensional)
!   節点変位ベクトルと各要素の応力を出力する．
!
!   Directory Structure
!   - inputFile
!       - inp_concentrated-f_quadrilateral.txt
!       - inp_concentrated-f_triangle.txt
!       - param.inc
!   - lecture07 (curdir)
!       - calc_stress.f90
!       - a.out
!       - rst.txt
!       - element_coordinates.txt
!       - kmat.txt
!       - Kmat_diag.txt
!       - dmat.txt
!       - bmat.txt
!       - bound.txt
! ---------------------------------------------------------------------

program FEM_triangle_element
    implicit none

    integer, parameter :: DIM = 2

    integer :: MAXNODE, MAXELEM, MAXBC
    include '../inputFile/param.inc' ! MAXNODE, MAXELEM, MAXBCをinclude

    integer :: NNODE, NDOF, NELEM, LNODS(3, MAXELEM), lelem, i, N_BC_GIVEN, inode
        ! NNODE: 節点数，NDOF: 総自由度数，NELEM: 要素数, LNODS: 要素を構成する接点のリスト, lelem: 対象要素の番号
    double precision :: X(MAXNODE), Y(MAXNODE), B(3, 6), D(3, 3), TK(2*MAXNODE, 2*MAXNODE), E, VNU, THICK, P(2*MAXNODE), I_BC_GIVEN(MAXBC), V_BC_GIVEN(MAXBC)
        ! (X, Y): 節点座標, B: ひずみ-変位mat, D: 応力-ひずみmat, TK: 全体剛性mat, E: ヤング率, VNU: ポアソン比, THICK: 厚さ
    double precision :: U(2*MAXNODE), STRAIN(3, MAXELEM), STRESS(3, MAXELEM), XE_C(2)
        ! U: 節点変位ベクトル, STRAIN: 要素のひずみ(epsilon_xx, epsilon_yy, gamma_xy), STRESS: 要素の応力(sigma_xx, sigma_yy, tau_xy)

    open(2, file = 'element_coordinates.txt') ! 節点座標が読み込めてるか確認するための書き出し先ファイル
    open(3, file = '../inputFile/inp_concentrated-f_triangle.txt') ! 節点座標が記述してあるファイル
    open(4, file = 'bmat.txt') ! bmatを記述するファイル
    open(5, file = 'dmat.txt') ! dmatを記述するファイル
    open(6, file = 'kmat.txt') ! 要素剛性マトリックスを記述するファイル
    open(7, file = 'Kmat_diag_P.txt') ! 全体剛性マトリックスの対角項と，B.C.を考慮した場合の右辺ベクトルPを記述するファイル
    open(8, file = 'rst.txt') ! 節点変位ベクトルと各要素の応力
    open(9, file = 'elem_yposi_ystress.csv') ! 中央の要素のy座標（重心位置），y方向垂直応力

    call inp(NNODE, NDOF, NELEM, X, Y, THICK, LNODS, E, VNU, DIM, P, N_BC_GIVEN, I_BC_GIVEN, V_BC_GIVEN) ! 座標系，パラメータのinput
    close(3)

    call dmat(D, E, VNU) ! dmat 作成

    TK(1:NDOF, 1:NDOF) = 0.D0 ! 全体剛性マトリックスを初期化

    do lelem = 1, nelem
        call bmat(lelem, B, X, Y, LNODS) ! bmat 作成
        call mergekmat(TK, lelem, LNODS, B, X, Y, D, THICK, MAXNODE)
    end do

    call bound(TK, P, NDOF, N_BC_GIVEN, I_BC_GIVEN, V_BC_GIVEN, MAXNODE)
    ! do i = 1, NDOF
    !     write(7, '(10F12.4)') TK(i, 1:NDOF)
    ! end do

    write(7, *) 'TK(i, i), P(i)'
    do i = 1, NDOF
        write(7, '(2F12.4)') TK(i, i), P(i)
    end do

    call gauss2(TK(1:NDOF, 1:NDOF), U, P, NDOF) ! 変位を求める
    ! write(8, *) TK
    write(8, *) 'U(i)'
    do inode = 1, NNODE
        write(8, '(I5, 2E15.7)') inode, U(2*inode-1), U(2*inode) ! 節点変位を記述
    end do

    do lelem = 1, NELEM
        call bmat(lelem, B, X, Y, LNODS)
        call calc_stress(lelem, STRAIN, STRESS, B, D, U, LNODS) ! ひずみと応力を計算
    end do

    write(8, *) 'STRESS'
    do lelem = 1, NELEM
        write(8, '(I5, 3E15.7)') lelem, STRESS(1:3, lelem) ! 各要素の応力を記述
    end do

    do lelem = 91, 100 ! 中央の要素をピックアップ
        call xe_c_posi(lelem, LNODS, X, Y, U, XE_C) ! 中央の要素のy座標（重心位置）をXE_Cに格納
        write(9, *) XE_C(2), ',', STRESS(2, lelem) ! y座標とy方向垂直応力を記述
    end do
    

    close(2)
    close(4)
    close(5)
    close(6)
    close(7)
    close(8)
    close(9)

    stop
end

! ---------------------------------------------------------------------
!   analytical model input
! ---------------------------------------------------------------------
subroutine inp(nnode, ndof, nelem, x, y, thick, lnods, e, vnu, dim, p, n_bc_given, i_bc_given, v_bc_given)
    implicit none

    integer :: nnode, ndof, nelem, lnods(3, *), dim, n_bc_given, i_bc_given(*), i_load_given, num, idof
    double precision :: thick, x(*), y(*), e, vnu, p(*), v_bc_given(*)

    integer :: inode, lelem
        ! for do文


    ! 機械的特性 ----------
    e = 200D3 ! ヤング率, [MPa]=[N/mm^2], `D3`は10^3を表している
    vnu = 0.4D0 ! ポアソン比, [-]

    thick = 1.D0

    ! inputファイルから読み込み ----------
    read(3, *) nnode, nelem ! 節点数，要素数の読み込み

    do inode = 1, nnode
        read(3, *) x(inode), y(inode) ! 各節点の座標の読み込み
        write(2, *) x(inode), y(inode) ! rst.txtに記述
    end do

    ndof = dim * nnode ! 総自由度数（未知変位成分の数）


    do lelem = 1, nelem
        read(3, *) lnods(1:3, lelem) ! 要素を構成する節点のリスト（必ず反時計回り）
    end do

    ! 境界条件の設定 ----------
    ! 変位境界条件 (x方向)
    num = 0 ! 境界条件の番号，インクリメントする
    idof = 1 ! x方向
    do inode = 1, nnode
        if((x(inode) .eq. 0D0) .and. (y(inode) .eq. 0D0)) then
            num = num + 1
            i_bc_given(num) = inode ! 節点番号
            v_bc_given(num) = 0.D0 ! 強制変位量
            i_bc_given(num) = (i_bc_given(num) - 1) * 2 + idof
        end if
    end do

    ! 変位境界条件（y方向)
    idof = 2 ! y方向
    do inode = 1, nnode
        if(y(inode) .eq. 0D0) then
            num = num + 1
            i_bc_given(num) = inode ! 節点番号
            v_bc_given(num) = 0.D0 ! 強制変位量
            i_bc_given(num) = (i_bc_given(num) - 1) * 2 + idof
        end if 
    end do

    n_bc_given = num

    ! 力学的境界条件
    p(1:ndof) = 0.D0 ! 初期化（力の作用しない自由表面では，節点力はゼロ）
    idof = 2 ! y方向
    do inode = 1, nnode
        if((x(inode) .eq. 0D0) .and. (y(inode) .eq. 50D0)) then
            i_load_given = (inode - 1) * 2 + idof ! 自由度番号
            p(i_load_given) = -10D0 * thick ! 集中荷重: 10N
        end if
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
    Ae = Ax2 / 2 ! 要素の面積の2倍

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

! ---------------------------------------------------------------------
!   consider Boundary Condition 
!       - 
! ---------------------------------------------------------------------
subroutine bound(tk, p, ndof, n_bc_given, i_bc_given, v_bc_given, maxnode)
    implicit none

    integer :: ndof, n_bc_given, i_bc_given(*), maxnode
        ! ndof: 総自由度数, n_bc_given: 変位境界条件の数, i_bc_given: 変位境界条件を設定する自由度, maxnode: 最大メモリ？
    double precision :: tk(2*maxnode, 2*maxnode), p(*), v_bc_given(*)
        ! tk: 全体剛性マトリックス, p: 境界条件考慮した時の右辺ベクトル, v_bc_given: 変位境界条件の既知変位の値

    integer :: i, j, K
        ! for do文，jは境界条件を設定する自由度番号を代入する

    ! 元の右辺ベクトルp(*)から系数行列のj列のbarU_j倍をっし引く
    do i = 1, n_bc_given ! i: 境界条件のインデックス
        j = i_bc_given(i) ! 境界条件(i)を設定する自由度番号
        do k = 1, ndof
            p(k) = p(k) - v_bc_given(j) * tk(k, j)
        end do
    end do

    ! 右辺ベクトルp(*)のj行に既知変位の値barU_jを代入する．（上書き）
    do i = 1, n_bc_given
        j = i_bc_given(i) ! 境界条件(i)を設定する自由度番号
        p(j) = v_bc_given(j)
    end do

    ! 系数行列のj行とj列をすべて0にする．対角項(i, j)を1にする
    do i = 1, n_bc_given
        j = i_bc_given(i) ! 境界条件(i)を設定する自由度番号
        do k = 1, ndof
            tk(k, j) = 0
            tk(j, k) = 0
        end do
        tk(j,j) = 1
    end do

end

! ---------------------------------------------------------------------
!   LU decomposition and Gaussian elimination
! ---------------------------------------------------------------------
subroutine gauss2(a, x, b, N)

    implicit none
    integer :: N, i, j, k
    double precision :: a(N, *), x(*), b(*), sum_au, sum_al, sum_ad, y(N), temp, z(N)

    
    ! 下三角行列 al, 対角行列 ad, 上三角行列 au 成分を左上から右下に向かって逐次的に求めていく ---------------
    j = 1
    ! ad(j, j) = a(j, j)

    j = 2
    a(j-1, j) = a(j-1, j) / a(j-1, j-1)
    a(j, j-1) = a(j, j-1) / a(j-1, j-1)
    a(j, j) = a(j, j) - a(j, j-1) * a(j-1, j-1) * a(j-1, j)

    do j = 3, N ! n: 行列[A]の次元
        a(1, j) = a(1, j) / a(1, 1)
        a(j, 1) = a(j, 1) / a(1, 1)

        do i = 2, j-1
            sum_au = 0.D0 ! 初期化
            do k = 1, i-1
                sum_au = sum_au + a(i, k) * a(k, k) * a(k, j)
                ! write(*, *) a(k,k)
            end do
            
            sum_al = 0.D0 ! 初期化
            do k = 1, i-1
                sum_al = sum_al + a(j, k) * a(k, k) * a(k, i)
            end do
            
            a(i, j) = (a(i, j) - sum_au) / a(i, i)
            a(j, i) = (a(j, i) - sum_al) / a(i, i)
        end do

        sum_ad = 0.D0 ! 初期化
        do k = 1, j-1
            sum_ad = sum_ad + a(j, k) * a(k, k) * a(k, j)
        end do
        a(j, j) = a(j, j) - sum_ad
    end do

    ! 未知ベクトルを求めていく ---------------
    y(1:N) = 0D0 ! 初期化
    ! 前進代入
    y(1) = b(1) ! 1行目
    do i = 2, N
        temp = 0.D0 ! 初期化
        do j = 1, i-1
            temp = temp + a(i, j) * y(j)
        end do  
        y(i) = b(i) - temp
    end do
    
    ! 対角項
    do i = 1, N
        z(i) = y(i) / a(i, i)
    end do

    ! 後退消去
    x(N) = z(N) ! n行目
    do j = N, 2, -1 ! 逆順のインクリメント．3つ目がstep size
        do k = 1, j-1
            z(k) = z(k) - a(k, j) * x(j)
        end do
        x(j-1) = z(j-1)
    end do


    return
end

! ---------------------------------------------------------------------
!   calculate strain and stress
! ---------------------------------------------------------------------
subroutine calc_stress(lelem, strain, stress, b, d, u, lnods)
    implicit none

    integer :: lelem, lnods(3, *), m, inode
    double precision :: strain(3, *), stress(3, *), b(3, 6), d(3, 3), u(*), ue(6)
    ! 対象要素の節点変位を持ってくる
    ue(1:6) = 0.D0 ! 初期化
    do m = 1, 3
        inode = lnods(m, lelem)
        ue(2*m-1) = u(2*inode-1)
        ue(2*m) = u(2*inode)
    end do
    ! ひずみ-変位関係式でひずみを計算
    strain(1:3, lelem) = matmul(b, ue)
    ! 応力-ひずみ間形式で応力を計算
    stress(1:3, lelem) = matmul(d, matmul(b, ue))
    
    return 
end

! ---------------------------------------------------------------------
!   
! ---------------------------------------------------------------------
subroutine xe_c_posi(lelem, lnods, x, y, u, xe_c)
    implicit none

    integer :: lelem, lnods(3, *), m, inode
    double precision :: x(*), y(*), u(*), xe(2, 3), xe_c(2)
    ! 配列xeに対象要素の変形後の座標を格納
    do m = 1, 3
        inode = lnods(m, lelem)
        xe(1, m) = x(inode) + u(2*inode-1)
        xe(2, m) = y(inode) + u(2*inode)
    end do
    ! 重心の座標
    xe_c(1) = (xe(1, 1) + xe(1, 2) + xe(1, 3)) / 3
    xe_c(2) = (xe(2, 1) + xe(2, 2) + xe(2, 3)) / 3

    return 
end
